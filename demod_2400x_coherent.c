// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// demod_2400x_coherent.c: coherent second pass for DEMOD_2400X (experimental)
//
// For a frame the magnitude decoder found, interpolate the complex I/Q at the
// chip centres (2 MHz, at the 12 MHz position it chose), estimate the carrier
// frequency offset and phase from the preamble and the first 56 bits, de-rotate
// and slice the real part. One dimension of noise instead of two, and no
// envelope detection loss.
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#include "dump1090.h"

#ifdef DEMOD_2400X

#include <math.h>
#include <complex.h>

#include "demod_2400x_coherent.h"
#include "dsp/helpers/demod2400x_dsp.h"
#include "dsp/helpers/demod2400x_coeffs.h"

#define PREAMBLE_CHIPS 16
#define DATA_CHIPS (2 * LONG_FRAME)
#define CHIPS (PREAMBLE_CHIPS + DATA_CHIPS)
#define EST_CHIPS (PREAMBLE_CHIPS + 2 * SHORT_FRAME)       // preamble + first 56 bits

// Largest carrier offset searched: +-62.5 kHz (about +-57 ppm at 1090 MHz), in
// radians per chip (0.5 us). Delay-and-multiply over 16 chips is unambiguous over
// exactly this range.
#define MAX_OFFSET_HZ 62500.0
#define CHIP_RATE 2e6
#define WMAX (2 * M_PI * MAX_OFFSET_HZ / CHIP_RATE)

float coh_last_freq_hz, coh_last_amp;     // experiment: last estimate, for diagnostics

const struct coherent_config *coherentConfig(void)
{
    static struct coherent_config cfg;
    static int init;
    if (!init) {
        const char *s;
        init = 1;
        cfg.enabled = (s = getenv("DEMOD_COHERENT")) ? atoi(s) : 0;
        cfg.estimator = (s = getenv("COH_EST")) ? atoi(s) : 0;
        cfg.refine = (s = getenv("COH_REFINE")) ? atoi(s) : 1;
        cfg.gate = (s = getenv("COH_GATE")) ? (unsigned) atoi(s) : 0;
        cfg.diag = (s = getenv("COH_DIAG")) ? atoi(s) : 0;
    }
    return &cfg;
}

// Complex filter output at 12 MHz position t: same polyphase table as upsample12
static inline float complex interpChip(const int16_t *iq, int t)
{
    const int n = t / 5, c = t % 5;
    int32_t si = 0, sq = 0;
    for (int k = 0; k < 8; k++) {
        si += upsample_coeffs[c][k] * iq[2 * (n + k)];
        sq += upsample_coeffs[c][k] * iq[2 * (n + k) + 1];
    }
    return (float) si + I * (float) sq;
}

int coherentChips(const int16_t *iq, unsigned iq_len, int pos, int32_t (*chips)[2], int nchips)
{
    if (pos < 0 || (unsigned) ((pos + 6 * (nchips - 1)) / 5 + 8) > iq_len)
        return -1;
    for (int k = 0; k < nchips; k++) {
        const int t = pos + 6 * k, n = t / 5, c = t % 5;
        int32_t si = 0, sq = 0;
        for (int j = 0; j < 8; j++) {
            si += upsample_coeffs[c][j] * iq[2 * (n + j)];
            sq += upsample_coeffs[c][j] * iq[2 * (n + j) + 1];
        }
        chips[k][0] = si;
        chips[k][1] = sq;
    }
    return 0;
}

// Frequency (radians per chip) by delay-and-multiply: on chips 16 apart (8 us)
// have the same phase apart from the rotation
static float estimateDelayMultiply(const float complex *x, int n)
{
    float complex r = 0;
    for (int k = 0; k + 16 < n; k++)
        r += x[k + 16] * conjf(x[k]);
    return cargf(r) / 16;
}

// Frequency by searching |sum x[k] e^{-jwk}| over +-WMAX, then a parabola
// through the best point and its neighbours
static float estimateDFT(const float complex *x, int n)
{
    enum { STEPS = 64 };
    const float step = 2 * WMAX / STEPS;
    float p[STEPS + 1];
    int best = 0;
    for (int s = 0; s <= STEPS; s++) {
        const float w = -WMAX + s * step;
        const float complex rot = cexpf(-I * w);
        float complex ph = 1, sum = 0;
        for (int k = 0; k < n; k++) {
            sum += x[k] * ph;
            ph *= rot;
        }
        p[s] = crealf(sum) * crealf(sum) + cimagf(sum) * cimagf(sum);
        if (p[s] > p[best])
            best = s;
    }
    float w = -WMAX + best * step;
    if (best > 0 && best < STEPS) {
        const float den = p[best - 1] - 2 * p[best] + p[best + 1];
        if (den < 0)
            w += 0.5f * step * (p[best - 1] - p[best + 1]) / den;
    }
    return w;
}

// Phase of the "on" chips after removing the rotation w
static float estimatePhase(const float complex *x, const uint8_t *on, int n, float w)
{
    float complex z = 0;
    for (int k = 0; k < n; k++)
        if (on[k])
            z += x[k] * cexpf(-I * w * k);
    return cargf(z);
}

// Slice bits from the de-rotated real parts; marks the on chips, returns the quality
static unsigned sliceBits(const float *y, int first_bit, int nbits, uint8_t *msg, uint8_t *on)
{
    float qual = 0, pow = 0;
    for (int i = first_bit; i < first_bit + nbits; i++) {
        const float a = y[PREAMBLE_CHIPS + 2 * i], b = y[PREAMBLE_CHIPS + 2 * i + 1];
        const int bit = a > b;
        if (bit)
            msg[i >> 3] |= (uint8_t) (0x80 >> (i & 7));
        else
            msg[i >> 3] &= (uint8_t) ~(0x80 >> (i & 7));
        on[PREAMBLE_CHIPS + 2 * i] = (uint8_t) bit;
        on[PREAMBLE_CHIPS + 2 * i + 1] = (uint8_t) !bit;
        qual += fabsf(a - b);
        pow += fabsf(a) + fabsf(b);
    }
    return (unsigned) (256 * qual / (pow + 1));
}

int coherentDecode(const int16_t *iq, unsigned iq_len, int pos, const uint8_t *mag_msg,
                   uint8_t *out, unsigned *q1, unsigned *q2)
{
    const struct coherent_config *cfg = coherentConfig();
    float complex x[CHIPS];
    float y[CHIPS];
    uint8_t on[CHIPS] = { 0 };

    if (pos < 0 || (unsigned) ((pos + 6 * (CHIPS - 1)) / 5 + 8) > iq_len)
        return -1;

    for (int k = 0; k < CHIPS; k++)
        x[k] = interpChip(iq, pos + 6 * k);

    // Known on chips: preamble pulses, and the first 56 bits as the magnitude decoder saw them
    on[0] = on[2] = on[7] = on[9] = 1;
    for (int i = 0; i < SHORT_FRAME; i++) {
        const int bit = (mag_msg[i >> 3] >> (7 - (i & 7))) & 1;
        on[PREAMBLE_CHIPS + 2 * i + !bit] = 1;
    }

    float w = cfg->estimator ? estimateDFT(x, EST_CHIPS) : estimateDelayMultiply(x, EST_CHIPS);
    float phi = estimatePhase(x, on, EST_CHIPS, w);

    for (int pass = 0; pass <= cfg->refine; pass++) {
        for (int k = 0; k < CHIPS; k++)
            y[k] = crealf(x[k] * cexpf(-I * (w * k + phi)));
        memset(out, 0, MODES_LONG_MSG_BYTES);
        *q1 = sliceBits(y, 0, SHORT_FRAME, out, on);
        *q2 = sliceBits(y, SHORT_FRAME, SHORT_FRAME, out, on);
        if (pass == cfg->refine)
            break;

        // Decision-directed refinement: weighted straight-line fit to the phase of
        // the on chips over the whole frame
        float sw = 0, sk = 0, st = 0, skk = 0, skt = 0;
        for (int k = 0; k < CHIPS; k++) {
            if (!on[k])
                continue;
            const float complex u = x[k] * cexpf(-I * (w * k + phi));
            const float wt = cabsf(u), th = cargf(u);
            sw += wt; sk += wt * k; st += wt * th; skk += wt * k * k; skt += wt * k * th;
        }
        const float den = sw * skk - sk * sk;
        if (sw <= 0 || den <= 0)
            break;
        const float b = (sw * skt - sk * st) / den;
        const float a = (st - b * sk) / sw;
        w += b;
        phi += a;
    }
    coh_last_freq_hz = (float) (w * CHIP_RATE / (2 * M_PI));
    float amp = 0; int non = 0;
    for (int k = 0; k < EST_CHIPS; k++) if (on[k]) { amp += cabsf(x[k]); non++; }
    coh_last_amp = non ? amp / non : 0;
    return 0;
}

#endif /* DEMOD_2400X */
