// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// demod_2400x_coherent.c: coherent second pass for DEMOD_2400X (experimental)
//
// For a frame the magnitude decoder found, interpolate the complex I/Q at the
// chip centres (2 MHz, at the 12 MHz position it chose), find the carrier
// frequency and channel (phase and amplitude) with a correlator over the first
// 56 bits, and slice the real part of all 112 bits with a fixed rotation.
//
// Every PPM bit has exactly one pulse, so the pair sum of its two chips is a
// carrier sample whatever the data. The frequency search runs on 7 block sums
// of 8 pair sums (one per byte) over 33 bins across +-62.5 kHz, refined to 1/8
// of a bin with a parabola; one exact correlation over the 56 pair sums then
// gives the channel h and the coherence rho2 = |h|^2 / (56 sum |c|^2).
// Transmitters without a stable carrier (8-PSK phase overlay, free-running or
// chirping oscillators) give a low rho2 and are left to the magnitude decoder.
//
// Integer arithmetic. The interpolation multiplies 16-bit samples by 16-bit taps
// into 32-bit sums; the chips are then normalised to 13 bits and everything else
// is 16-bit, with Q15 multiplies mulq15(a, b) = (a*b + 2^14) >> 15 (NEON vqrdmulh,
// AVX2 mulhrs). Only sums over many values (h, energies, qualities) are 32-bit,
// and the final gate comparison is 64-bit.
// Two stages: the second 56 bits are only interpolated if the gate passes.
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#include "dump1090.h"

#ifdef DEMOD_2400X

#include <math.h>

#include "demod_2400x_coherent.h"
#include "dsp/helpers/demod2400x_dsp.h"
#include "dsp/helpers/demod2400x_coeffs.h"

#define PREAMBLE_CHIPS 16
#define EST_BITS SHORT_FRAME                   // estimate over the first 56 bits
#define BLOCK 8                                // pair sums per block (one byte)
#define NBLOCKS (EST_BITS / BLOCK)
#define BINS 33                                // across +-62.5 kHz, ~3.9 kHz apart
#define FINE 8                                 // fine frequency steps per bin
#define FINE_BINS ((BINS - 1) * FINE + 1)
#define CHIP_BITS 13                           // chips normalised to |x| < 2^13
#define BLOCK_BITS 11                          // block sums normalised to |B| < 2^11

// Largest carrier offset searched: +-62.5 kHz (about +-57 ppm at 1090 MHz), in
// radians per chip (0.5 us); blocks of 8 bits (125 kHz block rate) are
// unambiguous over exactly this range
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
        cfg.gate = (s = getenv("COH_GATE")) ? (unsigned) atoi(s) : 96;
        cfg.rho2 = (s = getenv("COH_RHO2")) ? (float) atof(s) : 0.6f;
        cfg.diag = (s = getenv("COH_DIAG")) ? atoi(s) : 0;
    }
    return &cfg;
}

int coherentChips(const int16_t *iq, unsigned iq_len, int pos, int32_t (*chips)[2], int nchips)
{
    if (pos < 0 || (unsigned) ((pos + 6 * (nchips - 1)) / 5 + 8) > iq_len)
        return -1;
    for (int k = 0; k < nchips; k++) {
        const int t = pos + 6 * k, n = t / 5, c = t % 5;
        int32_t si = 0, sq = 0;                   // |s| <= 2048 * UPSAMPLE_ABS_SUM < 2^19
        for (int j = 0; j < 8; j++) {
            si += upsample_coeffs[c][j] * iq[2 * (n + j)];
            sq += upsample_coeffs[c][j] * iq[2 * (n + j) + 1];
        }
        chips[k][0] = si;
        chips[k][1] = sq;
    }
    return 0;
}

// 16-bit complex value / Q15 phasor (components within +-32767)
typedef struct { int16_t re, im; } c16;

// Q15 multiply, rounded: (a*b + 2^14) >> 15 (NEON vqrdmulh, AVX2 _mm256_mulhrs_epi16)
static inline int16_t mulq15(int16_t a, int16_t b)
{
    return (int16_t) (((int32_t) a * b + (1 << 14)) >> 15);
}

// x * p for a Q15 phasor p; |result| <= |x| (plus rounding), so 16 bits if |x| < 2^15 / sqrt(2)
static inline c16 rot16(c16 x, c16 p)
{
    c16 r = { (int16_t) (mulq15(x.re, p.re) - mulq15(x.im, p.im)),
              (int16_t) (mulq15(x.re, p.im) + mulq15(x.im, p.re)) };
    return r;
}

// Q15 phasors, built once:
//   bin_rot[j][m]: e^{-j w_j t_m}, w_j the coarse bin frequency, t_m the block centre (chips)
//   fine_chip[f], fine_bit[f]: e^{-j w_f} and e^{-j 2 w_f} for the fine frequencies
static c16 bin_rot[BINS][NBLOCKS];
static c16 fine_chip[FINE_BINS], fine_bit[FINE_BINS];

static c16 q15(double angle)
{
    c16 r = { (int16_t) lrint(32767 * cos(angle)), (int16_t) lrint(32767 * sin(angle)) };
    return r;
}

static void initTables(void)
{
    static int done;
    if (done)
        return;
    done = 1;
    for (int j = 0; j < BINS; j++) {
        const double w = -WMAX + 2 * WMAX * j / (BINS - 1);
        for (int m = 0; m < NBLOCKS; m++)
            bin_rot[j][m] = q15(-w * 2 * (m * BLOCK + (BLOCK - 1) / 2.0));
    }
    for (int f = 0; f < FINE_BINS; f++) {
        const double w = -WMAX + 2 * WMAX * f / (FINE_BINS - 1);
        fine_chip[f] = q15(-w);
        fine_bit[f] = q15(-2 * w);
    }
}

static int bitsNeeded(uint32_t v)
{
    return v ? 32 - __builtin_clz(v) : 0;
}

static inline uint32_t iabs32(int32_t v)
{
    return v < 0 ? (uint32_t) -v : (uint32_t) v;
}

// Shift that brings max |v| below 2^bits
static int scaleShift(uint32_t vmax, int bits)
{
    const int n = bitsNeeded(vmax);
    return n > bits ? n - bits : 0;
}

// Chips (int32) to 16 bits with a shared shift, clamped to |x| < 2^CHIP_BITS
static void toChips16(const int32_t (*ch)[2], c16 *out, int n, int shift)
{
    const int32_t lim = (1 << CHIP_BITS) - 1;
    for (int k = 0; k < n; k++) {
        int32_t re = ch[k][0] >> shift, im = ch[k][1] >> shift;
        out[k].re = (int16_t) (re > lim ? lim : re < -lim ? -lim : re);
        out[k].im = (int16_t) (im > lim ? lim : im < -lim ? -lim : im);
    }
}

int coherentDecode(const int16_t *iq, unsigned iq_len, int pos,
                   uint8_t *out, unsigned *q1, unsigned *q2, float *rho2_out)
{
    const struct coherent_config *cfg = coherentConfig();
    int32_t ch32[2 * EST_BITS][2];
    c16 x[2 * LONG_FRAME];                         // data chips: bit i at 2i, 2i+1

    initTables();
    *rho2_out = 0;

    // Stage 1: the first 56 bits, normalised to CHIP_BITS from their largest value
    if (coherentChips(iq, iq_len, pos + 6 * PREAMBLE_CHIPS, ch32, 2 * EST_BITS) < 0)
        return -1;
    uint32_t cmax = 0;
    for (int k = 0; k < 2 * EST_BITS; k++) {
        if (iabs32(ch32[k][0]) > cmax) cmax = iabs32(ch32[k][0]);
        if (iabs32(ch32[k][1]) > cmax) cmax = iabs32(ch32[k][1]);
    }
    const int cs = scaleShift(cmax, CHIP_BITS);
    toChips16(ch32, x, 2 * EST_BITS, cs);

    // Block sums of the pair sums (one per byte): 16 chips < 2^13 -> < 2^17 (int32),
    // then normalised to BLOCK_BITS
    int32_t bre[NBLOCKS], bim[NBLOCKS];
    uint32_t bmax = 0;
    for (int m = 0; m < NBLOCKS; m++) {
        int32_t re = 0, im = 0;
        for (int k = 2 * m * BLOCK; k < 2 * (m + 1) * BLOCK; k++) {
            re += x[k].re;
            im += x[k].im;
        }
        bre[m] = re; bim[m] = im;
        if (iabs32(re) > bmax) bmax = iabs32(re);
        if (iabs32(im) > bmax) bmax = iabs32(im);
    }
    const int bs = scaleShift(bmax, BLOCK_BITS);
    c16 blk[NBLOCKS];
    for (int m = 0; m < NBLOCKS; m++)
        blk[m] = (c16) { (int16_t) (bre[m] >> bs), (int16_t) (bim[m] >> bs) };

    // Coarse search: |sum_m B_m e^{-j w t_m}|^2 per bin; each rotated block < 2^11 * sqrt(2),
    // 7 of them < 20300 (int16); power < 2^30 (int32)
    int32_t p[BINS];
    int best = 0;
    for (int j = 0; j < BINS; j++) {
        int16_t re = 0, im = 0;
        for (int m = 0; m < NBLOCKS; m++) {
            const c16 r = rot16(blk[m], bin_rot[j][m]);
            re = (int16_t) (re + r.re);
            im = (int16_t) (im + r.im);
        }
        p[j] = re * re + im * im;
        if (p[j] > p[best])
            best = j;
    }

    // Parabola through the best bin and its neighbours, to 1/8 of a bin (powers >> 4: < 2^26)
    int f = best * FINE;
    if (best > 0 && best < BINS - 1) {
        const int32_t pm = p[best - 1] >> 4, p0 = p[best] >> 4, pp = p[best + 1] >> 4;
        const int32_t den = 2 * (2 * p0 - pm - pp);                  // >= 0 at the peak, < 2^29
        if (den > 0) {
            const int32_t num = FINE * (pp - pm);                    // |num| < 2^29
            int32_t d = (num >= 0 ? num + den / 2 : num - den / 2) / den;
            if (d > FINE / 2) d = FINE / 2;
            if (d < -FINE / 2) d = -FINE / 2;
            f += d;
        }
    }
    const c16 rc = fine_chip[f], rb = fine_bit[f];

    // Exact correlation over the 56 pair sums: c_i = x_a ph_i + x_b (rc ph_i), each rotated chip
    // < 2^13 * sqrt(2), so c_i < 23170 (int16); h < 56 * 23170 < 2^21 (int32);
    // sum |c_i >> 4|^2 < 56 * 2 * 1448^2 < 2^28 (uint32)
    int32_t hre = 0, him = 0;
    uint32_t e2 = 0;
    c16 ph = { 32767, 0 };
    for (int i = 0; i < EST_BITS; i++) {
        const c16 a = rot16(x[2 * i], ph), b = rot16(x[2 * i + 1], rot16(rc, ph));
        const int16_t cre = (int16_t) (a.re + b.re), cim = (int16_t) (a.im + b.im);
        hre += cre; him += cim;
        e2 += (uint32_t) ((cre >> 4) * (cre >> 4)) + (uint32_t) ((cim >> 4) * (cim >> 4));
        ph = rot16(ph, rb);
    }

    // Gate: rho2 = |h|^2 / (56 sum |c|^2). h >> 8 < 2^13, |h>>8|^2 < 2^27;
    // |h|^2 = 2^16 |h>>8|^2 and sum |c|^2 = 2^8 e2, so rho2 = 2^8 |h>>8|^2 / (56 e2)
    const int32_t hr8 = hre >> 8, hi8 = him >> 8;
    const uint32_t h2 = (uint32_t) (hr8 * hr8) + (uint32_t) (hi8 * hi8);
    *rho2_out = e2 ? 256.0f * h2 / (EST_BITS * (float) e2) : 0;
    const uint32_t thr = (uint32_t) lrintf(cfg->rho2 * 1024);       // rho2 in 1/1024
    if (!e2 || (uint64_t) h2 * 256 * 1024 < (uint64_t) thr * EST_BITS * e2)
        return 0;

    // Stage 2: the second 56 bits, with the same shift (clamped)
    int32_t ch2[2 * EST_BITS][2];
    if (coherentChips(iq, iq_len, pos + 6 * (PREAMBLE_CHIPS + 2 * EST_BITS), ch2, 2 * EST_BITS) < 0)
        return -1;
    toChips16(ch2, x + 2 * EST_BITS, 2 * EST_BITS, cs);

    // Slice all 112 bits: y = Re(x ph conj(g)), g = h scaled to 15 bits; with x ph < 2^13.5,
    // y = mulq15(re, g.re) + mulq15(im, g.im) < 2^14 (int16)
    const int gs = scaleShift(iabs32(hre) > iabs32(him) ? iabs32(hre) : iabs32(him), 15);
    const int32_t glim = 32767;
    int32_t gre = hre >> gs, gim = him >> gs;
    if (gre > glim) gre = glim;
    if (gim > glim) gim = glim;
    if (gre < -glim) gre = -glim;
    if (gim < -glim) gim = -glim;
    uint32_t qual[2] = { 0, 0 }, pw[2] = { 0, 0 };                  // 56 * 2^15 < 2^21
    memset(out, 0, MODES_LONG_MSG_BYTES);
    ph = (c16) { 32767, 0 };
    for (int i = 0; i < LONG_FRAME; i++) {
        const c16 a = rot16(x[2 * i], ph), b = rot16(x[2 * i + 1], rot16(rc, ph));
        const int16_t ya = (int16_t) (mulq15(a.re, (int16_t) gre) + mulq15(a.im, (int16_t) gim));
        const int16_t yb = (int16_t) (mulq15(b.re, (int16_t) gre) + mulq15(b.im, (int16_t) gim));
        if (ya > yb)
            out[i >> 3] |= (uint8_t) (0x80 >> (i & 7));
        qual[i >= SHORT_FRAME] += iabs32(ya - yb);
        pw[i >= SHORT_FRAME] += iabs32(ya) + iabs32(yb);
        ph = rot16(ph, rb);
    }
    *q1 = 256 * qual[0] / (pw[0] + 1);                              // 256 * 2^21 < 2^29
    *q2 = 256 * qual[1] / (pw[1] + 1);

    coh_last_freq_hz = (float) ((-WMAX + 2 * WMAX * f / (FINE_BINS - 1)) * CHIP_RATE / (2 * M_PI));
    coh_last_amp = 0;
    return 0;
}

#endif /* DEMOD_2400X */
