// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// demod_2400x.c: 2.4MHz Mode S demodulator working on I/Q samples
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.
//
// This file is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "dump1090.h"

#ifdef DEMOD_2400X

#include <assert.h>

#include "dsp/helpers/demod2400x_dsp.h"
#include "dsp/helpers/demod2400x_coeffs.h"
#ifdef DEMOD_2400X_COHERENT
// Coherent second pass (demod_2400x_coherent.c): frames whose magnitude bits fail
// the CRC are decoded again coherently. Statistics are printed at exit with COHERENT_STATS
// (settings in demod_2400x_coherent.h).
#include "demod_2400x_coherent.h"

#include <time.h>

// Coherent second pass statistics (experiment)
static struct {
    unsigned long attempts, coherent, rescued, rescued_known, rescued_unknown;
    unsigned long mixed[2];                 // rescued with one half from the magnitude bits
    unsigned long rescued_long;             // rescued 112-bit frames
    unsigned long diag_frames, diag_mag_err[4], diag_coh_err[4], diag_coh_worse, diag_coh_better;
    double seconds;
    unsigned long q_attempts[17], q_coherent[17], q_rescued[17];    // by magnitude quality / 16
} coh_stats;

// Only timed for the COHERENT_STATS report: reading the clock for every frame costs CPU
static double cohNow(void)
{
    if (!COHERENT_STATS)
        return 0;
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static void cohReport(void)
{
    fprintf(stderr, "coherent: quality gate %u, rho2 gate %.2f\n", COHERENT_GATE, COHERENT_RHO2 / 1024.0);
    fprintf(stderr, "coherent: %lu attempts, %lu coherent (rho2), %lu rescued (%lu known aircraft, %lu unknown), %.0f ms\n",
            coh_stats.attempts, coh_stats.coherent, coh_stats.rescued, coh_stats.rescued_known, coh_stats.rescued_unknown,
            coh_stats.seconds * 1e3);
    fprintf(stderr, "coherent: rescued %lu long frames; %lu with the magnitude second half, %lu with the magnitude first half\n",
            coh_stats.rescued_long, coh_stats.mixed[0], coh_stats.mixed[1]);
    fprintf(stderr, "coherent by magnitude quality: quality attempts coherent rescued\n");
    for (int b = 0; b < 17; b++)
        if (coh_stats.q_attempts[b])
            fprintf(stderr, "  %3d-%3d %8lu %7lu %6lu\n", 16 * b, 16 * b + 15,
                    coh_stats.q_attempts[b], coh_stats.q_coherent[b], coh_stats.q_rescued[b]);
    if (coh_stats.diag_frames)
        fprintf(stderr, "coherent diag: %lu accepted frames; bit errors 0/1/2/3+: magnitude %lu/%lu/%lu/%lu, "
                "coherent %lu/%lu/%lu/%lu; coherent better %lu, worse %lu\n",
                coh_stats.diag_frames,
                coh_stats.diag_mag_err[0], coh_stats.diag_mag_err[1], coh_stats.diag_mag_err[2], coh_stats.diag_mag_err[3],
                coh_stats.diag_coh_err[0], coh_stats.diag_coh_err[1], coh_stats.diag_coh_err[2], coh_stats.diag_coh_err[3],
                coh_stats.diag_coh_better, coh_stats.diag_coh_worse);
}

static int scoreIsKnown(int score)
{
    switch (score) {
    case SR_UNRELIABLE_KNOWN: case SR_DF18_2ERROR_KNOWN: case SR_DF17_2ERROR_KNOWN:
    case SR_DF11_IID_1ERROR_KNOWN: case SR_DF11_ACQ_1ERROR_KNOWN: case SR_DF11_IID_KNOWN:
    case SR_DF18_1ERROR_KNOWN: case SR_DF17_1ERROR_KNOWN: case SR_DF11_ACQ_KNOWN:
    case SR_DF18_KNOWN: case SR_DF17_KNOWN:
        return 1;
    default:
        return 0;
    }
}

static unsigned bitErrors(const uint8_t *a, const uint8_t *b, int bits)
{
    unsigned n = 0;
    for (int i = 0; i < bits / 8; i++)
        n += (unsigned) __builtin_popcount(a[i] ^ b[i]);
    return n;
}
#endif /* DEMOD_2400X_COHERENT */

// New samples per pass. Each FIFO buffer (about 131k samples) is demodulated in
// pieces this size, so the 12 MHz data (5 bytes per input sample) stays in the
// L2 cache between the upsample, the preamble search and the decodes. Only this
// demodulator: the FIFO, the SDR drivers and the stock demodulator are unchanged.
#define DEMOD_BLOCK 32768

// The metric index is this many 12 MHz samples (3 input samples, the FIR delay)
// before the true start of the preamble
#define UPSAMPLE_DELAY 15

// The 12 MHz samples are 8-bit log powers (demod_log8()), whatever USE_MAGNITUDE says
#ifndef PEAMBLE_THRES
#define PEAMBLE_THRES 86      // (12 L(minp) - sum L(gaps)) / 2, see demod2400x_dsp.h
#endif
#ifndef DATA_THRES
#define DATA_THRES    65      // decode_bits() quality, see demod2400x_dsp.h
#endif
#define DATA2_THRES   ((4*DATA_THRES)/5)


// Static data
struct staticdata_t {

    int usize;
    uint8_t *udata;

    int msize;
    uint8_t *metric;

    unsigned prev_mlen;     // new samples in the previous buffer (0 = no previous buffer)
    int next_start;         // first peak position to search in this buffer (skips the rest of a
                            // frame decoded at the end of the previous buffer)

    // Adaptive gain input, as demod_2400.c: 2.4 MHz magnitudes (0..65535) of this piece
    // (from getMagnitude(), shared with the stock demodulators) and the first sample not
    // yet passed to adaptive_update()
    uint16_t *amag;         // NULL when adaptive gain is off
    unsigned adaptive_end;
};

// One preamble candidate
struct candidate {
    int pos;                            // metric index = 12 MHz sample index from mag->data[0]
    unsigned quality;                   // data quality of the first 56 bits (0..176, decode_bits())
    unsigned quality2;                  // the same for the second 56 bits (0..256 if from the coherent pass)
    uint8_t msg[MODES_LONG_MSG_BYTES];  // demodulated frame
};


#ifdef DEMOD_2400X_COHERENT
// Sample dump for offline algorithm development (COHERENT_SAMPLES). One record per
// frame, little-endian:
//   int32 accepted (1: CRC-checked frame, 0: candidate that failed)
//   int32 pos (12 MHz position chosen by the magnitude decoder)
//   uint32 quality, quality2 (magnitude, first and second 56 bits)
//   int32 score; uint32 addr (0 if not accepted)
//   uint8 mag_msg[14] (magnitude bits), uint8 good_msg[14] (corrected bits if accepted)
//   int32 chips[240][2]: complex filter outputs (I, Q) at pos + 6k, before rounding;
//   k = 0..15 preamble, 16 + 2i and 17 + 2i the two chips of bit i
// Each chip is checked against the log power the decoder used (same rounding and log).
static struct { FILE *f; unsigned gate; unsigned long written, accepted, rejected, checked, mismatches; } dump;

static void dumpReport(void)
{
    fprintf(stderr, "samples: %lu records (%lu accepted, %lu failed candidates); "
            "%lu chips checked against udata, %lu mismatches\n",
            dump.written, dump.accepted, dump.rejected, dump.checked, dump.mismatches);
    if (dump.f)
        fclose(dump.f);
}

static void dumpFrame(const int16_t *iq, unsigned iq_len, const struct candidate *c, int accepted,
                      int score, uint32_t addr, const uint8_t *good_msg, const uint8_t *udata, int nt)
{
    static int init;
    if (!init) {
        const char *s = COHERENT_SAMPLES;
        init = 1;
        if (s && (dump.f = fopen(s, "wb")) == NULL)
            perror(s);
        dump.gate = COHERENT_SAMPLES_GATE;
        if (dump.f)
            atexit(dumpReport);
    }
    if (!dump.f || (!accepted && c->quality < dump.gate))
        return;

    int32_t chips[COHERENT_CHIPS][2];
    if (coherentChips(iq, iq_len, c->pos, chips, COHERENT_CHIPS) < 0)
        return;

    // Same rounding and log as upsample12: must equal what the magnitude decoder saw
    const int stride = NU_STRIDE(nt);
    for (int k = 0; k < COHERENT_CHIPS; k++) {
        const int t = c->pos + 6 * k;
        const int32_t si = (chips[k][0] + UPSAMPLE_ROUND) >> UPSAMPLE_SHIFT;
        const int32_t sq = (chips[k][1] + UPSAMPLE_ROUND) >> UPSAMPLE_SHIFT;
        const unsigned m = demod_log8((uint32_t) (si * si) + (uint32_t) (sq * sq));
        dump.checked++;
        if (t < nt && m != udata[(t % 6) * stride + t / 6])
            dump.mismatches++;
    }

    const int32_t hdr[6] = { accepted, c->pos, (int32_t) c->quality, (int32_t) c->quality2, score, (int32_t) addr };
    static const uint8_t zeros[14];
    fwrite(hdr, sizeof(hdr), 1, dump.f);
    fwrite(c->msg, 1, 14, dump.f);
    fwrite(good_msg ? good_msg : zeros, 1, 14, dump.f);
    fwrite(chips, sizeof(chips), 1, dump.f);
    dump.written++;
    if (accepted) dump.accepted++; else dump.rejected++;
}
#endif /* DEMOD_2400X_COHERENT */

// Allocate holding buffers
static int allocBuffer(int nsize, int *osize, void **bptr) {

   void *nptr;
   if (nsize > *osize) {

       if ((nptr = realloc(*bptr, nsize)) == NULL) {

          free(*bptr);
          *bptr = NULL;
          *osize = 0;
          return -1;
       }

       *bptr = nptr;
       *osize = nsize;
   }

   return 0;
}


// Allocate buffers
static int allocateBuffers(struct staticdata_t *stdata, int totalLength) {

   int nt = NU_SAMP(totalLength);
   int stride = NU_STRIDE(nt);
   int mt = METRIC_SIZE(nt);

   if (allocBuffer(6*stride*sizeof(uint8_t), &stdata->usize, (void **)&stdata->udata) < 0) return -1;
   if (allocBuffer(mt*sizeof(uint8_t), &stdata->msize, (void **)&stdata->metric) < 0) return -1;

   return 0;
}


//#define CHECK_METRIC
#ifdef CHECK_METRIC
static void checkMetric(const uint8_t *metric, int mlen) {

    unsigned max_noise = 0;
    unsigned sum_noise = 0;
    for (int i = 0; i < mlen>>1; i++) {

        if (metric[i] > max_noise)
            max_noise = metric[i];

        sum_noise += metric[i];
    }
    unsigned peak1 = 0, peak2 = 0;
    int pos1 = 0, pos2 = 0;

    /* Find largest peak */
    for (int i = mlen >> 1; i < mlen; i++) {
        if (metric[i] > peak1) {
            peak1 = metric[i];
            pos1 = i;
        }
    }

    /* Find largest peak sufficiently far away */
    #define SEP 6
    for (int i = mlen >> 1; i < mlen; i++) {
        if ((i > pos1 ? i - pos1 : pos1 - i) < SEP)
            continue;
        if (metric[i] > peak2) {
            peak2 = metric[i];
            pos2 = i;
        }
    }

    fprintf(stderr, "noise=%u avg=%u p1=%u@%d p2=%u@%d\n"
            "nr=%.2f fr=%.2f\n",
            max_noise,
            sum_noise/(unsigned)(mlen>>1),
            peak1, pos1,
            peak2, pos2,
            (double)peak1/max_noise,
            (double)peak1/peak2);
}
#endif


// Decode 56 bits starting at metric position 'pos'; returns a quality score (0..176)
static inline __attribute__((always_inline))
unsigned decodeShortAdbsData(int pos, const uint8_t * restrict udata, int nt, uint8_t * restrict msg)
{
    unsigned quality;
    starch_decode_bits(&udata[(pos % 6) * NU_STRIDE(nt) + pos / 6 + 16], SHORT_FRAME / 8, msg, &quality);
    return quality;
}


// Is metric position 'p' a preamble peak?
static inline __attribute__((always_inline)) int isPeak(const uint8_t *metric, int p)
{
    unsigned m0 = metric[p - 1], m1 = metric[p], m2 = metric[p + 1];

    return m1 >= PEAMBLE_THRES && m1 >= m0 && m1 > m2;
}


// Does any of the 8 metric bytes at 'm' reach PEAMBLE_THRES? (exact: bytes below
// 128 get 128 - THRES added without carrying into the next byte; bytes from 128 up
// have their top bit set already)
_Static_assert(PEAMBLE_THRES >= 1 && PEAMBLE_THRES <= 128, "anyAboveThres8() needs 1 <= PEAMBLE_THRES <= 128");
static inline __attribute__((always_inline)) int anyAboveThres8(const uint8_t *m)
{
    const uint64_t ones = 0x0101010101010101ULL;
    uint64_t x;
    memcpy(&x, m, sizeof(x));
    uint64_t y = (x & (0x7f * ones)) + (128 - PEAMBLE_THRES) * ones;
    return ((x | y) & (0x80 * ones)) != 0;
}


// Try the peak and its two neighbours, keep the one whose first 56 bits look best
static void bestPhase(int peak, const uint8_t *udata, int nt, struct candidate *c)
{
    uint8_t msg[MODES_SHORT_MSG_BYTES];

    for (int p = peak - 1; p <= peak + 1; p++) {
        unsigned q = decodeShortAdbsData(p, udata, nt, msg);
        if (p == peak - 1 || q > c->quality) {
            c->pos = p;
            c->quality = q;
            memcpy(c->msg, msg, sizeof(msg));
        }
    }
}


// Demodulate the second 56 bits whatever the DF says: scoreModesMessage() can
// correct the DF field, including the bit that selects 56 or 112 bits, so the
// frame length is only known after scoring
static void demodFrame(struct candidate *c, const uint8_t *udata, int nt)
{
    c->quality2 = decodeShortAdbsData(c->pos + 12 * SHORT_FRAME, udata, nt, c->msg + MODES_SHORT_MSG_BYTES);
}


// Signal power over the frame's data bits, as a fraction of full scale
// (same scale as mag->mean_power). *energy gets the total, for the noise estimate.
static double signalPower(const int16_t *iq, int pos, int bits, double *energy)
{
    int start = (pos + UPSAMPLE_DELAY) / 5 + 19;    // skip the 8 us preamble (19.2 samples)
    int len = bits * 12 / 5;                        // 2.4 samples per bit
    uint64_t sum = 0;

    for (int k = start; k < start + len; k++) {
        int32_t I = iq[2 * k], Q = iq[2 * k + 1];
        uint32_t p = (uint32_t) (I * I + Q * Q);
        sum += p < 2048 * 2048 ? p : 2048 * 2048;    // full scale, as the stock magnitude clamps
    }

    *energy = sum / 2048.0 / 2048.0;
    Modes.stats_current.signal_power_sum += *energy;
    Modes.stats_current.signal_power_count += len;

    return *energy / len;
}


// Arrival time: the 12 MHz sample within +-2 of 'pos' where the frame's data chips
// (from the corrected message) correlate best with the upsampled signal, i.e. the
// chip centres. bestPhase() picks 'pos' by bit quality, which is flat across about
// 3 samples on a clean pulse, so 'pos' alone is a poor time estimate for MLAT.
// The correlation is on magnitudes, 2^(L/16), so weak noisy chips don't count as much as
// strong ones (once per accepted frame, so a table is fine)
static int fineTiming(int pos, const uint8_t *msg, int bits, const uint8_t *udata, int nt)
{
    static const int order[5] = { 0, -1, 1, -2, 2 };   // ties go to the nearest to 'pos'
    static uint16_t mag[256];
    if (!mag[255])
        for (int l = 0; l < 256; l++)
            mag[l] = (uint16_t) lrint(exp2(l / 16.0));
    int stride = NU_STRIDE(nt);
    int best = pos;
    int32_t best_corr = INT32_MIN;

    if (pos < 2 || pos + 2 + 6 * (16 + 2 * bits) > nt)
        return pos;

    for (int k = 0; k < 5; k++) {
        int t = pos + order[k];
        const uint8_t *u = &udata[(t % 6) * stride + t / 6 + 16];     // first data chip
        int32_t corr = 0;

        for (int i = 0; i < bits; i++) {
            int d = mag[u[2 * i]] - mag[u[2 * i + 1]];                  // first chip - second
            corr += ((msg[i >> 3] >> (7 - (i & 7))) & 1) ? d : -d;
        }
        if (corr > best_corr) {
            best_corr = corr;
            best = t;
        }
    }
    return best;
}


#ifdef DEMOD_2400X_COHERENT
// Experiment (sic-dev branch): build switches, on by default here
#ifndef COHERENT_FIRST
#define COHERENT_FIRST 1        // 1: coherent bits first for every candidate past the gate, magnitude bits the fallback
#endif
#ifndef COHERENT_SIC
#define COHERENT_SIC 1          // 1: subtract each accepted coherent frame and search its span again (float prototype)
#endif

// What the canceller needs from emitFrame(): the accepted frame and its carrier frequency
static float sic_f;             // coherent kernel's frequency estimate (Hz), if sic_fok
static int sic_fok;
static uint8_t sic_lastmsg[MODES_LONG_MSG_BYTES];
static int sic_lastbits;
static int sicIsDup(int pos, const uint8_t *msg, int bits);
static struct { unsigned long cancelled, nocarrier, bounds, found, dups, accepted; } sic_stats;
#endif

// Score, decode and pass one frame on to the rest of dump1090.
// Returns its length in bits if it was accepted, else 0; *energy gets its signal energy.
static int emitFrame(const struct mag_buf *mag, const struct candidate *c,
                     const int16_t *iq, struct staticdata_t *st, int nt, double *energy)
{
    static const struct modesMessage zeroMessage;
    struct modesMessage mm;

    int score = scoreModesMessage(c->msg);
#ifdef DEMOD_2400X_COHERENT
    sic_fok = 0;
#endif

#ifdef DEMOD_2400X_COHERENT
    const struct candidate *mag_c = c;      // the magnitude decoder's candidate
    struct candidate cc;                    // its coherent version, if that decodes
    int repaired = 0;                       // bits flipped by soft-decision repair
    struct candidate coh_c;                 // the coherent bits before any repair, if coherent
    int coh_ok = 0;
    uint16_t rel[MODES_LONG_MSG_BITS];      // their reliabilities
    uint16_t mrel[MODES_LONG_MSG_BITS];     // the magnitude bits' (quality >= COHERENT_REPAIR_MAG_GATE)
    struct repairBits crb, mrb;             // what the repairs work out from each, shared by them

    // Coherent second pass for frames the magnitude bits don't decode. Only worth
    // trying for candidates whose bits already look fairly clean (quality gate);
    // the kernel itself stops early when the carrier isn't coherent (rho2 gate).
    // COHERENT_FIRST: the coherent bits first, the magnitude bits are the fallback (and the
    // repairs only run when both fail); otherwise only for frames the magnitude bits don't decode
    const int mag_ok = COHERENT_FIRST && score >= SR_ACCEPT_THRESHOLD;
    if (Modes.coherent && (COHERENT_FIRST || score < SR_ACCEPT_THRESHOLD) && c->quality >= COHERENT_GATE) {
        static int registered;
        if (!registered) {
            registered = 1;
            if (COHERENT_STATS)
                atexit(cohReport);
        }
        const double t0 = cohNow();
        unsigned q1, rho2, rho2b;
        cc = *c;
        const int qb = c->quality / 16 < 16 ? c->quality / 16 : 16;
        Modes.stats_current.demod_coherent_tried++;
        if (coherentDecode(iq, mag->validLength, c->pos, cc.msg, rel, &q1, &cc.quality2, &rho2, &rho2b) == 0) {
            coh_stats.attempts++;
            coh_stats.q_attempts[qb]++;
        }
        else
            rho2 = 0;
        if (rho2 >= COHERENT_RHO2) {
            extern float coh_last_freq_hz;
            sic_f = coh_last_freq_hz;
            sic_fok = 1;
        }
        if (rho2 >= COHERENT_RHO2) {
            Modes.stats_current.demod_coherent_coherent++;
            coh_stats.coherent++;
            coh_stats.q_coherent[qb]++;
            coh_c = cc;
            coh_ok = 1;
            crb.rel = rel;
            crb.ready = 0;
            int s2 = scoreModesMessage(cc.msg);
            if (COHERENT_MIX && s2 < SR_ACCEPT_THRESHOLD) {
                // Each half from the other decoder in turn: the coherent first half with the
                // magnitude second, then the magnitude first with the coherent second
                for (int h = 0; h < 2 && s2 < SR_ACCEPT_THRESHOLD; h++) {
                    struct candidate mix = cc;
                    const int off = h ? 0 : MODES_SHORT_MSG_BYTES;
                    memcpy(mix.msg + off, c->msg + off, MODES_SHORT_MSG_BYTES);
                    if (!h)
                        mix.quality2 = c->quality2;
                    const int sm = scoreModesMessage(mix.msg);
                    if (sm >= SR_ACCEPT_THRESHOLD) {
                        cc = mix;
                        s2 = sm;
                        coh_stats.mixed[h]++;
                    }
                }
            }
            // Soft-decision repair: the least reliable bits, DF17/18 only, and only for an
            // aircraft already seen (the CRC is all that checks the flipped bits)
            if (Modes.coherent_repair && !mag_ok && s2 < SR_ACCEPT_THRESHOLD) {
                struct candidate rep = cc;
                const int flips = softRepair(rep.msg, &crb);
                if (flips) {
                    const int sr = scoreModesMessage(rep.msg);
                    if (sr >= SR_ACCEPT_THRESHOLD && scoreIsKnown(sr)) {
                        cc = rep;
                        s2 = sr;
                        repaired = flips;
                        Modes.stats_current.demod_coherent_repaired++;
                    }
                }
            }
            if (s2 >= SR_ACCEPT_THRESHOLD) {
                Modes.stats_current.demod_coherent_rescued++;
                coh_stats.rescued++;
                if (cc.msg[0] & 0x80)
                    coh_stats.rescued_long++;
                coh_stats.q_rescued[qb]++;
                if (scoreIsKnown(s2))
                    coh_stats.rescued_known++;
                else
                    coh_stats.rescued_unknown++;
                c = &cc;
                score = s2;
            }
        }
        coh_stats.seconds += cohNow() - t0;
    }

    // The same repair from the magnitude bits, whose reliability is the log power
    // contrast of each bit's two chips: for transmitters without a steady carrier,
    // and coherent frames the coherent repair couldn't fix. Only fairly clean frames:
    // below the gate most candidates are noise, and the search costs CPU
    if (Modes.coherent && Modes.coherent_repair && score < SR_ACCEPT_THRESHOLD &&
        mag_c->quality >= COHERENT_REPAIR_MAG_GATE) {
        const uint8_t *chip = &st->udata[(mag_c->pos % 6) * NU_STRIDE(nt) + mag_c->pos / 6 + 16];
        for (int i = 0; i < MODES_LONG_MSG_BITS; i++) {
            const int d = chip[2 * i + 1] - chip[2 * i];
            mrel[i] = (uint16_t) (d < 0 ? -d : d);
        }
        mrb.rel = mrel;
        mrb.ready = 0;
        struct candidate rep = *mag_c;
        const int flips = softRepair(rep.msg, &mrb);
        if (flips) {
            const int sr = scoreModesMessage(rep.msg);
            if (sr >= SR_ACCEPT_THRESHOLD && scoreIsKnown(sr)) {
                cc = rep;
                c = &cc;
                score = sr;
                repaired = flips;
                Modes.stats_current.demod_magnitude_repaired++;
            }
        }
    }

    // Last, assume the frame is a DF17 from a tracked aircraft: its DF and address
    // bits are known, which leaves the repair for the rest (coherent bits first)
    for (int p = 0; p < 2 && Modes.coherent && Modes.coherent_repair && score < SR_ACCEPT_THRESHOLD; p++) {
        if (p == 0 ? !coh_ok : mag_c->quality < COHERENT_REPAIR_MAG_GATE)
            continue;
        struct candidate rep = p == 0 ? coh_c : *mag_c;
        // Fewer flips from the magnitude bits: their reliabilities are rougher, and a
        // smaller search makes a chance match of the CRC less likely
        const int changed = aprioriRepair(rep.msg, p == 0 ? &crb : &mrb,
                                          p == 0 ? COHERENT_REPAIR_FLIPS : COHERENT_APRIORI_MAG_FLIPS);
        if (changed) {
            const int sr = scoreModesMessage(rep.msg);
            if (sr >= SR_ACCEPT_THRESHOLD && scoreIsKnown(sr)) {
                cc = rep;
                c = &cc;
                score = sr;
                repaired = changed;
                Modes.stats_current.demod_apriori_repaired++;
            }
        }
    }
#endif /* DEMOD_2400X_COHERENT */

    if (score < SR_ACCEPT_THRESHOLD) {
#ifdef DEMOD_2400X_COHERENT
        dumpFrame(iq, mag->validLength, c, 0, score, 0, NULL, st->udata, nt);
#endif
        if (score >= SR_UNKNOWN_THRESHOLD)
            Modes.stats_current.demod_rejected_unknown_icao++;
        else
            Modes.stats_current.demod_rejected_bad++;
        return 0;
    }

    mm = zeroMessage;
    mm.score = score;

    if (decodeModesMessage(&mm, c->msg) < 0) {
        Modes.stats_current.demod_rejected_bad++;
        return 0;
    }
#ifdef DEMOD_2400X_COHERENT
    // What's left of a frame just subtracted can decode again: drop it
    if (COHERENT_SIC && sicIsDup(c->pos, mm.msg, mm.msgbits)) {
        sic_stats.dups++;
        return 0;
    }
#endif
#ifdef DEMOD_2400X_COHERENT
    // A soft-decision repair is a correction the CRC alone vouches for: not reliable,
    // like a message whose bit errors were repaired
    if (repaired)
        mm.reliable = 0;
#endif

    // A long frame without a CRC to check (DF16/20/21/24: address/parity only)
    // also needs a clean second half
    if (mm.msgbits == MODES_LONG_MSG_BITS && mm.msgtype != 17 && mm.msgtype != 18 &&
        c->quality2 < DATA2_THRES) {
        Modes.stats_current.demod_rejected_bad++;
        return 0;
    }
    Modes.stats_current.demod_accepted[mm.correctedbits]++;

#ifdef DEMOD_2400X_COHERENT
    if (c == mag_c)
        dumpFrame(iq, mag->validLength, c, 1, score, mm.addr, mm.msg, st->udata, nt);

    // Diagnostic: bit errors of the magnitude and coherent decoders on frames whose
    // correct bits are known (accepted by the magnitude bits alone)
    if (COHERENT_DIAG && c == mag_c) {
        static int registered;
        if (!registered) {
            registered = 1;
            atexit(cohReport);
        }
        uint8_t cmsg[MODES_LONG_MSG_BYTES];
        uint16_t crel[MODES_LONG_MSG_BITS];
        unsigned q1, q2, rho2, rho2b;
        if (coherentDecode(iq, mag->validLength, c->pos, cmsg, crel, &q1, &q2, &rho2, &rho2b) == 0 && rho2 >= COHERENT_RHO2) {
            const unsigned em = bitErrors(c->msg, mm.msg, mm.msgbits), ec = bitErrors(cmsg, mm.msg, mm.msgbits);
            coh_stats.diag_frames++;
            coh_stats.diag_mag_err[em < 3 ? em : 3]++;
            coh_stats.diag_coh_err[ec < 3 ? ec : 3]++;
            if (COHERENT_DUMP) {
                extern float coh_last_freq_hz, coh_last_amp;
                fprintf(stderr, "D %.0f %.0f %u %u %u %d %06x\n", coh_last_freq_hz, coh_last_amp, em, ec,
                        c->quality, mm.msgtype, mm.addr);
            }
            if (ec < em) coh_stats.diag_coh_better++;
            if (ec > em) coh_stats.diag_coh_worse++;
        }
    }
#endif /* DEMOD_2400X_COHERENT */

    // Timestamp at the end of bit 56, even for 112-bit frames (Beast / Radarcape convention)
    int tpos = fineTiming(c->pos, mm.msg, mm.msgbits, st->udata, nt);
    mm.timestampMsg = mag->sampleTimestamp + tpos + UPSAMPLE_DELAY + (8 + 56) * 12;
    mm.sysTimestampMsg = mag->sysTimestamp + receiveclock_ms_elapsed(mag->sampleTimestamp, mm.timestampMsg);

    mm.signalLevel = signalPower(iq, c->pos, mm.msgbits, energy);
    if (mm.signalLevel > Modes.stats_current.peak_signal_power)
        Modes.stats_current.peak_signal_power = mm.signalLevel;
    if (mm.signalLevel > 0.50119)
        Modes.stats_current.strong_signal_count++;      // above -3 dBFS

    // Adaptive gain: the samples since the last frame, then this frame's (preamble included)
    if (st->amag) {
        unsigned j = (unsigned) (c->pos + UPSAMPLE_DELAY) / 5;
        if (j > st->adaptive_end)
            adaptive_update(&st->amag[st->adaptive_end], j - st->adaptive_end, NULL);
        st->adaptive_end = j + (mm.msgbits + 8) * 12 / 5;
        adaptive_update(&st->amag[j], st->adaptive_end - j, &mm);
    }

#ifdef DEMOD_2400X_COHERENT
    memcpy(sic_lastmsg, mm.msg, MODES_LONG_MSG_BYTES);
    sic_lastbits = mm.msgbits;
    sic_stats.accepted++;
#endif
    useModesMessage(&mm);
    return mm.msgbits;
}


#ifdef DEMOD_2400X_COHERENT
// ---- SIC prototype (float; experiment only) ----------------------------------
// After an accepted frame with a coherent carrier: rebuild it at 12 MHz from its
// bits, the coherent kernel's frequency and a least-squares fit, subtract it from
// the complex 12 MHz signal, rewrite udata and the preamble metric over its span,
// and search again from up to one long frame before it.
#define SIC_PAD 12
#define SIC_MAXN (6 * (16 + 2 * MODES_LONG_MSG_BITS) + 2 * SIC_PAD + 16)
#define SIC_NCOL 5

static const float sic_h[19] = { -0.022f, -0.056f, 0.0f, 0.110f, 0.267f, 0.458f, 0.652f, 0.780f, 0.930f, 1.0f,
                                 0.974f, 0.857f, 0.684f, 0.484f, 0.267f, 0.078f, -0.057f, -0.127f, -0.127f };
static struct { int pos; int bytes; uint8_t msg[MODES_LONG_MSG_BYTES]; } sic_ring[16];
static unsigned sic_ring_n;
static int sic_floor, sic_skipto;    // search restarts: never before sic_floor; end of the usual skip

static void sicReport(void)
{
    fprintf(stderr, "sic: %lu accepted, %lu subtracted, %lu without a coherent carrier, %lu out of bounds, "
            "%lu found inside or before a decoded frame, %lu duplicates dropped\n",
            sic_stats.accepted, sic_stats.cancelled, sic_stats.nocarrier, sic_stats.bounds,
            sic_stats.found, sic_stats.dups);
}

static int sicIsDup(int pos, const uint8_t *msg, int bits)
{
    unsigned n = sic_ring_n < 16 ? sic_ring_n : 16;
    for (unsigned i = 0; i < n; i++) {
        int d = sic_ring[i].pos - pos;
        if (d >= -24 && d <= 24 && sic_ring[i].bytes == bits / 8 && !memcmp(sic_ring[i].msg, msg, bits / 8))
            return 1;
    }
    return 0;
}

// Solve the n x n complex system G x = b (row-major, re/im arrays) by Gaussian elimination
static int sicSolve(int n, float *gr, float *gi, float *br, float *bi)
{
    for (int k = 0; k < n; k++) {
        int piv = k;
        float best = gr[k * n + k] * gr[k * n + k] + gi[k * n + k] * gi[k * n + k];
        for (int r = k + 1; r < n; r++) {
            float v = gr[r * n + k] * gr[r * n + k] + gi[r * n + k] * gi[r * n + k];
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-20f) return -1;
        if (piv != k) {
            for (int j = 0; j < n; j++) {
                float t = gr[k * n + j]; gr[k * n + j] = gr[piv * n + j]; gr[piv * n + j] = t;
                t = gi[k * n + j]; gi[k * n + j] = gi[piv * n + j]; gi[piv * n + j] = t;
            }
            float t = br[k]; br[k] = br[piv]; br[piv] = t;
            t = bi[k]; bi[k] = bi[piv]; bi[piv] = t;
        }
        for (int r = k + 1; r < n; r++) {
            // f = G[r][k] / G[k][k]
            float dr = gr[k * n + k], di = gi[k * n + k], den = dr * dr + di * di;
            float fr = (gr[r * n + k] * dr + gi[r * n + k] * di) / den;
            float fi = (gi[r * n + k] * dr - gr[r * n + k] * di) / den;
            for (int j = k; j < n; j++) {
                gr[r * n + j] -= fr * gr[k * n + j] - fi * gi[k * n + j];
                gi[r * n + j] -= fr * gi[k * n + j] + fi * gr[k * n + j];
            }
            br[r] -= fr * br[k] - fi * bi[k];
            bi[r] -= fr * bi[k] + fi * br[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        float sr = br[k], si = bi[k];
        for (int j = k + 1; j < n; j++) {
            sr -= gr[k * n + j] * br[j] - gi[k * n + j] * bi[j];
            si -= gr[k * n + j] * bi[j] + gi[k * n + j] * br[j];
        }
        float dr = gr[k * n + k], di = gi[k * n + k], den = dr * dr + di * di;
        br[k] = (sr * dr + si * di) / den;
        bi[k] = (si * dr - sr * di) / den;
    }
    return 0;
}

// Subtract the last accepted frame (sic_lastmsg) found at 12 MHz position 'pos'.
// Returns 1 if subtracted.
static int sicCancel(const struct mag_buf *mag, struct staticdata_t *st, int nt, int mlen, int pos)
{
    const int16_t *iq = (const int16_t *) mag->data;
    const int bits = sic_lastbits;
    const int nchip = 16 + 2 * bits;
    uint8_t a[16 + 2 * MODES_LONG_MSG_BITS];
    static float zr[SIC_MAXN], zi[SIC_MAXN], m[SIC_MAXN];

    memset(a, 0, 16);
    a[0] = a[2] = a[7] = a[9] = 1;
    for (int i = 0; i < bits; i++) {
        int b = (sic_lastmsg[i >> 3] >> (7 - (i & 7))) & 1;
        a[16 + 2 * i] = (uint8_t) b;
        a[17 + 2 * i] = (uint8_t) !b;
    }

    const int t0 = pos - 3 - SIC_PAD, t1 = pos + 3 + 6 * (nchip - 1) + SIC_PAD;
    if (t0 < 0 || t1 >= nt || t1 >= mlen + 90 || (unsigned) (t1 / 5 + 9) > mag->validLength) {
        sic_stats.bounds++;
        return 0;
    }
    const int n = t1 - t0 + 1;
    for (int t = t0; t <= t1; t++) {
        int32_t ch[1][2];
        if (coherentChips(iq, mag->validLength, t, ch, 1) < 0) { sic_stats.bounds++; return 0; }
        zr[t - t0] = (float) ch[0][0];
        zi[t - t0] = (float) ch[0][1];
    }

    // Position: the best of pos-3 .. pos+3 by chip contrast
    int best = pos;
    float bc = -1e30f;
    for (int dp = -3; dp <= 3; dp++) {
        float s = 0;
        for (int k = 0; k < nchip; k++) {
            int x = pos + dp + 6 * k - t0;
            float v = sqrtf(zr[x] * zr[x] + zi[x] * zi[x]);
            s += a[k] ? v : -v;
        }
        if (s > bc) { bc = s; best = dp; }
    }
    pos += best;

    // Frequency: the coherent kernel's estimate, sign and +-2 kHz checked on the pulse chips
    float bw = 0, bp = -1;
    for (int sg = -1; sg <= 1; sg += 2) {
        for (int df = -2000; df <= 2000; df += 250) {
            const float w = 2 * (float) M_PI * (sg * sic_f + df) / 12e6f;
            float sr = 0, si = 0, c = 1, s = 0;
            const float c6 = cosf(6 * w), s6 = sinf(6 * w);
            for (int k = 0; k < nchip; k++) {
                if (a[k]) {
                    int x = pos + 6 * k - t0;
                    sr += zr[x] * c + zi[x] * s;      // z * e^{-jw t}
                    si += zi[x] * c - zr[x] * s;
                }
                const float nc = c * c6 - s * s6, ns = c * s6 + s * c6;
                c = nc; s = ns;
            }
            if (sr * sr + si * si > bp) { bp = sr * sr + si * si; bw = w; }
        }
    }

    // Envelope and its derivative
    memset(m, 0, sizeof(float) * (size_t) n);
    for (int k = 0; k < nchip; k++) {
        if (!a[k]) continue;
        for (int d = -9; d <= 9; d++) {
            int x = pos + 6 * k + d - t0;
            if (x >= 0 && x < n) m[x] += sic_h[d + 9];
        }
    }
    float mmax = 0;
    for (int x = 0; x < n; x++) if (m[x] > mmax) mmax = m[x];
    if (mmax <= 0) return 0;

    // Derotate once: y = z e^{-jw(t-pos)} (rotation updated step by step, renormalised)
    static float yr[SIC_MAXN], yi[SIC_MAXN], v[SIC_NCOL][SIC_MAXN];
    {
        const float cr = cosf(bw), ci = sinf(bw);
        float er = cosf(bw * (t0 - pos)), ei = sinf(bw * (t0 - pos));
        for (int x = 0; x < n; x++) {
            yr[x] = zr[x] * er + zi[x] * ei;
            yi[x] = zi[x] * er - zr[x] * ei;
            const float nr = er * cr - ei * ci, ni = er * ci + ei * cr;
            er = nr; ei = ni;
            if ((x & 63) == 63) { const float g = 1.5f - 0.5f * (er * er + ei * ei); er *= g; ei *= g; }
        }
    }
    // Real columns: m, m s, m mn, m mn^2, dm (s: time across the frame, mn = m / max)
    for (int x = 0; x < n; x++) {
        const float s = (x - n / 2.0f) / n, mn = m[x] / mmax;
        v[0][x] = m[x];
        v[1][x] = m[x] * s;
        v[2][x] = m[x] * mn;
        v[3][x] = m[x] * mn * mn;
        v[4][x] = (x > 0 && x < n - 1) ? (m[x + 1] - m[x - 1]) * 0.5f : 0;
    }
    float gr[SIC_NCOL * SIC_NCOL], gi[SIC_NCOL * SIC_NCOL] = { 0 }, br[SIC_NCOL], bi[SIC_NCOL];
    for (int i = 0; i < SIC_NCOL; i++) {
        for (int j = i; j < SIC_NCOL; j++) {
            float sum = 0;
            for (int x = 0; x < n; x++) sum += v[i][x] * v[j][x];
            gr[i * SIC_NCOL + j] = gr[j * SIC_NCOL + i] = sum;
        }
        float sr = 0, si = 0;
        for (int x = 0; x < n; x++) { sr += v[i][x] * yr[x]; si += v[i][x] * yi[x]; }
        br[i] = sr;
        bi[i] = si;
    }
    if (sicSolve(SIC_NCOL, gr, gi, br, bi) < 0)
        return 0;

    // Subtract (model rotated back), then the same rounding and log as upsample12
    const int stride = NU_STRIDE(nt);
    {
        const float cr = cosf(bw), ci = sinf(bw);
        float er = cosf(bw * (t0 - pos)), ei = sinf(bw * (t0 - pos));
        for (int x = 0; x < n; x++) {
            float hr = 0, hi = 0;
            for (int j = 0; j < SIC_NCOL; j++) { hr += br[j] * v[j][x]; hi += bi[j] * v[j][x]; }
            const float mr = hr * er - hi * ei, mi = hr * ei + hi * er;
            const int32_t ri = (int32_t) lrintf(zr[x] - mr), qi = (int32_t) lrintf(zi[x] - mi);
            const int32_t si = (ri + UPSAMPLE_ROUND) >> UPSAMPLE_SHIFT, sq = (qi + UPSAMPLE_ROUND) >> UPSAMPLE_SHIFT;
            const int t = t0 + x;
            st->udata[(t % 6) * stride + t / 6] = demod_log8((uint32_t) (si * si) + (uint32_t) (sq * sq));
            const float nr = er * cr - ei * ci, ni = er * ci + ei * cr;
            er = nr; ei = ni;
            if ((x & 63) == 63) { const float g = 1.5f - 0.5f * (er * er + ei * ei); er *= g; ei *= g; }
        }
    }

    // The preamble metric over the span, as locate_preamble (metric[p] uses udata p .. p + 90)
    int p0 = t0 - 90 > 0 ? t0 - 90 : 0, p1 = t1 < mlen - 1 ? t1 : mlen - 1;
    for (int p = p0; p <= p1; p++) {
        const uint8_t *u = st->udata + (p % 6) * stride + p / 6;
        unsigned minp = MIN(MIN(u[0], u[2]), MIN(u[7], u[9]));
        unsigned gap = u[1] + u[3] + u[4] + u[5] + u[6] + u[8] + u[10] + u[11] + u[12] + u[13] + u[14] + u[15];
        unsigned v = 12 * minp > gap ? MIN((12 * minp - gap) >> 1, 255u) : 0;
        st->metric[p] = v < PEAMBLE_THRES ? 0 : (uint8_t) v;
    }

    unsigned r = sic_ring_n++ % 16;
    sic_ring[r].pos = pos;
    sic_ring[r].bytes = bits / 8;
    memcpy(sic_ring[r].msg, sic_lastmsg, (size_t) (bits / 8));
    sic_stats.cancelled++;
    return 1;
}
// ---- end SIC ------------------------------------------------------------------
#endif /* DEMOD_2400X_COHERENT */


// Search the metric for preambles and decode them. Only peaks at positions
// 1..limit belong to this buffer; later ones are in the next buffer's overlap
// and are searched there, so nothing is decoded twice.
// Returns the total signal energy of the accepted frames.
static double searchFrames(const struct mag_buf *mag, struct staticdata_t *st, int mlen, int nt)
{
    int limit = 5 * (int) (mag->validLength - mag->overlap);
    int last = MIN(limit, mlen - 2);
    double energy_sum = 0;
    int p;

#ifdef DEMOD_2400X_COHERENT
    static int sic_registered;
    if (COHERENT_SIC && !sic_registered) {
        sic_registered = 1;
        atexit(sicReport);
    }
    sic_ring_n = 0;
    sic_floor = MAX(1, st->next_start);
    sic_skipto = 0;
#endif
    for (p = MAX(1, st->next_start); p <= last; p++) {

        // Most of the metric is noise below the threshold: skip it 8 at a time
        while (p + 7 <= last && !anyAboveThres8(st->metric + p))
            p += 8;
        if (p > last)
            break;

        if (!isPeak(st->metric, p))
            continue;
        Modes.stats_current.demod_preambles++;

        struct candidate c;
        bestPhase(p, st->udata, nt, &c);
        if (c.quality < DATA_THRES)
            continue;

        demodFrame(&c, st->udata, nt);

        double energy;
        int bits = emitFrame(mag, &c, (const int16_t *) mag->data, st, nt, &energy);
        if (bits) {
            energy_sum += energy;
#ifdef DEMOD_2400X_COHERENT
            // Subtract the frame, then search again from up to one long frame before it:
            // frames that start under it or end under it. Restarts never go back past
            // sic_floor, which only moves forward, so the loop ends.
            if (COHERENT_SIC && Modes.coherent) {
                if (c.pos < sic_skipto)
                    sic_stats.found++;      // the usual skip would never have tried it
                if (c.pos + bits * 12 > sic_skipto)
                    sic_skipto = c.pos + bits * 12;
                if (!sic_fok)
                    sic_stats.nocarrier++;
                else if (sicCancel(mag, st, nt, mlen, c.pos)) {
                    int back = c.pos - 12 * MODES_LONG_MSG_BITS;
                    if (back < sic_floor)
                        back = sic_floor;
                    if (c.pos + 3 > sic_floor)
                        sic_floor = c.pos + 3;
                    p = back - 1;
                    continue;
                }
                if (c.pos + 3 > sic_floor)
                    sic_floor = c.pos + 3;
            }
#endif
            // Skip to 8 bits before the end of the frame, as demod_2400.c does
            p = c.pos + bits * 12 - 1;
        }
    }

    // A skip that runs past 'limit' carries on into the next buffer,
    // whose position 0 is this buffer's 'limit'
    st->next_start = MAX(1, p - limit);
    return energy_sum;
}


// Demodulate 2400
//#define CHECK_USAMP
//#define CHECK_POW
// mag->data contains DC-removed signed 12-bit I/Q sample pairs (-2048..2047): [overlap | new samples]
static void demodulateBlock(struct mag_buf *mag, struct staticdata_t *st, uint16_t *amag)
{
   unsigned overlap = mag->overlap;
   unsigned dlen = mag->validLength;            // whole buffer, overlap included
   unsigned newlen = dlen - overlap;            // new samples in this buffer

   assert(overlap >= OVERLAP_WINDOW);

   // No usable skip from the previous buffer
   if (!st->prev_mlen || (mag->flags & MAGBUF_DISCONTINUOUS))
       st->next_start = 0;

   st->prev_mlen = newlen;

   // Upsample the whole buffer
   int ulen = starch_upsample12((const int16_t *) mag->data, (int) dlen, st->udata);

#ifdef CHECK_POW
   int max = 0;
   for (int i = 0; i < 6 * NU_STRIDE(ulen); i++) {
       int mm = st->udata[i];
       if (mm > max) max = mm;
   }
   fprintf(stderr, "Max log power = %d\n", max);
#endif
#ifdef CHECK_USAMP
   int stride = NU_STRIDE(ulen);
   for (int i = 0; i < ulen; i++) {
        int phase = i % 6;
        int pos   = i / 6;
        fprintf(stderr, "%4d %u\n", i, st->udata[phase * stride + pos]);
   }
#endif

   // Locate possible preambles
   int mlen = starch_locate_preamble(st->udata, ulen, PEAMBLE_THRES, st->metric);   // metric < PEAMBLE_THRES is 0
#ifdef CHECK_METRIC
   checkMetric(st->metric, mlen);
#endif

   // Magnitudes for adaptive gain, on the stock scale
   st->amag = amag;
   if (amag) {
       if (st->adaptive_end > newlen || (mag->flags & MAGBUF_DISCONTINUOUS))
           st->adaptive_end = 0;
   }

   // Decode the frames this buffer owns
   double signal_energy = searchFrames(mag, st, mlen, ulen);

   // Adaptive gain: the rest of the new samples; a frame running into the overlap
   // carries on into the next buffer, as demod_2400.c
   if (st->amag) {
       if (st->adaptive_end < newlen) {
           adaptive_update(&st->amag[st->adaptive_end], newlen - st->adaptive_end, NULL);
           st->adaptive_end = 0;
       } else {
           st->adaptive_end -= newlen;
       }
   }

   // Noise power = everything in the new samples that wasn't a decoded frame
   Modes.stats_current.noise_power_sum += mag->mean_power * newlen - signal_energy;
   Modes.stats_current.noise_power_count += newlen;
}

void demodulate2400x(struct mag_buf *mag)
{
   static struct staticdata_t stdata;

   const unsigned overlap = mag->overlap;
   const unsigned newlen = mag->validLength - overlap;
   if (allocateBuffers(&stdata, MIN(mag->totalLength, overlap + DEMOD_BLOCK)) < 0)
       return;

   // Magnitudes of the whole buffer for adaptive gain (made once, see getMagnitude())
   uint16_t *amag = NULL;
   if ((Modes.adaptive_burst_control || Modes.adaptive_range_control) && (amag = getMagnitude(mag)) == NULL)
       return;

   // Each piece: the previous piece's last 'overlap' samples, then up to
   // DEMOD_BLOCK new ones, just as if the FIFO buffers were that size
   for (unsigned start = 0; start < newlen; start += DEMOD_BLOCK) {
       struct mag_buf piece = *mag;
       piece.data = magbuf_sample(mag, start);
       piece.validLength = overlap + MIN(DEMOD_BLOCK, newlen - start);
       piece.totalLength = MIN(mag->totalLength - start, overlap + DEMOD_BLOCK);
       piece.sampleTimestamp = mag->sampleTimestamp + start * 12e6 / Modes.sample_rate;
       piece.sysTimestamp = mag->sysTimestamp + receiveclock_ms_elapsed(mag->sampleTimestamp, piece.sampleTimestamp);
       if (start)
           piece.flags &= ~MAGBUF_DISCONTINUOUS;
       demodulateBlock(&piece, &stdata, amag ? amag + start : NULL);
   }
}

#endif /* DEMOD_2400X */
