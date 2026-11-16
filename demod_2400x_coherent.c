// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// demod_2400x_coherent.c: coherent second pass for DEMOD_2400X (experimental)
//
// For a frame the magnitude decoder found but couldn't decode, decode it again
// coherently: the coherent_decode starch kernel (dsp/impl/demod2400x.*.c, method
// in dsp/helpers/demod2400x_coherent.h) estimates the carrier frequency and
// channel from the first 56 bits and slices the real part of all 112 bits.
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#include "dump1090.h"

#if defined(DEMOD_2400X) && defined(DEMOD_2400X_COHERENT)

#include "demod_2400x_coherent.h"
#include "dsp/helpers/demod2400x_coherent.h"

float coh_last_freq_hz, coh_last_amp;     // experiment: last estimate, for diagnostics

// For the chip-sample dump: chips (int32 I, Q) at pos + 6k, k = 0 .. nchips-1
int coherentChips(const int16_t *iq, unsigned iq_len, int pos, int32_t (*chips)[2], int nchips)
{
    if (pos < 0 || (unsigned) ((pos + 6 * (nchips - 1)) / 5 + 8) > iq_len)
        return -1;
    for (int k = 0; k < nchips; k++)
        coh_chip(iq, pos + 6 * k, &chips[k][0], &chips[k][1]);
    return 0;
}

// Bounds check, then the starch kernel; the rho2 gate (first half) lets the kernel
// stop early on frames that aren't coherent. rho2b: coherence of the second half
int coherentDecode(const int16_t *iq, unsigned iq_len, int pos, uint8_t *out, uint16_t *rel,
                   unsigned *q1, unsigned *q2, unsigned *rho2, unsigned *rho2b)
{
    uint32_t info[COH_INFO_N];

    if (pos < 0 || coh_iq_needed(pos) > iq_len)
        return -1;

    starch_coherent_decode(iq, pos, COHERENT_RHO2, out, info, rel);
    *rho2 = info[COH_INFO_RHO2];
    *rho2b = info[COH_INFO_RHO2B];
    *q1 = info[COH_INFO_Q1];
    *q2 = info[COH_INFO_Q2];
    coh_last_freq_hz = (float) (COH_WBIN_HZ * ((double) info[COH_INFO_WIDE] - (COH_WBINS - 1) / 2.0) +
                                COH_BIN_HZ * ((double) info[COH_INFO_FINE] / COH_FINE - (COH_BINS - 1) / 2.0));
    coh_last_amp = 0;
    return 0;
}

// CRC syndrome of a 112-bit frame that differs from a valid one only in bit b
static uint32_t bitSyndrome(int b)
{
    uint8_t m[MODES_LONG_MSG_BYTES] = { 0 };
    m[b >> 3] = (uint8_t) (0x80 >> (b & 7));
    return modesChecksum(m, MODES_LONG_MSG_BITS);
}

// Whether the DF bits in 'flip' (bit 4 = the frame's first bit) could all be among
// the COHERENT_REPAIR_BITS least reliable: each must at least be among those of the
// first half, where the raw reliabilities compare. A quick test before any scaling
static int dfFlippable(const uint16_t *rel, unsigned flip)
{
    if (__builtin_popcount(flip) > COHERENT_REPAIR_FLIPS)
        return 0;
    for (int b = 0; b < 5; b++) {
        if (!(flip & (0x10u >> b)))
            continue;
        int weaker = 0;
        for (int i = 0; i < 56; i++)
            if (rel[i] < rel[b] && ++weaker >= COHERENT_REPAIR_BITS)
                return 0;
    }
    return 1;
}

int coherentRepair(uint8_t *msg, const uint16_t *rel)
{
    // Most frames can't become a DF17/18 by flipping only unreliable bits
    const unsigned df = msg[0] >> 3;
    if (!dfFlippable(rel, df ^ 17) && !dfFlippable(rel, df ^ 18))
        return 0;

    // Each half has its own channel scaling, so a bit's reliability is |y_a - y_b|
    // over its half's sum: rel * (2^32 / sum) >> 16, one division per half
    uint32_t sum[2] = { 0, 0 };
    for (int i = 0; i < MODES_LONG_MSG_BITS; i++)
        sum[i / 56] += rel[i];
    if (!sum[0] || !sum[1])
        return 0;
    const uint64_t inv[2] = { 0xffffffffu / sum[0], 0xffffffffu / sum[1] };

    // The least reliable bits, in order
    int bit[COHERENT_REPAIR_BITS];
    uint32_t cost[COHERENT_REPAIR_BITS];
    int n = 0;
    for (int i = 0; i < MODES_LONG_MSG_BITS; i++) {
        const uint32_t c = (uint32_t) ((rel[i] * inv[i / 56]) >> 16);
        if (n == COHERENT_REPAIR_BITS && c >= cost[n - 1])
            continue;
        int j = n < COHERENT_REPAIR_BITS ? n++ : n - 1;
        for (; j > 0 && cost[j - 1] > c; j--) {
            cost[j] = cost[j - 1];
            bit[j] = bit[j - 1];
        }
        cost[j] = c;
        bit[j] = i;
    }

    // The same test exactly: the DF bits to change must all be among them
    unsigned weak_df = 0;                   // DF bits (bit 4 = first) among them
    for (int j = 0; j < n; j++)
        if (bit[j] < 5)
            weak_df |= 0x10u >> bit[j];
    if (((df ^ 17) & ~weak_df) && ((df ^ 18) & ~weak_df))
        return 0;

    uint32_t syn[COHERENT_REPAIR_BITS];
    for (int j = 0; j < n; j++)
        syn[j] = bitSyndrome(bit[j]);

    // Every combination of 1 .. COHERENT_REPAIR_FLIPS of them (indices in increasing
    // order); the syndrome is linear, so flipping bits XORs their syndromes in
    const uint32_t s0 = modesChecksum(msg, MODES_LONG_MSG_BITS);
    int best[COHERENT_REPAIR_FLIPS], nbest = 0;
    uint32_t best_cost = UINT32_MAX;
    int idx[COHERENT_REPAIR_FLIPS];
    for (int k = 1; k <= COHERENT_REPAIR_FLIPS && k <= n; k++) {
        for (int j = 0; j < k; j++)
            idx[j] = j;
        for (;;) {
            uint32_t s = s0, c = 0;
            unsigned d = df;
            for (int j = 0; j < k; j++) {
                s ^= syn[idx[j]];
                c += cost[idx[j]];
                if (bit[idx[j]] < 5)
                    d ^= 0x10u >> bit[idx[j]];
            }
            if (s == 0 && c < best_cost && (d == 17 || d == 18)) {
                best_cost = c;
                nbest = k;
                for (int j = 0; j < k; j++)
                    best[j] = bit[idx[j]];
            }
            // Next combination
            int j = k - 1;
            while (j >= 0 && idx[j] == n - k + j)
                j--;
            if (j < 0)
                break;
            idx[j]++;
            for (int m = j + 1; m < k; m++)
                idx[m] = idx[m - 1] + 1;
        }
    }

    for (int j = 0; j < nbest; j++)
        msg[best[j] >> 3] ^= (uint8_t) (0x80 >> (best[j] & 7));
    return nbest;
}

#endif /* DEMOD_2400X && DEMOD_2400X_COHERENT */
