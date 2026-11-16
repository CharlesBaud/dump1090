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

// Repair costs are in units of REPAIR_UNIT per average bit (see repairCosts())
#define REPAIR_UNIT (65536u / 56)

// Whether the DF bits in 'flip' (bit 4 = the frame's first bit) could all be among
// the frame's 'weak' least reliable bits: each must at least be among those of the
// first half, where the raw reliabilities compare. A quick test before any scaling
static int dfFlippable(const uint16_t *rel, unsigned flip, int weak)
{
    for (int b = 0; b < 5; b++) {
        if (!(flip & (0x10u >> b)))
            continue;
        int weaker = 0;
        for (int i = 0; i < 56; i++)
            if (rel[i] < rel[b] && ++weaker >= weak)
                return 0;
    }
    return 1;
}

// Each half has its own channel scaling, so a bit's cost (how unlikely it is to be
// wrong) is its reliability over its half's sum: rel * (2^32 / sum) >> 16, one
// division per half, REPAIR_UNIT for an average bit. Returns 0 if a half wasn't sliced
static int repairCosts(const uint16_t *rel, uint32_t *c)
{
    uint32_t sum[2] = { 0, 0 };
    for (int i = 0; i < MODES_LONG_MSG_BITS; i++)
        sum[i / 56] += rel[i];
    if (!sum[0] || !sum[1])
        return 0;
    const uint64_t inv[2] = { 0xffffffffu / sum[0], 0xffffffffu / sum[1] };
    for (int i = 0; i < MODES_LONG_MSG_BITS; i++)
        c[i] = (uint32_t) ((rel[i] * inv[i / 56]) >> 16);
    return 1;
}

// The n cheapest bits, in order, skipping those set in 'fixed' (bit 31 = the frame's
// first bit, covering bits 0..31). One pass, keeping a sorted list. Returns how many
static int weakestBits(const uint32_t *c, uint32_t fixed, int n, int *bit, uint32_t *cost)
{
    int k = 0;
    for (int i = 0; i < MODES_LONG_MSG_BITS; i++) {
        if (i < 32 && (fixed & (0x80000000u >> i)))
            continue;
        if (k == n && c[i] >= cost[k - 1])
            continue;
        int j = k < n ? k++ : k - 1;
        for (; j > 0 && cost[j - 1] > c[i]; j--) {
            cost[j] = cost[j - 1];
            bit[j] = bit[j - 1];
        }
        cost[j] = c[i];
        bit[j] = i;
    }
    return k;
}

// The cheapest set of up to COHERENT_REPAIR_FLIPS of the n bits whose syndromes XOR
// to s0 (that is, that makes the CRC correct), as a mask over bit[]; with check_df,
// only sets that leave a DF17/18 (df: the frame's DF). *best gets its cost and *second
// the next cheapest valid set's (UINT32_MAX if none). Returns 0 if there is none.
//
// A set of up to 4 is two disjoint sets of up to 2 whose syndromes XOR to s0. The
// sets (the empty one too) go in a small hash table by syndrome, and each is looked
// up against s0 XOR its own: about 80 lookups instead of trying every combination of
// up to 4 of 12 bits (793)
static unsigned flipSearch(uint32_t s0, const int *bit, const uint32_t *cost, int n,
                           int check_df, unsigned df, uint32_t *best, uint32_t *second)
{
    _Static_assert(COHERENT_REPAIR_FLIPS <= 4 && COHERENT_REPAIR_BITS <= 16, "flipSearch() pairs sets of at most 2 bits");
    uint32_t syn[COHERENT_REPAIR_BITS];
    for (int j = 0; j < n; j++)
        syn[j] = bitSyndrome(bit[j]);

    struct { uint32_t syn, cost; unsigned mask; } set[1 + COHERENT_REPAIR_BITS * (COHERENT_REPAIR_BITS + 1) / 2];
    int ns = 0;
    set[ns].syn = 0; set[ns].cost = 0; set[ns++].mask = 0;
    for (int i = 0; i < n; i++) {
        set[ns].syn = syn[i]; set[ns].cost = cost[i]; set[ns++].mask = 1u << i;
        for (int j = i + 1; j < n; j++) {
            set[ns].syn = syn[i] ^ syn[j]; set[ns].cost = cost[i] + cost[j]; set[ns++].mask = (1u << i) | (1u << j);
        }
    }
    uint8_t slot[256];                      // set index + 1 by the syndrome's low byte, linear probing
    memset(slot, 0, sizeof(slot));
    for (int k = 0; k < ns; k++) {
        unsigned h = set[k].syn & 255;
        while (slot[h])
            h = (h + 1) & 255;
        slot[h] = (uint8_t) (k + 1);
    }

    unsigned best_mask = 0;
    *best = *second = UINT32_MAX;
    for (int k = 0; k < ns; k++) {
        const uint32_t target = s0 ^ set[k].syn;
        for (unsigned h = target & 255; slot[h]; h = (h + 1) & 255) {
            const int m = slot[h] - 1;
            if (set[m].syn != target || (set[m].mask & set[k].mask))
                continue;
            const unsigned mask = set[m].mask | set[k].mask;
            const uint32_t c = set[m].cost + set[k].cost;
            if (!mask || mask == best_mask || c >= *second || __builtin_popcount(mask) > COHERENT_REPAIR_FLIPS)
                continue;
            if (check_df) {
                unsigned d = df;            // DF17/18 after the flips
                for (int j = 0; j < n; j++)
                    if ((mask & (1u << j)) && bit[j] < 5)
                        d ^= 0x10u >> bit[j];
                if (d != 17 && d != 18)
                    continue;
            }
            if (c < *best) {
                *second = *best;
                *best = c;
                best_mask = mask;
            } else {
                *second = c;
            }
        }
    }
    return best_mask;
}

// Confident enough: cheap, and clearly cheaper than any other valid repair
static int repairConfident(uint32_t best, uint32_t second)
{
    return best <= COHERENT_REPAIR_MAX_COST * REPAIR_UNIT &&
        (second == UINT32_MAX || second - best >= COHERENT_REPAIR_MIN_MARGIN * REPAIR_UNIT);
}

static int applyFlips(uint8_t *msg, unsigned mask, const int *bit, int n)
{
    int flips = 0;
    for (int j = 0; j < n; j++) {
        if (mask & (1u << j)) {
            msg[bit[j] >> 3] ^= (uint8_t) (0x80 >> (bit[j] & 7));
            flips++;
        }
    }
    return flips;
}

int softRepair(uint8_t *msg, const uint16_t *rel)
{
    // Most frames can't become a DF17/18 by flipping only unreliable bits
    const unsigned df = msg[0] >> 3;
    if ((__builtin_popcount(df ^ 17) > COHERENT_REPAIR_FLIPS || !dfFlippable(rel, df ^ 17, COHERENT_REPAIR_BITS)) &&
        (__builtin_popcount(df ^ 18) > COHERENT_REPAIR_FLIPS || !dfFlippable(rel, df ^ 18, COHERENT_REPAIR_BITS)))
        return 0;

    uint32_t c[MODES_LONG_MSG_BITS];
    if (!repairCosts(rel, c))
        return 0;
    int bit[COHERENT_REPAIR_BITS];
    uint32_t cost[COHERENT_REPAIR_BITS];
    const int n = weakestBits(c, 0, COHERENT_REPAIR_BITS, bit, cost);

    // The same test exactly: the DF bits to change must all be among them
    unsigned weak_df = 0;                   // DF bits (bit 4 = first) among them
    for (int j = 0; j < n; j++)
        if (bit[j] < 5)
            weak_df |= 0x10u >> bit[j];
    if (((df ^ 17) & ~weak_df) && ((df ^ 18) & ~weak_df))
        return 0;

    uint32_t best, second;
    const unsigned mask = flipSearch(modesChecksum(msg, MODES_LONG_MSG_BITS), bit, cost, n, 1, df, &best, &second);
    if (!mask || !repairConfident(best, second))
        return 0;
    return applyFlips(msg, mask, bit, n);
}

// DF and address bits of a DF17 (bits 0..4 and 8..31; the capability bits 5..7 vary)
#define APRIORI_FIXED 0xf8ffffffu

int aprioriRepair(uint8_t *msg, const uint16_t *rel)
{
    // The DF bits that would have to become 17 must be weak (a quick test first)
    if (!dfFlippable(rel, (unsigned) (msg[0] >> 3) ^ 17, COHERENT_APRIORI_WEAK))
        return 0;

    // A bit may only be changed if it is among the frame's COHERENT_APRIORI_WEAK
    // weakest: a frame from an untracked aircraft would otherwise be taken for a
    // tracked one a few address bits away, changing bits that were received right
    uint32_t c[MODES_LONG_MSG_BITS];
    if (!repairCosts(rel, c))
        return 0;
    int wbit[COHERENT_APRIORI_WEAK];
    uint32_t wcost[COHERENT_APRIORI_WEAK];
    const int nw = weakestBits(c, 0, COHERENT_APRIORI_WEAK, wbit, wcost);
    uint32_t weak = 0;                      // among the bits 0..31
    for (int j = 0; j < nw; j++)
        if (wbit[j] < 32)
            weak |= 0x80000000u >> wbit[j];

    // Reliable tracked aircraft, recently heard in a CRC-checked message, whose DF17
    // differs from these bits only in weak ones
    const uint32_t w = ((uint32_t) msg[0] << 24) | ((uint32_t) msg[1] << 16) | ((uint32_t) msg[2] << 8) | msg[3];
    uint32_t hyp[COHERENT_APRIORI_CANDIDATES];
    int nh = 0;
    for (const struct aircraft *a = Modes.aircrafts; a && nh < COHERENT_APRIORI_CANDIDATES; a = a->next) {
        const uint32_t h = (17u << 27) | a->addr;
        if (!((w ^ h) & APRIORI_FIXED & ~weak) && !(a->addr & MODES_NON_ICAO_ADDRESS) && a->reliable &&
            icaoFilterTest(a->addr))
            hyp[nh++] = h;
    }
    if (!nh)
        return 0;

    // The rest: the weakest bits outside the fixed ones
    int bit[COHERENT_REPAIR_BITS];
    uint32_t cost[COHERENT_REPAIR_BITS];
    const int n = weakestBits(c, APRIORI_FIXED, COHERENT_REPAIR_BITS, bit, cost);

    uint32_t best = UINT32_MAX, second = UINT32_MAX, best_h = 0;
    unsigned best_mask = 0;
    for (int k = 0; k < nh; k++) {
        const uint32_t d = (w ^ hyp[k]) & APRIORI_FIXED;
        uint32_t cfix = 0;
        for (int i = 0; i < 32; i++)
            if (d & (0x80000000u >> i))
                cfix += c[i];
        uint8_t m[MODES_LONG_MSG_BYTES];
        memcpy(m, msg, sizeof(m));
        const uint32_t v = w ^ d;
        m[0] = (uint8_t) (v >> 24); m[1] = (uint8_t) (v >> 16); m[2] = (uint8_t) (v >> 8); m[3] = (uint8_t) v;
        uint32_t b1, b2;
        const unsigned mask = flipSearch(modesChecksum(m, MODES_LONG_MSG_BITS), bit, cost, n, 0, 0, &b1, &b2);
        if (b1 == UINT32_MAX)
            continue;
        b1 += cfix;
        b2 = b2 == UINT32_MAX ? b2 : b2 + cfix;
        if (b1 < best) {
            second = best < b2 ? best : b2;
            best = b1;
            best_h = hyp[k];
            best_mask = mask;
        } else if (b1 < second) {
            second = b1;
        }
    }
    if (best == UINT32_MAX || !repairConfident(best, second))
        return 0;

    const uint32_t d = (w ^ best_h) & APRIORI_FIXED;
    const uint32_t v = w ^ d;
    msg[0] = (uint8_t) (v >> 24); msg[1] = (uint8_t) (v >> 16); msg[2] = (uint8_t) (v >> 8); msg[3] = (uint8_t) v;
    return __builtin_popcount(d) + applyFlips(msg, best_mask, bit, n);
}

#endif /* DEMOD_2400X && DEMOD_2400X_COHERENT */
