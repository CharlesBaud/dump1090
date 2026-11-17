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
        int weaker = 0;                     // counting them all is quicker than stopping early: no branches
        for (int i = 0; i < 56; i++)
            weaker += rel[i] < rel[b];
        if (weaker >= weak)
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

// None of the n cheapest bits costs more than this: each costs no more than the
// 112 - n + 1 bits from it up, and those add up to at most the 2 * 65536 all the
// costs add up to (repairCosts())
static uint32_t costCeiling(int n)
{
    return 2u * 65536 / (uint32_t) (MODES_LONG_MSG_BITS - n + 1);
}

// A bit and its cost (at most 65535) in one value, cost << 8 | bit, as the
// weakest_bits kernel gives them (cheapest first: see dsp/impl/weakest_bits.generic.c)
#define KEY_BIT(key)       ((int) ((key) & 255))
#define KEY_COST(key)      ((uint32_t) (key) >> 8)

// Every set of up to 2 of the n bits (the empty one too), by syndrome, in a small hash
// table: a valid repair of up to COHERENT_REPAIR_FLIPS bits is two disjoint sets
// whose syndromes XOR to the frame's. Built once per frame, then looked up for each
// syndrome to repair (the tracked-aircraft repair tries several)
struct flipTable {
    const uint32_t *key;                    // the bits, as cost << 8 | bit
    int n, ns;
    struct { uint32_t syn, cost; unsigned mask; } set[1 + COHERENT_REPAIR_BITS * (COHERENT_REPAIR_BITS + 1) / 2];
    uint8_t slot[256];                      // set index + 1 by the syndrome's low byte, linear probing
};

static void flipTableBuild(struct flipTable *t, const uint32_t *key, int n)
{
    _Static_assert(COHERENT_REPAIR_FLIPS <= 4 && COHERENT_REPAIR_BITS <= 16, "flipSearch() pairs sets of at most 2 bits");
    uint32_t syn[COHERENT_REPAIR_BITS], cost[COHERENT_REPAIR_BITS];
    for (int j = 0; j < n; j++) {
        syn[j] = modesBitSyndrome(KEY_BIT(key[j]));
        cost[j] = KEY_COST(key[j]);
    }

    int ns = 0;
    t->set[ns].syn = 0; t->set[ns].cost = 0; t->set[ns++].mask = 0;
    for (int i = 0; i < n; i++) {
        t->set[ns].syn = syn[i]; t->set[ns].cost = cost[i]; t->set[ns++].mask = 1u << i;
        for (int j = i + 1; j < n; j++) {
            t->set[ns].syn = syn[i] ^ syn[j]; t->set[ns].cost = cost[i] + cost[j]; t->set[ns++].mask = (1u << i) | (1u << j);
        }
    }
    memset(t->slot, 0, sizeof(t->slot));
    for (int k = 0; k < ns; k++) {
        unsigned h = t->set[k].syn & 255;
        while (t->slot[h])
            h = (h + 1) & 255;
        t->slot[h] = (uint8_t) (k + 1);
    }
    t->key = key;
    t->n = n;
    t->ns = ns;
}

// The cheapest valid set of up to max_flips bits for syndrome s0 (one that makes
// the CRC correct), as a mask over the table's bits; with check_df, only sets that leave a DF17/18 (df: the
// frame's DF). *best gets its cost and *second the next cheapest valid set's
// (UINT32_MAX if none). Returns 0 if there is none. About 80 lookups instead of
// trying every combination of up to 4 of 12 bits (793)
static unsigned flipSearch(const struct flipTable *t, uint32_t s0, int max_flips, int check_df, unsigned df,
                           uint32_t *best, uint32_t *second)
{
    unsigned best_mask = 0;
    *best = *second = UINT32_MAX;
    for (int k = 0; k < t->ns; k++) {
        const uint32_t target = s0 ^ t->set[k].syn;
        for (unsigned h = target & 255; t->slot[h]; h = (h + 1) & 255) {
            const int m = t->slot[h] - 1;
            if (t->set[m].syn != target || (t->set[m].mask & t->set[k].mask))
                continue;
            const unsigned mask = t->set[m].mask | t->set[k].mask;
            const uint32_t c = t->set[m].cost + t->set[k].cost;
            if (!mask || mask == best_mask || c >= *second || __builtin_popcount(mask) > max_flips)
                continue;
            if (check_df) {
                unsigned d = df;            // DF17/18 after the flips
                for (int j = 0; j < t->n; j++)
                    if ((mask & (1u << j)) && KEY_BIT(t->key[j]) < 5)
                        d ^= 0x10u >> KEY_BIT(t->key[j]);
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

static int applyFlips(uint8_t *msg, unsigned mask, const uint32_t *key, int n)
{
    int flips = 0;
    for (int j = 0; j < n; j++) {
        if (mask & (1u << j)) {
            const int b = KEY_BIT(key[j]);
            msg[b >> 3] ^= (uint8_t) (0x80 >> (b & 7));
            flips++;
        }
    }
    return flips;
}

// The frame's costs, the first time a repair needs them
static int repairCostsReady(struct repairBits *rb)
{
    if (!rb->ready) {
        rb->ready = repairCosts(rb->rel, rb->c) ? 1 : -1;
        rb->n = 0;
    }
    return rb->ready > 0;
}

// ... and its 'weakest' cheapest bits, in order (the list only grows: its first bits
// stay the same)
static int repairPrepare(struct repairBits *rb, int weakest)
{
    if (!repairCostsReady(rb))
        return 0;
    if (rb->n < weakest) {
        unsigned n;
        starch_weakest_bits(rb->c, (unsigned) weakest, costCeiling(weakest), rb->key, &n);
        rb->n = (int) n;
    }
    return 1;
}

int softRepair(uint8_t *msg, struct repairBits *rb)
{
    // Most frames can't become a DF17/18 by flipping only unreliable bits
    const unsigned df = msg[0] >> 3;
    if ((__builtin_popcount(df ^ 17) > COHERENT_REPAIR_FLIPS || !dfFlippable(rb->rel, df ^ 17, COHERENT_REPAIR_BITS)) &&
        (__builtin_popcount(df ^ 18) > COHERENT_REPAIR_FLIPS || !dfFlippable(rb->rel, df ^ 18, COHERENT_REPAIR_BITS)))
        return 0;

    if (!repairPrepare(rb, COHERENT_REPAIR_BITS))
        return 0;
    const int n = COHERENT_REPAIR_BITS;

    // The same test exactly: the DF bits to change must all be among them
    unsigned weak_df = 0;                   // DF bits (bit 4 = first) among them
    for (int j = 0; j < n; j++)
        if (KEY_BIT(rb->key[j]) < 5)
            weak_df |= 0x10u >> KEY_BIT(rb->key[j]);
    if (((df ^ 17) & ~weak_df) && ((df ^ 18) & ~weak_df))
        return 0;

    struct flipTable t;
    flipTableBuild(&t, rb->key, n);
    uint32_t best, second;
    const unsigned mask = flipSearch(&t, modesChecksum(msg, MODES_LONG_MSG_BITS), COHERENT_REPAIR_FLIPS, 1, df, &best, &second);
    if (!mask || !repairConfident(best, second))
        return 0;
    return applyFlips(msg, mask, rb->key, n);
}

// DF and address bits of a DF17 (bits 0..4 and 8..31; the capability bits 5..7 vary)
#define APRIORI_FIXED 0xf8ffffffu

// Reliable tracked aircraft, recently heard in a CRC-checked message, whose DF17
// differs from the frame's first 32 bits w only in bits set in 'weak'
static int aprioriCandidates(uint32_t w, uint32_t weak, uint32_t *hyp)
{
    int nh = 0;
    for (const struct aircraft *a = Modes.aircrafts; a && nh < COHERENT_APRIORI_CANDIDATES; a = a->next) {
        const uint32_t h = (17u << 27) | a->addr;
        if (!((w ^ h) & APRIORI_FIXED & ~weak) && !(a->addr & MODES_NON_ICAO_ADDRESS) && a->reliable &&
            icaoFilterTest(a->addr))
            hyp[nh++] = h;
    }
    return nh;
}

int aprioriRepair(uint8_t *msg, struct repairBits *rb, int max_flips)
{
    // The DF bits that would have to become 17 must be weak (a quick test first)
    if (!dfFlippable(rb->rel, (unsigned) (msg[0] >> 3) ^ 17, COHERENT_APRIORI_WEAK))
        return 0;

    // A bit may only be changed if it is among the frame's COHERENT_APRIORI_WEAK
    // weakest: a frame from an untracked aircraft would otherwise be taken for a
    // tracked one a few address bits away, changing bits that were received right.
    // Most frames have no such aircraft, so first a quick test without sorting: the
    // bits that cost no more than any of the weakest can (costCeiling())
    if (!repairCostsReady(rb))
        return 0;
    const uint32_t w = ((uint32_t) msg[0] << 24) | ((uint32_t) msg[1] << 16) | ((uint32_t) msg[2] << 8) | msg[3];
    const uint32_t ceiling = costCeiling(COHERENT_APRIORI_WEAK);
    uint32_t weak = 0;                      // among the bits 0..31
    for (int i = 0; i < 32; i++)
        if (rb->c[i] <= ceiling)
            weak |= 0x80000000u >> i;
    uint32_t hyp[COHERENT_APRIORI_CANDIDATES];
    if (!aprioriCandidates(w, weak, hyp))
        return 0;

    // Now exactly
    _Static_assert(COHERENT_APRIORI_WEAK >= COHERENT_REPAIR_BITS, "the weakest-bit list serves both repairs");
    repairPrepare(rb, COHERENT_APRIORI_WEAK);
    weak = 0;
    for (int j = 0; j < rb->n; j++)
        if (KEY_BIT(rb->key[j]) < 32)
            weak |= 0x80000000u >> KEY_BIT(rb->key[j]);
    const int nh = aprioriCandidates(w, weak, hyp);
    if (!nh)
        return 0;

    // The rest: the weakest bits outside the fixed ones. The weakest-bit list has them
    // in order, unless too many of its bits are fixed ones
    uint32_t key[COHERENT_REPAIR_BITS];
    int n = 0;
    for (int j = 0; j < rb->n && n < COHERENT_REPAIR_BITS; j++) {
        const int b = KEY_BIT(rb->key[j]);
        if (b >= 32 || !(APRIORI_FIXED & (0x80000000u >> b)))
            key[n++] = rb->key[j];
    }
    if (n < COHERENT_REPAIR_BITS) {
        uint32_t c[MODES_LONG_MSG_BITS];    // the fixed bits out of reach
        memcpy(c, rb->c, sizeof(c));
        for (int i = 0; i < 32; i++)
            if (APRIORI_FIXED & (0x80000000u >> i))
                c[i] = UINT32_MAX;
        unsigned nk;
        starch_weakest_bits(c, COHERENT_REPAIR_BITS, 65535, key, &nk);
        n = (int) nk;
    }
    struct flipTable t;
    flipTableBuild(&t, key, n);

    // Each aircraft: the syndrome with its DF and address bits set is the frame's
    // XOR the changed bits' syndromes (the CRC is linear)
    const uint32_t s = modesChecksum(msg, MODES_LONG_MSG_BITS);
    uint32_t best = UINT32_MAX, second = UINT32_MAX, best_h = 0;
    unsigned best_mask = 0;
    for (int k = 0; k < nh; k++) {
        uint32_t cfix = 0, s0 = s;
        for (uint32_t d = (w ^ hyp[k]) & APRIORI_FIXED; d; ) {
            const int i = __builtin_clz(d);  // the first bit left to change
            d ^= 0x80000000u >> i;
            cfix += rb->c[i];
            s0 ^= modesBitSyndrome(i);
        }
        uint32_t b1, b2;
        const unsigned mask = flipSearch(&t, s0, max_flips, 0, 0, &b1, &b2);
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
    return __builtin_popcount(d) + applyFlips(msg, best_mask, key, n);
}

#endif /* DEMOD_2400X && DEMOD_2400X_COHERENT */
