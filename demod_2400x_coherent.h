// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// demod_2400x_coherent.h: coherent second pass for DEMOD_2400X
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#ifndef DEMOD_2400X_COHERENT_H
#define DEMOD_2400X_COHERENT_H

#include <stdint.h>

// Settings (each can be overridden at build time, e.g. make CPPFLAGS=-DCOHERENT_STATS=1)
#ifndef COHERENT_GATE
#define COHERENT_GATE   96u     // minimum magnitude quality of the first 56 bits to try
#endif
#ifndef COHERENT_RHO2
#define COHERENT_RHO2   614     // minimum carrier coherence, in 1/1024 (0.6)
#endif
#ifndef COHERENT_MIX
#define COHERENT_MIX    0       // 1: if the coherent frame fails, also try it with each half from the magnitude bits
#endif
#ifndef COHERENT_STATS
#define COHERENT_STATS  0       // 1: print statistics at exit
#endif
#ifndef COHERENT_DIAG
#define COHERENT_DIAG   0       // 1: also decode accepted frames, to count bit errors
#endif
#ifndef COHERENT_DUMP
#define COHERENT_DUMP   0       // 1: with COHERENT_DIAG, one line per frame (frequency, errors)
#endif
#ifndef COHERENT_REPAIR_BITS
#define COHERENT_REPAIR_BITS  12      // soft-decision repair: the least reliable bits it may flip
#endif
#ifndef COHERENT_REPAIR_FLIPS
#define COHERENT_REPAIR_FLIPS 4       // and at most this many of them at once
#endif
#ifndef COHERENT_REPAIR_MAX_COST
#define COHERENT_REPAIR_MAX_COST   3u  // confidence: at most this many average bits' worth changed
#endif
#ifndef COHERENT_REPAIR_MIN_MARGIN
#define COHERENT_REPAIR_MIN_MARGIN 1u  // ... and the next-best valid repair this much costlier
#endif
#ifndef COHERENT_APRIORI_WEAK
#define COHERENT_APRIORI_WEAK      24  // tracked-aircraft repair: DF and address bits it changes must be among the frame's this many weakest
#endif
#ifndef COHERENT_APRIORI_CANDIDATES
#define COHERENT_APRIORI_CANDIDATES 8  // ... trying at most this many aircraft
#endif
#ifndef COHERENT_APRIORI_MAG_FLIPS
#define COHERENT_APRIORI_MAG_FLIPS 3  // ... flipping at most this many other bits when repairing the magnitude bits
#endif
#ifndef COHERENT_REPAIR_MAG_GATE
#define COHERENT_REPAIR_MAG_GATE 128u // minimum magnitude quality to repair from the magnitude bits
#endif
#ifndef COHERENT_SAMPLES
#define COHERENT_SAMPLES      NULL    // file name for the chip-sample dump, e.g. "chips.bin"
#endif
#ifndef COHERENT_SAMPLES_GATE
#define COHERENT_SAMPLES_GATE 257u    // also dump failed candidates of at least this quality (257: none)
#endif

// Coherent decode of the frame whose preamble starts at 12 MHz position 'pos'
// (the position the magnitude decoder chose). iq: 12-bit I/Q pairs, iq_len
// samples. out: 14 bytes; rel: 112 bit reliabilities |y_a - y_b| (0 where not
// sliced); q1, q2: quality of the two 56-bit halves (0..256); rho2, rho2b: carrier
// coherence of each half in 1/1024 (~1024 for a stable carrier); below
// COHERENT_RHO2 (first half) the bits aren't sliced.
// Returns 0, or -1 if the frame runs past iq_len.
int coherentDecode(const int16_t *iq, unsigned iq_len, int pos, uint8_t *out, uint16_t *rel,
                   unsigned *q1, unsigned *q2, unsigned *rho2, unsigned *rho2b);

// One frame's bits as the repairs see them, worked out once and shared by
// softRepair() and aprioriRepair(): set rel and ready = 0; the rest is filled in
// when a repair first needs it
struct repairBits {
    const uint16_t *rel;                    // each bit's reliability
    int ready;                              // 0: not yet, 1: costs worked out, -1: a half wasn't sliced
    uint32_t c[MODES_LONG_MSG_BITS];        // each bit's cost: reliability over its half's average
    int n;                                  // the n cheapest bits so far, in order (up to COHERENT_APRIORI_WEAK),
    uint32_t key[COHERENT_APRIORI_WEAK];    // each as its cost << 8 | the bit
};

// Soft-decision repair of bits that fail the CRC, from either decoder: rel holds
// each bit's reliability (coherent |y_a - y_b|, or the magnitude chips' log power
// contrast). Flips up to COHERENT_REPAIR_FLIPS of the COHERENT_REPAIR_BITS least
// reliable bits and keeps the cheapest combination (sum of the flipped
// reliabilities) that makes a DF17/18 with a correct CRC, if confident: its cost at
// most COHERENT_REPAIR_MAX_COST average bits, and the next-best valid combination at
// least COHERENT_REPAIR_MIN_MARGIN costlier. Returns the number of bits flipped, with
// msg repaired, or 0.
int softRepair(uint8_t *msg, struct repairBits *rb);

// Soft-decision repair assuming the frame is a DF17 from a tracked aircraft: for the
// reliable aircraft whose DF17 DF and address bits differ from the frame's only in
// bits among its COHERENT_APRIORI_WEAK weakest, those bits are set to the aircraft's,
// then up to max_flips of the COHERENT_REPAIR_BITS weakest other bits are flipped.
// The cheapest valid result over all of them, if confident. Returns the number of
// bits changed, with msg repaired, or 0.
int aprioriRepair(uint8_t *msg, struct repairBits *rb, int max_flips);

// Complex filter outputs (int32 I, Q, before rounding) at the chip centres
// pos + 6k, k = 0 .. nchips-1: 16 preamble chips then 2 per data bit.
// Returns 0, or -1 if they run past iq_len.
#define COHERENT_CHIPS (16 + 2 * 112)
int coherentChips(const int16_t *iq, unsigned iq_len, int pos, int32_t (*chips)[2], int nchips);

#endif
