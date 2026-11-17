// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// weakest_bits.neon.c: NEON version of weakest_bits.generic.c (the same result,
// bit for bit). All 112 keys sit in 14 vectors of 8; each round takes the
// smallest key and marks it taken, with no branches.
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#ifndef __ARM_NEON
#  error Neon support not available
#endif

#include <arm_neon.h>

#ifndef WB_BITS
#define WB_BITS   112           // bits in a frame
#define WB_MAX_N  32            // at most this many wanted
#define WB_OUT    0xffffu       // key of a bit left out
#endif
#define WB_VECS   (WB_BITS / 8)

// Smallest lane (ARMv7 has no across-vector min: pairwise steps instead)
static inline uint16_t wb_min_lane(uint16x8_t v)
{
#ifdef __aarch64__
    return vminvq_u16(v);
#else
    uint16x4_t m = vmin_u16(vget_low_u16(v), vget_high_u16(v));
    m = vpmin_u16(m, m);
    m = vpmin_u16(m, m);
    return vget_lane_u16(m, 0);
#endif
}

void STARCH_IMPL(weakest_bits, neon) (const uint32_t * restrict c, unsigned n, uint32_t ceiling,
                                      uint32_t * restrict key, unsigned * restrict count)
{
    static const uint16_t lane[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    const uint32x4_t ceil = vdupq_n_u32(ceiling);
    const uint16x8_t cap = vdupq_n_u16(511);
    uint16x8_t idx = vld1q_u16(lane);
    uint16x8_t k[WB_VECS];

    // Keys: min(cost * 3 / 8, 511) << 7 | bit, all ones if above the ceiling (where
    // the multiply may overflow)
    for (int v = 0; v < WB_VECS; v++) {
        const uint32x4_t a = vld1q_u32(c + 8 * v), b = vld1q_u32(c + 8 * v + 4);
        uint16x8_t q = vcombine_u16(vqshrn_n_u32(vmulq_n_u32(a, 3), 3), vqshrn_n_u32(vmulq_n_u32(b, 3), 3));
        q = vorrq_u16(vshlq_n_u16(vminq_u16(q, cap), 7), idx);
        const uint16x8_t out = vcombine_u16(vmovn_u32(vcgtq_u32(a, ceil)), vmovn_u32(vcgtq_u32(b, ceil)));
        k[v] = vorrq_u16(q, out);
        idx = vaddq_u16(idx, vdupq_n_u16(8));
    }

    if (n > WB_MAX_N)
        n = WB_MAX_N;
    unsigned got = 0;
    for (; got < n; got++) {
        // The smallest: a tree of mins, then across the lanes
        uint16x8_t m[WB_VECS / 2];
        for (int v = 0; v < WB_VECS / 2; v++)
            m[v] = vminq_u16(k[2 * v], k[2 * v + 1]);
        const uint16x8_t m01 = vminq_u16(vminq_u16(m[0], m[1]), vminq_u16(m[2], m[3]));
        const uint16x8_t m23 = vminq_u16(vminq_u16(m[4], m[5]), m[6]);
        const uint16_t x = wb_min_lane(vminq_u16(m01, m23));
        if (x == WB_OUT)
            break;
        const unsigned bit = x & 127;
        key[got] = c[bit] << 8 | bit;

        // Taken: its lane (the only one equal to it) becomes all ones
        const uint16x8_t xv = vdupq_n_u16(x);
        for (int v = 0; v < WB_VECS; v++)
            k[v] = vorrq_u16(k[v], vceqq_u16(k[v], xv));
    }
    *count = got;
}
