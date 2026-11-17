// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// weakest_bits.generic.c: the n cheapest of a frame's 112 bit costs, for the
// soft-decision repair (demod_2400x_coherent.c). Plain C reference.
//
// Bits costing more than 'ceiling' are left out. The rest are ordered by a 16-bit
// key, min(cost * 3 / 8, 511) << 7 | bit, so equal keys can't happen and the bit
// breaks near-ties (* 3 / 8: the ceilings the repairs use, 1297 and 1472, come to
// 486 and 552, so 9 bits are nearly all used); each comes out as cost << 8 | bit, with its full cost
// (at most 65535). *count gets how many: n (at most WB_MAX_N), or fewer if fewer
// bits are within the ceiling.
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#ifndef WB_BITS
#define WB_BITS   112           // bits in a frame
#define WB_MAX_N  32            // at most this many wanted
#define WB_OUT    0xffffu       // key of a bit left out
#endif

// One pass keeping a sorted list: the quickest way without vectors
void STARCH_IMPL(weakest_bits, generic) (const uint32_t * restrict c, unsigned n, uint32_t ceiling,
                                         uint32_t * restrict key, unsigned * restrict count)
{
    uint16_t list[WB_MAX_N];
    unsigned k = 0;

    if (n > WB_MAX_N)
        n = WB_MAX_N;
    for (unsigned i = 0; i < WB_BITS && n; i++) {
        if (c[i] > ceiling)
            continue;
        const uint32_t q = c[i] * 3 >> 3;  // c[i] <= ceiling <= 65535: no overflow
        const uint16_t x = (uint16_t) ((q < 511 ? q : 511) << 7 | i);
        if (k == n && x >= list[k - 1])
            continue;
        unsigned j = k < n ? k++ : k - 1;
        for (; j > 0 && list[j - 1] > x; j--)
            list[j] = list[j - 1];
        list[j] = x;
    }
    for (unsigned j = 0; j < k; j++) {
        const unsigned bit = list[j] & 127;
        key[j] = c[bit] << 8 | bit;
    }
    *count = k;
}
