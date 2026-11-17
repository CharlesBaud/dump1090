// Part of dump1090, a Mode S message decoder for RTLSDR devices.
//
// weakest_bits.benchmark.c: starch benchmark for weakest_bits. Every
// implementation is checked to be bit-exact against the plain C version
// (generic_generic) before it is timed.
//
// This file is free software: you may copy, redistribute and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 2 of the License, or (at your
// option) any later version.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define WB_BENCH_BITS   112
#define WB_BENCH_MAX    32
#define WB_BENCH_POISON 0xA5

// Costs like a frame's (an average bit costs about 1170), with near-ties that
// fall into one key, equal costs, costs over 511 * 4, and bits out of reach
static void weakest_bits_make(uint32_t *c, unsigned seed)
{
    srand(seed);
    for (int i = 0; i < WB_BENCH_BITS; ++i)
        c[i] = (uint32_t) (rand() % 2400);
    c[17] = c[5];
    c[90] = c[5] + 1;
    c[33] = c[34] = 3;
    c[60] = 2100;
    c[61] = 65535;
    c[2] = c[9] = UINT32_MAX;
}

// The key the kernel orders by (see weakest_bits.generic.c)
static unsigned weakest_bits_key(const uint32_t *c, unsigned bit)
{
    const uint32_t q = c[bit] * 3 >> 3;
    return (q < 511 ? q : 511) << 7 | bit;
}

void STARCH_BENCHMARK(weakest_bits) (void)
{
    uint32_t c[WB_BENCH_BITS], key[WB_BENCH_MAX];
    unsigned count;

    memset(key, WB_BENCH_POISON, sizeof(key));
    memset(&count, WB_BENCH_POISON, sizeof(count));
    for (unsigned seed = 1; seed <= 3; ++seed) {
        weakest_bits_make(c, seed);
        STARCH_BENCHMARK_RUN( weakest_bits, c, 12, 1297, key, &count );
        STARCH_BENCHMARK_RUN( weakest_bits, c, 24, 1472, key, &count );
        STARCH_BENCHMARK_RUN( weakest_bits, c, 12, 65535, key, &count );
        STARCH_BENCHMARK_RUN( weakest_bits, c, 32, 200, key, &count );   // fewer within the ceiling
    }
}

bool STARCH_BENCHMARK_VERIFY(weakest_bits) (const uint32_t *c, unsigned n, uint32_t ceiling, uint32_t *key, unsigned *count)
{
    uint32_t ref[WB_BENCH_MAX];
    unsigned ref_count;
    bool okay = true;

    memset(ref, WB_BENCH_POISON, sizeof(ref));
    starch_weakest_bits_generic_generic(c, n, ceiling, ref, &ref_count);

    // The reference itself: in order of the key, each within the ceiling, none missed
    for (unsigned j = 0; j < ref_count && okay; ++j) {
        const unsigned bit = ref[j] & 255;
        if (bit >= WB_BENCH_BITS || c[bit] > ceiling || (ref[j] >> 8) != c[bit] ||
            (j && weakest_bits_key(c, bit) <= weakest_bits_key(c, ref[j - 1] & 255))) {
            fprintf(stderr, "reference check failed: entry %u\n", j);
            okay = false;
        }
    }
    unsigned within = 0;
    for (int i = 0; i < WB_BENCH_BITS; ++i)
        within += c[i] <= ceiling;
    if (ref_count != (within < n ? within : n)) {
        fprintf(stderr, "reference check failed: count %u, %u within the ceiling\n", ref_count, within);
        okay = false;
    }

    if (*count != ref_count) {
        fprintf(stderr, "verification failed: count %u, expected %u\n", *count, ref_count);
        okay = false;
    }
    // Entries past the count are don't-care (and may hold an earlier timed run's)
    for (unsigned j = 0; j < ref_count && okay; ++j) {
        if (key[j] != ref[j]) {
            fprintf(stderr, "verification failed: entry %u: %08x, expected %08x\n", j, key[j], ref[j]);
            okay = false;
        }
    }

    memset(key, WB_BENCH_POISON, WB_BENCH_MAX * sizeof(*key));
    memset(count, WB_BENCH_POISON, sizeof(*count));
    return okay;
}
