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

// Settings
#define COHERENT_GATE   96u     // minimum magnitude quality of the first 56 bits to try
#define COHERENT_RHO2   614     // minimum carrier coherence, in 1/1024 (0.6)
#define COHERENT_DIAG   0       // 1: also decode accepted frames, to count bit errors

// Coherent decode of the frame whose preamble starts at 12 MHz position 'pos'
// (the position the magnitude decoder chose). iq: 12-bit I/Q pairs, iq_len
// samples. out: 14 bytes; q1, q2: quality of the two 56-bit halves (0..256);
// rho2: carrier coherence over the first 56 bits in 1/1024 (~1024 for a stable
// carrier); below COHERENT_RHO2 the bits aren't sliced.
// Returns 0, or -1 if the frame runs past iq_len.
int coherentDecode(const int16_t *iq, unsigned iq_len, int pos,
                   uint8_t *out, unsigned *q1, unsigned *q2, unsigned *rho2);

// Complex filter outputs (int32 I, Q, before rounding) at the chip centres
// pos + 6k, k = 0 .. nchips-1: 16 preamble chips then 2 per data bit.
// Returns 0, or -1 if they run past iq_len.
#define COHERENT_CHIPS (16 + 2 * 112)
int coherentChips(const int16_t *iq, unsigned iq_len, int pos, int32_t (*chips)[2], int nchips);

#endif
