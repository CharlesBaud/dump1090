# DEMOD_2400X demodulator

This branch (`demod-2400x-dev`, based on upstream `dev`) adds a new Mode S
demodulator for 2.4 MHz SDR input (`demod_2400x.c`, enabled with
`-DDEMOD_2400X -DUSE_MAGNITUDE` in the Makefile; the build reports version
`11.2x` unless `DUMP1090_VERSION` is set). A coherent second pass
(`-DDEMOD_2400X_COHERENT`) rescues frames that fail the CRC; see
[COHERENT-demod2400x.md](COHERENT-demod2400x.md). Both are on by default and can
be switched off when running:

* `--no-demod2400x`: the stock demodulator, exactly as upstream (stock
  converters, magnitudes in the FIFO; output byte-identical to upstream `dev`);
* `--no-coherent`: DEMOD_2400X without the coherent pass.

In outline:

* DC removal and conversion to 12-bit I/Q in the input converter, so the
  FIFO keeps I/Q samples for later processing.
* A 5-phase polyphase FIR upsamples to 12 MHz (least-squares fractional-delay
  coefficients, all integer arithmetic), then stores an 8-bit log of the power
  (8 × log2, from the bits of a float conversion: no square root).
* A preamble metric at every 12 MHz position, and the 112 bits are decoded
  from the best of the neighbouring positions.
* The arrival time is refined by correlating the decoded chips with the
  signal, so MLAT timestamps are unbiased across sub-sample phases.
* Adaptive gain and signal level statistics behave as in the stock
  demodulator.

The heavy kernels (`iqdata_conv_*`, `upsample12`, `locate_preamble`,
`decode_bits`) are starch kernels with generic, ARM NEON (ARMv7 and AArch64)
and x86 AVX2 versions. All versions give bit-identical results, which
`starch-benchmark -V` checks, including on full-range worst-case input.

## Installing over a dump1090-fa package

On a system with the dump1090-fa package (PiAware, or the FlightAware .deb),
build with `make`, then `sudo sh install-2400x-dev.sh`. It copies the build to
`/usr/local/bin/dump1090-2400x` and adds a systemd override so the dump1090-fa
service runs it with the same settings (`/etc/default/dump1090-fa`, PiAware
location, user, JSON directory); the packaged files are not changed. Run it
again after a new build; `sudo sh uninstall-2400x-dev.sh` goes back to the
packaged decoder installed by apt-get.

## Results

A 60 s capture from a Raspberry Pi 4 with a FlightAware Pro Stick Plus
(gain 49.6 dB), with an outdoor antenna on a 3 m pole and 12 m of LMR400,
decoded by stock dump1090-fa 11.1 and by this branch with the same options
(`--fix`, receiver position set, `--throttle`). "+ coherent" is DEMOD_2400X with the
coherent second pass; changes are against stock. Development also used an SDRplay
recording and synthetic frames; see
[NOTES-demod2400x.md](NOTES-demod2400x.md#test-data).

| | Stock | DEMOD_2400X | Change | + coherent | Change |
|---|---|---|---|---|---|
| Usable messages | 73,146 | 84,545 | +16% | 93,827 | +28% |
| DF17 (ADS-B) | 10,352 | 11,830 | +14% | 13,365 | +29% |
| Global CPR positions | 3,510 | 3,977 | +13% | 4,538 | +29% |
| Aircraft | 136 | 139 | +3 | 150 | +14 |

Positions and tracks are timed by the wall clock, so files are decoded with
`--throttle` (at the speed they were captured). Without it the minute is
decoded in a few seconds, aircraft appear to move many times faster, some
positions fail the speed check, and the counts vary from run to run and
with CPU speed. With `--throttle` every figure repeats exactly.

The gain is concentrated in weak signals. Messages by signal strength (from the
signal byte of the Beast output; the first fraction of a second, before the
connection, is missed by all three alike):

| RSSI (dBFS) | Stock | DEMOD_2400X | Change | + coherent | Change |
|---|---|---|---|---|---|
| below −27 | 144 | 304 | +111% | 1,035 | +619% |
| −27 to −24 | 2,732 | 4,724 | +73% | 6,749 | +147% |
| −24 to −21 | 6,454 | 8,588 | +33% | 9,927 | +54% |
| −21 to −18 | 9,587 | 11,423 | +19% | 12,642 | +32% |
| −18 to −15 | 12,013 | 13,558 | +13% | 14,731 | +23% |
| −15 to −10 | 17,956 | 19,675 | +10% | 21,184 | +18% |
| −10 to −5 | 9,902 | 10,639 | +7% | 11,337 | +14% |
| above −5 | 14,018 | 15,272 | +9% | 15,802 | +13% |

Position messages by distance from the receiver:

| Distance | Stock | DEMOD_2400X | Change | + coherent | Change | Aircraft (stock / DEMOD_2400X / + coherent) |
|---|---|---|---|---|---|---|
| 0–50 km | 737 | 818 | +11% | 883 | +20% | 12 / 13 / 13 |
| 50–100 km | 1,444 | 1,583 | +10% | 1,750 | +21% | 30 / 31 / 33 |
| 100–150 km | 893 | 1,014 | +14% | 1,173 | +31% | 33 / 33 / 33 |
| 150–200 km | 438 | 527 | +20% | 626 | +43% | 20 / 23 / 24 |
| 200–250 km | 121 | 145 | +20% | 196 | +62% | 10 / 13 / 15 |
| 250–300 km | 16 | 25 | +56% | 48 | +200% | 2 / 4 / 4 |
| over 300 km | 37 | 49 | +32% | 62 | +68% | 5 / 5 / 5 |

Weak synthetic frames (2,000 DF17 frames with a random ±100 kHz carrier
offset, decoded out of 2,000): 689 (stock), 1,407 (DEMOD_2400X), 1,592
(+ coherent).
The generator and checker are in [sim/](sim/).

## Timing (MLAT)

On synthetic frames at known 12 MHz ticks, covering all 5 sub-sample
phases:

| | Stock | DEMOD_2400X |
|---|---|---|
| Bias across phases (clean) | 0 to +2 ticks | none |
| Jitter, noisy frames (σ) | 0.48 ticks | 0 (every frame on the right tick) |

One tick is 83.3 ns. Running live on the Pi 4, mlat-client stays synchronized
with about 260–270 nearby receivers, as with stock.

## CPU

Demodulation time for the 60 s capture (`--fix`):

| | Stock | DEMOD_2400X | + coherent |
|---|---|---|---|
| Raspberry Pi 5 (Cortex-A76) | 2.1 s | 3.0 s (gcc 14), 2.8 s (clang 19) | 3.7 s (gcc), 3.6 s (clang) |
| Raspberry Pi 4 (Cortex-A72) | 4.0 s | 5.8 s | 7.5 s |
| Raspberry Pi 3B (Cortex-A53, 32-bit) | 8.6 s | 18.4 s | 23.5 s |

Measured with `--no-demod2400x`, `--no-coherent` and the default. The Pi 4
runs were on three cores while the live receiver was also running on the same
Pi, so they may be slightly high. The Pi 3B times are for the first 20 s of
the capture, times 3; on a Pi 3B the default uses about 39% of one core for
demodulation (see
[NOTES-demod2400x.md](NOTES-demod2400x.md#raspberry-pi-3b)). A Pi 3B needs a
heatsink: if it overheats and throttles to 600 MHz it can't keep up and drops
samples (`vcgencmd get_throttled` shows 0x0 if it has never throttled), and
`--no-coherent` or `--no-demod2400x` reduce the load. Live on the Pi 4,
with about 105,000 messages a minute, demodulation takes about 25% of one
core. These results come from one receiver and a one-minute capture; longer
side-by-side runs are in progress.

## Compilers

Built with the same options as the GitHub workflow
(`make RTLSDR=yes BLADERF=yes HACKRF=yes LIMESDR=yes SOAPYSDR=yes all`,
then `make test`):

| Compiler | Target | Warnings | `make test` |
|---|---|---|---|
| gcc 14.2 | AArch64 (Pi 5) | 0 | pass (74 validations) |
| clang 19.1 | AArch64 (Pi 5) | 0 | pass (74 validations) |
| gcc 14.2 (cross) | ARMv7 | 0 | `starch-benchmark -V` passes (74), run natively on AArch64 |
| gcc 14.2 (cross) | x86-64 AVX2 | 0 | pass under qemu-x86_64 (68 validations) |
| clang 19.1 (cross) | x86-64 AVX2 | 0 | pass under qemu-x86_64 (68 validations) |
| gcc (native) | x86-64 AVX2 | 0 | `starch-benchmark -V` passes on real hardware |

The gcc and clang AArch64 builds give byte-identical decoder output,
timestamps included. On x86, clang with `-ffast-math` compiles every AVX2
divide and square root to a reciprocal estimate (`vrcpps`/`vrsqrtps`), so the
kernels correct the estimate in both directions to stay exact. QEMU computes
these estimates exactly, so x86 correctness also needs a run on real
hardware. All the AVX2 kernels, including the mirrored upsampler,
`decode_bits` and `coherent_decode`, pass `starch-benchmark -V` on real x86-64
hardware, also after the change to 8-bit log power, which made the AVX2
`upsample12` 39% faster (1.43 → 0.87 ms per call) and `locate_preamble` 47%
faster (0.38 → 0.20 ms).
