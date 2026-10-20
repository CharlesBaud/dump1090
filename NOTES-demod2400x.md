# DEMOD_2400X development notes

What was built, why, and what was tried and dropped. Results are in
[README-demod2400x.md](README-demod2400x.md). Times are for the 60 s Pi 4
RTL-SDR capture decoded on a Raspberry Pi 5 unless stated otherwise.

## Test data

* **SDRplay RSP1B recording**: 2.4 MHz SC16 I/Q, 17.8 s of signal. Used for
  the early development and speed work. The stock demodulator crashes on
  RSP1B input; that is a separate issue, not fixed on this branch.
* **RTL-SDR recording**: 60 s at 2.4 MHz, UC8, from the live receiver: a
  Raspberry Pi 4 with the FlightAware Pro Stick Plus (blue, internal LNA and
  SAW filter) at gain 49.6 dB, outdoor antenna on a 3 m pole, 12 m of
  LMR400. Used for the stock comparison, the later tuning and CPU figures.
* **Synthetic frames** from `sim/gen_timing.c` (based on `gen_data.c`):
  DF17 frames at known 12 MHz ticks with noise and a random carrier offset,
  for sensitivity and MLAT timestamp accuracy.
* **Live**: the same Pi 4 receiver feeding FlightAware, including MLAT.

## Signal path

1. **Input conversion** (`iqdata_conv_uc8/sc16/sc16q11`, `convert.c`).
   Converts the SDR samples to 12-bit signed I/Q and removes DC in the same
   pass. The DC filter holds its correction for blocks of 16 samples
   (`IQ_DC_BLOCK`) and updates a running mean with a time constant of 4096
   samples, so SIMD versions need no serial loop. The filter state lives in
   `convert.c` (`converter_state`) and is passed in and out (`dc_in`,
   `dc_out`), never static, which also keeps benchmark verification
   repeatable. Mean level and power are computed after DC removal. The FIFO
   holds I/Q (`magbuf_sample()`), so later stages, such as a coherent
   demodulator, can reuse it.
2. **Upsampling** (`upsample12`). A 5-phase, 8-tap polyphase FIR takes
   2.4 MHz I/Q to 12 MHz, rounds by `UPSAMPLE_SHIFT` and takes the exact
   integer magnitude, `floor(sqrt(I² + Q²))`. Output is stored as 6 rows, one
   per 12 MHz tick within a 0.5 µs chip.
3. **Preamble search** (`locate_preamble`). The metric
   `512·minp / (2·minp + gap + 1)` at every 12 MHz position, from the
   smallest of the 4 preamble pulses and the sum of the 12 gaps. Below
   `PEAMBLE_THRES` it is set to 0 without dividing, which skips most noise.
4. **Candidates** (`demod_2400x.c`). Peaks are scanned 8 metric bytes at a
   time (a SWAR test against the threshold). For each peak, the first 56 bits
   are decoded at the peak and its two neighbours, and the best quality wins
   (`bestPhase`). The second 56 bits are always decoded too, because CRC
   correction can change the DF field, including the bit that selects a
   56- or 112-bit frame. A long frame other than DF17/18 is rejected if its
   second half quality is low.
5. **Bit decoding** (`decode_bits`). Bit = 1 if the first chip of the pair is
   larger. Quality = 256·Σ|a−b| / Σ(a+b).
6. **Output** (`emitFrame`). Scored and decoded by the stock code. Then:
   * Fine arrival time: the decoded chips are correlated with the signal at
     `pos−2 … pos+2` and the best offset is used for the MLAT timestamp.
   * Signal level from the I/Q, with each sample's power capped at full
     scale as the stock magnitude kernels do.
   * Adaptive gain gets 2.4 MHz magnitudes and each decoded message, fed the
     same way as the stock demodulator (only when adaptive gain is on).

## Decisions and why

| Change | Why | Effect |
|---|---|---|
| DC removal in the converter | One pass over the data; FIFO keeps clean I/Q for later use | Fixed an overlap write-offset bug on the way; 3325 → 3332 messages on the SDRplay recording |
| Always decode 112 bits | CRC correction can fix the DF bits, so length is only known after scoring | No lost long frames |
| Skip the metric divide below threshold | `floor(n/d) ≥ T ⇔ n ≥ T·d`; peaks below T are never used | ~1570 → ~1490 ms (SDRplay recording) |
| AVX2 divide/sqrt corrected both ways | `-ffast-math` lets the compiler use reciprocal estimates; results were one too small on real Intel hardware | Validation passes on Intel |
| Fine arrival time | Bit quality is flat over ~3 ticks on a clean pulse, so `pos` was biased (0 to +2 ticks by phase) and jittery (σ 0.88) | Exact on synthetic frames; MLAT synced with ~260 receivers live |
| Least-squares coefficients | The windowed sinc had uneven delay and gain across phases | Weak synthetic frames 1205 → 1405 |
| All phases sum to 128 | Less rounding error; still fits the overflow bound | Delay error 0.046 → 0.028 samples |
| `DATA_THRES` 80 → 65 | More messages, no extra CPU or junk tracks | +0.2% messages |
| `UPSAMPLE_SHIFT` 5 → 4 | One more bit for weak signals; power still < 2³¹ | Results unchanged within noise |
| Mirror-symmetric upsampler (NEON, AVX2) | Phase 4−c is phase c reversed, so pairs share sums `x[k] ± x[7−k]`: half the multiplies | `upsample12` 2.48 → 1.54 ms per call |
| `decode_bits` starch kernel | Called 5 times per candidate, ~2M candidates a minute | 82 → 26 ns per call; demod 4.8 → 4.3 s |
| Adaptive gain fed | It was not called at all by the new demodulator | Same gain decisions as stock |
| Signal power capped at full scale | Clipped frames reported up to +1.1 dBFS | RSSI matches stock (mean +0.01 dB, σ 0.2 dB) |
| Demodulate in 32768-sample pieces (`DEMOD_BLOCK`) | The 12 MHz data for a whole FIFO buffer is 1.3 MB; a piece's fits in L2 between upsample, preamble search and decodes. Done inside `demodulate2400x()`, as the drivers need the FIFO buffers to stay the size of their transfers | Output byte-identical; CPU −4.4% (RTL-SDR), −8.4% (SDRplay) on the Pi 5 |

### Overflow bound

Input is 12-bit (|x| ≤ 2048). Each filter output is at most
`2048·UPSAMPLE_ABS_SUM >> UPSAMPLE_SHIFT`, and the power I² + Q² must stay
below 2³¹ for signed 32-bit arithmetic and the exact square root. With
`UPSAMPLE_ABS_SUM` 244 (the largest Σ|tap|, printed by `sim/ls_coeffs.py`) and
shift 4, the worst case is 1.95·10⁹. A `_Static_assert` in
`demod2400x_coeffs.h` enforces it, and the benchmark runs worst-case input.

### Bit-exactness

Generic is the reference; NEON and AVX2 must match it exactly, which
`starch-benchmark -V` checks. That means exact methods only:

* Square root: float sqrt, then integer checks in both directions.
* Metric divide: float divide, then integer checks in both directions.
* The mirror trick is exact integer arithmetic. Generic keeps the direct form.

## Tried and dropped

* **2.4 MHz pre-detector gate** (skip upsampling where nothing is
  happening): real frames occur at 0.75–1.25 times the mean magnitude, so a
  gate cheap enough to help misses frames.
* **Sharpened upsampler** (`coeffs.py` parameter `a`): amplifies noise more
  than it sharpens pulses; weak frames 1205 → 616 at `a = 0.1`.
* **Lower `PEAMBLE_THRES`**: 50 → 40 gives +1.2% messages for +64% CPU;
  at 30, false aircraft appear.
* **Smaller FIFO buffers** (32768 × 48 instead of 131072 × 12): same gain
  as `DEMOD_BLOCK` with `--ifile`, but the RTL-SDR, HackRF, LimeSDR and
  bladeRF drivers put each transfer (up to 131072 samples) into one buffer
  and drop what doesn't fit, and the stock demodulator would change too.
* **Upsampling only around detected preambles**: the 6-phase metric needs
  the full 12 MHz data, and a 2.4 MHz four-pulse pre-test passes 42% of
  positions (95% of 64-sample blocks); upsampling on demand for the decodes
  alone is break-even.
* **Other preamble metrics** (fewer candidates for the same messages, to
  save CPU). Tested on a 60 s RTL-SDR and a 60 s SDRplay recording:
  * Noise-floor bias in the denominator (25th percentile of the buffer):
    +0.4% messages at a third fewer candidates on both, but needs a floor
    estimate and tuning constants.
  * Sum of the 4 pulses instead of the smallest, or the second smallest:
    fewer messages at any threshold.
  * Gap sum less its largest 1 or 2 chips (stray pulses from other replies
    only add energy): up to +3.3% messages on the SDRplay, none on the
    RTL-SDR, and the detector costs ~30% more. Fewer candidates save little
    (~0.4 µs each), so no CPU gain without losing messages.
* **Soft-decision bit flipping** (flip the weakest bits when the CRC fails):
  `--fix` already repairs any single bit, so this only helps 2-bit errors,
  at a false-decode risk. Not pursued.

## Testing

* `make test`: `cprtests` and `starch-benchmark -V`. Each kernel is checked
  against generic, including full-range worst-case input for `upsample12`
  and a 64-bit check of the generic upsampler.
* Cross-builds: ARMv7 (`ARCH=armv7l CC=arm-linux-gnueabihf-gcc
  LDFLAGS=-static`, runs natively on AArch64) and x86
  (`ARCH=x86_64 CC=x86_64-linux-gnu-gcc`, run with `qemu-x86_64 -cpu max`,
  correctness only; QEMU computes reciprocal estimates exactly, so it can't
  reproduce `-ffast-math` problems).
* Real data: stock and DEMOD_2400X on the same capture, compared message by
  message (contents, RSSI, timestamps).
* Synthetic data: `sim/gen_timing.c`. Frames need a carrier offset: at
  exactly 0 Hz the DC filter removes part of the signal.
* Adaptive gain with `--ifile` needs a gain table; for testing, a fake R820T
  table was patched into a scratch copy of `sdr.c` (not committed).

## For reviewers

### Where to start

* `demod_2400x.c`: the demodulator (buffer handling, peak search, scoring,
  fine timing, signal level) and the hook for the coherent pass.
* `dsp/impl/demod2400x.generic.c`: the plain C reference of every kernel.
  The NEON and AVX2 files must give exactly the same results.
* `demod_2400x_coherent.c/.h` and `dsp/helpers/demod2400x_coherent.h`: the
  coherent pass, described in [COHERENT-demod2400x.md](COHERENT-demod2400x.md).

### Changes to existing upstream files

Everything is inside `#ifdef DEMOD_2400X` (or `DEMOD_2400X_COHERENT`) unless
stated.

| File | Change | Why |
|---|---|---|
| `Makefile` | Adds `-DDEMOD_2400X -DUSE_MAGNITUDE -DDEMOD_2400X_COHERENT`, the new objects, default version `11.2x` | Builds the new demodulator in |
| `dump1090.c` | Options `--no-demod2400x` and `--no-coherent` (both on by default); `fifo_create()` with I/Q or magnitudes; `demodulate2400x()` or `demodulate2400()` | Chooses the demodulator when running |
| `dump1090.h` | `Modes.demod_2400x`, `Modes.coherent`; includes `demod_2400x.h` | |
| `fifo.c/h` | Not under `#ifdef`: `fifo_create()` takes `dtype_iq`, and the buffers and overlap copy are twice as large when it is 1; `magbuf_sample()` addresses sample n (`fifo_values_per_sample` is 1 for magnitudes, as before) | The FIFO carries interleaved I/Q for DEMOD_2400X, so it (and the coherent pass) can use the phase |
| `sdr_*.c` | Not under `#ifdef`, one line each: write at `magbuf_sample(outbuf, overlap)` (the same address as before with magnitudes) | Same, for every SDR driver and `--ifile` |
| `convert.c` | I/Q converters (`convert_iq_*`, calling the `iqdata_conv_*` kernels: 12-bit signed I/Q with DC removed) beside the stock ones; `init_converter()` picks by `Modes.demod_2400x` | DC is removed once, in the converter, so the FIFO holds clean I/Q |
| `demod_2400.c` | `demodulate2400AC()` computes magnitudes from the I/Q first (`getMagnitude()`) when the FIFO holds I/Q | Mode A/C still uses the stock demodulator, which expects magnitudes |
| `stats.c/h`, `net_io.c` | Three coherent counters; `--stats` lines and `"coherent":{...}` in `stats.json` `local`, only while the coherent pass runs | Live monitoring of the coherent pass. Consumers that ignore unknown fields are unaffected |
| `dsp/starchgen.py` | Seven new starch functions | The kernels |
| `wisdom.*` | Kernel order for the new functions | starch picks the fastest version |
| `README.md` | A short section pointing at these documents | |

`--no-demod2400x` gives byte-identical output to upstream `dev`'s stock build on
the 60 s capture (checked with `--raw`), and so does a build without the defines,
which then doesn't have the two options.

### Generated and derived files

* `dsp/generated/*` is produced by `python3 dsp/starchgen.py .`; review
  `dsp/starchgen.py` and `dsp/impl/*` instead.
* `dsp/helpers/demod2400x_coeffs.h` (upsampler taps) comes from
  `sim/ls_coeffs.py`; `dsp/helpers/demod2400x_coh_tables.h` (phasors) from
  `sim/coh_tables.py`.

### How to check it

* `make test`: CPR tests, and `starch-benchmark -V`, which checks every
  NEON and AVX2 kernel against generic (including worst-case input).
* Decode a recording with stock and with this branch and compare (see
  [Testing](#testing)); `--throttle` for positions and tracks.

### Open points

* DEMOD_2400X and the coherent pass are on by default (`--no-demod2400x`,
  `--no-coherent` to switch off). Upstream may prefer stock as the default.
* The default version string is `11.2x`; upstream's is `unknown` (packages
  set it).
* `--dcfilter` is still not supported (as upstream), although the converters
  now remove DC anyway.
* Mode A/C (`--modeac`) with DEMOD_2400X goes through `getMagnitude()`, which
  keeps a static buffer; it has had little testing.
* Statics: `getMagnitude()`'s buffer, the coherent statistics and sample-dump
  state in `demod_2400x.c`, and two diagnostic globals in
  `demod_2400x_coherent.c`. All are used from the demodulator thread only.
* The stock demodulator crashes on SDRplay RSP1B input; not addressed here.

## Raspberry Pi 3B

A Pi 3 Model B (Rev 1.2, Cortex-A53 at 1.2 GHz) running Volumio, a 32-bit OS,
so the ARMv7 NEON kernels. A static ARMv7 build of the dev branch, decoding
the first 20 s of the Pi 4 capture (`--fix`); the messages are identical to
the Pi 5's.

| Run | Messages | Demodulation per 20 s | One core |
|---|---|---|---|
| `--no-demod2400x` (stock) | 24,345 | 2.8 s | 14% |
| `--no-coherent` | 28,051 | 11.0 s | 55% |
| default (DEMOD_2400X + coherent) | 30,971 | 12.6 s | 63% |

Kernel times on the Pi 3B (NEON) against the Pi 5 (aarch64 NEON):

| Kernel | Pi 3B | Pi 5 |
|---|---|---|
| `upsample12` | 12.8 ms | 1.54 ms |
| `locate_preamble` | 15.6 / 10.5 ms | 1.64 / 1.22 ms |
| `coherent_decode` (full / stopped at the gate) | 11.8 / 5.5 µs | 1.80 / 0.81 µs |

So a Pi 3B can run it live, at roughly 65–70% of one core once USB input and
the rest are added (on the Pi 4 they add about 4 points). The upsampler and
preamble search, about 80% of it, cost the same whatever the traffic. The
risk is thermal throttling: at 600 MHz it could not keep up, so a heatsink
is advisable. With `--no-coherent` it is about 55%.

The static `dump1090` for ARMv7 needs `-lncurses`, which the cross-compiler
doesn't have; ncurses is only used by interactive mode, so the test build
linked a stub library (not committed).

## Gotchas

* `dsp/generated/*.o` don't depend on `dsp/helpers/*.h`. After a header edit,
  run `rm dsp/generated/*.o`, or generic keeps old code and verification
  compares against stale results.
* After adding a starch function or changing a wisdom file, rerun
  `python3 dsp/starchgen.py .`; the wisdom order is compiled into
  `dispatcher.c`.
* After changing the version string, run `make clean`: objects don't depend
  on the Makefile.
* starch finds `STARCH_IMPL` by regex, so implementations can't be generated
  by macros.
* starch doesn't pass a kernel's return value to its verifier, so
  `decode_bits` returns its quality through a pointer.

## Possible next steps

* Coherent demodulation as a fallback when magnitude demodulation fails,
  reusing the I/Q kept in the FIFO: done (`-DDEMOD_2400X_COHERENT`), see
  [COHERENT-demod2400x.md](COHERENT-demod2400x.md).
* An early DF check before decoding the second half and scoring.
