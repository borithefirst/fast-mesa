# fast-mesa

A performance fork of [MESA](https://github.com/MESAHub/mesa) r25.12.1 (Modules for Experiments in Stellar
Astrophysics). It runs the same physics with the same inlists and needs **about half the
instructions**. With the optional LTO build it is **2.36× faster** in CPU cycles across a 9-segment benchmark
of stellar evolution, and 1.4–3.8× depending on the phase. The accuracy requirement: every change is either **bit-identical** to stock MESA, or moves results
by far less than MESA's own timestep error.

This is an independent, unofficial modification of MESA, not endorsed by the MESA developers. It is based
on the `25.12.1` release tag, and all changes are on the `fast-mesa` branch. As MESA is, it is distributed under
the GNU LGPL v3 (see `LICENSE`). For science, cite MESA as usual (`CITATIONS.bib`).

## Results

The benchmark suite has 9 segments of stellar evolution, each 40–50 steps from a saved model, run on 1 thread.
"Error" is the difference from a reference run at a quarter of the timestep, i.e. MESA's own timestep error.

| segment | stock instr. | fast-mesa instr. | change | cycles speedup, `MESA_LTO=1` ¹ | log L error, stock | log L error, fast-mesa | fast-mesa − stock, log L |
|---|---|---|---|---|---|---|---|
| 1 Msun main sequence | 2.87e11 | 1.26e11 | −56% | 2.57× | 5.1e-4 | 5.1e-4 | 7.5e-11 |
| 1 Msun subgiant | 1.20e12 | 4.24e11 | −65% | 3.09× | 1.7e-3 | 1.7e-3 | 4.1e-8 |
| 1 Msun red giant branch | 1.03e12 | 5.58e11 | −46% | 2.15× | 1.9e-3 | 1.9e-3 | 3.6e-5 |
| 1 Msun RGB tip | 1.34e12 | 7.76e11 | −42% | 1.85× | 1.6e-3 | 1.6e-3 | 3.5e-4 |
| 1 Msun core He burning | 1.18e12 | 7.22e11 | −39% | 1.85× | 3.0e-7 | 3.0e-7 | 3.3e-11 |
| 1 Msun AGB | 2.32e12 | 7.53e11 | −68% | 3.77× | 3.2e-2 | 3.2e-2 | 2.7e-6 |
| 1 Msun WD cooling | 9.01e11 | 6.55e11 | −27% | 1.58× | — ² | — ² | 7.2e-9 |
| 15 Msun main sequence | 1.63e11 | 1.31e11 | −20% | 1.37× | 1.4e-2 | 1.4e-2 | 7.0e-11 |
| 300 Msun | 4.98e11 | 1.82e11 | −64% | 3.39× | 1.6e-2 | 1.6e-2 | 1.2e-9 |
| **total** | 8.92e12 | 4.33e12 | **−52%** | **2.36×** | | | |

¹ CPU cycles on a quiet machine for the LTO build with every speedup except `MESA_EPSG_LIN`, which was
measured separately. The default (non-LTO) build is measured here in instructions. They understate the time
saved by `MESA_FAST_CRLIBM`, which shortens latency more than it removes instructions.

² The WD segment's timestep is not set by `time_delta_coeff`, so the reference run is identical to stock.

**Checks:**
- **Accuracy:**
  - fast-mesa's error with respect to the quarter-timestep reference matches stock to the digits shown. The one
    exception is log Teff on the RGB: 1.21e-4 → 1.47e-4.
  - The largest deviation from stock, 3.5e-4 dex (RGB tip), is 5× smaller than stock's own timestep error
    there.
- **Bit-identity:**
  - With all six switches set to `0`, every history file is bit-identical to stock MESA 25.12.1 in all 9
    segments.
  - The bit-identical switches (memo, lazy Brunt, fastcr, LTO) were each verified bit-identical on all 9
    segments.

**Whole life** (1 Msun, pre-MS → WD, ~13.5k models, 8 threads), with memo + Skye partials + residual floor:
- **Speed:** −20% instructions.
- **Stages up to the end of core He burning:** the stage-end properties move by less than they do when stock
  MESA is simply run with a 1% smaller timestep.
- **After the AGB:** thermal pulses amplify any perturbation. The final WD mass changes by 0.09% (5e-4 Msun),
  the same order as the 1%-timestep run (2e-4 Msun).

## What changed

Each change is a separate commit, so `git log 25.12.1..fast-mesa` is the full list. The speedups are **on by
default**. Each one turns off with its environment variable set to `0`, and with all six off the results are
bit-identical to MESA 25.12.1.

| switch (env var) | module | what it does | bit-identical? |
|---|---|---|---|
| `MESA_EOS_MEMO` | star | Reuses EOS results that star recomputes with identical inputs (bitwise key). | yes |
| `MESA_LAZY_BRUNT` | star | Skips the Brunt B evaluation in `set_vars` when nothing reads it within the step. | yes |
| `MESA_FAST_CRLIBM` | math | Fast correctly-rounded front end for crlibm `exp_rd`/`log_rz`: double-double evaluation, a rounding decision from the exact residual, crlibm fallback near rounding boundaries, and a 4-entry per-thread cache. The returned bits are crlibm's by construction (verified on 80M arguments). Needs a CPU with FMA; otherwise it silently forwards to crlibm. | yes |
| `MESA_SKYE_DXA` | eos + star | Analytic composition partials from Skye, instead of finite-difference EOS re-evaluations. Newton converges in fewer iterations. | no: the converged solution differs within the solver tolerance |
| `MESA_RESID_FLOOR` | star | Floor-aware energy residual norm. The ULP of L made the gold2 tolerances unreachable in some phases, costing ~10 Newton iterations per step. Default factor 1. | no: it changes when Newton stops |
| `MESA_EPSG_LIN` | star | Linearizes the eps_grav composition term when \|dX\| over a step is below the threshold. Default 1e-8. | no: tested at ≤ 8e-6 dex |

The effect of the non-bit-identical changes on the tracks is in the Results section.

Also included:
- A backport of upstream PR #920 (math_lib without `IEEE_ARITHMETIC`): bit-identical, −25…−28%.
- **Optional link-time optimization.** Build MESA and your run directories with `MESA_LTO=1` in the environment.
  It is bit-identical and gives a further −7…−16%, but every run-directory link then takes ~2 minutes.
- **Experimental, off by default:**
  - `MESA_SCREEN_DXA=1|2`: composition partials of the screening factors in the net.
  - `MESA_PREDICTOR`: a Newton predictor.
  - `MESA_SOLVER_LOG=<file>`: a per-step solver log.

## Install

Same as MESA: install the [MESA SDK](http://user.astro.wisc.edu/~townsend/static.php?ref=mesasdk). Then:

```bash
git clone -b fast-mesa https://github.com/borithefirst/fast-mesa.git
cd fast-mesa
git lfs pull          # MESA's data archives are stored with Git LFS
export MESA_DIR=$PWD
./install
```

Run directories work unchanged. To compare with stock MESA in the same build:

```bash
export MESA_EOS_MEMO=0 MESA_SKYE_DXA=0 MESA_LAZY_BRUNT=0 MESA_FAST_CRLIBM=0 MESA_RESID_FLOOR=0 MESA_EPSG_LIN=0
```

## How it was measured

MESA r25.12.1 with SDK 24.7.1 on WSL2 Ubuntu 24.04, AMD Ryzen 7 7800X3D. The benchmark suite has 9 segments,
each 40–50 steps from a saved model, run on 1 thread and counted with `perf stat` (cycles and instructions):
- 1 Msun: main sequence, subgiant, RGB, RGB tip, core He burning, AGB, WD cooling;
- 15 Msun main sequence;
- 300 Msun.

Accuracy is measured against a reference run of every segment at a quarter of the timestep (`time_delta_coeff
= 0.25`). A whole-life 1 Msun run (pre-MS → WD, ~13.5k models) is compared with stock and with stock at a 1%
smaller timestep.
