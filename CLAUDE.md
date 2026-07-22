# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

This is a research fork of **CLASS** (Cosmic Linear Anisotropy Solving System), the Boltzmann
code for cosmological background/perturbation/power-spectrum computations
(upstream: https://github.com/lesgourg/class_public). The fork adds a **GBH** ("Generalized
Boltzmann Hierarchy") treatment of massive neutrinos / ncdm species, developed on top of vanilla
CLASS. See `README_GBH.md` for the running design notes, TODO list, and file-tagging convention
for this feature — read it before touching any GBH or ncdm code.

The GBH extension implements the physics from two papers by Caio Nascimento, kept locally in
`GBH_papers/`: `Nascimento_2021_Generalized_Boltzmann_hierarchy_massive_neutrinos.pdf`
(arXiv:2104.00703, the GBH hierarchy itself) and
`Nascimento_2023_fluid_approximation_massive_neutrinos.pdf` (arXiv:2303.09580, the fluid
approximation behind the `/*ncdm_caio_fa*/`-tagged code). Consult these for the derivations
behind any GBH/ncdm math before changing it.

The current branch (`ncdm_horizon_test`) is adding `pba->ncdm_horizon`, a comoving free-streaming
horizon computed for ncdm/GBH species during background integration.

## Build

```bash
make class      # build the standalone C binary `./class`
make             # also builds libclass.a and the `classy` python wrapper (pip install .)
make -j          # parallel build
make clean       # remove build/ dir, libclass.a, python/classy.c, python/build
```

Compiler/optimization/OpenMP flags are configured at the top of `Makefile` (default compiler
`gcc`, `CPP = g++ --std=c++11 ...`). Object files land in `build/`; `vpath` maps `source/`,
`tools/`, `main/`, `test/` for compilation.

Run the code against an `.ini` file (and optionally a `.pre` precision file):

```bash
./class explanatory.ini              # reference input file documenting every parameter
./class default.ini
./class default_gbh.ini              # example input enabling the GBH extension
./class test.ini cl_permille.pre     # ini + precision-override file
```

## Tests

C-level unit binaries (one per module) build via Makefile targets and run standalone, no arguments:

```bash
make test_background && ./test_background
make test_thermodynamics && ./test_thermodynamics
make test_perturbations && ./test_perturbations
make test_fourier && ./test_fourier
make test_transfer && ./test_transfer
make test_harmonic && ./test_harmonic
make test_hyperspherical && ./test_hyperspherical
make test_loops && ./test_loops           # and test_loops_omp for the OpenMP variant
```

Python wrapper tests use `nosetests` against many cosmological scenarios (`python/test_class.py`).
Key env vars it reads:

```bash
cd python
TEST_LEVEL=1 nosetests test_class.py            # 0-3, higher = more scenarios/scans
COMPARE_OUTPUT_GAUGE=1 nosetests -a test_scenario test_class.py   # sync vs Newtonian gauge diff
COMPARE_OUTPUT_REF=1 nosetests test_class.py    # compare against a reference `classyref`
```

The wrapper must be rebuilt (`make` or `make classy`) after any C-side struct/field change before
python tests will reflect it.

## Architecture

CLASS is a strict pipeline of modules, each with a `struct` holding its inputs/outputs, wired
together in `main/class.c` (and mirrored by `python/classy.pyx` for the wrapper). Every module
follows the same `<module>_init(...)` / `<module>_free(...)` contract and only reads structs
produced by earlier stages — never later ones:

```
input -> background -> thermodynamics -> perturbations -> primordial
      -> fourier -> transfer -> harmonic -> lensing -> distortions -> output
```

- `include/*.h` declares each module's struct and public functions; `source/*.c` implements them.
  `include/class.h` fixes the include order and is the map of the whole pipeline.
- `tools/` + matching `include/*.h` hold generic numerics used across modules: quadrature,
  hermite/spline interpolation (`arrays.c`), sparse solver, the `ndf15`/`rkck` ODE evolvers,
  the hyperspherical Bessel code, and the input-file parser.
- `external/` vendors optional physics engines selected at compile time in the Makefile
  (HyRec2020/RecfastCLASS for recombination, Halofit/HMcode for nonlinear P(k), heating/injection,
  distortions, GBH tables). `external/GBH/table_*.dat` are precomputed `w_n(x)` and `rho(x)`
  tables consumed when `gbh_use_table = 1`.
- Error handling is macro-based, not exceptions: `class_call(function, err_msg, err_out)`,
  `class_test(condition, err_msg, ...)`, `class_alloc`/`class_calloc`/`class_realloc` (see
  `include/common.h`). Every function returns `_SUCCESS_`/`_FAILURE_` and writes into a
  per-struct `ErrorMsg error_message` field on failure — follow this pattern in new code rather
  than introducing ad hoc error paths.
- Struct fields are allocated/indexed through `class_define_index(pba->index_bg_X, condition,
  index_bg, size)`-style macros, gated on `_TRUE_`/`_FALSE_` flags (e.g. `pba->has_gbh`,
  `pba->has_ncdm`). When adding a new physical species/quantity, follow this same
  index-then-conditionally-allocate pattern rather than unconditional allocation.
- `cpp/` provides a small C++ wrapper (`ClassEngine`) around the C API for external C++ callers;
  `python/` (`classy.pyx`/`cclassy.pxd`) is the Cython wrapper that `pip install .` builds into
  the `classy` python module.

### GBH / ncdm extension conventions (this fork)

- GBH-specific background code is bracketed with `/*GBH_bg_start*/ ... /*GBH_bg_end*/` comments
  (one-line additions tagged `//GBH_bg`); GBH-specific perturbation-related background code with
  `/*GBH_pt_start*/ ... /*GBH_pt_end*/` (`//GBH_pt`). These appear across `background.c`,
  `background.h`, `precisions.h`, `input.c`, `perturbations.c`. Grep for these tags before
  editing to find every place a change needs to be mirrored.
- `/*ncdm_caio_fa*/` tags mark code implementing "Caio's" fluid-approximation scheme for ncdm
  (an alternative closure for the ncdm hierarchy at high `x`), distinct from the GBH hierarchy
  itself — don't conflate the two when editing fluid-approximation logic.
- All GBH computation must stay conditioned on `pba->has_gbh` (set from `Omega0_gbh`/GBH input
  parameters in `input.c`); `README_GBH.md`'s TODO list tracks known-fragile spots (table
  interpolation bounds, `n_max_gbh` vs. table size checks, shooting for `rho0_gbh`, etc.) — check
  it for context before "fixing" something that looks like a bug there.
- Precomputed tables (`external/GBH/table_rho.dat`, `table_omegas.dat`) are keyed on
  `x = m*a/T` and read via `gbh_use_table`/`gbh_w_table` input parameters (`n_max_gbh_table`
  columns); `default_gbh.ini` documents these.

## Input files

`.ini` files (`explanatory.ini` is the canonical, fully-commented reference; `default.ini` and
`default_gbh.ini` are trimmed starting points) drive all runtime physics/output options. `.pre`
files (`cl_ref.pre`, `cl_permille.pre`) override precision-only parameters and are passed as a
second CLI argument. New input parameters get parsed in `source/input.c` and declared in the
relevant module's struct in `include/*.h`.
