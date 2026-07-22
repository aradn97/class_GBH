# Task 7 plan — Redefine GBH mass/degeneracy fields in `background.h`

Scope: `include/background.h` only. No `.c` consumers are touched here (those are
Tasks #8–10). The struct will compile standalone, but the rest of the codebase
will **not** compile again until Task #8 (input.c parsing) and Tasks #9–10
(background.c/perturbations.c consumers) land, since they still reference the
old scalar `pba->M_gbh`/`pba->N_gbh`. That's expected and tracked, not a bug in
this step.

## Current state (`include/background.h:105-115`)

```c
/*GBH_bg_start*/
int last_index_gbh_rho;
int last_index_gbh_w;
int gbh_use_table;
int n_max_gbh;
int gbh_init_condition_integrate;
double M_gbh;             /**< mass of gbh species in eV */
double N_gbh;             /**< number of gbh species */
double T0_gbh;            /**< present temperature of gbh / T_cmb */
int gbh_fluid_approximation;
/*GBH_bg_end*/
```

Two problems being fixed here, beyond the rename itself:
- `M_gbh`'s doc comment says "mass ... in eV" — that's already wrong today.
  `input.c` actually stores the *dimensionless* `m/T` ratio in it
  (`M_gbh = param1/_k_B_*_eV_/T0_gbh/T_cmb`); the eV value itself is never kept
  around. The new `m_gbh_in_eV` field is what will actually hold eV.
- `N_gbh`'s doc comment ("number of gbh species") describes what we *want* it
  to mean now, but today it's read as a scalar degeneracy-like multiplier —
  this task makes the comment true.

## Target state

Modeled directly on the ncdm block a few lines above it
(`background.h:83-98`), reusing the exact same field-pair idiom
(`array` + `array_default`):

```c
/*GBH_bg_start*/
int last_index_gbh_rho;
int last_index_gbh_w;
int gbh_use_table;
int n_max_gbh;
int gbh_init_condition_integrate;
int N_gbh;                          /**< number of distinguishable gbh species. Default 0 (off),
                                          mirrors N_ncdm exactly: the GBH input-parsing block in
                                          input.c only runs if N_gbh > 0. */
double * m_gbh_in_eV;                /**< list of gbh species masses in eV, one per species */
double * M_gbh;                      /**< vector of dimensionless mass ratios m_gbh/T_gbh per species,
                                          inferred from m_gbh_in_eV (mirrors M_ncdm) */
double * deg_gbh, deg_gbh_default;   /**< vector of degeneracy parameters per gbh species, and its
                                          default value (mirrors deg_ncdm/deg_ncdm_default, default 1.) */
double T0_gbh;                       /**< present temperature of gbh / T_cmb (stays a single shared
                                          value across species for now — not part of this task's scope) */
int gbh_fluid_approximation;
/*GBH_bg_end*/
```

Field-by-field mapping:

| Old | New | Type change | Notes |
|---|---|---|---|
| `double N_gbh` | `int N_gbh` | scalar → scalar, `double`→`int` | Becomes a real species counter, default `0`, exactly like `N_ncdm` (`background.h:83`) |
| `double M_gbh` | `double * M_gbh` | scalar → array | Now holds the dimensionless ratio *per species*; mirrors `M_ncdm` (`background.h:90`) |
| *(none)* | `double * m_gbh_in_eV` | new array | User-facing eV mass per species; mirrors `m_ncdm_in_eV` (`background.h:91`) |
| *(none)* | `double * deg_gbh` + `double deg_gbh_default` | new array + scalar | Degeneracy per species, default `1.0`; mirrors `deg_ncdm`/`deg_ncdm_default` (`background.h:98`) |
| `double T0_gbh` | unchanged | — | Out of scope for Task 7 per the earlier scoping discussion — stays a single shared value |

`last_index_gbh_rho`, `last_index_gbh_w`, `gbh_use_table`, `n_max_gbh`,
`gbh_init_condition_integrate`, `gbh_fluid_approximation` are untouched by
this task — they're single-species-table/precision machinery covered by the
later (not-yet-tracked) phase 2+ work.

## Edit sequencing

1. Read `include/background.h:105-115` (already done above) to confirm no
   other code in the same header block references these fields by their old
   names/types (it doesn't — they're declaration-only here).
2. Single `Edit` replacing the block above, old → new.
3. Grep the whole repo for `pba->N_gbh`, `pba->M_gbh` to produce the exact
   consumer list Tasks #8–10 need (this becomes their starting point, not
   something to fix in this task).
4. No build attempt yet — the repo is expected to fail to compile until
   Task #8 lands, since `source/input.c` still assigns into the old scalar
   fields. This task's own "done" condition is just: header compiles in
   isolation (no self-referential type errors) and every other GBH field
   comment/ordering stays consistent with the ncdm block's style.

## Out of scope / explicitly deferred

- `T0_gbh` staying shared (not split per-species) — confirmed out of scope
  earlier in this conversation.
- Any change to `.c` files — Tasks #8 (input.c), #9 (background.c), #10
  (perturbations.c).
- `Omega0_gbh` staying a single scalar (not `Omega0_gbh[]` +
  `Omega0_gbh_tot` like ncdm) — not part of the three fields the user asked
  to restructure; would be part of the larger not-yet-tracked phase 2+
  background-index-block work.
