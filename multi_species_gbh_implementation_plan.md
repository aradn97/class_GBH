# Multi-species GBH implementation plan (Tasks #13–19)

**Prerequisite:** Tasks #7–12 (see `task7_plan_gbh_species_fields.md`) must be
complete first. Those tasks restructure `N_gbh`/`deg_gbh`/`m_gbh_in_eV`/`M_gbh`
into properly-typed per-species fields, but every consumer still hardcodes
index `[0]` — single-species behavior is preserved exactly, nothing loops yet.
**This document is what makes `N_gbh > 1` actually evolve multiple species.**

All file:line references below are current as of commit `ba743dc5`
(`ncdm_horizon_test` branch, after the `(.)horizon_gbh` column-title fix was
pulled). **Tasks #7–12 will shift most of these line numbers** — re-grep the
relevant symbol before editing rather than trusting the numbers blindly; the
surrounding code and control flow described here won't change, just the exact
line count.

Throughout, "the ncdm precedent" means: CLASS's existing `N_ncdm`-species loop
pattern is the template to mirror. Where GBH has no ncdm precedent to copy
(the horizon/truncation-order machinery), that's called out explicitly.

---

## Task #13 — Per-species horizon ODE + truncation-order arrays

**Files:** `include/background.h`, `source/background.c`,
`include/perturbations.h`, `source/perturbations.c`

### Why this has to come first

`pba->gbh_horizon` (a single scalar today) directly drives the perturbation
truncation order (`ppv->n_max_gbh`/`l_max_gbh`) for *every* k-mode
(`perturbations.c:3919`, `x = k * pba->gbh_horizon`). Everything downstream —
state-vector sizing (#15), the hierarchy RHS (#17) — depends on each species
having its own truncation order, which depends on each species having its own
horizon. Get this wrong and every later task has to be redone.

**No ncdm precedent exists for this piece.** Even ncdm's own analogous horizon
(`pba->ncdm_horizon`) hardcodes `pba->M_ncdm[0]` (background.c:3427) and is
never consumed downstream — it's a dead diagnostic. GBH's horizon is
load-bearing, so this is new design, not a mirror of existing code.

### Current state — `include/background.h`

```
151:  double gbh_horizon; /**< conformal free streaming length of GBH species in Mpc */
...
210:  int index_bg_w_min1_gbh;
211:  int index_bg_rho_gbh;
212:  int index_bg_w_gbh;
213:  int index_bg_app_horizon_gbh;     /**< approximate comoving horizon distance of species in GBH*/   //GBH_pt
...
294:  int index_bi_app_horizon_gbh;  /**< {C} approximate horizon distance in Mpc for gbh */
```

### Target state — `include/background.h`

```c
double * gbh_horizon; /**< conformal free streaming length of each GBH species, in Mpc. Size N_gbh. */
...
int index_bg_app_horizon_gbh1;  /**< base of an N_gbh-sized block; species k at index_bg_app_horizon_gbh1+k */
...
int index_bi_app_horizon_gbh1;  /**< {C} base of an N_gbh-sized block of approximate horizon distances, one per species */
```

(Naming convention `_gbh1` for the block base mirrors `index_bg_rho_ncdm1` —
see Task #14 for the same idiom applied to `rho`/`w_min1`.)

### Current state — `source/background.c`

Index reservation (`background_indices`, line 1243):
```c
class_define_index(pba->index_bg_app_horizon_gbh,pba->has_gbh,index_bg,1);
```
line 1354:
```c
class_define_index(pba->index_bi_app_horizon_gbh,pba->has_gbh,index_bi,1);
```

Initial condition (line 3042, inside the background ODE IC-setting function):
```c
pvecback_integration[pba->index_bi_app_horizon_gbh] = pvecback_integration[pba->index_bi_tau]; //equal to tau initially
```

The ODE itself (`background_derivs`, line 3422):
```c
if (pba->has_gbh == _TRUE_){
  dy[pba->index_bi_app_horizon_gbh] = 1. / (a * H * sqrt(1. + pow(pba->M_gbh * a / 3., 2.)));
}
```
Recall from the earlier physics review: this is a faithful ODE form of the
2021 paper's Eq. (A6) `T := ∫dτ(q/ε)` at `q=3T₀`, converted via CLASS's
standard `dτ/d(ln a)=1/(aH)`. **Preserve this derivation exactly per species**
— just substitute `M_gbh` (scalar) → `M_gbh[k]` (per-species) in the same
formula.

Retrieval at end of background solve (line 2614):
```c
if (pba->has_gbh == _TRUE_) { //GBH_pt
  pba->gbh_horizon = pvecback_integration[pba->index_bi_app_horizon_gbh];
}
```

Output columns (`background_output_titles` line 3194, `background_output_data`
line 3285) — single column each, need to become per-species repeated columns
(same fix pattern as Task #14's output-column work; do both in the same pass
if convenient since they're adjacent code).

### Target state — `source/background.c`

- Index reservation: `class_define_index(pba->index_bg_app_horizon_gbh1,pba->has_gbh,index_bg,pba->N_gbh);` and the `index_bi` equivalent — straight copy of the ncdm block-reservation idiom (`background.c:1233`).
- IC-setting loop: `for(k=0;k<pba->N_gbh;k++) pvecback_integration[pba->index_bi_app_horizon_gbh1+k] = pvecback_integration[pba->index_bi_tau];`
- `background_derivs`: `for(k=0;k<pba->N_gbh;k++) dy[pba->index_bi_app_horizon_gbh1+k] = 1./(a*H*sqrt(1.+pow(pba->M_gbh[k]*a/3.,2.)));`
- Retrieval: `for(k=0;k<pba->N_gbh;k++) pba->gbh_horizon[k] = pvecback_integration[pba->index_bi_app_horizon_gbh1+k];` — remember `pba->gbh_horizon` must be `class_alloc`'d to size `N_gbh` first (in `background_gbh_init` or nearby, alongside the other per-species array allocations from Task #14).
- Output columns: loop like ncdm's `background.c:3268` pattern (`for(n=0;n<pba->N_gbh;n++){ sprintf(tmp,"(.)app_horizon_gbh_[%d]",n); class_store_columntitle(...); }` and matching `class_store_double` loop).

### Current state — `include/perturbations.h`

```
495:  int n_max_gbh;
496:  int l_max_gbh;
```
(fields of `struct perturbations_vector`)

### Target state — `include/perturbations.h`

```c
int N_gbh;          /**< number of distinct gbh species (copied from pba->N_gbh for convenience, mirrors ppv->N_ncdm) */
int * n_max_gbh;     /**< truncation order in velocity-moment index n, one per species. Size N_gbh. */
int * l_max_gbh;     /**< truncation order in multipole l, one per species. Size N_gbh. */
```
(mirrors `ppv->N_ncdm`/`ppv->l_max_ncdm`/`ppv->q_size_ncdm` exactly, declared
at `perturbations.h:521-523`)

### Current state — `source/perturbations.c:3908-3944` (`perturbations_vector_init`)

```c
if (pba->has_gbh == _TRUE_) {
  x = k * pba->gbh_horizon;
  class_test(ppr->gbh_nl_max_method!=0. && ppr->gbh_nl_max_method!=1., ppt->error_message,
              "gbh_nl_max_method should be either 0 or 1.");
  if(ppr->gbh_nl_max_method==0){
    ppv->n_max_gbh = std::min(static_cast<int>(ceil(pow(x,1.6)/5.)+3.),static_cast<int>(ceil(pow(ppr->gbh_FA_trigger,1.6)/5.)+3.));
    ppv->l_max_gbh = std::min(static_cast<int>(ceil(x/2.)+2.),static_cast<int>(ceil(ppr->gbh_FA_trigger/2.)+2.));
    if(ppv->l_max_gbh<3){ ppv->l_max_gbh = 3; }
    if(ppv->n_max_gbh<1){ ppv->n_max_gbh = 1; }
  }
  else if(ppr->gbh_nl_max_method==1){
    ppv->n_max_gbh = static_cast<int>(ceil(pow(ppr->gbh_FA_trigger,1.6)/5.)+3.);
    ppv->l_max_gbh = static_cast<int>(ceil(ppr->gbh_FA_trigger/2.)+2.);
  }
  class_test(ppv->n_max_gbh + ppv->l_max_gbh > pba->n_max_gbh, ppt->error_message, "...");
}
```

### Target state

Wrap the whole block in `for(k=0;k<pba->N_gbh;k++){ ... }`, replacing every
`x`/`ppv->n_max_gbh`/`ppv->l_max_gbh` with `x_k`/`ppv->n_max_gbh[k]`/
`ppv->l_max_gbh[k]`, and `x_k = k_mode * pba->gbh_horizon[k]` (careful: `k` is
already used as the wavenumber variable name in this function — **use a
different loop variable name, e.g. `species_k`, to avoid shadowing** the
existing `k` (wavenumber) parameter; this is a real footgun since the
existing code already uses `k` for wavenumber throughout
`perturbations_vector_init`).

`ppv->n_max_gbh`/`l_max_gbh` need `class_alloc`'d to size `N_gbh` before the
loop (mirrors `ppv->l_max_ncdm`/`q_size_ncdm` allocation at
`perturbations.c:4091-4092`). The final `class_test` bound-check against
`pba->n_max_gbh` (the background-table size limit, still a single scalar per
the Task #7-12 scoping decision) should become a check against
`max_k(n_max_gbh[k]+l_max_gbh[k])` or, more conservatively, just check each
species individually against the same shared `pba->n_max_gbh` bound.

---

## Task #14 — Convert background-vector indices to species blocks

**Files:** `include/background.h`, `source/background.c`

### Current state — `include/background.h:210-212`

```c
int index_bg_w_min1_gbh;    /**< P_{-1} of species in GBH, from quadrature integration*/   //GBH_bg
int index_bg_rho_gbh;       /**< density of species in GBH*/   //GBH_bg
int index_bg_w_gbh;         /**< equation of state of species in GBH*/   //GBH_bg
```

### Target state

```c
int index_bg_w_min1_gbh1;   /**< base of N_gbh-sized block; species k at index_bg_w_min1_gbh1+k */
int index_bg_rho_gbh1;      /**< base of N_gbh-sized block; species k at index_bg_rho_gbh1+k */
int index_bg_w_gbh1;        /**< base of an N_gbh*(n_max_gbh+1)-sized 2-D block (species-major):
                                  species k, moment n at index_bg_w_gbh1 + k*(n_max_gbh+1) + n */
```

The first two are direct copies of the `index_bg_rho_ncdm1` idiom
(`background.c:1233`, block size `pba->N_gbh`). **`index_bg_w_gbh1` has no
ncdm precedent** — ncdm's per-species background quantities
(`rho_ncdm1`/`p_ncdm1`/`pseudo_p_ncdm1`) are all scalars-per-species, never
their own sub-array, because ncdm doesn't have GBH's extra moment axis. Two
layout choices for the 2-D block:
  - **Species-major (recommended):** `index_bg_w_gbh1 + k*(n_max_gbh+1) + n`. Requires all species to share the same `n_max_gbh` background-table size (already true today — `pba->n_max_gbh` stays a single shared scalar per the Task #7-12 scoping decision), which keeps the stride constant and this simple.
  - Moment-major: `index_bg_w_gbh1 + n*N_gbh + k`. No clear advantage here since consumers always loop by species outermost; don't use this.

### Current state — `source/background.c:1237-1239`

```c
class_define_index(pba->index_bg_w_min1_gbh,pba->has_gbh,index_bg,1); //GBH_bg
class_define_index(pba->index_bg_rho_gbh,pba->has_gbh,index_bg,1); //GBH_bg
class_define_index(pba->index_bg_w_gbh,pba->has_gbh,index_bg,pba->n_max_gbh+1); //GBH_bg
```

### Target state

```c
class_define_index(pba->index_bg_w_min1_gbh1,pba->has_gbh,index_bg,pba->N_gbh);
class_define_index(pba->index_bg_rho_gbh1,pba->has_gbh,index_bg,pba->N_gbh);
class_define_index(pba->index_bg_w_gbh1,pba->has_gbh,index_bg,pba->N_gbh*(pba->n_max_gbh+1));
```

### Current state — `background_functions`, `source/background.c:552-669`

Single `if(pba->has_gbh==_TRUE_){...}` block (both the `gbh_use_table==1`
interpolation path, lines ~557-624, and the quadrature fallback, lines
~625-656) writing into `pvecback[pba->index_bg_rho_gbh]`,
`pvecback[pba->index_bg_w_min1_gbh]`, `pvecback[pba->index_bg_w_gbh+n]` for
one species, then accumulating into `rho_tot`/`p_tot`/`rho_r`/`rho_m` at
lines 658-666.

### Target state

Wrap in `for(species_k=0;species_k<pba->N_gbh;species_k++){...}`, mirroring
the ncdm loop at `background.c:502-547` exactly in structure:
```c
for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {
  class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm], pba->w_ncdm_bg[n_ncdm], pba->q_size_ncdm_bg[n_ncdm],
                                       pba->M_ncdm[n_ncdm], pba->factor_ncdm[n_ncdm], 1./a-1., ...), ...);
  pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
  rho_tot += rho_ncdm;
  ...
}
```
GBH's version becomes: `current_x_value = pba->M_gbh[species_k]*a`, table
lookups use `pba->x_gbh_bg`/`pba->w_gbh_bg` etc — **these stay single shared
tables for now** (per the Task #7-12 scoping note that `T0_gbh`/table paths
are out of scope), so multiple species with different masses all interpolate
against the *same* table, just at different `x` values — this is physically
valid only if all GBH species share the same underlying phase-space shape
(Fermi-Dirac, same `T0_gbh`), which is the current assumption. If species
need genuinely different temperatures/statistics later, the table
infrastructure itself needs per-species tables (out of scope here, flag it
as a known limitation in code comments).
`pvecback[pba->index_bg_rho_gbh1+species_k] = rho_gbh_k`, and the
`index_bg_w_gbh1` 2-D block write becomes
`pvecback[pba->index_bg_w_gbh1 + species_k*(pba->n_max_gbh+1) + n] = ...`.
Accumulate `rho_tot += rho_gbh_k` etc. inside the loop (sum over species,
same as ncdm).

### `background_gbh_init` (`source/background.c:1819-2093`) and `background_gbh_momenta` (`source/background.c:2121-2190`)

`background_gbh_init` currently reads one `w_table`/`rho_table` file pair
(paths from `ppr->gbh_w_file`/`gbh_rho_file`, single fixed strings, not even
`.ini`-configurable) and builds one quadrature grid
(`pba->q_gbh_bg`/`weights_gbh_bg`). Since the table itself stays shared
(previous paragraph), this function's *file-loading* part doesn't need a
species loop — but any per-species derived quantity it computes (if any) does.
`background_gbh_momenta` takes scalar `M`/`factor`/`qvec`/`wvec`/`qsize`
arguments — this function signature can stay as-is and simply be **called
once per species inside the `background_functions` loop above**, passing
`pba->M_gbh[species_k]` and `pba->factor_gbh` (still shared, since
`factor_gbh` derives from `T0_gbh` which is shared) — no signature change
needed, just more call sites.

### Output columns — `background_output_titles`/`background_output_data`

Lines 3182-3195 / 3274-3285. Follow the exact ncdm per-species column loop
pattern already in the same file (`background.c:3268`-area, e.g.
`for(n=0;n<pba->N_ncdm;n++){ sprintf(tmp,"(.)rho_ncdm[%d]",n); class_store_columntitle(...); }`).
Do this for `rho_gbh`, `w_min1_gbh`, the `w_gbh[n]` moment columns (nested
loop: species outer, moment inner), and `app_horizon_gbh` (from Task #13) in
the same pass since they're adjacent code blocks.

---

## Task #15 — Species-loop perturbation state-vector index allocation

**File:** `source/perturbations.c`

### Current state (`perturbations_vector_init`, lines 4116-4127)

```c
if (pba->has_gbh == _TRUE_) {
  if (ppw->approx[ppw->index_gbh_fa] == (int)gbh_fa_off) {
    class_define_index(ppv->index_pt_Delta_gbh,_TRUE_,index_pt,ppv->n_max_gbh);
    class_define_index(ppv->index_pt_Sigma_gbh,_TRUE_,index_pt,(ppv->n_max_gbh)*(ppv->l_max_gbh));
  }
  else{
    class_define_index(ppv->index_pt_delta_gbh,_TRUE_,index_pt,1);
    class_define_index(ppv->index_pt_theta_gbh,_TRUE_,index_pt,1);
    class_define_index(ppv->index_pt_sigma_gbh,_TRUE_,index_pt,1);
  }
}
```

Note this already reads `ppw->approx[ppw->index_gbh_fa]` as a **single**
flag deciding the mode for the (one) species — Task #18 makes this
per-species, so by the time Task #15 lands, expect this `if/else` to already
be inside (or about to become) a per-species loop where each species
independently checks its own FA-vs-exact flag. Sequence Task #15 assuming
Task #18's per-species `index_gbh_fa` block already exists (they're tightly
coupled — read Task #18 before starting #15's implementation, even though
they're numbered/blocked sequentially the other way; consider doing the
`index_gbh_fa` block-conversion piece of #18 first if it simplifies #15).

### The ncdm precedent (`perturbations.c:4085-4110`)

```c
if (pba->has_ncdm == _TRUE_) {
  ppv->index_pt_psi0_ncdm1 = index_pt; /* single base index for the whole multi-species block */
  ppv->N_ncdm = pba->N_ncdm;
  class_alloc(ppv->l_max_ncdm,ppv->N_ncdm*sizeof(double),ppt->error_message);
  class_alloc(ppv->q_size_ncdm,ppv->N_ncdm*sizeof(double),ppt->error_message);

  for (n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++){
    ...
    ppv->l_max_ncdm[n_ncdm] = ppr->l_max_ncdm;
    ppv->q_size_ncdm[n_ncdm] = pba->q_size_ncdm[n_ncdm];
    index_pt += (ppv->l_max_ncdm[n_ncdm]+1)*ppv->q_size_ncdm[n_ncdm];  /* manual running offset */
  }
}
```
`class_define_index` is deliberately **not used** here because each species'
block size varies — same reason GBH needs this pattern once `n_max_gbh[k]`/
`l_max_gbh[k]` (Task #13) differ per species.

### Target state

```c
if (pba->has_gbh == _TRUE_) {
  ppv->N_gbh = pba->N_gbh;
  class_alloc(ppv->index_pt_delta_gbh, ppv->N_gbh*sizeof(int), ppt->error_message); /* or Delta_gbh, per species FA state */
  class_alloc(ppv->index_pt_theta_gbh, ppv->N_gbh*sizeof(int), ppt->error_message);
  class_alloc(ppv->index_pt_sigma_gbh, ppv->N_gbh*sizeof(int), ppt->error_message);
  class_alloc(ppv->index_pt_Delta_gbh, ppv->N_gbh*sizeof(int), ppt->error_message);
  class_alloc(ppv->index_pt_Sigma_gbh, ppv->N_gbh*sizeof(int), ppt->error_message);

  for (species_k = 0; species_k < pba->N_gbh; species_k++){
    if (ppw->approx[ppw->index_gbh_fa+species_k] == (int)gbh_fa_off) {  /* per-species FA flag, see Task #18 */
      ppv->index_pt_Delta_gbh[species_k] = index_pt;
      index_pt += ppv->n_max_gbh[species_k];
      ppv->index_pt_Sigma_gbh[species_k] = index_pt;
      index_pt += ppv->n_max_gbh[species_k]*ppv->l_max_gbh[species_k];
    }
    else{
      ppv->index_pt_delta_gbh[species_k] = index_pt; index_pt += 1;
      ppv->index_pt_theta_gbh[species_k] = index_pt; index_pt += 1;
      ppv->index_pt_sigma_gbh[species_k] = index_pt; index_pt += 1;
    }
  }
}
```

**This changes `index_pt_delta_gbh` etc. from `int` to `int *`** in
`struct perturbations_vector` (`perturbations.h:494-502`) — update the
header alongside this.

### Every consumption site must change

Every place that currently reads `ppv->index_pt_Delta_gbh` (a bare int) as
`y[ppv->index_pt_Delta_gbh+n]` must become
`y[ppv->index_pt_Delta_gbh[species_k]+n]` inside a species loop. This
includes the approximation-switch state-copy blocks at (approximately, will
shift) `perturbations.c:4339-4954` and beyond — there are **many** of these
copy blocks (one per pair of approximation schemes being switched between);
grep `index_pt_Delta_gbh\|index_pt_Sigma_gbh\|index_pt_delta_gbh\|index_pt_theta_gbh\|index_pt_sigma_gbh`
across `perturbations.c` after Task #13/14 land to get the authoritative,
re-numbered list before starting — there were roughly 15-20 distinct call
sites as of this writing (see the grep output embedded in this repo's
conversation history / re-run `grep -n "index_pt_.*_gbh" source/perturbations.c`).

---

## Task #16 — `ppw` GBH workspace arrays + loop stress-energy accumulation

**Files:** `include/perturbations.h`, `source/perturbations.c`

### Current state — `include/perturbations.h`

No `ppw->delta_gbh`/`theta_gbh`/`sigma_gbh` fields exist at all (unlike ncdm's
`ppw->delta_ncdm`/`theta_ncdm`/`shear_ncdm` at `perturbations.h:611-613`).
GBH code reads the state vector directly at each use site using local C
scalars.

### Target state — `include/perturbations.h`

```c
double * delta_gbh; /**< relative density perturbation of each gbh species. Size N_gbh. */
double * theta_gbh; /**< velocity divergence theta of each gbh species. Size N_gbh. */
double * sigma_gbh; /**< shear of each gbh species. Size N_gbh. */
```
Allocate these in `perturbations_workspace_init` (near `perturbations.c:2805`
where `ppw->index_gbh_fa` is set up) to size `pba->N_gbh`.

### Current state — `perturbations_total_stress_energy`, `source/perturbations.c:7734-7768`

```c
double rho_gbh,p_gbh,rho_plus_p_gbh,pseudo_p_gbh,w_gbh,ca2_gbh,cg2_gbh,delta_gbh,theta_gbh,sigma_gbh; // line 7382, local scalars

if (pba->has_gbh == _TRUE_) {
  rho_gbh = ppw->pvecback[pba->index_bg_rho_gbh];
  p_gbh = ppw->pvecback[pba->index_bg_w_gbh+1]*rho_gbh;
  rho_plus_p_gbh = rho_gbh + p_gbh;

  if (ppw->approx[ppw->index_gbh_fa] == (int)gbh_fa_on){
    delta_gbh = y[ppw->pv->index_pt_delta_gbh];
    theta_gbh = y[ppw->pv->index_pt_theta_gbh];
    w_gbh = p_gbh/rho_gbh;
    pseudo_p_gbh = ppw->pvecback[pba->index_bg_w_gbh+2]*rho_gbh;
    sigma_gbh = y[ppw->pv->index_pt_sigma_gbh];
    cg2_gbh = w_gbh*(1.0-1.0/(3.0+3.0*w_gbh)*(3.0*w_gbh-2.0+pseudo_p_gbh/p_gbh));
    ppw->delta_p += cg2_gbh*rho_gbh*delta_gbh;
  }
  else{
    delta_gbh = 3.*y[ppw->pv->index_pt_Delta_gbh];
    theta_gbh = k*y[ppw->pv->index_pt_Sigma_gbh];
    sigma_gbh = y[ppw->pv->index_pt_Sigma_gbh+1];
    ppw->delta_p += rho_gbh*y[ppw->pv->index_pt_Delta_gbh+1];
  }
  ppw->delta_rho += rho_gbh*delta_gbh;
  ppw->rho_plus_p_theta += rho_plus_p_gbh*theta_gbh;
  ppw->rho_plus_p_shear += rho_plus_p_gbh*sigma_gbh;
  ppw->rho_plus_p_tot += rho_plus_p_gbh;
  if (ppt->has_source_delta_m == _TRUE_) {
    delta_rho_m += rho_gbh*delta_gbh;
    rho_m += rho_gbh;
  }
  if ((ppt->has_source_delta_m == _TRUE_) || (ppt->has_source_theta_m == _TRUE_)) {
    rho_plus_p_theta_m += rho_plus_p_gbh*theta_gbh;
    rho_plus_p_m += rho_plus_p_gbh;
  }
}
```

### Target state

Wrap the whole block in `for(species_k=0;species_k<pba->N_gbh;species_k++)`,
replacing `pba->index_bg_rho_gbh` → `pba->index_bg_rho_gbh1+species_k` (Task
#14's block), `ppw->index_gbh_fa` → `ppw->index_gbh_fa+species_k` (Task
#18's per-species block), `index_pt_delta_gbh`/etc → the `[species_k]`-indexed
versions from Task #15, and cache the result into
`ppw->delta_gbh[species_k] = delta_gbh;` etc. before accumulating into the
shared totals (`ppw->delta_rho += ...` stays a plain `+=` accumulator summed
across species, same as the ncdm loop does).

Directly mirrors ncdm's own loop over the same function, `perturbations.c:7646-7728`.

---

## Task #17 — Loop the Boltzmann hierarchy RHS in `perturbations_derivs` over species

**File:** `source/perturbations.c`, inside `perturbations_derivs` (starts
line 9461; the GBH-specific block was previously located around lines
10396-10534 pre-pull — re-grep `/*GBH_pt_start*/` near the end of that
function after Tasks #13-16 land, since line numbers will have shifted).

### Current state (structure, both branches)

```c
if(pba->has_gbh == _TRUE_){
  int N = pv->n_max_gbh+pv->l_max_gbh+1,ll,n;

  if (ppw->approx[ppw->index_gbh_fa] == (int)gbh_fa_on) {
    /* fluid-approximation delta'/theta'/sigma' equations, one species */
    rho_gbh = pvecback[pba->index_bg_rho_gbh];
    w_gbh = pvecback[pba->index_bg_w_gbh+1];
    ...
    dy[pv->index_pt_delta_gbh] = ...;
    dy[pv->index_pt_theta_gbh] = ...;
    dy[pv->index_pt_sigma_gbh] = ...;
  }
  else{
    /* exact GBH hierarchy, one species */
    double * w = &pvecback[pba->index_bg_w_gbh];
    for(n=0;n<pv->n_max_gbh;n++){
      ...
      dy[pv->index_pt_Delta_gbh+n] = ...;
      for(ll=1;ll<pv->l_max_gbh+1;ll++){
        ...
        dy[pv->index_pt_Sigma_gbh+n*pv->l_max_gbh+ll-1] = ...;
      }
    }
  }
}
```

(Full verbatim equations were already verified against the 2021/2023 papers
in the earlier physics review in this conversation — this task is purely
about adding the species loop, **not** re-deriving or changing any of the
physics/algebra inside it.)

### Target state

```c
if(pba->has_gbh == _TRUE_){
  for (species_k=0; species_k<pba->N_gbh; species_k++){
    int N = pv->n_max_gbh[species_k]+pv->l_max_gbh[species_k]+1,ll,n;

    if (ppw->approx[ppw->index_gbh_fa+species_k] == (int)gbh_fa_on) {
      rho_gbh = pvecback[pba->index_bg_rho_gbh1+species_k];
      w_gbh = pvecback[pba->index_bg_w_gbh1 + species_k*(pba->n_max_gbh+1) + 1];
      ...
      dy[pv->index_pt_delta_gbh[species_k]] = ...;
      dy[pv->index_pt_theta_gbh[species_k]] = ...;
      dy[pv->index_pt_sigma_gbh[species_k]] = ...;
    }
    else{
      double * w = &pvecback[pba->index_bg_w_gbh1 + species_k*(pba->n_max_gbh+1)];
      for(n=0;n<pv->n_max_gbh[species_k];n++){
        ...
        dy[pv->index_pt_Delta_gbh[species_k]+n] = ...;
        for(ll=1;ll<pv->l_max_gbh[species_k]+1;ll++){
          ...
          dy[pv->index_pt_Sigma_gbh[species_k]+n*pv->l_max_gbh[species_k]+ll-1] = ...;
        }
      }
    }
  }
}
```

**Every internal reference to `n_max_gbh`/`l_max_gbh`/`index_bg_rho_gbh`/
`index_bg_w_gbh`/`index_pt_*_gbh` inside the (long, ~140-line) equation block
needs the `[species_k]` or `+species_k`/`1+species_k*(...)` treatment** —
this is mechanical but extensive; the safest approach is a careful
find-and-replace pass followed by a diff review against the pre-loop version
with `species_k` fixed to `0`, to confirm the `N_gbh=1` case reduces to
byte-identical RHS values (good regression check before Task #19's full
verification).

---

## Task #18 — Loop initial conditions, FA switch, and source functions over species

**File:** `source/perturbations.c`

### (a) Initial conditions — `perturbations.c:6322-6404`

Two sub-paths, both currently single-species:
- `gbh_init_condition_integrate==_FALSE_` (default, algebraic ur-analogue,
  lines 6329-6337): loops over `n` (moments) using `pba->index_bg_w_gbh`
  already — just needs the outer species loop and per-species `w` pointer
  (same `index_bg_w_gbh1+species_k*(...)` substitution as Task #17).
- `gbh_init_condition_integrate==_TRUE_` (momentum-space quadrature IC,
  lines 6339-6400): uses `pba->M_gbh`, `pba->q_gbh`, `pba->weights_gbh`,
  `pba->dlnf0_dlnq_gbh`, `pba->factor_gbh` directly (all currently scalar/flat
  — after Task #14 `pba->M_gbh` is `M_gbh[k]`; the momentum grids
  `q_gbh`/`weights_gbh`/`dlnf0_dlnq_gbh` stay shared per the table-sharing
  design note in Task #14, so only `M_gbh[species_k]` needs substituting here,
  not a full `q_gbh[k][...]` double-indirection like ncdm's `q_ncdm[n_ncdm][index_q]`).
  Per README_GBH.md TODO #73 this quadrature IC path was already marked for
  eventual removal — worth flagging to the user whether to bother
  multi-species-ifying it at all, or just multi-species-ify the (default,
  faster) algebraic path and leave the quadrature IC path
  single-species/deprecated with a clear error if `N_gbh>1` and
  `gbh_init_condition_integrate=1` are both set.

### (b) FA switch — `perturbations.c:6854-6864`

```c
if (pba->has_gbh == _TRUE_) {
  if (k > ppr->gbh_FA_trigger/ppw->pvecback[pba->index_bg_app_horizon_gbh]) {
    ppw->approx[ppw->index_gbh_fa] = (int)gbh_fa_on;
  }
  else {
    ppw->approx[ppw->index_gbh_fa] = (int)gbh_fa_off;
  }
}
```
`ppw->index_gbh_fa` (`perturbations.c:2805`,
`class_define_index(ppw->index_gbh_fa,pba->has_gbh,index_ap,1);`) needs to
become an `N_gbh`-sized block:
`class_define_index(ppw->index_gbh_fa1,pba->has_gbh,index_ap,pba->N_gbh);`,
and the switch becomes a per-species loop comparing each species against its
own `pvecback[pba->index_bg_app_horizon_gbh1+species_k]` (Task #13's
per-species horizon). **Unlike ncdm's global `index_ap_ncdmfa`** (which stays
a single flag even for multi-species ncdm, since ncdm's FA trigger
(`tau/tau_k`) doesn't depend on species mass) — GBH's trigger is inherently
mass-dependent, so this one genuinely needs to diverge from the ncdm
precedent, as flagged when this task list was first proposed.

### (c) Source functions — `perturbations.c:1429,1444` (index reservation) and `8556-8570`/`8684-8698` (consumption)

```c
class_define_index(ppt->index_tp_delta_gbh,  ppt->has_source_delta_gbh, index_type,1);//GBH_pt
class_define_index(ppt->index_tp_theta_gbh,  ppt->has_source_theta_gbh, index_type,1);//GBH_pt
```
→
```c
class_define_index(ppt->index_tp_delta_gbh1,  ppt->has_source_delta_gbh, index_type,pba->N_gbh);
class_define_index(ppt->index_tp_theta_gbh1,  ppt->has_source_theta_gbh, index_type,pba->N_gbh);
```
mirroring `index_tp_delta_ncdm1` (`perturbations.c:1428`). Consumption sites
(currently single `_set_source_(ppt->index_tp_delta_gbh) = ...` statements)
become loops over `index_tp` from `index_tp_delta_gbh1` to
`index_tp_delta_gbh1+pba->N_gbh`, mirroring ncdm's loop at
`perturbations.c:8548-8553`. Also update the two `class_store_double`/
`class_store_double_or_default` output-transfer-function call sites at
`perturbations.c:480,510,527` similarly (these read `tk[ppt->index_tp_delta_gbh]`
today — become a species loop or an explicit per-species column, matching
however `perturbations_output` already handles `N_ncdm` columns for the ncdm
transfer functions in the same file, for consistency).

---

## Task #19 — End-to-end multi-species build and physics verification

No code changes — this is the acceptance test for Tasks #13-18 together.

**Build:** `make class` with a `.ini` setting `N_gbh=2`, `m_gbh=0.02,0.05`
(two distinct masses, comma-separated per the Task #7-12 `class_read_list_of_doubles_or_default`
parsing), `deg_gbh` left at default (1,1).

**Checks, in order of how much they isolate the bug if they fail:**
1. **Background only** (`./class` with `output=` empty or just background):
   confirm `rho_gbh` total at `z=0` (or any `z`) from the `N_gbh=2` run equals
   the sum of two independent `N_gbh=1` runs (one per mass), to floating-point
   tolerance. This isolates bugs to Task #13/#14 (background side) if it
   fails, before even touching perturbations.
2. **Perturbations, single k, background sources only**: with `N_gbh=1`, m
   matching Task #12's baseline exactly, confirm **zero regression** —
   identical `background.dat`/`perturbations.dat`/Cl/Pk output to the
   pre-Task-#13 baseline. This is the most important check: it confirms the
   whole species-loop refactor is a no-op for the case everyone was already
   relying on.
3. **Perturbations, `N_gbh=2`**: compare Cl/Pk against a manually-summed
   reference (run each mass through `N_gbh=1` separately, since GBH species
   are non-interacting except through gravity — verify the combined metric
   perturbations from the `N_gbh=2` run match what you'd get feeding both
   species' stress-energy into the same Einstein equations by hand, or at
   minimum confirm total matter power spectrum shifts in the expected
   direction/magnitude versus either single-species run alone).
4. Spot-check `gbh_horizon[k]`/`n_max_gbh[k]`/`l_max_gbh[k]` in the two
   species differ sensibly (heavier species → smaller horizon → smaller
   truncation order needed at fixed k, per the 2021 paper's `l_max~x/2`,
   `n_max~x^1.6/5` scaling, `x=k·horizon`).
