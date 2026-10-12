# Why the ncdm DFA was slower than the default fluid approximation, and what was done about it

This note records two separate performance problems found in the `ncdm` dispersive
fluid approximation (`ncdm_fluid_approximation = 4`, `ncdmfa_caio`, tagged
`/*ncdm_caio_fa*/` throughout the source), and the fixes applied to each. It also
records the traps hit along the way, because two of them are easy to fall into
again.

See item 77 of `README_GBH.md` for the original decision to compute `P_{-1}` only
when the fluid approximation is caio; this note is the follow-up to that.

## Background: what was observed

Timing runs comparing the standard `ncdm` module under CLASS's own fluid
approximation (`ncdm_fluid_approximation = 2`, `ncdmfa_CLASS`, "CFA") against the
same module under the DFA showed the DFA consistently slower. Measured within a
single batch (20 cosmologies x 20 repetitions, both gauges, both evolvers), the
massive-neutrino overhead ratio was

| | evolver 0 (`rk`) | evolver 1 (`ndf15`) |
|---|---|---|
| absolute DFA/CFA wall time | 1.136 - 1.143 | 1.093 - 1.110 |
| overhead ratio `Dt_DFA / Dt_CFA` | 1.38 - 1.43 | 1.26 - 1.28 |

flat across four decades in neutrino mass. This mattered because the DFA is
recommended for use inside the standard `ncdm` module independently of the GBH, so
its cost is a cost everyone pays.

The trigger is **not** the explanation. `ncdmfa_CLASS` and `ncdmfa_caio` switch at
exactly the same time: `perturbations.c`, in `perturbations_approximations`,

```c
if ((tau/tau_k > ppr->ncdm_fluid_trigger_tau_over_tau_k) &&
    (ppr->ncdm_fluid_approximation != ncdmfa_none)) {
  ppw->approx[ppw->index_ap_ncdmfa] = (int)ncdmfa_on;
}
```

with `tau_k = 1/k` and the trigger defaulting to 31. The method selector appears
there only as `!= ncdmfa_none`, i.e. "is any fluid approximation enabled". Which
method is chosen is consulted much later, in `perturbations_derivs`, purely to pick
the closure formulas. So both approximations run the exact hierarchy over the same
time interval, and the extra cost is entirely in what the DFA needs that the CFA
does not.

Two things fell out of that, treated separately below.

---

## Problem 1: `P_{-1}` was recomputed on every shooting iteration

### What the DFA needs

The DFA's effective sound speed needs `lambda = P_{-1}/rho`, where `P_{-1}` is an
extra moment of the background distribution. It is computed in
`background_ncdm_momenta` (`background.c`) as one additional term in the loop that
already accumulates `rho`, `p` and `pseudo_p`:

```c
if (p_min1!=NULL) *p_min1 += pow(epsilon,3.)/3.*wvec[index_q];
```

and the call site gates it on the fluid approximation being caio, so a CFA run pays
nothing. On its own this costs about **2 ms per background evaluation** for a
one-species setup at default precision, a ~26% increase to the background module.
That alone is not a problem.

### Why it became one

The timing runs passed `100*theta_s` rather than `h`, so CLASS shoots for `h`.
`input_try_unknown_parameters` runs a fresh `background_init` (and, for a
`theta_s` target, `thermodynamics_init`) on every trial evaluation, and shooting
takes **6 function evaluations** for this target. With the final real run that is
**7 background computations**, each paying the extra momentum integral.

Measured, isolating the background stage:

| | no shooting | with shooting |
|---|---|---|
| CFA | 8.1 ms | 154.5 ms |
| DFA | 10.2 ms | 170.0 ms |
| **DFA excess** | **2.1 ms** | **15.5 ms** |

15.5 / 2.1 = 7.4, matching the 6 trials plus one real run almost exactly. Nothing
in `background_init` consumes `P_{-1}` — it is written in one place
(`background.c`, storing into `index_bg_P_min1_ncdm1`) and read in exactly one
other (`perturbations.c`, inside the caio branch of the fluid approximation), which
a shooting trial never reaches, since the trial stops at `cs_background` or
`cs_thermodynamics`.

### The fix

A dedicated flag, `pba->compute_P_min1` (`include/background.h`):

- defaults to `_TRUE_`, set in `input_read_parameters` alongside
  `pba->ncdm_fluid_approximation` and in `input_default_params`;
- forced `_FALSE_` for the trial background in `input_try_unknown_parameters`,
  immediately beside the existing `ba.background_verbose = 0` (the trial works on
  its own local `struct background`, so the real run is untouched);
- checked at the `background_ncdm_momenta` call site, which now reads
  `((pba->ncdm_fluid_approximation == ncdmfa_caio) && (pba->compute_P_min1 == _TRUE_))`.

When the flag is `_FALSE_` the slot is filled with `0.` rather than left unwritten,
so the column is never splined from uninitialised memory.

### Why a flag rather than clearing `ncdm_fluid_approximation`

The first attempt simply set `ba.ncdm_fluid_approximation = ncdmfa_CLASS` on the
trial. That works for the gate, but **`index_bg_P_min1_ncdm1` is also defined
conditionally on caio** in `background_indices`:

```c
if (pba->ncdm_fluid_approximation == ncdmfa_caio){
  class_define_index(pba->index_bg_P_min1_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
}
```

so clearing the flag changes the number of columns in the trial's background table,
not just whether one of them is filled. The separate flag leaves the layout alone.

### Verification

Background + thermodynamics only, so shooting dominates the runtime; 21
interleaved repetitions, medians:

| binary | CFA | DFA | excess | DFA/CFA |
|---|---|---|---|---|
| before | 169.9 ms | 184.8 ms | **14.9 ms** | 1.0878 |
| after | 170.3 ms | 172.1 ms | **1.8 ms** | 1.0107 |

**88% of the excess removed**, and the 1.8 ms residual is the one remaining real
background evaluation. The two CFA baselines agree to 0.4 ms, so the noise floor on
this measurement is tight.

---

## Problem 2: background-only quantities rebuilt for every k at every step

### What was happening

In `perturbations_derivs`, the caio branch read:

```c
lambda_ncdm = pvecback[pba->index_bg_P_min1_ncdm1+n_ncdm] / rho_ncdm_bg;
c2_asp = (1. + w_ncdm) / (1. + lambda_ncdm) / 3.;
Omega_m = pba->Omega0_m/pow(a,3.)*pow(pba->H0/(a_prime_over_a/a),2.);
k_fs = sqrt(3. / 2. * Omega_m) * a_prime_over_a / sqrt(c2_asp);
ceff2_ncdm = ca2_ncdm + (c2_asp - ca2_ncdm) * exp(-4. / 3. * k_fs / k);
```

Only the last line depends on `k`. `lambda`, `c2_asp`, `Omega_m` and `k_fs` are
functions of the scale factor alone, yet all four were rebuilt for every mode at
every right-hand-side evaluation, for the whole fluid phase. By comparison the CFA
branch is just `ceff2_ncdm = ca2_ncdm`.

Benchmarked in isolation, the caio expression costs **6.8 ns per call** against
**0.87 ns** for the CFA one — two `pow`, two `sqrt` and one `exp` versus nothing.

### The fix

The useful observation is that the free-streaming scale simplifies exactly.
Substituting `Omega_m(a) = Omega0_m a^-3 (H0/H)^2` and `a_prime_over_a = aH`:

```
k_fs = sqrt(3/2 Omega_m) * aH / sqrt(c_asp^2)
     = sqrt(3/2 Omega0_m) * a^-3/2 * (H0/H) * aH / sqrt(c_asp^2)
     = sqrt(3/2 Omega0_m) * H0 / sqrt(a * c_asp^2)
```

`H` cancels completely, and with it both `pow` calls. The constant
`sqrt(3/2 Omega0_m)*H0` is stored once as `pba->kfs_prefactor`, and the
perturbation code becomes

```c
k_fs = pba->kfs_prefactor / sqrt(a*c2_asp);
```

one divide and one `sqrt` where there were two `pow` and two `sqrt`. `c2_asp`
itself is left computed inline exactly as before, so it is bit-for-bit unchanged.

### Where `kfs_prefactor` must be set, and why not earlier

It is assigned in `background_solve`, on the line immediately after

```c
pba->Omega0_m = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_m];
```

This placement is forced. **`pba->Omega0_m` is a derived quantity read off the last
row of the completed background table**, so it is not available inside
`background_functions` while the table is being built. An earlier version of this
fix precomputed `k_fs` into the background table from within `background_functions`
and silently used an unset `Omega0_m`. Anything that needs `Omega0_m` must be
computed after `background_solve` has filled the table.

An intermediate version also tried storing `c2_asp` itself as a background table
column. That was abandoned: `c2_asp` is proportional to `a^-2` at late times and
spans many orders of magnitude across the table, so interpolating it directly is
less well conditioned than interpolating `rho`, `p` and `P_{-1}` separately and
forming the ratio in `perturbations_derivs`, which is what the code does now.

### Verification

Correctness is exact. See the next section.

Timing, full pipeline with `h` supplied directly so shooting is out of the picture,
21 interleaved repetitions:

| binary | CFA | DFA | excess |
|---|---|---|---|
| before | 828.1 ms | 845.5 ms | 17.4 ms |
| after | 808.4 ms | 818.9 ms | 10.5 ms |

The direction is right, but **this number is not resolved**: the two CFA baselines
should be identical and differ by 19.7 ms, so the noise floor on that machine
exceeded the 7 ms effect being measured. The host was under heavy memory pressure
at the time (~60 MB free RAM, swap nearly exhausted), which also made `classy`
abort on any shooting run, so the whole verification was done with the standalone
binary. **Problem 2's speed-up should be re-measured on a quiet machine before any
figure is quoted.**

---

## Correctness of both fixes together

Compared against a binary built from the unpatched source, same `.ini`, same
cosmology:

- `ncdm_fluid_approximation = 2` (CFA): **bitwise identical** `P(k)`.
- `ncdm_fluid_approximation = 4` (DFA): `P(k)` differs by `5.95e-08`.

That residual is not an error introduced by the fix; it is an error removed. Before
the fix, the value of `h` found by shooting depended on which fluid approximation
the perturbation module was going to use later:

| | shot `h` |
|---|---|
| unpatched, CFA | 0.67059701989124 |
| unpatched, DFA | 0.67059700591162 |
| **patched, DFA** | **0.67059701989124** |

The 2e-8 discrepancy in the unpatched code is a floating-point artifact of
accumulating one extra term in the momentum loop during the trial backgrounds.
`P_{-1}` has no influence on `theta_s`, so the shot `h` should not depend on it,
and after the fix it does not — the patched DFA run shoots to exactly the CFA
value. The `5.95e-08` in `P(k)` is that corrected `h` propagating forward.

The algebraic identity behind problem 2's fix was checked at runtime rather than
only on paper: the old and new `k_fs` expressions were evaluated side by side
during an actual run and agreed to 12 significant figures.

---

## Two traps worth remembering

**`pba->Omega0_m` does not exist until the background table is finished.** It is
read off the last row of `background_table` at the end of `background_solve`.
Using it inside `background_functions` gives whatever it happened to hold before,
with no warning. This cost one full debug cycle.

**CLASS does not overwrite output files, it auto-increments them.** Running with
`root = foo_` repeatedly produces `foo_00_pk.dat`, `foo_01_pk.dat`, and so on. A
comparison script that globs and takes the first match will therefore keep reading
the *oldest* output. This produced a phantom `2.758e-03` discrepancy that survived
two complete rewrites of fix 2 before the stale file was noticed. Delete previous
outputs explicitly in any A/B harness.

---

## Practical notes

- Which input parameter you supply changes the measured DFA penalty by roughly a
  factor of two, because problem 1 only bites when CLASS shoots. A run given `h`
  directly never pays it. Any quoted timing should state which was used.
- The GBH fluid approximation has the same per-step pattern as problem 2 in
  `perturbations_derivs` (its own `c2_asp` / `Omega_m` / `k_fs` block), and the
  same identity applies there. That has **not** been changed, and is the obvious
  next candidate.
- `classy` must be rebuilt (`make && pip install .`) before any of this is visible
  from Python; the verification above used the standalone `./class` binary only.
