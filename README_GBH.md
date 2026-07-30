***************************************************************************************
Input parameters for GBH and their default values:

N_gbh                          # Number of distinct gbh species. Default 0 (gbh off), mirrors N_ncdm exactly:
                                the whole GBH sector is disabled unless N_gbh>0 is set explicitly.
m_gbh                          # Mass of each gbh species, in eV. One value if N_gbh=1, or N_gbh comma-separated
                                values if N_gbh>1 (mirrors m_ncdm). Required whenever N_gbh>0.
deg_gbh                        # Degeneracy of each gbh species. One value, or N_gbh comma-separated values if
                                N_gbh>1 (mirrors deg_ncdm). Default 1 per species.
T_gbh                          # present temperature of gbh / T_cmb. Default 0.71611. Shared across all gbh
                                species (not yet per-species).

gbh_use_table                 # Whether to use gbh w_n precomputed table for neutrinos background (1) or use quadrature to
                                perform integrals (0). Default is 1
                                ***Case 0 is not implemented yet.
n_max_gbh                     # Number of velocity moments for the background. Default 20. Note that n_max_gbh+1 
                                quantities will be kept tracked of.
n_max_gbh_table               # Number of columns in the input w table                                

gbh_w_table                   # directory of the file containing pre-computed table of w_n's

***************************************************************************************
Note on multi-species status (as of the N_gbh/deg_gbh/m_gbh restructuring):
N_gbh/deg_gbh/m_gbh are now real per-species input fields (pba->N_gbh, pba->deg_gbh[],
pba->m_gbh_in_eV[]/pba->M_gbh[]), same idiom as ncdm. However, only species index [0] is
actually consumed by background.c/perturbations.c so far -- setting N_gbh>1 will build and
run but will NOT evolve additional species yet. See multi_species_gbh_implementation_plan.md
for the remaining work (per-species background-vector indices, perturbation state-vector
indices, Boltzmann hierarchy loop, etc.) before N_gbh>1 is physically meaningful.
















*************************************************************************************** 
Flags used for implementing the Generalized Boltzmann Hierarchy are explained here







1) Between /*GBH_bg_start*/ and /*GBH_bg_end*/: All the places where the background hierarchy in powers of velocity is used. 
   For one-line codes: //GBH_bg
    The files include:
        precisions.h,   background.c,   input.c,    background.h,      perturbations.c

2) Between /*GBH_pt_start*/ and /*GBH_pt_end*/.
    For one-liners: //GBH_pt
    The files include:
        background.c    background.h

3) /*ncdm_caio_fa*/ determines the places where Caio's FA affects ncdm code.














***************************************************************************************
TODO:
Notes of caution during coding.

1. pba->has_gbh has been added. make sure to turn it on in all required circumstances
2. the above requires Omega0_gbh to be non-zero. introduce it properly.
DONE 3. fix the initial condition for rho at background
4. Be careful with the order of background columns in saved file
5. make sure that all gbh codes are within the condition if has_gbh=true.
DONE 6. the condition if (pba->gbh_use_table == 1) in parsing the table address must be after pba->gbh_use_table is assigned. double check.
7. consider removing struct background_gbh_tabulated_w {} as a separate structure
DONE 8. put default value for pba->n_max_gbh
9. make sure all computations are within pba->has_gbh condition
10. put a check for n_max_gbh comparing to the file
11. implement pba->M_gbh and x=ma/T0 in background_functions
DONE 12. T_gbh is once used in rho_ini and once in background_functions. check consistency
DONE 13. Do we need shooting for rho0_gbh? maybe use the quadrature in ncdm once at the start...
14. Initialization of last_index_gbh might be wrong
DONE 15. Does parser_read_list_of_strings automatically allocate memory to pba->gbh_table_address?
DONE (not necessary) 16. Consider adding a second pba->gbh_use_table just for rho.
DONE 17. for x>100 find rho
DONE 18. be careful with units of rho
DONE (input_read_parameters_species now has a class_test rejecting m_gbh_in_eV[n]<=0 for every species; an omitted `m_gbh` falls back to the 0.0 default and is caught by the same test) 19. In input.c, put some error condition if m_gbh is not given
20. why do we have ppr->gbh_use_table?
21. is x<1.e-3 handled correctly for rho?
22. is x<1.e-3 handled correctly for w_n's?
23. I'm introducing x=0 by hand. will the interpolation between that and x<1.e-3 give sensible results?
24. for gbh section in initial_conditions function, why don't you use the interpolation of the table that's being done in background_functions?
DONE 25. I'm adding has_gbh_pt to be able to turn-off perturbations from input file. later remove it and replace it with has_gbh.
26. define ppv->index_pt_delta_gbh,ppv->index_pt_Delta_gbh,ppv->n_max_gbh,ppw->index_gbh_fa in perturbations.h file and ppr->n_max_gbh in precisions.h. make sure ppv->n_max_gbh>1 so that we have the pressure density
DONE (obsolete note) 27. For the moment, for each k, I either use GBH at every tau or not at all. The case of switching to FA during evolution should be discussed later.
DONE 28. double check with Caio that the x column is the same for the two tables.
DONE (equation-by-equation check against the two GBH papers confirms the n_max*l_max sizing, i.e. n=0...n_max-1, is correct as-is; no (n+1) needed) 29. dimension of Sigma_gbh is (ppv->n_max_gbh)*(ppv->l_max_gbh). dont we want (n+1)l?
DONE (because it directly enters in the equation for sigma') 30. Why isn't this zero in the initial conditions? l3_ur = ktau_three*2./7./(12.*fracnu+45.)* ppr->curvature_ini;
31. Maybe implement the FA at very low k, x<1?
32. I'm using p_gbh = ppw->pvecback[pba->index_bg_P_gbh+1]; in perturbations. what if n_max is chosen to be 1?
33. NOTE: we are using cg2_gbh = w_gbh*(1.0-1.0/(3.0+3.0*w_gbh)*(3.0*w_gbh-2.0+pseudo_p_gbh/p_gbh));
        ppw->delta_p += cg2_gbh*rho_gbh*delta_gbh;  as the pressure contribution in fluid approximation (x>30).
34. Why is class adding all of the ncdm energy density to delta_rho_m? what if it's relativistic?
    if (ppt->has_source_delta_m == _TRUE_) {
            delta_rho_m += rho_gbh*delta_gbh; // contribution to delta rho_matter
            rho_m += rho_gbh;
        }
35. TBC: gravitational wave contribution of gbh species
DONE 36. add has_source_delta_gbh, ppt->index_tp_delta_gbh (and for theta) everywhere needed, not just in perturbations_sources
37. think about what you want to print in pt to outputs
38. /** TODO: use c_eff^2 (which is different from c_a^2 in DFA) and don't use 0.*/

DONE 39. sync gauge eqs?
DONE 40. be careful that n_max at perturbations doesn't exceed n_max of tables. we need at least pv->l_max_gbh+pv->n_max_gbh+1 columns of w
DONE (cause we have no high-precision setting) 41. We have no analogue for ncdmfa_none
42. We don't have class_define_index(ppt->index_tp_delta_gbh,   ppt->has_source_delta_gbh,  index_type,1); why
43. make sure in precision.h, evolver is ndf15
44. if you change the definition of Omega_m in neutrino horizon, change both k_fs_gbh def in perturbations and horizon_gbh in background (the latter is integrated)
DONE 45. Resolve the ***Pending*** issues in perturbations.c
DONE (discussed with Caio: Omega_m~Omega0_m/a^3(H0/H)^2) 46. in dy[pba->index_bi_horizon_gbh] = 2. * _PI_ * c_asp / (a * H * sqrt(3. / 2. * pvecback[pba->index_bg_Omega_M])); I'm using index_bg_Omega_M which includes relativistic neutrinos as well as non-relativistic. If we use index_bg_Omega_m we get zero in the denominator.
Also Omega_m = pvecback[pba->index_bg_Omega_M]; in perturbations.
47. GBH chemical potential is assumed to be zero, gbh degeneracy parameter assumed to be one (equivalent to deg_ncdm)
DONE 48. use only one 0.71611. now you're using it both in input and background modules
DONE 49. use quadrature rho_gbh at initial condition
DONE 50. remove rho_M as the matter density that includes neutrinos even if non-rel
51. for(n=0;n<N;n++)
          {
            w[n] = pvecback[pba->index_bg_P_gbh+n]/pvecback[pba->index_bg_P_gbh]/3.;
          }
          should be modified for speed-up, not to calculate it each time.

DONE 52. I had made the deadly mistake of defining ppv->n_max_gbh = 16;
      ppv->l_max_gbh = 8;
      only when fluid approx is off, but then I was using these numbers to define N and w[N] in
      perturbations_derivs outside the condition that checks whether we're doing GBH or fluid approx.
      This was giving huge memory issues.

DONE 53. free streaming scale at the moment does not include Omega_m. Fix it.
54. Don't get other inputs for gbh if N_gbh=0. put a condition in input.c
DONE 55. remove index_bg_lambda_gbh and index_bg_Pminus1_gbh. They are not used in perturbations.c.
56. extra stuff I've ordered to print: //GBH_pt_print. You can remove them at the end.
DONE 57. if the code is slow, look at the perturbations initial condition for gbh, where I integrate in momentum space. find a way to speed it up.
58. For the case pba->gbh_use_table == 1, also provide a table for P_{-1}, so we wouldnt need quadrature in the background
DONE 59. change 1000. in interpolation (both in input.c and background.c) to some precision parameter, that also needs to be input together with the file
DONE 60. for sync gauge, change FA
DONE 61. for sync gauge, change trunc scheme
DONE 62. for sync gauge, change init conditions
63. check whether in presence of fld, caio fluid approximation with sync gauge still works fine (finds correct sigma_ncdm or gbh)
DONE 64. no need to even compute quadrature weights if gbh_init_condition_integrate is FALSE.
DONE 65. save w_n's instead of P_n's?
66. make sure you don't redundantly read precision parameters in input.c
67. remove n_max_gbh as an input parameter if use_table is 1. It should be computed based on what we need for perturbations.
Warn in the paper or to the user how large the table needs to be
68. remove ppr->gbh_nl_max_method==1
DONE (background_gbh_init now derives pba->n_max_gbh_table from the actual table file's column count, rather than trusting the input.c placeholder default of 31, and warns + falls back to gbh_use_table=0/quadrature if it's too small for the requested n_max_gbh; verified this triggers correctly with an oversized gbh_FA_trigger. The stale default of 31 is now only ever reached, harmlessly, when gbh_use_table=0 and the table is never read) 69. The error n_max_gbh_table<n_max_gbh is not being triggered when n_max_gbh_table is not given in the input and the default 31 is being used. either remove n_max_gbh totally or put error if the default of n_max_gbh_table is smaller than n_max_gbh
DONE 70. set the default in input.c to pba->gbh_use_table = 1, and put the table somewhere so that user does not have to provide it.
DONE 71. move w_n table to external
DONE 72. do not hard code 1000 as the max x value of files
73. Remove the computation of gbh pert init conditions by quadrature.
DONE 74. ncdm caio FA: update, removing the explicit shear formula--> then remove index_bg_horizon_ncdm1 from background as well
DONE 75. remove extra precision params use_formula and gbh_FA_caio_transition_rate
76. /** - --> TODO: gauge transformation of delta, deltaP/rho (?) and theta using -= 3aH(1+w_ncdm) alpha for delta. */
this is now done for ncdm. do it for gbh
77. compute w_-1 at background only if gbh_fa is caio.
DONE 78. dont call class_call(background_ncdm_momenta(pba->q_gbh_bg, ...)) in the initial condition checker of background for gbh.
79. OPEN (not fixed this round -- documenting only): evolver=1 (ndf15) produces silently unphysical
    GBH results (sigma8 off by many orders of magnitude, d_gbh(k) blown up to ~1e16 for some k) for
    m_gbh above roughly 1.5-2 eV, well before the mass (~10 eV) where it actually raises an error.
    evolver=0 (rk) does not show this (see item 80 below for a related but separate rk gap that was
    fixed). gets worse under a *tighter* tol_perturbations_integration, consistent with a genuine
    ndf15/GBH-hierarchy ill-conditioning rather than a borderline-precision artifact. Until fixed,
    do not trust gbh_fluid_approximation results with evolver=1 above ~1.5-2 eV; use evolver=0.

    Direct investigation this round (m_gbh=5 eV, single species, k_output_values isolating one k):
    the actual `evolver_ndf15` crash ("Step size too small") reproduces on the *smallest* k in the
    internal source-sampling grid (~8e-6 Mpc^-1 for this setup) -- a mode that stays in the exact
    (non-fluid) hierarchy for the entire run (n_max=4, l_max=3, never crosses gbh_FA_trigger), not
    a large-k/small-scale mode as previously assumed. A temporary debug instrumentation was added to
    `perturbations_derivs` (guarded behind `GBH_DEBUG_FILE`/`GBH_DEBUG_FA_FILE`/`GBH_DEBUG_K` env
    vars, currently uncommitted in source/perturbations.c) to dump every Delta_gbh[n]/Sigma_gbh[n,l]
    value at every RHS evaluation for both evolvers at this k. Result: the dumped state vector is
    completely smooth and finite all the way to today for *both* evolvers -- no NaN, no blow-up, no
    discontinuity -- even though ndf15 fails on this exact k/tau range. This directly rules out,
    by computation rather than assertion:
      - the n_max-boundary closure ratios (`delta_next`, `sigma_ns`/`sigma_sn`) -- computed directly
        from real background output, both stay O(1e-6)-to-O(1), smooth, no blow-up, even deep in the
        w_n(x) table-extrapolation regime.
      - the fluid-approximation's background-derived coefficients (`ca2_gbh`, `c2_asp`, `k_fs`,
        `w_min1_gbh`/lambda_gbh) -- also directly computed from background.dat, all smooth, no sign
        flip, no division-by-zero.
      - the l=1 Newtonian-gauge metric-coupling term (`delta_l1`, proportional to k, so it *vanishes*
        as k->0 rather than diverging).
    Also ruled out a generic ndf15-at-high-mass issue unrelated to GBH: the identical cosmology /
    mass / evolver=1 / k-grid run with plain `ncdm` (no GBH at all) integrates fine
    (sigma8=0.557, sane) -- so the failure is specific to the GBH exact hierarchy's implementation,
    not massive species + ndf15 in general.
    Follow-up (this round): added a second temporary debug hook, this time inside
    `evolver_ndf15` itself (also guarded by env var, `GBH_JAC_DUMP_FILE`, currently uncommitted in
    tools/evolver_ndf15.c), that dumps the dense finite-difference Jacobian (`jac.dfdy`, always
    populated regardless of the sparse/dense internal representation) at the exact moment the
    `absh <= hmin` failure fires. For the same failing case (m=5eV, k~8e-6/Mpc, n_max=4, l_max=3):
    the Jacobian at failure (neq=203) is genuinely near-singular, condition number ~1e16 (machine
    precision), with 3 columns (state indices 1, 3, 5 -- photon delta_g, shear_g, and the l=4
    photon multipole, in the post-tight-coupling variable ordering) numerically zero or near-zero.
    This is suggestive but not fully conclusive on its own, since some Jacobian ill-conditioning at
    such a tiny k may be generic to any species combination (the relevant couplings scale as k^2,
    tiny for k~8e-6) -- this was not cross-checked against a comparably-conditioned ncdm Jacobian.

    A more decisive test: forcing `gbh_FA_trigger` down to a tiny value (0.0001) so that *every* k,
    including this problematic one, immediately uses the fluid approximation instead of the exact
    hierarchy -- with everything else unchanged (still evolver=1, still m=5 eV) -- makes the run
    succeed cleanly (sigma8=0.536, sane). This is the most direct evidence yet: the bug is
    specifically in the **exact (non-fluid) hierarchy's RHS at low truncation order**
    (small n_max/l_max, e.g. n_max=4/l_max=3), not in the fluid approximation, not a generic
    small-k-with-ndf15 issue, and not shared with plain ncdm. The likely next step is to audit the
    exact-hierarchy equations in `perturbations_derivs` (source/perturbations.c, the
    `else{//GBH eqs}` block) specifically for correctness -- not just numerical conditioning -- at
    very small n_max/l_max (e.g. an off-by-one or a fallback formula, such as the `ll==1` branch's
    `sigma_np = 3/(1+w[1])*delta_next` shortcut, or the interaction between the simultaneous
    n_max-boundary and l_max-boundary closures when n_max and l_max are both this small) rather
    than assuming the equations are correct and only mis-integrated.

    Code audit performed (this round): traced every array index used in the exact-hierarchy block
    of perturbations_derivs (the else{//GBH eqs} block) by hand for n_max=4, l_max=3, including
    the "corner" case where the n_max-boundary and l_max-boundary closures fire simultaneously
    (n=n_max-1 AND ll=l_max at once) and the ll==1 fallback (sigma_np = 3/(1+w[1])*delta_next). No
    out-of-bounds access or obviously wrong index/fallback was found. Note the existing code
    comment: "DO NOT use l_max=2, because the truncation scheme...is gauge invariant only for
    l_max>2" -- there is already a floor l_max>=3 for exactly this reason, and our failing case
    sits right at that floor (l_max=3), i.e. at the theoretical minimum the author already flagged
    as marginal.

    Decisive follow-up test distinguishing "too-coarse truncation" from "genuine dynamics": reran
    the same m=5eV/evolver=1 case with gbh_nl_max_method=1 (uniform n_max=19/l_max=10 for every k,
    well above the small-k floor). This does *not* crash -- but produces the previously-reported
    silent corruption instead (sigma8~1e11), from a *different*, mid-range k (~0.097/Mpc, not the
    ~8e-6/Mpc mode above). So raising the truncation order does not fix the underlying issue, it
    just changes which failure mode you get and which k triggers it -- ruling out "n_max/l_max too
    small" as the root cause.
    Dumping the state vector for that k (same debug hook) shows something new: Delta_gbh[n=0]
    under ndf15 grows smoothly through a sign flip around tau~6390 Mpc (a~0.236, from ~-58 to ~+62)
    and then runs away to >1e16 by the end of the integration. Under rk, at the exact same k and
    the exact same tau range, Delta_gbh[n=0] stays smooth, monotonic, same sign throughout, an
    order of magnitude smaller, with no runaway -- i.e. the two evolvers disagree in *sign*, not
    just precision. Counting RHS evaluations in that same tau window: rk takes ~349, ndf15 only
    ~46 (including repeated Jacobian-only calls that don't advance tau) -- roughly 7-8x fewer
    genuine time steps. This is the signature of a stiff or rapidly-varying feature in the exact
    hierarchy's dynamics that ndf15's adaptive step-size control under-resolves (accepts too large
    a step based on its local error estimate, aliases through real structure, and the resulting
    error compounds into a spurious runaway) -- a classic numerical under-resolution artifact, not
    a formula/indexing bug (which the audit above did not find). rk's step size, tied directly to
    physical timescales via perturbations_timescale rather than a measured local error, happens to
    stay fine enough to track this feature by construction.

    Current status: the failure is a genuine numerical-resolution problem in ndf15's handling of
    the GBH exact hierarchy's fast/stiff dynamics (manifesting as either a hard stop at low
    truncation order, or silent runaway at high truncation order, depending on which k/truncation
    combination is involved) -- not a coding bug in the closure formulas or index arithmetic, which
    were audited directly and found correct for the cases checked. A structural fix (e.g. a tighter
    default tol_perturbations_integration specifically for GBH-active runs, or capping ndf15's
    maximum step size directly rather than relying on its own error estimate) has not been
    implemented. Until then, evolver=0 (rk) remains the only evolver validated for GBH at these
    masses.
DONE 80. perturbations_timescale (used only by evolver=0/rk to pick its step size) checked
    pba->has_ncdm but not pba->has_gbh when deciding whether to resolve tau_k=1/k in the step-size
    estimate -- meaning a GBH-only run's rk steps could ignore k entirely once radiation streaming
    turned on, even while a GBH species was still in its exact (non-fluid) hierarchy needing
    k-resolved steps. Fixed by adding the missing pba->has_gbh check.