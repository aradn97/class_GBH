***************************************************************************************
Input parameters for GBH and their default values:

gbh_use_table                 # Whether to use gbh w_n precomputed table for neutrinos background (1) or use quadrature to
                                perform integrals (0). Default is 1
                                ***Case 0 is not implemented yet.
n_max_gbh                     # Number of velocity moments for the background. Default 20. Note that n_max_gbh+1 
                                quantities will be kept tracked of.
n_max_gbh_table               # Number of columns in the input w table                                

gbh_w_table                   # directory of the file containing pre-computed table of w_n's
















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
3. fix the initial condition for rho at background
4. Be careful with the order of background columns in saved file
5. make sure that all gbh codes are within the condition if has_gbh=true.
6. the condition if (pba->gbh_use_table == 1) in parsing the table address must be after pba->gbh_use_table is assigned. double check.
7. consider removing struct background_gbh_tabulated_w {} as a separate structure
DONE 8. put default value for pba->n_max_gbh
9. make sure all computations are within pba->has_gbh condition
10. put a check for n_max_gbh comparing to the file
11. implement pba->M_gbh and x=ma/T0 in background_functions
12. T_gbh is once used in rho_ini and once in background_functions. check consistency
13. Do we need shooting for rho0_gbh? maybe use the quadrature in ncdm once at the start...
14. Initialization of last_index_gbh might be wrong
15. Does parser_read_list_of_strings automatically allocate memory to pba->gbh_table_address?
16. Consider adding a second pba->gbh_use_table just for rho.
17. for x>100 find rho
18. be careful with units of rho
19. In input.c, put some error condition if m_gbh is not given
20. why do we have ppr->gbh_use_table?
21. is x<1.e-3 handled correctly for rho?
22. is x<1.e-3 handled correctly for w_n's?
23. I'm introducing x=0 by hand. will the interpolation between that and x<1.e-3 give sensible results?
24. for gbh section in initial_conditions function, why don't you use the interpolation of the table that's being done in background_functions?
25. I'm adding has_gbh_pt to be able to turn-off perturbations from input file. later remove it and replace it with has_gbh.
26. define ppv->index_pt_delta_gbh,ppv->index_pt_Delta_gbh,ppv->n_max_gbh,ppw->index_gbh_fa in perturbations.h file and ppr->n_max_gbh in precisions.h. make sure ppv->n_max_gbh>1 so that we have the pressure density
27. For the moment, for each k, I either use GBH at every tau or not at all. The case of switching to FA during evolution should be discussed later.
28. double check with Caio that the x column is the same for the two tables.
29. dimension of Sigma_gbh is (ppv->n_max_gbh)*(ppv->l_max_gbh). dont we want (n+1)l?
30. Why isn't this zero in the initial conditions? l3_ur = ktau_three*2./7./(12.*fracnu+45.)* ppr->curvature_ini;
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
38. work out gauge transformations to sync gauge in print variables: if (pba->has_gbh_pt == _TRUE_) {
        
        /** - --> TODO: gauge transformation of delta, deltaP/rho (?) and theta using -= 3aH(1+w_ncdm) alpha for delta. */
        
      }
39. sync gauge eqs?
DONE 40. be careful that n_max at perturbations doesn't exceed n_max of tables. we need at least pv->l_max_gbh+pv->n_max_gbh+1 columns of w
41. We have no analogue for ncdmfa_none
42. We don't have class_define_index(ppt->index_tp_delta_gbh,   ppt->has_source_delta_gbh,  index_type,1); why
43. make sure in precision.h, evolver is ndf15
44. if you change the definition of Omega_m in neutrino horizon, change both k_fs_gbh def in perturbations and horizon_gbh in background (the latter is integrated)
45. Resolve the ***Pending*** issues in perturbations.c
46. in dy[pba->index_bi_horizon_gbh] = 2. * _PI_ * c_asp / (a * H * sqrt(3. / 2. * pvecback[pba->index_bg_Omega_M])); I'm using index_bg_Omega_M which includes relativistic neutrinos as well as non-relativistic. If we use index_bg_Omega_m we get zero in the denominator.
Also Omega_m = pvecback[pba->index_bg_Omega_M]; in perturbations.
47. GBH chemical potential is assumed to be zero, gbh degeneracy parameter assumed to be one (equivalent to deg_ncdm)
48. use only one 0.71611. now you're using it both in input and background modules
49. use quadrature rho_gbh at initial condition
50. remove rho_M as the matter density that includes neutrinos even if non-rel
51. for(n=0;n<N;n++)
          {
            w[n] = pvecback[pba->index_bg_P_gbh+n]/pvecback[pba->index_bg_P_gbh]/3.;
          }
          should be modified for speed-up, not to calculate it each time.

52. I had made the deadly mistake of defining ppv->n_max_gbh = 16;
      ppv->l_max_gbh = 8;
      only when fluid approx is off, but then I was using these numbers to define N and w[N] in
      perturbations_derivs outside the condition that checks whether we're doing GBH or fluid approx.
      This was giving huge memory issues.

53. free streaming scale at the moment does not include Omega_m. Fix it.
54. Don't get other inputs for gbh if N_gbh=0. put a condition in input.c
