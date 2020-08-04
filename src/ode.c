#include "pihm.h"

int Ode(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
{
    int             i;
    double         *y;
    double         *dy;
    pihm_struct     pihm;
    elem_struct    *elem;
    river_struct   *river;

    y = NV_DATA(CV_Y);
    dy = NV_DATA(CV_Ydot);
    pihm = (pihm_struct)pihm_data;

    elem = &pihm->elem[0];
    river = &pihm->river[0];

    /*
     * Initialization of RHS of ODEs
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < NumStateVar(); i++)
    {
        dy[i] = 0.0;
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem[i].ws.surf  = MAX(y[SURF(i)], 0.0);
        elem[i].ws.unsat = MAX(y[UNSAT(i)], 0.0);
        elem[i].ws.gw    = MAX(y[GW(i)], 0.0);

#if defined(_DGW_)
        elem[i].ws.unsat_geol = MAX(y[FBRUNSAT(i)], 0.0);
        elem[i].ws.gw_geol    = MAX(y[FBRGW(i)], 0.0);
#endif

#if defined(_BGC_) && !defined(_LUMPEDBGC_)
        elem[i].ns.sminn = MAX(y[SOLUTE_SOIL(i, 0)], 0.0);
#endif

#if defined(_CYCLES_)
        elem[i].ps.no3 = MAX(y[SOLUTE_SOIL(i, NO3)], 0.0);
        elem[i].ps.nh4 = MAX(y[SOLUTE_SOIL(i, NH4)], 0.0);
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            elem[i].chms.tot_mol[k] = MAX(y[SOLUTE_SOIL(i, k)], 0.0);
# if defined(_DGW_)
            elem[i].chms.tot_mol[k] = MAX(y[SOLUTE_GEOL(i, k)], 0.0);
# endif
        }
#endif
    }

#if defined(_BGC_) && defined(_LUMPEDBGC_)
    elem[LUMPEDBGC].ns.sminn = MAX(y[LUMPEDBGC_SMINN], 0.0);
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].ws.stage = MAX(y[RIVER(i)], 0.0);

#if defined(_BGC_) && !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
        river[i].ns.streamn = MAX(y[SOLUTE_RIVER(i, 0)], 0.0);
#endif

#if defined(_CYCLES_)
        river[i].ns.no3 = MAX(y[SOLUTE_RIVER(i, NO3)], 0.0);
        river[i].ns.nh4 = MAX(y[SOLUTE_RIVER(i, NH4)], 0.0);
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            river[i].chms.tot_mol[k] = MAX(y[SOLUTE_RIVER(i, k)], 0.0);
        }
#endif
    }

    /*
     * PIHM Hydrology fluxes
     */
    Hydrol(&pihm->ctrl, pihm->elem, pihm->river);

#if _OBSOLETE_
#if defined(_BGC_)
    /*
     * Nitrogen transport fluxes
     */
# if defined(_LUMPEDBGC_)
    NLeachingLumped(pihm->elem, pihm->river);
# elif defined(_LEACHING_)
    NLeaching(pihm->elem);
# else
    NTransport(pihm->elem, pihm->river);
# endif
#endif
#endif

    /*
     * Calculate solute concentrations
     */
#if defined(_BGC_)
    SoluteConc(pihm->elem, pihm->river);
#elif defined(_CYCLES_)
    SoluteConc(t + (realtype)pihm->ctrl.tout[0] -
        (realtype)pihm->ctrl.tout[pihm->ctrl.cstep], pihm->elem, pihm->river);
#elif defined(_RT_)
    SoluteConc(pihm->chemtbl, &pihm->rttbl, pihm->elem, pihm->river);
#endif


#if defined(_BGC_) || defined(_CYCLES_)
    SoluteTranspt(0.0, 0.0, 0.0, pihm->elem, pihm->river);
#elif defined(_RT_)
    SoluteTranspt(pihm->rttbl.diff_coef, pihm->rttbl.disp_coef,
        pihm->rttbl.cementation, pihm->elem, pihm->river);
#endif

    /*
     * Build RHS of ODEs
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        /*
         * Vertical water fluxes for surface and subsurface
         */
        dy[SURF(i)] += elem[i].wf.pcpdrp - elem[i].wf.infil -
            elem[i].wf.edir_surf;
        dy[UNSAT(i)] += elem[i].wf.infil - elem[i].wf.recharge -
            elem[i].wf.edir_unsat - elem[i].wf.ett_unsat;
        dy[GW(i)] += elem[i].wf.recharge - elem[i].wf.edir_gw - elem[i].wf.ett_gw;

#if defined(_DGW_)
        /*
         * Vertical water fluxes for fractured bedrock
         */
        dy[GW(i)] -= elem[i].wf.infil_geol;

        dy[FBRUNSAT(i)] += elem[i].wf.infil_geol - elem[i].wf.rechg_geol;
        dy[FBRGW(i)]    += elem[i].wf.rechg_geol;
# if defined(_LUMPED_)
        dy[FBRGW(i)] -= elem[i].wf.dgw_runoff;
# endif
#endif

        /*
         * Horizontal water fluxes
         */
        for (j = 0; j < NUM_EDGE; j++)
        {
            dy[SURF(i)]  -= elem[i].wf.overland[j] / elem[i].topo.area;
            dy[GW(i)]    -= elem[i].wf.subsurf[j] / elem[i].topo.area;
#if defined(_DGW_)
            dy[FBRGW(i)] -= elem[i].wf.dgw[j] / elem[i].topo.area;
#endif
        }

        dy[UNSAT(i)] /= elem[i].soil.porosity;
        dy[GW(i)]    /= elem[i].soil.porosity;
#if defined(_DGW_)
        dy[FBRUNSAT(i)] /= elem[i].geol.porosity;
        dy[FBRGW(i)]    /= elem[i].geol.porosity;
#endif

#if _OBSOLETE_
#if defined(_BGC_) && !defined(_LUMPEDBGC_)
# if !defined(_LEACHING_)
        /*
         * BGC N transport fluxes
         */
        dy[SURFN(i)] +=
            (elem[i].nf.ndep_to_sminn + elem[i].nf.nfix_to_sminn) / DAYINSEC -
            elem[i].nsol.infilflux;
        dy[SMINN(i)] += elem[i].nsol.infilflux + elem[i].nsol.snksrc;
# else
        dy[SMINN(i)] +=
            (elem[i].nf.ndep_to_sminn + elem[i].nf.nfix_to_sminn) / DAYINSEC +
            elem[i].nsol.snksrc;
# endif

        for (j = 0; j < NUM_EDGE; j++)
        {
# if !defined(_LEACHING_)
            dy[SURFN(i)] -= elem[i].nsol.ovlflux[j] / elem[i].topo.area;
# endif
            dy[SMINN(i)] -= elem[i].nsol.subflux[j] / elem[i].topo.area;
        }
#endif
#endif

#if defined(_CYCLES_) || defined(_BGC_) || defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
# if defined(_CYCLES_)
            dy[SOLUTE_SOIL(i, k)] += elem[i].solute[k].infil +
                Profile(elem[i].ps.nlayers, elem[i].solute[k].snksrc);
# else
            dy[SOLUTE_SOIL(i, k)] += elem[i].solute[k].infil +
                elem[i].solute[k].snksrc;
# endif

# if defined(_DGW_)
            dy[SOLUTE_SOIL(i, k)] -= elem[i].solute[k].infil_geol;

            dy[SOLUTE_GEOL(i, k)] += elem[i].solute[k].infil_geol +
                elem[i].solute[k].snksrc_geol;
#  if defined(_LUMPED_)
            dy[SOLUTE_GEOL(i, k)] -= elem[i].solute[k].dgw_leach;
#  endif
# endif

            for (j = 0; j < NUM_EDGE; j++)
            {
                dy[SOLUTE_SOIL(i, k)] -= elem[i].solute[k].subflux[j] /
                    elem[i].topo.area;
# if defined(_DGW_)
                dy[SOLUTE_GEOL(i, k)] -= elem[i].solute[k].dgwflux[j] /
                    elem[i].topo.area;
# endif
            }
        }
#endif
    }

#if defined(_BGC_) && defined(_LUMPEDBGC_)
    elem_struct    *elem;

    elem = &pihm->elem[LUMPEDBGC];

    dy[LUMPEDBGC_SMINN] +=
        (elem->nf.ndep_to_sminn + elem->nf.nfix_to_sminn) / DAYINSEC +
        elem->nsol.snksrc;
#endif

    /*
     * ODEs for river segments
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;

        for (j = 0; j < NUM_RIVFLX; j++)
        {
            /* Note the limitation due to
             * d(v) / dt = a * dy / dt + y * da / dt
             * for cs other than rectangle */
            dy[RIVER(i)] -= river[i].wf.rivflow[j] / river[i].topo.area;
        }

#if _OBSOLETE_
#if defined(_BGC_) && !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
        for (j = 0; j <= 6; j++)
        {
            dy[STREAMN(i)] -= river[i].nsol.flux[j] / river[i].topo.area;
        }
#endif
#endif

#if defined(_CYCLES_) || defined(_BGC_) || defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            for (j = 0; j < NUM_RIVFLX; j++)
            {
                dy[SOLUTE_RIVER(i, k)] -= river[i].solute[k].flux[j] /
                    river[i].topo.area;
            }
        }
#endif
    }

    return 0;
}

int NumStateVar(void)
{
    /*
     * Return number of state variables
     */
    int             nsv;

    nsv = 3 * nelem + nriver;

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
    nsv += nsolute * (nelem + nriver);
#endif

#if TEMP_DISABLED
#if defined(_BGC_)
# if defined(_LUMPEDBGC_)
    nsv += 1;
# else
    nsv += 2 * nelem + 2 * nriver;
# endif
#endif
#endif

#if defined(_DGW_)
    nsv += 2 * nelem;
# if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
    nsv += nsolute * nelem;
# endif
#endif

    return nsv;
}

void SetCVodeParam(pihm_struct pihm, void *cvode_mem, SUNLinearSolver *sun_ls,
    N_Vector CV_Y)
{
    int             cv_flag;
    static int      reset;
    N_Vector        abstol;
#if defined(_BGC_) || defined(_CYCLES_)
    const double    TRANSP_TOL = 1.0E-5;
#elif defined(_RT_)
    const double    TRANSP_TOL = 1.0E-8;
#endif

    pihm->ctrl.maxstep = pihm->ctrl.stepsize;

    if (reset)
    {
        /* When model spins-up and recycles forcing, use CVodeReInit to reset
         * solver time, which does not allocates memory */
        cv_flag = CVodeReInit(cvode_mem, 0.0, CV_Y);
        CheckCVodeFlag(cv_flag);
    }
    else
    {
        cv_flag = CVodeInit(cvode_mem, Ode, 0.0, CV_Y);
        CheckCVodeFlag(cv_flag);
        reset = 1;

        *sun_ls = SUNLinSol_SPGMR(CV_Y, PREC_NONE, 0);

        /* Attach the linear solver */
        CVodeSetLinearSolver(cvode_mem, *sun_ls, NULL);

        /* When BGC, Cycles, or RT module is turned on, both water storage and
         * transport variables are in the CVODE vector. A vector of absolute
         * tolerances is needed to specify different absolute tolerances for
         * water storage variables and transport variables */
        abstol = N_VNew(NumStateVar());
#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
        SetAbsTolArray(pihm->ctrl.abstol, TRANSP_TOL, abstol);
#else
        SetAbsTolArray(pihm->ctrl.abstol, abstol);
#endif

        cv_flag = CVodeSVtolerances(cvode_mem, (realtype)pihm->ctrl.reltol,
            abstol);
        CheckCVodeFlag(cv_flag);

        N_VDestroy(abstol);

        /* Specifies PIHM data block and attaches it to the main cvode memory
         * block */
        cv_flag = CVodeSetUserData(cvode_mem, pihm);
        CheckCVodeFlag(cv_flag);

        /* Specifies the initial step size */
        cv_flag = CVodeSetInitStep(cvode_mem, (realtype)pihm->ctrl.initstep);
        CheckCVodeFlag(cv_flag);

        /* Indicates if the BDF stability limit detection algorithm should be
         * used */
        cv_flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
        CheckCVodeFlag(cv_flag);

        /* Specifies an upper bound on the magnitude of the step size */
        cv_flag = CVodeSetMaxStep(cvode_mem, (realtype)pihm->ctrl.maxstep);
        CheckCVodeFlag(cv_flag);

        /* Specifies the maximum number of steps to be taken by the solver in
         * its attempt to reach the next output time */
        cv_flag = CVodeSetMaxNumSteps(cvode_mem, pihm->ctrl.stepsize * 10);
        CheckCVodeFlag(cv_flag);
    }
}

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
void SetAbsTolArray(double hydrol_tol, double transp_tol, N_Vector abstol)
#else
void SetAbsTolArray(double hydrol_tol, N_Vector abstol)
#endif
{
    int             i;
    int             num_hydrol_var;

    num_hydrol_var = 3 * nelem + nriver;

#if defined(_DGW_)
    num_hydrol_var += 2 * nelem;
#endif

    /* Set absolute errors for hydrologic state variables */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < num_hydrol_var; i++)
    {
        NV_Ith(abstol, i) = (realtype)hydrol_tol;
    }

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
    /* Set absolute errors for solute state variables */
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = num_hydrol_var; i < NumStateVar(); i++)
    {
        NV_Ith(abstol, i) = (realtype)transp_tol;
    }
#endif
}

void SolveCVode(double cputime, const ctrl_struct *ctrl, int *t,
    void *cvode_mem, N_Vector CV_Y)
{
    realtype        solvert;
    realtype        tout;
    pihm_t_struct   pihm_time;
    int             starttime;
    int             nextptr;
    int             cv_flag;
    double          progress;

    starttime = ctrl->starttime;
    nextptr = ctrl->tout[ctrl->cstep + 1];

    tout = (realtype)(nextptr - starttime);

    /* Specifies the value of the independent variable t past which the solution
     * is not to proceed */
    cv_flag = CVodeSetStopTime(cvode_mem, tout);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVode(cvode_mem, tout, CV_Y, &solvert, CV_NORMAL);
    CheckCVodeFlag(cv_flag);

    *t = roundi(solvert) + starttime;

    pihm_time = PIHMTime(*t);

    progress = ((double)ctrl->cstep + 1.0) / (double)ctrl->nstep;

    if (ctrl->cstep == 0)
    {
        pihm_printf(VL_NORMAL, "\n");
    }

    if (debug_mode)
    {
        pihm_printf(VL_NORMAL, "\033[1A\rStep = %s (t = %d)\n",
            pihm_time.str, *t);
        ProgressBar(progress);
    }
    else if (spinup_mode)
    {
        if (pihm_time.t % DAYINSEC == 0)
        {
            pihm_printf(VL_NORMAL, "\033[1A\rStep = %s\n", pihm_time.str);
            ProgressBar(progress);
        }
    }
    else if (pihm_time.t % 3600 == 0)
    {
        pihm_printf(VL_NORMAL, "\033[1A\rStep = %s (cputime %8.2f s)\n",
            pihm_time.str, cputime);
        ProgressBar(progress);
    }
}

void AdjCVodeMaxStep(void *cvode_mem, ctrl_struct *ctrl)
{
    /* Variable CVODE max step (to reduce oscillations) */
    long int        nst;
    long int        ncfn;
    long int        nni;
    static long int nst0;
    static long int ncfn0;
    static long int nni0;
    int             cv_flag;
    double          nsteps;
    double          nfails;
    double          niters;

    /* Gets the cumulative number of internal steps taken by the solver (total
     * so far) */
    cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
    CheckCVodeFlag(cv_flag);

    /* Gets the number of nonlinear convergence failures that have occurred */
    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    CheckCVodeFlag(cv_flag);

    /* Gets the number of nonlinear iterations performed */
    cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    CheckCVodeFlag(cv_flag);

    nsteps = (double)(nst - nst0);
    nfails = (double)(ncfn - ncfn0) / nsteps;
    niters = (double)(nni - nni0) / nsteps;

    ctrl->maxstep /= (nfails > ctrl->nncfn || niters >= ctrl->nnimax) ?
         ctrl->decr : 1.0;

    ctrl->maxstep *= (nfails == 0.0 && niters <= ctrl->nnimin) ?
        ctrl->incr : 1.0;

    ctrl->maxstep = MIN(ctrl->maxstep, ctrl->stepsize);
    ctrl->maxstep = MAX(ctrl->maxstep, ctrl->stmin);

    /* Updates the upper bound on the magnitude of the step size */
    cv_flag = CVodeSetMaxStep(cvode_mem, (realtype)ctrl->maxstep);
    CheckCVodeFlag(cv_flag);

    nst0  = nst;
    ncfn0 = ncfn;
    nni0  = nni;
}
