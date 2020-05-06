#include "pihm.h"

int Ode(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
{
    int             i;
    double         *y;
    double         *dy;
    pihm_struct     pihm;

    y = NV_DATA(CV_Y);
    dy = NV_DATA(CV_Ydot);
    pihm = (pihm_struct)pihm_data;

    /*
     * Initialization of RHS of ODEs
     */
    for (i = 0; i < NumStateVar(); i++)
    {
        dy[i] = 0.0;
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct    *elem;

        elem = &pihm->elem[i];

        elem->ws.surf = (y[SURF(i)] >= 0.0) ? y[SURF(i)] : 0.0;
        elem->ws.unsat = (y[UNSAT(i)] >= 0.0) ? y[UNSAT(i)] : 0.0;
        elem->ws.gw = (y[GW(i)] >= 0.0) ? y[GW(i)] : 0.0;

#if defined(_FBR_)
        elem->ws.fbr_unsat = MAX(y[FBRUNSAT(i)], 0.0);
        elem->ws.fbr_gw = MAX(y[FBRGW(i)], 0.0);
#endif

#if defined(_BGC_) && !defined(_LUMPED_)
        elem->ns.surfn = (y[SURFN(i)] >= 0.0) ? y[SURFN(i)] : 0.0;
        elem->ns.sminn = (y[SMINN(i)] >= 0.0) ? y[SMINN(i)] : 0.0;
#endif

#if defined(_CYCLES_OBSOLETE_)
        elem->np.no3 = (y[NO3(i)] >= 0.0) ? y[NO3(i)] : 0.0;
        elem->np.nh4 = (y[NH4(i)] >= 0.0) ? y[NH4(i)] : 0.0;
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            elem->chms.t_mole[k] = MAX(y[SOLUTE_SOIL(i, k)], 0.0);
# if defined(_FBR_)
            elem->chms.t_mole[k] = MAX(y[SOLUTE_GEOL(i, k)], 0.0);
# endif
        }
#endif
    }

#if defined(_BGC_) && defined(_LUMPED_)
    pihm->elem[LUMPED].ns.sminn = (y[LUMPED_SMINN] >= 0.0) ?
        y[LUMPED_SMINN] : 0.0;
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *river;

        river = &pihm->river[i];

        river->ws.stage = MAX(y[RIVER(i)], 0.0);

#if defined(_BGC_) && !defined(_LUMPED_) && !defined(_LEACHING_)
        river->ns.streamn = (y[STREAMN(i)] >= 0.0) ? y[STREAMN(i)] : 0.0;
        river->ns.sminn = (y[RIVBEDN(i)] >= 0.0) ? y[RIVBEDN(i)] : 0.0;
#endif

        river->wf.rivflow[UP_CHANL2CHANL] = 0.0;

#if defined(_CYCLES_OBSOLETE_)
        river->ns.streamno3 = (y[STREAMNO3(i)] > 0.0) ? y[STREAMNO3(i)] : 0.0;
        river->ns.bedno3 = (y[RIVBEDNO3(i)] > 0.0) ? y[RIVBEDNO3(i)] : 0.0;

        river->ns.streamnh4 = (y[STREAMNH4(i)] > 0.0) ? y[STREAMNH4(i)] : 0.0;
        river->ns.bednh4 = (y[RIVBEDNH4(i)] > 0.0) ? y[RIVBEDNH4(i)] : 0.0;
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            river->chms.t_mole[k] = MAX(y[SOLUTE_RIVER(i, k)], 0.0);
        }
#endif
    }

    /*
     * PIHM Hydrology fluxes
     */
    Hydrol(pihm->elem, pihm->river, &pihm->ctrl);

#if defined(_BGC_)
    /*
     * Nitrogen transport fluxes
     */
# if defined(_LUMPED_)
    NLeachingLumped(pihm->elem, pihm->river);
# elif defined(_LEACHING_)
    NLeaching(pihm->elem);
# else
    NTransport(pihm->elem, pihm->river);
# endif
#endif

#if defined(_CYCLES_OBSOLETE_)
    /*
     * NO3 and NH4 transport fluxes
     */
    NTransport(t + (realtype)pihm->ctrl.tout[0] -
        (realtype)pihm->ctrl.tout[pihm->ctrl.cstep], pihm->elem, pihm->river);
#endif

#if defined(_RT_)
    /*
     * Aqueous species transport fluxes
     */
    SoluteConc(pihm->chemtbl, &pihm->rttbl, pihm->elem, pihm->river);
#endif

#if defined(_BGC_) || defined(_CYCLES_OBSOLETE_)
    SoluteTransp(0.0, 0.0, 0.0, pihm->elem, pihm->river);
#elif defined(_RT_)
    SoluteTranspt(pihm->rttbl.DiffCoe, pihm->rttbl.DispCoe,
        pihm->rttbl.Cementation, pihm->elem, pihm->river);
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
        elem_struct    *elem;

        elem = &pihm->elem[i];

        /*
         * Vertical water fluxes for surface and subsurface
         */
        dy[SURF(i)] += elem->wf.pcpdrp - elem->wf.infil - elem->wf.edir_surf;
        dy[UNSAT(i)] += elem->wf.infil - elem->wf.rechg - elem->wf.edir_unsat -
            elem->wf.ett_unsat;
        dy[GW(i)] += elem->wf.rechg - elem->wf.edir_gw - elem->wf.ett_gw;

#if defined(_FBR_)
        /*
         * Vertical water fluxes for fractured bedrock
         */
        dy[GW(i)] -= elem->wf.fbr_infil;

        dy[FBRUNSAT(i)] += elem->wf.fbr_infil - elem->wf.fbr_rechg;
        dy[FBRGW(i)] += elem->wf.fbr_rechg;
# if defined(_TGM_)
        dy[FBRGW(i)] -= elem->wf.fbr_discharge;
# endif
#endif

        /*
         * Horizontal water fluxes
         */
        for (j = 0; j < NUM_EDGE; j++)
        {
            dy[SURF(i)] -= elem->wf.ovlflow[j] / elem->topo.area;
            dy[GW(i)] -= elem->wf.subsurf[j] / elem->topo.area;
#if defined(_FBR_)
            dy[FBRGW(i)] -= elem->wf.fbrflow[j] / elem->topo.area;
#endif
        }

        dy[UNSAT(i)] /= elem->soil.porosity;
        dy[GW(i)] /= elem->soil.porosity;
#if defined(_FBR_)
        dy[FBRUNSAT(i)] /= elem->geol.porosity;
        dy[FBRGW(i)] /= elem->geol.porosity;
#endif

#if defined(_BGC_) && !defined(_LUMPED_)
# if !defined(_LEACHING_)
        /*
         * BGC N transport fluxes
         */
        dy[SURFN(i)] +=
            (elem->nf.ndep_to_sminn + elem->nf.nfix_to_sminn) / DAYINSEC -
            elem->nsol.infilflux;
        dy[SMINN(i)] += elem->nsol.infilflux + elem->nsol.snksrc;
# else
        dy[SMINN(i)] +=
            (elem->nf.ndep_to_sminn + elem->nf.nfix_to_sminn) / DAYINSEC +
            elem->nsol.snksrc;
# endif

        for (j = 0; j < NUM_EDGE; j++)
        {
# if !defined(_LEACHING_)
            dy[SURFN(i)] -= elem->nsol.ovlflux[j] / elem->topo.area;
# endif
            dy[SMINN(i)] -= elem->nsol.subflux[j] / elem->topo.area;
        }
#endif

#if defined(_CYCLES_OBSOLETE_)
        /*
         * Cycles NO3 and NH4 transport fluxes
         */
        dy[NO3(i)] += elem->no3sol.snksrc / DAYINSEC;
        dy[NH4(i)] += elem->nh4sol.snksrc / DAYINSEC;

        for (j = 0; j < NUM_EDGE; j++)
        {
            dy[NO3(i)] -= elem->no3sol.flux[j] / elem->topo.area;
            dy[NH4(i)] -= elem->nh4sol.flux[j] / elem->topo.area;
        }
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            dy[SOLUTE_SOIL(i, k)] += (elem->solute[k].infil +
                elem->chmf.react[k]) / elem->topo.area;
# if defined(_FBR_)
            dy[SOLUTE_SOIL(i, k)] -= elem->solute[k].fbr_infil /
                elem->topo.area;

            dy[SOLUTE_GEOL(i, k)] += (elem->solute[k].fbr_infil +
                elem->chmf.react_geol[k]) / elem->topo.area;
#  if defined(_TGM_)
            dy[SOLUTE_GEOL(i, k)] -= elem->solute[k].fbr_discharge /
                elem->topo.area;
#  endif
# endif

            for (j = 0; j < NUM_EDGE; j++)
            {
                dy[SOLUTE_SOIL(i, k)] -= elem->solute[k].subflux[j] /
                    elem->topo.area;
# if defined(_FBR_)
                dy[SOLUTE_GEOL(i, k)] -= elem->solute[k].fbrflow[j] /
                    elem->topo.area;
# endif
            }
        }
#endif
    }

#if defined(_BGC_) && defined(_LUMPED_)
    elem_struct    *elem;

    elem = &pihm->elem[LUMPED];

    dy[LUMPED_SMINN] +=
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
        river_struct   *river;

        river = &(pihm->river[i]);

        for (j = 0; j < NUM_RIVFLX; j++)
        {
            /* Note the limitation due to
             * d(v) / dt = a * dy / dt + y * da / dt
             * for cs other than rectangle */
            dy[RIVER(i)] -= river->wf.rivflow[j] / river->topo.area;
        }

#if defined(_BGC_) && !defined(_LUMPED_) && !defined(_LEACHING_)
        for (j = 0; j <= 6; j++)
        {
            dy[STREAMN(i)] -= river->nsol.flux[j] / river->topo.area;
        }

        dy[RIVBEDN(i)] += -river->nsol.flux[LEFT_AQUIF2AQUIF] -
            river->nsol.flux[RIGHT_AQUIF2AQUIF] -
            river->nsol.flux[DOWN_AQUIF2AQUIF] -
            river->nsol.flux[UP_AQUIF2AQUIF] + river->nsol.flux[CHANL_LKG];

        dy[RIVBEDN(i)] /= river->topo.area;
#endif

#if defined(_CYCLES_OBSOLETE_)
        for (j = 0; j <= 6; j++)
        {
            dy[STREAMNO3(i)] -= river->no3sol.flux[j] / river->topo.area;
            dy[STREAMNH4(i)] -= river->nh4sol.flux[j] / river->topo.area;
        }

        dy[RIVBEDNO3(i)] += -river->no3sol.flux[LEFT_AQUIF2AQUIF] -
            river->no3sol.flux[RIGHT_AQUIF2AQUIF] -
            river->no3sol.flux[DOWN_AQUIF2AQUIF] -
            river->no3sol.flux[UP_AQUIF2AQUIF] + river->no3sol.flux[CHANL_LKG];
        dy[RIVBEDNH4(i)] += -river->nh4sol.flux[LEFT_AQUIF2AQUIF] -
            river->nh4sol.flux[RIGHT_AQUIF2AQUIF] -
            river->nh4sol.flux[DOWN_AQUIF2AQUIF] -
            river->nh4sol.flux[UP_AQUIF2AQUIF] + river->nh4sol.flux[CHANL_LKG];

        dy[RIVBEDNO3(i)] /= river->topo.area;
        dy[RIVBEDNH4(i)] /= river->topo.area;
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            for (j = 0; j < NUM_RIVFLX; j++)
            {
                dy[SOLUTE_RIVER(i, k)] -= river->solute[k].flux[j] /
                    river->topo.area;
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

#if defined(_RT_)
    nsv += nsolute * (nelem + nriver);
#endif

#if defined(_BGC_)
# if defined(_LUMPED_)
    nsv += 1;
# else
    nsv += 2 * nelem + 2 * nriver;
# endif
#endif

#if defined(_CYCLES_OBSOLETE_)
    nsv += 2 * nelem + 4 * nriver;
#endif

#if defined(_FBR_)
    nsv += 2 * nelem;
#endif

#if defined(_FBR_) && defined(_RT_)
    nsv += nsolute * nelem;
#endif

    return nsv;
}
void SetCVodeParam(pihm_struct pihm, void *cvode_mem, SUNLinearSolver *sun_ls,
    N_Vector CV_Y)
{
    int             cv_flag;
    static int      reset;
    N_Vector        abstol;
#if defined(_BGC_) || defined(_CYCLES_OBSOLETE_)
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
#if defined(_BGC_) || defined(_CYCLES_OBSOLETE_) || defined(_RT_)
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

#if defined(_BGC_) || defined(_CYCLES_OBSOLETE_) || defined(_RT_)
void SetAbsTolArray(double hydrol_tol, double transp_tol, N_Vector abstol)
#else
void SetAbsTolArray(double hydrol_tol, N_Vector abstol)
#endif
{
    int             i;
    int             num_hydrol_var;

    num_hydrol_var = 3 * nelem + nriver;

#if defined(_FBR_)
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

#if defined(_BGC_) || defined(_CYCLES_OBSOLETE_) || defined(_RT_)
    /* Set absolute errors for nitrogen state variables */
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = num_hydrol_var; i < NumStateVar(); i++)
    {
        NV_Ith(abstol, i) = (realtype)transp_tol;
    }
#endif
}

void SolveCVode(const ctrl_struct *ctrl, double cputime, int *t,
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
        PIHMprintf(VL_NORMAL, "\n");
    }

    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "\033[1A\rStep = %s (t = %d)\n",
            pihm_time.str, *t);
        ProgressBar(progress);
    }
    else if (spinup_mode)
    {
        if (pihm_time.t % DAYINSEC == 0)
        {
            PIHMprintf(VL_NORMAL, "\033[1A\rStep = %s\n", pihm_time.str);
            ProgressBar(progress);
        }
    }
    else if (pihm_time.t % 3600 == 0)
    {
        PIHMprintf(VL_NORMAL, "\033[1A\rStep = %s (cputime %8.2f s)\n",
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

    if (nfails > ctrl->nncfn || niters >= ctrl->nnimax)
    {
        ctrl->maxstep /= ctrl->decr;
    }

    if (nfails == 0.0 && niters <= ctrl->nnimin)
    {
        ctrl->maxstep *= ctrl->incr;
    }

    ctrl->maxstep = (ctrl->maxstep < ctrl->stepsize) ?
        ctrl->maxstep : ctrl->stepsize;
    ctrl->maxstep = (ctrl->maxstep > ctrl->stmin) ?
        ctrl->maxstep : ctrl->stmin;

    /* Updates the upper bound on the magnitude of the step size */
    cv_flag = CVodeSetMaxStep(cvode_mem, (realtype)ctrl->maxstep);
    CheckCVodeFlag(cv_flag);

    nst0 = nst;
    ncfn0 = ncfn;
    nni0 = nni;
}
