#include "pihm.h"

int ODE(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
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
        elem->ws.fbr_unsat = (y[FBRUNSAT(i)] >= 0.0) ? y[FBRUNSAT(i)] : 0.0;
        elem->ws.fbr_gw = (y[FBRGW(i)] >= 0.0) ? y[FBRGW(i)] : 0.0;
#endif

#if defined(_BGC_) && !defined(_LUMPED_)
        elem->ns.surfn = (y[SURFN(i)] >= 0.0) ? y[SURFN(i)] : 0.0;
        elem->ns.sminn = (y[SMINN(i)] >= 0.0) ? y[SMINN(i)] : 0.0;
#endif

#if defined(_CYCLES_)
        elem->np.no3 = (y[NO3(i)] >= 0.0) ? y[NO3(i)] : 0.0;
        elem->np.nh4 = (y[NH4(i)] >= 0.0) ? y[NH4(i)] : 0.0;
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            elem->chms_unsat.t_mole[k] =
                (y[UNSAT_MOLE(i, k)] >= 0.0) ? y[UNSAT_MOLE(i, k)] : 0.0;
            elem->chms_gw.t_mole[k] =
                (y[GW_MOLE(i, k)] >= 0.0) ? y[GW_MOLE(i, k)] : 0.0;
# if defined(_FBR_)
            elem->chms_fbrunsat.t_mole[k] = MAX(y[FBRUNSAT_MOLE(i, k)], 0.0);
            elem->chms_fbrgw.t_mole[k] = MAX(y[FBRGW_MOLE(i, k)], 0.0);
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

        river->ws.stage = (y[RIVSTG(i)] >= 0.0) ? y[RIVSTG(i)] : 0.0;
        river->ws.gw = (y[RIVGW(i)] >= 0.0) ? y[RIVGW(i)] : 0.0;

#if defined(_BGC_) && !defined(_LUMPED_) && !defined(_LEACHING_)
        river->ns.streamn = (y[STREAMN(i)] >= 0.0) ? y[STREAMN(i)] : 0.0;
        river->ns.sminn = (y[RIVBEDN(i)] >= 0.0) ? y[RIVBEDN(i)] : 0.0;
#endif

        river->wf.rivflow[UP_CHANL2CHANL] = 0.0;
        river->wf.rivflow[UP_AQUIF2AQUIF] = 0.0;

#if defined(_CYCLES_)
        river->ns.streamno3 = (y[STREAMNO3(i)] > 0.0) ? y[STREAMNO3(i)] : 0.0;
        river->ns.bedno3 = (y[RIVBEDNO3(i)] > 0.0) ? y[RIVBEDNO3(i)] : 0.0;

        river->ns.streamnh4 = (y[STREAMNH4(i)] > 0.0) ? y[STREAMNH4(i)] : 0.0;
        river->ns.bednh4 = (y[RIVBEDNH4(i)] > 0.0) ? y[RIVBEDNH4(i)] : 0.0;
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            river->chms_stream.t_mole[k] =
                (y[STREAM_MOLE(i, k)] >= 0.0) ? y[STREAM_MOLE(i, k)] : 0.0;
            river->chms_rivbed.t_mole[k] =
                (y[RIVBED_MOLE(i, k)] >= 0.0) ? y[RIVBED_MOLE(i, k)] : 0.0;
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

#if defined(_CYCLES_)
    /*
     * NO3 and NH4 transport fluxes
     */
    NTransport(t + (realtype)pihm->ctrl.tout[0] -
        (realtype)pihm->ctrl.tout[pihm->ctrl.cstep], pihm->elem, pihm->river);
#endif

#if defined(_RT_)
    /*
     * Aquous species transport fluxes
     */
    Transport(pihm->chemtbl, &pihm->rttbl, pihm->elem, pihm->river);
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

        /* Check NAN errors for dy */
        CheckDy(dy[SURF(i)], "element", "surface water", i + 1, (double)t);
        CheckDy(dy[UNSAT(i)], "element", "unsat water", i + 1, (double)t);
        CheckDy(dy[GW(i)], "element", "groundwater", i + 1, (double)t);
#if defined(_FBR_)
        CheckDy(dy[FBRUNSAT(i)], "element", "fbr unsat", i + 1, (double)t);
        CheckDy(dy[FBRGW(i)], "element", "fbr groundwater", i + 1, (double)t);
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

        /* Check NAN errors for dy */
        CheckDy(dy[SURFN(i)], "element", "surface N", i + 1, (double)t);
        CheckDy(dy[SMINN(i)], "element", "soil mineral N", i + 1, (double)t);
#endif

#if defined(_CYCLES_)
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

        /* Check NAN errors for dy */
        CheckDy(dy[NO3(i)], "element", "NO3", i + 1, (double)t);
        CheckDy(dy[NH4(i)], "element", "NH4", i + 1, (double)t);
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            dy[UNSAT_MOLE(i, k)] += elem->chmf.infil[k] - elem->chmf.rechg[k] +
                elem->chmf.react_unsat[k];
            dy[GW_MOLE(i, k)] += elem->chmf.rechg[k] + elem->chmf.react_gw[k];
# if defined(_FBR_)
            dy[GW_MOLE(i, k)] -= elem->chmf.fbr_infil[k];

            dy[FBRUNSAT_MOLE(i, k)] += elem->chmf.fbr_infil[k] -
                elem->chmf.fbr_rechg[k] + elem->chmf.react_fbrunsat[k];
            dy[FBRGW_MOLE(i, k)] +=
                elem->chmf.fbr_rechg[k] + elem->chmf.react_fbrgw[k];
#  if defined(_TGM_)
            dy[FBRGW_MOLE(i, k)] -= elem->chmf.fbr_discharge[k];
#  endif
# endif

            for (j = 0; j < NUM_EDGE; j++)
            {
                dy[UNSAT_MOLE(i, k)] -= elem->chmf.unsatflux[j][k];
                dy[GW_MOLE(i, k)] -= elem->chmf.subflux[j][k];
# if defined(_FBR_)
                dy[FBRUNSAT_MOLE(i, k)] -= elem->chmf.fbr_unsatflux[j][k];
                dy[FBRGW_MOLE(i, k)] -= elem->chmf.fbrflow[j][k];
# endif
            }

            CheckDy(dy[UNSAT_MOLE(i, k)], "element", "unsat chem", i + 1,
                (double)t);
            CheckDy(dy[GW_MOLE(i, k)], "element", "gw chem", i + 1, (double)t);
# if defined(_FBR_)
            CheckDy(dy[FBRUNSAT_MOLE(i, k)], "element", "fbr unsat chem", i + 1,
                (double)t);
            CheckDy(dy[FBRGW_MOLE(i, k)], "element", "fbr gw chem", i + 1,
                (double)t);
# endif
        }
#endif
    }

#if defined(_BGC_) && defined(_LUMPED_)
    elem_struct    *elem;

    elem = &pihm->elem[LUMPED];

    dy[LUMPED_SMINN] +=
        (elem->nf.ndep_to_sminn + elem->nf.nfix_to_sminn) / DAYINSEC +
        elem->nsol.snksrc;

    /* Check NAN errors for dy */
    CheckDy(dy[LUMPED_SMINN], "lumped", "soil mineral N", LUMPED + 1,
        (double)t);
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

        for (j = 0; j <= 6; j++)
        {
            /* Note the limitation due to
             * d(v) / dt = a * dy / dt + y * da / dt
             * for cs other than rectangle */
            dy[RIVSTG(i)] -= river->wf.rivflow[j] / river->topo.area;
        }

#if defined(_FBR_) && defined(_TGM_)
        dy[RIVSTG(i)] -= (river->wf.rivflow[LEFT_FBR2CHANL] +
            river->wf.rivflow[RIGHT_FBR2CHANL]) / river->topo.area;
#endif

        dy[RIVGW(i)] += -river->wf.rivflow[LEFT_AQUIF2AQUIF] -
            river->wf.rivflow[RIGHT_AQUIF2AQUIF] -
            river->wf.rivflow[DOWN_AQUIF2AQUIF] -
            river->wf.rivflow[UP_AQUIF2AQUIF] + river->wf.rivflow[CHANL_LKG];

        dy[RIVGW(i)] /= river->matl.porosity * river->topo.area;

        /* Check NAN errors for dy */
        CheckDy(dy[RIVSTG(i)], "river", "stage", i + 1, (double)t);
        CheckDy(dy[RIVGW(i)], "river", "groundwater", i + 1, (double)t);

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

        /* Check NAN errors for dy */
        CheckDy(dy[STREAMN(i)], "river", "stream N", i + 1, (double)t);
        CheckDy(dy[RIVBEDN(i)], "river", "bed mineral N", i + 1, (double)t);
#endif

#if defined(_CYCLES_)
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

        /* Check NAN errors for dy */
        CheckDy(dy[STREAMNO3(i)], "river", "stream NO3", i + 1, (double)t);
        CheckDy(dy[RIVBEDNO3(i)], "river", "bed NO3", i + 1, (double)t);
        CheckDy(dy[STREAMNH4(i)], "river", "stream NH4", i + 1, (double)t);
        CheckDy(dy[RIVBEDNH4(i)], "river", "bed NH4", i + 1, (double)t);
#endif

#if defined(_RT_)
        int             k;

        for (k = 0; k < NumSpc; k++)
        {
            for (j = 0; j <= 6; j++)
            {
                dy[STREAM_MOLE(i, k)] -= river->chmf.flux[j][k];
            }

# if defined(_FBR_) && defined(_TGM_)
            dy[STREAM_MOLE(i, k)] -= river->chmf.flux[LEFT_FBR2CHANL][k] +
                river->chmf.flux[RIGHT_FBR2CHANL][k];
# endif

            dy[RIVBED_MOLE(i, k)] += -river->chmf.flux[LEFT_AQUIF2AQUIF][k] -
                river->chmf.flux[RIGHT_AQUIF2AQUIF][k] -
                river->chmf.flux[DOWN_AQUIF2AQUIF][k] -
                river->chmf.flux[UP_AQUIF2AQUIF][k] +
                river->chmf.flux[CHANL_LKG][k];

            CheckDy(dy[STREAM_MOLE(i, k)], "river", "stream chem", i + 1, (double)t);
            CheckDy(dy[RIVBED_MOLE(i, k)], "river", "bed chem", i + 1, (double)t);
        }
#endif
    }

    return 0;
}

void CheckDy(double dy, const char *type, const char *varname, int ind,
    double t)
{
    if (isnan(dy))
    {
        PIHMprintf(VL_ERROR,
            "Error: NAN error for %s %d (%s) at %lf\n", type, ind, varname, t);
        PIHMexit(EXIT_FAILURE);
    }
}

int NumStateVar(void)
{
    /*
     * Return number of state variables
     */
    int             nsv;

    nsv = 3 * nelem + 2 * nriver;

#if defined(_RT_)
    nsv += NumSpc * (2 * nelem + 2 * nriver);
#endif

#if defined(_BGC_)
# if defined(_LUMPED_)
    nsv += 1;
# else
    nsv += 2 * nelem + 2 * nriver;
# endif
#endif

#if defined(_CYCLES_)
    nsv += 2 * nelem + 4 * nriver;
#endif

#if defined(_FBR_)
    nsv += 2 * nelem;
#endif

#if defined(_FBR_) && defined(_RT_)
    nsv += NumSpc * 2 * nelem;
#endif

    return nsv;
}
void SetCVodeParam(pihm_struct pihm, void *cvode_mem, SUNLinearSolver *sun_ls,
    N_Vector CV_Y)
{
    int             cv_flag;
    static int      reset;
#if defined(_BGC_) || defined(_CYCLES_)
    N_Vector        abstol;
    const double    SMINN_TOL = 1.0E-5;
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
        cv_flag = CVodeInit(cvode_mem, ODE, 0.0, CV_Y);
        CheckCVodeFlag(cv_flag);
        reset = 1;

        *sun_ls = SUNLinSol_SPGMR(CV_Y, PREC_NONE, 0);

        /* Attach the linear solver */
        CVodeSetLinearSolver(cvode_mem, *sun_ls, NULL);

#if defined(_BGC_) || defined(_CYCLES_)
        /* When BGC module is turned on, both water storage and nitrogen storage
         * variables are in the CVODE vector. A vector of absolute tolerances is
         * needed to specify different absolute tolerances for water storage
         * variables and nitrogen storage variables */
        abstol = N_VNew(NumStateVar());
        SetAbsTol(pihm->ctrl.abstol, SMINN_TOL, abstol);

        cv_flag = CVodeSVtolerances(cvode_mem, (realtype)pihm->ctrl.reltol,
                abstol);
        CheckCVodeFlag(cv_flag);

        N_VDestroy(abstol);
#else
        cv_flag = CVodeSStolerances(cvode_mem, (realtype)pihm->ctrl.reltol,
                (realtype)pihm->ctrl.abstol);
        CheckCVodeFlag(cv_flag);
#endif

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

void SetAbsTol(double hydrol_tol, double sminn_tol, N_Vector abstol)
{
    int             i;

    /* Set absolute errors for hydrologic state variables */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < 3 * nelem + 2 * nriver; i++)
    {
        NV_Ith(abstol, i) = (realtype)hydrol_tol;
    }

    /* Set absolute errors for nitrogen state variables */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 3 * nelem + 2 * nriver; i < NumStateVar(); i++)
    {
        NV_Ith(abstol, i) = (realtype)sminn_tol;
    }
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
    int             progress;

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

    progress = (int)(((double)ctrl->cstep + 1) / (double)ctrl->nstep * 100.0);

    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "\r Step = %s (%d)", pihm_time.str, *t);
        ProgressBar(progress);
    }
    else if (spinup_mode)
    {
        if (pihm_time.t % DAYINSEC == 0)
        {
            PIHMprintf(VL_NORMAL, "\r Step = %s", pihm_time.str);
            ProgressBar(progress);
        }
    }
    else if (pihm_time.t % 3600 == 0)
    {
        PIHMprintf(VL_NORMAL,
            "\r Step = %s (cputime %.2f)", pihm_time.str, cputime);
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
