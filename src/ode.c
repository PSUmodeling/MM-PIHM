#include "pihm.h"

int Ode(sunrealtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
{
    int             i;
    double         *y;
    double         *dy;
    pihm_struct    *pihm;
    elem_struct    *elem;
    river_struct   *river;

    y = NV_DATA(CV_Y);
    dy = NV_DATA(CV_Ydot);
    pihm = (pihm_struct *)pihm_data;

    elem = &pihm->elem[0];
    river = &pihm->river[0];

    // Initialization of RHS of ODEs
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
        elem[i].ws.surf = MAX(y[SURF(i)], 0.0);
        elem[i].ws.unsat = MAX(y[UNSAT(i)], 0.0);
        elem[i].ws.gw = MAX(y[GW(i)], 0.0);

#if defined(_BGC_)
        elem[i].ns.sminn = MAX(y[SOLUTE_SOIL(i, 0)], 0.0);
#endif

#if defined(_CYCLES_)
        elem[i].ps.no3 = MAX(y[SOLUTE_SOIL(i, NO3)], 0.0);
        elem[i].ps.nh4 = MAX(y[SOLUTE_SOIL(i, NH4)], 0.0);
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].ws.stage = MAX(y[RIVER(i)], 0.0);

#if defined(_BGC_)
        river[i].ns.streamn = MAX(y[SOLUTE_RIVER(i, 0)], 0.0);
#endif

#if defined(_CYCLES_)
        river[i].ns.no3 = MAX(y[SOLUTE_RIVER(i, NO3)], 0.0);
        river[i].ns.nh4 = MAX(y[SOLUTE_RIVER(i, NH4)], 0.0);
#endif
    }

    // PIHM Hydrology fluxes
    Hydrol(&pihm->ctrl, pihm->elem, pihm->river);

    // Calculate solute concentrations
#if defined(_TRANSPORT_) && !defined(_CYCLES_)
    SoluteConc(pihm->elem, pihm->river);
#elif defined(_CYCLES_)
    SoluteConc(t + (sunrealtype)pihm->ctrl.tout[0] - (sunrealtype)pihm->ctrl.tout[pihm->ctrl.cstep], pihm->elem, pihm->river);
#endif

#if defined(_TRANSPORT_)
    SoluteTranspt(pihm->elem, pihm->river);
#endif

    // Build RHS of ODEs
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        // Vertical water fluxes for surface and subsurface
        dy[SURF(i)] += elem[i].wf.pcpdrp - elem[i].wf.infil - elem[i].wf.edir_surf;
        dy[UNSAT(i)] += elem[i].wf.infil - elem[i].wf.recharge - elem[i].wf.edir_unsat - elem[i].wf.ett_unsat;
        dy[GW(i)] += elem[i].wf.recharge - elem[i].wf.edir_gw - elem[i].wf.ett_gw;

        // Horizontal water fluxes
        for (j = 0; j < NUM_EDGE; j++)
        {
            dy[SURF(i)] -= elem[i].wf.overland[j] / elem[i].topo.area;
            dy[GW(i)] -= elem[i].wf.subsurf[j] / elem[i].topo.area;
        }

        dy[UNSAT(i)] /= elem[i].soil.porosity;
        dy[GW(i)] /= elem[i].soil.porosity;

#if defined(_TRANSPORT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
# if defined(_CYCLES_)
            dy[SOLUTE_SOIL(i, k)] += elem[i].solute[k].infil + Profile(elem[i].ps.nlayers, elem[i].solute[k].snksrc);
# else
            dy[SOLUTE_SOIL(i, k)] += elem[i].solute[k].infil + elem[i].solute[k].snksrc;
# endif

            for (j = 0; j < NUM_EDGE; j++)
            {
                dy[SOLUTE_SOIL(i, k)] -= elem[i].solute[k].subflux[j] / elem[i].topo.area;
            }
        }
#endif
    }

    // ODEs for river segments
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;

        for (j = 0; j < NUM_RIVFLX; j++)
        {
            // Note the limitation due to d(v) / dt = a * dy / dt + y * da / dt for cs other than rectangle
            dy[RIVER(i)] -= river[i].wf.rivflow[j] / river[i].topo.area;
        }

#if defined(_TRANSPORT_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            for (j = 0; j < NUM_RIVFLX; j++)
            {
                dy[SOLUTE_RIVER(i, k)] -= river[i].solute[k].flux[j] / river[i].topo.area;
            }
        }
#endif
    }

    return 0;
}

int NumStateVar(void)
{
    // Return number of state variables
    int             nsv;

    nsv = 3 * nelem + nriver;

#if defined(_TRANSPORT_)
    nsv += nsolute * (nelem + nriver);
#endif

    return nsv;
}

void SetCVodeParam(pihm_struct *pihm, cvode_struct *cvode)
{
    int             cv_flag;
    static int      reset;
    N_Vector        abstol;
#if defined(_TRANSPORT_)
    const double    TRANSP_TOL = 1.0E-5;
#endif

    pihm->ctrl.maxstep = pihm->ctrl.stepsize;

    if (reset)
    {
        // When model spins-up and recycles forcing, use CVodeReInit to reset solver time, which does not allocates
        // memory
        cv_flag = CVodeReInit(cvode->memory_block, 0.0, cvode->CV_Y);
        CheckCVodeFlag(cv_flag);
    }
    else
    {
        cv_flag = CVodeInit(cvode->memory_block, Ode, 0.0, cvode->CV_Y);
        CheckCVodeFlag(cv_flag);
        reset = 1;

        cvode->sun_ls = SUNLinSol_SPGMR(cvode->CV_Y, SUN_PREC_NONE, 0, cvode->sunctx);

        // Attach the linear solver
        CVodeSetLinearSolver(cvode->memory_block, cvode->sun_ls, NULL);

        // When BGC or Cycles module is turned on, both water storage and transport variables are in the CVODE vector. A
        // vector of absolute tolerances is needed to specify different absolute tolerances for water storage variables
        // and transport variables
        abstol = N_VNew(NumStateVar(), cvode->sunctx);
#if defined(_TRANSPORT_)
        SetAbsTolArray(pihm->ctrl.abstol, TRANSP_TOL, abstol);
#else
        SetAbsTolArray(pihm->ctrl.abstol, abstol);
#endif

        cv_flag = CVodeSVtolerances(cvode->memory_block, (sunrealtype)pihm->ctrl.reltol, abstol);
        CheckCVodeFlag(cv_flag);

        N_VDestroy(abstol);

        // Specifies PIHM data block and attaches it to the main cvode memory block
        cv_flag = CVodeSetUserData(cvode->memory_block, cvode->user_data);
        CheckCVodeFlag(cv_flag);

        // Specifies the initial step size
        cv_flag = CVodeSetInitStep(cvode->memory_block, (sunrealtype)pihm->ctrl.initstep);
        CheckCVodeFlag(cv_flag);

        // Indicates if the BDF stability limit detection algorithm should be used
        cv_flag = CVodeSetStabLimDet(cvode->memory_block, SUNTRUE);
        CheckCVodeFlag(cv_flag);

        // Specifies an upper bound on the magnitude of the step size
        cv_flag = CVodeSetMaxStep(cvode->memory_block, (sunrealtype)pihm->ctrl.maxstep);
        CheckCVodeFlag(cv_flag);

        // Specifies the maximum number of steps to be taken by the solver in its attempt to reach the next output time
        cv_flag = CVodeSetMaxNumSteps(cvode->memory_block, pihm->ctrl.stepsize * 10);
        CheckCVodeFlag(cv_flag);
    }
}

#if defined(_TRANSPORT_)
void SetAbsTolArray(double hydrol_tol, double transp_tol, N_Vector abstol)
#else
void SetAbsTolArray(double hydrol_tol, N_Vector abstol)
#endif
{
    int             i;
    int             num_hydrol_var;

    num_hydrol_var = 3 * nelem + nriver;

    // Set absolute errors for hydrologic state variables
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < num_hydrol_var; i++)
    {
        NV_Ith(abstol, i) = (sunrealtype)hydrol_tol;
    }

#if defined(_TRANSPORT_)
    // Set absolute errors for solute state variables
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = num_hydrol_var; i < NumStateVar(); i++)
    {
        NV_Ith(abstol, i) = (sunrealtype)transp_tol;
    }
#endif
}

void SolveCVode(double cputime, const ctrl_struct *ctrl, int *t, cvode_struct *cvode)
{
    sunrealtype    solvert;
    sunrealtype    tout;
    pihm_t_struct   pihm_time;
    int             starttime;
    int             nextptr;
    int             cv_flag;
    double          progress;

    starttime = ctrl->starttime;
    nextptr = ctrl->tout[ctrl->cstep + 1];

    tout = (sunrealtype)(nextptr - starttime);

    // Specifies the value of the independent variable t past which the solution is not to proceed
    cv_flag = CVodeSetStopTime(cvode->memory_block, tout);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVode(cvode->memory_block, tout, cvode->CV_Y, &solvert, CV_NORMAL);
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
        pihm_printf(VL_NORMAL, "\033[1A\rStep = %s (t = %d)\n", pihm_time.str, *t);
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
        pihm_printf(VL_NORMAL, "\033[1A\rStep = %s (cputime %8.2f s)\n", pihm_time.str, cputime);
        ProgressBar(progress);
    }
}

void AdjCVodeMaxStep(cvode_struct *cvode, ctrl_struct *ctrl)
{
    // Variable CVODE max step (to reduce oscillations)
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

    // Gets the cumulative number of internal steps taken by the solver (total so far)
    cv_flag = CVodeGetNumSteps(cvode->memory_block, &nst);
    CheckCVodeFlag(cv_flag);

    // Gets the number of nonlinear convergence failures that have occurred
    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode->memory_block, &ncfn);
    CheckCVodeFlag(cv_flag);

    // Gets the number of nonlinear iterations performed
    cv_flag = CVodeGetNumNonlinSolvIters(cvode->memory_block, &nni);
    CheckCVodeFlag(cv_flag);

    nsteps = (double)(nst - nst0);
    nfails = (double)(ncfn - ncfn0) / nsteps;
    niters = (double)(nni - nni0) / nsteps;

    ctrl->maxstep /= (nfails > ctrl->nncfn || niters >= ctrl->nnimax) ? ctrl->decr : 1.0;

    ctrl->maxstep *= (nfails == 0.0 && niters <= ctrl->nnimin) ? ctrl->incr : 1.0;

    ctrl->maxstep = MIN(ctrl->maxstep, ctrl->stepsize);
    ctrl->maxstep = MAX(ctrl->maxstep, ctrl->stmin);

    // Updates the upper bound on the magnitude of the step size
    cv_flag = CVodeSetMaxStep(cvode->memory_block, (sunrealtype)ctrl->maxstep);
    CheckCVodeFlag(cv_flag);

    nst0  = nst;
    ncfn0 = ncfn;
    nni0  = nni;
}
