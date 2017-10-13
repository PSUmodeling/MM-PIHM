#include "pihm.h"

int ODE(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
{
    int             i;
    double         *y;
    double         *dy;
    double          dt;
    pihm_struct     pihm;

    y = NV_DATA(CV_Y);
    dy = NV_DATA(CV_Ydot);
    pihm = (pihm_struct)pihm_data;

    dt = (double)pihm->ctrl.stepsize;

    /*
     * Initialization of temporary state variables
     */
    for (i = 0; i < NumStateVar(); i++)
    {
        dy[i] = 0.0;
    }

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct    *elem;

        elem = &pihm->elem[i];

        elem->ws.surf = (y[SURF(i)] >= 0.0) ? y[SURF(i)] : 0.0;
        elem->ws.unsat = (y[UNSAT(i)] >= 0.0) ? y[UNSAT(i)] : 0.0;
        elem->ws.gw = (y[GW(i)] >= 0.0) ? y[GW(i)] : 0.0;

#ifdef _FBR_
        elem->ws.fbr_unsat = (y[FBRUNSAT(i)] >= 0.0) ? y[FBRUNSAT(i)] : 0.0;
        elem->ws.fbr_gw = (y[FBRGW(i)] >= 0.0) ? y[FBRGW(i)] : 0.0;
#endif

#ifdef _BGC_
        elem->ns.surfn = (y[SURFN(i)] >= 0.0) ? y[SURFN(i)] : 0.0;
        elem->ns.sminn = (y[SMINN(i)] >= 0.0) ? y[SMINN(i)] : 0.0;
#endif
    }

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *river;
        river = &pihm->river[i];

        river->ws.stage = (y[RIVSTG(i)] >= 0.0) ? y[RIVSTG(i)] : 0.0;
        river->ws.gw = (y[RIVGW(i)] >= 0.0) ? y[RIVGW(i)] : 0.0;

#ifdef _BGC_
        river->ns.streamn = (y[STREAMN(i)] >= 0.0) ? y[STREAMN(i)] : 0.0;
        river->ns.sminn = (y[RIVBEDN(i)] >= 0.0) ? y[RIVBEDN(i)] : 0.0;
#endif

        river->wf.rivflow[UP_CHANL2CHANL] = 0.0;
        river->wf.rivflow[UP_AQUIF2AQUIF] = 0.0;
    }

    /*
     * PIHM Hydrology fluxes
     */
    Hydrol(pihm->elem, pihm->river, &pihm->ctrl);

#ifdef _BGC_
    /*
     * Nitrogen transport fluxes
     */
    NTransport(pihm->elem, pihm->river);
#endif

    /*
     * Build RHS of ODEs
     */
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        elem_struct    *elem;

        elem = &pihm->elem[i];

        dy[SURF(i)] += elem->wf.pcpdrp - elem->wf.infil - elem->wf.edir_surf;
        dy[UNSAT(i)] += elem->wf.infil - elem->wf.rechg - elem->wf.edir_unsat -
            elem->wf.ett_unsat;
        dy[GW(i)] += elem->wf.rechg - elem->wf.edir_gw - elem->wf.ett_gw;

#ifdef _FBR_
        dy[GW(i)] -= elem->wf.fbr_infil;

        dy[FBRUNSAT(i)] += elem->wf.fbr_infil - elem->wf.fbr_rechg;
        dy[FBRGW(i)] += elem->wf.fbr_rechg;
#endif

        for (j = 0; j < NUM_EDGE; j++)
        {
            dy[SURF(i)] -= elem->wf.ovlflow[j] / elem->topo.area;
            dy[GW(i)] -= elem->wf.subsurf[j] / elem->topo.area;
#ifdef _FBR_
            dy[FBRGW(i)] -= elem->wf.fbrflow[j] / elem->topo.area;
#endif
        }

        dy[UNSAT(i)] /= elem->soil.porosity;
        dy[GW(i)] /= elem->soil.porosity;
#ifdef _FBR_
        dy[FBRUNSAT(i)] /= elem->geol.porosity;
        dy[FBRGW(i)] /= elem->geol.porosity;
#endif

        /* Check NAN errors for dy */
        CheckDy(dy[SURF(i)], "element", "surface water", i + 1, (double)t);
        CheckDy(dy[UNSAT(i)], "element", "unsat water", i + 1, (double)t);
        CheckDy(dy[GW(i)], "element", "groundwater", i + 1, (double)t);
#ifdef _FBR_
        CheckDy(dy[FBRUNSAT(i)], "element", "fbr unsat water", i + 1,
            (double)t);
        CheckDy(dy[GW(i)], "element", "fbr groundwater", i + 1, (double)t);
#endif

#ifdef _BGC_
        dy[SURFN(i)] +=
            (elem->nf.ndep_to_sminn + elem->nf.nfix_to_sminn) / DAYINSEC -
            elem->nsol.infilflux;
        dy[SMINN(i)] += elem->nsol.infilflux + elem->nsol.snksrc;

        for (j = 0; j < NUM_EDGE; j++)
        {
            dy[SURFN(i)] -= elem->nsol.ovlflux[j] / elem->topo.area;
            dy[SMINN(i)] -= elem->nsol.subflux[j] / elem->topo.area;
        }

        /* Check NAN errors for dy */
        CheckDy(dy[SURFN(i)], "element", "surface N", i + 1, (double)t);
        CheckDy(dy[SMINN(i)], "element", "soil mineral N", i + 1, (double)t);
#endif
    }

#ifdef _OPENMP
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

        dy[RIVGW(i)] += -river->wf.rivflow[LEFT_AQUIF2AQUIF] -
            river->wf.rivflow[RIGHT_AQUIF2AQUIF] -
            river->wf.rivflow[DOWN_AQUIF2AQUIF] -
            river->wf.rivflow[UP_AQUIF2AQUIF] + river->wf.rivflow[CHANL_LKG];

        dy[RIVGW(i)] /= river->matl.porosity * river->topo.area;

        /* Check NAN errors for dy */
        CheckDy(dy[RIVSTG(i)], "river", "stage", i + 1, (double)t);
        CheckDy(dy[RIVGW(i)], "river", "groundwater", i + 1, (double)t);

#ifdef _BGC_
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
    }

    return 0;
}

void CheckDy(double dy, const char *type, const char *varname, int ind, double t)
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

#ifdef _BGC_
    nsv = 5 * nelem + 4 * nriver;
#elif _FBR_
    nsv = 5 * nelem + 2 * nriver;
#else
    nsv = 3 * nelem + 2 * nriver;
#endif

    return nsv;
}
void SetCVodeParam(pihm_struct pihm, void *cvode_mem, N_Vector CV_Y)
{
    int             flag;

    pihm->ctrl.maxstep = pihm->ctrl.stepsize;

    flag = CVodeInit(cvode_mem, ODE, 0.0, CV_Y);
    flag = CVodeSStolerances(cvode_mem, (realtype)pihm->ctrl.reltol,
        (realtype)pihm->ctrl.abstol);
    flag = CVodeSetUserData(cvode_mem, pihm);
    flag = CVodeSetInitStep(cvode_mem, (realtype)pihm->ctrl.initstep);
    flag = CVodeSetStabLimDet(cvode_mem, TRUE);
    flag = CVodeSetMaxStep(cvode_mem, (realtype)pihm->ctrl.maxstep);
    flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
}

void SolveCVode(int starttime, int *t, int nextptr, double cputime,
    void *cvode_mem, N_Vector CV_Y)
{
    realtype        solvert;
    realtype        tout;
    pihm_t_struct   pihm_time;
    int             flag;

    tout = (realtype)(nextptr - starttime);

    flag = CVodeSetStopTime(cvode_mem, tout);
    flag = CVode(cvode_mem, tout, CV_Y, &solvert, CV_NORMAL);

    *t = (int)round(solvert) + starttime;

    pihm_time = PIHMTime(*t);

    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, " Step = %s (%d)\n", pihm_time.str, *t);
    }
    else if (spinup_mode)
    {
        if (pihm_time.t % DAYINSEC == 0)
        {
            PIHMprintf(VL_NORMAL, " Step = %s\n", pihm_time.str);
        }
    }
    else if (pihm_time.t % 3600 == 0)
    {
        PIHMprintf(VL_NORMAL,
            " Step = %s (cputime %f)\n", pihm_time.str, cputime);
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
    int             flag;
    double          nsteps;
    double          nfails;
    double          niters;

    flag = CVodeGetNumSteps(cvode_mem, &nst);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);

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

    flag = CVodeSetMaxStep(cvode_mem, (realtype)ctrl->maxstep);

    nst0 = nst;
    ncfn0 = ncfn;
    nni0 = nni;
}
