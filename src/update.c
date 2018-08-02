#include "pihm.h"

void Summary(elem_struct *elem, river_struct *river, N_Vector CV_Y,
    double stepsize)
{
    double         *y;
    int             i;

    y = NV_DATA(CV_Y);

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          subrunoff;

        elem[i].ws.surf = y[SURF(i)];
        elem[i].ws.unsat = y[UNSAT(i)];
        elem[i].ws.gw = y[GW(i)];

#if defined(_FBR_)
        elem[i].ws.fbr_unsat = y[FBRUNSAT(i)];
        elem[i].ws.fbr_gw = y[FBRGW(i)];
#endif

#if defined(_FBR_)
        MassBalance(&elem[i].ws, &elem[i].ws0, &elem[i].wf, &subrunoff,
            &elem[i].soil, &elem[i].geol, elem[i].topo.area, stepsize);
#else
        MassBalance(&elem[i].ws, &elem[i].ws0, &elem[i].wf, &subrunoff,
            &elem[i].soil, elem[i].topo.area, stepsize);
#endif

#if defined(_NOAH_)
        elem[i].wf.runoff2 = subrunoff;
# if defined(_FBR_)
        elem[i].wf.runoff2 += elem[i].wf.fbr_infil;
# endif

        elem[i].ps.nwtbl = FindWaterTable(elem[i].ps.sldpth, elem[i].ps.nsoil,
            elem[i].ws.gw, elem[i].ps.satdpth);

        CalcLatFlx(&elem[i].ps, &elem[i].wf);
#endif

        elem[i].ws0 = elem[i].ws;

#if defined(_BGC_) && !defined(_LUMPED_)
        elem[i].ns.surfn = (y[SURFN(i)] > 0.0) ? y[SURFN(i)] : 0.0;
        elem[i].ns.sminn = (y[SMINN(i)] > 0.0) ? y[SMINN(i)] : 0.0;

        elem[i].ns.nleached_snk += (elem[i].nt.surfn0 + elem[i].nt.sminn0) -
            (elem[i].ns.surfn + elem[i].ns.sminn) +
            elem[i].nf.ndep_to_sminn / DAYINSEC * stepsize +
            elem[i].nf.nfix_to_sminn / DAYINSEC * stepsize +
            elem[i].nsol.snksrc * stepsize;

        elem[i].nt.surfn0 = elem[i].ns.surfn;
        elem[i].nt.sminn0 = elem[i].ns.sminn;
#endif
    }

#if defined(_BGC_) && defined(_LUMPED_)
    elem[LUMPED].ns.sminn = (y[LUMPED_SMINN] > 0.0) ? y[LUMPED_SMINN] : 0.0;

    elem[LUMPED].ns.nleached_snk += (elem[LUMPED].nt.sminn0 -
        elem[LUMPED].ns.sminn) +
        elem[LUMPED].nf.ndep_to_sminn / DAYINSEC * stepsize +
        elem[LUMPED].nf.nfix_to_sminn / DAYINSEC * stepsize +
        elem[LUMPED].nsol.snksrc * stepsize;

    elem[LUMPED].nt.sminn0 = elem[LUMPED].ns.sminn;
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].ws.stage = y[RIVSTG(i)];
        river[i].ws.gw = y[RIVGW(i)];

        river[i].ws0 = river[i].ws;
    }
}

#if defined(_FBR_)
void MassBalance(const wstate_struct *ws, const wstate_struct *ws0,
    wflux_struct *wf, double *subrunoff, const soil_struct *soil,
    const geol_struct *geol, double area, double stepsize)
#else
void MassBalance(const wstate_struct *ws, const wstate_struct *ws0,
    wflux_struct *wf, double *subrunoff, const soil_struct *soil,
    double area, double stepsize)
#endif
{
    int             j;
    double          soilw0, soilw1;
#if defined(_FBR_)
    double          fbrw0, fbrw1;
    double          fbrrunoff;
#endif

    /*
     * Calculate infiltration based on mass conservation
     */
    soilw0 = ws0->gw + ws0->unsat;
    soilw0 = (soilw0 > soil->depth) ? soil->depth : soilw0;
    soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

    soilw1 = ws->gw + ws->unsat;
    soilw1 = (soilw1 > soil->depth) ? soil->depth : soilw1;
    soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;

    /* Subsurface runoff rate */
    *subrunoff = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        *subrunoff += wf->subsurf[j] / area;
    }

    wf->infil = (soilw1 - soilw0) * soil->porosity / stepsize + *subrunoff +
        wf->edir_unsat + wf->edir_gw + wf->ett_unsat + wf->ett_gw;

#if defined(_FBR_)
    fbrw0 = ws0->fbr_gw + ws0->fbr_unsat;
    fbrw1 = ws->fbr_gw + ws->fbr_unsat;

    fbrrunoff = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        fbrrunoff += wf->fbrflow[j] / area;
    }

    wf->fbr_infil = (fbrw1 - fbrw0) * geol->porosity / stepsize + fbrrunoff;

    wf->infil += wf->fbr_infil;
#endif

    if (wf->infil < 0.0)
    {
        *subrunoff -= wf->infil;
        wf->infil = 0.0;
    }
}
