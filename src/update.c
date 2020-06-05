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

        elem[i].ps.nwtbl = FindWaterTable(elem[i].ps.soil_depth, elem[i].ps.nlayers,
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

#if defined(_CYCLES_)
        elem[i].ps.no3 = MAX(y[SOLUTE_SOIL(i, NO3)], 0.0);
        elem[i].ps.nh4 = MAX(y[SOLUTE_SOIL(i, NH4)], 0.0);

        UpdateNProfile(stepsize, &elem[i].soil, &elem[i].ws, &elem[i].ns,
            elem[i].solute, elem[i].ns.no3, elem[i].ns.nh4, &elem[i].ps);

        elem[i].ns0 = elem[i].ns;
        elem[i].ps.no3_prev = elem[i].ps.no3;
        elem[i].ps.nh4_prev = elem[i].ps.nh4;
#endif

#if defined(_RT_)
        int             k;
        double          storage;

        storage = (elem[i].ws.gw + elem[i].ws.unsat) * elem[i].soil.porosity +
            elem[i].soil.smcmin * elem[i].soil.depth;

        for (k = 0; k < nsolute; k++)
        {
            elem[i].chms.tot_mol[k] = MAX(y[SOLUTE_SOIL(i, k)], 0.0);

            /* Calculate concentrations */
            elem[i].chms.tot_conc[k] = elem[i].chms.tot_mol[k] / storage;
            elem[i].chms.tot_conc[k] = MAX(elem[i].chms.tot_conc[k], ZERO_CONC);
        }

# if defined(_FBR_)
        storage = (elem[i].ws.fbr_unsat + elem[i].ws.fbr_gw) *
            elem[i].geol.porosity + elem[i].geol.smcmin * elem[i].geol.depth;
        storage = MAX(storage, 0.0);

        for (k = 0; k < nsolute; k++)
        {
            elem[i].chms_geol.tot_mol[k] = MAX(y[SOLUTE_GEOL(i, k)], 0.0);

            /* Calculate concentrations */
            elem[i].chms_geol.tot_conc[k] = (storage > 0.0) ?
                elem[i].chms_geol.tot_mol[k] / storage : 0.0;
            elem[i].chms_geol.tot_conc[k] =
                MAX(elem[i].chms_geol.tot_conc[k], ZERO_CONC);
        }
# endif
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
        river[i].ws.stage = y[RIVER(i)];

        river[i].ws0 = river[i].ws;

#if defined(_CYCLES_)
        river[i].ns.no3 = MAX(y[SOLUTE_RIVER(i, NO3)], 0.0);
        river[i].ns.nh4 = MAX(y[SOLUTE_RIVER(i, NH4)], 0.0);
#endif

#if defined(_RT_)
        int             k;
        double          storage;

        storage = MAX(river[i].ws.stage, 0.0);

        for (k = 0; k < nsolute; k++)
        {
            river[i].chms.tot_mol[k] = MAX(y[SOLUTE_RIVER(i, k)], 0.0);

            /* Calculate concentrations */
            river[i].chms.tot_conc[k] = (storage > DEPTHR) ?
                river[i].chms.tot_mol[k] / storage : ZERO_CONC;
            river[i].chms.tot_conc[k] = MAX(river[i].chms.tot_conc[k], ZERO_CONC);
        }
#endif
    }
}

#if defined(_FBR_)
void MassBalance(const wstate_struct *ws, const wstate_struct *ws0,
    wflux_struct *wf, double *subrunoff, const soil_struct *soil,
    const soil_struct *geol, double area, double stepsize)
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
     * The MassBalance subroutine adjusts model infiltration rate and recharge
     * rate based on mass balance in the soil column. This subroutine also
     * calculates an equivalent infiltration rate, which is used as the boundary
     * condition in the Noah soil moisture calculation. Unlike the infiltration
     * rate, the equivalent infiltration rate is calculated without considering
     * oversturation, and is based on the mass balance of the whole soil column,
     * instead of just the unsaturated zone.
     */

    /* Adjust recharge */
    soilw0 = ws0->gw;
    soilw1 = ws->gw;

    *subrunoff = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        *subrunoff += wf->subsurf[j] / area;
    }

    wf->rechg = (soilw1 - soilw0) * soil->porosity / stepsize + *subrunoff +
        wf->edir_gw + wf->ett_gw;

    /* Adjust infiltration */
    soilw0 = ws0->unsat;
    soilw1 = ws->unsat;

    wf->infil = (soilw1 - soilw0) * soil->porosity / stepsize + wf->rechg +
        wf->edir_unsat + wf->ett_unsat;

#if defined(_FBR_)
    /* Adjust bedrock recharge */
    fbrw0 = ws0->fbr_gw;
    fbrw1 = ws->fbr_gw;

    fbrrunoff = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        fbrrunoff += wf->fbrflow[j] / area;
    }

    wf->fbr_rechg = (fbrw1 - fbrw0) * geol->porosity / stepsize + fbrrunoff;

    /* Adjust bedrock infiltration */
    fbrw0 = ws0->fbr_unsat;
    fbrw1 = ws->fbr_unsat;

    wf->fbr_infil = (fbrw1 - fbrw0) * geol->porosity / stepsize + wf->fbr_rechg;

    /* Further adjust soil infiltration and recharge rate to take into account
     * bedrock leakage */
    wf->rechg += wf->fbr_infil;
    wf->infil += wf->fbr_infil;
#endif

    /*
     * Calculate equivalent infiltration based on mass conservation
     */
    soilw0 = ws0->gw + ws0->unsat;
    soilw0 = (soilw0 > soil->depth) ? soil->depth : soilw0;
    soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

    soilw1 = ws->gw + ws->unsat;
    soilw1 = (soilw1 > soil->depth) ? soil->depth : soilw1;
    soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;

    wf->eqv_infil = (soilw1 - soilw0) * soil->porosity / stepsize + *subrunoff +
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

    wf->eqv_infil += wf->fbr_infil;
#endif

    if (wf->eqv_infil < 0.0)
    {
        *subrunoff -= wf->eqv_infil;
        wf->eqv_infil = 0.0;
    }
}
