#include "pihm.h"

void Transport(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[])
{
    int             i;

    /*
     * Calculate chemical concentrations
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k, kk;
        double          strg_gw;
        double          strg_unsat;

        strg_gw = GWStrg(elem[i].soil.depth, elem[i].soil.smcmax,
            elem[i].soil.smcmin, elem[i].ws.gw);
        strg_unsat = UnsatWaterStrg(elem[i].soil.depth, elem[i].soil.smcmax,
            elem[i].soil.smcmin, elem[i].ws.gw, elem[i].ws.unsat);

        for (k = 0; k < NumSpc; k++)
        {
            /* Initialize chemical fluxes */
            elem[i].chmf.infil[k] = 0.0;
            elem[i].chmf.rechg[k] = 0.0;

            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].chmf.unsatflux[j][k] = 0.0;
                elem[i].chmf.subflux[j][k] = 0.0;
            }

            /* Calculate concentrations */
            elem[i].chms_unsat.t_conc[k] = (strg_unsat > DEPTHR) ?
                elem[i].chms_unsat.t_mole[k] / strg_unsat : 0.0;

            elem[i].chms_gw.t_conc[k] = (strg_gw > DEPTHR) ?
                elem[i].chms_gw.t_mole[k] / strg_gw : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->NumSsc; kk++)
                {
                    if ((rttbl->Totalconc[k][kk + rttbl->NumStc] != 0) &&
                        (chemtbl[kk + rttbl->NumStc].itype != AQUEOUS))
                    {
                        elem[i].chms_gw.t_conc[k] -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms_gw.s_conc[kk];
                        elem[i].chms_unsat.t_conc[k] -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms_unsat.s_conc[kk];
                    }
                }
            }

            elem[i].chms_unsat.t_conc[k] =
                MAX(elem[i].chms_unsat.t_conc[k], 0.0);
            elem[i].chms_gw.t_conc[k] =
                MAX(elem[i].chms_gw.t_conc[k], 0.0);
        }

#if defined(_FBR_)
        strg_gw = GWStrg(elem[i].geol.depth, elem[i].geol.smcmax,
            elem[i].geol.smcmin, elem[i].ws.fbr_gw);
        strg_unsat = UnsatWaterStrg(elem[i].geol.depth, elem[i].geol.smcmax,
            elem[i].geol.smcmin, elem[i].ws.fbr_gw, elem[i].ws.fbr_unsat);

        for (k = 0; k < NumSpc; k++)
        {
            /* Initialize chemical fluxes */
            elem[i].chmf.fbr_infil[k] = 0.0;
            elem[i].chmf.fbr_rechg[k] = 0.0;

            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].chmf.fbr_unsatflux[j][k] = 0.0;
                elem[i].chmf.fbrflow[j][k] = 0.0;
            }

            /* Calculate concentrations */
            elem[i].chms_fbrunsat.t_conc[k] = (strg_unsat > DEPTHR) ?
                elem[i].chms_fbrunsat.t_mole[k] / strg_unsat : 0.0;

            elem[i].chms_fbrgw.t_conc[k] = (strg_gw > DEPTHR) ?
                elem[i].chms_fbrgw.t_mole[k] / strg_gw : 0.0;

            if (chemtbl[k].mtype == MIXED_MA)
            {
                for (kk = 0; kk < rttbl->NumSsc; kk++)
                {
                    if ((rttbl->Totalconc[k][kk + rttbl->NumStc] != 0) &&
                        (chemtbl[kk + rttbl->NumStc].itype != AQUEOUS))
                    {
                        elem[i].chms_fbrgw.t_conc[k] -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms_fbrgw.s_conc[kk];
                        elem[i].chms_fbrunsat.t_conc[k] -=
                            rttbl->Totalconc[k][kk + rttbl->NumStc] *
                            elem[i].chms_fbrunsat.s_conc[kk];
                    }
                }
            }

            elem[i].chms_fbrunsat.t_conc[k] =
                MAX(elem[i].chms_fbrunsat.t_conc[k], 0.0);
            elem[i].chms_fbrgw.t_conc[k] =
                MAX(elem[i].chms_fbrgw.t_conc[k], 0.0);
        }
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j, k;
        double          strg_rivbed;
        double          strg_stream;

        strg_rivbed = RivBedStrg(&river[i].matl, &river[i].ws);
        strg_stream = river[i].ws.stage;

        for (k = 0; k < NumSpc; k++)
        {
            /* Initialize chemical fluxes */
            for (j = 0; j < NUM_RIVFLX; j++)
            {
                river[i].chmf.flux[j][k] = 0.0;
            }

            /* Calculate concentrations */
            river[i].chms_stream.t_conc[k] = (strg_stream > DEPTHR) ?
                river[i].chms_stream.t_mole[k] / strg_stream : 0.0;
            river[i].chms_stream.t_conc[k] =
                MAX(river[i].chms_stream.t_conc[k], 0.0);

            river[i].chms_rivbed.t_conc[k] = (strg_rivbed > DEPTHR) ?
                river[i].chms_rivbed.t_mole[k] / strg_rivbed : 0.0;
            river[i].chms_rivbed.t_conc[k] =
                MAX(river[i].chms_rivbed.t_conc[k], 0.0);
        }
    }

    /*
     * Calculate chemical fluxes
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k;
        elem_struct    *nabr;

        for (k = 0; k < NumSpc; k++)
        {
            /* Infiltration */
            elem[i].chmf.infil[k] = elem[i].wf.infil * elem[i].topo.area *
                ((elem[i].wf.infil > 0.0) ?
                elem[i].prcps.t_conc[k] * rttbl->Condensation : 0.0);

            /* Interface between unsaturated zone and groundwater */
            elem[i].chmf.rechg[k] = AdvDiffDisp(chemtbl[k].DiffCoe,
                chemtbl[k].DispCoe, rttbl->Cementation,
                elem[i].chms_unsat.t_conc[k], elem[i].chms_gw.t_conc[k],
                elem[i].soil.smcmax, 0.5 * elem[i].soil.depth,
                elem[i].topo.area, elem[i].wf.rechg * elem[i].topo.area);

            /* Element to element */
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] == 0)
                {
                    elem[i].chmf.unsatflux[j][k] = 0.0;
                    elem[i].chmf.subflux[j][k] = 0.0;
                }
                else if (elem[i].nabr_river[j] == 0)
                {
                    nabr = &elem[elem[i].nabr[j] - 1];

                    /* Unsaturated zone diffusion */
                    elem[i].chmf.unsatflux[j][k] =
                        AdvDiffDisp(chemtbl[k].DiffCoe, chemtbl[k].DispCoe,
                        rttbl->Cementation, elem[i].chms_unsat.t_conc[k],
                        nabr->chms_unsat.t_conc[k],
                        0.5 * elem[i].soil.smcmax + 0.5 * nabr->soil.smcmax,
                        elem[i].topo.nabrdist[j],
                        0.5 * MAX(elem[i].soil.depth - elem[i].ws.gw, 0.0) +
                        0.5 * MAX(nabr->soil.depth - nabr->ws.gw, 0.0), 0.0);

                    /* Groundwater advection, diffusion, and dispersion */
                    elem[i].chmf.subflux[j][k] =
                        AdvDiffDisp(chemtbl[k].DiffCoe, chemtbl[k].DispCoe,
                        rttbl->Cementation, elem[i].chms_gw.t_conc[k],
                        nabr->chms_gw.t_conc[k],
                        0.5 * elem[i].soil.smcmax + 0.5 * nabr->soil.smcmax,
                        elem[i].topo.nabrdist[j],
                        0.5 * MAX(elem[i].ws.gw, 0.0) +
                        0.5 * MAX(nabr->ws.gw, 0.0),
                        elem[i].wf.subsurf[j]);
                }
                else
                {
                    /* No flux between unsaturated zone and river.
                     * River-groundwater interactions are calculated later */
                    elem[i].chmf.unsatflux[j][k] = 0.0;
                }
            }   /* End of element to element */

#if defined(_FBR_)
            /* Bedrock infiltration */
            elem[i].chmf.fbr_infil[k] = elem[i].wf.fbr_infil *
                elem[i].topo.area * ((elem[i].wf.fbr_infil > 0.0) ?
                elem[i].chms_gw.t_conc[k] : elem[i].chms_fbrgw.t_conc[k]);

            /* Interface between unsaturated bedrock and deep groundwater */
            elem[i].chmf.fbr_rechg[k] = AdvDiffDisp(chemtbl[k].DiffCoe,
                chemtbl[k].DispCoe, rttbl->Cementation,
                elem[i].chms_fbrunsat.t_conc[k], elem[i].chms_fbrgw.t_conc[k],
                elem[i].geol.smcmax, 0.5 * elem[i].geol.depth,
                elem[i].topo.area, elem[i].wf.fbr_rechg * elem[i].topo.area);

# if defined(_TGM_)
            /* Fractured bedrock discharge to river.
             * Note that FBR discharge is always non-negative */
            elem[i].chmf.fbr_discharge[k] = elem[i].wf.fbr_discharge *
                elem[i].topo.area * elem[i].chms_fbrgw.t_conc[k];
# endif

            /* Element to element */
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] == 0)
                {
                    elem[i].chmf.fbr_unsatflux[j][k] = 0.0;
                    /* Diffusion and dispersion are ignored for boundary fluxes
                     */
                    elem[i].chmf.fbrflow[j][k] =
                        (elem[i].attrib.fbrbc_type[j] == 0) ?
                        0.0 : elem[i].wf.fbrflow[j] *
                        ((elem[i].wf.fbrflow[j] > 0.0) ?
                        elem[i].chms_fbrgw.t_conc[k] :
                        elem[i].fbr_bc.conc[j][k]);
                }
                else
                {
                    nabr = &elem[elem[i].nabr[j] - 1];

                    /* Unsaturated zone diffusion */
                    elem[i].chmf.fbr_unsatflux[j][k] =
                        AdvDiffDisp(chemtbl[k].DiffCoe, chemtbl[k].DispCoe,
                        rttbl->Cementation, elem[i].chms_fbrunsat.t_conc[k],
                        nabr->chms_fbrunsat.t_conc[k],
                        0.5 * elem[i].geol.smcmax + 0.5 * nabr->geol.smcmax,
                        elem[i].topo.nabrdist[j],
                        0.5 * MAX(elem[i].geol.depth - elem[i].ws.fbr_gw, 0.0) +
                        0.5 * MAX(nabr->geol.depth - nabr->ws.fbr_gw, 0.0),
                        0.0);

                    /* Groundwater advection, diffusion, and dispersion */
                    elem[i].chmf.fbrflow[j][k] =
                        AdvDiffDisp(chemtbl[k].DiffCoe, chemtbl[k].DispCoe,
                        rttbl->Cementation, elem[i].chms_fbrgw.t_conc[k],
                        nabr->chms_fbrgw.t_conc[k],
                        0.5 * elem[i].geol.smcmax + 0.5 * nabr->geol.smcmax,
                        elem[i].topo.nabrdist[j],
                        0.5 * MAX(elem[i].ws.fbr_gw, 0.0) +
                        0.5 * MAX(nabr->ws.fbr_gw, 0.0),
                        elem[i].wf.fbrflow[j]);
                }
            }   /* End of element to element */
#endif
        }   /* End of species loop */
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        elem_struct    *left;
        elem_struct    *right;
        int             j, k;

        for (k = 0; k < NumSpc; k++)
        {
            /* Downstream and upstream */
            if (river[i].down > 0)
            {
                down = &river[river[i].down - 1];

                /* Stream */
                river[i].chmf.flux[DOWN_CHANL2CHANL][k] =
                    river[i].wf.rivflow[DOWN_CHANL2CHANL] *
                    ((river[i].wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
                    river[i].chms_stream.t_conc[k] :
                    down->chms_stream.t_conc[k]);

                /* Bed */
                river[i].chmf.flux[DOWN_AQUIF2AQUIF][k] =
                    river[i].wf.rivflow[DOWN_AQUIF2AQUIF] *
                    ((river[i].wf.rivflow[DOWN_AQUIF2AQUIF] > 0.0) ?
                    river[i].chms_rivbed.t_conc[k] :
                    down->chms_rivbed.t_conc[k]);
            }
            else
            {
                river[i].chmf.flux[DOWN_CHANL2CHANL][k] =
                    river[i].wf.rivflow[DOWN_CHANL2CHANL] *
                    river[i].chms_stream.t_conc[k];

                river[i].chmf.flux[DOWN_AQUIF2AQUIF][k] = 0.0;
            }

            /* Left and right banks */
            left = &elem[river[i].leftele - 1];
            right = &elem[river[i].rightele - 1];

            if (river[i].leftele > 0)
            {
                river[i].chmf.flux[LEFT_SURF2CHANL][k] =
                    river[i].wf.rivflow[LEFT_SURF2CHANL] *
                    ((river[i].wf.rivflow[LEFT_SURF2CHANL] > 0.0) ?
                    river[i].chms_stream.t_conc[k] :
                    left->prcps.t_conc[k] * rttbl->Condensation);

                river[i].chmf.flux[LEFT_AQUIF2CHANL][k] =
                    river[i].wf.rivflow[LEFT_AQUIF2CHANL] *
                    ((river[i].wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
                    river[i].chms_stream.t_conc[k] :
                    left->chms_gw.t_conc[k]);

                river[i].chmf.flux[LEFT_AQUIF2AQUIF][k] =
                    river[i].wf.rivflow[LEFT_AQUIF2AQUIF] *
                    ((river[i].wf.rivflow[LEFT_AQUIF2AQUIF] > 0.0) ?
                    river[i].chms_rivbed.t_conc[k] :
                    left->chms_gw.t_conc[k]);

#if defined(_FBR_) && defined(_TGM_)
                river[i].chmf.flux[LEFT_FBR2CHANL][k] =
                    river[i].wf.rivflow[LEFT_FBR2CHANL] *
                    left->chms_fbrgw.t_conc[k];
#endif

                for (j = 0; j < NUM_EDGE; j++)
                {
                    if (left->nabr_river[j] == i + 1)
                    {
                        left->chmf.subflux[j][k] =
                            -(river[i].chmf.flux[LEFT_AQUIF2CHANL][k] +
                            river[i].chmf.flux[LEFT_AQUIF2AQUIF][k]);
                        break;
                    }
                }

            }

            if (river[i].rightele > 0)
            {
                river[i].chmf.flux[RIGHT_SURF2CHANL][k] =
                    river[i].wf.rivflow[RIGHT_SURF2CHANL] *
                    ((river[i].wf.rivflow[RIGHT_SURF2CHANL] > 0.0) ?
                    river[i].chms_stream.t_conc[k] :
                    right->prcps.t_conc[k] * rttbl->Condensation);

                river[i].chmf.flux[RIGHT_AQUIF2CHANL][k] =
                    river[i].wf.rivflow[RIGHT_AQUIF2CHANL] *
                    ((river[i].wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
                    river[i].chms_stream.t_conc[k] :
                    right->chms_gw.t_conc[k]);

                river[i].chmf.flux[RIGHT_AQUIF2AQUIF][k] =
                    river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] *
                    ((river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] > 0.0) ?
                    river[i].chms_rivbed.t_conc[k] :
                    right->chms_gw.t_conc[k]);

#if defined(_FBR_) && defined(_TGM_)
                river[i].chmf.flux[RIGHT_FBR2CHANL][k] =
                    river[i].wf.rivflow[RIGHT_FBR2CHANL] *
                    right->chms_fbrgw.t_conc[k];
#endif

                for (j = 0; j < NUM_EDGE; j++)
                {
                    if (right->nabr_river[j] == i + 1)
                    {
                        right->chmf.subflux[j][k] =
                            -(river[i].chmf.flux[RIGHT_AQUIF2CHANL][k] +
                            river[i].chmf.flux[RIGHT_AQUIF2AQUIF][k]);
                        break;
                    }
                }
            }

            river[i].chmf.flux[CHANL_LKG][k] = river[i].wf.rivflow[CHANL_LKG] *
                ((river[i].wf.rivflow[CHANL_LKG] > 0.0) ?
                river[i].chms_stream.t_conc[k] :
                river[i].chms_rivbed.t_conc[k]);
        }
    }

    /*
     * Accumulate to get in-flow for down segments
     */
    for (i = 0; i < nriver; i++)
    {
        int             k;
        river_struct   *down;

        for (k = 0; k < NumSpc; k++)
        {
            if (river[i].down > 0)
            {
                down = &river[river[i].down - 1];

                down->chmf.flux[UP_CHANL2CHANL][k] -=
                    river[i].chmf.flux[DOWN_CHANL2CHANL][k];

                down->chmf.flux[UP_AQUIF2AQUIF][k] -=
                    river[i].chmf.flux[DOWN_AQUIF2AQUIF][k];
            }
        }
    }
}

double AdvDiffDisp(double DiffCoe, double DispCoe, double cementation,
    double conc_up, double conc_down, double porosity, double distance,
    double area, double wflux)
{
    /*
     * Calculate the total of advection, diffusion and dispersion
     */
    double          inv_dist;
    double          diff_conc;
    double          diff_flux, disp_flux;

    inv_dist = 1.0 / distance;

    /* Difference in concentration (mol kg-1 water) */
    diff_conc = conc_up - conc_down;

    /* Diffusion flux, effective diffusion coefficient  */
    diff_flux = DiffCoe * area * pow(porosity, cementation) *
        inv_dist * diff_conc;

    /* Longitudinal dispersion */
    disp_flux = fabs(wflux) * DispCoe * inv_dist * diff_conc;

    return wflux * ((wflux > 0.0) ? conc_up : conc_down) +
        diff_flux + disp_flux;
}
