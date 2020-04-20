#include "pihm.h"

void SoluteTransp(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k, kk;
        double          storage;

        for (k = 0; k < NumSpc; k++)
        {
            /* Initialize chemical fluxes */
            elem[i].solute[k].infil = 0.0;

            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].subflux[j] = 0.0;
            }
        }

#if defined(_FBR_)
        for (k = 0; k < NumSpc; k++)
        {
            /* Initialize chemical fluxes */
            elem[i].solute[k].fbr_infil = 0.0;

            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].fbrflow[j] = 0.0;
            }
        }
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j, k;

        for (k = 0; k < NumSpc; k++)
        {
            /* Initialize chemical fluxes */
            for (j = 0; j < NUM_RIVFLX; j++)
            {
                river[i].solute[k].flux[j] = 0.0;
            }
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
            elem[i].solute[k].infil = elem[i].wf.infil * elem[i].topo.area *
                ((elem[i].wf.infil > 0.0) ?
                elem[i].prcps.t_conc[k] * rttbl->Condensation : 0.0);

            /* Element to element */
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] == 0)
                {
                    elem[i].solute[k].subflux[j] = 0.0;
                }
                else
                {
                    double          wflux;

                    nabr = &elem[elem[i].nabr[j] - 1];

                    if (elem[i].nabr_river[j] == 0)
                    {
                        wflux = elem[i].wf.subsurf[j];
                    }
                    else
                    {
                        river_struct           *river_ptr;

                        river_ptr = &river[elem[i].nabr_river[j] - 1];

                        wflux = elem[i].wf.subsurf[j] +
                            ((elem[i].ind == river_ptr->leftele) ?
                            river_ptr->wf.rivflow[LEFT_AQUIF2CHANL] :
                            river_ptr->wf.rivflow[RIGHT_AQUIF2CHANL]);
                    }

                    /* Advection, diffusion, and dispersion between triangular
                     * elements */
                    elem[i].solute[k].subflux[j] = AdvDiffDisp(chemtbl[k].DiffCoe,
                        chemtbl[k].DispCoe, rttbl->Cementation,
                        elem[i].chms.t_conc[k], nabr->chms.t_conc[k],
                        0.5 * (elem[i].soil.smcmax + nabr->soil.smcmax),
                        elem[i].topo.nabrdist[j],
                        0.5 * (elem[i].soil.depth + nabr->soil.depth), wflux);
                }
            }   /* End of element to element */

#if defined(_FBR_)
            /* Bedrock infiltration */
            elem[i].solute[k].fbr_infil = elem[i].wf.fbr_infil *
                elem[i].topo.area * ((elem[i].wf.fbr_infil > 0.0) ?
                elem[i].chms.t_conc[k] : elem[i].chms_geol.t_conc[k]);

# if defined(_TGM_)
            /* Fractured bedrock discharge to river.
             * Note that FBR discharge is always non-negative */
            elem[i].solute[k].fbr_discharge = elem[i].wf.fbr_discharge *
                elem[i].topo.area * elem[i].chms_fbrgw.t_conc[k];
# endif

            /* Element to element */
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] == 0)
                {
                    /* Diffusion and dispersion are ignored for boundary fluxes
                     */
                    elem[i].solute[k].fbrflow[j] =
                        (elem[i].attrib.fbrbc_type[j] == 0) ?
                        0.0 : elem[i].wf.fbrflow[j] *
                        ((elem[i].wf.fbrflow[j] > 0.0) ?
                        elem[i].chms_geol.t_conc[k] :
                        elem[i].fbr_bc.conc[j][k]);
                }
                else
                {
                    nabr = &elem[elem[i].nabr[j] - 1];

                    /* Groundwater advection, diffusion, and dispersion */
                    elem[i].solute[k].fbrflow[j] =
                        AdvDiffDisp(chemtbl[k].DiffCoe, chemtbl[k].DispCoe,
                        rttbl->Cementation, elem[i].chms_geol.t_conc[k],
                        nabr->chms_geol.t_conc[k],
                        0.5 * (elem[i].geol.smcmax + nabr->geol.smcmax),
                        elem[i].topo.nabrdist[j],
                        0.5 * (elem[i].geol.depth + nabr->geol.depth),
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
                river[i].solute[k].flux[DOWN_CHANL2CHANL] =
                    river[i].wf.rivflow[DOWN_CHANL2CHANL] *
                    ((river[i].wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
                    river[i].chms.t_conc[k] : down->chms.t_conc[k]);
            }
            else
            {
                river[i].solute[k].flux[DOWN_CHANL2CHANL] =
                    river[i].wf.rivflow[DOWN_CHANL2CHANL] *
                    river[i].chms.t_conc[k];
            }

            /* Left and right banks */
            left = &elem[river[i].leftele - 1];
            right = &elem[river[i].rightele - 1];

            if (river[i].leftele > 0)
            {
                river[i].solute[k].flux[LEFT_SURF2CHANL] =
                    river[i].wf.rivflow[LEFT_SURF2CHANL] *
                    ((river[i].wf.rivflow[LEFT_SURF2CHANL] > 0.0) ?
                    river[i].chms.t_conc[k] :
                    left->prcps.t_conc[k] * rttbl->Condensation);

                river[i].solute[k].flux[LEFT_AQUIF2CHANL] =
                    river[i].wf.rivflow[LEFT_AQUIF2CHANL] *
                    ((river[i].wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
                    river[i].chms.t_conc[k] : left->chms.t_conc[k]);

#if defined(_FBR_) && defined(_TGM_)
                river[i].solute[k].flux[LEFT_FBR2CHANL] =
                    river[i].wf.rivflow[LEFT_FBR2CHANL] *
                    left->chms_fbrgw.t_conc[k];
#endif

                for (j = 0; j < NUM_EDGE; j++)
                {
                    if (left->nabr_river[j] == i + 1)
                    {
                        left->solute[k].subflux[j] -=
                            river[i].solute[k].flux[LEFT_AQUIF2CHANL];
                        break;
                    }
                }

            }

            if (river[i].rightele > 0)
            {
                river[i].solute[k].flux[RIGHT_SURF2CHANL] =
                    river[i].wf.rivflow[RIGHT_SURF2CHANL] *
                    ((river[i].wf.rivflow[RIGHT_SURF2CHANL] > 0.0) ?
                    river[i].chms.t_conc[k] :
                    right->prcps.t_conc[k] * rttbl->Condensation);

                river[i].solute[k].flux[RIGHT_AQUIF2CHANL] =
                    river[i].wf.rivflow[RIGHT_AQUIF2CHANL] *
                    ((river[i].wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
                    river[i].chms.t_conc[k] :
                    right->chms.t_conc[k]);

#if defined(_FBR_) && defined(_TGM_)
                river[i].solute[k].flux[RIGHT_FBR2CHANL] =
                    river[i].wf.rivflow[RIGHT_FBR2CHANL] *
                    right->chms_fbrgw.t_conc[k];
#endif

                for (j = 0; j < NUM_EDGE; j++)
                {
                    if (right->nabr_river[j] == i + 1)
                    {
                        right->solute[k].subflux[j] -=
                            river[i].solute[k].flux[RIGHT_AQUIF2CHANL];
                        break;
                    }
                }
            }
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

                down->solute[k].flux[UP_CHANL2CHANL] -=
                    river[i].solute[k].flux[DOWN_CHANL2CHANL];
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
