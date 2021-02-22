#include "pihm.h"

void SoluteTranspt(double diff_coef, double disp_coef, double cementation, elem_struct elem[], river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k;

        for (k = 0; k < nsolute; k++)
        {
            // Initialize chemical fluxes
            elem[i].solute[k].infil = 0.0;

            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].subflux[j] = 0.0;
            }
        }

#if defined(_DGW_)
        for (k = 0; k < nsolute; k++)
        {
            // Initialize chemical fluxes
            elem[i].solute[k].infil_geol = 0.0;
#if defined(_LUMPED_)
            elem[i].solute[k].dgw_leach = 0.0;
#endif

            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].dgwflux[j] = 0.0;
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

        for (k = 0; k < nsolute; k++)
        {
            // Initialize chemical fluxes
            for (j = 0; j < NUM_RIVFLX; j++)
            {
                river[i].solute[k].flux[j] = 0.0;
            }
        }
    }

    // Calculate chemical fluxes
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k;
        elem_struct    *nabr;

        for (k = 0; k < nsolute; k++)
        {
            // Infiltration
            elem[i].solute[k].infil = elem[i].wf.infil * ((elem[i].wf.infil > 0.0) ? elem[i].solute[k].conc_surf : 0.0);

            // Element to element
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

                        wflux = elem[i].wf.subsurf[j] + ((elem[i].ind == river_ptr->left) ?
                            river_ptr->wf.rivflow[AQUIFER_LEFT] : river_ptr->wf.rivflow[AQUIFER_RIGHT]);
                    }

                    // Advection, diffusion, and dispersion between triangular elements
                    elem[i].solute[k].subflux[j] = AdvDiffDisp(diff_coef, disp_coef, cementation,
                        elem[i].solute[k].conc, nabr->solute[k].conc, 0.5 * (elem[i].soil.smcmax + nabr->soil.smcmax),
                        elem[i].topo.dist_nabr[j], 0.5 * (elem[i].soil.depth + nabr->soil.depth), wflux);
                }
            }   // End of element to element

#if defined(_DGW_)
            // Bedrock infiltration
            elem[i].solute[k].infil_geol = elem[i].wf.infil_geol * ((elem[i].wf.infil_geol > 0.0) ?
                elem[i].solute[k].conc : elem[i].solute[k].conc_geol);

# if defined(_LUMPED_)
            // Deep zone leaching to river. Note that deep zone runoff is always non-negative
            elem[i].solute[k].dgw_leach = elem[i].wf.dgw_runoff * elem[i].solute[k].conc_geol;
# endif

            // Element to element
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] == 0)
                {
                    // Diffusion and dispersion are ignored for boundary fluxes
                    elem[i].solute[k].dgwflux[j] = (elem[i].attrib.bc_geol[j] == 0) ?
                        0.0 : elem[i].wf.dgw[j] * ((elem[i].wf.dgw[j] > 0.0) ?
                        elem[i].solute[k].conc_geol : elem[i].bc_geol.conc[j][k]);
                }
                else
                {
                    nabr = &elem[elem[i].nabr[j] - 1];

                    // Groundwater advection, diffusion, and dispersion
                    elem[i].solute[k].dgwflux[j] = AdvDiffDisp(diff_coef, disp_coef, cementation,
                        elem[i].solute[k].conc_geol, nabr->solute[k].conc_geol,
                        0.5 * (elem[i].geol.smcmax + nabr->geol.smcmax), elem[i].topo.dist_nabr[j],
                        0.5 * (elem[i].geol.depth + nabr->geol.depth), elem[i].wf.dgw[j]);
                }
            }   // End of element to element
#endif
        }   // End of species loop
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            // Downstream and upstream
            if (river[i].down > 0)
            {
                down = &river[river[i].down - 1];

                // Stream
                river[i].solute[k].flux[DOWNSTREAM] = river[i].wf.rivflow[DOWNSTREAM] *
                    ((river[i].wf.rivflow[DOWNSTREAM] > 0.0) ? river[i].solute[k].conc : down->solute[k].conc);
            }
            else
            {
                river[i].solute[k].flux[DOWNSTREAM] = river[i].wf.rivflow[DOWNSTREAM] * river[i].solute[k].conc;
            }
        }

            // Left and right banks
        if (river[i].left > 0)
        {
#if defined(_DGW_) && defined(_LUMPED_)
            RiverElemSoluteFlow(SURF_LEFT, AQUIFER_LEFT, DGW_LEFT, &elem[river[i].left - 1], &river[i]);
#else
            RiverElemSoluteFlow(SURF_LEFT, AQUIFER_LEFT, &elem[river[i].left - 1], &river[i]);
#endif
        }

        if (river[i].right > 0)
        {
#if defined(_DGW_) && defined(_LUMPED_)
            RiverElemSoluteFlow(SURF_RIGHT, AQUIFER_RIGHT, DGW_RIGHT, &elem[river[i].right - 1], &river[i]);
#else
            RiverElemSoluteFlow(SURF_RIGHT, AQUIFER_RIGHT, &elem[river[i].right - 1], &river[i]);
#endif
        }
    }

    // Accumulate to get in-flow for down segments
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            if (river[i].down > 0)
            {
                river[river[i].down - 1].solute[k].flux[UPSTREAM] -= river[i].solute[k].flux[DOWNSTREAM];
            }
        }
    }
}

#if defined(_DGW_) && defined(_LUMPED_)
void RiverElemSoluteFlow(int surf_to_chanl, int aquif_to_chanl, int dgw_to_chanl, elem_struct *bank,
    river_struct *river)
#else
void RiverElemSoluteFlow(int surf_to_chanl, int aquif_to_chanl, elem_struct *bank, river_struct *river)
#endif
{
    int             j, k;

    for (k = 0; k < nsolute; k++)
    {
        river->solute[k].flux[surf_to_chanl] = river->wf.rivflow[surf_to_chanl] *
            ((river->wf.rivflow[surf_to_chanl] > 0.0) ? river->solute[k].conc : bank->solute[k].conc_surf);

        river->solute[k].flux[aquif_to_chanl] = river->wf.rivflow[aquif_to_chanl] *
            ((river->wf.rivflow[aquif_to_chanl] > 0.0) ? river->solute[k].conc : bank->solute[k].conc);

#if defined(_DGW_) && defined(_LUMPED_)
        river->solute[k].flux[dgw_to_chanl] = river->wf.rivflow[dgw_to_chanl] * bank->solute[k].conc_geol;
#endif

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (bank->nabr_river[j] == river->ind)
            {
                bank->solute[k].subflux[j] -= river->solute[k].flux[aquif_to_chanl];
                break;
            }
        }
    }
}

double AdvDiffDisp(double diff_coef, double disp_coef, double cementation, double conc_up, double conc_down,
    double porosity, double distance, double area, double wflux)
{
    // Calculate the total of advection, diffusion and dispersion
    double          inv_dist;
    double          diff_conc;
    double          diff_flux, disp_flux;

    inv_dist = 1.0 / distance;

    // Difference in concentration (mol kg-1 water)
    diff_conc = conc_up - conc_down;

    // Diffusion flux, effective diffusion coefficient
    diff_flux = diff_coef * area * pow(porosity, cementation) * inv_dist * diff_conc;

    // Longitudinal dispersion
    disp_flux = fabs(wflux) * disp_coef * inv_dist * diff_conc;

    return wflux * ((wflux > 0.0) ? conc_up : conc_down) + diff_flux + disp_flux;
}
