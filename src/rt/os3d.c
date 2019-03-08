/******************************************************************************
* This subroutine is used to calculate the advection diffusion dispersion
* It uses a similar OS3D scheme as detailed in Crunchflow user's manual.
*
* If you have any questions, concerns, suggestions, please contact me at
* the following address:
*     Developer: Chen Bao <baochen.d.s@gmail.com>
*     Version  : pre-alpha
*     Date     : June, 2013
*****************************************************************************/
#include "pihm.h"

#define max(a,b) ((a)>(b) ? (a):(b))
#define min(a,b) ((a)<(b) ? (a):(b))
#define EPSILON 1.0E-20

void OS3D(double stepsize, Chem_Data CD)
{
    /* Input stepsize in the unit of second */
    double        **dconc = (double **)malloc(CD->NumOsv * sizeof(double *));
    int             i;
    double          unit_c;
    double        **tconc;

    unit_c = 1.0 / 1440.0;

    tconc = (double **)malloc(CD->NumOsv * sizeof(double *));

    /* Initalize the allocated array */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumOsv; i++)
    {
        int             j;

        tconc[i] = (double *)malloc(CD->NumSpc * sizeof(double));
        dconc[i] = (double *)malloc(CD->NumSpc * sizeof(double));

        for (j = 0; j < CD->NumSpc; j++)
        {
            tconc[i][j] = CD->Vcele[i].t_conc[j];
            dconc[i][j] = 0.0;
        }
    }

    for (i = 0; i < CD->NumFac; i++)
    {
        int             j;

        if (CD->Flux[i].BC != NO_FLOW)
        {
            for (j = 0; j < CD->NumSpc; j++)
            {
                dconc[CD->Flux[i].nodeup - 1][j] +=
                    Dconc(&CD->Flux[i], CD->Vcele, CD->chemtype,
                        CD->Cementation, CD->TVDFlg, j);
            }
        }
    }

    /* Local time step part */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumOsv; i++)
    {
        int             j, k;
        double          diff_conc;
        double          timelps;
        double          adpstep;
        double         *tmpconc = (double *)malloc(CD->NumSpc * sizeof(double));

        adpstep = CD->Vcele[i].rt_step;

        if (CD->CptFlg == 1)
        {
            if ((CD->Vcele[i].rt_step < stepsize / 60) &&
                (CD->Vcele[i].height_t > 1.0E-3) &&
                (CD->Vcele[i].height_o > 1.0E-3))
            {
                /* Use its intrinsic smaller step for small/fast flowing cells
                 *  ~= slow cells (in term of time marching). */
                if (CD->Vcele[i].type == UNSAT_VOL ||
                    CD->Vcele[i].type == GW_VOL)
                {
                    timelps = 0.0;

                    while (timelps < stepsize / 60)
                    {
                        adpstep = (adpstep > stepsize / 60 - timelps) ?
                            stepsize / 60 - timelps : adpstep;

                        diff_conc = 0.0;
                        for (j = 0; j < CD->NumSpc; j++)
                        {
                            tmpconc[j] = dconc[i][j] * adpstep +
                                CD->Vcele[i].t_conc[j] *
                                (CD->Vcele[i].porosity * 0.5 *
                                (CD->Vcele[i].vol_o + CD->Vcele[i].vol));

                            if ((tmpconc[j] < 0.0) &&
                                (strcmp(CD->chemtype[j].ChemName, "'H+'")))
                            {
                                fprintf(stderr, "Local time stepping check\n");
                                fprintf(stderr,
                                    "negative concentration change at species %s !\n",
                                    CD->chemtype[j].ChemName);
                                fprintf(stderr, "Change from fluxes: %8.4g\n",
                                    dconc[i][j] * adpstep);
                                fprintf(stderr, "Original mass: %8.4g\n",
                                    CD->Vcele[i].t_conc[j] *
                                    (CD->Vcele[i].porosity * 0.5 *
                                        (CD->Vcele[i].vol_o +
                                            CD->Vcele[i].vol)));
                                fprintf(stderr,
                                    "Old Conc: %8.4g\t New Conc: %8.4g\n",
                                    CD->Vcele[i].t_conc[j], tmpconc[j]);
                                ReportError(CD->Vcele[i], CD);
                                CD->Vcele[i].illness++;
                            }

                            tmpconc[j] =
                                tmpconc[j] / (CD->Vcele[i].porosity * 0.5 *
                                (CD->Vcele[i].vol + CD->Vcele[i].vol_o));

                            if (CD->Vcele[i].illness < 20)
                            {
                                diff_conc =
                                    max(fabs(CD->Vcele[i].t_conc[j] -
                                    tmpconc[j]), diff_conc);
                                CD->Vcele[i].t_conc[j] = tmpconc[j];
                            }
                        }

                        if (diff_conc > 1.0E-6) /* Which lead to the change in
                                                 * the flux of mass between
                                                 * cells */
                        {
                            for (j = 0; j < CD->NumSpc; j++)
                                tmpconc[j] = 0.0;

                            for (k = 0; k < CD->NumFac; k++)
                            {
                                if (CD->Flux[k].BC != NO_FLOW)
                                {
                                    if (CD->Flux[k].nodeup == CD->Vcele[i].index)
                                    {
                                        for (j = 0; j < CD->NumSpc; j++)
                                        {
                                            tmpconc[j] += Dconc(&CD->Flux[k],
                                                CD->Vcele, CD->chemtype,
                                                CD->Cementation, CD->TVDFlg, j);
                                        }
                                    }
                                }
                            }

                            for (j = 0; j < CD->NumSpc; j++)
                            {
                                dconc[i][j] = tmpconc[j];
                            }
                        }

                        timelps += adpstep;
                        if (timelps >= stepsize / 60)
                            break;
                    }
                }
            }
        }
        /* CD->Vcele[i].rt_step >= stepsize */

        /* For blocks with very small content, we just skip it */
        if ((CD->Vcele[i].height_t > 1.0E-3) &&
            (CD->Vcele[i].height_o > 1.0E-3))
        {
            if (CD->CptFlg)
            {
                if (CD->Vcele[i].rt_step < stepsize / 60)
                {
                    continue;
                }
            }

            for (j = 0; j < CD->NumSpc; j++)
            {
                tmpconc[j] =
                    dconc[i][j] * stepsize / 60 +
                    CD->Vcele[i].t_conc[j] * (CD->Vcele[i].porosity *
                    CD->Vcele[i].vol_o);

                tmpconc[j] = tmpconc[j] /
                    (CD->Vcele[i].porosity * CD->Vcele[i].vol);
            }
            if (CD->Vcele[i].illness < 20)
                for (j = 0; j < CD->NumSpc; j++)
                {

                    if ((tmpconc[j] < 0.0) &&
                        (strcmp(CD->chemtype[j].ChemName, "'H+'")) &&
                        (CD->Vcele[i].type == GW_VOL ||
                        CD->Vcele[i].type == UNSAT_VOL))
                    {
                        fprintf(stderr,
                            "negative concentration change at species %s !\n",
                            CD->chemtype[j].ChemName);
                        fprintf(stderr, "Change from fluxes: %8.4g\t",
                            dconc[i][j] * stepsize / 60);
                        fprintf(stderr, "Original mass: %8.4g\n",
                            CD->Vcele[i].t_conc[j] *
                            (CD->Vcele[i].porosity * CD->Vcele[i].vol_o));
                        fprintf(stderr,
                            "New mass: %8.4g\t New Volume: %8.4g\t Old Conc: %8.4g\t New Conc: %8.4g\t Timestep: %8.4g\n",
                            dconc[i][j] * stepsize / 60 +
                            CD->Vcele[i].t_conc[j] *
                            (CD->Vcele[i].porosity * CD->Vcele[i].vol_o),
                            CD->Vcele[i].porosity * CD->Vcele[i].height_t *
                            CD->Vcele[i].area, CD->Vcele[i].t_conc[j],
                            tmpconc[j], CD->Vcele[i].rt_step);
                        ReportError(CD->Vcele[i], CD);
                        CD->Vcele[i].illness++;
                    }
                    tconc[i][j] = tmpconc[j];
                }
        }

        free(tmpconc);
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumOsv; i++)
    {
        int             j;

        for (j = 0; j < CD->NumSpc; j++)
        {
            CD->Vcele[i].t_conc[j] = tconc[i][j];
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumOsv; i++)
    {
        free(tconc[i]);
        free(dconc[i]);
    }
    free(tconc);
    free(dconc);
}

double Dconc(const face *Flux, const vol_conc Vcele[], const species chemtype[],
    double cementation, int TVDFlg, int spc_ind)
{
    int             node_1, node_2, node_3, node_4;
    int             node_5_trib;
    double          flux_t, flux_t_trib;
    double          distance;
    double          velocity;
    double          area;
    double          inv_dist;
    double          diff_conc;
    double          diff_flux, disp_flux;
    double          temp_dconc;
    double          temp_dconc_trib;
    double          temp_conc;
    double          temp_conc_trib;
    double          r_, beta_;
    double          unit_c = 1.0 / 1440;

    node_1 = Flux->nodeup - 1;
    node_2 = Flux->nodelo - 1;
    node_3 = Flux->nodeuu - 1;
    node_4 = Flux->nodell - 1;
    node_5_trib = Flux->node_trib - 1;
    flux_t = -Flux->flux;
    flux_t_trib = -Flux->flux_trib;
    distance = Flux->distance;
    velocity = -Flux->velocity;
    area = Flux->s_area;
    inv_dist = 1.0 / distance;
    /* Difference in concentration, * in M/kg w */
    diff_conc = Vcele[node_2].t_conc[spc_ind] - Vcele[node_1].t_conc[spc_ind];
    diff_flux = 0.0;
    disp_flux = 0.0;

    if (Flux->BC == DISPERSION)    /* Only for soil flow (not river flow)*/
    {
        diff_flux = chemtype[spc_ind].DiffCoe *
            pow(Vcele[node_1].porosity, cementation);
        /* Diffusion flux, effective diffusion coefficient  */
        if (velocity < 0.0)
        {
            disp_flux = velocity * chemtype[spc_ind].DispCoe;
        }
        else
        {
            disp_flux = -velocity * chemtype[spc_ind].DispCoe;
        }
        /* Longitudinal dispersion */
        diff_flux = -diff_flux * 86400 * inv_dist * diff_conc * area;
        /* Diffusion is in the opposite direction of conc gradient */
        disp_flux = disp_flux * 86400 * inv_dist * diff_conc * area;
    }

    /* Use temp_conc to store the concentration at the surfaces
     * Use temp_dconc to store the concentration changes at the cell */
    temp_dconc = 0.0;
    temp_dconc_trib = 0.0;
    temp_conc = 0.0;
    temp_conc_trib = 0.0;

    if (TVDFlg == 0)
    {
        temp_conc = (flux_t > 0) ?
            Vcele[node_2].t_conc[spc_ind] : Vcele[node_1].t_conc[spc_ind];

        /* Add tributary */
        temp_conc_trib = (flux_t_trib > 0) ?
            Vcele[node_5_trib].t_conc[spc_ind] : 0.0;
    }
    else if (TVDFlg == 1)
    {
        if (flux_t > 0)
        {
            if (node_4 > 0)
            {
                r_ = (Vcele[node_2].t_conc[spc_ind] -
                    Vcele[node_4].t_conc[spc_ind] + EPSILON) /
                    (Vcele[node_1].t_conc[spc_ind] -
                    Vcele[node_2].t_conc[spc_ind] + EPSILON);
                beta_ = max(0, min(min(2, 2 * r_), (2 + r_) / 3));
                temp_conc = Vcele[node_2].t_conc[spc_ind] + beta_ *
                    (Vcele[node_1].t_conc[spc_ind] -
                    Vcele[node_2].t_conc[spc_ind]) * 0.5;
            }
            else
            {
                temp_conc = Vcele[node_2].t_conc[spc_ind];
            }
        }
        else
        {
            if (node_3 > 0)
            {
                r_ = (Vcele[node_1].t_conc[spc_ind] -
                    Vcele[node_3].t_conc[spc_ind] +
                    EPSILON) / (Vcele[node_2].t_conc[spc_ind] -
                    Vcele[node_1].t_conc[spc_ind] + EPSILON);
                beta_ = max(0, min(min(2, 2 * r_), (2 + r_) / 3));
                temp_conc =
                    Vcele[node_1].t_conc[spc_ind] +
                    beta_ * (Vcele[node_2].t_conc[spc_ind] -
                    Vcele[node_1].t_conc[spc_ind]) * 0.5;
            }
            else
            {
                temp_conc = Vcele[node_1].t_conc[spc_ind];
            }
        }

        /* Add tributary */
        temp_conc_trib = (flux_t_trib > 0) ?
            Vcele[node_5_trib].t_conc[spc_ind] : 0.0;
    }
    else
    {
        temp_conc = (flux_t > 0) ?
            Vcele[node_2].t_conc[spc_ind] : Vcele[node_1].t_conc[spc_ind];

        /* Add tributary */
        temp_conc_trib = (flux_t_trib > 0) ?
            Vcele[node_5_trib].t_conc[spc_ind] : 0.0;
    }

    /* Advective flux */
    temp_dconc += temp_conc * flux_t + temp_conc_trib * flux_t_trib;

    if (Flux->BC == DISPERSION)
    {
        temp_dconc -= diff_flux + disp_flux;
    }

    temp_dconc *= unit_c;

    return temp_dconc;
}
