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

void OS3D(realtype t, realtype stepsize, Chem_Data CD)
{
    /* Input t and stepsize in the unit of minute */
    double        **dconc = (double **)malloc(CD->NumOsv * sizeof(double *));
    int             i, j, k, abnormalflg;
    double          temp_dconc, diff_conc, unit_c, timelps, invavg, adpstep;
    double         *tmpconc = (double *)malloc(CD->NumSpc * sizeof(double));

    abnormalflg = 0;
    unit_c = 1.0 / 1440.0;

    /* Initalize the allocated array */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumOsv; i++)
    {
        int             j;
        dconc[i] = (double *)malloc(CD->NumSpc * sizeof(double));
        for (j = 0; j < CD->NumSpc; j++)
        {
            dconc[i][j] = 0.0;
        }
    }

    for (i = 0; i < CD->NumFac; i++)
    {
        CD->Flux[i].q = 0.0;

        for (j = 0; j < CD->NumSpc; j++)
        {
            temp_dconc = Dconc(&CD->Flux[i], CD->Vcele, CD->chemtype,
                CD->Cementation, CD->TVDFlg, j);

            CD->Flux[i].q = temp_dconc;
            dconc[CD->Flux[i].nodeup - 1][j] += temp_dconc;
        }
    }

    /* Local time step part */
    for (i = 0; i < CD->NumOsv; i++)
    {
        if (CD->CptFlg == 1)
        {
            if ((CD->Vcele[i].rt_step < stepsize) &&
                (CD->Vcele[i].height_t > 1.0E-3) &&
                (CD->Vcele[i].height_o > 1.0E-3))
            {
                /* Use its intrinsic smaller step for small/fast flowing cells
                 *  ~= slow cells (in term of time marching). */
                if (i < 2 * CD->NumEle + CD->NumRiv - CD->RivOff)
                {
                    timelps = t;
                    invavg = 1.0 / stepsize;
                    adpstep = CD->Vcele[i].rt_step;

                    while (timelps < t + stepsize)
                    {
                        adpstep = (adpstep > t + stepsize - timelps) ?
                            t + stepsize - timelps : adpstep;

                        diff_conc = 0.0;
                        for (j = 0; j < CD->NumSpc; j++)
                        {
                            tmpconc[j] = dconc[i][j] * adpstep +
                                CD->Vcele[i].t_conc[j] *
                                (CD->Vcele[i].porosity * 0.5 *
                                (CD->Vcele[i].vol_o + CD->Vcele[i].vol));

                            if (CD->PrpFlg)
                            {
                                if (CD->Vcele[i].q > 0.0)
                                {
                                    if (strcmp(CD->chemtype[j].ChemName, "'DOC'") == 0)
                                    {
                                        tmpconc[j] +=
                                            CD->Precipitation.t_conc[j] *
                                            CD->Vcele[i].q * adpstep * unit_c *
                                            CD->Condensation * CD->CalPrcpconc;
                                    }
                                    else
                                    {
                                        tmpconc[j] +=
                                            CD->Precipitation.t_conc[j] *
                                            CD->Vcele[i].q * adpstep * unit_c *
                                            CD->Condensation;
                                    }
                                }
                                tmpconc[j] +=
                                    CD->Precipitation.t_conc[j] *
                                    CD->Vcele[i].q * adpstep * unit_c *
                                    CD->Condensation;
                                if (CD->Vcele[i].q < 0.0)
                                {
                                    tmpconc[j] += 0.0;  /* n_0 design */
                                }
                            }
                            if ((tmpconc[j] < 0.0) &&
                                (strcmp(CD->chemtype[j].ChemName, "'H+'")))
                            {
                                fprintf(stderr, "Local time stepping check\n");
                                fprintf(stderr,
                                    "negative concentration change at species %s !\n",
                                    CD->chemtype[j].ChemName);
                                fprintf(stderr, "Change from fluxes: %8.4g\n",
                                    dconc[i][j] * adpstep);
                                fprintf(stderr,
                                    "Change from precipitation: %8.4g\n",
                                    CD->Precipitation.t_conc[j] *
                                    CD->Vcele[i].q * adpstep * unit_c *
                                    CD->Condensation);
                                fprintf(stderr, "Original mass: %8.4g\n",
                                    CD->Vcele[i].t_conc[j] *
                                    (CD->Vcele[i].porosity * 0.5 *
                                        (CD->Vcele[i].vol_o +
                                            CD->Vcele[i].vol)));
                                fprintf(stderr,
                                    "Old Conc: %8.4g\t New Conc: %8.4g\n",
                                    CD->Vcele[i].t_conc[j], tmpconc[j]);
                                ReportError(CD->Vcele[i], CD);
                                abnormalflg = 1;
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
                                if (CD->Flux[k].nodeup == CD->Vcele[i].index)
                                {
                                    for (j = 0; j < CD->NumSpc; j++)
                                    {
                                        temp_dconc = Dconc(&CD->Flux[k],
                                            CD->Vcele, CD->chemtype,
                                            CD->Cementation, CD->TVDFlg, j);

                                        tmpconc[j] += temp_dconc;
                                    }

                                }
                            }

                            for (j = 0; j < CD->NumSpc; j++)
                            {
                                dconc[i][j] = tmpconc[j];
                            }
                        }

                        timelps += adpstep;
                        if (timelps >= t + stepsize)
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
                if (CD->Vcele[i].rt_step < stepsize)
                {
                    continue;
                }
                if ((i >= 2 * CD->NumEle) && (i < 2 * CD->NumEle + CD->NumRiv))
                {
                    /* Treated in the above section */
                    continue;
                }
            }

            for (j = 0; j < CD->NumSpc; j++)
            {
                tmpconc[j] =
                    dconc[i][j] * stepsize +
                    CD->Vcele[i].t_conc[j] * (CD->Vcele[i].porosity *
                    CD->Vcele[i].vol_o);
                /* Need consider the change of concentration at the unsat
                 * zone from precipitation */
                if (CD->PrpFlg)
                {
                    if (CD->Vcele[i].q > 0.0)
                    {
                        if (strcmp(CD->chemtype[j].ChemName, "'DOC'") == 0)
                        {
                            tmpconc[j] += CD->Precipitation.t_conc[j] *
                                CD->Vcele[i].q * adpstep * unit_c *
                                CD->Condensation * CD->CalPrcpconc;
                        }
                        else
                        {
                            tmpconc[j] += CD->Precipitation.t_conc[j] *
                                CD->Vcele[i].q * adpstep * unit_c *
                                CD->Condensation;
                        }
                    }
                    if (CD->Vcele[i].q < 0.0)
                    {
                        tmpconc[j] += CD->Vcele[i].t_conc[j] *
                            CD->Vcele[i].q * stepsize * unit_c;
                    }
                }
                tmpconc[j] = tmpconc[j] /
                    (CD->Vcele[i].porosity * CD->Vcele[i].vol);
            }
            if (CD->Vcele[i].illness < 20)
                for (j = 0; j < CD->NumSpc; j++)
                {

                    if ((tmpconc[j] < 0.0) &&
                        (strcmp(CD->chemtype[j].ChemName, "'H+'")) &&
                        (i < CD->NumEle * 2))
                    {
                        fprintf(stderr,
                            "negative concentration change at species %s !\n",
                            CD->chemtype[j].ChemName);
                        fprintf(stderr, "Change from fluxes: %8.4g\t",
                            dconc[i][j] * stepsize);
                        fprintf(stderr, "Original mass: %8.4g\n",
                            CD->Vcele[i].t_conc[j] *
                            (CD->Vcele[i].porosity * CD->Vcele[i].vol_o));
                        fprintf(stderr,
                            "New mass: %8.4g\t New Volume: %8.4g\t Old Conc: %8.4g\t New Conc: %8.4g\t Timestep: %8.4g\n",
                            dconc[i][j] * stepsize +
                            CD->Vcele[i].t_conc[j] *
                            (CD->Vcele[i].porosity * CD->Vcele[i].vol_o),
                            CD->Vcele[i].porosity * CD->Vcele[i].height_t *
                            CD->Vcele[i].area, CD->Vcele[i].t_conc[j],
                            tmpconc[j], CD->Vcele[i].rt_step);
                        ReportError(CD->Vcele[i], CD);
                        CD->Vcele[i].illness++;
                        abnormalflg = 1;
                    }
                    CD->Vcele[i].t_conc[j] = tmpconc[j];
                }
        }
    }

    for (i = 0; i < CD->NumOsv; i++)
    {
        free(dconc[i]);
    }
    free(dconc);
    free(tmpconc);
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
        diff_flux = -diff_flux * inv_dist * diff_conc * area;
        /* Diffusion is in the opposite direction of conc gradient */
        disp_flux = disp_flux * inv_dist * diff_conc * area;
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
