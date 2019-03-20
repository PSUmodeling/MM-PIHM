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

void OS3D(double stepsize, const chemtbl_struct chemtbl[],
    const rttbl_struct *rttbl, Chem_Data CD)
{
    int             i;

    /* Initalize the allocated array */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumVol; i++)
    {
        int             j;

        if (CD->Vcele[i].type != VIRTUAL_VOL)
        {
            for (j = 0; j < NumSpc; j++)
            {
                CD->Vcele[i].transp_flux[j] = 0.0;
            }
        }
    }

    for (i = 0; i < CD->NumFac; i++)
    {
        int             j;

        if (CD->Flux[i].BC != NO_FLOW)
        {
            for (j = 0; j < NumSpc; j++)
            {
                CD->Vcele[CD->Flux[i].nodeup - 1].transp_flux[j] +=
                    Dconc(&CD->Flux[i], CD->Vcele, chemtbl,
                        rttbl->Cementation, 0, j);
            }
        }
    }
}

double Dconc(const face *Flux, const vol_conc Vcele[], const chemtbl_struct chemtbl[],
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
        diff_flux = chemtbl[spc_ind].DiffCoe *
            pow(Vcele[node_1].porosity, cementation);
        /* Diffusion flux, effective diffusion coefficient  */
        if (velocity < 0.0)
        {
            disp_flux = velocity * chemtbl[spc_ind].DispCoe;
        }
        else
        {
            disp_flux = -velocity * chemtbl[spc_ind].DispCoe;
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
    //else if (TVDFlg == 1)
    //{
    //    if (flux_t > 0)
    //    {
    //        if (node_4 > 0)
    //        {
    //            r_ = (Vcele[node_2].t_conc[spc_ind] -
    //                Vcele[node_4].t_conc[spc_ind] + EPSILON) /
    //                (Vcele[node_1].t_conc[spc_ind] -
    //                Vcele[node_2].t_conc[spc_ind] + EPSILON);
    //            beta_ = max(0, min(min(2, 2 * r_), (2 + r_) / 3));
    //            temp_conc = Vcele[node_2].t_conc[spc_ind] + beta_ *
    //                (Vcele[node_1].t_conc[spc_ind] -
    //                Vcele[node_2].t_conc[spc_ind]) * 0.5;
    //        }
    //        else
    //        {
    //            temp_conc = Vcele[node_2].t_conc[spc_ind];
    //        }
    //    }
    //    else
    //    {
    //        if (node_3 > 0)
    //        {
    //            r_ = (Vcele[node_1].t_conc[spc_ind] -
    //                Vcele[node_3].t_conc[spc_ind] +
    //                EPSILON) / (Vcele[node_2].t_conc[spc_ind] -
    //                Vcele[node_1].t_conc[spc_ind] + EPSILON);
    //            beta_ = max(0, min(min(2, 2 * r_), (2 + r_) / 3));
    //            temp_conc =
    //                Vcele[node_1].t_conc[spc_ind] +
    //                beta_ * (Vcele[node_2].t_conc[spc_ind] -
    //                Vcele[node_1].t_conc[spc_ind]) * 0.5;
    //        }
    //        else
    //        {
    //            temp_conc = Vcele[node_1].t_conc[spc_ind];
    //        }
    //    }

    //    /* Add tributary */
    //    temp_conc_trib = (flux_t_trib > 0) ?
    //        Vcele[node_5_trib].t_conc[spc_ind] : 0.0;
    //}

    /* Advective flux */
    temp_dconc += temp_conc * flux_t + temp_conc_trib * flux_t_trib;

    if (Flux->BC == DISPERSION)
    {
        temp_dconc -= diff_flux + disp_flux;
    }

    return temp_dconc;
}
