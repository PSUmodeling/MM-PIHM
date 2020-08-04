#include "pihm.h"

#define TOL        1E-7
#define SKIP_JACOB 1

void Reaction(double stepsize, const chemtbl_struct chemtbl[],
    const kintbl_struct kintbl[], const rttbl_struct *rttbl, elem_struct elem[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          storage;
        double          satn;
        double          ftemp;
        int             k;

        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;

        satn = (elem[i].ws.unsat + elem[i].ws.gw) / elem[i].soil.depth;
        satn = MAX(satn, SATMIN);
        satn = MIN(satn, 1.0);

        for (k = 0; k < rttbl->num_spc; k++)
        {
            elem[i].chmf.react[k] = 0.0;
        }

        if (satn > 1.0E-2)
        {
            int             klayer;
            double          avg_stc = 0.0;

            for (klayer = 0; klayer < elem[i].ps.nlayers; klayer++)
            {
                avg_stc += elem[i].es.stc[klayer] * elem[i].ps.soil_depth[klayer];
            }
            avg_stc /= elem[i].soil.depth;

            ftemp = SoilTempFactor(avg_stc);

            ReactControl(chemtbl, kintbl, rttbl, stepsize, satn, ftemp,
                &elem[i].chms, elem[i].chmf.react);
        }

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].snksrc = elem[i].chmf.react[k] * storage;
        }

#if defined(_DGW_)
        storage = (elem[i].ws.unsat_geol + elem[i].ws.gw_geol) *
            elem[i].geol.porosity + elem[i].geol.depth * elem[i].geol.smcmin;

        satn = (elem[i].ws.unsat_geol + elem[i].ws.gw_geol) / elem[i].geol.depth;
        satn = MAX(satn, SATMIN);
        satn = MIN(satn, 1.0);

        for (k = 0; k < rttbl->num_spc; k++)
        {
            elem[i].chmf.react_geol[k] = 0.0;
        }

        /*
         * Deep zone
         */
        if (satn > 1.0E-2)
        {
            ftemp = SoilTempFactor(elem[i].ps.tbot);
            ReactControl(chemtbl, kintbl, rttbl, stepsize, satn, ftemp,
                &elem[i].chms_geol, elem[i].chmf.react_geol);
        }

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].snksrc_geol = elem[i].chmf.react_geol[k] *
                storage;
        }
#endif
    }
}

int _React(double stepsize, const chemtbl_struct chemtbl[],
    const kintbl_struct kintbl[], const rttbl_struct *rttbl, double satn,
    double ftemp, chmstate_struct *chms)
{
    int             i, j, k;
    int             control;
    int             num_spe = rttbl->num_stc + rttbl->num_ssc;
    int             min_pos;
    int             pivot_flg;
    int             mn, in;
    double          monodterm = 1.0, inhibterm = 1.0;
    double          residue[MAXSPS];
    double          residue_t[MAXSPS];
    double          tmpconc[MAXSPS];
    double          tot_conc[MAXSPS];
    double          area[MAXSPS];
    double          error[MAXSPS];
    double          gamma[MAXSPS];
    double          Keq[MAXSPS];
    double          Rate_pre[MAXSPS];
    double          IAP[MAXSPS];
    double          dependency[MAXSPS];
    double          Rate_spe[MAXSPS];
    double          Rate_spet[MAXSPS];
    const int       SUFEFF = 1;
    double          tmpval;
    double          tmpprb;
    double          inv_sat;
    double          I;
    double          Iroot;
    double          tmpKeq;
    double          adh;
    double          bdh;
    double          bdt;
    double          maxerror = 1.0;
    double          surf_ratio;
    double          tot_cec;
    double          tmpprb_inv;
    realtype      **jcb;
    sunindextype    p[MAXSPS];
    realtype        x_[MAXSPS];

    control = 0;
    tmpprb = 1.0E-2;
    tmpprb_inv = 1.0 / tmpprb;
    inv_sat = 1.0 / satn;

    for (i = 0; i < rttbl->num_min; i++)
    {
        area[i] = chms->ssa[i + rttbl->num_stc - rttbl->num_min] *
            chms->prim_conc[i + rttbl->num_stc - rttbl->num_min] *
            chemtbl[i + rttbl->num_stc - rttbl->num_min].molar_mass;
    }

    if (SUFEFF)
    {
        if (satn < 1.0)
        {
            //surf_ratio = 1.0;  /* # 1 function */
            surf_ratio = exp(satn) - 1.0;    /* # 3 function */
            //surf_ratio = 1.0 - pow(exp(-1.0/satn), 0.6); /* # 4 function */
            for (i = 0; i < rttbl->num_min; i++)
            {
                area[i] *= surf_ratio;
            }
        }
    }   /* Lichtner's 2 third law if SUF_EFF is turned on. */

    for (j = 0; j < rttbl->num_stc; j++)
    {
        Rate_spe[j] = 0.0;
    }

    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (kintbl[i].type == TST)
        {
            IAP[i] = 0.0;
            for (j = 0; j < rttbl->num_stc; j++)
            {
                IAP[i] += log10(chms->prim_actv[j]) *
                    rttbl->dep_kin[min_pos][j];
            }
            IAP[i] = pow(10, IAP[i]);
            tmpKeq = pow(10, rttbl->keq_kin[min_pos]);
            dependency[i] = 1.0;
            for (k = 0; k < kintbl[i].ndep; k++)
                dependency[i] *=
                    pow(chms->prim_actv[kintbl[i].dep_index[k]],
                    kintbl[i].dep_power[k]);
            /* Calculate the predicted rate depending on the type of rate law */
            Rate_pre[i] = area[min_pos] * (pow(10, kintbl[i].rate)) *
                dependency[i] * (1 - (IAP[i] / tmpKeq));
            /* Rate_pre: rate per reaction (mol L-1 water s-1)
             * area: m2 L-1 water
             * rate: mol m-2 s-1
             * dependency: dimensionless */
        }
        else if (kintbl[i].type == MONOD)
        {
            monodterm = 1.0;    /* re-set for new species */
            inhibterm = 1.0;    /*re-set for new species */

            /* Calculate rate */
            for (mn = 0; mn < kintbl[i].nmonod; mn++)
            {
                monodterm *= chms->prim_conc[kintbl[i].monod_index[mn]] /
                    (chms->prim_conc[kintbl[i].monod_index[mn]] +
                    kintbl[i].monod_para[mn]);
            }

            for (in = 0; in < kintbl[i].ninhib; in++)
            {
                inhibterm *= kintbl[i].inhib_para[in] /
                    (kintbl[i].inhib_para[in] +
                    chms->prim_conc[kintbl[i].inhib_index[in]]);
            }

            /* Based on CrunchTope */
            Rate_pre[i] =
                area[min_pos] * pow(10, kintbl[i].rate) * monodterm * ftemp;
        }

        for (j = 0; j < rttbl->num_stc; j++)
        {
            Rate_spe[j] += Rate_pre[i] * rttbl->dep_kin[min_pos][j];
        }
    }

    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;
        if (Rate_pre[i] < 0.0)
        {
            /* Mineral cutoff when mineral is disappearing */
            if (chms->prim_conc[min_pos + rttbl->num_stc - rttbl->num_min] < 1.0E-8)
            {
                area[min_pos] = 0.0;
            }
        }
    }

    for (i = 0; i < rttbl->num_spc; i++)
    {
        if (chemtbl[i].itype == AQUEOUS)
        {
            /* Aqueous species, saturation term for aqueous volume */
            Rate_spe[i] *= inv_sat;
        }
    }

    jcb = newDenseMat(rttbl->num_stc - rttbl->num_min,
        rttbl->num_stc - rttbl->num_min);

    if (rttbl->tmp_coup == 0)
    {
        for (i = 0; i < rttbl->num_ssc; i++)
        {
            Keq[i] = rttbl->keq[i];
        }
    }

    adh = rttbl->adh;
    bdh = rttbl->bdh;
    bdt = rttbl->bdt;

    for (i = 0; i < rttbl->num_stc; i++)
    {
        tmpconc[i] = log10(chms->prim_conc[i]);
    }
    for (i = 0; i < rttbl->num_ssc; i++)
    {
        tmpconc[i + rttbl->num_stc] = log10(chms->sec_conc[i]);
    }
    tot_cec = 0.0;
    for (i = 0; i < num_spe; i++)
    {
        tot_cec += (chemtbl[i].itype == CATION_ECHG) ?
            pow(10, tmpconc[i]) : 0.0;
    }

    I = 0;
    for (i = 0; i < num_spe; i++)
    {
        I += 0.5 * pow(10, tmpconc[i]) * chemtbl[i].charge * chemtbl[i].charge;
    }
    Iroot = sqrt(I);
    for (i = 0; i < num_spe; i++)
    {
        switch (chemtbl[i].itype)
        {
            case AQUEOUS:
                gamma[i] =
                    (-adh * chemtbl[i].charge * chemtbl[i].charge * Iroot) /
                    (1.0 + bdh * chemtbl[i].size_fac * Iroot) + bdt * I;
                break;
            case ADSORPTION:
                gamma[i] = log10(satn);
                break;
            case CATION_ECHG:
                gamma[i] = -log10(tot_cec);
                break;
            case MINERAL:
                gamma[i] = -tmpconc[i];
                break;
        }
    }

    while (maxerror > TOL)
    {
        for (i = 0; i < rttbl->num_ssc; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->num_sdc; j++)
            {
                tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
            }
            tmpval -= Keq[i] + gamma[i + rttbl->num_stc];
            tmpconc[i + rttbl->num_stc] = tmpval;
        }

        for (j = 0; j < rttbl->num_stc; j++)
        {
            Rate_spet[j] = 0.0;
        }

        for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
        {
            min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

            if (kintbl[i].type == TST)
            {
                IAP[i] = 0.0;
                for (j = 0; j < rttbl->num_stc; j++)
                {
                    if (chemtbl[j].itype != MINERAL)
                    {
                        IAP[i] += (tmpconc[j] + gamma[j]) *
                            rttbl->dep_kin[min_pos][j];
                    }
                }
                IAP[i] = pow(10, IAP[i]);
                tmpKeq = pow(10, rttbl->keq_kin[min_pos]);
                /*
                 * if ( IAP[i] < tmpKeq)
                 * rct_drct[i] = 1.0;
                 * if ( IAP[i] > tmpKeq)
                 * rct_drct[i] = -1.0;
                 * if ( IAP[i] == tmpKeq)
                 * rct_drct[i] = 0.0;
                 */
                dependency[i] = 0.0;
                for (k = 0; k < kintbl[i].ndep; k++)
                {
                    dependency[i] += (tmpconc[kintbl[i].dep_index[k]] +
                        gamma[kintbl[i].dep_index[k]]) *
                        kintbl[i].dep_power[k];
                }
                dependency[i] = pow(10, dependency[i]);
                /* Calculate predicted rate depending on type of rate law */
                Rate_pre[i] = area[min_pos] * (pow(10, kintbl[i].rate)) *
                    dependency[i] * (1.0 - (IAP[i] / tmpKeq));
                /* Rate_pre: in mol / L water / s
                 * area: m2/L water
                 * rate: mol/m2/s
                 * dependency: dimensionless; */
            }
            else if (kintbl[i].type == MONOD)
            {
                monodterm = 1.0;
                inhibterm = 1.0;

                /* Calculate rate */
                for (mn = 0; mn < kintbl[i].nmonod; mn++)
                {
                    monodterm *= chms->prim_conc[kintbl[i].monod_index[mn]] /
                        (chms->prim_conc[kintbl[i].monod_index[mn]] +
                        kintbl[i].monod_para[mn]);
                }

                for (in = 0; in < kintbl[i].ninhib; in++)
                {
                    inhibterm *= kintbl[i].inhib_para[in] /
                        (kintbl[i].inhib_para[in] +
                        chms->prim_conc[kintbl[i].inhib_index[in]]);
                }

                /* Based on CrunchTope */
                Rate_pre[i] =
                    area[min_pos] * pow(10, kintbl[i].rate) * monodterm * ftemp;
            }

            for (j = 0; j < rttbl->num_stc; j++)
            {
                Rate_spet[j] += Rate_pre[i] * rttbl->dep_kin[min_pos][j];
            }
            /* Adjust the unit of the calculated rate. Note that for mineral,
             * the unit of rate and the unit of concentration are mol/L porous
             * media. For the aqueous species, the unit of the rate and the unit
             * of the concentration are mol/L pm and mol/L water respectively.*/
        }

        for (i = 0; i < rttbl->num_spc; i++)
        {
            Rate_spet[i] *= (chemtbl[i].itype == AQUEOUS) ? inv_sat : 1.0;
        }

        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
            {
                tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
            }
            tot_conc[i] = tmpval;
            residue[i] = tmpval - (chms->tot_conc[i] +
                (Rate_spe[i] + Rate_spet[i]) * stepsize * 0.5);
        }
        if (control % SKIP_JACOB == 0)
        {
            for (k = 0; k < rttbl->num_stc - rttbl->num_min; k++)
            {
                tmpconc[k] += tmpprb;
                for (i = 0; i < rttbl->num_ssc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_sdc; j++)
                    {
                        tmpval +=
                            (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
                    }
                    tmpval -= Keq[i] + gamma[i + rttbl->num_stc];
                    tmpconc[i + rttbl->num_stc] = tmpval;
                }
                for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                    {
                        tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                    }
                    residue_t[i] = tmpval - (chms->tot_conc[i] +
                        (Rate_spe[i] + Rate_spet[i]) * stepsize * 0.5);
                    jcb[k][i] = (residue_t[i] - residue[i]) * tmpprb_inv;
                }
                tmpconc[k] -= tmpprb;
            }
        }
        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            x_[i] = -residue[i];
        }

        pivot_flg = denseGETRF(jcb, rttbl->num_stc - rttbl->num_min,
            rttbl->num_stc - rttbl->num_min, p);
        if (pivot_flg != 0)
        {
            return 1;
        }

        denseGETRS(jcb, rttbl->num_stc - rttbl->num_min, p, x_);

        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            if (fabs(x_[i]) < 0.3)
            {
                tmpconc[i] += x_[i];
            }
            else
            {
                tmpconc[i] += (x_[i] < 0) ? -0.3 : 0.3;
            }
            error[i] = residue[i] / tot_conc[i];
        }
        maxerror = fabs(error[0]);
        for (i = 1; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            maxerror = MAX(fabs(error[i]), maxerror);
        }
        control++;
        if (control > 10)
        {
            return 1;
        }
    }

    destroyMat(jcb);

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->num_sdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
        }
        tmpval -= Keq[i] + gamma[i + rttbl->num_stc];
        tmpconc[i + rttbl->num_stc] = tmpval;
    }

    for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
        {
            tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
        }
        tot_conc[i] = tmpval;
        residue[i] = tmpval - chms->tot_conc[i];
        error[i] = residue[i] / tot_conc[i];
    }
    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        if (i < rttbl->num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms->tot_conc[i] +=
                    (Rate_spe[i] + Rate_spet[i]) * stepsize * 0.5;
                chms->prim_actv[i] = 1.0;
                chms->prim_conc[i] = chms->tot_conc[i];
            }
            else
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = pow(10, (tmpconc[i] + gamma[i]));
                chms->tot_conc[i] = tot_conc[i];
            }
        }
        else
        {
            chms->sec_conc[i - rttbl->num_stc] = pow(10, tmpconc[i]);
#if TEMP_DISABLED
            chms->s_actv[i - rttbl->num_stc] = pow(10, (tmpconc[i] + gamma[i]));
#endif
        }
    }

    return 0;
}

void ReactControl(const chemtbl_struct chemtbl[], const kintbl_struct kintbl[],
    const rttbl_struct *rttbl, double stepsize, double satn, double ftemp,
    chmstate_struct *chms, double react_flux[])
{
    double          t_conc0[MAXSPS];
    double          substep;
    double          step_counter = 0.0;
    int             flag;
    int             k;

    for (k = 0; k < rttbl->num_spc; k++)
    {
        t_conc0[k] = chms->tot_conc[k];
    }

    substep = stepsize;

    while (1.0 - step_counter / stepsize > 1.0E-10 && substep > 30.0)
    {
        flag = _React(substep, chemtbl, kintbl, rttbl, satn, ftemp, chms);

        if (flag == 0)
        {
            /* Reaction passed with current step */
            step_counter += substep;
        }
        else
        {
            substep *= 0.5;
        }
    }

    if (roundi(step_counter) == roundi(stepsize))
    {
        if (roundi(substep) != roundi(stepsize))
        {
            pihm_printf(VL_VERBOSE,
                " Reaction passed with minimum step of %f s.\n", substep);
        }
    }
    else
    {
        pihm_printf(VL_VERBOSE, " Reaction failed.\n");
    }

    for (k = 0; k < rttbl->num_spc; k++)
    {
        react_flux[k] =
            (chms->tot_conc[k] - t_conc0[k]) / stepsize;
    }
}

double SoilTempFactor(double stc)
{
    return pow(2.0, (stc - TFREEZ - 20.0) / 10.0);
}
