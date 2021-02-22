#include "pihm.h"

#define TOL                     1E-7
#define SKIP_JACOB              1

void Reaction(double stepsize, const chemtbl_struct chemtbl[], const kintbl_struct kintbl[], const rttbl_struct *rttbl,
    elem_struct elem[])
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

        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity + elem[i].soil.depth * elem[i].soil.smcmin;

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

            ftemp = SoilTempFactor(rttbl->q10, avg_stc);

            ReactControl(chemtbl, kintbl, rttbl, stepsize, satn, ftemp, &elem[i].chms, elem[i].chmf.react);
        }

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].snksrc = elem[i].chmf.react[k] * storage;
        }

#if defined(_DGW_)
        storage = (elem[i].ws.unsat_geol + elem[i].ws.gw_geol) * elem[i].geol.porosity +
            elem[i].geol.depth * elem[i].geol.smcmin;

        satn = (elem[i].ws.unsat_geol + elem[i].ws.gw_geol) / elem[i].geol.depth;
        satn = MAX(satn, SATMIN);
        satn = MIN(satn, 1.0);

        for (k = 0; k < rttbl->num_spc; k++)
        {
            elem[i].chmf.react_geol[k] = 0.0;
        }

        // Deep zone
        if (satn > 1.0E-2)
        {
            ftemp = SoilTempFactor(rttbl->q10, elem[i].ps.tbot);
            ReactControl(chemtbl, kintbl, rttbl, stepsize, satn, ftemp, &elem[i].chms_geol, elem[i].chmf.react_geol);
        }

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].snksrc_geol = elem[i].chmf.react_geol[k] * storage;
        }
#endif
    }
}

int _React(double stepsize, const chemtbl_struct chemtbl[], const kintbl_struct kintbl[], const rttbl_struct *rttbl,
    double satn, double ftemp, chmstate_struct *chms)
{
    int             i, j, k;
    int             kmonod, kinhib;
    int             control;
    int             min_pos;
    int             pivot_flg;
    double          monodterm = 1.0, inhibterm = 1.0;
    double          residue[MAXSPS];
    double          residue_t[MAXSPS];
    double          tmpconc[MAXSPS];
    double          tot_conc[MAXSPS];
    double          area[MAXSPS];
    double          gamma[MAXSPS];
    double          keq[MAXSPS];
    double          rate_pre[MAXSPS];
    double          iap[MAXSPS];
    double          dependency[MAXSPS];
    double          rate_spe[MAXSPS];
    double          rate_spet[MAXSPS];
    double          tmpval;
    double          inv_sat;
    double          imat;
    double          iroot;
    double          temp_keq;
    double          adh;
    double          bdh;
    double          bdt;
    double          max_error;
    double          surf_ratio;
    double          tot_cec;
    realtype      **jcb;
    sunindextype    p[MAXSPS];
    realtype        x[MAXSPS];
    const int       SUFEFF = 1;
    const double    TMPPRB = 1.0E-2;
    const double    TMPPRB_INV = 1.0 / TMPPRB;

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
            surf_ratio = (satn <= rttbl->sw_thld) ?
                pow(satn / rttbl->sw_thld, rttbl->sw_exp) : pow((1.0 - satn) / (1.0 - rttbl->sw_thld), rttbl->sw_exp);

            for (i = 0; i < rttbl->num_min; i++)
            {
                area[i] *= surf_ratio;
            }
        }
    }   // Lichtner's 2 third law if SUF_EFF is turned on

    for (i = 0; i < rttbl->num_stc; i++)
    {
        rate_spe[i] = 0.0;
    }

    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (kintbl[i].type == TST)
        {
            iap[i] = 0.0;
            for (j = 0; j < rttbl->num_stc; j++)
            {
                iap[i] += log10(chms->prim_actv[j]) * rttbl->dep_kin[i][j];
            }
            iap[i] = pow(10, iap[i]);

            temp_keq = pow(10, rttbl->keq_kin[i]);

            dependency[i] = 1.0;
            for (k = 0; k < kintbl[i].ndep; k++)
            {
                dependency[i] *= pow(chms->prim_actv[kintbl[i].dep_index[k]], kintbl[i].dep_power[k]);
            }

            // Calculate the predicted rate depending on the type of rate law
            //   rate_pre: rate per reaction (mol L-1 water s-1)
            //   area: m2 L-1 water
            //   rate: mol m-2 s-1
            //   dependency: dimensionless
            rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * dependency[i] * (1.0 - iap[i] / temp_keq);
        }
        else if (kintbl[i].type == MONOD)
        {
            monodterm = 1.0;    // re-set for new species
            inhibterm = 1.0;    //re-set for new species

            // Calculate rate
            for (kmonod = 0; kmonod < kintbl[i].nmonod; kmonod++)
            {
                monodterm *=
                    chms->prim_conc[kintbl[i].monod_index[kmonod]] / (chms->prim_conc[kintbl[i].monod_index[kmonod]] +
                    kintbl[i].monod_para[kmonod]);
            }

            for (kinhib = 0; kinhib < kintbl[i].ninhib; kinhib++)
            {
                inhibterm *= kintbl[i].inhib_para[kinhib] / (kintbl[i].inhib_para[kinhib] +
                    chms->prim_conc[kintbl[i].inhib_index[kinhib]]);
            }

            // Based on CrunchTope
            rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm * ftemp;
        }

        for (j = 0; j < rttbl->num_stc; j++)
        {
            rate_spe[j] += rate_pre[i] * rttbl->dep_kin[i][j];
        }
    }

    for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
    {
        min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

        if (rate_pre[i] < 0.0)
        {
            // Mineral cutoff when mineral is disappearing
            area[min_pos] = (chms->prim_conc[kintbl[i].position] < 1.0E-8) ? 0.0 : area[min_pos];
        }
    }

    for (i = 0; i < rttbl->num_spc; i++)
    {
        // Aqueous species, saturation term for aqueous volume
        rate_spe[i] *= (chemtbl[i].itype == AQUEOUS) ? inv_sat : 1.0;
    }

    jcb = newDenseMat(rttbl->num_stc - rttbl->num_min, rttbl->num_stc - rttbl->num_min);

    if (rttbl->tmp_coup == 0)
    {
        for (i = 0; i < rttbl->num_ssc; i++)
        {
            keq[i] = rttbl->keq[i];
        }
    }

    adh = rttbl->adh;
    bdh = rttbl->bdh;
    bdt = rttbl->bdt;

    tot_cec = 0.0;
    imat = 0;
    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        tmpconc[i] = (i < rttbl->num_stc) ? log10(chms->prim_conc[i]) : log10(chms->sec_conc[i - rttbl->num_stc]);

        tot_cec += (chemtbl[i].itype == CATION_ECHG) ? pow(10, tmpconc[i]) : 0.0;

        imat += 0.5 * pow(10, tmpconc[i]) * chemtbl[i].charge * chemtbl[i].charge;
    }
    iroot = sqrt(imat);

    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        switch (chemtbl[i].itype)
        {
            case AQUEOUS:
                gamma[i] = (-adh * chemtbl[i].charge * chemtbl[i].charge * iroot) /
                    (1.0 + bdh * chemtbl[i].size_fac * iroot) + bdt * imat;
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

    control = 0;
    max_error = 0.0;
    do
    {
        for (i = 0; i < rttbl->num_ssc; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->num_sdc; j++)
            {
                tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
            }
            tmpval -= keq[i] + gamma[i + rttbl->num_stc];
            tmpconc[i + rttbl->num_stc] = tmpval;
        }

        for (j = 0; j < rttbl->num_stc; j++)
        {
            rate_spet[j] = 0.0;
        }

        for (i = 0; i < rttbl->num_mkr + rttbl->num_akr; i++)
        {
            min_pos = kintbl[i].position - rttbl->num_stc + rttbl->num_min;

            if (kintbl[i].type == TST)
            {
                iap[i] = 0.0;
                for (j = 0; j < rttbl->num_stc; j++)
                {
                    iap[i] += (chemtbl[j].itype != MINERAL) ? (tmpconc[j] + gamma[j]) * rttbl->dep_kin[i][j] : 0.0;
                }
                iap[i] = pow(10, iap[i]);
                temp_keq = pow(10, rttbl->keq_kin[i]);

#if NOT_YET_IMPLEMENTED
                if (iap[i] < temp_keq)
                {
                    rct_drct[i] = 1.0;
                }
                else if (iap[i] > temp_keq)
                {
                    rct_drct[i] = -1.0;
                }
                else
                {
                    rct_drct[i] = 0.0;
                }
#endif

                dependency[i] = 0.0;
                for (k = 0; k < kintbl[i].ndep; k++)
                {
                    dependency[i] += (tmpconc[kintbl[i].dep_index[k]] +
                        gamma[kintbl[i].dep_index[k]]) * kintbl[i].dep_power[k];
                }
                dependency[i] = pow(10, dependency[i]);

                // Calculate predicted rate depending on type of rate law
                // rate_pre: in mol / L water / s
                // area: m2/L water
                // rate: mol/m2/s
                // dependency: dimensionless
                rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * dependency[i] * (1.0 - (iap[i] / temp_keq));
            }
            else if (kintbl[i].type == MONOD)
            {
                monodterm = 1.0;
                inhibterm = 1.0;

                // Calculate rate
                for (kmonod = 0; kmonod < kintbl[i].nmonod; kmonod++)
                {
                    monodterm *= chms->prim_conc[kintbl[i].monod_index[kmonod]] /
                        (chms->prim_conc[kintbl[i].monod_index[kmonod]] +
                        kintbl[i].monod_para[kmonod]);
                }

                for (kinhib = 0; kinhib < kintbl[i].ninhib; kinhib++)
                {
                    inhibterm *= kintbl[i].inhib_para[kinhib] / (kintbl[i].inhib_para[kinhib] +
                        chms->prim_conc[kintbl[i].inhib_index[kinhib]]);
                }

                // Based on CrunchTope
                rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm * ftemp;
            }

            for (j = 0; j < rttbl->num_stc; j++)
            {
                rate_spet[j] += rate_pre[i] * rttbl->dep_kin[i][j];
            }
            // Adjust the unit of the calculated rate. Note that for mineral, the unit of rate and the unit of
            // concentration are mol/L porous media. For the aqueous species, the unit of the rate and the unit of the
            // concentration are mol/L pm and mol/L water respectively.
        }

        for (i = 0; i < rttbl->num_spc; i++)
        {
            rate_spet[i] *= (chemtbl[i].itype == AQUEOUS) ? inv_sat : 1.0;
        }

        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
            {
                tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
            }
            tot_conc[i] = tmpval;
            residue[i] = tmpval - (chms->tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
        }

        if (control % SKIP_JACOB == 0)
        {
            for (k = 0; k < rttbl->num_stc - rttbl->num_min; k++)
            {
                tmpconc[k] += TMPPRB;
                for (i = 0; i < rttbl->num_ssc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_sdc; j++)
                    {
                        tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
                    }
                    tmpval -= keq[i] + gamma[i + rttbl->num_stc];
                    tmpconc[i + rttbl->num_stc] = tmpval;
                }
                for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                    {
                        tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                    }
                    residue_t[i] = tmpval - (chms->tot_conc[i] + (rate_spe[i] + rate_spet[i]) * stepsize * 0.5);
                    jcb[k][i] = (residue_t[i] - residue[i]) * TMPPRB_INV;
                }
                tmpconc[k] -= TMPPRB;
            }
        }
        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            x[i] = -residue[i];
        }

        pivot_flg = denseGETRF(jcb, rttbl->num_stc - rttbl->num_min, rttbl->num_stc - rttbl->num_min, p);
        if (pivot_flg != 0)
        {
            return 1;
        }

        denseGETRS(jcb, rttbl->num_stc - rttbl->num_min, p, x);

        max_error = 0.0;
        for (i = 0; i < rttbl->num_stc - rttbl->num_min; i++)
        {
            if (fabs(x[i]) < 0.3)
            {
                tmpconc[i] += x[i];
            }
            else
            {
                tmpconc[i] += (x[i] < 0) ? -0.3 : 0.3;
            }

            max_error = MAX(fabs(residue[i] / tot_conc[i]), max_error);
        }

        control++;
        if (control > 10)
        {
            return 1;
        }
    } while (max_error > TOL);

    destroyMat(jcb);

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->num_sdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl->dep_mtx[i][j];
        }
        tmpval -= keq[i] + gamma[i + rttbl->num_stc];
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
    }

    for (i = 0; i < rttbl->num_stc + rttbl->num_ssc; i++)
    {
        if (i < rttbl->num_stc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms->tot_conc[i] += (rate_spe[i] + rate_spet[i]) * stepsize * 0.5;
                chms->prim_actv[i] = 1.0;
                chms->prim_conc[i] = chms->tot_conc[i];
            }
            else
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = pow(10, tmpconc[i] + gamma[i]);
                chms->tot_conc[i] = tot_conc[i];
            }
        }
        else
        {
            chms->sec_conc[i - rttbl->num_stc] = pow(10, tmpconc[i]);
#if TEMP_DISABLED
            chms->s_actv[i - rttbl->num_stc] = pow(10, tmpconc[i] + gamma[i]);
#endif
        }
    }

    return 0;
}

void ReactControl(const chemtbl_struct chemtbl[], const kintbl_struct kintbl[], const rttbl_struct *rttbl,
    double stepsize, double satn, double ftemp, chmstate_struct *chms, double react_flux[])
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
            // Reaction passed with current step
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
            pihm_printf(VL_VERBOSE, " Reaction passed with minimum step of %f s.\n", substep);
        }
    }
    else
    {
        pihm_printf(VL_VERBOSE, " Reaction failed.\n");
    }

    for (k = 0; k < rttbl->num_spc; k++)
    {
        react_flux[k] = (chms->tot_conc[k] - t_conc0[k]) / stepsize;
    }
}

double SoilTempFactor(double q10, double stc)
{
    return pow(q10, (stc - TFREEZ - 20.0) / 10.0);
}
