#include "pihm.h"

#define TOL        1E-7

void Speciation(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < rttbl->num_spc; k++)
        {
            river[i].chms.prim_conc[k] = river[i].chms.tot_conc[k];
        }

        if (river[i].ws.stage > DEPTHR)
        {
            _Speciation(chemtbl, rttbl, 0, &river[i].chms);
        }
    }
}

int _Speciation(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    int speciation_flg, chmstate_struct *chms)
{
    /* if speciation flg = 1, pH is defined
     * if speciation flg = 0, all defined value is total concentration */
    int             i, j, k;
    int             num_spe = rttbl->num_stc + rttbl->num_ssc;
    double          tmpval;
    double          tmpprb = 1E-2;
    double          I;
    double          Iroot;
    double          residue[MAXSPS];
    double          residue_t[MAXSPS];
    double          tmpconc[MAXSPS];
    double          tot_conc[MAXSPS];
    double          error[MAXSPS];
    double          gamma[MAXSPS];
    double          Keq[MAXSPS];
    double          current_totconc[MAXSPS];
    double          adh;
    double          bdh;
    double          bdt;
    realtype      **jcb;

    for (k = 0; k < MAXSPS; k++)
    {
        residue[k] = 0.0;
        residue_t[k] = 0.0;
        tmpconc[k] = 0.0;
        tot_conc[k] = 0.0;
        error[k] = 0.0;
        gamma[k] = 0.0;
        Keq[k] = 0.0;
        current_totconc[k] = 0.0;
    }

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
        /* Using log10 conc as the primary unknowns. Works better because
         * negative numbers are not a problem. */
        tmpconc[i] = log10(chms->prim_conc[i]);
    }

    for (i = 0; i < rttbl->num_ssc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->num_sdc; j++)
        {
            tmpval += tmpconc[j] * rttbl->dep_mtx[i][j];
        }
        tmpval -= Keq[i];
        tmpconc[i + rttbl->num_stc] = tmpval;
    }

    for (i = 0; i < rttbl->num_stc; i++)
    {
        current_totconc[i] = chms->tot_conc[i];
    }

    if (speciation_flg == 1)
    {
        /* pH is defined, total concentration is calculated from the activity of
         * H. Dependency is the same but the total concentration for H need not
         * be solved */
        jcb = newDenseMat(rttbl->num_stc - 1, rttbl->num_stc - 1);
        sunindextype    p[MAXSPS];
        realtype        x_[MAXSPS];
        double          maxerror = 1;

        while (maxerror > TOL)
        {
            if (rttbl->actv_mode == 1)
            {
                I = 0;
                /* Calculate the ionic strength in this block */
                for (i = 0; i < num_spe; i++)
                {
                    I += 0.5 * pow(10, tmpconc[i]) *
                        chemtbl[i].charge * chemtbl[i].charge;
                }
                Iroot = sqrt(I);
                for (i = 0; i < num_spe; i++)
                {
                    if (chemtbl[i].itype == MINERAL)
                    {
                        gamma[i] = -tmpconc[i];
                    }
                    /* aqueous species in the unit of mol/L, however the solids
                     * are in the unit of mol/L porous media
                     * the activity of solid is 1, the log10 of activity is 0.
                     * by assigning gamma[minerals] to negative of the
                     * tmpconc[minerals], we ensured the log 10 of activity of
                     * solids are 0*/
                    else
                    {
                        gamma[i] = (-adh * chemtbl[i].charge *
                            chemtbl[i].charge * Iroot) /
                            (1 + bdh * chemtbl[i].size_fac * Iroot) + bdt * I;
                        }
                    if (strcmp(chemtbl[i].name, "'H+'") == 0)
                    {
                        tmpconc[i] = log10(chms->prim_actv[i]) - gamma[i];
                    }
                }
                /* gamma stores log10gamma[i] */
            }

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
            for (i = 0; i < rttbl->num_stc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < num_spe; j++)
                {
                    tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                }
                tot_conc[i] = tmpval;
                if (strcmp(chemtbl[i].name, "'H+'") == 0)
                {
                    chms->tot_conc[i] = tot_conc[i];
                }
                residue[i] = tmpval - chms->tot_conc[i];
                /* update the total concentration of H+ for later stage RT at
                 * initialization */
            }
            int             row, col;
            col = 0;
            for (k = 0; k < rttbl->num_stc; k++)
            {
                if (strcmp(chemtbl[k].name, "'H+'") != 0)
                {
                    tmpconc[k] += tmpprb;
                    for (i = 0; i < rttbl->num_ssc; i++)
                    {
                        tmpval = 0.0;
                        for (j = 0; j < rttbl->num_sdc; j++)
                        {
                            tmpval += (tmpconc[j] + gamma[j]) *
                                rttbl->dep_mtx[i][j];
                        }
                        tmpval -= Keq[i] + gamma[i + rttbl->num_stc];
                        tmpconc[i + rttbl->num_stc] = tmpval;
                    }
                    row = 0;
                    for (i = 0; i < rttbl->num_stc; i++)
                    {
                        if (strcmp(chemtbl[i].name, "'H+'") != 0)
                        {
                            tmpval = 0.0;
                            for (j = 0; j < rttbl->num_stc + rttbl->num_ssc;
                                j++)
                            {
                                tmpval += rttbl->conc_contrib[i][j] *
                                    pow(10, tmpconc[j]);
                            }
                            residue_t[i] = tmpval - chms->tot_conc[i];
                            jcb[col][row] =
                                (residue_t[i] - residue[i]) / (tmpprb);
                            row++;
                        }
                    }
                    tmpconc[k] -= tmpprb;
                    col++;
                }
            }
            row = 0;
            for (i = 0; i < rttbl->num_stc; i++)
                if (strcmp(chemtbl[i].name, "'H+'") != 0)
                    x_[row++] = -residue[i];
            if (denseGETRF(jcb, rttbl->num_stc - 1, rttbl->num_stc - 1, p) != 0)
            {
                pihm_printf(VL_ERROR, "Speciation error.\n");
                pihm_exit(EXIT_FAILURE);
            }
            denseGETRS(jcb, rttbl->num_stc - 1, p, x_);

            row = 0;
            for (i = 0; i < rttbl->num_stc; i++)
            {
                if (strcmp(chemtbl[i].name, "'H+'") != 0)
                    tmpconc[i] += x_[row++];
                error[i] = residue[i] / tot_conc[i];
            }
            maxerror = fabs(error[0]);
            for (i = 1; i < rttbl->num_stc; i++)
                if (fabs(error[i]) > maxerror)
                    maxerror = fabs(error[i]);
        }
    }
    else
    {
        jcb = newDenseMat(rttbl->num_stc, rttbl->num_stc);
        sunindextype    p[MAXSPS];
        realtype        x_[MAXSPS];
        double          maxerror = 1;

        while (maxerror > TOL)
        {
            if (rttbl->actv_mode == 1)
            {
                I = 0.0;
                /* Calculate the ionic strength in this block */
                for (i = 0; i < num_spe; i++)
                {
                    I += 0.5 * pow(10, tmpconc[i]) *
                        chemtbl[i].charge * chemtbl[i].charge;
                }
                Iroot = sqrt(I);
                for (i = 0; i < num_spe; i++)
                {
                    if (chemtbl[i].itype == MINERAL)
                        gamma[i] = -tmpconc[i];
                    else
                        gamma[i] =
                            (-adh * chemtbl[i].charge *
                            chemtbl[i].charge * Iroot) /
                            (1.0 + bdh * chemtbl[i].size_fac * Iroot) + bdt * I;
                }
            }
            /* gamma stores log10gamma[i]. */
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
            for (i = 0; i < rttbl->num_stc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                {
                    tmpval += rttbl->conc_contrib[i][j] * pow(10, tmpconc[j]);
                }
                tot_conc[i] = tmpval;
                residue[i] = tmpval - chms->tot_conc[i];
            }

            for (k = 0; k < rttbl->num_stc; k++)
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
                for (i = 0; i < rttbl->num_stc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->num_stc + rttbl->num_ssc; j++)
                    {
                        tmpval += rttbl->conc_contrib[i][j] *
                            pow(10, tmpconc[j]);
                    }
                    residue_t[i] = tmpval - chms->tot_conc[i];
                    jcb[k][i] = (residue_t[i] - residue[i]) / tmpprb;
                }
                tmpconc[k] -= tmpprb;
            }
            for (i = 0; i < rttbl->num_stc; i++)
            {
                x_[i] = -residue[i];
            }
            if (denseGETRF(jcb, rttbl->num_stc, rttbl->num_stc, p) != 0)
            {
                pihm_printf(VL_ERROR, "Speciation error.\n");
                pihm_exit(EXIT_FAILURE);
            }
            denseGETRS(jcb, rttbl->num_stc, p, x_);
            for (i = 0; i < rttbl->num_stc; i++)
            {
                tmpconc[i] += x_[i];
                error[i] = residue[i] / tot_conc[i];
            }
            maxerror = fabs(error[0]);
            for (i = 1; i < rttbl->num_stc; i++)
            {
                maxerror = MAX(fabs(error[i]), maxerror);
            }
        }
    }
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
    for (i = 0; i < rttbl->num_stc; i++)
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
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = 1.0;
            }
            else
            {
                chms->prim_conc[i] = pow(10, tmpconc[i]);
                chms->prim_actv[i] = pow(10, (tmpconc[i] + gamma[i]));
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
    destroyMat(jcb);

    return (0);
}
