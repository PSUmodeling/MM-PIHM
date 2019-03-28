#include "pihm.h"

#define TOL        1E-7

void Speciation(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    elem_struct elem[], river_struct river[])
{
    int             i;

    if (rttbl->RecFlg == KIN_REACTION)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nriver; i++)
        {
            int             k;

            for (k = 0; k < rttbl->NumStc; k++)
            {
                river[i].chms_stream.p_conc[k] = (chemtbl[k].itype == MINERAL) ?
                    river[i].chms_stream.t_conc[k] :
                    river[i].chms_stream.t_conc[k] * 0.1;

                river[i].chms_rivbed.p_conc[k] = (chemtbl[k].itype == MINERAL) ?
                    river[i].chms_rivbed.t_conc[k] :
                    river[i].chms_rivbed.t_conc[k] * 0.1;
            }

            _Speciation(chemtbl, rttbl, 0, &river[i].chms_stream);

            _Speciation(chemtbl, rttbl, 0, &river[i].chms_rivbed);
        }
    }
    else
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            _Speciation(chemtbl, rttbl, 0, &elem[i].chms_unsat);

            _Speciation(chemtbl, rttbl, 0, &elem[i].chms_gw);
        }
    }
}

int _Speciation(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    int speciation_flg, chmstate_struct *chms)
{
    /* if speciation flg = 1, pH is defined
     * if speciation flg = 0, all defined value is total concentration */
    int             i, j, k;
    int             num_spe = rttbl->NumStc + rttbl->NumSsc;
    double          tmpval;
    double          tmpprb = 1E-2;
    double          I;
    double          Iroot;
    double          residue[MAXSPS];
    double          residue_t[MAXSPS];
    double          tmpconc[MAXSPS];
    double          totconc[MAXSPS];
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
        totconc[k] = 0.0;
        error[k] = 0.0;
        gamma[k] = 0.0;
        Keq[k] = 0.0;
        current_totconc[k] = 0.0;
    }

    if (rttbl->TEMcpl == 0)
    {
        for (i = 0; i < rttbl->NumSsc; i++)
        {
            Keq[i] = rttbl->Keq[i];
        }
    }

    adh = rttbl->adh;
    bdh = rttbl->bdh;
    bdt = rttbl->bdt;

    for (i = 0; i < rttbl->NumStc; i++)
    {
        /* Using log10 conc as the primary unknowns. Works better because
         * negative numbers are not a problem. */
        tmpconc[i] = log10(chms->p_conc[i]);
    }

    for (i = 0; i < rttbl->NumSsc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->NumSdc; j++)
        {
            tmpval += tmpconc[j] * rttbl->Dependency[i][j];
        }
        tmpval -= Keq[i];
        tmpconc[i + rttbl->NumStc] = tmpval;
    }

    for (i = 0; i < rttbl->NumStc; i++)
    {
        current_totconc[i] = chms->t_conc[i];
    }

    if (speciation_flg == 1)
    {
        /* pH is defined, total concentration is calculated from the activity of
         * H. Dependency is the same but the total concentration for H need not
         * be solved */
        jcb = newDenseMat(rttbl->NumStc - 1, rttbl->NumStc - 1);
        long int       *p;
        realtype       *x_;
        double          maxerror = 1;

        p = (long int *)malloc((rttbl->NumStc - 1) * sizeof(long int));
        x_ = (realtype *)malloc((rttbl->NumStc - 1) * sizeof(realtype));

        while (maxerror > TOL)
        {
            if (rttbl->ACTmod == 1)
            {
                I = 0;
                /* Calculate the ionic strength in this block */
                for (i = 0; i < num_spe; i++)
                {
                    I += 0.5 * pow(10, tmpconc[i]) *
                        chemtbl[i].Charge * chemtbl[i].Charge;
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
                        gamma[i] = (-adh * chemtbl[i].Charge *
                            chemtbl[i].Charge * Iroot) /
                            (1 + bdh * chemtbl[i].SizeF * Iroot) + bdt * I;
                        }
                    if (strcmp(chemtbl[i].ChemName, "'H+'") == 0)
                    {
                        tmpconc[i] = log10(chms->p_actv[i]) - gamma[i];
                    }
                }
                /* gamma stores log10gamma[i] */
            }

            for (i = 0; i < rttbl->NumSsc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < rttbl->NumSdc; j++)
                {
                    tmpval += (tmpconc[j] + gamma[j]) * rttbl->Dependency[i][j];
                }
                tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                tmpconc[i + rttbl->NumStc] = tmpval;
            }
            for (i = 0; i < rttbl->NumStc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < num_spe; j++)
                {
                    tmpval += rttbl->Totalconc[i][j] * pow(10, tmpconc[j]);
                }
                totconc[i] = tmpval;
                if (strcmp(chemtbl[i].ChemName, "'H+'") == 0)
                {
                    chms->t_conc[i] = totconc[i];
                }
                residue[i] = tmpval - chms->t_conc[i];
                /* update the total concentration of H+ for later stage RT at
                 * initialization */
            }
            int             row, col;
            col = 0;
            for (k = 0; k < rttbl->NumStc; k++)
            {
                if (strcmp(chemtbl[k].ChemName, "'H+'") != 0)
                {
                    tmpconc[k] += tmpprb;
                    for (i = 0; i < rttbl->NumSsc; i++)
                    {
                        tmpval = 0.0;
                        for (j = 0; j < rttbl->NumSdc; j++)
                        {
                            tmpval += (tmpconc[j] + gamma[j]) *
                                rttbl->Dependency[i][j];
                        }
                        tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                        tmpconc[i + rttbl->NumStc] = tmpval;
                    }
                    row = 0;
                    for (i = 0; i < rttbl->NumStc; i++)
                    {
                        if (strcmp(chemtbl[i].ChemName, "'H+'") != 0)
                        {
                            tmpval = 0.0;
                            for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
                            {
                                tmpval += rttbl->Totalconc[i][j] *
                                    pow(10, tmpconc[j]);
                            }
                            residue_t[i] = tmpval - chms->t_conc[i];
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
            for (i = 0; i < rttbl->NumStc; i++)
                if (strcmp(chemtbl[i].ChemName, "'H+'") != 0)
                    x_[row++] = -residue[i];
            if (denseGETRF(jcb, rttbl->NumStc - 1, rttbl->NumStc - 1, p) != 0)
            {
                PIHMprintf(VL_ERROR, "Speciation error.\n");
                PIHMexit(EXIT_FAILURE);
            }
            denseGETRS(jcb, rttbl->NumStc - 1, p, x_);

            row = 0;
            for (i = 0; i < rttbl->NumStc; i++)
            {
                if (strcmp(chemtbl[i].ChemName, "'H+'") != 0)
                    tmpconc[i] += x_[row++];
                error[i] = residue[i] / totconc[i];
            }
            maxerror = fabs(error[0]);
            for (i = 1; i < rttbl->NumStc; i++)
                if (fabs(error[i]) > maxerror)
                    maxerror = fabs(error[i]);
        }

        free(p);
        free(x_);
    }
    else
    {
        jcb = newDenseMat(rttbl->NumStc, rttbl->NumStc);
        long int       *p;
        realtype       *x_;
        double          maxerror = 1;

        p = (long int *)malloc(rttbl->NumStc * sizeof(long int));
        x_ = (realtype *)malloc(rttbl->NumStc * sizeof(realtype));

        while (maxerror > TOL)
        {
            if (rttbl->ACTmod == 1)
            {
                I = 0.0;
                /* Calculate the ionic strength in this block */
                for (i = 0; i < num_spe; i++)
                {
                    I += 0.5 * pow(10, tmpconc[i]) *
                        chemtbl[i].Charge * chemtbl[i].Charge;
                }
                Iroot = sqrt(I);
                for (i = 0; i < num_spe; i++)
                {
                    if (chemtbl[i].itype == MINERAL)
                        gamma[i] = -tmpconc[i];
                    else
                        gamma[i] =
                            (-adh * chemtbl[i].Charge *
                            chemtbl[i].Charge * Iroot) /
                            (1.0 + bdh * chemtbl[i].SizeF * Iroot) + bdt * I;
                }
            }
            /* gamma stores log10gamma[i]. */
            for (i = 0; i < rttbl->NumSsc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < rttbl->NumSdc; j++)
                {
                    tmpval += (tmpconc[j] + gamma[j]) * rttbl->Dependency[i][j];
                }
                tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                tmpconc[i + rttbl->NumStc] = tmpval;
            }
            for (i = 0; i < rttbl->NumStc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
                {
                    tmpval += rttbl->Totalconc[i][j] * pow(10, tmpconc[j]);
                }
                totconc[i] = tmpval;
                residue[i] = tmpval - chms->t_conc[i];
            }

            for (k = 0; k < rttbl->NumStc; k++)
            {
                tmpconc[k] += tmpprb;
                for (i = 0; i < rttbl->NumSsc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->NumSdc; j++)
                    {
                        tmpval +=
                            (tmpconc[j] + gamma[j]) * rttbl->Dependency[i][j];
                    }
                    tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                    tmpconc[i + rttbl->NumStc] = tmpval;
                }
                for (i = 0; i < rttbl->NumStc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
                    {
                        tmpval += rttbl->Totalconc[i][j] * pow(10, tmpconc[j]);
                    }
                    residue_t[i] = tmpval - chms->t_conc[i];
                    jcb[k][i] = (residue_t[i] - residue[i]) / tmpprb;
                }
                tmpconc[k] -= tmpprb;
            }
            for (i = 0; i < rttbl->NumStc; i++)
            {
                x_[i] = -residue[i];
            }
            if (denseGETRF(jcb, rttbl->NumStc, rttbl->NumStc, p) != 0)
            {
                PIHMprintf(VL_ERROR, "Speciation error.\n");
                PIHMexit(EXIT_FAILURE);
            }
            denseGETRS(jcb, rttbl->NumStc, p, x_);
            for (i = 0; i < rttbl->NumStc; i++)
            {
                tmpconc[i] += x_[i];
                error[i] = residue[i] / totconc[i];
            }
            maxerror = fabs(error[0]);
            for (i = 1; i < rttbl->NumStc; i++)
            {
                maxerror = MAX(fabs(error[i]), maxerror);
            }
        }

        free(p);
        free(x_);
    }
    for (i = 0; i < rttbl->NumSsc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->NumSdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * rttbl->Dependency[i][j];
        }
        tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
        tmpconc[i + rttbl->NumStc] = tmpval;
    }
    for (i = 0; i < rttbl->NumStc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
        {
            tmpval += rttbl->Totalconc[i][j] * pow(10, tmpconc[j]);
        }
        totconc[i] = tmpval;
        residue[i] = tmpval - chms->t_conc[i];
        error[i] = residue[i] / totconc[i];
    }
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        if (i < rttbl->NumStc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                chms->p_conc[i] = pow(10, tmpconc[i]);
                chms->p_actv[i] = 1.0;
            }
            else
            {
                chms->p_conc[i] = pow(10, tmpconc[i]);
                chms->p_actv[i] = pow(10, (tmpconc[i] + gamma[i]));
            }
        }
        else
        {
            chms->s_conc[i - rttbl->NumStc] = pow(10, tmpconc[i]);
#if TEMP_DISABLED
            chms->s_actv[i - rttbl->NumStc] = pow(10, (tmpconc[i] + gamma[i]));
#endif
        }
    }
    destroyMat(jcb);

    return (0);
}
