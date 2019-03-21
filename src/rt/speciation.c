#include "pihm.h"

#define TOL        1E-7

int Speciation(chemtbl_struct chemtbl[], rttbl_struct *rttbl, Chem_Data CD, int cell,
    int speciation_flg)
{
    /* if speciation flg = 1, pH is defined
     * if speciation flg = 0, all defined value is total concentration */
    int             i, j, k, num_spe = rttbl->NumStc + rttbl->NumSsc;
    double         *residue, *residue_t, *tmpconc, *totconc;
    double          tmpval, tmpprb = 1E-2, I, Iroot;
    double         *error, *gamma, *Keq, *current_totconc, adh, bdh, bdt;
    realtype      **jcb;

    residue = (double *)calloc(rttbl->NumStc, sizeof(double));
    residue_t = (double *)calloc(rttbl->NumStc, sizeof(double));
    tmpconc = (double *)calloc(rttbl->NumStc + rttbl->NumSsc, sizeof(double));
    totconc = (double *)calloc(rttbl->NumStc, sizeof(double));
    error = (double *)calloc(rttbl->NumStc, sizeof(double));
    gamma = (double *)calloc(num_spe, sizeof(double));
    Keq = (double *)calloc(rttbl->NumSsc, sizeof(double));
    current_totconc = (double *)calloc(rttbl->NumStc, sizeof(double));

    if (rttbl->TEMcpl == 0)
    {
        for (i = 0; i < rttbl->NumSsc; i++)
            Keq[i] = rttbl->Keq[i];
    }

    adh = rttbl->adh;
    bdh = rttbl->bdh;
    bdt = rttbl->bdt;

    for (i = 0; i < rttbl->NumStc; i++)
    {
        /* Using log10 conc as the primary unknowns. Works better because
         * negative numbers are not a problem. */
        tmpconc[i] = log10(CD->Vcele[cell].p_conc[i]);
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

    /* Initial speciation to get secondary species, no activity corrections */
    for (i = 0; i < num_spe; i++)
        gamma[i] = 0;

    for (i = 0; i < rttbl->NumStc; i++)
    {
        current_totconc[i] = CD->Vcele[cell].t_conc[i];
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
                    I += 0.5 * pow(10, tmpconc[i]) *
                        chemtbl[i].Charge * chemtbl[i].Charge;
                Iroot = sqrt(I);
                for (i = 0; i < num_spe; i++)
                {
                    if (chemtbl[i].itype == MINERAL)
                        gamma[i] = -tmpconc[i];
                    /* aqueous species in the unit of mol/L, however the solids
                     * are in the unit of mol/L porous media
                     * the activity of solid is 1, the log10 of activity is 0.
                     * by assigning gamma[minerals] to negative of the
                     * tmpconc[minerals], we ensured the log 10 of activity of
                     * solids are 0*/
                    else
                        gamma[i] =
                            (-adh * chemtbl[i].Charge * chemtbl[i].Charge * Iroot) / (1 +
                            bdh * chemtbl[i].SizeF * Iroot) + bdt * I;
                    if (strcmp(chemtbl[i].ChemName, "'H+'") == 0)
                    {
                        tmpconc[i] =
                            log10(CD->Vcele[cell].p_actv[i]) - gamma[i];
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
                    CD->Vcele[cell].t_conc[i] = totconc[i];
                residue[i] = tmpval - CD->Vcele[cell].t_conc[i];
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
                            tmpval +=
                                (tmpconc[j] + gamma[j]) * rttbl->Dependency[i][j];
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
                                tmpval +=
                                    rttbl->Totalconc[i][j] * pow(10, tmpconc[j]);
                            }
                            residue_t[i] = tmpval - CD->Vcele[cell].t_conc[i];
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
                ReportError(chemtbl, rttbl, CD->Vcele[cell], CD);
                return (1);
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
                    I += 0.5 * pow(10,
                        tmpconc[i]) * chemtbl[i].Charge * chemtbl[i].Charge;
                }
                Iroot = sqrt(I);
                for (i = 0; i < num_spe; i++)
                {
                    if (chemtbl[i].itype == MINERAL)
                        gamma[i] = -tmpconc[i];
                    else
                        gamma[i] =
                            (-adh * chemtbl[i].Charge * chemtbl[i].Charge * Iroot) / (1 +
                            bdh * chemtbl[i].SizeF * Iroot) + bdt * I;
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
                residue[i] = tmpval - CD->Vcele[cell].t_conc[i];
            }

            for (k = 0; k < rttbl->NumStc; k++)
            {
                tmpconc[k] += tmpprb;
                for (i = 0; i < rttbl->NumSsc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->NumSdc; j++)
                        tmpval +=
                            (tmpconc[j] + gamma[j]) * rttbl->Dependency[i][j];
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
                    residue_t[i] = tmpval - CD->Vcele[cell].t_conc[i];
                    jcb[k][i] = (residue_t[i] - residue[i]) / (tmpprb);
                }
                tmpconc[k] -= tmpprb;
            }
            for (i = 0; i < rttbl->NumStc; i++)
                x_[i] = -residue[i];
            if (denseGETRF(jcb, rttbl->NumStc, rttbl->NumStc, p) != 0)
            {
                ReportError(chemtbl, rttbl, CD->Vcele[cell], CD);
                return (1);
            }
            denseGETRS(jcb, rttbl->NumStc, p, x_);
            for (i = 0; i < rttbl->NumStc; i++)
            {
                tmpconc[i] += x_[i];
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
        residue[i] = tmpval - CD->Vcele[cell].t_conc[i];
        error[i] = residue[i] / totconc[i];
    }
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        if (i < rttbl->NumStc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                CD->Vcele[cell].p_conc[i] = pow(10, tmpconc[i]);
                CD->Vcele[cell].p_actv[i] = 1.0;
            }
            else
            {
                CD->Vcele[cell].p_conc[i] = pow(10, tmpconc[i]);
                CD->Vcele[cell].p_actv[i] = pow(10, (tmpconc[i] + gamma[i]));
            }
        }
        else
        {
            CD->Vcele[cell].s_conc[i - rttbl->NumStc] = pow(10, tmpconc[i]);
#if TEMP_DISABLED
            CD->Vcele[cell].s_actv[i - rttbl->NumStc] =
                pow(10, (tmpconc[i] + gamma[i]));
#endif
        }
    }
    destroyMat(jcb);

    free(residue);
    free(residue_t);
    free(tmpconc);
    free(totconc);
    free(error);
    free(gamma);
    free(Keq);
    free(current_totconc);

    return (0);
}

