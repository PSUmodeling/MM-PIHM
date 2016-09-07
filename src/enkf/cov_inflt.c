#include "pihm.h"

void CovInflt (enkf_struct ens, enkf_struct ens0)
{
    int             i, j, k;
    int             ne;

    double         *xp0, *xp;

    double        **x, *x0;
    double          average0, average;
    double          param_min, param_max;
    double          param_std;
    double          c1, c2, c;

    //double TotalWaterVolume[ens->ne];
    //double MassCoeff[ens->ne];

    ne = ens->ne;

    xp0 = (double *)malloc (sizeof (double) * ne);
    xp = (double *)malloc (sizeof (double) * ne);

    x = (double **)malloc (ne * sizeof (double));
    x0 = (double *)malloc (ne * sizeof (double));

    PIHMprintf (VL_NORMAL, "\n*****Parameters********\n");

    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].update == 1 && ens->update_param == 1)
        {
            PIHMprintf (VL_NORMAL, "%s\n", ens->param[i].name);
            /* Make pointers point to the parameter that needs to be updated */
            for (j = 0; j < ne; j++)
            {
                x[j] = &ens->member[j].param[i];
                x0[j] = ens0->member[j].param[i];
            }

            /* Take log if the parameter is conductivity or Czil */
            if (ens->param[i].type == LOG_TYPE)
            {
                for (j = 0; j < ne; j++)
                {
                    *x[j] = log10 (*x[j]);
                    x0[j] = log10 (x0[j]);
                }
            }

            average0 = 0.0;
            average = 0.0;

            for (j = 0; j < ne; j++)
            {
                average0 += x0[j];
                average += *x[j];
            }

            average /= (double)ne;
            average0 /= (double)ne;

            for (j = 0; j < ne; j++)
            {
                xp0[j] = x0[j] - average0;
                xp[j] = *x[j] - average;
            }

            if (average < ens->param[i].max - 0.25 * ens->param[i].init_std &&
                average > ens->param[i].min + 0.25 * ens->param[i].init_std)
            {
                /* Calculate new perturbations and standard deviation */
                param_std = 0.0;

                for (j = 0; j < ne; j++)
                {
                    xp[j] =
                        (1.0 - ens->weight) * xp[j] + ens->weight * xp0[j];
                    param_std += xp[j] * xp[j];
                }
                param_std = sqrt (param_std / ((double)ne - 1.0));

                /* Coavariance inflation */
                param_min = 999.0;
                param_max = -999.0;

                for (j = 0; j < ne; j++)
                {
                    if (param_std < 0.25 * ens->param[i].init_std)
                    {
                        xp[j] =
                            0.25 * ens->param[i].init_std / param_std * xp[j];
                    }
                    param_min =
                        (param_min <
                        average + xp[j]) ? param_min : (average + xp[j]);
                    param_max =
                        (param_max >
                        average + xp[j]) ? param_max : (average + xp[j]);
                }

                c1 = (ens->param[i].max - 1.0E-6 - average) / (param_max -
                    average);
                c2 = (average - ens->param[i].min - 1.0E-6) / (average -
                    param_min);
                c = (c1 < c2) ? c1 : c2;
                c = (c < 1.0) ? c : 1.0;

                for (j = 0; j < ne; j++)
                {
                    *x[j] = average + c * xp[j];
                    if (ens->param[i].type == LOG_TYPE)
                    {
                        *x[j] = pow (10.0, *x[j]);
                    }
                    PIHMprintf (VL_NORMAL, "%lf\t", *x[j]);
                }
                if (ens->param[i].type == LOG_TYPE)
                {
                    average = pow (10.0, average);
                }
                PIHMprintf (VL_NORMAL, "mean: %lf\n", average);
            }
            else
            {
                PIHMprintf (VL_NORMAL,
                    "EnKF analysis %lf is out of range. "
                    "Parameter is not updated\n", average);
                average = average0;
                for (j = 0; j < ne; j++)
                {
                    *x[j] = x0[j];
                    if (ens->param[i].type == LOG_TYPE)
                    {
                        *x[j] = pow (10.0, x0[j]);
                    }
                }

                if (ens->param[i].type == LOG_TYPE)
                {
                    average = pow (10.0, average);
                }
            }
        }
    }

    if (ens->update_var == 1)
    {
        for (i = 0; i < MAXVAR; i++)
        {
            if (ens->var[i].dim > 0)
            {
                for (k = 0; k < ens->var[i].dim; k++)
                {
                    for (j = 0; j < ne; j++)
                    {
                        x[j] = &ens->member[j].var[i][k];
                        x0[j] = ens0->member[j].var[i][k];
                    }

                    average = 0.0;
                    average0 = 0.0;

                    for (j = 0; j < ne; j++)
                    {
                        average += *x[j];
                        average0 += x0[j];
                    }
                    average /= (double)ne;
                    average0 /= (double)ne;

                    for (j = 0; j < ne; j++)
                    {
                        xp0[j] = x0[j] - average0;
                        xp[j] = *x[j] - average;
                        *x[j] =
                            average + (1.0 - ens->weight) * xp[j] +
                            ens->weight * xp0[j];
                        if ((int)ens->var[i].max != BADVAL)
                        {
                            *x[j] =
                                (*x[j] <
                                ens->var[i].max) ? *x[j] : ens->var[i].max;
                        }

                        if ((int)ens->var[i].min != BADVAL)
                        {
                            *x[j] =
                                (*x[j] >
                                ens->var[i].min) ? *x[j] : ens->var[i].min;
                        }
                    }
                }
            }
        }
    }

    free (xp0);
    free (xp);
    free (x);
    free (x0);

    PIHMprintf (VL_NORMAL, "Inflation done!\n");
}
