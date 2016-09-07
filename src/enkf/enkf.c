#include "pihm.h"

void EnKF (double *xa, double obs, double obs_err, double *xf, int ne)
{
    int             i;
    double          y_hxm;
    double          km;
    double         *hxa;
    double         *xp;
    double          var;
    double          fac;
    double          s;
    double          beta;
    double          fac2;

    hxa = (double *)malloc (ne * sizeof (double));
    xp = (double *)malloc (ne * sizeof (double));

    for (i = 0; i < ne; i++)
    {
        xa[i] -= xa[ne];
    }

    y_hxm = obs - xf[ne];

    var = 0;

    for (i = 0; i < ne; i++)
    {
        hxa[i] = xf[i] - xf[ne];
        var += hxa[i] * hxa[i];
    }

    fac = 1.0 / ((double)ne - 1.0);
    s = fac * var + obs_err * obs_err;
    beta = 1.0 / (1.0 + sqrt ((obs_err * obs_err) / s));
    fac2 = fac / s;

    km = 0.0;

    for (i = 0; i < ne; i++)
    {
        xp[i] = xa[i];
        km += xp[i] * hxa[i];
    }

    km *= fac2;

    xa[ne] += km * y_hxm;

    km *= beta;

    for (i = 0; i < ne; i++)
    {
        xa[i] += km * (0.0 - hxa[i]);
        xa[i] += xa[ne];
    }

    free (hxa);
    free (xp);
}
