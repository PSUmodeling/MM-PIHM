#ifndef FORCING_HEADER
#define FORCING_HEADER

typedef struct TSD_type
{
    char            name[15];
    int             index;
    int             length;     /* length of time series */
    int             iCounter;   /* interpolation counter */
    double          TSFactor;
    double        **TS;         /* 2D time series data */
} TSD;

double monthly_lai (double t, int LC_type);
double monthly_rl (double t, int LC_type);
double monthly_mf (double t);
double Interpolation (TSD * Data, double t);
void MultiInterpolation (TSD * Data, double t, double *forcing, int num_forcing);
void        update (double, void *);

#endif
