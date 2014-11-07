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

double      Interpolation (TSD *, double);
void        update (double, void *);

#endif

