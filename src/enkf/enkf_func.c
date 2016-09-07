#include "pihm.h"

void PauseParal (int id)
{
    int             ii = 0;
    char            hostname[256];

    gethostname (hostname, sizeof (hostname));

    PIHMprintf (VL_NORMAL, "PID %d (%d) on %s ready for attach\n",
        getpid (), id, hostname);

    while (0 == ii)
    {
        sleep (5);
    }
}

void JobHandout (int starttime, int endtime, int startmode,
    ensmbr_struct *member, double *param, int ne, int total_jobs)
{
    int             dest;
    int             msg[3];
    int             i, j;

    msg[0] = starttime;
    msg[1] = endtime;
    msg[2] = startmode;

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXPARAM; j++)
        {
            param[i * MAXPARAM + j] = member[i].param[j];
        }
    }

    for (dest = 1; dest < total_jobs + 1; dest++)
    {
        MPI_Send (msg, 3, MPI_INT, dest, CYCLE_TAG, MPI_COMM_WORLD);
        MPI_Send (param, ne * MAXPARAM, MPI_DOUBLE, dest, PARAM_TAG,
            MPI_COMM_WORLD);
    }
}

void JobRecv (int *starttime, int *endtime, int *startmode, double *param,
    int ne)
{
    int             msg[3];
    MPI_Status      status;

    MPI_Recv (msg, 3, MPI_INT, 0, CYCLE_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv (param, ne * MAXPARAM, MPI_DOUBLE, 0, PARAM_TAG, MPI_COMM_WORLD,
        &status);

    *starttime = msg[0];
    *endtime = msg[1];
    *startmode = msg[2];
}

void JobHandIn (int total_jobs)
{
    int             ierr;
    int             success;
    int             received = 0;
    int             source;

    MPI_Status      status;

    while (received < total_jobs)
    {
        ierr =
            MPI_Recv (&success, 1, MPI_INT, MPI_ANY_SOURCE, SUCCESS_TAG,
            MPI_COMM_WORLD, &status);
        received++;
        source = status.MPI_SOURCE;
    }
}

void FreeEns (enkf_struct ens)
{
    int             i, j;

    for (i = 0; i < ens->nobs; i++)
    {
        free (ens->obs[i].var_ind);
        free (ens->obs[i].weight);
        for (j = 0; j < ens->obs[i].dim; j++)
        {
            free (ens->obs[i].k[j]);
            free (ens->obs[i].b[j]);
        }
        free (ens->obs[i].k);
        free (ens->obs[i].b);
    }

    free (ens->obs);

    for (i = 0; i < ens->ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                free (ens->member[i].var[j]);
            }
        }
    }

    free (ens->member);
}

void MapVar (var_struct *var, int numele, int numriv)
{
    int             i, k;
    int             n = 0;

    for (i = 0; i < MAXVAR; i++)
    {
        var[i].dim = 0;
    }

    for (i = 0; i < MAXVAR; i++)
    {
        switch (i)
        {
            case SURF_CTRL:
                strcpy (var[n].name, "surf");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case UNSAT_CTRL:
                strcpy (var[n].name, "unsat");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case GW_CTRL:
                strcpy (var[n].name, "gw");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case RIVSTG_CTRL:
                strcpy (var[n].name, "stage");
                var[n].dim = numriv;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case RIVGW_CTRL:
                strcpy (var[n].name, "rivgw");
                var[n].dim = numriv;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case SNOW_CTRL:
                strcpy (var[n].name, "snow");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = BADVAL;
                n++;
                break;
            case CMC_CTRL:
                strcpy (var[n].name, "is");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = BADVAL;
                n++;
                break;
            case INFIL_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RECHARGE_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case EC_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ETT_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case EDIR_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX0_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX1_CTRL:
                strcpy (var[n].name, "rivflx1");
                var[n].dim = numriv;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case RIVFLX2_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX3_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX4_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX5_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX6_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX7_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX8_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX9_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX10_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SUBFLX_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SURFFLX_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case T1_CTRL:
                strcpy (var[n].name, "t1");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case STC_CTRL:
                for (k = 0; k < MAXLYR; k++)
                {
                    sprintf (var[n].name, "stc%d", k);
                    var[n].dim = numele;
                    var[n].min = BADVAL;
                    var[n].max = BADVAL;
                    n++;
                }
                break;
            case SMC_CTRL:
                for (k = 0; k < MAXLYR; k++)
                {
                    sprintf (var[n].name, "smc%d", k);
                    var[n].dim = numele;
                    var[n].min = BADVAL;
                    var[n].max = BADVAL;
                    n++;
                }
                break;
            case SH2O_CTRL:
                for (k = 0; k < MAXLYR; k++)
                {
                    sprintf (var[n].name, "swc%d", k);
                    var[n].dim = numele;
                    var[n].min = BADVAL;
                    var[n].max = BADVAL;
                    n++;
                }
                break;
            case SNOWH_CTRL:
                strcpy (var[n].name, "snowh");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = BADVAL;
                n++;
                break;
            case ALBEDO_CTRL:
                strcpy (var[n].name, "albedo");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = 1.0;
                n++;
                break;
            case LE_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SH_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case G_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ETP_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ESNOW_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ROOTW_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SOILM_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SOLAR_CTRL:
                var[n].dim = 0;
                n++;
                break;
            default:
                break;
        }

    }
}

void Calib2Mbr (calib_struct cal, double *param)
{
    param[KSATH] = cal.ksath;
    param[KSATV] = cal.ksatv;
    param[KINF] = cal.kinfv;
    param[KMACH] = cal.kmach;
    param[KMACV] = cal.kmacv;
    param[DINF] = cal.dinf;
    param[RZD] = cal.rzd;
    param[DMAC] = cal.dmac;
    param[POROSITY] = cal.porosity;
    param[ALPHA] = cal.alpha;
    param[BETA] = cal.beta;
    param[AREAFV] = cal.areafv;
    param[AREAFH] = cal.areafh;
    param[VEGFRAC] = cal.vegfrac;
    param[ALBEDO] = cal.albedo;
    param[ROUGH] = cal.rough;
    param[PRCP] = cal.prcp;
    param[SFCTMP] = cal.sfctmp;
    param[EC] = cal.ec;
    param[ETT] = cal.ett;
    param[EDIR] = cal.edir;
    param[RIVROUGH] = cal.rivrough;
    param[RIVKSATH] = cal.rivksath;
    param[RIVKSATV] = cal.rivksatv;
    param[RIVBEDTHICK] = cal.rivbedthick;
    param[RIVDEPTH] = cal.rivdepth;
    param[RIVSHPCOEFF] = cal.rivshpcoeff;
#ifdef _NOAH_
    param[THETAREF] = cal.smcref;
    param[THETAW] = cal.smcwlt;
    param[RSMIN] = cal.rsmin;
    param[DRIP] = cal.drip;
    param[INTCP] = cal.cmcmax;
    param[CZIL] = cal.czil;
    param[FXEXP] = cal.fxexp;
    param[CFACTR] = cal.cfactr;
    param[RGL] = cal.rgl;
    param[HS] = cal.hs;
#endif
}

void Mbr2Cal (calib_struct *cal, const double *param)
{
    cal->ksath = param[KSATH];
    cal->ksatv = param[KSATV];
    cal->kinfv = param[KINF];
    cal->kmach = param[KMACH];
    cal->kmacv = param[KMACV];
    cal->dinf = param[DINF];
    cal->rzd = param[RZD];
    cal->dmac = param[DMAC];
    cal->porosity = param[POROSITY];
    cal->alpha = param[ALPHA];
    cal->beta = param[BETA];
    cal->areafv = param[AREAFV];
    cal->areafh = param[AREAFH];
    cal->vegfrac = param[VEGFRAC];
    cal->albedo = param[ALBEDO];
    cal->rough = param[ROUGH];
    cal->prcp = param[PRCP];
    cal->sfctmp = param[SFCTMP];
    cal->ec = param[EC];
    cal->ett = param[ETT];
    cal->edir = param[EDIR];
    cal->rivrough = param[RIVROUGH];
    cal->rivksath = param[RIVKSATH];
    cal->rivksatv = param[RIVKSATV];
    cal->rivbedthick = param[RIVBEDTHICK];
    cal->rivdepth = param[RIVDEPTH];
    cal->rivshpcoeff = param[RIVSHPCOEFF];
#ifdef _NOAH_
    cal->smcref = param[THETAREF];
    cal->smcwlt = param[THETAW];
    cal->rsmin = param[RSMIN];
    cal->drip = param[DRIP];
    cal->cmcmax = param[INTCP];
    cal->czil = param[CZIL];
    cal->fxexp = param[FXEXP];
    cal->cfactr = param[CFACTR];
    cal->rgl = param[RGL];
    cal->hs = param[HS];
#endif
}

double Randn ()
{
    double          temp1, temp2;
    double          x;

    temp1 = (double)(rand () % MAXINT + 1) / MAXINT;
    temp2 = (double)(rand () % MAXINT + 1) / MAXINT;

    x = sqrt (-2.0 * log (temp1)) * cos (2.0 * PI * temp2);

    return (x);
}

void GenRandNum (int ne, int nparam, double **randnum, double lower,
    double upper)
{
    int             i, j, k;
    double          corr[MAXPARAM][MAXPARAM];
    double          mean[MAXPARAM];
    double          std;
    int             corr_flag;
    int             std_flag = 1;
    double          s1, s2;
    double          max = -999.0;

    if (nparam > 1)
    {
        corr_flag = 1;
    }
    else
    {
        corr_flag = 0;
    }

    srand (time (NULL));

    do
    {
        max = -999.0;

        for (j = 0; j < nparam; j++)
        {
            std_flag = 1;

            while (std_flag)
            {
                mean[j] = 0.0;
                std = 0.0;

                for (i = 0; i < ne; i++)
                {
                    randnum[j][i] = Randn ();
                    mean[j] += randnum[j][i];
                }

                mean[j] /= (double)ne;

                for (i = 0; i < ne; i++)
                {
                    std +=
                        (randnum[j][i] - mean[j]) * (randnum[j][i] - mean[j]);
                }
                std = sqrt (std / ((double)ne - 1.0));

                std_flag = 0;

                for (i = 0; i < ne; i++)
                {
                    randnum[j][i] = (randnum[j][i] - mean[j]) / std * 1.0;

                    if (randnum[j][i] <= lower || randnum[j][i] >= upper)
                    {
                        std_flag = 1;
                        break;
                    }
                }
            }
        }

        for (i = 0; i < nparam; i++)
        {
            for (j = 0; j < nparam; j++)
            {
                corr[i][j] = 0.0;
                s1 = 0.0;
                s2 = 0.0;

                for (k = 0; k < ne; k++)
                {
                    corr[i][j] =
                        (randnum[i][k] - mean[i]) * (randnum[j][k] -
                        mean[j]) + corr[i][j];
                    s1 = s1 + (randnum[i][k] - mean[i]) * (randnum[i][k] -
                        mean[i]);
                    s2 = s2 + (randnum[j][k] - mean[j]) * (randnum[j][k] -
                        mean[j]);
                }

                corr[i][j] = corr[i][j] / sqrt (s1 * s2);

                if (fabs (corr[i][j]) > max && i != j)
                {
                    max = corr[i][j];
                }
            }
        }

        if (max < CORRMAX)
        {
            corr_flag = 0;
        }
    } while (corr_flag);
}
