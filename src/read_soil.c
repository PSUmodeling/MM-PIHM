#include "pihm.h"

void ReadSoil(const char *filename, soiltbl_struct *soiltbl)
{
    FILE           *soil_file;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             texture;
    const int       TOPSOIL = 1;
    const int       SUBSOIL = 0;
    int             ptf_used = 0;
    int             lno = 0;

    soil_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    /* Start reading soil file */
    NextLine(soil_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMSOIL", &soiltbl->number, 'i', filename, lno);

    soiltbl->silt = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->clay = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->om = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->bd = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->kinfv = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->ksatv = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->ksath = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->smcmax = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->smcmin = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->qtz = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->alpha = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->beta = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->areafh = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->areafv = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->dmac = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->smcref = (double *)malloc(soiltbl->number * sizeof(double));
    soiltbl->smcwlt = (double *)malloc(soiltbl->number * sizeof(double));

    /* Skip header line */
    NextLine(soil_file, cmdstr, &lno);

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine(soil_file, cmdstr, &lno);
        match = sscanf(cmdstr,
            "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index, &soiltbl->silt[i], &soiltbl->clay[i], &soiltbl->om[i],
            &soiltbl->bd[i],
            &soiltbl->kinfv[i], &soiltbl->ksatv[i], &soiltbl->ksath[i],
            &soiltbl->smcmax[i], &soiltbl->smcmin[i],
            &soiltbl->alpha[i], &soiltbl->beta[i],
            &soiltbl->areafh[i], &soiltbl->areafv[i],
            &soiltbl->dmac[i], &soiltbl->qtz[i]);

        if (match != 16 || i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading properties of the %dth soil type.\n", i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }

        /* Fill in missing organic matter and bulk density values */
        soiltbl->om[i] = (soiltbl->om[i] > 0.0) ? soiltbl->om[i] : 2.5;
        soiltbl->bd[i] = (soiltbl->bd[i] > 0.0) ? soiltbl->bd[i] : 1.3;

        /* Fill missing hydraulic properties using PTFs */
        if (soiltbl->kinfv[i] < 0.0)
        {
            soiltbl->kinfv[i] = PtfKv(soiltbl->silt[i], soiltbl->clay[i],
                soiltbl->om[i], soiltbl->bd[i], TOPSOIL);
            ptf_used = 1;
        }
        if (soiltbl->ksatv[i] < 0.0)
        {
            soiltbl->ksatv[i] = PtfKv(soiltbl->silt[i], soiltbl->clay[i],
                soiltbl->om[i], soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->ksath[i] < 0.0)
        {
            soiltbl->ksath[i] = 10.0 * soiltbl->ksatv[i];
            ptf_used = 1;
        }
        if (soiltbl->smcmax[i] < 0.0)
        {
            soiltbl->smcmax[i] = PtfThetas(soiltbl->silt[i], soiltbl->clay[i],
                soiltbl->om[i], soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->smcmin[i] < 0.0)
        {
            soiltbl->smcmin[i] = PtfThetar(soiltbl->silt[i], soiltbl->clay[i]);
            ptf_used = 1;
        }
        if (soiltbl->alpha[i] < 0.0)
        {
            soiltbl->alpha[i] = PtfAlpha(soiltbl->silt[i], soiltbl->clay[i],
                soiltbl->om[i], soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->beta[i] < 0.0)
        {
            soiltbl->beta[i] = PtfBeta(soiltbl->silt[i], soiltbl->clay[i],
                soiltbl->om[i], soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->qtz[i] < 0.0)
        {
            texture = SoilTex(soiltbl->silt[i], soiltbl->clay[i]);
            soiltbl->qtz[i] = Qtz(texture);
            ptf_used = 1;
        }

        /* Calculate field capacity and wilting point */
        soiltbl->smcref[i] = FieldCapacity(soiltbl->beta[i], soiltbl->ksatv[i],
            soiltbl->smcmax[i], soiltbl->smcmin[i]);
        soiltbl->smcwlt[i] = WiltingPoint(soiltbl->smcmax[i],
            soiltbl->smcmin[i], soiltbl->alpha[i], soiltbl->beta[i]);
    }

    NextLine(soil_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "DINF", &soiltbl->dinf, 'd', filename, lno);

    NextLine(soil_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACV_RO", &soiltbl->kmacv_ro, 'd', filename, lno);

    NextLine(soil_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACH_RO", &soiltbl->kmach_ro, 'd', filename, lno);

    if (ptf_used)
    {
        pihm_printf(VL_NORMAL, "\nA priori soil parameter values (uncalibrated) "
                              "estimated from pedotransfer functions:\n");
        pihm_printf(VL_NORMAL,
            "%-7s\t%-15s\t%-15s\t%-15s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\n",
            "TYPE", "KINFV", "KSATV", "KSATH", "SMCMAX", "SMCMIN", "ALPHA",
            "BETA", "QTZ");
        pihm_printf(VL_NORMAL,
            "%-7s\t%-15s\t%-15s\t%-15s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\n",
            "-", "m s-1", "m s-1", "m s-1", "m3 m-3", "m3 m-3", "m-1",
            "-", "-");
        for (i = 0; i < soiltbl->number; i++)
        {
            pihm_printf(VL_NORMAL,
                "%-7d\t%-15.3le\t%-15.3le\t%-15.3le\t%-7.3lf\t%-7.3lf\t"
                "%-7.3lf\t%-7.3lf\t%-7.3lf\n",
                i + 1, soiltbl->kinfv[i], soiltbl->ksatv[i],
                soiltbl->ksath[i], soiltbl->smcmax[i], soiltbl->smcmin[i],
                soiltbl->alpha[i], soiltbl->beta[i], soiltbl->qtz[i]);
        }
    }

    fclose(soil_file);
}
