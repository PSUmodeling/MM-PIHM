#include "pihm.h"

#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80

void Lookup(FILE *database, const calib_struct *cal, chemtbl_struct chemtbl[],
    kintbl_struct kintbl[], rttbl_struct *rttbl)
{
    double          tmpval[WORDS_LINE];
    char            line[LINE_WIDTH];
    int             i, j, k, l;
    int             ind;
    int             keq_position = BADVAL;
    int             total_temp_points;
    int             mn, in;
    char          **tmpstr = (char **)malloc(WORDS_LINE * sizeof(char *));
    char            cmdstr[MAXSTRING];
    int             lno = 0;

    for (i = 0; i < WORDS_LINE; i++)
    {
        tmpstr[i] = (char *)malloc(WORD_WIDTH * sizeof(char));
        memset(tmpstr[i], 0, WORD_WIDTH);
    }

    /*
     * Find temperature point in database
     */
    NextLine(database, cmdstr, &lno);
    while(MatchWrappedKey(cmdstr, "'temperature points'") != 0)
    {
        NextLine(database, cmdstr, &lno);
    }
    ReadTempPoints(cmdstr, rttbl->Temperature, &total_temp_points,
        &keq_position);
    if (keq_position == BADVAL)
    {
        PIHMprintf(VL_ERROR, "Error reading temperature points "
            "in %s near Line %d", ".cdbs", lno);
        PIHMexit(EXIT_FAILURE);
    }

    /*
     * Read Debye Huckel parameters from database
     */
    rttbl->adh = BADVAL;
    rttbl->bdh = BADVAL;
    rttbl->bdt = BADVAL;

    while (MatchWrappedKey(cmdstr, "'Debye-Huckel adh'") != 0)
    {
        NextLine(database, cmdstr, &lno);
    }
    ReadDHParam(cmdstr, keq_position, &rttbl->adh);
    if (roundi(rttbl->adh) == BADVAL)
    {
        PIHMprintf(VL_ERROR, "Error reading Debye Huckel parameters "
            "in %s near Line %d", ".cdbs", lno);
        PIHMexit(EXIT_FAILURE);
    }

    while (MatchWrappedKey(cmdstr, "'Debye-Huckel bdh'") != 0)
    {
        NextLine(database, cmdstr, &lno);
    }
    ReadDHParam(cmdstr, keq_position, &rttbl->bdh);
    if (roundi(rttbl->bdh) == BADVAL)
    {
        PIHMprintf(VL_ERROR, "Error reading Debye Huckel parameters "
            "in %s near Line %d", ".cdbs", lno);
        PIHMexit(EXIT_FAILURE);
    }

    while (MatchWrappedKey(cmdstr, "'Debye-Huckel bdt'") != 0)
    {
        NextLine(database, cmdstr, &lno);
    }
    ReadDHParam(cmdstr, keq_position, &rttbl->bdt);
    if (roundi(rttbl->bdt) == BADVAL)
    {
        PIHMprintf(VL_ERROR, "Error reading Debye Huckel parameters "
            "in %s near Line %d", ".cdbs", lno);
        PIHMexit(EXIT_FAILURE);
    }

    PIHMprintf(VL_VERBOSE,
        " Debye-Huckel Parameters set to A=%6.4f; B=%6.4f; b=%6.4f\n\n",
        rttbl->adh, rttbl->bdh, rttbl->bdt);


    for (i = 0; i < MAXSPS; i++)
    {
        for (j = 0; j < MAXSPS; j++)
        {
            rttbl->Dependency[i][j] = 0.0;      /* NumSsc x NumSdc */
            rttbl->Dep_kinetic[i][j] = 0.0;     /* (NumMkr + NumAkr) x NumStc */
            rttbl->Dep_kinetic_all[i][j] = 0.0; /* NumMin x NumStc */
            rttbl->Totalconc[i][j] = 0.0;       /* NumStc x (NumStc + NumSsc) */
#if NOT_YET_IMPLEMENTED
            rttbl->Totalconck[i][j] = 0.0;      /* NumStc x (NumStc + NumSsc) */
#endif
        }
        /* Keqs of equilibrium: kinetic and kinetic all */
        rttbl->Keq[i] = 0.0;                    /* NumSsc */
        rttbl->KeqKinect[i] = 0.0;              /* NumMkr + NumAkr */
        rttbl->KeqKinect_all[i] = 0.0;          /* NumMin */
    }

    /*
     * Update species parameters
     */
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        if (strcmp(chemtbl[i].ChemName, "pH") == 0)
        {
            strcpy(chemtbl[i].ChemName, "H+");
        }
        wrap(chemtbl[i].ChemName);

        chemtbl[i].DiffCoe = rttbl->DiffCoe;
        chemtbl[i].DispCoe = rttbl->DispCoe;
        chemtbl[i].Charge = 0.0;
        chemtbl[i].SizeF = 1.0;
        chemtbl[i].MolarMass = BADVAL;
    }

    /*
     * Read parameters for primary species
     */
    while (MatchWrappedKey(cmdstr, "'End of primary'") != 0)
    {
        NextLine(database, cmdstr, &lno);
        ReadPrimary(cmdstr, rttbl->NumStc, chemtbl);
    }

    /*
     * Read parameters for secondary species
     */
    while (MatchWrappedKey(cmdstr, "'End of secondary'") != 0)
    {
        NextLine(database, cmdstr, &lno);
        ReadSecondary(cmdstr, total_temp_points, keq_position, chemtbl, rttbl);
    }

    /*
     * Read parameters for minerals
     */
    while (MatchWrappedKey(cmdstr, "'End of minerals'") != 0)
    {
        NextLine(database, cmdstr, &lno);
        ReadMinerals(cmdstr, total_temp_points, keq_position, chemtbl, rttbl);
    }

    /*
     * Build dependencies
     */
    for (i = 0; i < rttbl->NumMkr + rttbl->NumAkr; i++)
    {
        ind = kintbl[i].position - rttbl->NumStc + rttbl->NumMin;
        PIHMprintf(VL_VERBOSE,
            " Selecting the kinetic species %s from all possible species.\n\n",
            chemtbl[kintbl[i].position].ChemName);
        rttbl->KeqKinect[i] = rttbl->KeqKinect_all[ind];
        for (k = 0; k < rttbl->NumStc; k++)
        {
            rttbl->Dep_kinetic[i][k] = rttbl->Dep_kinetic_all[ind][k];
        }
    }

    while (strcmp(cmdstr, "End of surface complexation\r\n") != 0 &&
        strcmp(cmdstr, "End of surface complexation\n") != 0)
    {
        NextLine(database, cmdstr, &lno);
        ReadAdsorption(cmdstr, total_temp_points, keq_position, chemtbl, rttbl);
    }

    while (strcmp(cmdstr, "End of exchange\r\n") != 0 &&
        strcmp(cmdstr, "End of exchange\n") != 0)
    {
        NextLine(database, cmdstr, &lno);
        ReadCationEchg(cmdstr, cal->Xsorption, chemtbl, rttbl);
    }

    for (i = 0; i < MAXSPS; i++)
    {
        for (j = 0; j < MAXDEP; j++)
        {
            kintbl[i].dep_position[j] = 0;
            kintbl[i].monod_position[j] = 0;
            kintbl[i].inhib_position[j] = 0;
        }
    }

    for (i = 0; i < rttbl->NumMkr; i++)
    {
        FindLine(database, "BOF", &lno, ".cdbs");
        NextLine(database, cmdstr, &lno);
        while (strcmp(cmdstr, "Begin mineral kinetics\r\n") != 0 &&
            strcmp(cmdstr, "Begin mineral kinetics\n") != 0)
        {
            NextLine(database, cmdstr, &lno);
        }

        while (strcmp(cmdstr, "End of mineral kinetics\r\n") != 0 &&
            strcmp(cmdstr, "End of mineral kinetics\n") != 0)
        {
            ReadMinKin(database, rttbl->NumStc, cal->rate, &lno, cmdstr, chemtbl, &kintbl[i]);
            NextLine(database, cmdstr, &lno);
        }
    }

    for (i = 0; i < rttbl->NumStc; i++)
    {
        rttbl->Totalconc[i][i] = 1.0;
    }

    for (i = rttbl->NumStc; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        for (j = 0; j < rttbl->NumSdc; j++)
        {
            rttbl->Totalconc[j][i] += rttbl->Dependency[i - rttbl->NumStc][j];
        }
    }
    if (rttbl->NumSsc > 0)
    {
        PIHMprintf(VL_VERBOSE, " \n Dependency Matrix!\n");

        PIHMprintf(VL_VERBOSE, "%-15s", " ");
        for (i = 0; i < rttbl->NumSdc; i++)
        {
            PIHMprintf(VL_VERBOSE, "%-14s", chemtbl[i].ChemName);
        }
        PIHMprintf(VL_VERBOSE, "\n");

        for (i = 0; i < rttbl->NumSsc; i++)
        {
            PIHMprintf(VL_VERBOSE, " %-14s",
                chemtbl[i + rttbl->NumStc].ChemName);
            for (j = 0; j < rttbl->NumSdc; j++)
            {
                PIHMprintf(VL_VERBOSE, "%-14.2f", rttbl->Dependency[i][j]);
            }
            PIHMprintf(VL_VERBOSE, " %6.2f\n", rttbl->Keq[i]);
        }
    }

    PIHMprintf(VL_VERBOSE, " \n Total Concentration Matrix!\n");
    PIHMprintf(VL_VERBOSE, "%-18s", " ");
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        PIHMprintf(VL_VERBOSE, "%-14s", chemtbl[i].ChemName);
    }
    PIHMprintf(VL_VERBOSE, "\n");
    for (i = 0; i < rttbl->NumStc; i++)
    {
        PIHMprintf(VL_VERBOSE, " Sum%-14s", chemtbl[i].ChemName);
        for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
        {
            PIHMprintf(VL_VERBOSE, "%-14.2f", rttbl->Totalconc[i][j]);
        }
        PIHMprintf(VL_VERBOSE, "\n");
    }

    PIHMprintf(VL_VERBOSE, " \n Kinetic Mass Matrx!\n");
    PIHMprintf(VL_VERBOSE, "%-15s", " ");
    for (i = 0; i < rttbl->NumStc; i++)
    {
        PIHMprintf(VL_VERBOSE, "%-14s", chemtbl[i].ChemName);
    }
    PIHMprintf(VL_VERBOSE, "\n");
    for (j = 0; j < rttbl->NumMkr + rttbl->NumAkr; j++)
    {
        PIHMprintf(VL_VERBOSE, " %-14s",
            chemtbl[j + NumSpc + rttbl->NumAds + rttbl->NumCex].ChemName);
        for (i = 0; i < rttbl->NumStc; i++)
        {
            PIHMprintf(VL_VERBOSE, "%-14f", rttbl->Dep_kinetic[j][i]);
        }
        PIHMprintf(VL_VERBOSE, " Keq = %-6.2f\n", rttbl->KeqKinect[j]);
    }

#if NOT_YET_IMPLEMENTED
    /* Use calibration coefficient to produce new Keq values for
     * 1) CO2, 2) other kinetic reaction */
    double          Cal_PCO2 = 1.0;
    double          Cal_Keq = 1.0;
    for (i = 0; i < rttbl->NumAkr + rttbl->NumMkr; i++)
    {
        rttbl->KeqKinect[i] +=
            (!strcmp(chemtbl[i + NumSpc + rttbl->NumAds + rttbl->NumCex].ChemName,
            "'CO2(*g)'")) ? log10(Cal_PCO2) : log10(Cal_Keq);
    }

    PIHMprintf(VL_VERBOSE, "\n Kinetic Mass Matrx (calibrated Keq)! \n");
    PIHMprintf(VL_VERBOSE, "%-15s", " ");
    for (i = 0; i < rttbl->NumStc; i++)
        PIHMprintf(VL_VERBOSE, "%-14s", chemtbl[i].ChemName);
    PIHMprintf(VL_VERBOSE, "\n");
    for (j = 0; j < rttbl->NumMkr + rttbl->NumAkr; j++)
    {
        PIHMprintf(VL_VERBOSE, " %-14s",
            chemtbl[j + NumSpc + rttbl->NumAds + rttbl->NumCex].ChemName);
        for (i = 0; i < rttbl->NumStc; i++)
        {
            PIHMprintf(VL_VERBOSE, "%-14.2f", rttbl->Dep_kinetic[j][i]);
        }
        PIHMprintf(VL_VERBOSE, " Keq = %-6.2f\n", rttbl->KeqKinect[j]);
    }
    PIHMprintf(VL_VERBOSE, "\n");
#endif

    PIHMprintf(VL_VERBOSE,
        " \n Mass action species type determination "
        "(0: immobile, 1: mobile, 2: Mixed) \n");
    for (i = 0; i < NumSpc; i++)
    {
        chemtbl[i].mtype = (chemtbl[i].itype == AQUEOUS) ?
            MOBILE_MA : IMMOBILE_MA;

        for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
        {
            chemtbl[i].mtype = (rttbl->Totalconc[i][j] != 0 &&
                chemtbl[j].itype != chemtbl[i].mtype) ?
                MIXED_MA : chemtbl[i].mtype;
        }
        PIHMprintf(VL_VERBOSE, " %12s\t%10d\n",
            chemtbl[i].ChemName, chemtbl[i].mtype);
    }

    PIHMprintf(VL_VERBOSE,
        " \n Individual species type determination "
        "(1: aqueous, 2: adsorption, 3: ion exchange, 4: solid)\n");
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        PIHMprintf(VL_VERBOSE, " %12s\t%10d\n",
            chemtbl[i].ChemName, chemtbl[i].itype);
    }

    for (i = 0; i < WORDS_LINE; i++)
    {
        free(tmpstr[i]);
    }
    free(tmpstr);
}

void wrap(char *str)
{
    char            word[WORD_WIDTH];
    sprintf(word, "'%s'", str);
    strcpy(str, word);
}

int MatchWrappedKey(const char cmdstr[], const char key[])
{
    char            optstr[MAXSTRING];

    sscanf(cmdstr, "'%[^']'", optstr);
    wrap(optstr);

    return (strcmp(optstr, key) == 0) ? 0 : 1;
}

void ReadTempPoints(const char cmdstr[], double tmp, int *total_points,
    int *keq_position)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i;
    double          val;

    if (sscanf(cmdstr + bytes_consumed, "'%*[^']'%n", &bytes_now) != 0)
    {
        return;
    }
    bytes_consumed += bytes_now;

    if (sscanf(cmdstr + bytes_consumed, "%d%n", total_points, &bytes_now) != 1)
    {
        return;
    }
    bytes_consumed += bytes_now;

    for (i = 0; i < *total_points; i++)
    {
        if (sscanf(cmdstr + bytes_consumed, "%lf%n", &val, &bytes_now) != 1)
        {
            return;
        }
        bytes_consumed += bytes_now;

        if (fabs(val - tmp) < 1.0E-5)
        {
            *keq_position = i + 1;
            PIHMprintf(VL_VERBOSE,
                "\n Temperature point %.2f found in database (Pos. %d).\n\n",
                    val, *keq_position);
            return;
        }
    }
}

void ReadDHParam(const char cmdstr[], int tmp_position, double *param)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i;

    if (sscanf(cmdstr + bytes_consumed, "'%*[^']'%n", &bytes_now) != 0)
    {
        return;
    }
    bytes_consumed += bytes_now;

    for (i = 0; i < tmp_position; i++)
    {
        if (sscanf(cmdstr + bytes_consumed, "%lf%n", param, &bytes_now) != 1)
        {
            return;
        }
        bytes_consumed += bytes_now;
    }
}

void ReadPrimary(const char cmdstr[], int NumStc, chemtbl_struct chemtbl[])
{
    int             i;

    for (i = 0; i < NumStc; i++)
    {
        if (MatchWrappedKey(cmdstr, chemtbl[i].ChemName) == 0)
        {
            if (sscanf(cmdstr, "'%*[^']' %lf %lf %lf", &chemtbl[i].SizeF,
                &chemtbl[i].Charge, &chemtbl[i].MolarMass) != 3)
            {
                PIHMprintf(VL_ERROR,
                    "Error reading primary species parameters for %s\n",
                    chemtbl[i].ChemName);
                PIHMexit(EXIT_FAILURE);
            }
            PIHMprintf(VL_VERBOSE,
                " Primary species %s found in database.\n"
                " MolarMass = %6.4f\n\n",
                chemtbl[i].ChemName, chemtbl[i].MolarMass);
            break;
        }
    }
}

void ReadSecondary(const char cmdstr[], int npoints, int keq_position,
    chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k;
    int             ind;
    int             ndep;
    double          dep;
    char            chemn[MAXSTRING];

    sscanf(cmdstr + bytes_consumed, "'%[^']' %d%n", chemn, &ndep, &bytes_now);
    bytes_consumed += bytes_now;
    wrap(chemn);

    for (i = 0; i < rttbl->NumSsc; i++)
    {
        ind = i + rttbl->NumStc;

        if (MatchWrappedKey(chemn, chemtbl[ind].ChemName) == 0)
        {
            PIHMprintf(VL_VERBOSE,
                " Secondary species %s found in database!\n",
                chemtbl[ind].ChemName);
            PIHMprintf(VL_VERBOSE, " %s", cmdstr);
            chemtbl[ind].itype = AQUEOUS;

            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn,
                    &bytes_now);
                bytes_consumed += bytes_now;
                wrap(chemn);

                for (k = 0; k < rttbl->NumSdc; k++)
                {
                    if (MatchWrappedKey(chemtbl[k].ChemName, chemn) == 0)
                    {
                        rttbl->Dependency[i][k] = dep;
                    }
                }
            }

            for (j = 0; j < npoints; j++)
            {
                if (j + 1 == keq_position)
                {
                    sscanf(cmdstr + bytes_consumed, "%lf%n", &rttbl->Keq[i],
                        &bytes_now);
                    bytes_consumed += bytes_now;
                    PIHMprintf(VL_VERBOSE,
                        " Keq = %6.4f\n", rttbl->Keq[i]);
                }
                else
                {
                    sscanf(cmdstr + bytes_consumed, "%*f%n", &bytes_now);
                    bytes_consumed += bytes_now;
                }
            }

            sscanf(cmdstr + bytes_consumed, "%lf %lf %lf", &chemtbl[ind].SizeF,
                &chemtbl[ind].Charge, &chemtbl[ind].MolarMass);
            PIHMprintf(VL_VERBOSE,
                " MolarMass = %6.4f, Charge = %6.4f,"
                " SizeFactor = %6.4f\n\n",
                chemtbl[ind].MolarMass, chemtbl[ind].Charge,
                chemtbl[ind].SizeF);

            break;
        }
    }
}

void ReadMinerals(const char cmdstr[], int npoints, int keq_position,
    chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k, l;
    int             ind;
    int             ndep;
    double          dep;
    double          mvol;
    char            chemn[MAXSTRING];

    sscanf(cmdstr + bytes_consumed, "'%[^']' %lf %d%n", chemn, &mvol, &ndep,
        &bytes_now);
    bytes_consumed += bytes_now;
    wrap(chemn);

    for (i = 0; i < rttbl->NumMin; i++)
    {
        ind = i + NumSpc + rttbl->NumAds + rttbl->NumCex;

        if (MatchWrappedKey(chemn, chemtbl[ind].ChemName) == 0)
        {
            PIHMprintf(VL_VERBOSE, " Mineral %s found in database!\n",
                chemtbl[ind].ChemName);
            PIHMprintf(VL_VERBOSE, " %s", cmdstr);

            chemtbl[ind].MolarVolume = mvol;
            chemtbl[ind].itype = MINERAL;

            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn,
                    &bytes_now);
                bytes_consumed += bytes_now;
                wrap(chemn);

                for (k = 0; k < rttbl->NumStc + rttbl->NumSsc; k++)
                {
                    if (MatchWrappedKey(chemtbl[k].ChemName, chemn) == 0)
                    {
                        if (k < rttbl->NumStc)
                        {
                            rttbl->Dep_kinetic_all[i][k] = dep;
                        }
                        else
                        {
                            for (l = 0; l < NumSpc; l++)
                            {
                                rttbl->Dep_kinetic_all[i][l] += dep *
                                    rttbl-> Dependency[k - rttbl->NumStc][l];
                            }
                            rttbl->KeqKinect_all[i] += dep *
                                rttbl->Keq[k - rttbl->NumStc];
                        }

                        break;
                    }
                }
            }

            for (j = 0; j < npoints; j++)
            {
                if (j + 1 == keq_position)
                {
                    sscanf(cmdstr + bytes_consumed, "%lf%n",
                        &rttbl->KeqKinect_all[i], &bytes_now);
                    bytes_consumed += bytes_now;
                }
                else
                {
                    sscanf(cmdstr + bytes_consumed, "%*f%n", &bytes_now);
                    bytes_consumed += bytes_now;
                }
            }

            sscanf(cmdstr + bytes_consumed, "%lf", &chemtbl[ind].MolarMass);

            rttbl->Dep_kinetic_all[i][ind] = -1.0;

            PIHMprintf(VL_VERBOSE, " Keq = %6.4f\n",
                rttbl->KeqKinect_all[i]);
            chemtbl[ind].Charge = 0;
            PIHMprintf(VL_VERBOSE,
                " MolarMass = %6.4f, MolarVolume = %6.4f\n\n",
                chemtbl[ind].MolarMass, chemtbl[ind].MolarVolume);

            break;
        }
    }
}

void ReadAdsorption(const char cmdstr[], int npoints, int keq_position,
    chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k;
    int             ind;
    int             ndep;
    double          dep;
    char            chemn[MAXSTRING];

    sscanf(cmdstr + bytes_consumed, "'%[^']' %d%n", chemn, &ndep, &bytes_now);
    bytes_consumed += bytes_now;
    wrap(chemn);

    for (i = 0; i < rttbl->NumSsc; i++)
    {
        ind = i + rttbl->NumStc;

        if (strcmp(chemtbl[ind].ChemName, chemn) == 0)
        {
            PIHMprintf(VL_VERBOSE,
                " Secondary surface complexation %s found in database!\n",
                chemtbl[ind].ChemName);
            PIHMprintf(VL_VERBOSE, " %s", cmdstr);
            chemtbl[ind].itype = ADSORPTION;
            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn,
                    &bytes_now);
                bytes_consumed += bytes_now;
                wrap(chemn);

                for (k = 0; k < rttbl->NumSdc; k++)
                {
                    if (strcmp(chemtbl[k].ChemName, chemn) == 0)
                    {
                        rttbl->Dependency[i][k] = dep;
                        break;
                    }
                }
            }

            for (j = 0; j < npoints; j++)
            {
                if (j + 1 == keq_position)
                {
                    sscanf(cmdstr + bytes_consumed, "%lf%n", &rttbl->Keq[i],
                        &bytes_now);
                    bytes_consumed += bytes_now;
                    PIHMprintf(VL_VERBOSE,
                        " Keq = %6.4f\n", rttbl->Keq[i]);
                }
                else
                {
                    sscanf(cmdstr + bytes_consumed, "%*f%n", &bytes_now);
                    bytes_consumed += bytes_now;
                }
            }
        }
    }
}

void ReadCationEchg(const char cmdstr[], double calval,
    chemtbl_struct chemtbl[], rttbl_struct *rttbl)
{
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i, j, k;
    int             ind;
    int             ndep;
    double          dep;
    char            chemn[MAXSTRING];

    if (sscanf(cmdstr + bytes_consumed, "'%[^']' %d%n",
        chemn, &ndep, &bytes_now) != 2)
    {
        return;
    }
    bytes_consumed += bytes_now;
    wrap(chemn);

    for (i = 0; i < rttbl->NumSsc; i++)
    {
        ind = i + rttbl->NumStc;

        if (strcmp(chemtbl[ind].ChemName, chemn) == 0)
        {
            PIHMprintf(VL_VERBOSE,
                " Secondary ion exchange %s found in database!\n",
                chemtbl[ind].ChemName);
            PIHMprintf(VL_VERBOSE, " %s", cmdstr);
            chemtbl[ind].itype = CATION_ECHG;
            for (j = 0; j < ndep; j++)
            {
                sscanf(cmdstr + bytes_consumed, "%lf '%[^']'%n", &dep, chemn,
                    &bytes_now);
                bytes_consumed += bytes_now;
                wrap(chemn);

                for (k = 0; k < rttbl->NumSdc; k++)
                {
                    if (strcmp(chemtbl[k].ChemName, chemn) == 0)
                    {
                        rttbl->Dependency[i][k] = dep;
                        break;
                    }
                }
            }
            sscanf(cmdstr + bytes_consumed, "%lf", &rttbl->Keq[i]);
            PIHMprintf(VL_VERBOSE, " Keq = %6.4f \n", rttbl->Keq[i]);

            rttbl->Keq[i] += calval;
            PIHMprintf(VL_VERBOSE, " After calibration: Keq = %6.4f \n",
                rttbl->Keq[i]);
        }
    }
}

void ReadMinKin(FILE *database, int NumStc, double calval, int *lno,
    char cmdstr[], chemtbl_struct chemtbl[], kintbl_struct *kintbl)
{
    int             i, j, k;
    int             ind;
    int             ndep;
    int             bytes_now;
    int             bytes_consumed = 0;
    double          dep;
    char            chemn[MAXSTRING];
    char            optstr[MAXSTRING];
    char            label[MAXSTRING];

    sscanf(cmdstr, "%s", chemn);
    wrap(chemn);

    if (strcmp(chemtbl[kintbl->position].ChemName, chemn) == 0)
    {
        NextLine(database, cmdstr, lno);
        sscanf(cmdstr, "%*s = %s", label);

        if (strcmp(kintbl->Label, label) == 0)
        {
            PIHMprintf(VL_VERBOSE,
                " \n Mineral kinetics %s %s found in database!\n",
                chemtbl[kintbl->position].ChemName, kintbl->Label);

            /* Read type */
            NextLine(database, cmdstr, lno);
            sscanf(cmdstr, "%s = %s", optstr, label);
            if (strcmp(optstr, "label") == 0)
            {
                if (strcmp(label, "tst") == 0)
                {
                    kintbl->type = TST;
                }
                else if (strcmp(label, "PrecipitationOnly") == 0)
                {
                    kintbl->type = PRCP_ONLY;
                }
                else if (strcmp(label, "DissolutionOnly") == 0)
                {
                    kintbl->type = DISS_ONLY;
                }
                else if (strcmp(label, "monod") == 0)
                {
                    kintbl->type = MONOD;
                }
            }

            /* Read rate */
            NextLine(database, cmdstr, lno);
            sscanf(cmdstr, "%s = %lf", optstr, &kintbl->rate);
            if (strcmp(optstr, "rate(25C)") == 0)
            {
                PIHMprintf(VL_VERBOSE, " Rate is %f\n", kintbl->rate);

                kintbl->rate += calval;
                PIHMprintf(VL_VERBOSE,
                    " After calibration: Rate is %f, cal->Rate = %f \n",
                    kintbl->rate, calval);
            }

            /* Read activation */
            NextLine(database, cmdstr, lno);
            sscanf(cmdstr, "%s = %lf", optstr, &kintbl->actv);
            if (strcmp(optstr, "activation") == 0)
            {
                PIHMprintf(VL_VERBOSE, " Activation is %f\n", kintbl->actv);
            }

            /* Read dependence
             * Assume that all mineral kinetics only depend on one species*/
            NextLine(database, cmdstr, lno);
            if (sscanf(cmdstr, "%s : %s %lf", optstr, chemn, &dep) == 3)
            {
                if (strcmp(optstr, "dependence") == 0)
                {
                    wrap(chemn);
                    kintbl->dep_position[0] =
                        FindChem(chemn, chemtbl, NumStc);
                    kintbl->num_dep = 1;
                    kintbl->dep_power[0] = dep;
                    PIHMprintf(VL_VERBOSE, " Dependency: %s %f\n",
                        chemn, kintbl->dep_power[0]);
                }
            }
            else
            {
                if (strcmp(optstr, "dependence") == 0)
                {
                    kintbl->dep_position[0] = -1;
                    kintbl->num_dep = 0;
                    kintbl->dep_power[0] = 0.0;
                    PIHMprintf(VL_VERBOSE, " No dependency.\n");
                }
            }

            NextLine(database, cmdstr, lno);
            while (cmdstr[0] != '+')
            {
                sscanf(cmdstr, "%s", optstr);

                if (strcmp(optstr, "biomass") == 0)
                {
                    /* Biomass term */
                    sscanf(cmdstr, "%*s : %s", chemn);
                    wrap(chemn);
                    PIHMprintf(VL_VERBOSE, " Biomass species: %s \n",
                        chemn);
                    kintbl->biomass_position =
                        FindChem(chemn, chemtbl, NumStc);
                    PIHMprintf(VL_VERBOSE,
                        " Biomass species position: %d \n",
                        kintbl->biomass_position);
                }
                else if (strcmp(optstr, "num_monod") == 0)
                {
                    /* Monod term */
                    sscanf(cmdstr, "%*s = %d", &kintbl->num_monod);
                    PIHMprintf(VL_VERBOSE,
                        " Number of monod term: %d\n", kintbl->num_monod);
                    if (kintbl->num_monod > 0)
                    {
                        int             mn;

                        NextLine(database, cmdstr, lno);
                        sscanf(cmdstr, "%*s : %n", &bytes_now);
                        bytes_consumed += bytes_now;

                        for (mn = 0; mn < kintbl->num_monod; mn++)
                        {
                            sscanf(cmdstr + bytes_consumed, "%s %lf",
                                chemn, &kintbl->monod_para[mn]);
                            wrap(chemn);
                            kintbl->monod_position[mn] =
                                FindChem(chemn, chemtbl, NumStc);
                            if (kintbl->monod_position[mn] < 0)
                            {
                                PIHMprintf(VL_ERROR,
                                    "Error finding monod_terms in species "
                                    "table.\n");
                                PIHMexit(EXIT_FAILURE);
                            }
                            PIHMprintf(VL_VERBOSE,
                                " Monod term: %s %f\n",
                                chemtbl[kintbl->monod_position[mn]].ChemName,
                                kintbl->monod_para[mn]);
                        }
                    }
                }
                else if (strcmp(optstr, "num_inhibition"))
                {
                    /* Inhibition term */
                    sscanf(cmdstr, "%*s = %d", &kintbl->num_inhib);
                    PIHMprintf(VL_VERBOSE,
                        " Number of inhibition term: %d\n", kintbl->num_inhib);

                    if (kintbl->num_inhib > 0)
                    {
                        int             in;

                        NextLine(database, cmdstr, lno);
                        sscanf(cmdstr, "%*s : %n", &bytes_now);
                        bytes_consumed += bytes_now;

                        for (in = 0; in < kintbl->num_inhib; in++)
                        {
                            sscanf(cmdstr + bytes_consumed, "%s %lf",
                                chemn, &kintbl->inhib_para[in]);
                            wrap(chemn);
                            kintbl->inhib_position[in] =
                                FindChem(chemn, chemtbl, NumStc);
                            if (kintbl->inhib_position[in] < 0)
                            {
                                PIHMprintf(VL_ERROR,
                                    "Error finding inhibition term in "
                                    "species table.\n");
                                PIHMexit(EXIT_FAILURE);
                            }
                            PIHMprintf(VL_VERBOSE,
                                " Inhibition term: %s %f\n",
                                chemtbl[kintbl->inhib_position[in]].ChemName,
                                kintbl->inhib_para[in]);
                        }
                    }
                }

                NextLine(database, cmdstr, lno);
            }
        }
    }
}
