/******************************************************************************
* This subroutine is used to calculate the reactions
* It uses a similar OS3D scheme as detailed inCrunchflow user's manual.
*
* If you have any questions, concerns, suggestions, please contact me at
* the following address:
*     Developer: Chen Bao <baochen.d.s@gmail.com>
*     Version  : pre-alpha
*     Date     : June, 2013
*****************************************************************************/
#include "pihm.h"

#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80
#define TOL        1E-7
#define SKIP_JACOB 1
#define sqr(a)  (a)*(a)

void wrap(char *str)
{

    char            word[WORD_WIDTH];
    sprintf(word, "'%s'", str);
    strcpy(str, word);
}

void Unwrap(char *str, const char *str0)
{
    int             i, j = 0;

    for (i = 0; i < (int)strlen(str0); i++)
    {
        if (str0[i] != '\'')
        {
            str[j] = str0[i];
            j++;
        }
    }

    str[j] = '\0';
}

void ReportError(const chemtbl_struct chemtbl[], const rttbl_struct *rttbl, vol_conc cell, Chem_Data CD)
{
    /*
     * This subroutine checks all the important parameters within a cell and output it to the error log
     */
    /* Diagnotic purposes */
    int             i;
    fprintf(stderr, "Error found at cell %d\n", cell.index);
    fprintf(stderr,
        "Volumetric properties:\n\t sat: %8.6f\t, height_o: %8.6f\t, height_t: %8.6f\t, area: %8.6f, Error: %d\n",
        cell.sat, cell.height_o, cell.height_t, cell.area, cell.illness);
    fprintf(stderr,
        "Chemical total concentrations: Chemical  Tot_conc  Log(Tot_conc)\n");
    for (i = 0; i < rttbl->NumStc; i++)
    {
        fprintf(stderr, "\t%20s\t%8.4f\t%8.4f\t\n", chemtbl[i].ChemName,
            cell.t_conc[i], log10(cell.t_conc[i]));
    }
    fprintf(stderr, "\n");
}

int SpeciationType(FILE *database, char *Name)
{
    /* This subroutine is used to find out what the input species is.
     * 0) not found within database
     * 1) aqueous
     * 2) adsorption
     * 3) cation exchange
     * 4) mineral
     */
    double          tmpval[WORDS_LINE];
    int             i, return_val;
    char            line[LINE_WIDTH], word[WORD_WIDTH];

    if (strcmp(Name, "pH") == 0)
        return (1);

    char          **tmpstr = (char **)malloc(WORDS_LINE * sizeof(char *));
    for (i = 0; i < WORDS_LINE; i++)
        tmpstr[i] = (char *)malloc(WORD_WIDTH * sizeof(char));

    return_val = 0;

    sprintf(word, "'%s'", Name);
    rewind(database);
    fgets(line, LINE_WIDTH, database);
    while (keymatch(line, "'End of primary'", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            if (strcmp(word, tmpstr[0]) == 0)
            {
                for (i = 0; i < WORDS_LINE; i++)
                    free(tmpstr[i]);
                free(tmpstr);
                return (AQUEOUS);
            }
        }
        fgets(line, LINE_WIDTH, database);
    }
    fgets(line, LINE_WIDTH, database);
    while (keymatch(line, "'End of secondary'", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            if (strcmp(word, tmpstr[0]) == 0)
            {
                for (i = 0; i < WORDS_LINE; i++)
                    free(tmpstr[i]);
                free(tmpstr);
                return (5);
            }
        }
        fgets(line, LINE_WIDTH, database);
    }
    fgets(line, LINE_WIDTH, database);
    while (keymatch(line, "'End of minerals'", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            if (strcmp(word, tmpstr[0]) == 0)
            {
                for (i = 0; i < WORDS_LINE; i++)
                    free(tmpstr[i]);
                free(tmpstr);
                return (MINERAL);
            }
        }
        fgets(line, LINE_WIDTH, database);
    }
    fgets(line, LINE_WIDTH, database);
    while (strcmp(line, "End of surface complexation\r\n") != 1)
    {
        /* Notice that in crunchflow database, starting from surface
         * complexation, there is not apostrophe marks around blocking keywords
         */
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            if (strcmp(word, tmpstr[0]) == 0)
            {
                for (i = 0; i < WORDS_LINE; i++)
                    free(tmpstr[i]);
                free(tmpstr);
                return (ADSORPTION);
            }
        }
        fgets(line, LINE_WIDTH, database);
    }
    fgets(line, LINE_WIDTH, database);
    while (!feof(database))
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            if (strcmp(word, tmpstr[0]) == 0)
            {
                for (i = 0; i < WORDS_LINE; i++)
                    free(tmpstr[i]);
                free(tmpstr);
                return (CATION_ECHG);
            }
        }
        fgets(line, LINE_WIDTH, database);
    }

    for (i = 0; i < WORDS_LINE; i++)
        free(tmpstr[i]);
    free(tmpstr);

    return (0);
}

void Lookup(FILE *database, chemtbl_struct chemtbl[], kintbl_struct kintbl[],
    rttbl_struct *rttbl, Chem_Data CD)
{
    double          tmpval[WORDS_LINE];
    /* Kinetic reactions is currently only applicable to minerals */
    int             i, j, k, l, keq_position = 0, total_temp_points;
    int             mn, in;     // 08.19 Wei
    char            line[LINE_WIDTH], tmp[WORD_WIDTH];
    char          **tmpstr = (char **)malloc(WORDS_LINE * sizeof(char *));

    for (i = 0; i < WORDS_LINE; i++)
    {
        tmpstr[i] = (char *)malloc(WORD_WIDTH * sizeof(char));
        memset(tmpstr[i], 0, WORD_WIDTH);
    }

    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
        wrap(chemtbl[i].ChemName);
    rewind(database);
    fgets(line, LINE_WIDTH, database);
    while (keymatch(line, "'temperature points'", tmpval, tmpstr) != 1)
    {
        fgets(line, LINE_WIDTH, database);
    }
    total_temp_points = tmpval[0];
    for (i = 0; i < tmpval[0]; i++)
    {
        if (tmpval[i + 1] == rttbl->Temperature)
        {
            fprintf(stderr,
                "\n Temperature point %6.4f C found in database!\n\n",
                tmpval[i + 1]);
            keq_position = i + 1;
        }
    }
    while (keymatch(line, "'Debye-Huckel adh'", tmpval, tmpstr) != 1)
    {
        fgets(line, LINE_WIDTH, database);
    }
    CD->DH.adh = tmpval[keq_position - 1];
    while (keymatch(line, "'Debye-Huckel bdh'", tmpval, tmpstr) != 1)
    {
        fgets(line, LINE_WIDTH, database);
    }
    CD->DH.bdh = tmpval[keq_position - 1];
    while (keymatch(line, "'Debye-Huckel bdt'", tmpval, tmpstr) != 1)
    {
        fgets(line, LINE_WIDTH, database);
    }
    CD->DH.bdt = tmpval[keq_position - 1];

    fprintf(stderr,
        " Debye-Huckel Parameters set to A=%6.4f; B=%6.4f; b=%6.4f\n\n",
        CD->DH.adh, CD->DH.bdh, CD->DH.bdt);

    rewind(database);
    fgets(line, LINE_WIDTH, database);
    while (keymatch(line, "'End of primary'", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            for (i = 0; i < rttbl->NumStc; i++)
                if (strcmp(chemtbl[i].ChemName, tmpstr[0]) == 0)
                {
                    fprintf(stderr,
                        " Primary species %s found in database!\n MolarMass = %6.4f\n\n",
                        chemtbl[i].ChemName, tmpval[2]);
                    chemtbl[i].MolarMass = tmpval[2];
                    chemtbl[i].Charge = tmpval[1];
                    chemtbl[i].SizeF = tmpval[0];
                }
        }
        fgets(line, LINE_WIDTH, database);
    }
    while (keymatch(line, "'End of secondary'", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            for (i = rttbl->NumStc; i < rttbl->NumSsc + rttbl->NumStc; i++)
                if (strcmp(chemtbl[i].ChemName, tmpstr[0]) == 0)
                {
                    fprintf(stderr,
                        " Secondary species %s found in database!\n",
                        chemtbl[i].ChemName);
                    fprintf(stderr, " %s", line);
                    chemtbl[i].itype = AQUEOUS;
                    for (j = 0; j < WORDS_LINE; j++)
                    {
                        for (k = 0; k < rttbl->NumSdc; k++)
                            if (strcmp(chemtbl[k].ChemName,
                                    tmpstr[j]) == 0)
                                CD->Dependency[i - rttbl->NumStc][k] =
                                    atof(tmpstr[j - 1]);
                    }
                    CD->Keq[i - rttbl->NumStc] =
                        tmpval[(int)tmpval[0] + keq_position];
                    fprintf(stderr, " Keq = %6.4f\n", CD->Keq[i - rttbl->NumStc]);
                    chemtbl[i].MolarMass =
                        tmpval[(int)tmpval[0] + total_temp_points + 3];
                    chemtbl[i].Charge =
                        tmpval[(int)tmpval[0] + total_temp_points + 2];
                    chemtbl[i].SizeF =
                        tmpval[(int)tmpval[0] + total_temp_points + 1];
                    fprintf(stderr,
                        " MolarMass = %6.4f, Charge = %6.4f, SizeFactor = %6.4f\n\n",
                        chemtbl[i].MolarMass, chemtbl[i].Charge,
                        chemtbl[i].SizeF);
                }
        }
        fgets(line, LINE_WIDTH, database);
    }
    while (keymatch(line, "'End of minerals'", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            for (i = NumSpc + rttbl->NumAds + rttbl->NumCex; i < rttbl->NumStc; i++)
                if (strcmp(chemtbl[i].ChemName, tmpstr[0]) == 0)
                {
                    fprintf(stderr, " Mineral %s found in database!\n",
                        chemtbl[i].ChemName);
                    fprintf(stderr, " %s", line);
                    chemtbl[i].itype = MINERAL;
                    CD->KeqKinect_all[i - NumSpc - rttbl->NumAds -
                        rttbl->NumCex] = tmpval[(int)tmpval[1] + keq_position + 1];
                    for (j = 1; j < WORDS_LINE; j++)
                    {
                        for (k = 0; k < rttbl->NumStc + rttbl->NumSsc; k++)
                            if (strcmp(chemtbl[k].ChemName,
                                    tmpstr[j]) == 0)
                            {
                                if (k < rttbl->NumStc)
                                {
                                    CD->Dep_kinetic_all[i - rttbl->NumStc +
                                        rttbl->NumMin][k] = atof(tmpstr[j - 1]);
                                }
                                else
                                {
                                    for (l = 0; l < NumSpc; l++)    /* Number of primary species in the rt simulator      */
                                        CD->Dep_kinetic_all[i - rttbl->NumStc +
                                            rttbl->NumMin][l] +=
                                            atof(tmpstr[j -
                                                1]) * CD->Dependency[k -
                                            rttbl->NumStc][l];
                                    CD->KeqKinect_all[i - NumSpc -
                                        rttbl->NumAds - rttbl->NumCex] +=
                                        atof(tmpstr[j - 1]) * CD->Keq[k -
                                        rttbl->NumStc];
                                }
                            }
                    }
                    CD->Dep_kinetic_all[i - rttbl->NumStc + rttbl->NumMin][i] = -1.0;
                    fprintf(stderr, " Keq = %6.4f\n",
                        CD->KeqKinect_all[i - NumSpc - rttbl->NumAds -
                            rttbl->NumCex]);
                    chemtbl[i].MolarMass =
                        tmpval[(int)tmpval[1] + total_temp_points + 2];
                    chemtbl[i].MolarVolume = tmpval[0];
                    chemtbl[i].Charge = 0;
                    fprintf(stderr,
                        " MolarMass = %6.4f, MolarVolume = %6.4f\n\n",
                        chemtbl[i].MolarMass, chemtbl[i].MolarVolume);
                }
        }
        fgets(line, LINE_WIDTH, database);
    }

    for (i = 0; i < rttbl->NumMkr + rttbl->NumAkr; i++)
    {
        for (j = NumSpc + rttbl->NumAds + rttbl->NumCex; j < rttbl->NumStc; j++)
        {
            strcpy(tmp, kintbl[i].species);
            wrap(tmp);
            if (strcmp(tmp, chemtbl[j].ChemName) == 0)
            {
                fprintf(stderr,
                    " Selecting the kinetic species %s from all possible species.\n\n",
                    tmp);
                CD->KeqKinect[i] =
                    CD->KeqKinect_all[j - rttbl->NumStc + rttbl->NumMin];
                for (k = 0; k < rttbl->NumStc; k++)
                {
                    CD->Dep_kinetic[i][k] =
                        CD->Dep_kinetic_all[j - rttbl->NumStc + rttbl->NumMin][k];
                }
            }
        }
    }
    while (strcmp(line, "End of surface complexation\r\n") != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            for (i = rttbl->NumStc; i < rttbl->NumSsc + rttbl->NumStc; i++)
                if (strcmp(chemtbl[i].ChemName, tmpstr[0]) == 0)
                {
                    fprintf(stderr,
                        " Secondary surface complexation %s found in database!\n",
                        chemtbl[i].ChemName);
                    fprintf(stderr, " %s", line);
                    chemtbl[i].itype = ADSORPTION;
                    for (j = 0; j < WORDS_LINE; j++)
                    {
                        for (k = 0; k < rttbl->NumSdc; k++)
                            if (strcmp(chemtbl[k].ChemName,
                                    tmpstr[j]) == 0)
                                CD->Dependency[i - rttbl->NumStc][k] =
                                    atof(tmpstr[j - 1]);
                    }
                    CD->Keq[i - rttbl->NumStc] =
                        tmpval[(int)tmpval[0] + keq_position];
                    fprintf(stderr, " Keq = %6.4f\n", CD->Keq[i - rttbl->NumStc]);
                }
        }
        fgets(line, LINE_WIDTH, database);
    }
    while (!feof(database))
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            for (i = rttbl->NumStc; i < rttbl->NumSsc + rttbl->NumStc; i++)
                if (strcmp(chemtbl[i].ChemName, tmpstr[0]) == 0)
                {
                    fprintf(stderr,
                        " Secondary ion exchange %s found in database!\n",
                        chemtbl[i].ChemName);
                    fprintf(stderr, " %s", line);
                    chemtbl[i].itype = CATION_ECHG;
                    for (j = 0; j < WORDS_LINE; j++)
                    {
                        for (k = 0; k < rttbl->NumSdc; k++)
                            if (strcmp(chemtbl[k].ChemName,
                                    tmpstr[j]) == 0)
                                CD->Dependency[i - rttbl->NumStc][k] =
                                    atof(tmpstr[j - 1]);
                    }
                    CD->Keq[i - rttbl->NumStc] = tmpval[(int)tmpval[0] + 1];
                    fprintf(stderr, " Keq = %6.4f \n", CD->Keq[i - rttbl->NumStc]);

                    CD->Keq[i - rttbl->NumStc] =
                        tmpval[(int)tmpval[0] + 1] + CD->CalXsorption;
                    fprintf(stderr, " After calibration: Keq = %6.4f \n",
                        CD->Keq[i - rttbl->NumStc]);
                }
        }
        fgets(line, LINE_WIDTH, database);
    }

    for (i = 0; i < rttbl->NumMkr; i++)
    {
        rewind(database);
        fgets(line, LINE_WIDTH, database);
        while (strcmp(line, "Begin mineral kinetics\r\n") != 0)
        {
            fgets(line, LINE_WIDTH, database);
        }
        fgets(line, LINE_WIDTH, database);
        while (strcmp(line, "End of mineral kinetics\r\n") != 0)
        {
            if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
            {
                if (strcmp(kintbl[i].species, tmpstr[0]) == 0)
                {
                    fgets(line, LINE_WIDTH, database);
                    keymatch(line, "NULL", tmpval, tmpstr);
                    if (strcmp(kintbl[i].Label, tmpstr[2]) == 0)
                    {
                        fprintf(stderr,
                            " \n Mineral kinetics %s %s found in database!\n",
                            kintbl[i].species, kintbl[i].Label);
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[2], "tst") == 0)
                            kintbl[i].type = 1;
                        if (strcmp(tmpstr[2], "PrecipitationOnly") == 0)
                            kintbl[i].type = 2;
                        if (strcmp(tmpstr[2], "DissolutionOnly") == 0)
                            kintbl[i].type = 3;
                        if (strcmp(tmpstr[2], "monod") == 0)
                            kintbl[i].type = 4;
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "rate(25C)") == 0)
                        {
                            kintbl[i].rate = tmpval[0];
                            fprintf(stderr, " Rate is %f\n",
                                kintbl[i].rate);

                            kintbl[i].rate = tmpval[0] + CD->CalRate;
                            fprintf(stderr,
                                " After calibration: Rate is %f, CD->CalRate = %f \n",
                                kintbl[i].rate, CD->CalRate);
                        }
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "activation") == 0)
                        {
                            kintbl[i].actv = tmpval[0];
                            fprintf(stderr, " Activation is %f\n",
                                kintbl[i].actv);
                        }
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "dependence") == 0)
                        {
                            strcpy(kintbl[i].dep_species[0], tmpstr[2]);
                            wrap(kintbl[i].dep_species[0]);
                            kintbl[i].num_dep = 1;
                            /* Assume that all mineral kinetic only depend on
                             * one species !! */

                            /* Require further elaboration after this !! */
                            kintbl[i].dep_position[0] = -999;
                            for (k = 0; k < rttbl->NumStc; k++)
                            {
                                if (strcmp(kintbl[i].dep_species[0],
                                        chemtbl[k].ChemName) == 0)
                                    kintbl[i].dep_position[0] = k;
                            }
                            if (kintbl[i].dep_position[0] == -999)
                            {
                                kintbl[i].num_dep = 0;
                            }
                            kintbl[i].dep_power[0] = tmpval[0];
                            fprintf(stderr, " Dependency: %s %f\n",
                                kintbl[i].dep_species[0],
                                kintbl[i].dep_power[0]);
                        }

                        /* Biomass term */
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "biomass") == 0)
                        {
                            strcpy(kintbl[i].biomass_species, tmpstr[2]);
                            wrap(kintbl[i].biomass_species);
                            fprintf(stderr, " Biomass species: %s \n",
                                kintbl[i].biomass_species);
                            kintbl[i].biomass_position = -999;
                            for (k = 0; k < rttbl->NumStc; k++)
                            {
                                if (strcmp(kintbl[i].biomass_species,
                                        chemtbl[k].ChemName) == 0)
                                {
                                    kintbl[i].biomass_position = k;
                                    fprintf(stderr,
                                        " Biomass species position: %d \n",
                                        kintbl[i].biomass_position);
                                }
                            }
                        }

                        /* Monod term */
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "num_monod") == 0)
                        {
                            kintbl[i].num_monod = tmpval[0];
                            fprintf(stderr, " Number of monod term: %d\n",
                                kintbl[i].num_monod);
                        }
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "monod_terms") == 0)
                        {
                            for (mn = 0; mn < kintbl[i].num_monod; mn++)
                            {
                                /* Note tmpstr indexing not same with tmpval
                                 * indexing */
                                strcpy(kintbl[i].monod_species[mn], tmpstr[mn * 2 + 2]);
                                wrap(kintbl[i].monod_species[mn]);
                                kintbl[i].monod_position[mn] = -999;
                                for (k = 0; k < rttbl->NumStc; k++)
                                {
                                    if (strcmp(
                                        kintbl[i].monod_species[mn],
                                        chemtbl[k].ChemName) == 0)
                                        kintbl[i].monod_position[mn] = k;
                                }
                                kintbl[i].monod_para[mn] = tmpval[mn + 0];
                                fprintf(stderr, " Monod term: %s %f\n",
                                    kintbl[i].monod_species[mn],
                                    kintbl[i].monod_para[mn]);
                            }
                        }

                        /* Inhibition term */
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "num_inhibition") == 0)
                        {
                            kintbl[i].num_inhib = tmpval[0];
                            fprintf(stderr, " Number of inhibition term: %d\n",
                                kintbl[i].num_inhib);
                        }
                        fgets(line, LINE_WIDTH, database);
                        keymatch(line, "NULL", tmpval, tmpstr);
                        if (strcmp(tmpstr[0], "inhibition") == 0)
                        {
                            for (in = 0; in < kintbl[i].num_inhib; in++)
                            {
                                /* Note tmpstr indexing not same with tmpval */
                                strcpy(kintbl[i].inhib_species[in],
                                    tmpstr[in * 2 + 2]);
                                wrap(kintbl[i].inhib_species[in]);
                                kintbl[i].inhib_position[in] = -999;
                                for (k = 0; k < rttbl->NumStc; k++)
                                {
                                    if (strcmp(kintbl[i].inhib_species[in],
                                            chemtbl[k].ChemName) == 0)
                                        kintbl[i].inhib_position[in] = k;
                                }
                                kintbl[i].inhib_para[in] = tmpval[in + 0];
                                fprintf(stderr, " Inhibition term: %s %f\n",
                                    kintbl[i].inhib_species[in],
                                    kintbl[i].inhib_para[in]);
                            }
                        }
                        wrap(kintbl[i].species);
                    }
                }
            }
            fgets(line, LINE_WIDTH, database);
        }
    }

    for (i = 0; i < rttbl->NumMkr; i++)
        for (j = 0; j < rttbl->NumStc; j++)
        {
            if (strcmp(kintbl[i].species, chemtbl[j].ChemName) == 0)
            {
                kintbl[i].position = j;
                fprintf(stderr,
                    " \n Position_check (NumMkr[i] vs NumStc[j]) %s (%d), %s (%d)\n",
                    kintbl[i].species, i, chemtbl[j].ChemName, j);
            }
        }

    for (i = 0; i < rttbl->NumStc; i++)
        CD->Totalconc[i][i] = 1.0;

    for (i = rttbl->NumStc; i < rttbl->NumStc + rttbl->NumSsc; i++)
        for (j = 0; j < rttbl->NumSdc; j++)
            CD->Totalconc[j][i] += CD->Dependency[i - rttbl->NumStc][j];
    if (rttbl->NumSsc > 0)
    {
        fprintf(stderr, " \n Dependency Matrix!\n");

        fprintf(stderr, "%-15s", " ");
        for (i = 0; i < rttbl->NumSdc; i++)
            fprintf(stderr, "%-14s", chemtbl[i].ChemName);
        fprintf(stderr, "\n");

        for (i = 0; i < rttbl->NumSsc; i++)    /* Number of secondary speices in the simulator */
        {
            fprintf(stderr, " %-14s", chemtbl[i + rttbl->NumStc].ChemName);
            for (j = 0; j < rttbl->NumSdc; j++)    /* Number of independent species (others depending on these species) */
                fprintf(stderr, "%-14.2f", CD->Dependency[i][j]);
            fprintf(stderr, " %6.2f\n", CD->Keq[i]);
        }
    }

    fprintf(stderr, " \n Total Concentration Matrix!\n");
    fprintf(stderr, "%-18s", " ");
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
        fprintf(stderr, "%-14s", chemtbl[i].ChemName);
    fprintf(stderr, "\n");
    for (i = 0; i < rttbl->NumStc; i++)
    {
        fprintf(stderr, " Sum%-14s", chemtbl[i].ChemName);
        for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
            fprintf(stderr, "%-14.2f", CD->Totalconc[i][j]);
        fprintf(stderr, "\n");
    }

    fprintf(stderr, " \n Kinetic Mass Matrx!\n");
    fprintf(stderr, "%-15s", " ");
    for (i = 0; i < rttbl->NumStc; i++)
        fprintf(stderr, "%-14s", chemtbl[i].ChemName);
    fprintf(stderr, "\n");
    for (j = 0; j < rttbl->NumMkr + rttbl->NumAkr; j++)
    {
        fprintf(stderr, " %-14s",
            chemtbl[j + NumSpc + rttbl->NumAds + rttbl->NumCex].ChemName);
        for (i = 0; i < rttbl->NumStc; i++)
        {
            fprintf(stderr, "%-14f", CD->Dep_kinetic[j][i]);
        }
        fprintf(stderr, " Keq = %-6.2f\n", CD->KeqKinect[j]);
    }
    for (i = 0; i < WORDS_LINE; i++)
        free(tmpstr[i]);
    free(tmpstr);
}

int Speciation(chemtbl_struct chemtbl[], rttbl_struct *rttbl, Chem_Data CD, int cell)
{
    /* if speciation flg = 1, pH is defined
     * if speciation flg = 0, all defined value is total concentration */
    int             i, j, k, speciation_flg = CD->SPCFlg, num_spe =
        rttbl->NumStc + rttbl->NumSsc;
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
            Keq[i] = CD->Keq[i];
    }

    adh = CD->DH.adh;
    bdh = CD->DH.bdh;
    bdt = CD->DH.bdt;

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
            tmpval += tmpconc[j] * CD->Dependency[i][j];
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
                        sqr(chemtbl[i].Charge);
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
                            (-adh * sqr(chemtbl[i].Charge) * Iroot) / (1 +
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
                    tmpval += (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
                }
                tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                tmpconc[i + rttbl->NumStc] = tmpval;
            }
            for (i = 0; i < rttbl->NumStc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < num_spe; j++)
                {
                    tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
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
                                (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
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
                                    CD->Totalconc[i][j] * pow(10, tmpconc[j]);
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
                        tmpconc[i]) * sqr(chemtbl[i].Charge);
                }
                Iroot = sqrt(I);
                for (i = 0; i < num_spe; i++)
                {
                    if (chemtbl[i].itype == MINERAL)
                        gamma[i] = -tmpconc[i];
                    else
                        gamma[i] =
                            (-adh * sqr(chemtbl[i].Charge) * Iroot) / (1 +
                            bdh * chemtbl[i].SizeF * Iroot) + bdt * I;
                }
            }
            /* gamma stores log10gamma[i]. */
            for (i = 0; i < rttbl->NumSsc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < rttbl->NumSdc; j++)
                {
                    tmpval += (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
                }
                tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                tmpconc[i + rttbl->NumStc] = tmpval;
            }
            for (i = 0; i < rttbl->NumStc; i++)
            {
                tmpval = 0.0;
                for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
                {
                    tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
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
                            (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
                    tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                    tmpconc[i + rttbl->NumStc] = tmpval;
                }
                for (i = 0; i < rttbl->NumStc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
                    {
                        tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
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
            tmpval += (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
        }
        tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
        tmpconc[i + rttbl->NumStc] = tmpval;
    }
    for (i = 0; i < rttbl->NumStc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
        {
            tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
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


void React(double stepsize, const chemtbl_struct chemtbl[],
    const kintbl_struct kintbl[], const rttbl_struct *rttbl, ctrl_struct *ctrl,
    Chem_Data CD, vol_conc *Vcele)
{
    int             k, j;
    double          substep;

    if (Vcele->illness < 20)
    {
        if (_React(stepsize, chemtbl, kintbl, rttbl, CD, Vcele))
        {
            fprintf(stderr, "  ---> React failed at cell %12d.\t",
                    Vcele->index);

            substep = 0.5 * stepsize;
            k = 2;

            while ((j = _React(substep, chemtbl, kintbl, rttbl, CD, Vcele)))
            {
                substep = 0.5 * substep;
                k = 2 * k;
                if (substep < 0.5)
                    break;
            }

            if (j == 0)
            {
                fprintf(stderr,
                        " Reaction passed with step equals to %f (1/%d)\n",
                        substep, k);
                for (j = 1; j < k; j++)
                {
                    _React(substep, chemtbl, kintbl, rttbl, CD, Vcele);
                }
            }
        }
    }
}

int _React(double stepsize, const chemtbl_struct chemtbl[],
    const kintbl_struct kintbl[], const rttbl_struct *rttbl, Chem_Data CD, vol_conc *Vcele)
{
    if (Vcele->sat < 1.0E-2)
    {
        /* very dry, no reaction can take place */
        return (0);
    }
    int             i, j, k, control, num_spe =
        rttbl->NumStc + rttbl->NumSsc, min_pos, pivot_flg;
    int             mn, in;
    double          monodterm = 1.0, inhibterm = 1.0;
    int             stc = rttbl->NumStc, ssc = rttbl->NumSsc, nkr =
        rttbl->NumMkr + rttbl->NumAkr, smc = rttbl->NumMin;
    double         *residue, *residue_t, *tmpconc, *totconc, *area, *error,
        *gamma, *Keq, *Rate_pre, *IAP, *dependency, *Rate_spe, *Rate_spe_t,
        *Rate_spet;
    residue = (double *)malloc(stc * sizeof(double));
    residue_t = (double *)malloc(stc * sizeof(double));
    tmpconc = (double *)malloc(num_spe * sizeof(double));
    totconc = (double *)malloc(stc * sizeof(double));
    area = (double *)malloc(smc * sizeof(double));
    error = (double *)malloc(stc * sizeof(double));
    gamma = (double *)malloc(num_spe * sizeof(double));
    Keq = (double *)malloc(ssc * sizeof(double));
    Rate_pre = (double *)malloc(nkr * sizeof(double));
    IAP = (double *)malloc(nkr * sizeof(double));
    dependency = (double *)malloc(nkr * sizeof(double));
    Rate_spe = (double *)malloc(stc * sizeof(double));
    Rate_spe_t = (double *)malloc(stc * sizeof(double));
    Rate_spet = (double *)malloc(stc * sizeof(double));
    long int       *p =
        (long int *)malloc((rttbl->NumStc - rttbl->NumMin) * sizeof(long int));
    realtype       *x_ =
        (realtype *)malloc((rttbl->NumStc - rttbl->NumMin) * sizeof(realtype));
    double          tmpval, tmpprb, inv_sat, I, Iroot, tmpKeq, adh,
        bdh, bdt, maxerror = 1, surf_ratio, tot_cec, tmpprb_inv;
    realtype      **jcb;

    /* Build model data structure from pointer address */
    control = 0;
    tmpprb = 1.0E-2;
    tmpprb_inv = 1.0 / tmpprb;
    inv_sat = 1.0 / Vcele->sat;

    for (i = 0; i < rttbl->NumMin; i++)
    {
        area[i] = CD->CalSSA *
            Vcele->p_para[i + rttbl->NumStc - rttbl->NumMin] *
            Vcele->p_conc[i + rttbl->NumStc - rttbl->NumMin] *
            chemtbl[i + rttbl->NumStc - rttbl->NumMin].MolarMass;
    }

    if (CD->SUFEFF)
    {
        if (Vcele->sat < 1.0)
        {
            //surf_ratio = 1.0;  /* # 1 function */
            surf_ratio = exp(Vcele->sat) - 1.0;    /* # 3 function */
            //surf_ratio = 1.0 - pow(exp(-1.0/Vcele->sat), 0.6); /* # 4 function */
            for (i = 0; i < rttbl->NumMin; i++)
            {
                area[i] *= surf_ratio;
            }
        }
    }   /* Lichtner's 2 third law if SUF_EFF is turned on. */

    for (j = 0; j < rttbl->NumStc; j++)
    {
        Rate_spe[j] = 0.0;
    }

    for (i = 0; i < rttbl->NumMkr + rttbl->NumAkr; i++)
    {
        min_pos = kintbl[i].position - rttbl->NumStc + rttbl->NumMin;

        if (kintbl[i].type == 1)  /* TST rate */
        {
            IAP[i] = 0.0;
            for (j = 0; j < rttbl->NumStc; j++)
            {
                IAP[i] += log10(Vcele->p_actv[j]) *
                    CD->Dep_kinetic[min_pos][j];
            }
            IAP[i] = pow(10, IAP[i]);
            tmpKeq = pow(10, CD->KeqKinect[min_pos]);
            dependency[i] = 1.0;
            for (k = 0; k < kintbl[i].num_dep; k++)
                dependency[i] *=
                    pow(Vcele->p_actv[kintbl[i].dep_position[k]],
                    kintbl[i].dep_power[k]);
            /* Calculate the predicted rate depending on the type of rate law!  */
            Rate_pre[i] = area[min_pos] * (pow(10, kintbl[i].rate)) *
                dependency[i] * (1 - (IAP[i] / tmpKeq));
            /* Rate_pre: rate per reaction (mol (L water)-1 s-1)
             * area: m2/L water
             * rate: mol/m2/s
             * dependency: dimensionless */
        }
        else if (kintbl[i].type == 4) /* Monod rate */
        {
            monodterm = 1.0;    /* re-set for new species */
            inhibterm = 1.0;    /*re-set for new species */

            /* Calculate rate */
            for (mn = 0; mn < kintbl[i].num_monod; mn++)
            {
                monodterm *=
                    Vcele->p_conc[kintbl[i].monod_position[mn]] /
                    (Vcele->p_conc[kintbl[i].monod_position[mn]] +
                    kintbl[i].monod_para[mn]);
            }

            for (in = 0; in < kintbl[i].num_inhib; in++)
            {
                inhibterm *= kintbl[i].inhib_para[in] /
                    (kintbl[i].inhib_para[in] +
                    Vcele->p_conc[kintbl[i].inhib_position[in]]);
            }

            /* Based on CrunchTope */
            Rate_pre[i] = area[min_pos] * pow(10, kintbl[i].rate) * monodterm;
        }

        for (j = 0; j < rttbl->NumStc; j++)
        {
            Rate_spe[j] += Rate_pre[i] * CD->Dep_kinetic[min_pos][j];
        }
    }

    for (i = 0; i < rttbl->NumMkr + rttbl->NumAkr; i++)
    {
        min_pos = kintbl[i].position - rttbl->NumStc + rttbl->NumMin;
        if (Rate_pre[i] < 0.0)
        {
            if (Vcele->p_conc[min_pos + rttbl->NumStc - rttbl->NumMin] < 1.0E-8) /* mineral cutoff when mineral is disappearing */
                area[min_pos] = 0.0;
        }
    }

    for (i = 0; i < NumSpc; i++)
        if (chemtbl[i].itype == AQUEOUS) /* 01.21 aqueous species, saturation term for aqueous volume */
            Rate_spe[i] *= inv_sat;

    jcb = newDenseMat(rttbl->NumStc - rttbl->NumMin, rttbl->NumStc - rttbl->NumMin);

    if (rttbl->TEMcpl == 0)
    {
        for (i = 0; i < rttbl->NumSsc; i++)
            Keq[i] = CD->Keq[i];
    }

    adh = CD->DH.adh;
    bdh = CD->DH.bdh;
    bdt = CD->DH.bdt;

    for (i = 0; i < rttbl->NumStc; i++)
    {
        tmpconc[i] = log10(Vcele->p_conc[i]);
    }
    for (i = 0; i < rttbl->NumSsc; i++)
    {
        tmpconc[i + rttbl->NumStc] = log10(Vcele->s_conc[i]);
    }
    tot_cec = 0.0;
    for (i = 0; i < num_spe; i++)
    {
        if (chemtbl[i].itype == CATION_ECHG)
        {
            tot_cec += pow(10, tmpconc[i]);
        }
    }

    I = 0;
    for (i = 0; i < num_spe; i++)
    {
        I += 0.5 * pow(10, tmpconc[i]) * sqr(chemtbl[i].Charge);
    }
    Iroot = sqrt(I);
    for (i = 0; i < num_spe; i++)
    {
        switch (chemtbl[i].itype)
        {
            case AQUEOUS:
                gamma[i] =
                    (-adh * sqr(chemtbl[i].Charge) * Iroot) / (1 +
                    bdh * chemtbl[i].SizeF * Iroot) + bdt * I;
                break;
            case ADSORPTION:
                gamma[i] = log10(Vcele->sat);
                break;
            case CATION_ECHG:
                gamma[i] = -log10(tot_cec);
                break;
            case MINERAL:
                gamma[i] = -tmpconc[i];
                break;
        }
    }

    while (maxerror > TOL)
    {
        for (i = 0; i < rttbl->NumSsc; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->NumSdc; j++)
            {
                tmpval += (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
            }
            tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
            tmpconc[i + rttbl->NumStc] = tmpval;
        }

        for (j = 0; j < rttbl->NumStc; j++)
        {
            Rate_spet[j] = 0.0;
        }

        for (i = 0; i < rttbl->NumMkr + rttbl->NumAkr; i++)
        {
            min_pos = kintbl[i].position - rttbl->NumStc + rttbl->NumMin;

            if (kintbl[i].type == 1)  /* TST rate */
            {
                IAP[i] = 0.0;
                for (j = 0; j < rttbl->NumStc; j++)
                {
                    if (chemtbl[j].itype != MINERAL)
                    {
                        IAP[i] += (tmpconc[j] + gamma[j]) *
                            CD->Dep_kinetic[min_pos][j];
                    }
                }
                IAP[i] = pow(10, IAP[i]);
                tmpKeq = pow(10, CD->KeqKinect[min_pos]);
                /*
                 * if ( IAP[i] < tmpKeq)
                 * rct_drct[i] = 1.0;
                 * if ( IAP[i] > tmpKeq)
                 * rct_drct[i] = -1.0;
                 * if ( IAP[i] == tmpKeq)
                 * rct_drct[i] = 0.0;
                 */
                dependency[i] = 0.0;
                for (k = 0; k < kintbl[i].num_dep; k++)
                    dependency[i] +=
                        (tmpconc[kintbl[i].dep_position[k]] +
                        gamma[kintbl[i].dep_position[k]]) *
                        kintbl[i].dep_power[k];
                dependency[i] = pow(10, dependency[i]);
                /* Calculate the predicted rate depending on the type of rate law!  */
                Rate_pre[i] = area[min_pos] * (pow(10, kintbl[i].rate)) *
                    dependency[i] * (1 - (IAP[i] / tmpKeq));
                /* Rate_pre: in mol / L water / s
                 * area: m2/L water
                 * rate: mol/m2/s
                 * dependency: dimensionless;
                 */
            }
            else if (kintbl[i].type == 4)
            {
                monodterm = 1.0;
                inhibterm = 1.0;

                /* Calculate rate */
                for (mn = 0; mn < kintbl[i].num_monod; mn++)
                {
                    monodterm *=
                        Vcele->p_conc[kintbl[i].monod_position[mn]] /
                        (Vcele->p_conc[kintbl[i].monod_position[mn]] +
                        kintbl[i].monod_para[mn]);
                }

                for (in = 0; in < kintbl[i].num_inhib; in++)
                {
                    inhibterm *=
                        kintbl[i].inhib_para[in] /
                        (kintbl[i].inhib_para[in] +
                        Vcele->p_conc[kintbl[i].inhib_position[in]]);
                }

                /* Based on CrunchTope */
                Rate_pre[i] =
                    area[min_pos] * pow(10, kintbl[i].rate) * monodterm;
            }

            for (j = 0; j < rttbl->NumStc; j++)
            {
                Rate_spet[j] += Rate_pre[i] * CD->Dep_kinetic[min_pos][j];
            }
            /* Adjust the unit of the calcuated rate. Note that for mineral, the
             * unit of rate and the unit of concentration are mol/L porous media
             * For the aqueous species, the unit of the rate and the unit of the
             * concentration are mol/L pm and mol/L water respectively. */
        }

        for (i = 0; i < NumSpc; i++)
            if (chemtbl[i].itype == AQUEOUS)
                Rate_spet[i] *= inv_sat;

        for (i = 0; i < rttbl->NumStc - rttbl->NumMin; i++)
        {
            tmpval = 0.0;
            for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
            {
                tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
            }
            totconc[i] = tmpval;
            residue[i] = tmpval - (Vcele->t_conc[i] +
                (Rate_spe[i] + Rate_spet[i]) * stepsize * 0.5);
        }
        if (control % SKIP_JACOB == 0)  /* update jacobian every the other iteration */
        {
            for (k = 0; k < rttbl->NumStc - rttbl->NumMin; k++)
            {
                tmpconc[k] += tmpprb;
                for (i = 0; i < rttbl->NumSsc; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->NumSdc; j++)
                        tmpval +=
                            (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
                    tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
                    tmpconc[i + rttbl->NumStc] = tmpval;
                }
                for (i = 0; i < rttbl->NumStc - rttbl->NumMin; i++)
                {
                    tmpval = 0.0;
                    for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
                    {
                        tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
                    }
                    residue_t[i] = tmpval - (Vcele->t_conc[i] +
                        (Rate_spe[i] + Rate_spet[i]) * stepsize * 0.5);
                    jcb[k][i] = (residue_t[i] - residue[i]) * tmpprb_inv;
                }
                tmpconc[k] -= tmpprb;
            }
        }
        for (i = 0; i < rttbl->NumStc - rttbl->NumMin; i++)
            x_[i] = -residue[i];

        pivot_flg = denseGETRF(jcb, rttbl->NumStc - rttbl->NumMin, rttbl->NumStc - rttbl->NumMin, p);   // 09.17
        if (pivot_flg != 0)
        {
            Vcele->illness++;
            return (1);
        }

        denseGETRS(jcb, rttbl->NumStc - rttbl->NumMin, p, x_);

        for (i = 0; i < rttbl->NumStc - rttbl->NumMin; i++)
        {
            if (fabs(x_[i]) < 0.3)
                tmpconc[i] += x_[i];
            else
            {
                if (x_[i] < 0)
                    tmpconc[i] += -0.3;
                else
                    tmpconc[i] += 0.3;
            }
            error[i] = residue[i] / totconc[i];
        }
        maxerror = fabs(error[0]);
        for (i = 1; i < rttbl->NumStc - rttbl->NumMin; i++)
            if (fabs(error[i]) > maxerror)
                maxerror = fabs(error[i]);
        control++;
        if (control > 10)
            return (1);
    }

    destroyMat(jcb);

    for (i = 0; i < rttbl->NumSsc; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->NumSdc; j++)
        {
            tmpval += (tmpconc[j] + gamma[j]) * CD->Dependency[i][j];
        }
        tmpval -= Keq[i] + gamma[i + rttbl->NumStc];
        tmpconc[i + rttbl->NumStc] = tmpval;
    }

    for (i = 0; i < rttbl->NumStc - rttbl->NumMin; i++)
    {
        tmpval = 0.0;
        for (j = 0; j < rttbl->NumStc + rttbl->NumSsc; j++)
        {
            tmpval += CD->Totalconc[i][j] * pow(10, tmpconc[j]);
        }
        totconc[i] = tmpval;
        residue[i] = tmpval - Vcele->t_conc[i];
        error[i] = residue[i] / totconc[i];
    }
    for (i = 0; i < rttbl->NumStc + rttbl->NumSsc; i++)
    {
        if (i < rttbl->NumStc)
        {
            if (chemtbl[i].itype == MINERAL)
            {
                Vcele->t_conc[i] +=
                    (Rate_spe[i] + Rate_spet[i]) * stepsize * 0.5;
                Vcele->p_actv[i] = 1.0;
                Vcele->p_conc[i] = Vcele->t_conc[i];
            }
            else
            {
                Vcele->p_conc[i] = pow(10, tmpconc[i]);
                Vcele->p_actv[i] = pow(10, (tmpconc[i] + gamma[i]));
                Vcele->t_conc[i] = totconc[i];
            }
        }
        else
        {
            Vcele->s_conc[i - rttbl->NumStc] = pow(10, tmpconc[i]);
#if TEMP_DISABLED
            Vcele->s_actv[i - rttbl->NumStc] =
                pow(10, (tmpconc[i] + gamma[i]));
#endif
        }
    }

    free(residue);
    free(residue_t);
    free(tmpconc);
    free(totconc);
    free(area);
    free(error);
    free(gamma);
    free(Keq);
    free(Rate_pre);
    free(IAP);
    free(dependency);
    free(Rate_spe);
    free(Rate_spe_t);
    free(Rate_spet);
    free(p);
    free(x_);
    return (0);
}

int keymatch(const char *line, const char *keyword, double *value, char **strval)
{
    /* A very general and convinient way of reading datafile and input file */
    /* find keyword in line, assign the value after keyword to value array if there is any */
    /* store both numbers and strings in order for later use, buffer required */
    /* if is keyword not found return 0. If comments, return 2. Otherwise return 1 */
    int             i;

    for (i = 0; i < WORDS_LINE; i++)
        value[i] = 0.0;

    if ((line[0] == '!') || (line[0] == '#'))
    {
        /* assign a special flag for comments */
        return (2);
    }

    int             j, k;
    int             words_line = WORDS_LINE;
    int             keyfoundflag = 0;

    char          **words;
    words = (char **)malloc(WORDS_LINE * sizeof(char *));

    for (i = 0; i < WORDS_LINE; i++)
    {
        words[i] = (char *)malloc(WORD_WIDTH * sizeof(char));
        memset(words[i], 0, WORD_WIDTH);
    }
    i = j = k = 0;

    /* Partition the line into words */
    while (i < (int)strlen(line))
    {
        if (line[i] != 39)
        {
            while (line[i] != 9 && line[i] != 0 && line[i] != 10
                && line[i] != 32 && line[i] != 13)
            {
                words[k][j++] = line[i++];
                if (line[i] == 9 || line[i] == 32 || line[i] == 13)
                {
                    k++;
                    j = 0;
                }
            }
        }
        else
        {
            words[k][j++] = line[i++];
            while (line[i] != 39)
            {
                words[k][j++] = line[i++];
            }
            words[k++][j] = line[i++];
            j = 0;
        }
        i++;
    }

    words_line = k + 1;

    for (i = 0; i < words_line; i++)
        if (strcmp(words[i], keyword) == 0)
            keyfoundflag = 1;

    j = k = 0;
    for (i = 0; i < words_line; i++)
    {
        strcpy(strval[k++], words[i]);
        if (realcheck(words[i]) == 1)
            value[j++] = atof(words[i]);
    }

    for (i = 0; i < WORDS_LINE; i++)
        free(words[i]);
    free(words);
    return (keyfoundflag);

}

int realcheck(const char *words)
{
    int             flg = 1, i;
    if (((words[0] >= '0') && (words[0] <= '9')) ||
        (words[0] == '.') || (words[0] == '-') || (words[0] == '+'))
    {
        for (i = 0; i < (int)strlen(words); i++)
        {
            /* Ascii 10 is new line and 13 is carriage return */
            if ((words[i] > '9' || words[i] < '+') && (words[i] != 'E')
                && (words[i] != 'e') && (words[i] != 10) && (words[i] != 13))
            {
                flg = 0;
            }
        }
    }
    else
    {
        flg = 0;
    }
    return (flg);
}

