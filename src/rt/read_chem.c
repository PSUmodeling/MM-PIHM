#include "pihm.h"

#define ZERO   1E-20

void ReadChem(const char chem_filen[], const char cdbs_filen[],
    const pihm_struct pihm, Chem_Data CD)
{
    int             i, j;
    int             match;
    int             chem_ind;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    char            chemn[MAXSPS][MAXSTRING];
    int             p_type[MAXSPS];
    int             lno = 0;
    FILE           *chem_fp;
    FILE           *database;

    chem_fp = fopen(chem_filen, "r");
    CheckFile(chem_fp, chem_filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", chem_filen);

    database = fopen(cdbs_filen, "r");
    CheckFile(database, cdbs_filen);

    /* Default control variable if not found in input file */
    CD->TVDFlg = 1;
    CD->Cementation = 1.0;
    CD->ACTmod = 0;
    CD->DHEdel = 0;
    CD->TEMcpl = 0;
    CD->EffAds = 0;
    CD->RelMin = 0;
    CD->AvgScl = 10;
    CD->Condensation = 1.0;
    CD->NumBTC = 0;
    CD->NumPUMP = 0;
    CD->SUFEFF = 1;
    CD->CnntVelo = 0.01;

    /*
     * Runtime block
     */
    FindLine(chem_fp, "RUNTIME", &lno, chem_filen);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_TYPE", &CD->conc_init, 'i', chem_filen, lno);
    if (debug_mode)
    {
        switch (CD->conc_init)
        {
            case 0:
                PIHMprintf(VL_NORMAL,
                    "  Concentration will be initialized using .cini\n");
                break;
            case 1:
                PIHMprintf(VL_NORMAL,
                    "  Concentration will be initialized using .rtic\n");
                break;
            default:
                break;
        }
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TVD", temp_str, 'w', chem_filen, lno);
    if (strcasecmp(temp_str, "false") == 0)
    {
        CD->TVDFlg = 0;
    }
    else if (strcasecmp(temp_str, "true") == 0)
    {
        CD->TVDFlg = 1;
    }
    else
    {
        PIHMprintf(VL_ERROR, "  TVD flag input error.\n");
        PIHMexit(EXIT_FAILURE);
    }
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  Total variation diminishing set to %d %s.\n",
                CD->TVDFlg, temp_str);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ACTIVITY", &CD->ACTmod, 'i', chem_filen, lno);
    if (debug_mode)
    {
        /* 0 for unity activity coefficient and 1 for DH equation update */
        PIHMprintf(VL_NORMAL, "  Activity correction is set to %d. \n",
            CD->ACTmod);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ACT_COE_DELAY", &CD->DHEdel, 'i', chem_filen, lno);
    if (debug_mode)
    {
        /* 0 for delay and 1 for no delay (solving together) */
        PIHMprintf(VL_NORMAL,
            "  Activity coefficient update delay is set to %d. \n", CD->DHEdel);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "THERMO", &CD->TEMcpl, 'i', chem_filen, lno);
    if (debug_mode)
    {
        /* 0 for delay and 1 for no delay (solving together) */
        PIHMprintf(VL_NORMAL, "  Coupling of thermo modelling is set to %d. \n",
            CD->DHEdel);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RELMIN", &CD->RelMin, 'i', chem_filen, lno);
    if (debug_mode)
    {
        switch (CD->RelMin)
        {
            case 0:
                PIHMprintf(VL_NORMAL,
                    "  Using absolute mineral volume fraction. \n");
                break;
            case 1:
                PIHMprintf(VL_NORMAL,
                    "  Using relative mineral volume fraction. \n");
                break;
            default:
                break;
        }
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "EFFADS", &CD->EffAds, 'i', chem_filen, lno);
    if (debug_mode)
    {
        switch (CD->EffAds)
        {
            case 0:
                PIHMprintf(VL_NORMAL,
                    "  Using the normal adsorption model. \n");
                break;
            case 1:
                PIHMprintf(VL_NORMAL,
                    "  Using the coupled MIM and adsorption model. \n");
                break;
                /* under construction. */
            default:
                break;
        }
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TRANSPORT_ONLY", &CD->RecFlg, 'i', chem_filen, lno);
    if (debug_mode)
    {
        switch (CD->RecFlg)
        {
            case 0:
                PIHMprintf(VL_NORMAL, "  Transport only mode disabled.\n");
                break;
            case 1:
                PIHMprintf(VL_NORMAL, "  Transport only mode enabled. \n");
                break;
                /* under construction. */
            default:
                break;
        }
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRECIPITATION", &CD->PrpFlg, 'i', chem_filen, lno);
    if (debug_mode)
    {
        switch (CD->PrpFlg)
        {
            case 0:
                PIHMprintf(VL_NORMAL, "  No precipitation condition. \n");
                break;
            case 1:
                PIHMprintf(VL_NORMAL,
                    "  Precipitation condition is to be specified. \n");
                break;
            case 2:
                PIHMprintf(VL_NORMAL,
                    "  Precipitation condition is specified via file *.prep. \n");
                break;
                /* under construction. */
            default:
                break;
        }
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RT_DELAY", &CD->RT_delay, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL,
            "  Flux-PIHM-RT will start after running PIHM for %d days. \n",
            CD->RT_delay);
        /* under construction. */
    }
    CD->RT_delay *= DAYINSEC;

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CONDENSATION", &CD->Condensation, 'd', chem_filen,
        lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  The concentrations of infiltrating rainfall "
            "is set to be %f times of concentrations in precipitation. \n",
            CD->Condensation);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "AVGSCL", &CD->AvgScl, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL,
            "  Averaging window for asynchronous reaction %d. \n",
            CD->AvgScl);
        /* under construction. */
    }

    /*
     * Global block
     */
    FindLine(chem_fp, "GLOBAL", &lno, chem_filen);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "T_SPECIES", &CD->NumStc, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d chemical species specified. \n", CD->NumStc);
        /* H2O is always a primary species */
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "S_SPECIES", &CD->NumSsc, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d secondary species specified. \n",
            CD->NumSsc);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_SPECIES", &CD->NumMin, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d minerals specified. \n", CD->NumMin);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ADSORPTION", &CD->NumAds, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d surface complexation specified. \n",
            CD->NumAds);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CATION_EXCHANGE", &CD->NumCex, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d cation exchange specified. \n", CD->NumCex);
    }

    /* The number of species that are mobile, later used in the OS3D subroutine */
    NumSpc = CD->NumStc - (CD->NumMin + CD->NumAds + CD->NumCex);

    /* The number of species that others depend on */
    CD->NumSdc = CD->NumStc - CD->NumMin;

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MINERAL_KINETIC", &CD->NumMkr, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d mineral kinetic reaction(s) specified. \n",
            CD->NumMkr);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "AQUEOUS_KINETIC", &CD->NumAkr, 'i', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  %d aqueous kinetic reaction(s) specified. \n",
            CD->NumAkr);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DIFFUSION", &CD->DiffCoe, 'd', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  Diffusion coefficient = %g [cm2/s] \n",
            CD->DiffCoe);
    }
    CD->DiffCoe /= 1.0E4;       /* Convert from cm2 s-1 to m2 s-1 */

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DISPERSION", &CD->DispCoe, 'd', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  Dispersion coefficient = %2.2f [m] \n",
            CD->DispCoe);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CEMENTATION", &CD->Cementation, 'd', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  Cementation factor = %2.1f \n", CD->Cementation);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TEMPERATURE", &CD->Temperature, 'd', chem_filen, lno);
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, "  Temperature = %3.1f \n\n", CD->Temperature);
    }

    /*
     * Output block
     */
    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "OUTPUT", &CD->NumBTC, 'i', chem_filen, lno);
    CD->BTC_loc = (int *)malloc(CD->NumBTC * sizeof(int));
    for (i = 0; i < CD->NumBTC; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        sscanf(cmdstr, "%d", &CD->BTC_loc[i]);
    }
    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL,
            "  %d breakthrough points specified. \n", CD->NumBTC);
        PIHMprintf(VL_NORMAL, "  --");
        for (i = 0; i < CD->NumBTC; i++)
        {
            PIHMprintf(VL_NORMAL, " Grid %d ", CD->BTC_loc[i] + 1);
        }
        PIHMprintf(VL_NORMAL, "are breakthrough points.\n\n");
    }

    /* Pump block */
    CD->CalGwinflux = pihm->cal.gwinflux;

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PUMP", &CD->NumPUMP, 'i', chem_filen, lno);
    CD->pumps = (Pump *) malloc(CD->NumPUMP * sizeof(Pump));
    for (i = 0; i < CD->NumPUMP; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%d %s %lf %lf",
            &CD->pumps[i].Pump_Location,
            CD->pumps[i].Name_Species,
            &CD->pumps[i].Injection_rate,
            &CD->pumps[i].Injection_conc) != 4)
        {
            PIHMprintf(VL_NORMAL, "Error reading pump information.\n");
            PIHMexit(EXIT_FAILURE);
        }
        CD->pumps[i].flow_rate =
            CD->pumps[i].Injection_rate / CD->pumps[i].Injection_conc /
            365 * 1E-3;
        CD->pumps[i].Injection_rate *=
            CD->pumps[i].Injection_rate * CD->CalGwinflux;
        CD->pumps[i].flow_rate = CD->pumps[i].flow_rate * CD->CalGwinflux;

        if (debug_mode)
        {
            PIHMprintf(VL_NORMAL,
                "  -- Rate %g [moles/year] of '%s' (pos: %d) at Grid '%d' with a concentration of %g [M/L]. \n",
                CD->pumps[i].Injection_rate, CD->pumps[i].Name_Species,
                (CD->pumps[i].Position_Species + 1), CD->pumps[i].Pump_Location,
                CD->pumps[i].Injection_conc);
            PIHMprintf(VL_NORMAL, "  -- Flow rate is then %g [m3/d]. \n",
                CD->pumps[i].flow_rate);

        }
    }

    /*
     * Precipitation block
     *
     * The precipitation block is read twice. The first time chemical names are
     * read and then sort to order. The second time the concentrations and SSA's
     * are read.
     */
    /* Read the first time */
    FindLine(chem_fp, "PRECIPITATION_CONC", &lno, chem_filen);
    for (i = 0; i < CD->NumStc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        sscanf(cmdstr, "%s", chemn[i]);
        p_type[i] = SpeciationType(database, chemn[i]);
    }

    SortChem(chemn, p_type, CD->NumStc, CD->chemtype);

    if (debug_mode)
    {
        PIHMprintf(VL_NORMAL, " \n");
        PIHMprintf(VL_NORMAL, "  ---------------------------------\n");
        PIHMprintf(VL_NORMAL, "  The condition of precipitation is \n");
        PIHMprintf(VL_NORMAL, "  ---------------------------------\n");
    }

    CD->Precipitation.t_conc =
        (double *)malloc(CD->NumStc * sizeof(double));
    CD->Precipitation.p_conc =
        (double *)malloc(CD->NumStc * sizeof(double));
    CD->Precipitation.p_para =
        (double *)malloc(CD->NumStc * sizeof(double));
    CD->Precipitation.s_conc = NULL;
    for (i = 0; i < CD->NumStc; i++)
    {
        CD->Precipitation.t_conc[i] = ZERO;
        CD->Precipitation.p_conc[i] = ZERO;
    }

    FindLine(chem_fp, "BOF", &lno, chem_filen);
    FindLine(chem_fp, "PRECIPITATION_CONC", &lno, chem_filen);
    for (i = 0; i < CD->NumStc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        chem_ind = FindChem(chemn[i], CD->chemtype, CD->NumStc);

        if (CD->chemtype[chem_ind].itype == AQUEOUS)
        {
            sscanf(cmdstr, "%*s %lf", &CD->Precipitation.t_conc[chem_ind]);
            PIHMprintf(VL_NORMAL, "  %-28s %g \n",
                CD->chemtype[chem_ind].ChemName, CD->Precipitation.t_conc[chem_ind]);

            if (strcasecmp(CD->chemtype[chem_ind].ChemName, "pH") == 0)
            {
                /* Change the pH of precipitation into total concentraion of H
                 * Skip the speciation for rain */
                CD->Precipitation.t_conc[chem_ind] =
                    (CD->Precipitation.t_conc[chem_ind] < 7.0) ?
                    pow(10, -CD->Precipitation.t_conc[chem_ind]) :
                    -pow(10, CD->Precipitation.t_conc[chem_ind] - 14);
            }
        }
        if (CD->chemtype[chem_ind].itype == MINERAL)
        {
            sscanf(cmdstr, "%*s %lf %s %lf", &CD->Precipitation.t_conc[chem_ind],
                temp_str, &CD->Precipitation.p_para[chem_ind]);
            PIHMprintf(VL_NORMAL,
                "  mineral %-20s %6.4f \t specific surface area %6.4f\n",
                CD->chemtype[chem_ind].ChemName, CD->Precipitation.t_conc[chem_ind],
                CD->Precipitation.p_para[chem_ind]);
        }
        if (CD->chemtype[chem_ind].ChemName[0] == '>' ||
            CD->chemtype[chem_ind].itype == ADSORPTION)
        {
            /* adsorptive sites and species start with > */
            sscanf(cmdstr, "%*s %lf", &CD->Precipitation.t_conc[chem_ind]);
            CD->Precipitation.p_para[chem_ind] = 0;
            PIHMprintf(VL_NORMAL, " surface complex %s\t %6.4f\n",
                CD->chemtype[chem_ind].ChemName, CD->Precipitation.t_conc[chem_ind]);
        }
        if (CD->chemtype[chem_ind].itype == CATION_ECHG)
        {
            sscanf(cmdstr, "%*s %lf", &CD->Precipitation.t_conc[chem_ind]);
            CD->Precipitation.p_para[chem_ind] = 0;
            PIHMprintf(VL_NORMAL, " cation exchange %s\t %6.4f\n",
                CD->chemtype[chem_ind].ChemName, CD->Precipitation.t_conc[chem_ind]);
        }
    }

    /*
     * Secondary_species block
     */
    FindLine(chem_fp, "SECONDARY_SPECIES", &lno, chem_filen);
    for (i = 0; i < CD->NumSsc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s", CD->chemtype[CD->NumStc + i].ChemName) != 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading secondary_species in %s near Line %d.\n",
                chem_filen, lno);
        }
    }

    /* Minerals block */
    FindLine(chem_fp, "MINERALS", &lno, chem_filen);

    for (i = 0; i < MAXSPS; i++)
    {
        for (j = 0; j < MAXDEP; j++)
        {
            CD->kinetics[i].dep_position[j] = 0;
            CD->kinetics[i].monod_position[j] = 0;
            CD->kinetics[i].inhib_position[j] = 0;
        }
    }

    for (i = 0; i < CD->NumMkr; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%s %s %s",
            CD->kinetics[i].species, temp_str, CD->kinetics[i].Label);
        if (match != 3 || strcasecmp(temp_str, "-label") != 0)
        {
            PIHMprintf(VL_ERROR,
                "Error reading mineral information in %s near Line %d.\n",
                chem_filen, lno);
            PIHMexit(EXIT_FAILURE);
        }

        if (debug_mode)
        {
            PIHMprintf(VL_NORMAL,
                "  Kinetic reaction on '%s' is specified, label '%s'. \n",
                CD->kinetics[i].species, CD->kinetics[i].Label);
        }
    }

    fclose(chem_fp);
    fclose(database);
}
