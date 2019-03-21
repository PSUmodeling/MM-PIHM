#include "pihm.h"

#define ZERO   1E-20

void ReadChem(const char chem_filen[], const char cdbs_filen[],
    const pihm_struct pihm, chemtbl_struct chemtbl[], kintbl_struct kintbl[],
    rttbl_struct *rttbl, ctrl_struct *ctrl, Chem_Data CD)
{
    int             i, j;
    int             match;
    int             chem_ind;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    char            chemn[MAXSPS][MAXSTRING];
    int             p_type[MAXSPS];
    int             lno = 0;
    double          conc;
    FILE           *chem_fp;
    FILE           *db_fp;

    chem_fp = fopen(chem_filen, "r");
    CheckFile(chem_fp, chem_filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", chem_filen);

    db_fp = fopen(cdbs_filen, "r");
    CheckFile(db_fp, cdbs_filen);

    /* Default control variable if not found in input file */
    rttbl->Cementation = 1.0;
    rttbl->ACTmod = 0;
    rttbl->TEMcpl = 0;
    CD->EffAds = 0;
    rttbl->RelMin = 0;
    ctrl->AvgScl = 10;
    rttbl->Condensation = 1.0;
    CD->NumBTC = 0;
    CD->NumPUMP = 0;
    CD->SUFEFF = 1;
    CD->CnntVelo = 0.01;

    /*
     * Runtime block
     */
    PIHMprintf(VL_VERBOSE, "\n Runtime block\n");
    FindLine(chem_fp, "RUNTIME", &lno, chem_filen);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_TYPE", &ctrl->read_rt_restart, 'i',
        chem_filen, lno);
    switch (ctrl->read_rt_restart)
    {
        case 0:
            PIHMprintf(VL_VERBOSE,
                    "  Concentration will be initialized using .cini\n");
            break;
        case 1:
            PIHMprintf(VL_VERBOSE,
                    "  Concentration will be initialized using .rtic\n");
            break;
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ACTIVITY", &rttbl->ACTmod, 'i', chem_filen, lno);
    /* 0 for unity activity coefficient and 1 for DH equation update */
    PIHMprintf(VL_VERBOSE, "  Activity correction is set to %d. \n",
        rttbl->ACTmod);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "THERMO", &rttbl->TEMcpl, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  Coupling of thermo modelling is set to %d. \n",
        rttbl->TEMcpl);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RELMIN", &rttbl->RelMin, 'i', chem_filen, lno);
    switch (rttbl->RelMin)
    {
        case 0:
            PIHMprintf(VL_VERBOSE,
                "  Using absolute mineral volume fraction. \n");
            break;
        case 1:
            PIHMprintf(VL_VERBOSE,
                "  Using relative mineral volume fraction. \n");
            break;
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TRANSPORT_ONLY", &rttbl->RecFlg, 'i', chem_filen, lno);
    switch (rttbl->RecFlg)
    {
        case 0:
            PIHMprintf(VL_VERBOSE, "  Transport only mode disabled.\n");
            break;
        case 1:
            PIHMprintf(VL_VERBOSE, "  Transport only mode enabled. \n");
            break;
            /* under construction. */
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRECIPITATION", &ctrl->PrpFlg, 'i', chem_filen, lno);
    switch (ctrl->PrpFlg)
    {
        case 0:
            PIHMprintf(VL_VERBOSE, "  No precipitation condition. \n");
            break;
        case 1:
            PIHMprintf(VL_VERBOSE,
                "  Precipitation condition is to be specified. \n");
            break;
        case 2:
            PIHMprintf(VL_VERBOSE,
                "  Precipitation condition is specified via file *.prep. \n");
            break;
            /* under construction. */
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RT_DELAY", &ctrl->RT_delay, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE,
        "  Flux-PIHM-RT will start after running PIHM for %d seconds. \n",
        ctrl->RT_delay);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CONDENSATION", &rttbl->Condensation, 'd',
        chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  The concentrations of infiltrating rainfall "
        "is set to be %f times of concentrations in precipitation. \n",
        rttbl->Condensation);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "AVGSCL", &ctrl->AvgScl, 'i', chem_filen, lno);
    if (ctrl->AvgScl < ctrl->stepsize || ctrl->AvgScl % ctrl->stepsize > 0)
    {
        PIHMprintf(VL_ERROR,
            "Error: Reaction step size "
            "should be an integral multiple of model step size.\n");
        PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", chem_filen, lno);
        PIHMexit(EXIT_FAILURE);
    }
    else
    {
        PIHMprintf(VL_VERBOSE,
            "  Averaging window for asynchronous reaction %d seconds.\n",
            ctrl->AvgScl);
    }

    /*
     * Global block
     */
    PIHMprintf(VL_VERBOSE, "\n Global block\n");
    FindLine(chem_fp, "GLOBAL", &lno, chem_filen);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "T_SPECIES", &rttbl->NumStc, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d chemical species specified. \n", rttbl->NumStc);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "S_SPECIES", &rttbl->NumSsc, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d secondary species specified. \n", rttbl->NumSsc);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_SPECIES", &rttbl->NumMin, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d minerals specified. \n", rttbl->NumMin);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ADSORPTION", &rttbl->NumAds, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d surface complexation specified. \n", rttbl->NumAds);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CATION_EXCHANGE", &rttbl->NumCex, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d cation exchange specified. \n", rttbl->NumCex);

    /* The number of species that are mobile */
    NumSpc = rttbl->NumStc - (rttbl->NumMin + rttbl->NumAds + rttbl->NumCex);
    if (NumSpc <= 0)
    {
        PIHMprintf(VL_ERROR,
            "Error: number of total species should be larger than the sum of "
            " mineral, surface complexation, and cation exchange species.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /* The number of species that others depend on */
    rttbl->NumSdc = rttbl->NumStc - rttbl->NumMin;
    if (rttbl->NumSdc < 0)
    {
        PIHMprintf(VL_ERROR,
            "Error: number of total species should not be smaller than number "
            "of minerals.\n");
        PIHMexit(EXIT_FAILURE);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MINERAL_KINETIC", &rttbl->NumMkr, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d mineral kinetic reaction(s) specified. \n",
        rttbl->NumMkr);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "AQUEOUS_KINETIC", &rttbl->NumAkr, 'i', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  %d aqueous kinetic reaction(s) specified. \n",
        rttbl->NumAkr);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DIFFUSION", &rttbl->DiffCoe, 'd', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  Diffusion coefficient = %g cm2 s-1 \n", rttbl->DiffCoe);
    rttbl->DiffCoe /= 1.0E4;       /* Convert from cm2 s-1 to m2 s-1 */

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DISPERSION", &rttbl->DispCoe, 'd', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  Dispersion coefficient = %2.2f m \n", rttbl->DispCoe);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CEMENTATION", &rttbl->Cementation, 'd', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  Cementation factor = %2.1f \n", rttbl->Cementation);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TEMPERATURE", &rttbl->Temperature, 'd', chem_filen, lno);
    PIHMprintf(VL_VERBOSE, "  Temperature = %3.1f \n\n", rttbl->Temperature);

    /*
     * Primary species block
     */
    PIHMprintf(VL_VERBOSE, "\n Primary species block\n");
    FindLine(chem_fp, "PRIMARY_SPECIES", &lno, chem_filen);
    for (i = 0; i < rttbl->NumStc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s", chemn[i]) != 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading primary_species in %s near Line %d.\n",
                chem_filen, lno);
        }
        p_type[i] = SpeciationType(db_fp, chemn[i]);
    }

    SortChem(chemn, p_type, rttbl->NumStc, chemtbl);

    /*
     * Secondary_species block
     */
    PIHMprintf(VL_VERBOSE, "\n Secondary species block\n");
    FindLine(chem_fp, "SECONDARY_SPECIES", &lno, chem_filen);
    for (i = 0; i < rttbl->NumSsc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s", chemtbl[rttbl->NumStc + i].ChemName) != 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading secondary_species in %s near Line %d.\n",
                chem_filen, lno);
        }
    }

    /*
     * Minerals block
     */
    PIHMprintf(VL_VERBOSE, "\n Minerals block\n");
    FindLine(chem_fp, "MINERALS", &lno, chem_filen);
    for (i = 0; i < rttbl->NumMkr; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s %*s %s", temp_str, kintbl[i].Label) != 2)
        {
            PIHMprintf(VL_ERROR,
                "Error reading mineral information in %s near Line %d.\n",
                chem_filen, lno);
            PIHMexit(EXIT_FAILURE);
        }

        PIHMprintf(VL_VERBOSE,
            "  Kinetic reaction on '%s' is specified, label '%s'.\n",
            temp_str, kintbl[i].Label);

        kintbl[i].position = FindChem(temp_str, chemtbl, rttbl->NumStc);

        if (kintbl[i].position < 0)
        {
            PIHMprintf(VL_ERROR,
                "Error finding mineral %s in species table.\n", temp_str);
            PIHMexit(EXIT_FAILURE);
        }
        else
        {
            PIHMprintf(VL_VERBOSE,
                "  Position_check (NumMkr[i] vs NumStc[j]) (%d, %d)\n",
                i, kintbl[i].position);
        }
    }

    /*
     * Precipitation block
     */
    PIHMprintf(VL_VERBOSE, "\n Precipitation concentration block\n");
    FindLine(chem_fp, "PRECIPITATION_CONC", &lno, chem_filen);
    PIHMprintf(VL_VERBOSE, "  ---------------------------------\n");
    PIHMprintf(VL_VERBOSE, "  The condition of precipitation is \n");
    PIHMprintf(VL_VERBOSE, "  ---------------------------------\n");

    for (i = 0; i < rttbl->NumStc; i++)
    {
        rttbl->prcp_conc[i] = ZERO;
    }

    for (i = 0; i < rttbl->NumStc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s %lf", temp_str, &conc) != 2)
        {
            PIHMprintf(VL_ERROR, "Error reading precipitation concentration "
                "in %s near Line %d.\n", chem_filen, lno);
        }
        chem_ind = FindChem(temp_str, chemtbl, rttbl->NumStc);
        if (chem_ind < 0)
        {
            PIHMprintf(VL_ERROR, "Error finding chemical %s.\n", chemn[i]);
            PIHMexit(EXIT_FAILURE);
        }
        else
        {
            PIHMprintf(VL_NORMAL, "  %-28s %g \n",
                chemtbl[chem_ind].ChemName, conc);

            if (strcasecmp(chemtbl[chem_ind].ChemName, "pH") == 0)
            {
                /* Change the pH of precipitation into total concentraion of H
                 * Skip the speciation for rain */
                conc = (conc < 7.0) ?  pow(10, -conc) : -pow(10, conc - 14);
            }

            rttbl->prcp_conc[chem_ind] = conc;
        }
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


    fclose(chem_fp);
    fclose(db_fp);
}
