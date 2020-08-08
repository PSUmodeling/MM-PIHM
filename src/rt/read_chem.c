#include "pihm.h"

void ReadChem(const char chem_fn[], const char cdbs_fn[],
    chemtbl_struct chemtbl[], kintbl_struct kintbl[], rttbl_struct *rttbl,
    forc_struct *forc, ctrl_struct *ctrl)
{
    int             i;
    int             chem_ind;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    char            chemn[MAXSPS][MAXSTRING];
    int             p_type[MAXSPS];
    int             lno = 0;
    double          conc;
    FILE           *chem_fp;
    FILE           *db_fp;

    chem_fp = pihm_fopen(chem_fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", chem_fn);

    db_fp = pihm_fopen(cdbs_fn, "r");

    /*
     * Runtime block
     */
    pihm_printf(VL_VERBOSE, "\n Runtime block\n");
    FindLine(chem_fp, "RUNTIME", &lno, chem_fn);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_TYPE", 'i', chem_fn, lno, &ctrl->read_rt_restart);
    switch (ctrl->read_rt_restart)
    {
        case 0:
            pihm_printf(VL_VERBOSE,
                "  Concentrations will be initialized using .cini\n");
            break;
        case 1:
            pihm_printf(VL_VERBOSE,
                "  Concentrations will be initialized using .rtic\n");
            break;
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", 'i', chem_fn, lno, &ctrl->write_rt_restart);
    if (ctrl->write_rt_restart)
    {
        pihm_printf(VL_VERBOSE,
            "  Concentrations will be written into .rtic\n");
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ACTIVITY", 'i', chem_fn, lno, &rttbl->actv_mode);
    /* 0 for unity activity coefficient and 1 for DH equation update */
    pihm_printf(VL_VERBOSE, "  Activity correction is set to %d. \n",
        rttbl->actv_mode);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "THERMO", 'i', chem_fn, lno, &rttbl->tmp_coup);
    pihm_printf(VL_VERBOSE, "  Coupling of thermo modelling is set to %d. \n",
        rttbl->tmp_coup);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RELMIN", 'i', chem_fn, lno, &rttbl->rel_min);
    switch (rttbl->rel_min)
    {
        case 0:
            pihm_printf(VL_VERBOSE,
                "  Using absolute mineral volume fraction. \n");
            break;
        case 1:
            pihm_printf(VL_VERBOSE,
                "  Using relative mineral volume fraction. \n");
            break;
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TRANSPORT_ONLY", 'i', chem_fn, lno, &rttbl->transpt_flag);
    switch (rttbl->transpt_flag)
    {
        case KIN_REACTION:
            pihm_printf(VL_VERBOSE, "  Transport only mode disabled.\n");
            break;
        case TRANSPORT_ONLY:
            pihm_printf(VL_VERBOSE, "  Transport only mode enabled. \n");
            break;
            /* under construction. */
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRECIPITATION", 'i', chem_fn, lno, &forc->prcp_flag);
    switch (forc->prcp_flag)
    {
        case 0:
            pihm_printf(VL_VERBOSE, "  No precipitation condition. \n");
            break;
        case 1:
            pihm_printf(VL_VERBOSE,
                "  Precipitation condition is to be specified. \n");
            break;
        case 2:
            pihm_printf(VL_VERBOSE,
                "  Precipitation condition is specified via file *.prep. \n");
            break;
        default:
            break;
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RT_DELAY", 'i', chem_fn, lno, &ctrl->RT_delay);
    pihm_printf(VL_VERBOSE,
        "  Flux-PIHM-RT will start after running PIHM for %d seconds. \n",
        ctrl->RT_delay);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CONDENSATION", 'd', chem_fn, lno, &rttbl->cond);
    pihm_printf(VL_VERBOSE, "  The concentrations of infiltrating rainfall "
        "is set to be %f times of concentrations in precipitation. \n",
        rttbl->cond);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "AVGSCL", 'i', chem_fn, lno, &ctrl->AvgScl);
    if (ctrl->AvgScl < ctrl->stepsize || ctrl->AvgScl % ctrl->stepsize > 0)
    {
        pihm_printf(VL_ERROR,
            "Error: Reaction step size "
            "should be an integral multiple of model step size.\n");
        pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", chem_fn, lno);
        pihm_exit(EXIT_FAILURE);
    }
    else
    {
        pihm_printf(VL_VERBOSE,
            "  Averaging window for asynchronous reaction %d seconds.\n",
            ctrl->AvgScl);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ctrl->prtvrbl[CHEM_CTRL] = ReadPrintCtrl(cmdstr, "OUTINTVL", chem_fn, lno);
    if (ctrl->prtvrbl[CHEM_CTRL] > 0)
    {
        pihm_printf(VL_VERBOSE, "  Chemical concentration output interval "
            "is set to %d seconds.\n", ctrl->prtvrbl[CHEM_CTRL]);
    }
    else if (ctrl->prtvrbl[CHEM_CTRL] < 0)
    {
        switch (ctrl->prtvrbl[CHEM_CTRL])
        {
            case -1:
                pihm_printf(VL_VERBOSE,
                    "  Chemical concentration output is yearly.\n");
                break;
            case -2:
                pihm_printf(VL_VERBOSE,
                    "  Chemical concentration output is monthly.\n");
                break;
            case -3:
                pihm_printf(VL_VERBOSE,
                    "  Chemical concentration output is daily.\n");
                break;
            case -4:
                pihm_printf(VL_VERBOSE,
                    "  Chemical concentration output is hourly.\n");
                break;
        }
    }
    else
    {
        pihm_printf(VL_VERBOSE,
            "  Chemical concentration output is turned off.\n");
    }

    /*
     * Global block
     */
    pihm_printf(VL_VERBOSE, "\n Global block\n");
    FindLine(chem_fp, "GLOBAL", &lno, chem_fn);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "T_SPECIES", 'i', chem_fn, lno, &rttbl->num_stc);
    pihm_printf(VL_VERBOSE, "  %d chemical species specified. \n", rttbl->num_stc);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "S_SPECIES", 'i', chem_fn, lno, &rttbl->num_ssc);
    pihm_printf(VL_VERBOSE, "  %d secondary species specified. \n", rttbl->num_ssc);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_SPECIES", 'i', chem_fn, lno, &rttbl->num_min);
    pihm_printf(VL_VERBOSE, "  %d minerals specified. \n", rttbl->num_min);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ADSORPTION", 'i', chem_fn, lno, &rttbl->num_ads);
    pihm_printf(VL_VERBOSE, "  %d surface complexation specified. \n", rttbl->num_ads);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CATION_EXCHANGE", 'i', chem_fn, lno, &rttbl->num_cex);
    pihm_printf(VL_VERBOSE, "  %d cation exchange specified. \n", rttbl->num_cex);

    /* The number of species that are mobile */
    rttbl->num_spc = rttbl->num_stc - (rttbl->num_min + rttbl->num_ads + rttbl->num_cex);
    if (rttbl->num_spc <= 0)
    {
        pihm_printf(VL_ERROR,
            "Error: number of total species should be larger than the sum of "
            " mineral, surface complexation, and cation exchange species.\n");
        pihm_exit(EXIT_FAILURE);
    }
    nsolute = rttbl->num_spc;

    /* The number of species that others depend on */
    rttbl->num_sdc = rttbl->num_stc - rttbl->num_min;
    if (rttbl->num_sdc < 0)
    {
        pihm_printf(VL_ERROR,
            "Error: number of total species should not be smaller than number "
            "of minerals.\n");
        pihm_exit(EXIT_FAILURE);
    }

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MINERAL_KINETIC", 'i', chem_fn, lno, &rttbl->num_mkr);
    pihm_printf(VL_VERBOSE, "  %d mineral kinetic reaction(s) specified. \n",
        rttbl->num_mkr);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "AQUEOUS_KINETIC", 'i', chem_fn, lno, &rttbl->num_akr);
    pihm_printf(VL_VERBOSE, "  %d aqueous kinetic reaction(s) specified. \n",
        rttbl->num_akr);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DIFFUSION", 'd', chem_fn, lno, &rttbl->diff_coef);
    pihm_printf(VL_VERBOSE, "  Diffusion coefficient = %g cm2 s-1 \n", rttbl->diff_coef);
    rttbl->diff_coef /= 1.0E4;       /* Convert from cm2 s-1 to m2 s-1 */

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DISPERSION", 'd', chem_fn, lno, &rttbl->disp_coef);
    pihm_printf(VL_VERBOSE, "  Dispersion coefficient = %2.2f m \n", rttbl->disp_coef);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CEMENTATION", 'd', chem_fn, lno, &rttbl->cementation);
    pihm_printf(VL_VERBOSE, "  Cementation factor = %2.1f \n", rttbl->cementation);

    NextLine(chem_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TEMPERATURE", 'd', chem_fn, lno, &rttbl->tmp);
    pihm_printf(VL_VERBOSE, "  Temperature = %3.1f \n\n", rttbl->tmp);

    /*
     * Primary species block
     */
    pihm_printf(VL_VERBOSE, "\n Primary species block\n");
    FindLine(chem_fp, "PRIMARY_SPECIES", &lno, chem_fn);
    for (i = 0; i < rttbl->num_stc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s", chemn[i]) != 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading primary_species in %s near Line %d.\n",
                chem_fn, lno);
        }
        p_type[i] = SpeciesType(db_fp, chemn[i]);
    }

    SortChem(chemn, p_type, rttbl->num_stc, chemtbl);

    /*
     * Secondary_species block
     */
    pihm_printf(VL_VERBOSE, "\n Secondary species block\n");
    FindLine(chem_fp, "SECONDARY_SPECIES", &lno, chem_fn);
    for (i = 0; i < rttbl->num_ssc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s", chemtbl[rttbl->num_stc + i].name) != 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading secondary_species in %s near Line %d.\n",
                chem_fn, lno);
        }
    }

    /*
     * Minerals block
     */
    pihm_printf(VL_VERBOSE, "\n Minerals block\n");
    FindLine(chem_fp, "MINERALS", &lno, chem_fn);
    for (i = 0; i < rttbl->num_mkr; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s %*s %s", temp_str, kintbl[i].label) != 2)
        {
            pihm_printf(VL_ERROR,
                "Error reading mineral information in %s near Line %d.\n",
                chem_fn, lno);
            pihm_exit(EXIT_FAILURE);
        }

        pihm_printf(VL_VERBOSE,
            "  Kinetic reaction on '%s' is specified, label '%s'.\n",
            temp_str, kintbl[i].label);

        kintbl[i].position = FindChem(temp_str, chemtbl, rttbl->num_stc);

        if (kintbl[i].position < 0)
        {
            pihm_printf(VL_ERROR,
                "Error finding mineral %s in species table.\n", temp_str);
            pihm_exit(EXIT_FAILURE);
        }
        else
        {
            pihm_printf(VL_VERBOSE,
                "  Position_check (num_mkr[i] vs num_stc[j]) (%d, %d)\n",
                i, kintbl[i].position);
        }
    }

    /*
     * Precipitation block
     */
    pihm_printf(VL_VERBOSE, "\n Precipitation concentration block\n");
    FindLine(chem_fp, "PRECIPITATION_CONC", &lno, chem_fn);
    pihm_printf(VL_VERBOSE, "  ---------------------------------\n");
    pihm_printf(VL_VERBOSE, "  The condition of precipitation is \n");
    pihm_printf(VL_VERBOSE, "  ---------------------------------\n");

    for (i = 0; i < rttbl->num_stc; i++)
    {
        rttbl->prcp_conc[i] = 0.0;
    }

    for (i = 0; i < rttbl->num_stc; i++)
    {
        NextLine(chem_fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s %lf", temp_str, &conc) != 2)
        {
            pihm_printf(VL_ERROR, "Error reading precipitation concentration "
                "in %s near Line %d.\n", chem_fn, lno);
        }
        chem_ind = FindChem(temp_str, chemtbl, rttbl->num_stc);
        if (chem_ind < 0)
        {
            pihm_printf(VL_ERROR, "Error finding chemical %s.\n", chemn[i]);
            pihm_exit(EXIT_FAILURE);
        }
        else
        {
            pihm_printf(VL_VERBOSE, "  %-28s %g \n",
                chemtbl[chem_ind].name, conc);

            if (strcasecmp(chemtbl[chem_ind].name, "pH") == 0)
            {
                /* Change the pH of precipitation into total concentraion of H
                 * Skip the speciation for rain */
                conc = (conc < 7.0) ?  pow(10, -conc) : -pow(10, conc - 14);
            }

            rttbl->prcp_conc[chem_ind] = conc;
        }
    }

    fclose(chem_fp);
    fclose(db_fp);
}

int ParseLocation(const char str[], const char filen[], int lno)
{
    char            type_str[MAXSTRING];
    int             loc;

    if (strstr(str, ":") == 0)
    {
        pihm_printf(VL_ERROR,
            "Error: location format error in %s near Line %d.\n"
            "Please use \"UNSAT:XX\", \"GW:XX\", or \"RIVER:XX\", "
            "where XX is the index of grid.\n", filen, lno);
        pihm_exit(EXIT_FAILURE);
    }
    sscanf(str, "%[^:]:%d", type_str, &loc);
    if (strcasecmp(type_str, "UNSAT") == 0)
    {
        /* Do nothing */
    }
    else if (strcasecmp(type_str, "GW") == 0)
    {
        loc += nelem;
    }
    else if (strcasecmp(type_str, "RIVER") == 0)
    {
        loc *= -1;
    }
    else
    {
        pihm_printf(VL_ERROR,
            "Error: location format error in %s near Line %d.\n"
            "Please use \"UNSAT:XX\", \"GW:XX\", or \"RIVER:XX\", "
            "where XX is the index of grid.\n", filen, lno);
        pihm_exit(EXIT_FAILURE);
    }

    return loc;
}

int SpeciesType(FILE *fp, const char *chemn)
{
    /*
     * This subroutine is used to find out what the input species is.
     * 0) not found within database
     * 1) aqueous
     * 2) adsorption
     * 3) cation exchange
     * 4) mineral
     */
    int             return_val;
    char            tempn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             lno = 0;

    if (strcmp(chemn, "pH") == 0)
    {
        return AQUEOUS;
    }

    return_val = 0;

    sprintf(tempn, "'%s'", chemn);

    FindLine(fp, "BOF", &lno, ".cdbs");

    NextLine(fp, cmdstr, &lno);
    while (MatchWrappedKey(cmdstr, "'End of primary'") != 0)
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            return AQUEOUS;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (MatchWrappedKey(cmdstr, "'End of secondary'") != 0)
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            return 5;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (MatchWrappedKey(cmdstr, "'End of minerals'") != 0)
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            return MINERAL;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (strcmp(cmdstr, "End of surface complexation\r\n") != 0 &&
        strcmp(cmdstr, "End of surface complexation\n") != 0)
    {
        /* Notice that in CrunchFlow database, starting from surface
         * complexation, there is not apostrophe marks around blocking keywords
         */
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            return ADSORPTION;
        }
        NextLine(fp, cmdstr, &lno);
    }

    while (!feof(fp))
    {
        if (MatchWrappedKey(cmdstr, tempn) == 0)
        {
            return CATION_ECHG;
        }
        NextLine(fp, cmdstr, &lno);
    }

    return 0;
}

void SortChem(char chemn[MAXSPS][MAXSTRING], const int p_type[MAXSPS], int nsps,
    chemtbl_struct chemtbl[])
{
    int             i, j;
    int             temp;
    int             rank[MAXSPS];
    int             ranked_type[MAXSPS];

    for (i = 0; i < nsps; i++)
    {
        rank[i] = i;
        ranked_type[i] = p_type[i];
    }

    for (i = 0; i < nsps - 1; i++)
    {
        for (j = 0; j < nsps - i - 1; j++)
        {
            if (ranked_type[j] > ranked_type[j + 1])
            {
                temp = rank[j];
                rank[j] = rank[j + 1];
                rank[j + 1] = temp;

                temp = ranked_type[j];
                ranked_type[j] = ranked_type[j + 1];
                ranked_type[j + 1] = temp;
            }
        }
    }

    for (i = 0; i < nsps; i++)
    {
        strcpy(chemtbl[i].name, chemn[rank[i]]);
        chemtbl[i].itype = p_type[rank[i]];
    }
}
