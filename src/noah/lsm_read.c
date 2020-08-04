#include "pihm.h"

void ReadLsm(const char filename[], ctrl_struct *ctrl,
    siteinfo_struct *siteinfo, noahtbl_struct *noahtbl)
{
    int             kz;
    FILE           *lsm_file;
    int             bytes_now;
    int             bytes_consumed = 0;
    char            cmdstr[MAXSTRING];
    char            buffer[MAXSTRING];
    int             lno = 0;

    /*
     * Open *.lsm file
     */
    lsm_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    /*
     * Start reading lsm_file
     */
    FindLine(lsm_file, "BOF", &lno, filename);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "LATITUDE", 'd', filename, lno, &siteinfo->latitude);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "LONGITUDE", 'd', filename, lno, &siteinfo->longitude);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NSOIL", 'i', filename, lno, &ctrl->nlayers);
    if (ctrl->nlayers > MAXLYR - 1)
    {
        pihm_printf(VL_ERROR,
            "The number of soil layers should not be larger than %d.\n",
            MAXLYR - 1);
        pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
        pihm_exit(EXIT_FAILURE);
    }

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SLDPTH_DATA", 's', filename, lno, buffer);

    for (kz = 0; kz < ctrl->nlayers; kz++)
    {
        if (sscanf(buffer + bytes_consumed, "%lf%n", &ctrl->soil_depth[kz],
            &bytes_now) != 1)
        {
            pihm_printf(VL_ERROR, "Error reading soil layer depths.\n"
                "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
        bytes_consumed += bytes_now;
    }

#if defined(_CYCLES_)
    for (kz = ctrl->nlayers; kz < MAXLYR; kz++)
    {
        ctrl->soil_depth[kz] = ctrl->soil_depth[kz - 1];
    }
#endif

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "RAD_MODE_DATA", 'i', filename, lno, &ctrl->rad_mode);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SBETA_DATA", 'd', filename, lno, &noahtbl->sbeta);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "FXEXP_DATA", 'd', filename, lno, &noahtbl->fxexp);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "CSOIL_DATA", 'd', filename, lno, &noahtbl->csoil);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SALP_DATA", 'd', filename, lno, &noahtbl->salp);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "FRZK_DATA", 'd', filename, lno, &noahtbl->frzk);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "ZBOT_DATA", 'd', filename, lno, &noahtbl->zbot);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "TBOT_DATA", 'd', filename, lno, &noahtbl->tbot);
    siteinfo->tavg = noahtbl->tbot;

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "CZIL_DATA", 'd', filename, lno, &noahtbl->czil);

    NextLine(lsm_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "LVCOEF_DATA", 'd', filename, lno, &noahtbl->lvcoef);

    /* Output control */
    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[T1_CTRL] = ReadPrintCtrl(cmdstr, "T1", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[STC_CTRL] = ReadPrintCtrl(cmdstr, "STC", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SMC_CTRL] = ReadPrintCtrl(cmdstr, "SMC", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SH2O_CTRL] = ReadPrintCtrl(cmdstr, "SH2O", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SNOWH_CTRL] = ReadPrintCtrl(cmdstr, "SNOWH", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ALBEDO_CTRL] = ReadPrintCtrl(cmdstr, "ALBEDO", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[LE_CTRL] = ReadPrintCtrl(cmdstr, "LE", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SH_CTRL] = ReadPrintCtrl(cmdstr, "SH", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[G_CTRL] = ReadPrintCtrl(cmdstr, "G", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ETP_CTRL] = ReadPrintCtrl(cmdstr, "ETP", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ESNOW_CTRL] = ReadPrintCtrl(cmdstr, "ESNOW", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ROOTW_CTRL] = ReadPrintCtrl(cmdstr, "ROOTW", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SOILM_CTRL] = ReadPrintCtrl(cmdstr, "SOILM", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SOLAR_CTRL] = ReadPrintCtrl(cmdstr, "SOLAR", filename, lno);

    NextLine(lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[CH_CTRL] = ReadPrintCtrl(cmdstr, "CH", filename, lno);

    fclose(lsm_file);
}

void ReadRad(const char filename[], forc_struct *forc)
{
    int             i, j;
    FILE           *rad_file;
    int             index;
    char            cmdstr[MAXSTRING];
    int             lno = 0;

    rad_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    FindLine(rad_file, "BOF", &lno, filename);

    forc->nrad = CountOccurr(rad_file, "RAD_TS");

    if (forc->nrad != forc->nmeteo)
    {
        pihm_printf(VL_ERROR, "The number of radiation forcing time series "
            "should be the same as the number of\nmeteorological forcing time "
            "series.\nError in %s.\n", filename);
        pihm_exit(EXIT_FAILURE);
    }

    forc->rad = (tsdata_struct *)malloc(forc->nrad * sizeof(tsdata_struct));

    FindLine(rad_file, "BOF", &lno, filename);

    NextLine(rad_file, cmdstr, &lno);
    for (i = 0; i < forc->nrad; i++)
    {
        ReadKeyword(cmdstr, "RAD_TS", 'i', filename, lno, &index);

        if (i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading the %dth radiation forcing time series.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }

        /* Skip header lines */
        NextLine(rad_file, cmdstr, &lno);
        NextLine(rad_file, cmdstr, &lno);
        forc->rad[i].length = CountLine(rad_file, cmdstr, 1, "RAD_TS");
    }

    /* Rewind and read */
    FindLine(rad_file, "BOF", &lno, filename);
    for (i = 0; i < forc->nrad; i++)
    {
        /* Skip header lines */
        NextLine(rad_file, cmdstr, &lno);
        NextLine(rad_file, cmdstr, &lno);
        NextLine(rad_file, cmdstr, &lno);

        forc->rad[i].ftime = (int *)malloc(forc->rad[i].length * sizeof(int));
        forc->rad[i].data =
            (double **)malloc(forc->rad[i].length * sizeof(double *));
        for (j = 0; j < forc->rad[i].length; j++)
        {
            forc->rad[i].data[j] = (double *)malloc(2 * sizeof(double));
            NextLine(rad_file, cmdstr, &lno);
            ReadTs(cmdstr, 2, &forc->rad[i].ftime[j], &forc->rad[i].data[j][0]);
        }
    }

    fclose(rad_file);
}

void ReadGlacierIce(const char filename[], double iceh[])
{
    int             i;
    FILE           *ice_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    ice_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    NextLine(ice_file, cmdstr, &lno);
    for (i = 0; i < nelem; i++)
    {
        NextLine(ice_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf", &index, &iceh[i]);
        if (match != 2)
        {
            pihm_printf(VL_ERROR,
                "Error reading glacier ice depth of the %dth element.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    fclose(ice_file);
}
