#include "pihm.h"

void ReadLsm(const char fn[], ctrl_struct *ctrl, siteinfo_struct *siteinfo, noahtbl_struct *noahtbl)
{
    int             kz;
    FILE           *fp;
    int             bytes_now;
    int             bytes_consumed = 0;
    char            cmdstr[MAXSTRING];
    char            buffer[MAXSTRING];
    int             lno = 0;

    // Open *.lsm file
    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    // Start reading lsm_file
    FindLine(fp, "BOF", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "LATITUDE", 'd', fn, lno, &siteinfo->latitude);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "LONGITUDE", 'd', fn, lno, &siteinfo->longitude);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NSOIL", 'i', fn, lno, &ctrl->nlayers);
    if (ctrl->nlayers > MAXLYR - 1)
    {
        pihm_printf(VL_ERROR, "The number of soil layers should not be larger than %d.\n", MAXLYR - 1);
        pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
        pihm_exit(EXIT_FAILURE);
    }

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SLDPTH_DATA", 's', fn, lno, buffer);

    for (kz = 0; kz < ctrl->nlayers; kz++)
    {
        if (sscanf(buffer + bytes_consumed, "%lf%n", &ctrl->soil_depth[kz], &bytes_now) != 1)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
        bytes_consumed += bytes_now;
    }

#if defined(_CYCLES_)
    for (kz = ctrl->nlayers; kz < MAXLYR; kz++)
    {
        ctrl->soil_depth[kz] = ctrl->soil_depth[kz - 1];
    }
#endif

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RAD_MODE_DATA", 'i', fn, lno, &ctrl->rad_mode);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SBETA_DATA", 'd', fn, lno, &noahtbl->sbeta);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "FXEXP_DATA", 'd', fn, lno, &noahtbl->fxexp);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CSOIL_DATA", 'd', fn, lno, &noahtbl->csoil);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SALP_DATA", 'd', fn, lno, &noahtbl->salp);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "FRZK_DATA", 'd', fn, lno, &noahtbl->frzk);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ZBOT_DATA", 'd', fn, lno, &noahtbl->zbot);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TBOT_DATA", 'd', fn, lno, &noahtbl->tbot);
    siteinfo->tavg = noahtbl->tbot;

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CZIL_DATA", 'd', fn, lno, &noahtbl->czil);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "LVCOEF_DATA", 'd', fn, lno, &noahtbl->lvcoef);

    // Output control
    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[T1_CTRL] = ReadPrintCtrl(cmdstr, "T1", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[STC_CTRL] = ReadPrintCtrl(cmdstr, "STC", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SMC_CTRL] = ReadPrintCtrl(cmdstr, "SMC", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SH2O_CTRL] = ReadPrintCtrl(cmdstr, "SH2O", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SNOWH_CTRL] = ReadPrintCtrl(cmdstr, "SNOWH", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[ALBEDO_CTRL] = ReadPrintCtrl(cmdstr, "ALBEDO", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[LE_CTRL] = ReadPrintCtrl(cmdstr, "LE", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SH_CTRL] = ReadPrintCtrl(cmdstr, "SH", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[G_CTRL] = ReadPrintCtrl(cmdstr, "G", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[ETP_CTRL] = ReadPrintCtrl(cmdstr, "ETP", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[ESNOW_CTRL] = ReadPrintCtrl(cmdstr, "ESNOW", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[ROOTW_CTRL] = ReadPrintCtrl(cmdstr, "ROOTW", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SOILM_CTRL] = ReadPrintCtrl(cmdstr, "SOILM", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SOLAR_CTRL] = ReadPrintCtrl(cmdstr, "SOLAR", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[CH_CTRL] = ReadPrintCtrl(cmdstr, "CH", fn, lno);

    fclose(fp);
}

void ReadRad(const char fn[], forc_struct *forc)
{
    int             i, j;
    FILE           *fp;
    int             index;
    char            cmdstr[MAXSTRING];
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    FindLine(fp, "BOF", &lno, fn);

    forc->nrad = CountOccurr(fp, "RAD_TS");

    if (forc->nrad != forc->nmeteo)
    {
        pihm_printf(VL_ERROR, "The number of radiation forcing time series should be the same as the number of\n"
            "meteorological forcing time series.\nError in %s.\n", fn);
        pihm_exit(EXIT_FAILURE);
    }

    forc->rad = (tsdata_struct *)malloc(forc->nrad * sizeof(tsdata_struct));

    FindLine(fp, "BOF", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    for (i = 0; i < forc->nrad; i++)
    {
        ReadKeyword(cmdstr, "RAD_TS", 'i', fn, lno, &index);

        if (i != index - 1)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }

        // Check header lines
        NextLine(fp, cmdstr, &lno);
        if (!CheckHeader(cmdstr, 3, "TIME", "SDIR", "SDIF"))
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
        forc->rad[i].length = CountLine(fp, cmdstr, 1, "RAD_TS");
    }

    // Rewind and read
    FindLine(fp, "BOF", &lno, fn);
    for (i = 0; i < forc->nrad; i++)
    {
        // Skip header lines
        NextLine(fp, cmdstr, &lno);
        NextLine(fp, cmdstr, &lno);

        forc->rad[i].ftime = (int *)malloc(forc->rad[i].length * sizeof(int));
        forc->rad[i].data = (double **)malloc(forc->rad[i].length * sizeof(double *));
        for (j = 0; j < forc->rad[i].length; j++)
        {
            forc->rad[i].data[j] = (double *)malloc(2 * sizeof(double));
            NextLine(fp, cmdstr, &lno);
            if (!ReadTs(cmdstr, 2, &forc->rad[i].ftime[j], &forc->rad[i].data[j][0]))
            {
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }
        }
    }

    fclose(fp);
}

void ReadGlacierIce(const char fn[], double iceh[])
{
    int             i;
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    NextLine(fp, cmdstr, &lno);
    for (i = 0; i < nelem; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf", &index, &iceh[i]);
        if (match != 2)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    fclose(fp);
}
