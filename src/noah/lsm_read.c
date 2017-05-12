#include "pihm.h"

void ReadLsm (char *filename, double *latitude, double *longitude,
    ctrl_struct *ctrl, noahtbl_struct *noahtbl)
{
    int             i;
    FILE           *lsm_file;
    int             match;
    int             bytes_now;
    int             bytes_consumed = 0;
    char            cmdstr[MAXSTRING];
    char            buffer[MAXSTRING];
    int             lno = 0;

    /*
     * Open *.lsm file
     */
    lsm_file = fopen (filename, "r");
    CheckFile (lsm_file, filename);
    PIHMprintf (VL_VERBOSE, " Reading %s\n", filename);

    /*
     * Start reading lsm_file
     */
    FindLine (lsm_file, "BOF", &lno, filename);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "LATITUDE", latitude, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "LONGITUDE", longitude, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NSOIL", &ctrl->nsoil, 'i', filename, lno);
    if (ctrl->nsoil > MAXLYR - 1)
    {
        PIHMprintf (VL_ERROR,
            "The number of soil layers should not be larger than %d.\n",
            MAXLYR - 1);
        PIHMprintf (VL_ERROR,
            "Error in %s near Line %d.\n", filename, lno);
        PIHMexit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SLDPTH_DATA", buffer, 's', filename, lno);

    for (i = 0; i < ctrl->nsoil; i++)
    {
        match =
            sscanf (buffer + bytes_consumed, "%lf%n", &ctrl->sldpth[i],
            &bytes_now);
        if (match != 1)
        {
            PIHMprintf (VL_ERROR,
                "Error reading soil layer depths.\n");
            PIHMprintf (VL_ERROR,
                "Error in %s near Line %d.\n", filename, lno);
            PIHMexit (EXIT_FAILURE);
        }
        bytes_consumed += bytes_now;
    }

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RAD_MODE_DATA", &ctrl->rad_mode, 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SBETA_DATA", &noahtbl->sbeta, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "FXEXP_DATA", &noahtbl->fxexp, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "CSOIL_DATA", &noahtbl->csoil, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SALP_DATA", &noahtbl->salp, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "FRZK_DATA", &noahtbl->frzk, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ZBOT_DATA", &noahtbl->zbot, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "TBOT_DATA", &noahtbl->tbot, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "CZIL_DATA", &noahtbl->czil, 'd', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "LVCOEF_DATA", &noahtbl->lvcoef, 'd', filename, lno);

    /* Output control */
    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[T1_CTRL] = ReadPrtCtrl (cmdstr, "T1", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[STC_CTRL] = ReadPrtCtrl (cmdstr, "STC", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SMC_CTRL] = ReadPrtCtrl (cmdstr, "SMC", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SH2O_CTRL] = ReadPrtCtrl (cmdstr, "SH2O", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SNOWH_CTRL] = ReadPrtCtrl (cmdstr, "SNOWH", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ALBEDO_CTRL] = ReadPrtCtrl (cmdstr, "ALBEDO", filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[LE_CTRL] = ReadPrtCtrl (cmdstr, "LE", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SH_CTRL] = ReadPrtCtrl (cmdstr, "SH", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[G_CTRL] = ReadPrtCtrl (cmdstr, "G", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ETP_CTRL] = ReadPrtCtrl (cmdstr, "ETP", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ESNOW_CTRL] = ReadPrtCtrl (cmdstr, "ESNOW", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[ROOTW_CTRL] = ReadPrtCtrl (cmdstr, "ROOTW", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SOILM_CTRL] = ReadPrtCtrl (cmdstr, "SOILM", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[SOLAR_CTRL] = ReadPrtCtrl (cmdstr, "SOLAR", filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ctrl->prtvrbl[CH_CTRL] = ReadPrtCtrl (cmdstr, "CH", filename, lno);

    fclose (lsm_file);
}

void ReadRad (char *filename, forc_struct *forc)
{
    int             i, j;
    FILE           *rad_file;
    int             index;
    char            cmdstr[MAXSTRING];
    int             lno = 0;

    rad_file = fopen (filename, "r");
    CheckFile (rad_file, filename);
    PIHMprintf (VL_VERBOSE, " Reading %s\n", filename);

    FindLine (rad_file, "BOF", &lno, filename);

    forc->nrad = CountOccurance (rad_file, "RAD_TS");

    if (forc->nrad != forc->nmeteo)
    {
        PIHMprintf (VL_ERROR,
            "The number of radiation forcing time series should be the same "
            "as the number of meteorlogical forcing time series.\n");
        PIHMprintf (VL_ERROR, "Error in %s.\n", filename);
        PIHMexit (EXIT_FAILURE);
    }

    forc->rad = (tsdata_struct *)malloc (forc->nrad * sizeof (tsdata_struct));

    FindLine (rad_file, "BOF", &lno, filename);

    NextLine (rad_file, cmdstr, &lno);
    for (i = 0; i < forc->nrad; i++)
    {
        ReadKeyword (cmdstr, "RAD_TS", &index, 'i', filename, lno);

        if (i != index - 1)
        {
            PIHMprintf (VL_ERROR,
                "Error reading the %dth radiation forcing time series.\n",
                i + 1);
            PIHMprintf (VL_ERROR,
                "Error in %s near Line %d.\n", filename, lno);
            PIHMexit (EXIT_FAILURE);
        }

        /* Skip header lines */
        NextLine (rad_file, cmdstr, &lno);
        NextLine (rad_file, cmdstr, &lno);
        forc->rad[i].length = CountLine (rad_file, cmdstr, 1, "RAD_TS");
    }

    /* Rewind and read */
    FindLine (rad_file, "BOF", &lno, filename);
    for (i = 0; i < forc->nrad; i++)
    {
        /* Skip header lines */
        NextLine (rad_file, cmdstr, &lno);
        NextLine (rad_file, cmdstr, &lno);
        NextLine (rad_file, cmdstr, &lno);

        forc->rad[i].ftime =
            (int *)malloc (forc->rad[i].length * sizeof (int));
        forc->rad[i].data =
            (double **)malloc (forc->rad[i].length * sizeof (double *));
        for (j = 0; j < forc->rad[i].length; j++)
        {
            forc->rad[i].data[j] = (double *)malloc (2 * sizeof (double));
            NextLine (rad_file, cmdstr, &lno);
            ReadTS (cmdstr, &forc->rad[i].ftime[j],
                &forc->rad[i].data[j][0], 2);
        }
    }

    fclose (rad_file);
}
