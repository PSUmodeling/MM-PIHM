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
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

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
        fprintf (stderr, "Error reading %s.\n", filename);
        fprintf (stderr,
            "The number of soil layers should not be larger than %d.\n",
            MAXLYR - 1);
        PIHMExit (EXIT_FAILURE);
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
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr, "Please check SLDPTH_DATA.\n");
            PIHMExit (EXIT_FAILURE);
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
    ReadKeyword (cmdstr, "T1", &ctrl->prtvrbl[T1_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "STC", &ctrl->prtvrbl[STC_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SMC", &ctrl->prtvrbl[SMC_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SH2O", &ctrl->prtvrbl[SH2O_CTRL], 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SNOWH", &ctrl->prtvrbl[SNOWH_CTRL], 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ALBEDO", &ctrl->prtvrbl[ALBEDO_CTRL], 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "LE", &ctrl->prtvrbl[LE_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SH", &ctrl->prtvrbl[SH_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "G", &ctrl->prtvrbl[G_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ETP", &ctrl->prtvrbl[ETP_CTRL], 'i', filename, lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ESNOW", &ctrl->prtvrbl[ESNOW_CTRL], 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ROOTW", &ctrl->prtvrbl[ROOTW_CTRL], 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SOILM", &ctrl->prtvrbl[SOILM_CTRL], 'i', filename,
        lno);

    NextLine (lsm_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SOLAR", &ctrl->prtvrbl[SOLAR_CTRL], 'i', filename,
        lno);

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
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    FindLine (rad_file, "BOF", &lno, filename);

    forc->nrad = CountOccurance (rad_file, "RAD_TS");

    if (forc->nrad != forc->nmeteo)
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        fprintf (stderr,
            "The number of radiation forcing time series should be the same as the number of meteorlogical forcing time series.\n");
        PIHMExit (EXIT_FAILURE);
    }

    forc->rad = (tsdata_struct *)malloc (forc->nrad * sizeof (tsdata_struct));

    FindLine (rad_file, "BOF", &lno, filename);

    NextLine (rad_file, cmdstr, &lno);
    for (i = 0; i < forc->nrad; i++)
    {
        ReadKeyword (cmdstr, "RAD_TS", &index, 'i', filename, lno);

        if (i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr,
                "Cannot read information of the %dth forcing series.\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
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
