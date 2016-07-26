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

    /*
     * Open *.lsm file
     */
    lsm_file = fopen (filename, "r");

    if (NULL == lsm_file)
    {
        fprintf (stderr, "Error opening %s.\n", filename);
        PIHMExit (EXIT_FAILURE);
    }

    if (verbose_mode)
    {
        printf ("Reading %s.\n", filename);
    }

    /*
     * Start reading lsm_file
     */
    FindLine (lsm_file, "BOF");
    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "LATITUDE", latitude, 'd'))
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "LONGITUDE", longitude, 'd'))
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "NSOIL", &ctrl->nsoil, 'i'))
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        PIHMExit (EXIT_FAILURE);
    }
    if (ctrl->nsoil > MAXLYR - 1)
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        fprintf (stderr,
            "The number of soil layers should not be larger than %d.\n",
            MAXLYR - 1);
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SLDPTH_DATA", buffer, 's'))
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        PIHMExit (EXIT_FAILURE);
    }

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

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "RAD_MODE_DATA", &ctrl->rad_mode, 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SBETA_DATA", &noahtbl->sbeta, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "FXEXP_DATA", &noahtbl->fxexp, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "CSOIL_DATA", &noahtbl->csoil, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SALP_DATA", &noahtbl->salp, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "FRZK_DATA", &noahtbl->frzk, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "ZBOT_DATA", &noahtbl->zbot, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "TBOT_DATA", &noahtbl->tbot, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "CZIL_DATA", &noahtbl->czil, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "LVCOEF_DATA", &noahtbl->lvcoef, 'd'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    /* Output control */
    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "T1", &ctrl->prtvrbl[T1_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "STC", &ctrl->prtvrbl[STC_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SMC", &ctrl->prtvrbl[SMC_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SH2O", &ctrl->prtvrbl[SH2O_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SNOWH", &ctrl->prtvrbl[SNOWH_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "ALBEDO", &ctrl->prtvrbl[ALBEDO_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "LE", &ctrl->prtvrbl[LE_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SH", &ctrl->prtvrbl[SH_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "G", &ctrl->prtvrbl[G_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "ETP", &ctrl->prtvrbl[ETP_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "ESNOW", &ctrl->prtvrbl[ESNOW_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "ROOTW", &ctrl->prtvrbl[ROOTW_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SOILM", &ctrl->prtvrbl[SOILM_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    NextLine (lsm_file, cmdstr);
    if (!ReadKeyword (cmdstr, "SOLAR", &ctrl->prtvrbl[SOLAR_CTRL], 'i'))
    {
        fprintf (stderr, "Pleareadheck SLDPTH_DATA.\n");
        PIHMExit (EXIT_FAILURE);
    }

    fclose (lsm_file);
}

void ReadRad (char *filename, forc_struct *forc)
{
    int             i, j;
    FILE           *rad_file;
    int             index;
    char            cmdstr[MAXSTRING];

    rad_file = fopen (filename, "r");

    if (NULL == rad_file)
    {
        fprintf (stderr, "Error opening %s.\n", filename);
        PIHMExit (EXIT_FAILURE);
    }

    if (verbose_mode)
    {
        printf ("Reading %s.\n", filename);
    }

    FindLine (rad_file, "BOF");

    forc->nrad = CountOccurance (rad_file, "RAD_TS");

    if (forc->nrad != forc->nmeteo)
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        fprintf (stderr,
            "The number of radiation forcing time series should be the same as the number of meteorlogical forcing time series.\n");
        PIHMExit (EXIT_FAILURE);
    }

    forc->rad = (tsdata_struct *)malloc (forc->nrad * sizeof (tsdata_struct));

    FindLine (rad_file, "BOF");

    NextLine (rad_file, cmdstr);
    for (i = 0; i < forc->nrad; i++)
    {
        if (!ReadKeyword (cmdstr, "RAD_TS", &index, 'i'))
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            PIHMExit (EXIT_FAILURE);
        }

        if (i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr,
                "Cannot read information of the %dth forcing series.\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
        }

        /* Skip header lines */
        NextLine (rad_file, cmdstr);
        NextLine (rad_file, cmdstr);
        forc->rad[i].length = CountLine (rad_file, cmdstr, 1, "RAD_TS");
    }

    /* Rewind and read */
    FindLine (rad_file, "BOF");
    for (i = 0; i < forc->nrad; i++)
    {
        /* Skip header lines */
        NextLine (rad_file, cmdstr);
        NextLine (rad_file, cmdstr);
        NextLine (rad_file, cmdstr);

        forc->rad[i].ftime =
            (int *)malloc (forc->rad[i].length * sizeof (int));
        forc->rad[i].data =
            (double **)malloc (forc->rad[i].length * sizeof (double *));
        for (j = 0; j < forc->rad[i].length; j++)
        {
            forc->rad[i].data[j] = (double *)malloc (2 * sizeof (double));
            NextLine (rad_file, cmdstr);
            ReadTS (cmdstr, &forc->rad[i].ftime[j],
                &forc->rad[i].data[j][0], 2);
        }
    }

    fclose (rad_file);
}
