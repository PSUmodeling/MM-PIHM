/*****************************************************************************
 * File        : read_alloc.c
 * Function    : read parameter files
 *----------------------------------------------------------------------------
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0...........................
 * a) Addition of three new input files: file.calib, file.lc and file.geol
 * b) Declaration and allocation  of new variables for new process, shape
 *    representations  and calibration (e.g. ET, Infiltration, Macropore,
 *    Stormflow, Element beneath river, river shapes, river bed property,
 *    thresholds for root zone, infiltration and macropore depths, land cover
 *    attributes etc)                                                      
 ****************************************************************************/

#include "pihm.h"

void read_alloc (char *project, Model_Data DS, Control_Data  CS)
{
    char           *simulation;
    char           *token, *tempname;
    struct tm      *timeinfo;

    tempname = (char *)malloc ((strlen (project) + 1) * sizeof (char));
    strcpy (tempname, project);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        simulation = (char *)malloc ((strlen (token) + 1) * sizeof (char));
        strcpy (simulation, token);
    }
    else
    {
        simulation = (char *)malloc ((strlen (project) + 1)
           * sizeof (char));
        strcpy (simulation, project);
    }
    free (tempname);

    if (CS->Verbose)
        printf ("\nStart reading in input files:\n");

    ReadRiv (project, DS, CS);

    ReadMesh (project, DS, CS);

    ReadAtt (project, DS, CS);

    ReadSoil (project, DS, CS);

    ReadGeol (project, DS, CS);

    ReadLC (project, DS, CS);

    ReadForc (project, DS, CS);

    ReadIbc (project, DS, CS);

    ReadPara (project, DS, CS);

    ReadCalib (simulation, DS, CS);
}

void ReadRiv (char *simulation, Model_Data DS, Control_Data CS)
{
    int             i, j;
    char           *fn;
    FILE           *riv_file;   /* Pointer to .riv file */
    struct tm      *timeinfo;
    time_t          rawtime;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "riv");
    fn = (char *)malloc ((2 * strlen (simulation) + 12) * sizeof (char));
    sprintf (fn, "input/%s/%s.riv", simulation, simulation);
    riv_file = fopen (fn, "r");
    free (fn);

    if (riv_file == NULL)
    {
        printf ("\n Fatal Error: %s.riv is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading riv_file */
    fscanf (riv_file, "%d %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s", &DS->NumRiv);
    DS->Riv = (river_segment *) malloc (DS->NumRiv * sizeof (river_segment));

    for (i = 0; i < DS->NumRiv; i++)
    {
        fscanf (riv_file, "%d", &(DS->Riv[i].index));
        fscanf (riv_file, "%d %d", &(DS->Riv[i].FromNode), &(DS->Riv[i].ToNode));
        fscanf (riv_file, "%d", &(DS->Riv[i].down));
        fscanf (riv_file, "%d %d", &(DS->Riv[i].LeftEle), &(DS->Riv[i].RightEle));
        fscanf (riv_file, "%d %d", &(DS->Riv[i].shape), &(DS->Riv[i].material));
        fscanf (riv_file, "%d %d", &(DS->Riv[i].IC), &(DS->Riv[i].BC));
        fscanf (riv_file, "%d", &(DS->Riv[i].reservoir));
    }

    fscanf (riv_file, "%*s");
    fscanf (riv_file, "%d %*s %*s %*s", &DS->NumRivShape);
    DS->Riv_Shape = (river_shape *) malloc (DS->NumRivShape * sizeof (river_shape));

    for (i = 0; i < DS->NumRivShape; i++)
    {
        fscanf (riv_file, "%d", &DS->Riv_Shape[i].index);
        fscanf (riv_file, "%lf", &DS->Riv_Shape[i].depth);
        fscanf (riv_file, "%d %lf", &DS->Riv_Shape[i].interpOrd, &DS->Riv_Shape[i].coeff);
    }

    fscanf (riv_file, "%*s");
    fscanf (riv_file, "%d %*s %*s %*s %*s %*s", &DS->NumRivMaterial);
    DS->Riv_Mat = (river_material *) malloc (DS->NumRivMaterial * sizeof (river_material));

    for (i = 0; i < DS->NumRivMaterial; i++)
        fscanf (riv_file, "%d %lf %lf %lf %lf %lf", &DS->Riv_Mat[i].index, &DS->Riv_Mat[i].Rough, &DS->Riv_Mat[i].Cwr, &DS->Riv_Mat[i].KsatH, &DS->Riv_Mat[i].KsatV, &DS->Riv_Mat[i].bedThick);

    fscanf (riv_file, "%*s");
    fscanf (riv_file, "%d %*s", &DS->NumRivIC);
    DS->Riv_IC = (river_IC *) malloc (DS->NumRivIC * sizeof (river_IC));

    for (i = 0; i < DS->NumRivIC; i++)
        fscanf (riv_file, "%d %lf", &DS->Riv_IC[i].index, &DS->Riv_IC[i].value);

    fscanf (riv_file, "%*s");
    fscanf (riv_file, "%d", &DS->NumRivBC);
    DS->TSD_Riv = (TSD *) malloc (DS->NumRivBC * sizeof (TSD));

    for (i = 0; i < DS->NumRivBC; i++)
    {
        fscanf (riv_file, "%s %d %d", DS->TSD_Riv[i].name, &DS->TSD_Riv[i].index, &DS->TSD_Riv[i].length);

        DS->TSD_Riv[i].TS = (realtype **) malloc ((DS->TSD_Riv[i].length) * sizeof (realtype *));
        for (j = 0; j < DS->TSD_Riv[i].length; j++)
            DS->TSD_Riv[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));

        for (j = 0; j < DS->TSD_Riv[i].length; j++)
        {
            fscanf (riv_file, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_Riv[i].TS[j][1]);
            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            rawtime = timegm (timeinfo);
            DS->TSD_Riv[i].TS[j][0] = (realtype) rawtime;
        }
    }

    /* read in reservoir information */
    fscanf (riv_file, "%*s");
    fscanf (riv_file, "%d", &DS->NumRes);
    if (DS->NumRes > 0)
    {
        /* read in reservoir information */
    }

    fclose (riv_file);
}

void ReadMesh (char *simulation, Model_Data DS, Control_Data CS)
{
    FILE           *mesh_file;  /* Pointer to .mesh file */
    char           *fn;
    int             i;

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "mesh");
    fn = (char *)malloc ((2 * strlen (simulation) + 13) * sizeof (char));
    sprintf (fn, "input/%s/%s.mesh", simulation, simulation);
    mesh_file = fopen (fn, "r");
    free (fn);

    if (mesh_file == NULL)
    {
        printf ("\n  Fatal Error: %s.mesh is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading mesh_file */
    fscanf (mesh_file, "%d %*s %*s %*s %*s %*s %*s", &DS->NumEle);
    DS->Ele = (element *) malloc ((DS->NumEle + DS->NumRiv) * sizeof (element));

    /* read in elements information */
    for (i = 0; i < DS->NumEle; i++)
    {
        fscanf (mesh_file, "%d", &(DS->Ele[i].index));
        fscanf (mesh_file, "%d %d %d", &(DS->Ele[i].node[0]), &(DS->Ele[i].node[1]), &(DS->Ele[i].node[2]));
        fscanf (mesh_file, "%d %d %d", &(DS->Ele[i].nabr[0]), &(DS->Ele[i].nabr[1]), &(DS->Ele[i].nabr[2]));
    }

    /* read in nodes information */
    fscanf (mesh_file, "%d %*s %*s %*s %*s", &DS->NumNode);
    DS->Node = (nodes *) malloc (DS->NumNode * sizeof (nodes));
    for (i = 0; i < DS->NumNode; i++)
    {
        fscanf (mesh_file, "%d", &(DS->Node[i].index));
        fscanf (mesh_file, "%lf %lf", &(DS->Node[i].x), &(DS->Node[i].y));
        fscanf (mesh_file, "%lf %lf", &(DS->Node[i].zmin), &(DS->Node[i].zmax));
    }


    /* finish reading mesh_files */
    fclose (mesh_file);
}

void ReadAtt (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    int             i, j;
    FILE           *att_file;   /* Pointer to .att file */

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "att");
    fn = (char *)malloc ((2 * strlen (simulation) + 12) * sizeof (char));
    sprintf (fn, "input/%s/%s.att", simulation, simulation);
    att_file = fopen (fn, "r");
    free (fn);

    if (att_file == NULL)
    {
        printf ("\n  Fatal Error: %s.att is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading att_file */
    DS->Ele_IC = (element_IC *) malloc (DS->NumEle * sizeof (element_IC));
    for (i = 0; i < 16; i++)
        fscanf (att_file, "%*s");
    for (i = 0; i < DS->NumEle; i++)
    {
        fscanf (att_file, "%*d");
        fscanf (att_file, "%d %d %d", &(DS->Ele[i].soil), &(DS->Ele[i].geol), &(DS->Ele[i].LC));

        fscanf (att_file, "%lf %lf %lf %lf %lf", &(DS->Ele_IC[i].interception), &(DS->Ele_IC[i].snow), &(DS->Ele_IC[i].surf), &(DS->Ele_IC[i].unsat), &(DS->Ele_IC[i].sat));
        fscanf (att_file, "%d %d", &(DS->Ele[i].meteo), &(DS->Ele[i].LAI));
        fscanf (att_file, "%d", &(DS->Ele[i].source));
        for (j = 0; j < 3; j++)
            fscanf (att_file, "%d", &(DS->Ele[i].BC[j]));
        fscanf (att_file, "%d", &(DS->Ele[i].Macropore));
    }


    /* finish reading att_files */
    fclose (att_file);
}

void ReadSoil (char *simulation, Model_Data DS, Control_Data CS)
{
    FILE           *soil_file;  /* Pointer to .soil file */
    char           *fn;
    int             i;

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "soil");
    fn = (char *)malloc ((2 * strlen (simulation) + 13) * sizeof (char));
    sprintf (fn, "input/%s/%s.soil", simulation, simulation);
    soil_file = fopen (fn, "r");
    free (fn);

    if (soil_file == NULL)
    {
        printf ("\n  Fatal Error: %s.soil is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading soil_file */
    fscanf (soil_file, "%d %*s %*s %*s %*s %*s %*s %*s %*s %*s", &DS->NumSoil);
    DS->Soil = (soils *) malloc (DS->NumSoil * sizeof (soils));

    for (i = 0; i < DS->NumSoil; i++)
    {
        fscanf (soil_file, "%d", &(DS->Soil[i].index));
        /* Note: Soil KsatH and macKsatH is not used in model calculation
         * anywhere */
        fscanf (soil_file, "%lf", &(DS->Soil[i].KsatV));
        fscanf (soil_file, "%lf %lf %lf", &(DS->Soil[i].ThetaS), &(DS->Soil[i].ThetaR), &(DS->Soil[i].infD));
        fscanf (soil_file, "%lf %lf", &(DS->Soil[i].Alpha), &(DS->Soil[i].Beta));
        fscanf (soil_file, "%lf %lf", &(DS->Soil[i].hAreaF), &(DS->Soil[i].macKsatV));
        fscanf (soil_file, "%lf", &(DS->Soil[i].qtz));
    }

    fclose (soil_file);
}

void ReadGeol (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    FILE           *geol_file;  /* Pointer to .geol file */
    int             i;

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "geol");
    fn = (char *)malloc ((2 * strlen (simulation) + 13) * sizeof (char));
    sprintf (fn, "input/%s/%s.geol", simulation, simulation);
    geol_file = fopen (fn, "r");
    free (fn);

    if (geol_file == NULL)
    {
        printf ("\n  Fatal Error: %s.geol is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading geol_file */
    fscanf (geol_file, "%d %*s %*s %*s %*s %*s %*s %*s %*s %*s", &DS->NumGeol);
    DS->Geol = (geol *) malloc (DS->NumGeol * sizeof (geol));

    for (i = 0; i < DS->NumGeol; i++)
    {
        fscanf (geol_file, "%d", &(DS->Geol[i].index));
        /* Geol macKsatV is not used in model calculation anywhere */
        fscanf (geol_file, "%lf %lf", &(DS->Geol[i].KsatH), &(DS->Geol[i].KsatV));
        fscanf (geol_file, "%lf %lf", &(DS->Geol[i].ThetaS), &(DS->Geol[i].ThetaR));
        fscanf (geol_file, "%lf %lf", &(DS->Geol[i].Alpha), &(DS->Geol[i].Beta));
        fscanf (geol_file, "%lf %lf %lf", &(DS->Geol[i].vAreaF), &(DS->Geol[i].macKsatH), &(DS->Geol[i].macD));
    }

    fclose (geol_file);


}

void ReadLC (char *simulation, Model_Data DS, Control_Data CS)
{
    FILE           *lc_file;    /* Pointer to .lc file */
    int             i;

    if (CS->Verbose)
        printf ("  Reading vegprmt.tbl\n");
    lc_file = fopen ("input/vegprmt.tbl", "r");

    if (lc_file == NULL)
    {
        printf ("\n  Fatal Error: land cover file input/vegprmt.tbl is in use or does not exist!\n");
        exit (1);
    }

    /* start reading land cover file */
    fscanf (lc_file, "%d %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s", &DS->NumLC);

    DS->LandC = (LC *) malloc (DS->NumLC * sizeof (LC));

    for (i = 0; i < DS->NumLC; i++)
    {
        fscanf (lc_file, "%d", &(DS->LandC[i].index));
        fscanf (lc_file, "%lf", &(DS->LandC[i].VegFrac));
        fscanf (lc_file, "%lf", &(DS->LandC[i].RzD));
        fscanf (lc_file, "%lf", &(DS->LandC[i].Rmin));
        fscanf (lc_file, "%lf", &(DS->LandC[i].Rs_ref));
        fscanf (lc_file, "%lf", &(DS->LandC[i].h_s));
        fscanf (lc_file, "%lf", &(DS->LandC[i].snup));
        fscanf (lc_file, "%lf %lf", &(DS->LandC[i].LAImin), &(DS->LandC[i].LAImax));
        fscanf (lc_file, "%lf %lf", &(DS->LandC[i].Emiss_min), &(DS->LandC[i].Emiss_max));
        fscanf (lc_file, "%lf %lf", &(DS->LandC[i].Albedo_min), &(DS->LandC[i].Albedo_max));
        fscanf (lc_file, "%lf %lf", &(DS->LandC[i].z0_min), &(DS->LandC[i].z0_max));
        fscanf (lc_file, "%lf%*[^\n]", &(DS->LandC[i].Rough));
    }

    fscanf (lc_file, "%*s %lf", &(DS->Tref));
    fscanf (lc_file, "%*s %lf", &(DS->fx_canopy));
    fscanf (lc_file, "%*s %lf", &(DS->Rmax));
    fscanf (lc_file, "%*s %d", &(DS->bare));

    DS->ISFactor = (realtype *) malloc (DS->NumLC * sizeof (realtype));
    for ( i = 0; i < DS->NumLC; i++)
        DS->ISFactor[i] = 0.0002;
    fclose (lc_file);


}

void ReadForc (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    FILE           *forc_file;  /* Pointer to .forc file */
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             ind;
    int            *count;
    int             i, j;
    struct tm      *timeinfo;
    time_t          rawtime;
    int             read_lai;
    int             read_ss;
    FILE           *lai_file;
    int             NumForcing;
    char           *laifn;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "forc");
    fn = (char *)malloc ((2 * strlen (simulation) + 13) * sizeof (char));
    sprintf (fn, "input/%s/%s.forc", simulation, simulation);
    forc_file = fopen (fn, "r");
    free (fn);

    if (forc_file == NULL)
    {
        printf ("\n  Fatal Error: %s.forc is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading forc_file */

    NumForcing = 7;

    /*
     * Forcing TS:
     * 0: Precipitation;
     * 1: Surface temperature;
     * 2: Relative humidity;
     * 3: Surface wind speed;
     * 4: Downward solar radiation;
     * 5: Downward longwave radiation (or G in PIHM);
     * 6: Surface air pressure;
     */

    fscanf (forc_file, "%*s %d", &DS->NumTS);

    DS->TSD_meteo = (TSD *) malloc (DS->NumTS * sizeof (TSD));

    rewind(forc_file);          /* For safety reasons, rewind and skip top line */
    fgets (cmdstr, MAXSTRING, forc_file);

    fgets (cmdstr, MAXSTRING, forc_file);
    while (!feof (forc_file))
    {
        if (cmdstr[0] != '\n' && cmdstr[0] != '\0' && cmdstr[0] != '\t')
        {
            sscanf (cmdstr, "%s", optstr);
            if (strcasecmp ("METEO_TS", optstr) == 0)
            {
                sscanf (cmdstr, "%*s %d %*s %*f", &ind);
                DS->TSD_meteo[ind - 1].length = 0;
                count = &(DS->TSD_meteo[ind - 1].length);
            }
            else if (strcasecmp ("TIME", optstr) == 0)
            {
                /* Do nothing */
            }
            else if (strcasecmp ("TS", optstr) == 0)
            {
                /* Do nothing */
            }
            else
            {
                (*count)++;
            }
        }
        fgets (cmdstr, MAXSTRING, forc_file);
    }

    if (ind != DS->NumTS)
    {
        printf ("ERROR!\n");
        exit (1);
    }

    for (i = 0; i < DS->NumTS; i++)
    {
        DS->TSD_meteo[i].TS = (realtype **) malloc ((DS->TSD_meteo[i].length) * sizeof (realtype *));
        DS->TSD_meteo[i].iCounter = 0;
        for (j = 0; j < DS->TSD_meteo[i].length; j++)
            DS->TSD_meteo[i].TS[j] = (realtype *) malloc ((NumForcing + 1)* sizeof (realtype));
    }

    rewind(forc_file);
    fscanf (forc_file, "%*s %*d");

    for (i = 0; i < DS->NumTS; i++)
    {
        fscanf (forc_file, "%*s %*d %*s %lf", &(DS->TSD_meteo[i].TSFactor));
        fscanf (forc_file, "%*s %*s %*s %*s %*s %*s %*s %*s");
        fscanf (forc_file, "%*s %*s %*s %*s %*s %*s %*s %*s");

        for (j = 0; j < DS->TSD_meteo[i].length; j++)
        {
            fscanf (forc_file, "%d-%d-%d %d:%d %lf %lf %lf %lf %lf %lf %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_meteo[i].TS[j][1], &DS->TSD_meteo[i].TS[j][2], &DS->TSD_meteo[i].TS[j][3], &DS->TSD_meteo[i].TS[j][4], &DS->TSD_meteo[i].TS[j][5], &DS->TSD_meteo[i].TS[j][6], &DS->TSD_meteo[i].TS[j][7]);
            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            timeinfo->tm_sec = 0;
            rawtime = timegm (timeinfo);
            DS->TSD_meteo[i].TS[j][0] = (realtype) rawtime;
        }
    }

    fclose (forc_file);

    DS->windH = (realtype *) malloc (DS->NumTS * sizeof (realtype));
    for (i = 0; i < DS->NumTS; i++)
        DS->windH[i] = DS->TSD_meteo[i].TSFactor;

    read_lai = 0;
    read_ss = 0;

    for (i = 0; i < DS->NumEle; i++)
    {
        if (DS->Ele[i].LAI > 0)
            read_lai = 1;
        if (DS->Ele[i].source >0)
            read_ss = 1;
    }

    DS->NumLAI = 0;
    if (read_lai == 1)
    {
        if (CS->Verbose)
            printf ("  Reading %s.%s\n", simulation, "lai");
        laifn = (char *)malloc ((2 * strlen (simulation) + 12) * sizeof (char));
        sprintf (laifn, "input/%s/%s.lai", simulation, simulation);
        lai_file = fopen (laifn, "r");
        free (laifn);

        if (lai_file == NULL)
        {
            printf ("\n  Fatal Error: %s.lai is in use or does not exist!\n",
               simulation);
            exit (1);
        }

        /* start reading lai_file */
        fscanf (lai_file, "%*s %d", &(DS->NumLAI));

        DS->TSD_lai = (TSD *) malloc (DS->NumLAI * sizeof (TSD));

        rewind(lai_file);          /* For safety reasons, rewind and skip top line */
        fgets (cmdstr, MAXSTRING, lai_file);

        fgets (cmdstr, MAXSTRING, lai_file);
        while (!feof (lai_file))
        {
            if (cmdstr[0] != '\n' && cmdstr[0] != '\0' && cmdstr[0] != '\t')
            {
                sscanf (cmdstr, "%s", optstr);
                if (strcasecmp ("LAI_TS", optstr) == 0)
                {
                    sscanf (cmdstr, "%*s %d", &ind);
                    DS->TSD_lai[ind - 1].length = 0;
                    count = &(DS->TSD_lai[ind - 1].length);
                }
                else if (strcasecmp ("TIME", optstr) == 0)
                {
                    /* Do nothing */
                }
                else if (strcasecmp ("TS", optstr) == 0)
                {
                    /* Do nothing */
                }
                else
                    (*count)++;
            }
            fgets (cmdstr, MAXSTRING, lai_file);
        }

        if (ind != DS->NumLAI)
        {
            printf ("ERROR!\n");
            exit (1);
        }

        for (i = 0; i < DS->NumLAI; i++)
        {
            DS->TSD_lai[i].TS = (realtype **) malloc ((DS->TSD_lai[i].length) * sizeof (realtype *));
            for (j = 0; j < DS->TSD_lai[i].length; j++)
                DS->TSD_lai[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));
        }

        rewind(lai_file);

        fscanf (lai_file, "%*s %*d");

        for (i = 0; i < DS->NumLAI; i++)
        {
            fscanf (lai_file, "%*s %*d");
            fscanf (lai_file, "%*s %*s");
            fscanf (lai_file, "%*s %*s");
            for (j = 0; j < DS->TSD_lai[i].length; j++)
            {
                fscanf (lai_file, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_lai[i].TS[j][1]);
                timeinfo->tm_year = timeinfo->tm_year - 1900;
                timeinfo->tm_mon = timeinfo->tm_mon - 1;
                timeinfo->tm_sec = 0;
                rawtime = timegm (timeinfo);
                DS->TSD_lai[i].TS[j][0] = (realtype) rawtime;
            }
        }
        fclose (lai_file);
    }

    if (read_ss == 1)
    {
        /* Read .ss */
    }

}

void ReadIbc (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    int             i, j;
    time_t          rawtime;
    struct tm      *timeinfo;
    FILE           *ibc_file;   /* Pointer to .ibc file */

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "ibc");
    fn = (char *)malloc ((2 * strlen (simulation) + 12) * sizeof (char));
    sprintf (fn, "input/%s/%s.ibc", simulation, simulation);
    ibc_file = fopen (fn, "r");
    free (fn);

    if (ibc_file == NULL)
    {
        printf ("\n  Fatal Error: %s.ibc is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /*
     * start reading ibc_file 
     */
    fscanf (ibc_file, "%d %d", &DS->Num1BC, &DS->Num2BC);

    if (DS->Num1BC + DS->Num2BC > 0)
        DS->TSD_EleBC = (TSD *) malloc ((DS->Num1BC + DS->Num2BC) * sizeof (TSD));

    if (DS->Num1BC > 0)
    {
        /* For elements with Dirichilet Boundary Conditions */
        for (i = 0; i < DS->Num1BC; i++)
        {
            fscanf (ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index, &DS->TSD_EleBC[i].length);
            DS->TSD_EleBC[i].TS = (realtype **) malloc ((DS->TSD_EleBC[i].length) * sizeof (realtype *));
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
                DS->TSD_EleBC[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));

            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
            {
                fscanf (ibc_file, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_EleBC[i].TS[j][1]);
                timeinfo->tm_year = timeinfo->tm_year - 1900;
                timeinfo->tm_mon = timeinfo->tm_mon - 1;
                timeinfo->tm_sec = 0;
                rawtime = timegm (timeinfo);
                DS->TSD_EleBC[i].TS[j][0] = (realtype) rawtime;
            }
        }
    }

    if (DS->Num2BC > 0)
    {
        /* For elements with Neumann (non-natural) Boundary Conditions */
        for (i = DS->Num1BC; i < DS->Num1BC + DS->Num2BC; i++)
        {
            fscanf (ibc_file, "%s %d %d", DS->TSD_EleBC[i].name, &DS->TSD_EleBC[i].index, &DS->TSD_EleBC[i].length);

            DS->TSD_EleBC[i].TS = (realtype **) malloc ((DS->TSD_EleBC[i].length) * sizeof (realtype *));

            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
                DS->TSD_EleBC[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
            {
                fscanf (ibc_file, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_EleBC[i].TS[j][1]);
                timeinfo->tm_year = timeinfo->tm_year - 1900;
                timeinfo->tm_mon = timeinfo->tm_mon - 1;
                timeinfo->tm_sec = 0;
                rawtime = timegm (timeinfo);
                DS->TSD_EleBC[i].TS[j][0] = (realtype) rawtime;
            }
        }
    }
    fclose (ibc_file);
}

void ReadPara (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    FILE           *para_file;  /* Pointer to .para file */
    int             i, j;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    time_t          rawtime;
    struct tm      *timeinfo;

    int             NumTout;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    if (CS->Verbose)
        printf ("  Reading %s.%s\n", simulation, "para");
    fn = (char *)malloc ((2 * strlen (simulation) + 13) * sizeof (char));
    sprintf (fn, "input/%s/%s.para", simulation, simulation);
    para_file = fopen (fn, "r");
    free (fn);

    if (para_file == NULL)
    {
        printf ("\n  Fatal Error: %s.para is in use or does not exist!\n",
           simulation);
        exit (1);
    }

    /* start reading para_file */
    /* Set default values for parameters */
    CS->Ascii = 0;              /* YS */
    CS->Spinup = 0;             /* YS */
    CS->init_type = 0;
    DS->UnsatMode = 2;
    DS->SurfMode = 2;
    DS->RivMode = 2;
    CS->Solver = 2;
    CS->GSType = 1;
    CS->MaxK = 0;
    CS->delt = 0;
    CS->abstol = BADVAL;
    CS->reltol = BADVAL;
    CS->InitStep = BADVAL;
    CS->MaxStep = BADVAL;
    CS->ETStep = BADVAL;
    CS->StartTime = BADVAL;
    CS->EndTime = BADVAL;
    CS->outtype = BADVAL;
    CS->a = BADVAL;
    CS->b = BADVAL;

    CS->PrintGW = 0;
    CS->PrintSurf = 0;
    CS->PrintSnow = 0;
    CS->PrintRivStg = 0;
    CS->PrintRech = 0;
    CS->PrintIS = 0;
    CS->PrintUnsat = 0;
    for (j = 0; j < 3; j++)
        CS->PrintET[j] = 0;
    for (j = 0; j < 10; j++)
        CS->PrintRivFlx[j] = 0;

    /* Read through parameter file to find parameters */
    fgets (cmdstr, MAXSTRING, para_file);

    while (!feof (para_file))
    {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0')
        {
            sscanf (cmdstr, "%s", optstr);

            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#')
            {
                fgets (cmdstr, MAXSTRING, para_file);
                continue;
            }

            /* Get Model Parameters */
            if (strcasecmp ("INIT_MODE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->init_type);
            else if (strcasecmp ("ASCII_OUTPUT", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->Ascii);
            else if (strcasecmp ("SPINUP_MODE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->Spinup);
            else if (strcasecmp ("UNSAT_MODE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &DS->UnsatMode);
            else if (strcasecmp ("SAT_MODE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &DS->SurfMode);
            else if (strcasecmp ("RIV_MODE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &DS->RivMode);
            else if (strcasecmp ("SOLVER", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->Solver);
            else if (strcasecmp ("GSTYPE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->GSType);
            else if (strcasecmp ("MAXK", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->MaxK);
            else if (strcasecmp ("DELTA", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->delt);
            else if (strcasecmp ("ABSTOL", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->abstol);
            else if (strcasecmp ("RELTOL", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->reltol);
            else if (strcasecmp ("INIT_SOLVER_STEP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->InitStep);
            else if (strcasecmp ("MAX_SOLVER_STEP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->MaxStep);
            else if (strcasecmp ("LSM_STEP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->ETStep);
            else if (strcasecmp ("START", optstr) == 0)
            {
                sscanf (cmdstr, "%*s %d-%d-%d %d:%d", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min);
                timeinfo->tm_year = timeinfo->tm_year - 1900;
                timeinfo->tm_mon = timeinfo->tm_mon - 1;
                timeinfo->tm_sec = 0;
                rawtime = timegm (timeinfo);
                CS->StartTime = (realtype) rawtime;
            }
            else if (strcasecmp ("END", optstr) == 0)
            {
                sscanf (cmdstr, "%*s %d-%d-%d %d:%d", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min);
                timeinfo->tm_year = timeinfo->tm_year - 1900;
                timeinfo->tm_mon = timeinfo->tm_mon - 1;
                timeinfo->tm_sec = 0;
                rawtime = timegm (timeinfo);
                CS->EndTime = (realtype) rawtime;
            }
            else if (strcasecmp ("OUTPUT_TYPE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->outtype);
            else if (strcasecmp ("STEPSIZE_FACTOR", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->a);
            else if (strcasecmp ("MODEL_STEPSIZE", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->b);
            else if (strcasecmp ("GW", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintGW);
            else if (strcasecmp ("SURF", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintSurf);
            else if (strcasecmp ("SNOW", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintSnow);
            else if (strcasecmp ("RIVSTG", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivStg);
            else if (strcasecmp ("RECHARGE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRech);
            else if (strcasecmp ("CMC", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintIS);
            else if (strcasecmp ("UNSAT", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintUnsat);
            else if (strcasecmp ("EC", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintET[0]);
            else if (strcasecmp ("ETT", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintET[1]);
            else if (strcasecmp ("EDIR", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintET[2]);
            else if (strcasecmp ("RIVFLX0", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[0]);
            else if (strcasecmp ("RIVFLX1", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[1]);
            else if (strcasecmp ("RIVFLX2", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[2]);
            else if (strcasecmp ("RIVFLX3", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[3]);
            else if (strcasecmp ("RIVFLX4", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[4]);
            else if (strcasecmp ("RIVFLX5", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[5]);
            else if (strcasecmp ("RIVFLX6", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[6]);
            else if (strcasecmp ("RIVFLX7", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[7]);
            else if (strcasecmp ("RIVFLX8", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[8]);
            else if (strcasecmp ("RIVFLX9", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &CS->PrintRivFlx[9]);
            /* Unrecognized Parameter Flag */
            else
            {
                printf
                   ("\n  Parameter:%s cannot be recognized. Please see User's Manual for more details!\n",
                   optstr);
                exit (1);
            }
        }
        fgets (cmdstr, MAXSTRING, para_file);
    }

    fclose (para_file);

    if (CS->a != 1.0)
        NumTout = (int)(log (1. - (CS->EndTime - CS->StartTime) * (1. - CS->a) / CS->b) / log (CS->a));
    else
    {
        if ((CS->EndTime - CS->StartTime) / CS->b - ((int)(CS->EndTime - CS->StartTime) / CS->b) > 0)
            NumTout = (int)((CS->EndTime - CS->StartTime) / CS->b);
        else
            NumTout = (int)((CS->EndTime - CS->StartTime) / CS->b - 1);
    }

    CS->NumSteps = NumTout + 1;

    CS->Tout = (realtype *) malloc ((CS->NumSteps + 1) * sizeof (realtype));

    for (i = 0; i < CS->NumSteps + 1; i++)
    {
        if (i == 0)
            CS->Tout[i] = CS->StartTime;
        else
            CS->Tout[i] = CS->Tout[i - 1] + pow (CS->a, i) * CS->b;
    }

    if (CS->Tout[CS->NumSteps] < CS->EndTime)
        CS->Tout[CS->NumSteps] = CS->EndTime;
    if (CS->abstol == BADVAL)
    {
        printf ("\n  Fatal Error: Absolute Tolerance (ABSTOL) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->reltol == BADVAL)
    {
        printf ("\n  Fatal Error: Relative  Tolerance (RELTOL) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->InitStep == BADVAL)
    {
        printf ("\n  Fatal Error: Initial time-step (INIT_STEP) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->MaxStep == BADVAL)
    {
        printf ("\n  Fatal Error: Maximum time-step (MAX_STEP) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->StartTime == BADVAL)
    {
        printf ("\n  Fatal Error: Simulation start time (START yyyy-mm-dd hh:mm) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->EndTime == BADVAL)
    {
        printf ("\n  Fatal Error: Simulation end time (END yyyy-mm-dd hh:mm) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->ETStep == BADVAL)
    {
        printf ("\n  Fatal Error: Land surface model time-step (LSM_STEP) must be defined in .para file!\n");
        exit (1);
    }
    if (CS->outtype == BADVAL)
    {
        printf ("\n  Fatal Error: Output step-size type (OUTPUT_TYPE) must be defined in .para file!\n");
        exit (1);
    }
    if ((CS->outtype == 0) && (CS->a == BADVAL || CS->b == BADVAL))
    {
        printf ("\n  Fatal Error: Output step-size factor (A) and base step-size (B) must be defined in .para file!\n");
        exit (1);
    }
}

void ReadCalib (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    FILE           *global_calib;   /* Pointer to .calib file */

    if (CS->Verbose)
        printf ("  Reading calibration file\n");

    fn = (char *)malloc ((2 * strlen (simulation) + 14) * sizeof (char));
    sprintf (fn, "input/%s/%s.calib", simulation, simulation);
    global_calib = fopen (fn, "r");
    free (fn);

    if (global_calib == NULL)
    {
        printf ("\n  Fatal Error: %s.calib is in use or does not exist!\n", simulation);
        exit (1);
    }

    /* start reading calib_file */
    CS->Cal.KsatH = 1.0;
    CS->Cal.KsatV = 1.0;
    CS->Cal.infKsatV = 1.0;
    CS->Cal.macKsatH = 1.0;
    CS->Cal.macKsatV = 1.0;
    CS->Cal.infD = 1.0;
    CS->Cal.RzD = 1.0;
    CS->Cal.macD = 1.0;
    CS->Cal.Porosity = 1.0;
    CS->Cal.Alpha = 1.0;
    CS->Cal.Beta = 1.0;
    CS->Cal.vAreaF = 1.0;
    CS->Cal.hAreaF = 1.0;
    CS->Cal.VegFrac = 1.0;
    CS->Cal.Albedo = 1.0;
    CS->Cal.Rough = 1.0;
    CS->Cal.Prep = 1.0;
    CS->Cal.Temp = 1.0;
    DS->pcCal.Et0 = 1.0;
    DS->pcCal.Et1 = 1.0;
    DS->pcCal.Et2 = 1.0;
    CS->Cal.rivRough = 1.0;
    CS->Cal.rivKsatH = 1.0;
    CS->Cal.rivKsatV = 1.0;
    CS->Cal.rivbedThick = 1.0;
    CS->Cal.rivDepth = 1.0;
    CS->Cal.rivShapeCoeff = 1.0;

    CS->Cal.Rmin = 1.0;
    CS->Cal.ThetaRef = 1.0;
    CS->Cal.ThetaW = 1.0;
#ifdef _FLUX_PIHM_
    CS->Cal.TF = 1.0;
    CS->Cal.IS = 1.0;
    CS->Cal.Czil = 1.0;
    CS->Cal.fx_soil = 1.0;
    CS->Cal.fx_canopy = 1.0;
    CS->Cal.Rs_ref = 1.0;
    CS->Cal.h_s = 1.0;
#endif
#ifdef _RT_
    CS->Cal.PCO2 = 1.0;
    CS->Cal.Keq  = 1.0;
    CS->Cal.Site_den = 1.0;
    CS->Cal.SSA  = 1.0;
    CS->Cal.Prep_conc = 1.0;
#endif

    fgets (cmdstr, MAXSTRING, global_calib);

    while (!feof (global_calib))
    {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0')
        {
            sscanf (cmdstr, "%s", optstr);
            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#')
            {
                fgets (cmdstr, MAXSTRING, global_calib);
                continue;
            }
            /* Get calibration coefficients */
            if (strcasecmp ("KSATH", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.KsatH);
            else if (strcasecmp ("KSATV", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.KsatV);
            else if (strcasecmp ("KINF", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.infKsatV);
            else if (strcasecmp ("KMACSATH", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.macKsatH);
            else if (strcasecmp ("KMACSATV", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.macKsatV);
            else if (strcasecmp ("DINF", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.infD);
            else if (strcasecmp ("DROOT", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.RzD);
            else if (strcasecmp ("DMAC", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.macD);
            else if (strcasecmp ("POROSITY", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Porosity);
            else if (strcasecmp ("ALPHA", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Alpha);
            else if (strcasecmp ("BETA", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Beta);
            else if (strcasecmp ("MACVF", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.vAreaF);
            else if (strcasecmp ("MACHF", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.hAreaF);
            else if (strcasecmp ("VEGFRAC", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.VegFrac);
            else if (strcasecmp ("ALBEDO", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Albedo);
            else if (strcasecmp ("ROUGH", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Rough);
            else if (strcasecmp ("PRCP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Prep);
            else if (strcasecmp ("SFCTMP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Temp);
            else if (strcasecmp ("EC", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &DS->pcCal.Et0);
            else if (strcasecmp ("ETT", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &DS->pcCal.Et1);
            else if (strcasecmp ("EDIR", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &DS->pcCal.Et2);
            else if (strcasecmp ("ROUGH_RIV", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.rivRough);
            else if (strcasecmp ("KRIVH", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.rivKsatH);
            else if (strcasecmp ("KRIVV", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.rivKsatV);
            else if (strcasecmp ("BEDTHCK", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.rivbedThick);
            else if (strcasecmp ("RIV_DPTH", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.rivDepth);
            else if (strcasecmp ("RIV_WDTH", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.rivShapeCoeff);
            else if (strcasecmp ("RS", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Rmin);
            else if (strcasecmp ("WLTSMC", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.ThetaW);
            else if (strcasecmp ("REFSMC", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.ThetaRef);
#ifdef _FLUX_PIHM_
            else if (strcasecmp ("DRIP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.TF);
            else if (strcasecmp ("CMCMAX", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.IS);
            else if (strcasecmp ("CZIL", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Czil);
            else if (strcasecmp ("FXEXP", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.fx_soil);
            else if (strcasecmp ("CFACTR", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.fx_canopy);
            else if (strcasecmp ("RGL", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.Rs_ref);
            else if (strcasecmp ("HS", optstr) == 0)
                sscanf (cmdstr, "%*s %lf", &CS->Cal.h_s);
#endif
            /* Unrecognized Parameter Flag */
            else
            {
                printf
                   ("\n  Parameter: %s cannot be recognized. Please see User's Manual for more details!\n",
                   optstr);
                exit (1);
            }
        }
        fgets (cmdstr, MAXSTRING, global_calib);
    }

    /* finish reading calib file */
    fclose (global_calib);
}

void FreeData (Model_Data DS, Control_Data  CS)
{
    int             i, j, k;

    /*
     * free river
     */
    for (i = 0; i < DS->NumRivBC; i++)
    {
        for (j = 0; j < DS->TSD_Riv[i].length; j++)
            free (DS->TSD_Riv[i].TS[j]);
        free (DS->TSD_Riv[i].TS);
    }

    free (DS->Riv);
    free (DS->Riv_IC);
    free (DS->Riv_Shape);
    free (DS->Riv_Mat);
    free (DS->TSD_Riv);
    /*
     * free mesh
     */

    free (DS->Ele);
    free (DS->Node);
    /*
     * free att
     */
    free (DS->Ele_IC);
    /*
     * free soil
     */
    free (DS->Soil);
    /*
     * free geol
     */
    free (DS->Geol);
    /*
     * free lc
     */
    free (DS->LandC);

    for (j = 0; j < DS->NumLAI; j++)
    {
        for (k = 0; k < DS->TSD_lai[j].length; k++)
            free (DS->TSD_lai[j].TS[k]);
        free (DS->TSD_lai[j].TS);
    }
    free (DS->TSD_lai);

    /*
     * free forc
     */
    for (j = 0; j < DS->NumTS; j++)
    {
        for (k = 0; k < DS->TSD_meteo[j].length; k++)
            free (DS->TSD_meteo[j].TS[k]);
        free (DS->TSD_meteo[j].TS);
    }
    free (DS->TSD_meteo);

    free (DS->ISFactor);
    /*
     * free ibc
     */
    if (DS->Num1BC > 0)
        for (i = 0; i < DS->Num1BC; i++)
        {
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
                free (DS->TSD_EleBC[i].TS[j]);
            free (DS->TSD_EleBC[i].TS);
        }
    if (DS->Num2BC > 0)
        for (i = DS->Num1BC; i < DS->Num1BC + DS->Num2BC; i++)
        {
            for (j = 0; j < DS->TSD_EleBC[i].length; j++)
                free (DS->TSD_EleBC[i].TS[j]);
            free (DS->TSD_EleBC[i].TS);
        }


    if (DS->Num1BC + DS->Num2BC > 0)
        free (DS->TSD_EleBC);
    /*
     * free para
     */
    free (CS->Tout);
    /*
     * free initialize.c
     */
    for (i = 0; i < DS->NumEle; i++)
        free (DS->FluxSurf[i]);
    free (DS->FluxSurf);
    for (i = 0; i < DS->NumEle; i++)
        free (DS->FluxSub[i]);
    free (DS->FluxSub);
    for (i = 0; i < DS->NumEle; i++)
        free (DS->EleET[i]);
    free (DS->EleET);
    for (i = 0; i < DS->NumRiv; i++)
        free (DS->FluxRiv[i]);
    free (DS->FluxRiv);
    free (DS->EleNetPrep);
    free (DS->windH);
    free (DS->EleSurf);
    free (DS->EleGW);
    free (DS->EleUnsat);
    free (DS->EleMacAct);
    free (DS->RivStg);
    free (DS->ElePrep);
    free (DS->EleViR);
    free (DS->Recharge);
    free (DS->EleIS);
    free (DS->EleISmax);
    free (DS->EleISsnowmax);
    free (DS->EleSnow);
    free (DS->EleSnowGrnd);
    free (DS->EleSnowCanopy);
    free (DS->EleTF);
    free (DS->Albedo);
#ifdef _FLUX_PIHM_
    free (DS->SfcSat);
    free (DS->EleETsat);
    free (DS->EleFCR);
    free (DS->avg_inf);
    free (DS->avg_rech);
    for (i = 0; i < DS->NumEle; i++)
        free (DS->avg_subflux[i]);
    free (DS->avg_subflux);
#endif
    /*
     * free Print
     */
    for (i = 0; i < CS->NumPrint; i++)
    {
        free (CS->PCtrl[i].PrintVar);
        free (CS->PCtrl[i].buffer);
    }
    /*
     * free DummyY
     */
    free (DS->DummyY);
}
