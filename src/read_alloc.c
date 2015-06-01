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
    char            cmdstr[MAXSTRING];
    int             match;

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

    /*
     * Read river segment block
     */

    /* Read number of river segments */
    FindLine (riv_file, "BOF");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumRiv);
    if (match != 1)
    {
        printf ("Cannot read number of river segments!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    DS->Riv = (river_segment *) malloc (DS->NumRiv * sizeof (river_segment));

    /* Read river segment information */
    for (i = 0; i < DS->NumRiv; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d %d %d %d %d",
            &DS->Riv[i].index, &DS->Riv[i].FromNode, &DS->Riv[i].ToNode,
            &DS->Riv[i].down, &DS->Riv[i].LeftEle, &DS->Riv[i].RightEle,
            &DS->Riv[i].shape, &DS->Riv[i].material, &DS->Riv[i].IC,
            &DS->Riv[i].BC, &DS->Riv[i].reservoir);
        if (match != 11 || i != DS->Riv[i].index - 1)
        {
            printf ("Cannot read river segment information for the %dth segment!\n", i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river shape information
     */
    FindLine (riv_file, "SHAPE");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumRivShape);
    if (match != 1)
    {
        printf ("Cannot read number of river shapes!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    DS->Riv_Shape = (river_shape *) malloc (DS->NumRivShape * sizeof (river_shape));

    for (i = 0; i < DS->NumRivShape; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %d %lf",
            &DS->Riv_Shape[i].index, &DS->Riv_Shape[i].depth,
            &DS->Riv_Shape[i].interpOrd, &DS->Riv_Shape[i].coeff);
        if (match != 4 || i != DS->Riv_Shape[i].index -1 )
        {
            printf ("Cannot read river shape information for the %dth shape!\n", i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river material information
     */
    FindLine (riv_file, "MATERIAL");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumRivMaterial);
    if (match != 1)
    {
        printf ("Cannot read number of river materials!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    DS->Riv_Mat = (river_material *) malloc (DS->NumRivMaterial * sizeof (river_material));

    for (i = 0; i < DS->NumRivMaterial; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf",
            &DS->Riv_Mat[i].index, &DS->Riv_Mat[i].Rough,
            &DS->Riv_Mat[i].Cwr, &DS->Riv_Mat[i].KsatH,
            &DS->Riv_Mat[i].KsatV, &DS->Riv_Mat[i].bedThick);
        if (match != 6 || i != DS->Riv_Mat[i].index - 1)
        {
            printf ("Cannot read information for the %dth material!\n", i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river initial condition information
     */
    FindLine (riv_file, "IC");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumRivIC);
    if (match != 1)
    {
        printf ("Cannot read number of river materials!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    DS->Riv_IC = (river_IC *) malloc (DS->NumRivIC * sizeof (river_IC));

    /* Rewind and read river initial condition information */
    for (i = 0; i < DS->NumRivIC; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf", &DS->Riv_IC[i].index, &DS->Riv_IC[i].value);
        if (match != 2 || i != DS->Riv_IC[i].index - 1)
        {
            printf ("Cannot read information for the %dth initial condition!\n", i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river boundary condition block
     */
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &DS->NumRivBC);
    if (match != 1)
    {
        printf ("Cannot read number of river boundary conditions!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    DS->TSD_Riv = (TSD *) malloc (DS->NumRivBC * sizeof (TSD));

    for (i = 0; i < DS->NumRivBC; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &DS->TSD_Riv[i].index);
        if (match != 1 || i != DS->TSD_Riv[i].index - 1)
        {
            printf ("Cannot read information of the %dth river boudnary condition!\n", i);
            printf (".riv file format error!\n");
            exit (1);
        }
        NextLine (riv_file, cmdstr);
        NextLine (riv_file, cmdstr);
        DS->TSD_Riv[i].length = CountLine (riv_file, 2, "RIV_TS", "RES");
    }

    FindLine (riv_file, "BC");
    for (i = 0; i < DS->NumRivBC; i++)
    {
        NextLine (riv_file, cmdstr);
        NextLine (riv_file, cmdstr);
        NextLine (riv_file, cmdstr);

        DS->TSD_Riv[i].TS = (realtype **) malloc ((DS->TSD_Riv[i].length) * sizeof (realtype *));
        DS->TSD_Riv[i].iCounter = 0;
        for (j = 0; j < DS->TSD_Riv[i].length; j++)
        {
            DS->TSD_Riv[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));
            NextLine (riv_file, cmdstr);
            match = sscanf (cmdstr, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_Riv[i].TS[j][1]);
            if (match != 6)
            {
                printf (".riv file format error!\n");
                exit (1);
            }

            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            rawtime = timegm (timeinfo);
            DS->TSD_Riv[i].TS[j][0] = (realtype) rawtime;
        }
    }

    NextLine (riv_file, cmdstr);
    sscanf (cmdstr, "%*s %d", &DS->NumRes);
    /* Read Reservoir information */

    fclose (riv_file);
}

void ReadMesh (char *simulation, Model_Data DS, Control_Data CS)
{
    FILE           *mesh_file;  /* Pointer to .mesh file */
    char           *fn;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;

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
    NextLine (mesh_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumEle);
    if (match != 1)
    {
        printf ("Cannot read number of elements!\n");
        printf (".mesh file format error!\n");
        exit (1);
    }

    DS->Ele = (element *) malloc ((DS->NumEle + DS->NumRiv) * sizeof (element));

    /* read in elements information */
    for (i = 0; i < DS->NumEle; i++)
    {
        NextLine (mesh_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d", &DS->Ele[i].index,
            &DS->Ele[i].node[0], &DS->Ele[i].node[1], &DS->Ele[i].node[2],
            &DS->Ele[i].nabr[0], &DS->Ele[i].nabr[1], &DS->Ele[i].nabr[2]);
        if (match != 7 || i != DS->Ele[i].index - 1)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".mesh file format error!\n");
            exit (1);
        }
    }

    /* read in nodes information */
    NextLine (mesh_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumNode);
    if (match != 1)
    {
        printf ("Cannot read number of nodes!\n");
        printf (".mesh file format error!\n");
        exit (1);
    }

    DS->Node = (nodes *) malloc (DS->NumNode * sizeof (nodes));
    for (i = 0; i < DS->NumNode; i++)
    {
        NextLine (mesh_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf", &(DS->Node[i].index),
            &DS->Node[i].x, &DS->Node[i].y,
            &DS->Node[i].zmin, &DS->Node[i].zmax);
        if (match != 5 || i != DS->Node[i].index - 1)
        {
            printf ("Cannot read information of the %dth node!\n", i + 1);
            printf (".mesh file format error!\n");
            exit (1);
        }
    }

    /* finish reading mesh_files */
    fclose (mesh_file);
}

void ReadAtt (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    int             i;
    FILE           *att_file;   /* Pointer to .att file */
    char            cmdstr[MAXSTRING];
    int             match;

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

    NextLine (att_file, cmdstr);
    for (i = 0; i < DS->NumEle; i++)
    {
        NextLine (att_file, cmdstr);
        match = sscanf (cmdstr, "%*d %d %d %d %lf %lf %lf %lf %lf %d %d %d %d %d %d %d",
            &DS->Ele[i].soil, &DS->Ele[i].geol, &DS->Ele[i].LC,
            &DS->Ele_IC[i].interception, &DS->Ele_IC[i].snow, &DS->Ele_IC[i].surf,
            &DS->Ele_IC[i].unsat, &DS->Ele_IC[i].sat,
            &DS->Ele[i].meteo, &DS->Ele[i].LAI,
            &DS->Ele[i].source,
            &DS->Ele[i].BC[0], &DS->Ele[i].BC[1], &DS->Ele[i].BC[2],
            &DS->Ele[i].Macropore);
        if (match != 15)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".att file format error!\n");
            exit (1);
        }
    }

    /* finish reading att_files */
    fclose (att_file);
}

void ReadSoil (char *simulation, Model_Data DS, Control_Data CS)
{
    FILE           *soil_file;  /* Pointer to .soil file */
    char           *fn;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;

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
    NextLine (soil_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumSoil);
    if (match != 1)
    {
        printf ("Cannot read number of soil types!\n");
        printf (".soil file format error!\n");
        exit (1);
    }

    DS->Soil = (soils *) malloc (DS->NumSoil * sizeof (soils));

    for (i = 0; i < DS->NumSoil; i++)
    {
        NextLine (soil_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &DS->Soil[i].index,
            &DS->Soil[i].KsatV,
            &DS->Soil[i].ThetaS, &DS->Soil[i].ThetaR,
            &DS->Soil[i].infD,
            &DS->Soil[i].Alpha, &DS->Soil[i].Beta,
            &DS->Soil[i].hAreaF, &DS->Soil[i].macKsatV,
            &DS->Soil[i].qtz);
        if (match != 10 || i != DS->Soil[i].index - 1)
        {
            printf ("Cannot read information of the %dth soil type!\n", i + 1);
            printf (".soil file format error!\n");
            exit (1);
        }
    }

    fclose (soil_file);
}

void ReadGeol (char *simulation, Model_Data DS, Control_Data CS)
{
    char           *fn;
    FILE           *geol_file;  /* Pointer to .geol file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;

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
    NextLine (geol_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumGeol);
    if (match != 1)
    {
        printf ("Cannot read number of geology types!\n");
        printf (".geol file format error!\n");
        exit (1);
    }

    DS->Geol = (geol *) malloc (DS->NumGeol * sizeof (geol));

    for (i = 0; i < DS->NumGeol; i++)
    {
        NextLine (geol_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &DS->Geol[i].index,
            &DS->Geol[i].KsatH, &DS->Geol[i].KsatV,
            &DS->Geol[i].ThetaS, &DS->Geol[i].ThetaR,
            &DS->Geol[i].Alpha, &DS->Geol[i].Beta,
            &DS->Geol[i].vAreaF, &DS->Geol[i].macKsatH, &DS->Geol[i].macD);
        if (match != 10 || i != DS->Geol[i].index - 1)
        {
            printf ("Cannot read information of the %dth geology type!\n", i + 1);
            printf (".geol file format error!\n");
            exit (1);
        }
    }

    fclose (geol_file);
}

void ReadLC (char *simulation, Model_Data DS, Control_Data CS)
{
    FILE           *lc_file;    /* Pointer to .lc file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;

    if (CS->Verbose)
        printf ("  Reading vegprmt.tbl\n");
    lc_file = fopen ("input/vegprmt.tbl", "r");

    if (lc_file == NULL)
    {
        printf ("\n  Fatal Error: land cover file input/vegprmt.tbl is in use or does not exist!\n");
        exit (1);
    }

    /* start reading land cover file */
    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%d", &DS->NumLC);
    if (match != 1)
    {
        printf ("Cannot read number of landcover types!\n");
        printf ("Land cover file format error!\n");
        exit (1);
    }

    DS->LandC = (LC *) malloc (DS->NumLC * sizeof (LC));

    for (i = 0; i < DS->NumLC; i++)
    {
        NextLine (lc_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &DS->LandC[i].index,
            &DS->LandC[i].VegFrac, &DS->LandC[i].RzD, 
            &DS->LandC[i].Rmin, &DS->LandC[i].Rs_ref,
            &DS->LandC[i].h_s,
            &DS->LandC[i].snup, 
            &DS->LandC[i].LAImin, &DS->LandC[i].LAImax,
            &DS->LandC[i].Emiss_min, &DS->LandC[i].Emiss_max,
            &DS->LandC[i].Albedo_min, &DS->LandC[i].Albedo_max,
            &DS->LandC[i].z0_min, &DS->LandC[i].z0_max,
            &DS->LandC[i].Rough);
        if (match != 16 || i != DS->LandC[i].index - 1)
        {
            printf ("Cannot read information of the %dth landcover type!\n", i + 1);
            printf ("Landcover file format error!\n");
            exit (1);
        }
    }

    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %lf", &DS->Tref);
    if (match != 1)
    {
        printf ("Cannot read information of optimal temperature!\n");
        printf ("Landcover file format error!\n");
        exit (1);
    }
    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %lf", &DS->fx_canopy);
    if (match != 1)
    {
        printf ("Cannot read information of canopy evaporation rate!\n");
        printf ("Landcover file format error!\n");
        exit (1);
    }
    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %lf", &DS->Rmax);
    if (match != 1)
    {
        printf ("Cannot read information of canopy cuticular resistance!\n");
        printf ("Landcover file format error!\n");
        exit (1);
    }
    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &(DS->bare));
    if (match != 1)
    {
        printf ("Cannot read information of Bare soil type!\n");
        printf ("Landcover file format error!\n");
        exit (1);
    }

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
    int             match;

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

    FindLine (forc_file, "BOF");
    NextLine (forc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &DS->NumTS);
    if (match != 1)
    {
        printf ("Cannot read number of meteorological forcing time series!\n");
        printf (".forc file format error!\n");
        exit (1);
    }

    DS->TSD_meteo = (TSD *) malloc (DS->NumTS * sizeof (TSD));
    DS->windH = (realtype *) malloc (DS->NumTS * sizeof (realtype));

    for (i = 0; i < DS->NumTS; i++)
    {
        NextLine (forc_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d %*s %lf", &DS->TSD_meteo[i].index, &DS->windH[i]);
        if (match != 2 || i != DS->TSD_meteo[i].index - 1)
        {
            printf ("Cannot read information of the %dth forcing series!\n", i);
            printf (".forc file format error!\n");
            exit (1);
        }
        NextLine (forc_file, cmdstr);
        NextLine (forc_file, cmdstr);
        DS->TSD_meteo[i].length = CountLine (forc_file, 1, "METEO_TS");
    }

    FindLine (forc_file, "NUM_METEO_TS");
    for (i = 0; i < DS->NumTS; i++)
    {
        NextLine (forc_file, cmdstr);
        NextLine (forc_file, cmdstr);
        NextLine (forc_file, cmdstr);

        DS->TSD_meteo[i].TS = (realtype **) malloc ((DS->TSD_meteo[i].length) * sizeof (realtype *));
        DS->TSD_meteo[i].iCounter = 0;
        for (j = 0; j < DS->TSD_meteo[i].length; j++)
        {
            DS->TSD_meteo[i].TS[j] = (realtype *) malloc ((NumForcing + 1)* sizeof (realtype));
            NextLine (forc_file, cmdstr);
            match = sscanf (cmdstr, "%d-%d-%d %d:%d %lf %lf %lf %lf %lf %lf %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_meteo[i].TS[j][1], &DS->TSD_meteo[i].TS[j][2], &DS->TSD_meteo[i].TS[j][3], &DS->TSD_meteo[i].TS[j][4], &DS->TSD_meteo[i].TS[j][5], &DS->TSD_meteo[i].TS[j][6], &DS->TSD_meteo[i].TS[j][7]);
            if (match != 12)
            {
                printf (".forc file format error (Line %d)!\n", j + 1);
                exit (1);
            }
            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            timeinfo->tm_sec = 0;
            rawtime = timegm (timeinfo);
            DS->TSD_meteo[i].TS[j][0] = (realtype) rawtime;
        }
    }

    fclose (forc_file);

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
        FindLine (lai_file, "BOF");
        NextLine (lai_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &DS->NumLAI);
        if (match != 1)
        {
            printf ("Cannot read number of LAI time series!\n");
            printf (".lai file format error!\n");
            exit (1);
        }


        DS->TSD_lai = (TSD *) malloc (DS->NumLAI * sizeof (TSD));

        for (i = 0; i < DS->NumLAI; i++)
        {
            NextLine (lai_file, cmdstr);
            match = sscanf (cmdstr, "%*s %d", &DS->TSD_lai[i].index);
            if (match != 1 || i != DS->TSD_lai[i].index - 1)
            {
                printf ("Cannot read information of the %dth LAI series!\n", i);
                printf (".lai file format error!\n");
                exit (1);
            }
            NextLine (lai_file, cmdstr);
            NextLine (lai_file, cmdstr);
            DS->TSD_lai[i].length = CountLine (lai_file, 1, "LAI_TS");
        }

        FindLine (lai_file, "NUM_LAI_TS");
        for (i = 0; i < DS->NumLAI; i++)
        {
            NextLine (lai_file, cmdstr);
            NextLine (lai_file, cmdstr);
            NextLine (lai_file, cmdstr);

            DS->TSD_lai[i].TS = (realtype **) malloc ((DS->TSD_lai[i].length) * sizeof (realtype *));
            DS->TSD_lai[i].iCounter = 0;
            for (j = 0; j < DS->TSD_lai[i].length; j++)
            {
                DS->TSD_lai[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));
                NextLine (lai_file, cmdstr);
                match = sscanf (cmdstr, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_lai[i].TS[j][1]);
                if (match != 6)
                {
                    printf (".forc file format error!\n");
                    exit (1);
                }
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
    char            cmdstr[MAXSTRING];
    int             match;

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
    FindLine (ibc_file, "BOF");
    NextLine (ibc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &DS->NumBC);
    if (match != 1)
    {
        printf ("Cannot read number of boundary condition time series!\n");
        printf (".ibc file format error!\n");
        exit (1);
    }
        
    DS->TSD_EleBC = (TSD *) malloc (DS->NumBC * sizeof (TSD));

    for (i = 0; i < DS->Num1BC; i++)
    {
        NextLine (ibc_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &DS->TSD_EleBC[i].index);
        if (match != 1 || i != DS->TSD_EleBC[i].index - 1)
        {
            printf ("Cannot read information of the %dth boundary condition series!\n", i);
            printf (".ibc file format error!\n");
            exit (1);
        }
        NextLine (ibc_file, cmdstr);
        NextLine (ibc_file, cmdstr);

        DS->TSD_EleBC[i].length = CountLine (ibc_file, 1, "BC_TS");
    }

    FindLine (ibc_file, "NUM_BC_TS");
    for (i = 0; i < DS->NumBC; i++)
    {
        NextLine (ibc_file, cmdstr);
        NextLine (ibc_file, cmdstr);
        NextLine (ibc_file, cmdstr);

        DS->TSD_EleBC[i].TS = (realtype **) malloc ((DS->TSD_EleBC[i].length) * sizeof (realtype *));
        DS->TSD_EleBC[i].iCounter = 0;
        for (j = 0; j < DS->TSD_EleBC[i].length; j++)
        {
            DS->TSD_EleBC[i].TS[j] = (realtype *) malloc (2 * sizeof (realtype));
            NextLine (ibc_file, cmdstr);
            match = sscanf (cmdstr, "%d-%d-%d %d:%d %lf", &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour, &timeinfo->tm_min, &DS->TSD_EleBC[i].TS[j][1]);
            if (match != 6)
            {
                printf (".ibc file format error (Line %d)!\n", j + 1);
                exit (1);
            }
            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            timeinfo->tm_sec = 0;
            rawtime = timegm (timeinfo);
            DS->TSD_EleBC[i].TS[j][0] = (realtype) rawtime;
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
