
/******************************************************************************
 * File		:   lsm_func.c                                                    *
 * Function	:   LSM related functions                                         *
 * Model	:   Flux-PIHM                                                     *
 * Version	:   September, 2014 (based on PIHM 2.2)                           *
 * Developer of Flux-PIHM:	Yuning Shi (yshi@psu.edu)                         *
 * Developer of PIHM2.0:	Mukesh Kumar (muk139@psu.edu)                     * 
 * Developer of PIHM1.0:	Yizhong Qu   (quyizhong@gmail.com)                * 
 *----------------------------------------------------------------------------*
 * For questions or comments, please contact                                  *
 *      --> Yuning Shi (yshi@psu.edu)                                         *
 * This code is free for research purpose only.                               *
 * Please provide relevant references if you use this code in your research   *
 *  work                                                                      *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "../pihm.h"               /* Data Model and Variable Declarations     */
#include "../spa/spa.h"
#include "noah.h"
#include "flux_pihm.h"

void LSM_read (char *filename, LSM_STRUCT LSM)
{
    int             i, j;
    char           *fn;
    char           *projectname;
    char           *token, *tempname;
    FILE           *lsm_file;
    int             ensemble_mode;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];

    /*
     * Detect if model is running in ensemble mode
     */
    tempname = (char *)malloc ((strlen (filename) + 1) * sizeof (char));
    strcpy (tempname, filename);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        projectname = (char *)malloc ((strlen (token) + 1) * sizeof (char));
        strcpy (projectname, token);
        ensemble_mode = 1;
    }
    else
    {
        projectname = (char *)malloc ((strlen (filename) + 1)
           * sizeof (char));
        strcpy (projectname, filename);
        ensemble_mode = 0;
    }
    free (tempname);

    /*
     * Open *.lsm file
     */
    if (ensemble_mode == 0)
        printf ("\n 10) reading %s.%-5s ...\t\t", projectname, "lsm");
    fn = (char *)malloc ((strlen (projectname) + 11) * sizeof (char));
    sprintf (fn, "input/%s.lsm", projectname);
    lsm_file = fopen (fn, "r");
    free (fn);

    if (lsm_file == NULL)
    {
        printf ("\n  Fatal Error: %s.lsm is in use or does not exist!\n",
           projectname);
        exit (1);
    }
    free (projectname);

    /*
     * Start reading lsm_file
     */
    fscanf (lsm_file, "%*s %lf", &LSM->LATITUDE);
    fscanf (lsm_file, "%*s %lf", &LSM->LONGITUDE);
    fscanf (lsm_file, "%*s %d", &LSM->STD_NSOIL);
    fscanf (lsm_file, "%*s");
    LSM->STD_SLDPTH = (double *)malloc (LSM->STD_NSOIL * sizeof (double));
    for (i = 0; i < LSM->STD_NSOIL; i++)
        fscanf (lsm_file, "%lf", &LSM->STD_SLDPTH[i]);
    fscanf (lsm_file, "%*s %d", &LSM->RAD_MODE);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.SBETA_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.FXEXP_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.CSOIL_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.SALP_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.FRZK_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.ZBOT_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.TBOT_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.CZIL_DATA);
    fscanf (lsm_file, "%*s %lf", &LSM->GENPRMT.LVCOEF_DATA);

    LSM->PRINT_T1 = 0;
    LSM->PRINT_STC = 0;
    LSM->PRINT_SMC = 0;
    LSM->PRINT_SH2O = 0;
    LSM->PRINT_SNOWH = 0;
    LSM->PRINT_ALBEDO = 0;
    LSM->PRINT_LE = 0;
    LSM->PRINT_SH = 0;
    LSM->PRINT_G = 0;
    LSM->PRINT_ETP = 0;

    fgets (cmdstr, MAXSTRING, lsm_file);

    while (!feof (lsm_file))
    {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0')
        {
            sscanf (cmdstr, "%s", optstr);

            /*
             * Handle case of comment line in which '#' is indented 
             */
            if (optstr[0] == '#')
            {
                fgets (cmdstr, MAXSTRING, lsm_file);
                continue;
            }

            /*
             * Get Model Parameters 
             */
            if (strcasecmp ("T1", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_T1);
            else if (strcasecmp ("STC", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_STC);
            else if (strcasecmp ("SMC", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_SMC);
            else if (strcasecmp ("SH2O", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_SH2O);
            else if (strcasecmp ("SNOWH", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_SNOWH);
            else if (strcasecmp ("ALBEDO", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_ALBEDO);
            else if (strcasecmp ("LE", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_LE);
            else if (strcasecmp ("SH", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_SH);
            else if (strcasecmp ("G", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_G);
            else if (strcasecmp ("ETP", optstr) == 0)
                sscanf (cmdstr, "%*s %d", &LSM->PRINT_ETP);
            /*
             * Unrecognized Parameter Flag 
             */
            else
            {
                printf
                   ("\n  Parameter:%s cannot be recognized. Please see User's Manual for more details!\n",
                   optstr);
                exit (1);
            }
        }
        fgets (cmdstr, MAXSTRING, lsm_file);
    }

    fclose (lsm_file);
    if (ensemble_mode == 0)
        printf ("done.\n");
}

void
LSM_initialize (char *filename, Model_Data DS, Control_Data * CS,
   LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;
    int             i, j, k, KZ;
    int             icounter;
    double          ZSOIL[LSM->STD_NSOIL];
    double          AquiferDepth;
    double          a_x, a_y, b_x, b_y, c_x, c_y, distX, distY;
    double          a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax;
    double          vector1[3], vector2[3], normal_vector[3], vector[3], H, c,
       se, ce;
    int             nodes[2];
    double          x1, y1, z1, x2, y2, z2, xc, yc, zc;
    double          c1, c2, ce1, ce2, se1, se2, phi1, phi2;
    int             ind, ind1, ind2;
    char           *fn;
    FILE           *init_file;


    if (CS->Verbose)
    {
        printf ("\nLSM vegetation parameters\n");
        printf
           ("\tSHDFAC\tNROOT\tRS\tRGL\tHS\tSNUP\tLAIMAX\tLAIMIN\tEMISSMIN\tEMISSMAX\tALBEDOMIN\tALBEDOMAX\tZ0MIN\tZ0MAX\tCMCFACTR\n");
    }

    LSM->VEGTBL.LUCATS = DS->NumLC;
    for (i = 0; i < LSM->VEGTBL.LUCATS; i++)
    {
        LSM->VEGTBL.SHDTBL[i] =
           (double)(CS->Cal.VegFrac * DS->LandC[i].VegFrac);
        LSM->VEGTBL.NROTBL[i] =
           FindLayer (LSM, (double)(CS->Cal.RzD * DS->LandC[i].RzD));
        /*
         * RSTBL: convert from day m-1 to s m-1 
         */
        LSM->VEGTBL.RSTBL[i] =
           (double)(CS->Cal.Rmin * DS->LandC[i].Rmin * 24. * 3600.);
        /*
         * RGLTBL: convert from J day-1 m-2 to W m-2 
         */
        LSM->VEGTBL.RGLTBL[i] =
           (double)(CS->Cal.Rs_ref * DS->LandC[i].Rs_ref / 24. / 3600.);
        LSM->VEGTBL.HSTBL[i] = (double)(CS->Cal.h_s * DS->LandC[i].h_s);
        LSM->VEGTBL.SNUPTBL[i] = (double)DS->LandC[i].snup;
        LSM->VEGTBL.LAIMINTBL[i] = (double)(DS->LandC[i].LAImin);
        LSM->VEGTBL.LAIMAXTBL[i] = (double)(DS->LandC[i].LAImax);
        LSM->VEGTBL.EMISSMINTBL[i] = (double)(DS->LandC[i].Emiss_min);
        LSM->VEGTBL.EMISSMAXTBL[i] = (double)(DS->LandC[i].Emiss_max);
        LSM->VEGTBL.ALBEDOMINTBL[i] =
           (double)(CS->Cal.Albedo * DS->LandC[i].Albedo_min);
        LSM->VEGTBL.ALBEDOMAXTBL[i] =
           (double)(CS->Cal.Albedo * DS->LandC[i].Albedo_max);
        LSM->VEGTBL.Z0MINTBL[i] = (double)DS->LandC[i].z0_min;
        LSM->VEGTBL.Z0MAXTBL[i] = (double)DS->LandC[i].z0_max;
        LSM->VEGTBL.CMCFACTRTBL[i] = (double)(CS->Cal.IS * DS->ISFactor[i]);
        if (CS->Verbose)
        {
            printf ("%-3d\t", i + 1);
            printf ("%-3.2f\t", LSM->VEGTBL.SHDTBL[i]);
            printf ("%-2d\t", LSM->VEGTBL.NROTBL[i]);
            printf ("%-6.2f\t", LSM->VEGTBL.RSTBL[i]);
            printf ("%-6.2f\t", LSM->VEGTBL.RGLTBL[i]);
            printf ("%-5.2f\t", LSM->VEGTBL.HSTBL[i]);
            printf ("%-4.2f\t", LSM->VEGTBL.SNUPTBL[i]);
            printf ("%-4.2f\t", LSM->VEGTBL.LAIMINTBL[i]);
            printf ("%-4.2f\t", LSM->VEGTBL.LAIMAXTBL[i]);
            printf ("%-4.2f\t\t", LSM->VEGTBL.EMISSMINTBL[i]);
            printf ("%-4.2f\t\t", LSM->VEGTBL.EMISSMAXTBL[i]);
            printf ("%-4.2f\t\t", LSM->VEGTBL.ALBEDOMINTBL[i]);
            printf ("%-4.2f\t\t", LSM->VEGTBL.ALBEDOMAXTBL[i]);
            printf ("%-4.2f\t", LSM->VEGTBL.Z0MINTBL[i]);
            printf ("%-4.2f\t", LSM->VEGTBL.Z0MAXTBL[i]);
            printf ("%-6.5f\n", LSM->VEGTBL.CMCFACTRTBL[i]);
        }
    }

    /*
     * TOPT: convert from degree C to Kalvin 
     */
    LSM->VEGTBL.TOPT_DATA = 273.15 + (double)DS->Tref;
    LSM->VEGTBL.CFACTR_DATA = (double)(CS->Cal.fx_canopy * DS->fx_canopy);
    LSM->VEGTBL.RSMAX_DATA = DS->Rmax;
    LSM->VEGTBL.BARE = DS->bare;

    if (CS->Verbose)
    {
        printf ("TOPT_DATA\n%-.1f\n", LSM->VEGTBL.TOPT_DATA);
        printf ("CFACTR_DATA\n%-.2f\n", LSM->VEGTBL.CFACTR_DATA);
        printf ("RSMAX_DATA\n%-.1f\n", LSM->VEGTBL.RSMAX_DATA);
        printf ("BARE\n%-d\n", LSM->VEGTBL.BARE);
    }

    LSM->SOILTBL.SLCATS = DS->NumSoil;
    if (CS->Verbose)
    {
        printf ("\nLSM soil parameters\n");
        printf
           ("\tDRYSMC\tMAXSMC\tREFSMC\tSATDK\t\tWLTSMC\tQTZ\tALPHA\tBETA\tMINSMC\tMACKSAT\t\tAREAF\tNMACD\n");
    }
    for (i = 0; i < LSM->SOILTBL.SLCATS; i++)
    {
        LSM->SOILTBL.DRYSMC[i] =
           (double)(CS->Cal.ThetaW * (DS->Soil[i].ThetaW -
              DS->Soil[i].ThetaR) + DS->Soil[i].ThetaR);
        LSM->SOILTBL.MAXSMC[i] =
           (double)(CS->Cal.Porosity * (DS->Geol[i].ThetaS -
              DS->Geol[i].ThetaR) + DS->Geol[i].ThetaR);
        LSM->SOILTBL.REFSMC[i] =
           (double)(CS->Cal.ThetaRef * (DS->Soil[i].ThetaRef -
              DS->Soil[i].ThetaR) + DS->Geol[i].ThetaR);
        /*
         * SATDK: convert from m day-1 to m s-1 
         */
        LSM->SOILTBL.SATDK[i] =
           (double)(CS->Cal.KsatV * DS->Geol[i].KsatV / 24. / 3600.);
        LSM->SOILTBL.WLTSMC[i] =
           (double)(CS->Cal.ThetaW * (DS->Soil[i].ThetaW -
              DS->Soil[i].ThetaR) + DS->Soil[i].ThetaR);
        LSM->SOILTBL.QTZ[i] = (double)DS->Soil[i].qtz;

        LSM->SOILTBL.VGA[i] = (double)(CS->Cal.Alpha * DS->Soil[i].Alpha);
        LSM->SOILTBL.VGB[i] = (double)(CS->Cal.Beta * DS->Soil[i].Beta);
        LSM->SOILTBL.MINSMC[i] = (double)DS->Soil[i].ThetaR;
        LSM->SOILTBL.MACKSAT[i] =
           (double)(CS->Cal.macKsatV * DS->Soil[i].macKsatV / 24. / 3600.);
        LSM->SOILTBL.AREAF[i] = (double)(CS->Cal.hAreaF * DS->Soil[i].hAreaF);
        LSM->SOILTBL.NMACD[i] =
           FindLayer (LSM, (double)(CS->Cal.macD * DS->Geol[i].macD));
        if (CS->Verbose)
        {
            printf ("%-3d\t", i + 1);
            printf ("%-5.3f\t", LSM->SOILTBL.DRYSMC[i]);
            printf ("%-5.3f\t", LSM->SOILTBL.MAXSMC[i]);
            printf ("%-5.3f\t", LSM->SOILTBL.REFSMC[i]);
            printf ("%-5.3E\t", LSM->SOILTBL.SATDK[i]);
            printf ("%-5.3f\t", LSM->SOILTBL.WLTSMC[i]);
            printf ("%-4.2f\t", LSM->SOILTBL.QTZ[i]);
            printf ("%-6.3f\t", LSM->SOILTBL.VGA[i]);
            printf ("%-6.3f\t", LSM->SOILTBL.VGB[i]);
            printf ("%-5.3f\t", LSM->SOILTBL.MINSMC[i]);
            printf ("%-5.3E\t", LSM->SOILTBL.MACKSAT[i]);
            printf ("%-5.3f\t", LSM->SOILTBL.AREAF[i]);
            printf ("%-2d\n", LSM->SOILTBL.NMACD[i]);
        }
    }

    //  LSM->GENPRMT.SBETA_DATA = LSM->SBETA;
    LSM->GENPRMT.FXEXP_DATA =
       (double)CS->Cal.fx_soil * LSM->GENPRMT.FXEXP_DATA;
    //  LSM->GENPRMT.CSOIL_DATA = ;
    //  LSM->GENPRMT.SALP_DATA = DS->salp;
    //  LSM->GENPRMT.FRZK_DATA = DS->frzk;
    //  LSM->GENPRMT.ZBOT_DATA = DS->zbot;
    LSM->GENPRMT.CZIL_DATA = (double)(CS->Cal.Czil) * LSM->GENPRMT.CZIL_DATA;
    //  LSM->GENPRMT.LVCOEF_DATA = ;

    if (CS->Verbose)
    {
        printf ("\nLSM general parameters\n");
        printf ("SBETA_DATA\t%-4.2f\n", LSM->GENPRMT.SBETA_DATA);
        printf ("FXEXP_DATA\t%-4.2f\n", LSM->GENPRMT.FXEXP_DATA);
        printf ("CSOIL_DATA\t%-4.2f\n", LSM->GENPRMT.CSOIL_DATA);
        printf ("SALP_DATA\t%-4.2f\n", LSM->GENPRMT.SALP_DATA);
        printf ("FRZK_DATA\t%-4.2f\n", LSM->GENPRMT.FRZK_DATA);
        printf ("ZBOT_DATA\t%-4.2f\n", LSM->GENPRMT.ZBOT_DATA);
        printf ("CZIL_DATA\t%-4.2f\n", LSM->GENPRMT.CZIL_DATA);
        printf ("LVCOEF_DATA\t%-4.2f\n", LSM->GENPRMT.LVCOEF_DATA);
    }

    LSM->GRID = (GRID_TYPE *) malloc (DS->NumEle * sizeof (GRID_TYPE));

    ZSOIL[0] = LSM->STD_SLDPTH[0];
    for (j = 0; j < LSM->STD_NSOIL; j++)
        ZSOIL[j] = ZSOIL[j - 1] + LSM->STD_SLDPTH[j];

    //  LSM->SLDPTH = (double **)malloc(DS->NumEle * sizeof(double *));

    for (i = 0; i < DS->NumEle; i++)
    {
        NOAH = &(LSM->GRID[i]);
        /*
         * Initialize model grid soil depths
         */
        NOAH->SLDPTH =
           (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
        AquiferDepth = (double)(DS->Ele[i].zmax - DS->Ele[i].zmin);

        if (AquiferDepth <= ZSOIL[0])
        {
            NOAH->SLDPTH[0] = AquiferDepth;
            NOAH->NSOIL = 1;
            for (j = 1; j < LSM->STD_NSOIL + 1; j++)
                NOAH->SLDPTH[j] = BADVAL;
        }
        else if (AquiferDepth <= ZSOIL[LSM->STD_NSOIL - 1])
        {
            for (j = 1; j < LSM->STD_NSOIL + 1; j++)
            {
                if (AquiferDepth <= ZSOIL[j])
                {
                    for (k = 0; k < j; k++)
                        NOAH->SLDPTH[k] = LSM->STD_SLDPTH[k];
                    NOAH->SLDPTH[j] = AquiferDepth - ZSOIL[j - 1];
                    NOAH->NSOIL = j + 1;

                    /*
                     * The following calculations gurantee that each layer is
                     * thicker than the layer on top
                     */
                    if (NOAH->SLDPTH[j] < NOAH->SLDPTH[j - 1])
                    {
                        NOAH->SLDPTH[j - 1] =
                           NOAH->SLDPTH[j - 1] + NOAH->SLDPTH[j];
                        NOAH->SLDPTH[j] = BADVAL;
                        NOAH->NSOIL = NOAH->NSOIL - 1;
                    }
                    for (k = j + 1; k < LSM->STD_NSOIL + 1; k++)
                        NOAH->SLDPTH[k] = BADVAL;
                    break;
                }
            }
        }
        else
        {
            for (j = 0; j < LSM->STD_NSOIL; j++)
                NOAH->SLDPTH[j] = LSM->STD_SLDPTH[j];
            NOAH->SLDPTH[LSM->STD_NSOIL] =
               AquiferDepth - ZSOIL[LSM->STD_NSOIL - 1];
            NOAH->NSOIL = LSM->STD_NSOIL + 1;
            if (NOAH->SLDPTH[LSM->STD_NSOIL] <
               NOAH->SLDPTH[LSM->STD_NSOIL - 1])
            {
                NOAH->SLDPTH[LSM->STD_NSOIL - 1] =
                   NOAH->SLDPTH[LSM->STD_NSOIL - 1] +
                   NOAH->SLDPTH[LSM->STD_NSOIL];
                NOAH->SLDPTH[LSM->STD_NSOIL] = BADVAL;
                NOAH->NSOIL = NOAH->NSOIL - 1;
            }
        }

        /*
         * if (CS->Verbose)
         * {
         * printf("Ele %d %d layers, soil depth:\t", i, NOAH->NSOIL);
         * for (j = 0; j < LSM->STD_NSOIL + 1; j++)
         * printf("%lf\t", NOAH->SLDPTH[j]);
         * printf("\n");
         * }
         */
        /*
         * Initialize topographic radiation related parameters
         */
        if (LSM->RAD_MODE == 1)
        {
            a_x = (double)DS->Node[DS->Ele[i].node[0] - 1].x;
            b_x = (double)DS->Node[DS->Ele[i].node[1] - 1].x;
            c_x = (double)DS->Node[DS->Ele[i].node[2] - 1].x;
            a_y = (double)DS->Node[DS->Ele[i].node[0] - 1].y;
            b_y = (double)DS->Node[DS->Ele[i].node[1] - 1].y;
            c_y = (double)DS->Node[DS->Ele[i].node[2] - 1].y;

            a_zmin = (double)DS->Node[DS->Ele[i].node[0] - 1].zmin;
            b_zmin = (double)DS->Node[DS->Ele[i].node[1] - 1].zmin;
            c_zmin = (double)DS->Node[DS->Ele[i].node[2] - 1].zmin;
            a_zmax = (double)DS->Node[DS->Ele[i].node[0] - 1].zmax;
            b_zmax = (double)DS->Node[DS->Ele[i].node[1] - 1].zmax;
            c_zmax = (double)DS->Node[DS->Ele[i].node[2] - 1].zmax;

            vector1[0] = a_x - c_x;
            vector1[1] = a_y - c_y;
            vector1[2] = a_zmax - c_zmax;
            vector2[0] = b_x - c_x;
            vector2[1] = b_y - c_y;
            vector2[2] = b_zmax - c_zmax;

            normal_vector[0] =
               vector1[1] * vector2[2] - vector1[2] * vector2[1];
            normal_vector[1] =
               vector1[2] * vector2[0] - vector1[0] * vector2[2];
            normal_vector[2] =
               vector1[0] * vector2[1] - vector1[1] * vector2[0];

            if (normal_vector[2] < 0)
            {
                normal_vector[0] = -normal_vector[0];
                normal_vector[1] = -normal_vector[1];
                normal_vector[2] = -normal_vector[2];
            }

            c = sqrt (normal_vector[0] * normal_vector[0] +
               normal_vector[1] * normal_vector[1]);
            NOAH->SLOPE = atan (c / normal_vector[2]) * 180. / PI;
            ce = normal_vector[0] / c;
            se = normal_vector[1] / c;
            NOAH->ASPECT = acos (ce) * 180. / PI;
            if (se < 0)
                NOAH->ASPECT = 360. - NOAH->ASPECT;
            NOAH->ASPECT = mod (360. - NOAH->ASPECT + 270., 360.);
            NOAH->SVF = 0.;
            for (j = 0; j < 36; j++)
                NOAH->H_PHI[j] = 90.;

            for (j = 0; j < DS->NumEle; j++)
            {
                for (k = 0; k < 3; k++)
                {
                    switch (k)
                    {
                        case 0:
                            nodes[0] = 1;
                            nodes[1] = 2;
                            break;
                        case 1:
                            nodes[0] = 0;
                            nodes[1] = 2;
                            break;
                        case 2:
                            nodes[0] = 0;
                            nodes[1] = 1;
                            break;
                    }
                    x1 = (double)DS->Node[DS->Ele[j].node[nodes[0]] - 1].x;
                    y1 = (double)DS->Node[DS->Ele[j].node[nodes[0]] - 1].y;
                    z1 = (double)DS->Node[DS->Ele[j].node[nodes[0]] - 1].zmax;
                    x2 = (double)DS->Node[DS->Ele[j].node[nodes[1]] - 1].x;
                    y2 = (double)DS->Node[DS->Ele[j].node[nodes[1]] - 1].y;
                    z2 = (double)DS->Node[DS->Ele[j].node[nodes[1]] - 1].zmax;

                    xc = 0.5 * (x1 + x2);
                    yc = 0.5 * (y1 + y2);
                    zc = 0.5 * (z1 + z2);

                    vector[0] = xc - DS->Ele[i].x;
                    vector[1] = yc - DS->Ele[i].y;
                    vector[2] = zc - DS->Ele[i].zmax;
                    c = sqrt (vector[0] * vector[0] + vector[1] * vector[1]);
                    H = atan (c / vector[2]) * 180. / PI;
                    H = H < 0. ? 90. : H;

                    vector1[0] = x1 - (double)DS->Ele[i].x;
                    vector1[1] = y1 - (double)DS->Ele[i].y;
                    vector1[2] = z1 - (double)DS->Ele[i].zmax;
                    vector2[0] = x2 - (double)DS->Ele[i].x;
                    vector2[1] = y2 - (double)DS->Ele[i].y;
                    vector2[2] = z2 - (double)DS->Ele[i].zmax;

                    c1 =
                       sqrt (vector1[0] * vector1[0] +
                       vector1[1] * vector1[1]);
                    c2 =
                       sqrt (vector2[0] * vector2[0] +
                       vector2[1] * vector2[1]);
                    ce1 = vector1[0] / c1;
                    se1 = vector1[1] / c1;
                    phi1 = acos (ce1) * 180. / PI;
                    if (se1 < 0)
                        phi1 = 360. - phi1;
                    phi1 = mod (360. - phi1 + 270., 360.);

                    ce2 = vector2[0] / c2;
                    se2 = vector2[1] / c2;
                    phi2 = acos (ce2) * 180. / PI;
                    if (se2 < 0)
                        phi2 = 360. - phi2;
                    phi2 = mod (360. - phi2 + 270., 360.);

                    if (fabs (phi1 - phi2) > 180.)
                    {
                        ind1 = 0.;
                        ind2 = (int)floor ((phi1 < phi2 ? phi1 : phi2) / 10.);
                        for (ind = ind1; ind <= ind2; ind++)
                        {
                            if (H < NOAH->H_PHI[ind])
                                NOAH->H_PHI[ind] = H;
                        }
                        ind1 = (int)floor ((phi1 > phi2 ? phi1 : phi2) / 10.);
                        ind2 = 35;
                        for (ind = ind1; ind <= ind2; ind++)
                        {
                            if (H < NOAH->H_PHI[ind])
                                NOAH->H_PHI[ind] = H;
                        }
                    }
                    else
                    {
                        ind1 = (int)floor ((phi1 < phi2 ? phi1 : phi2) / 10.);
                        ind2 = (int)floor ((phi1 > phi2 ? phi1 : phi2) / 10.);
                        for (ind = ind1; ind <= ind2; ind++)
                        {
                            if (H < NOAH->H_PHI[ind])
                                NOAH->H_PHI[ind] = H;
                        }
                    }
                }
            }

            for (ind = 0; ind < 36; ind++)
            {
                NOAH->SVF =
                   NOAH->SVF +
                   0.5 / PI * (cos (NOAH->SLOPE * PI / 180.) *
                   (pow (sin (NOAH->H_PHI[ind] * PI / 180.),
                         2)) +
                   sin (NOAH->SLOPE * PI / 180.) * cos ((ind * 10. + 5. -
                         NOAH->ASPECT) * PI / 180.) * NOAH->H_PHI[ind] * PI /
                   180. -
                   sin (NOAH->H_PHI[ind] * PI / 180.) *
                   cos (NOAH->H_PHI[ind] * PI / 180.)) * 10. / 180. * PI;
            }
            if (CS->Verbose)
            {
                printf ("Ele %d: slope = %lf, aspect = %lf, svf = %lf\t", i,
                   NOAH->SLOPE, NOAH->ASPECT, NOAH->SVF);
                for (ind = 0; ind < 36; ind++)
                    printf ("%lf\t", NOAH->H_PHI[ind]);
                printf ("\n");
            }
        }

        NOAH->SNOTIME1 = 0.0;
        NOAH->RIBB = 0.0;

        NOAH->SHEAT = BADVAL;
        NOAH->ETA_KINEMATIC = BADVAL;
        NOAH->ETA = BADVAL;
        NOAH->FDOWN = BADVAL;
        NOAH->EC = BADVAL;
        NOAH->EDIR = BADVAL;
        NOAH->ETT = BADVAL;
        NOAH->ESNOW = BADVAL;
        NOAH->DRIP = BADVAL;
        NOAH->DEW = BADVAL;
        NOAH->BETA = BADVAL;
        NOAH->T1 = BADVAL;
        NOAH->SNOWH = BADVAL;
        NOAH->SNEQV = BADVAL;
        NOAH->ETP = BADVAL;
        NOAH->SSOIL = BADVAL;
        NOAH->FLX1 = BADVAL;
        NOAH->FLX2 = BADVAL;
        NOAH->FLX3 = BADVAL;
        NOAH->SNOMLT = BADVAL;
        NOAH->SNCOVR = BADVAL;
        NOAH->RUNOFF1 = BADVAL;
        NOAH->RUNOFF2 = BADVAL;
        NOAH->RUNOFF3 = BADVAL;
        NOAH->RC = BADVAL;
        NOAH->PC = BADVAL;
        NOAH->RCS = BADVAL;
        NOAH->RCT = BADVAL;
        NOAH->RCSOIL = BADVAL;
        NOAH->SOILW = BADVAL;
        NOAH->SOILM = BADVAL;
        NOAH->Q1 = BADVAL;
        NOAH->SMCWLT = BADVAL;
        NOAH->SMCDRY = BADVAL;
        NOAH->SMCREF = BADVAL;
        NOAH->SMCMAX = BADVAL;
        NOAH->RSMIN = BADVAL;
        NOAH->NROOT = BADVAL;

        NOAH->PCPDRP = 0.;
        NOAH->DRIP = 0.;

        NOAH->DT = BADVAL;

        /*
         * TBOT: convert from degree C to K 
         */
        NOAH->TBOT = NOAH->TBOT + 273.15;
        NOAH->VEGTYP = DS->Ele[i].LC;
        NOAH->SOILTYP = DS->Ele[i].soil;
        NOAH->SNOALB = 0.75;
        NOAH->ZLVL = 3.0;
        NOAH->ZLVL_WIND = DS->windH[DS->Ele[i].WindVel - 1];
        NOAH->ISURBAN = BADVAL;
        NOAH->SHDMIN = 0.01;
        NOAH->SHDMAX = 0.96;

        NOAH->USEMONALB = 0;
        NOAH->RDLAI2D = 0;
        NOAH->IZ0TLND = 1;

        NOAH->ET = (double *)malloc (NOAH->NSOIL * sizeof (double));
        NOAH->SMAV = (double *)malloc (NOAH->NSOIL * sizeof (double));
        NOAH->RTDIS = (double *)malloc (NOAH->NSOIL * sizeof (double));

        NOAH->COSZ = BADVAL;
        NOAH->PRCPRAIN = BADVAL;
        NOAH->SOLARDIRECT = BADVAL;

        NOAH->EMISSI = 0.96;
        NOAH->ALBEDO = 0.18;
        NOAH->Z0 = 0.1;

        NOAH->Z0BRD = BADVAL;

        NOAH->CZIL = LSM->GENPRMT.CZIL_DATA;

        NOAH->CH = 1.e-4;
        NOAH->CM = 1.e-4;
    }

    /*
     * Set initial conditions for land surface variables
     */


    if (CS->init_type < 3)      /* Relaxation */
    {
        for (i = 0; i < DS->NumEle; i++)
        {
            NOAH = &(LSM->GRID[i]);
            NOAH->STC =
               (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
            NOAH->SMC =
               (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
            NOAH->SH2O =
               (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
            /*
             * T1: convert from degree C to K 
             */
            NOAH->T1 =
               (double)Interpolation (&DS->Forcing[1][DS->Ele[i].temp - 1],
               CS->StartTime) + 273.15;
            NOAH->STC[0] =
               NOAH->T1 + (NOAH->T1 -
               NOAH->TBOT) / LSM->GENPRMT.ZBOT_DATA * NOAH->SLDPTH[0] * 0.5;
            for (j = 1; j < LSM->STD_NSOIL + 1; j++)
                if (NOAH->SLDPTH[j] > 0)
                    NOAH->STC[j] =
                       NOAH->STC[j - 1] + (NOAH->T1 -
                       NOAH->TBOT) / LSM->GENPRMT.ZBOT_DATA *
                       (NOAH->SLDPTH[j - 1] + NOAH->SLDPTH[j]) * 0.5;
                else
                    NOAH->STC[j] = BADVAL;
            for (j = 0; j < LSM->STD_NSOIL + 1; j++)
            {
                if (NOAH->SLDPTH[j] > 0)
                {
                    NOAH->SMC[j] = DS->Ele[i].ThetaS;
                    NOAH->SH2O[j] = DS->Ele[i].ThetaS;
                }
                else
                {
                    NOAH->SMC[j] = BADVAL;
                    NOAH->SH2O[j] = BADVAL;
                }
            }
            NOAH->SNOWH = 0.;
            NOAH->CMC = (double)DS->EleIS[i];
            NOAH->SNEQV = (double)DS->EleSnow[i];
        }
    }
    else
    {
        /*
         * Hot start mode
         */
        fn = (char *)malloc ((strlen (filename) + 15) * sizeof (char));
        sprintf (fn, "input/%s.lsminit", filename);
        init_file = fopen (fn, "r");
        free (fn);
        if (init_file == NULL)
        {
            printf
               ("\n Fatal Error: %s.lsminit is in use of does not exist!\n",
               filename);
            exit (1);
        }
        for (i = 0; i < DS->NumEle; i++)
        {
            NOAH = &(LSM->GRID[i]);
            NOAH->STC =
               (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
            NOAH->SMC =
               (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
            NOAH->SH2O =
               (double *)malloc ((LSM->STD_NSOIL + 1) * sizeof (double));
            fscanf (init_file, "%lf %lf", &NOAH->T1, &NOAH->SNOWH);
            for (j = 0; j < LSM->STD_NSOIL + 1; j++)
                fscanf (init_file, "%lf", &NOAH->STC[j]);
            for (j = 0; j < LSM->STD_NSOIL + 1; j++)
                fscanf (init_file, "%lf", &NOAH->SMC[j]);
            for (j = 0; j < LSM->STD_NSOIL + 1; j++)
                fscanf (init_file, "%lf", &NOAH->SH2O[j]);
            NOAH->CMC = (double)DS->EleIS[i];
            NOAH->SNEQV = (double)DS->EleSnow[i];
        }
    }

    //  for(i = 0; i < DS->NumEle; i++)
    //    {
    //      printf( "%lf\t%lf", LSM->GRID[i].T1, LSM->GRID[i].SNOWH);
    //        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
    //            printf( "\t%lf", LSM->GRID[i].STC[j]);
    //        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
    //            printf( "\t%lf", LSM->GRID[i].SMC[j]);
    //        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
    //            printf( "\t%lf", LSM->GRID[i].SH2O[j]);
    //        printf( "\n");
    //    }

}

void
LSM_initialize_output (char *filename, Model_Data DS, LSM_STRUCT LSM,
   char *outputdir)
{
    FILE           *Ofile;
    char           *ofn;
    int             i, j, ensemble_mode, icounter;

    if (strstr (filename, ".") != 0)
        ensemble_mode = 1;
    else
        ensemble_mode = 0;

    if (ensemble_mode == 0)
        printf ("\nInitializing LSM output files ...\t");

    icounter = 0;
    if (LSM->PRINT_T1 > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.Tsfc", outputdir, filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_T1;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].T1);
        icounter++;
    }
    if (LSM->PRINT_STC > 0)
    {
        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
        {
            sprintf (LSM->PCtrl[icounter].name, "%s%s.TSOIL%d", outputdir,
               filename, j);
            LSM->PCtrl[icounter].Interval = LSM->PRINT_STC;
            LSM->PCtrl[icounter].NumVar = DS->NumEle;
            LSM->PCtrl[icounter].PrintVar =
               (double **)malloc (LSM->PCtrl[icounter].NumVar *
               sizeof (double *));
            for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
                LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].STC[j]);
            icounter++;
        }
    }
    if (LSM->PRINT_SMC > 0)
    {
        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
        {
            sprintf (LSM->PCtrl[icounter].name, "%s%s.SM%d", outputdir,
               filename, j);
            LSM->PCtrl[icounter].Interval = LSM->PRINT_SMC;
            LSM->PCtrl[icounter].NumVar = DS->NumEle;
            LSM->PCtrl[icounter].PrintVar =
               (double **)malloc (LSM->PCtrl[icounter].NumVar *
               sizeof (double *));
            for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
                LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].SMC[j]);
            icounter++;
        }
    }
    if (LSM->PRINT_SH2O > 0)
    {
        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
        {
            sprintf (LSM->PCtrl[icounter].name, "%s%s.SW%d", outputdir,
               filename, j);
            LSM->PCtrl[icounter].Interval = LSM->PRINT_SH2O;
            LSM->PCtrl[icounter].NumVar = DS->NumEle;
            LSM->PCtrl[icounter].PrintVar =
               (double **)malloc (LSM->PCtrl[icounter].NumVar *
               sizeof (double *));
            for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
                LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].SH2O[j]);
            icounter++;
        }
    }
    if (LSM->PRINT_SNOWH > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.snowH", outputdir,
           filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_SNOWH;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].SNOWH);
        icounter++;
    }
    if (LSM->PRINT_ALBEDO > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.Albedo", outputdir,
           filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_ALBEDO;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].ALBEDO);
        icounter++;
    }
    if (LSM->PRINT_LE > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.LE", outputdir, filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_LE;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].ETA);
        icounter++;
    }
    if (LSM->PRINT_SH > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.SH", outputdir, filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_SH;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].SHEAT);
        icounter++;
    }
    if (LSM->PRINT_G > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.G", outputdir, filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_G;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].SSOIL);
        icounter++;
    }
    if (LSM->PRINT_ETP > 0)
    {
        sprintf (LSM->PCtrl[icounter].name, "%s%s.ETP", outputdir, filename);
        LSM->PCtrl[icounter].Interval = LSM->PRINT_ETP;
        LSM->PCtrl[icounter].NumVar = DS->NumEle;
        LSM->PCtrl[icounter].PrintVar =
           (double **)malloc (LSM->PCtrl[icounter].NumVar *
           sizeof (double *));
        for (i = 0; i < LSM->PCtrl[icounter].NumVar; i++)
            LSM->PCtrl[icounter].PrintVar[i] = &(LSM->GRID[i].ETP);
        icounter++;
    }
    LSM->NPRINT = icounter;

    for (i = 0; i < LSM->NPRINT; i++)
    {
        Ofile = fopen (LSM->PCtrl[i].name, "w");
        fclose (Ofile);
        LSM->PCtrl[i].buffer =
           (double *)calloc (LSM->PCtrl[i].NumVar, sizeof (double));
    }
    printf ("done.\n");
}

void LSM_FreeData (Model_Data DS, LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;

    int             i;
    for (i = 0; i < DS->NumEle; i++)
    {
        NOAH = &(LSM->GRID[i]);
        free (NOAH->STC);
        free (NOAH->SMC);
        free (NOAH->SH2O);
        free (NOAH->SLDPTH);
        free (NOAH->ET);
        free (NOAH->SMAV);
        free (NOAH->RTDIS);
    }
    for (i = 0; i < LSM->NPRINT; i++)
    {
        free (LSM->PCtrl[i].PrintVar);
        free (LSM->PCtrl[i].buffer);
    }
    free (LSM->GRID);
    free (LSM->STD_SLDPTH);
}

void LSM_PrintInit (Model_Data DS, LSM_STRUCT LSM, char *filename)
{
    FILE           *init_file;
    char           *init_name;
    int             i, j;
    init_name = (char *)malloc ((strlen (filename) + 15) * sizeof (char));
    sprintf (init_name, "input/%s.lsminit", filename);
    init_file = fopen (init_name, "w");
    free (init_name);

    for (i = 0; i < DS->NumEle; i++)
    {
        fprintf (init_file, "%lf\t%lf", LSM->GRID[i].T1, LSM->GRID[i].SNOWH);
        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
            fprintf (init_file, "\t%lf", LSM->GRID[i].STC[j]);
        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
            fprintf (init_file, "\t%lf", LSM->GRID[i].SMC[j]);
        for (j = 0; j < LSM->STD_NSOIL + 1; j++)
            fprintf (init_file, "\t%lf", LSM->GRID[i].SH2O[j]);
        fprintf (init_file, "\n");
    }

    fclose (init_file);
}

int FindLayer (LSM_STRUCT LSM, double depth)
{
    int             layer;
    int             j = 0, ind = 0;
    double          dsum = 0;

    if (depth <= 0)
        layer = 0;
    else
    {
        while (dsum < depth)
        {
            if (LSM->STD_SLDPTH[j] < 0)
                break;
            dsum = dsum + LSM->STD_SLDPTH[j];
            ind = j;
            j++;
        }
        layer = ind + 1;
        layer = layer > LSM->STD_NSOIL + 1 ? (LSM->STD_NSOIL + 1) : layer;
    }
    return layer;
}

double mod (double a, double N)
{
    return a - N * floor (a / N);
}
