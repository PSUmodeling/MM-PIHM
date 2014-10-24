/*****************************************************************************
 * File		: pihm.c
 * Function	: Main program file
 * Developer of PIHM 2.2    :	Xuan Yu	    (xxy113@psu.edu)
 * Developer of PIHM 2.0    :	Mukesh Kumar	(muk139@psu.edu)
 * Developer of PIHM 1.0    :	Yizhong Qu	(quyizhong@gmail.com)
 * Version                  :   September, 2014
 *----------------------------------------------------------------------------
 * This code is free for research purpose only.
 * Please provide relevant references if you use this code in your research
 * 	work
 *----------------------------------------------------------------------------
 * DEVELOPMENT RELATED REFERENCES:
 * Flux-PIHM:
 *  a) Shi, Y. et al., 2013, "Development of a coupled land surface
 *      hydrologic model and evaluation at a critical zone observatory",
 *      Journal of Hydrometeorology, 14, 1401--1420.
 *  b) Shi, Y. et al., 2014, "Evaluation of the parameter sensitivity of 
 *	a coupled land surface hydrologic model. Journal of Hydrometeorology,
 *	15, 279--299.
 * PIHM V2:
 *  a) Kumar, M., 2008, "Development and Implementation of a Multiscale, 
 *	Multiprocess Hydrologic Model". PhD Thesis, Penn State University
 * PIHM V1:
 *  a) Qu, Y., 2005, "An Integrated hydrologic model for multiprocec
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU
 *  b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation
 *	for multiprocess watershed simulation". Water Resources Research
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "pihm.h"               /* Data Model and Variable Declarations */

#ifdef _FLUX_PIHM_
#include "noah/noah.h"
#endif
#ifdef _BGC_
#include "bgc/bgc.h"
#endif
/*
 * Main Function
 */
int main (int argc, char *argv[])
{
    Model_Data      mData;      /* Model Data */
    Control_Data    cData;      /* Solver Control Data */
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             flag;       /* flag to test return value */
    FILE           *iproj;      /* Project File */
    FILE           *coupling;
    int             N;          /* Problem size */
    int             i, j, k;    /* loop index */
    realtype        t;          /* simulation time */
    struct tm      *timestamp;
    time_t         *rawtime;
    realtype        NextPtr, StepSize;  /* stress period & step size */
    realtype        cvode_val;
    long int        cvode_int;
    char           *filename, *outputdir, str[11];
    char            system_cmd[1024];
#ifdef _FLUX_PIHM_
    LSM_STRUCT      LSM;
#endif

#ifdef _BGC_
    bgc_struct      BGCM;
#endif

    system ("clear");

    rawtime = (time_t *) malloc (sizeof (time_t));

    if (0 == (mkdir ("output", 0755)))
        printf (" Output directory was created.\n\n");

    /*
     * Project Input Name
     */
    if (argc != 2)
    {
        iproj = fopen ("input/projectName.txt", "r");
        if (iproj == NULL)
        {
            printf ("\tUsage ./pihm project_name\n");
            printf ("\t         OR              \n");
            printf ("\tUsage ./pihm, and have a file in the current directory named projectName.txt with the project name in it\n");
            fclose (iproj);
            exit (0);
        }
        else
        {
            filename = (char *)malloc (20 * sizeof (char));
            fscanf (iproj, "%s", filename);
            fclose (iproj);
        }
    }
    else
    {
        /*
         * get user specified file name in command line
         */
        filename = (char *)malloc ((strlen (argv[1]) + 1)* sizeof (char));
        strcpy (filename, argv[1]);
    }

    time (rawtime);
    timestamp = localtime (rawtime);

    printf ("\t\t########  #### ##     ## ##     ##\n");
    printf ("\t\t##     ##  ##  ##     ## ###   ###\n"); 
    printf ("\t\t##     ##  ##  ##     ## #### ####\n");
    printf ("\t\t########   ##  ######### ## ### ##\n");
    printf ("\t\t##         ##  ##     ## ##     ##\n");
    printf ("\t\t##         ##  ##     ## ##     ##\n"); 
    printf ("\t\t##        #### ##     ## ##     ##\n");
    printf ("\n\t    The Penn State Integrated Hydrologic Model\n");

#ifdef _FLUX_PIHM_
    printf ("\n\t    * Land surface module turned on.\n");
#endif
#ifdef _BGC_
    printf ("\n\t    * Biogeochemistry module turned on.\n");
#endif

    /*
     * Create output directory based on projectname and time
     */
    sprintf (str, "%2.2d%2.2d%2.2d%2.2d%2.2d", timestamp->tm_year + 1900 - 2000, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
    outputdir = (char *)malloc ((strlen (filename) + 20) * sizeof (char));
    sprintf (outputdir, "output/%s.%s/", filename, str);
    printf ("\nOutput directory: %s\n", outputdir);
    mkdir (outputdir, 0755);

    /*
     * Save input files into output directory
     */
    sprintf (system_cmd, "cp input/%s.para %s/%s.para.bak", filename, outputdir, filename);
    system (system_cmd);
    sprintf (system_cmd, "cp input/%s.calib %s/%s.calib.bak", filename, outputdir, filename);
    system (system_cmd);

    /*
     * Allocate memory for model data structure
     */
    mData = (Model_Data) malloc (sizeof *mData);
#ifdef _FLUX_PIHM_
    LSM = (LSM_STRUCT) malloc (sizeof *LSM);
#endif
#ifdef _BGC_
    BGCM = (bgc_struct) malloc (sizeof *BGCM);
#endif

    read_alloc (filename, mData, &cData);
#ifdef _FLUX_PIHM_
    LSM_read (filename, LSM);
#endif
#ifdef _BGC_
    BGC_read (filename, BGCM);
#endif

    //if(mData->UnsatMode ==1)
    //{    
    //}
    if (mData->UnsatMode == 2)
    {
        /* problem size */
        N = 3 * mData->NumEle + 2 * mData->NumRiv;
        mData->DummyY = (realtype *) malloc ((3 * mData->NumEle + 2 * mData->NumRiv) * sizeof (realtype));
    }
    /*
     * initial state variable depending on machine 
     */
    CV_Y = N_VNew_Serial (N);

    /*
     * initialize mode data structure 
     */
    initialize (filename, mData, &cData, CV_Y);
#ifdef _FLUX_PIHM_
    LSM_initialize (filename, mData, &cData, LSM);
#endif

    /*
     * initialize output files and structures 
     */
    initialize_output (filename, mData, &cData, outputdir);
#ifdef _FLUX_PIHM_
    LSM_initialize_output (filename, mData, LSM, outputdir);
#endif

    printf ("\n\nSolving ODE system ... \n\n");

    /*
     * allocate memory for solver 
     */
    cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        printf ("Fatal error: CVodeMalloc failed. \n");
        return (1);
    }

    flag = CVodeSetFdata (cvode_mem, mData);
    flag = CVodeSetInitStep (cvode_mem, cData.InitStep);
    flag = CVodeSetStabLimDet (cvode_mem, TRUE);
    flag = CVodeSetMaxStep (cvode_mem, cData.MaxStep);
    flag = CVodeMalloc (cvode_mem, f, cData.StartTime, CV_Y, CV_SS, cData.reltol, &cData.abstol);
    flag = CVSpgmr (cvode_mem, PREC_NONE, 0);
    //  flag = CVSpgmrSetGSType(cvode_mem, MODIFIED_GS);

    /*
     * set start time 
     */
    t = cData.StartTime;

    /*
     * start solver in loops 
     */
    for (i = 0; i < cData.NumSteps; i++)
    {
        /*
         * inner loops to next output points with ET step size control
         */
        while (t < cData.Tout[i + 1])
        {
            if (t + cData.ETStep >= cData.Tout[i + 1])
                NextPtr = cData.Tout[i + 1];
            else
                NextPtr = t + cData.ETStep;
            StepSize = NextPtr - t;

            mData->dt = StepSize;

            if ((int)t % (int)cData.ETStep == 0)
            {
#ifdef _FLUX_PIHM_
                /*
                 * calculate surface energy balance
                 */
                PIHM2Noah (t, cData.ETStep, mData, LSM, coupling);
                Noah2PIHM (mData, LSM);
#else
                /*
                 * calculate Interception Storage and ET
                 */
                is_sm_et (t, cData.ETStep, mData, CV_Y);
#endif
            }

#ifdef COUPLE_I
            t = NextPtr;
#else
            flag = CVodeSetMaxNumSteps(cvode_mem, (long int)(StepSize* 10));         /* Added to adatpt to larger time step. YS */
            flag = CVode (cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);
            flag = CVodeGetCurrentTime(cvode_mem, &cvode_val);
#endif
            *rawtime = (int)t;
            timestamp = gmtime (rawtime);
//            if ((int)*rawtime % 3600 == 0)
                printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
            summary (mData, CV_Y, t - StepSize, StepSize);
            update (t, mData);
        }

        /*
         * Print outputs 
         */
        for (j = 0; j < cData.NumPrint; j++)
            PrintData (cData.PCtrl[j], t, StepSize, cData.Ascii);
//        printf("%lf\n", cData.PCtrl[0].
#ifdef _FLUX_PIHM_
        for (j = 0; j < LSM->NPRINT; j++)
            PrintData (LSM->PCtrl[j], t, StepSize, cData.Ascii);
#endif
    }
    if (cData.Spinup)
    {
        PrintInit (mData, filename);
#ifdef _FLUX_PIHM_
        LSM_PrintInit (mData, LSM, filename);
#endif
    }

    printf ("\n Done. \n");
#ifdef COUPLE_I
    fclose (coupling);
#endif
#ifdef COUPLE_O
    fclose (coupling);
#endif
    /*
     * Free memory
     */
    N_VDestroy_Serial (CV_Y);

    /*
     * Free integrator memory
     */
    CVodeFree (&cvode_mem);

    free (outputdir);
    free (filename);
    free (rawtime);
    FreeData (mData, &cData);
#ifdef _FLUX_PIHM_
    LSM_FreeData (mData, LSM);
    free (LSM);
#endif
    free (mData);
    return 0;
}
