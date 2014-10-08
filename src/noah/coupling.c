/*****************************************************************************
 * File		: coupling.c                                                 *
 * Function	: Coupling between PIHM and Noah LSM                         *
 * Model	: Flux-PIHM	            	                             *
 * Version	: August, 2014 (based on PIHM 2.2)			     *
 * Developer of Flux-PIHM   :	Yuning Shi      (yshi@psu.edu)               *
 * Developer of PIHM2.0     :	Mukesh Kumar    (muk139@psu.edu)	     * 
 * Developer of PIHM1.0     :	Yizhong Qu      (quyizhong@gmail.com)	     * 
 *---------------------------------------------------------------------------*
 * For questions or comments, please contact                                 *
 *      --> Yuning Shi (yshi@psu.edu)                                        *
 * This code is free for research purpose only.                              *
 * Please provide relevant references if you use this code in your research  *
 *  work                                                                     *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "../pihm.h"
#include "noah.h"
#include "flux_pihm.h"
#include "../spa/spa.h"

void PIHM2Noah (realtype t, realtype stepsize, Model_Data DS, LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;

    double          TRESH = 0.95, A2 = 17.67, A3 = 273.15, A4 = 29.65, T0 =
       273.16, ELWV = 2.501e6, A23M4, E0 = 611.0, RV = 461.0, EPSILON = 0.622;
    double          E;
    double          SVP, SVP1 = 611.2, SVP2 = 17.67, SVP3 = 29.65, SVPT0 =
       273.15;
    double          T1V, TH2V, T2V, RHO, ES, RH;
    double          ZSOIL[LSM->STD_NSOIL + 1];
    int             i, j, KZ;
    int             spa_result;

    double          Soldown, Sdir, Sdif, gvf;
    double          incidence, azimuth180;
    time_t         *rawtime;
    struct tm      *timestamp;

    spa_data        spa;
    rawtime = (time_t *) malloc (sizeof (time_t));
    *rawtime = (int)t *60;
    timestamp = gmtime (rawtime);
    free (rawtime);

    spa.year = timestamp->tm_year + 1900;
    spa.month = timestamp->tm_mon + 1;
    spa.day = timestamp->tm_mday;
    spa.hour = timestamp->tm_hour;
    spa.minute = timestamp->tm_min;
    spa.second = timestamp->tm_sec;
    spa.timezone = 0;
    spa.delta_t = 67;
    spa.delta_ut1 = 0;
    spa.atmos_refract = 0.5667;

    spa.longitude = LSM->LONGITUDE;
    spa.latitude = LSM->LATITUDE;
    spa.elevation = 0.;
    for (i = 0; i < DS->NumEle; i++)
        spa.elevation = spa.elevation + (double)DS->Ele[i].zmax;
    spa.elevation = spa.elevation / (double)DS->NumEle;
    /*
     * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
     */
    spa.pressure =
       1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293., 5.26);
    spa.temperature = LSM->GENPRMT.TBOT_DATA;

    spa.function = SPA_ZA;
    spa_result = spa_calculate (&spa);

    if (spa_result != 0)
    {
        printf ("SPA Error Code: %d\n", spa_result);
        exit (0);
    }
    spa.azimuth180 = mod ((360. + spa.azimuth180), 360.);

    for (i = 0; i < DS->NumEle; i++)
    {
        NOAH = &(LSM->GRID[i]);

        /*
         * Read time step
         */

        NOAH->DT = (double)(stepsize * 60.);    /* DT: convert from d to s */

        /*
         * Read forcing
         */

        NOAH->SFCSPD = Interpolation (&DS->Forcing[3][DS->Ele[i].WindVel - 1], t) / 24. / 3600.;    /* SFCSPD: convert from m day-1 to m s-1 */
        NOAH->SFCTMP = Interpolation (&DS->Forcing[1][DS->Ele[i].temp - 1], t) + 273.15;    /* SFCTMP: convert from degree C to Kalvain */
        RH = Interpolation (&DS->Forcing[2][DS->Ele[i].humidity - 1], t) * 100.;    /* RH: convert from 100% to % */
        NOAH->SFCPRS = Interpolation (&DS->Forcing[6][DS->Ele[i].pressure - 1], t) / 100.;  /* SFCPRS: convert from m day-1 to m s-1 */
        NOAH->LONGWAVE = Interpolation (&DS->Forcing[5][DS->Ele[i].Ldown - 1], t) / 24. / 3600.;    /* LONGWAVE: convert from J day-1 m-2 to W m-2 */
        NOAH->PRCP = Interpolation (&DS->Forcing[0][DS->Ele[i].prep - 1], t) * 1000. / 24. / 3600.; /* PRCP: convert from m day-1 to kg m-2 s-1 */

        /*
         * Calculate solar radiation
         */
        if (LSM->RAD_MODE > 0)
        {
            Sdir =
               (double)Interpolation (&DS->Forcing[11][DS->Ele[i].Sdown - 1],
               t);
            Sdif =
               (double)Interpolation (&DS->Forcing[12][DS->Ele[i].Sdown - 1],
               t);

            if (spa.zenith > NOAH->H_PHI[(int)floor (spa.azimuth180 / 10.)])
                Sdir = 0.;
            incidence =
               180. / PI * acos (cos (spa.zenith * PI / 180.) *
               cos (NOAH->SLOPE * PI / 180.) +
               sin (spa.zenith * PI / 180.) * sin (NOAH->SLOPE * PI / 180.) *
               cos ((spa.azimuth180 - NOAH->ASPECT) * PI / 180.));
            incidence = incidence > 90. ? 90. : incidence;
            gvf = (1. + cos (NOAH->SLOPE * PI / 180.)) / 2. - NOAH->SVF;
            gvf = gvf < 0. ? 0. : gvf;
            Soldown =
               Sdir * cos (incidence * PI / 180.) + NOAH->SVF * Sdif +
               0.2 * gvf * (Sdir * cos (spa.zenith * PI / 180.) + Sdif);
        }
        else
        {
            if (DS->NumTS[11] > 0 && DS->NumTS[12] > 0)
            {
                Sdir =
                   (double)Interpolation (&DS->Forcing[11][DS->Ele[i].Sdown -
                      1], t);
                Sdif =
                   (double)Interpolation (&DS->Forcing[12][DS->Ele[i].Sdown -
                      1], t);
                Soldown = Sdir * cos (spa.zenith * PI / 180.);
                Soldown = Soldown < 0. ? 0. : Soldown;
                Soldown = Soldown + Sdif;
            }
            else
                Soldown =
                   Interpolation (&DS->Forcing[4][DS->Ele[i].Sdown - 1], t);
        }
        NOAH->SOLDN = Soldown / 24. / 3600.;    /* SOLDN: convert from J day-1 m-2 to W m-2 */

        NOAH->SFCPRS = NOAH->SFCPRS * 1.e2;

        /*
         * Initiate LSM variables
         */

        RH = RH / 100.0;

        SVP =
           SVP1 * exp (SVP2 * (NOAH->SFCTMP - SVPT0) / (NOAH->SFCTMP - SVP3));
        E = RH * SVP;

        NOAH->Q2 = (0.622 * E) / (NOAH->SFCPRS - (1.0 - 0.622) * E);

        if (NOAH->PRCP > 0 && NOAH->SFCTMP < 273.15)
            NOAH->FFROZP = 1.0;
        else
            NOAH->FFROZP = 0.0;

        NOAH->TH2 = NOAH->SFCTMP + (0.0098 * NOAH->ZLVL);
        T1V = NOAH->T1 * (1.0 + 0.61 * NOAH->Q2);
        TH2V = NOAH->TH2 * (1.0 + 0.61 * NOAH->Q2);
        T2V = NOAH->SFCTMP * (1.0 + 0.61 * NOAH->Q2);
        RHO = NOAH->SFCPRS / (RD * T2V);

        A23M4 = A2 * (A3 - A4);

        ES = E0 * exp (ELWV / RV * (1. / A3 - 1. / NOAH->SFCTMP));
        NOAH->Q2SAT = EPSILON * ES / (NOAH->SFCPRS - (1 - EPSILON) * ES);

        NOAH->DQSDT2 = NOAH->Q2SAT * A23M4 / pow (NOAH->SFCTMP - A4, 2);

        if (NOAH->USEMONALB)
            NOAH->ALB = 0.18;
        else
            NOAH->ALB = BADVAL;

        if (NOAH->RDLAI2D)
            NOAH->XLAI = 2.0;
        else
            NOAH->XLAI = BADVAL;

        NOAH->SHDFAC = DS->LandC[DS->Ele[i].LC - 1].VegFrac;

        NOAH->XLAI = Interpolation (&DS->Forcing[7][DS->Ele[i].LC - 1], t);

        if (NOAH->Q1 == BADVAL)
            NOAH->Q1 = NOAH->Q2;

        SFCDIF_off (&(NOAH->ZLVL), &(NOAH->ZLVL_WIND), &(NOAH->Z0), &T1V,
           &TH2V, &(NOAH->SFCSPD), &(NOAH->CZIL), &(NOAH->CM), &(NOAH->CH),
           &(NOAH->VEGTYP), &(NOAH->ISURBAN), &(NOAH->IZ0TLND));

        NOAH->SOLNET = NOAH->SOLDN * (1.0 - NOAH->ALBEDO);
        NOAH->LWDN = NOAH->LONGWAVE * NOAH->EMISSI;

        ZSOIL[0] = -NOAH->SLDPTH[0];
        for (KZ = 1; KZ < NOAH->NSOIL; KZ++)
            ZSOIL[KZ] = -NOAH->SLDPTH[KZ] + ZSOIL[KZ - 1];

        NOAH->RUNOFF2 = 0;
        for (j = 0; j < 3; j++)
            NOAH->RUNOFF2 = NOAH->RUNOFF2 +
               (double)(DS->FluxSub[i][j] / DS->Ele[i].area / 24. / 3600.);    /* RUNOFF2: convert from m d-1 to m s-1 */
        NOAH->INF = (double)(DS->EleViR[i] / 24. / 3600.);  /* Infiltration: convert from m d-1 to m s-1 */
        NOAH->NWTBL =
           FindLayer (LSM,
           (double)(DS->Ele[i].zmax - DS->Ele[i].zmin - DS->EleGW[i]));
        if (NOAH->NWTBL > NOAH->NSOIL)
            NOAH->NWTBL = NOAH->NSOIL;

        REDPRM (NOAH, LSM, ZSOIL);
#ifdef _DEBUG_
        if (i == CHECKELE * 1000)
        {
            printf ("\nForcing %d\n", i);
            printf ("DT (s)\t%f\n", NOAH->DT);
            printf ("SFCSPD (m s-1)\t%f\n", NOAH->SFCSPD);
            printf ("SFCTMP (K)\t%f\n", NOAH->SFCTMP);
            printf ("RH (%)\t%f\n", RH);
            printf ("SFCPRS (hPa)\t%f\n", NOAH->SFCPRS);
            printf ("SOLDN (W m-2)\t%f\n", NOAH->SOLDN);
            printf ("LONGWAVE (W m-1)\t%f\n", NOAH->LONGWAVE);
            printf ("PRCP (kg m-2 s-1)\t%f\n", NOAH->PRCP);
            printf ("SFCSPD (m-2 m-2)\t%f\n", NOAH->XLAI);
            printf ("NWTBL (m)\t%d\n", NOAH->NWTBL);
            printf ("INF (m s-1)\t%f\n", NOAH->INF);
            printf ("RUNOFF2 (m s-1)\t%f\n", NOAH->RUNOFF2);
        }
#endif

        SFLX (NOAH);

        NOAH->DRIP = 1.e3 * NOAH->DRIP / NOAH->DT;  /* Convert DRIP from m/timestep to kg m{-2} s{-1} (mm/s) */
    }
}

void Noah2PIHM (Model_Data DS, LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;
    double          ETsat;
    double          FCR, ACRT, DICE, SUM;
    int             IALP1;
    int             J, JJ, K;
    int             CVFRZ = 3;

    int             i, j;
    for (i = 0; i < DS->NumEle; i++)
    {
        NOAH = &(LSM->GRID[i]);

        DS->EleNetPrep[i] = (realtype) NOAH->PCPDRP * 24. * 3600.;  /* EleNetPrep: convert from m s-1 to m day-1 */
        /*
         * EleET: convert from W m-2 to m day-1 
         */
        DS->EleET[i][0] = (realtype) NOAH->EC / LVH2O / 1000. * 24. * 3600.;
        DS->EleET[i][1] = (realtype) NOAH->ETT / LVH2O / 1000. * 24. * 3600.;
        DS->EleET[i][2] = (realtype) NOAH->EDIR / LVH2O / 1000. * 24. * 3600.;

        /*
         * Calculate transpiration from saturated zone
         */
        DS->EleETsat[i] = 0;
        ETsat = 0;
        if (NOAH->ETT > 0)
        {
            if (NOAH->NWTBL <= NOAH->NROOT)
            {
                for (j = (NOAH->NWTBL <= 0 ? 0 : NOAH->NWTBL - 1);
                   j < NOAH->NROOT; j++)
                    ETsat = ETsat + NOAH->ET[j];
                DS->EleETsat[i] = (realtype) (ETsat / NOAH->ETT);
                DS->EleETsat[i] = DS->EleETsat[i] > 1. ? 1. : DS->EleETsat[i];
                DS->EleETsat[i] = DS->EleETsat[i] < 0. ? 0. : DS->EleETsat[i];
            }
        }
        DS->EleSnow[i] = (realtype) NOAH->SNEQV;
        DS->EleIS[i] = (realtype) NOAH->CMC;

        /*
         * Calculate surface saturation ratio for PIHM infiltration
         */
        DS->SfcSat[i] =
           (NOAH->SH2O[0] - NOAH->SMCMIN) / (NOAH->SMCMAX - NOAH->SMCMIN);

        /*
         * Calculate infiltration reduction factor due to frozen soil
         */
        DS->EleFCR[i] = 1.;
        DICE = 0.;
        for (j = 0; j < NOAH->NSOIL; j++)
            DICE = DICE + NOAH->SLDPTH[j] * (NOAH->SMC[j] - NOAH->SH2O[j]);
        FCR = 1.;
        if (DICE > 1.e-2)
        {
            ACRT = (double)CVFRZ *NOAH->FRZX / DICE;
            SUM = 1.;
            IALP1 = CVFRZ - 1;
            for (J = 1; J < IALP1 + 1; J++)
            {
                K = 1;
                for (JJ = J + 1; JJ < IALP1; JJ++)
                    K = K * JJ;
                SUM = SUM + pow (ACRT, (double)(CVFRZ - J)) / (double)K;
            }
            FCR = 1. - exp (-ACRT) * SUM;
        }
        DS->EleFCR[i] = (realtype) FCR;
    }
}
