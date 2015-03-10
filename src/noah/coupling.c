/*****************************************************************************
 * File		: coupling.c
 * Function	: Coupling between PIHM and Noah LSM
 * Version	: August, 2014
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "../pihm.h"
#include "noah.h"
//#include "flux_pihm.h"
#include "../spa/spa.h"

void PIHM2Noah (realtype t, realtype stepsize, Model_Data PIHM, LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;

    double          TRESH = 0.95, A2 = 17.67, A3 = 273.15, A4 = 29.65, T0 = 273.16, ELWV = 2.501e6, A23M4, E0 = 611.0, RV = 461.0, EPSILON = 0.622;
    double          E;
    double          SVP, SVP1 = 611.2, SVP2 = 17.67, SVP3 = 29.65, SVPT0 = 273.15;
    double          T1V, TH2V, T2V, RHO, ES, RH;
//    double          ZSOIL[LSM->STD_NSOIL + 1];
    int             i, j, KZ;
    int             spa_result;

    double          Soldown, Sdir, Sdif, gvf;
    double          incidence, azimuth180;
    time_t         *rawtime;
    struct tm      *timestamp;
    double          metarr[7];

    spa_data        spa;
    rawtime = (time_t *) malloc (sizeof (time_t));
    *rawtime = (int)t;
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
    for (i = 0; i < PIHM->NumEle; i++)
        spa.elevation = spa.elevation + (double)PIHM->Ele[i].zmax;
    spa.elevation = spa.elevation / (double)PIHM->NumEle;
    /*
     * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
     */
    spa.pressure = 1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293., 5.26);
    spa.temperature = LSM->GENPRMT.TBOT_DATA;

    spa.function = SPA_ZA;
    spa_result = spa_calculate (&spa);

    if (spa_result != 0)
    {
        printf ("SPA Error Code: %d\n", spa_result);
        exit (0);
    }
    spa.azimuth180 = mod ((360. + spa.azimuth180), 360.);

    for (i = 0; i < PIHM->NumEle; i++)
//    for (i = 0; i < 1; i++)
    {
        NOAH = &(LSM->GRID[i]);

        /* Read time step */
        NOAH->DT = (double)stepsize;    /* DT: convert from d to s */

        MultiInterpolation (&PIHM->TSD_meteo[PIHM->Ele[i].meteo - 1], t, &metarr[0], 7);
        /* Read forcing */
//        NOAH->SFCSPD = Interpolation (&PIHM->Forcing[SFCSPD_TS][PIHM->Ele[i].WindVel - 1], t); /// 24. / 3600.;    /* SFCSPD: convert from m day-1 to m s-1 */
//        NOAH->SFCTMP = Interpolation (&PIHM->Forcing[SFCTMP_TS][PIHM->Ele[i].temp - 1], t);// + 273.15;    /* SFCTMP: convert from degree C to Kalvain */
//        RH = Interpolation (&PIHM->Forcing[RH_TS][PIHM->Ele[i].humidity - 1], t);// * 100.;    /* RH: convert from 100% to % */
//        NOAH->SFCPRS = Interpolation (&PIHM->Forcing[PRES_TS][PIHM->Ele[i].pressure - 1], t) / 100.;  /* SFCPRS: convert from Pa to hPa */
//        NOAH->LONGWAVE = Interpolation (&PIHM->Forcing[LONGWAVE_TS][PIHM->Ele[i].Ldown - 1], t);// / 24. / 3600.;    /* LONGWAVE: convert from J day-1 m-2 to W m-2 */
//        NOAH->PRCP = Interpolation (&PIHM->Forcing[PRCP_TS][PIHM->Ele[i].prep - 1], t);// * 1000. / 24. / 3600.; /* PRCP: convert from m day-1 to kg m-2 s-1 */
        NOAH->SFCSPD = metarr[SFCSPD_TS];
        NOAH->SFCTMP = metarr[SFCTMP_TS];
        RH = metarr[RH_TS];
        NOAH->SFCPRS = metarr[PRES_TS] / 100.;
        NOAH->LONGWAVE = metarr[LONGWAVE_TS];
        NOAH->PRCP = metarr[PRCP_TS];

        /* Calculate solar radiation */
        if (LSM->RAD_MODE > 0)
        {
//            Sdir = (double)Interpolation (&LSM->Forcing[SOLAR_DIR_TS][PIHM->Ele[i].Sdown - 1], t);
//            Sdif = (double)Interpolation (&LSM->Forcing[SOLAR_DIF_TS][PIHM->Ele[i].Sdown - 1], t);
            MultiInterpolation (&LSM->TSD_rad[PIHM->Ele[i].meteo - 1], t, &metarr[0], 2);
            Sdir = metarr[SOLAR_DIR_TS];
            Sdif = metarr[SOLAR_DIF_TS];

            Soldown = topo_radiation (Sdir, Sdif, spa.zenith, spa.azimuth180, NOAH->SLOPE, NOAH->ASPECT, NOAH->H_PHI, NOAH->SVF);
//            if (spa.zenith > NOAH->H_PHI[(int)floor (spa.azimuth180 / 10.)])
//                Sdir = 0.;
//            incidence = 180. / PI * acos (cos (spa.zenith * PI / 180.) * cos (NOAH->SLOPE * PI / 180.) + sin (spa.zenith * PI / 180.) * sin (NOAH->SLOPE * PI / 180.) * cos ((spa.azimuth180 - NOAH->ASPECT) * PI / 180.));
//            incidence = incidence > 90. ? 90. : incidence;
//            gvf = (1. + cos (NOAH->SLOPE * PI / 180.)) / 2. - NOAH->SVF;
//            gvf = gvf < 0. ? 0. : gvf;
//            Soldown = Sdir * cos (incidence * PI / 180.) + NOAH->SVF * Sdif + 0.2 * gvf * (Sdir * cos (spa.zenith * PI / 180.) + Sdif);
        }
        else
        {
//            if (LSM->NumTS[SOLAR_DIR_TS] > 0 && LSM->NumTS[SOLAR_DIF_TS] > 0)
//            {
//                Sdir = (double)Interpolation (&LSM->Forcing[SOLAR_DIR_TS][PIHM->Ele[i].Sdown - 1], t);
//                Sdif = (double)Interpolation (&LSM->Forcing[SOLAR_DIF_TS][PIHM->Ele[i].Sdown - 1], t);
//                Soldown = Sdir * cos (spa.zenith * PI / 180.);
//                Soldown = Soldown < 0. ? 0. : Soldown;
//                Soldown = Soldown + Sdif;
//            }
//            else
//                Soldown = Interpolation (&PIHM->Forcing[SOLAR_TS][PIHM->Ele[i].Sdown - 1], t);
                Soldown = metarr[SOLAR_TS];
        }
        NOAH->SOLDN = Soldown;// / 24. / 3600.;    /* SOLDN: convert from J day-1 m-2 to W m-2 */

        NOAH->SFCPRS = NOAH->SFCPRS * 1.e2;

        /*
         * Initiate LSM variables
         */

        RH = RH / 100.0;

        SVP = SVP1 * exp (SVP2 * (NOAH->SFCTMP - SVPT0) / (NOAH->SFCTMP - SVP3));
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

        NOAH->SHDFAC = PIHM->LandC[PIHM->Ele[i].LC - 1].VegFrac;

        if (PIHM->Ele[i].LAI > 0)
        {
            NOAH->XLAI = Interpolation(&PIHM->TSD_lai[PIHM->Ele[i].LAI - 1], t);
        }
        else
            NOAH->XLAI = monthly_lai (t, PIHM->Ele[i].LC);
        //NOAH->XLAI = Interpolation (&PIHM->[LAI_TS][PIHM->Ele[i].LC - 1], t);

        if (NOAH->Q1 == BADVAL)
            NOAH->Q1 = NOAH->Q2;

        SFCDIF_off (&(NOAH->ZLVL), &(NOAH->ZLVL_WIND), &(NOAH->Z0), &T1V, &TH2V, &(NOAH->SFCSPD), &(NOAH->CZIL), &(NOAH->CM), &(NOAH->CH), &(NOAH->VEGTYP), &(NOAH->ISURBAN), &(NOAH->IZ0TLND));
//        SFCDIF_off (&(NOAH->ZLVL), &(NOAH->Z0), &(NOAH->Z0), &T1V, &TH2V, &(NOAH->SFCSPD), &(NOAH->CZIL), &(NOAH->CM), &(NOAH->CH), &(NOAH->VEGTYP), &(NOAH->ISURBAN), &(NOAH->IZ0TLND));

        NOAH->SOLNET = NOAH->SOLDN * (1.0 - NOAH->ALBEDO);
        NOAH->LWDN = NOAH->LONGWAVE * NOAH->EMISSI;

//        ZSOIL[0] = -NOAH->SLDPTH[0];
//        for (KZ = 1; KZ < NOAH->NSOIL; KZ++)
//            ZSOIL[KZ] = -NOAH->SLDPTH[KZ] + ZSOIL[KZ - 1];

        NOAH->RUNOFF2 = 0;
        for (j = 0; j < 3; j++)
            NOAH->RUNOFF2 = NOAH->RUNOFF2 + (double)(PIHM->FluxSub[i][j] / PIHM->Ele[i].area );/// 24. / 3600.);    /* RUNOFF2: convert from m d-1 to m s-1 */
        NOAH->INF = (double)(PIHM->EleViR[i] );/// 24. / 3600.);  /* Infiltration: convert from m d-1 to m s-1 */
        NOAH->NWTBL = FindLayer (LSM, (double)(PIHM->Ele[i].zmax - PIHM->Ele[i].zmin - PIHM->EleGW[i]));
        if (NOAH->NWTBL > NOAH->NSOIL)
            NOAH->NWTBL = NOAH->NSOIL;

//        REDPRM (NOAH, LSM, ZSOIL);
#ifdef _DEBUG_
        if (i == 0)
        {
            printf ("\nForcing %d\n", i);
            printf ("DT (s)\t%f\n", NOAH->DT);
            printf ("SFCSPD (m s-1)\t%f\n", NOAH->SFCSPD);
            printf ("SFCTMP (K)\t%f\n", NOAH->SFCTMP);
            printf ("RH (100%)\t%f\n", RH);
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

//        if (i == 0) printf("SOLAR = %lf, T = %lf, Q = %lf, SOIL = %lf\n", NOAH->RCS, NOAH->RCT, NOAH->RCQ, NOAH->RCSOIL);

        NOAH->DRIP = 1.e3 * NOAH->DRIP / NOAH->DT;  /* Convert DRIP from m/timestep to kg m{-2} s{-1} (mm/s) */
    }
}

void Noah2PIHM (Model_Data PIHM, LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;
    double          ETsat;
    double          FCR, ACRT, DICE, SUM;
    int             IALP1;
    int             J, JJ, K;
    int             CVFRZ = 3;

    int             i, j;
    for (i = 0; i < PIHM->NumEle; i++)
    {
        NOAH = &(LSM->GRID[i]);
	PIHM->Ele[i].temp   = (realtype) NOAH->STC[2];
        PIHM->EleNetPrep[i] = (realtype) NOAH->PCPDRP;// * 1000.;// * 24. * 3600.;  /* EleNetPrep: convert from m s-1 to m day-1 */
        /*
         * EleET: convert from W m-2 to m s-1 
         */
        PIHM->EleET[i][0] = (realtype) NOAH->EC / LVH2O / 1000. ;//* 24. * 3600.;
        PIHM->EleET[i][1] = (realtype) NOAH->ETT / LVH2O / 1000. ;//* 24. * 3600.;
        PIHM->EleET[i][2] = (realtype) NOAH->EDIR / LVH2O / 1000. ;//* 24. * 3600.;

        /* Calculate transpiration from saturated zone */
        PIHM->EleETsat[i] = 0;
        ETsat = 0;
        if (NOAH->ETT > 0)
        {
            if (NOAH->NWTBL <= NOAH->NROOT)
            {
                for (j = (NOAH->NWTBL <= 0 ? 0 : NOAH->NWTBL - 1); j < NOAH->NROOT; j++)
                    ETsat = ETsat + NOAH->ET[j];
                PIHM->EleETsat[i] = (realtype) (ETsat / NOAH->ETT);
                PIHM->EleETsat[i] = PIHM->EleETsat[i] > 1. ? 1. : PIHM->EleETsat[i];
                PIHM->EleETsat[i] = PIHM->EleETsat[i] < 0. ? 0. : PIHM->EleETsat[i];
            }
        }
        PIHM->EleSnow[i] = (realtype) NOAH->SNEQV;
        PIHM->EleIS[i] = (realtype) NOAH->CMC;

        /* Calculate surface saturation ratio for PIHM infiltration */
        PIHM->SfcSat[i] = (NOAH->SH2O[0] - NOAH->SMCMIN) / (NOAH->SMCMAX - NOAH->SMCMIN);
        PIHM->SfcSat[i] = PIHM->SfcSat[i] > 0. ? PIHM->SfcSat[i] : 0.;
        PIHM->SfcSat[i] = PIHM->SfcSat[i] < 1. ? PIHM->SfcSat[i] : 1.;

        /* Calculate infiltration reduction factor due to frozen soil */
        PIHM->EleFCR[i] = 1.;
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
        PIHM->EleFCR[i] = (realtype) FCR;
#ifdef _DEBUG_
        if (i==0)
            printf("Prcp = %lf m d-1, ET = %lf m d-1 (%lf), %lf m d-1 (%lf), %lf m d-1 (%lf)\n", PIHM->EleNetPrep[i] * 24 * 3600., PIHM->EleET[i][0]* 24. * 3600.0, NOAH->EC, PIHM->EleET[i][1]*24. * 3600.0, NOAH->ETT, PIHM->EleET[i][2]* 24. * 3600.0, NOAH->EDIR);
#endif
    }
}
