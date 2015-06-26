/*****************************************************************************
 * File		: coupling.c
 * Function	: Coupling between PIHM and Noah LSM
 * Version	: August, 2014
 ****************************************************************************/

#include "pihm.h"
#include "noah.h"
#include "spa.h"

void PIHM2Noah (realtype t, realtype stepsize, Model_Data PIHM, LSM_STRUCT LSM)
{
    GRID_TYPE      *NOAH;

    double          A2 = 17.67, A3 = 273.15, A4 = 29.65, ELWV = 2.501e6, A23M4, E0 = 611.0, RV = 461.0, EPSILON = 0.622;
    double          E;
    double          SVP, SVP1 = 611.2, SVP2 = 17.67, SVP3 = 29.65, SVPT0 = 273.15;
    double          T1V, TH2V, T2V, RHO, ES, RH;
    int             i, j;
    int             spa_result;

    double          Soldown, Sdir, Sdif;
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
    spa.elevation = 0.0;
    for (i = 0; i < PIHM->NumEle; i++)
        spa.elevation = spa.elevation + (double)PIHM->Ele[i].zmax;
    spa.elevation = spa.elevation / (double)PIHM->NumEle;
    /* Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) */
    spa.pressure = 1013.25 * pow ((293.0 - 0.0065 * spa.elevation) / 293.0, 5.26);
    spa.temperature = LSM->GENPRMT.TBOT_DATA;

    spa.function = SPA_ZA;
    spa_result = spa_calculate (&spa);

    if (spa_result != 0)
    {
        printf ("SPA Error Code: %d\n", spa_result);
        exit (1);
    }
    spa.azimuth180 = mod ((360.0 + spa.azimuth180), 360.0);

    for (i = 0; i < PIHM->NumEle; i++)
    {
        NOAH = &(LSM->GRID[i]);

        /* Read time step */
        NOAH->DT = (double)stepsize;

        /* Read forcing */
        MultiInterpolation (&PIHM->TSD_meteo[PIHM->Ele[i].meteo - 1], t, &metarr[0], 7);

        NOAH->SFCSPD = metarr[SFCSPD_TS];
        NOAH->SFCTMP = metarr[SFCTMP_TS];
        RH = metarr[RH_TS];
        NOAH->SFCPRS = metarr[PRES_TS] / 100.0;
        NOAH->LONGWAVE = metarr[LONGWAVE_TS];
        NOAH->PRCP = metarr[PRCP_TS];

        /* Calculate solar radiation */
        if (LSM->RAD_MODE > 0)
        {
            MultiInterpolation (&LSM->TSD_rad[PIHM->Ele[i].meteo - 1], t, &metarr[0], 2);
            Sdir = metarr[SOLAR_DIR_TS];
            Sdif = metarr[SOLAR_DIF_TS];

            Soldown = topo_radiation (Sdir, Sdif, spa.zenith, spa.azimuth180, NOAH->SLOPE, NOAH->ASPECT, NOAH->H_PHI, NOAH->SVF);
        }
        else
            Soldown = metarr[SOLAR_TS];

        NOAH->SOLDN = Soldown;

        NOAH->SFCPRS = NOAH->SFCPRS * 1.0e2;

        /* Initiate LSM variables */

        RH = RH / 100.0;

        SVP = SVP1 * exp (SVP2 * (NOAH->SFCTMP - SVPT0) / (NOAH->SFCTMP - SVP3));
        E = RH * SVP;

        NOAH->Q2 = (0.622 * E) / (NOAH->SFCPRS - (1.0 - 0.622) * E);

        if (NOAH->PRCP > 0.0 && NOAH->SFCTMP < 273.15)
            NOAH->FFROZP = 1.0;
        else
            NOAH->FFROZP = 0.0;

        NOAH->TH2 = NOAH->SFCTMP + (0.0098 * NOAH->ZLVL);
        T1V = NOAH->T1 * (1.0 + 0.61 * NOAH->Q2);
        TH2V = NOAH->TH2 * (1.0 + 0.61 * NOAH->Q2);
        T2V = NOAH->SFCTMP * (1.0 + 0.61 * NOAH->Q2);
        RHO = NOAH->SFCPRS / (RD * T2V);

        A23M4 = A2 * (A3 - A4);

        ES = E0 * exp (ELWV / RV * (1.0 / A3 - 1.0 / NOAH->SFCTMP));
        NOAH->Q2SAT = EPSILON * ES / (NOAH->SFCPRS - (1.0 - EPSILON) * ES);

        NOAH->DQSDT2 = NOAH->Q2SAT * A23M4 / pow (NOAH->SFCTMP - A4, 2);

        if (NOAH->USEMONALB)
            NOAH->ALB = 0.18;
        else
            NOAH->ALB = BADVAL;

        if (NOAH->RDLAI2D)
            NOAH->XLAI = 2.0;
        //else
        //    NOAH->XLAI = BADVAL;

        NOAH->SHDFAC = PIHM->LandC[PIHM->Ele[i].LC - 1].VegFrac;

#ifndef _BGC_
        if (PIHM->Ele[i].LAI > 0)
            NOAH->XLAI = Interpolation(&PIHM->TSD_lai[PIHM->Ele[i].LAI - 1], t);
        else
            NOAH->XLAI = monthly_lai (t, PIHM->Ele[i].LC);

        NOAH->CMCMAX = PIHM->ISFactor[PIHM->Ele[i].LC - 1] * NOAH->XLAI;
#endif
        if (NOAH->Q1 == BADVAL)
            NOAH->Q1 = NOAH->Q2;

        SFCDIF_off (&(NOAH->ZLVL), &(NOAH->ZLVL_WIND), &(NOAH->Z0), &T1V, &TH2V, &(NOAH->SFCSPD), &(NOAH->CZIL), &(NOAH->CM), &(NOAH->CH), &(NOAH->VEGTYP), &(NOAH->ISURBAN), &(NOAH->IZ0TLND));

        NOAH->SOLNET = NOAH->SOLDN * (1.0 - NOAH->ALBEDO);
        NOAH->LWDN = NOAH->LONGWAVE * NOAH->EMISSI;

        NOAH->RUNOFF2 = 0.0;
        for (j = 0; j < 3; j++)
            NOAH->RUNOFF2 = NOAH->RUNOFF2 + (double)(PIHM->avg_subflux[i][j] / PIHM->Ele[i].area );
        NOAH->INF = (double)(PIHM->avg_inf[i]);
        NOAH->NWTBL = FindLayer (LSM, (double)(PIHM->Ele[i].zmax - PIHM->Ele[i].zmin - PIHM->EleGW[i]));
        if (NOAH->NWTBL > NOAH->NSOIL)
            NOAH->NWTBL = NOAH->NSOIL;

        NOAH->MAC_STATUS = PIHM->EleMacAct[i];

        SFLX (NOAH);

        NOAH->DRIP = 1.0e3 * NOAH->DRIP / NOAH->DT;  /* Convert DRIP from m/s to kg m{-2} s{-1} (mm/s) */
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
#ifdef _RT_
	PIHM->Ele[i].temp   = (realtype) NOAH->STC[2];
#endif
        PIHM->EleNetPrep[i] = (realtype) NOAH->PCPDRP;
        /*
         * EleET: convert from W m-2 to m s-1 
         */
        PIHM->EleET[i][0] = (realtype) NOAH->EC / LVH2O / 1000.0 ;
        PIHM->EleET[i][1] = (realtype) NOAH->ETT / LVH2O / 1000.0 ;
        PIHM->EleET[i][2] = (realtype) NOAH->EDIR / LVH2O / 1000.0 ;

        /* Calculate transpiration from saturated zone */
        PIHM->EleETsat[i] = 0.0;
        ETsat = 0.0;
        if (NOAH->ETT > 0.0)
        {
            if (NOAH->NWTBL <= NOAH->NROOT)
            {
                for (j = (NOAH->NWTBL <= 0 ? 0 : NOAH->NWTBL - 1); j < NOAH->NROOT; j++)
                    ETsat = ETsat + NOAH->ET[j];
                PIHM->EleETsat[i] = (realtype) (ETsat / NOAH->ETT);
                PIHM->EleETsat[i] = PIHM->EleETsat[i] > 1.0 ? 1.0 : PIHM->EleETsat[i];
                PIHM->EleETsat[i] = PIHM->EleETsat[i] < 0.0 ? 0.0 : PIHM->EleETsat[i];
            }
        }
        PIHM->EleSnow[i] = (realtype) NOAH->SNEQV;
        PIHM->EleIS[i] = (realtype) NOAH->CMC;

        /* Calculate surface saturation ratio for PIHM infiltration */
        PIHM->SfcSat[i] = (NOAH->SH2O[0] - NOAH->SMCMIN) / (NOAH->SMCMAX - NOAH->SMCMIN);
        PIHM->SfcSat[i] = PIHM->SfcSat[i] > 0.0 ? PIHM->SfcSat[i] : 0.0;
        PIHM->SfcSat[i] = PIHM->SfcSat[i] < 1.0 ? PIHM->SfcSat[i] : 1.0;

        /* Calculate infiltration reduction factor due to frozen soil */
        PIHM->EleFCR[i] = 1.0;
        DICE = 0.0;
        for (j = 0; j < NOAH->NSOIL; j++)
            DICE = DICE + NOAH->SLDPTH[j] * (NOAH->SMC[j] - NOAH->SH2O[j]);
        FCR = 1.0;
        if (DICE > 1.0e-2)
        {
            ACRT = (double)CVFRZ *NOAH->FRZX / DICE;
            SUM = 1.0;
            IALP1 = CVFRZ - 1;
            for (J = 1; J < IALP1 + 1; J++)
            {
                K = 1;
                for (JJ = J + 1; JJ < IALP1; JJ++)
                    K = K * JJ;
                SUM = SUM + pow (ACRT, (double)(CVFRZ - J)) / (double)K;
            }
            FCR = 1.0 - exp (-ACRT) * SUM;
        }
        PIHM->EleFCR[i] = (realtype) FCR;
    }
}
