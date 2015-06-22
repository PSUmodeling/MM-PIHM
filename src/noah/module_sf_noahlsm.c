#include "noah.h"

void SFLX (GRID_TYPE * NOAH)
{
/*----------------------------------------------------------------------
* SUBROUTINE SFLX - UNIFIED NOAHLSM VERSION 1.0 JULY 2007
* ----------------------------------------------------------------------
* SUB-DRIVER FOR "Noah LSM" FAMILY OF PHYSICS SUBROUTINES FOR A
* SOIL/VEG/SNOWPACK LAND-SURFACE MODEL TO UPDATE SOIL MOISTURE, SOIL
* ICE, SOIL TEMPERATURE, SKIN TEMPERATURE, SNOWPACK WATER CONTENT,
* SNOWDEPTH, AND ALL TERMS OF THE SURFACE ENERGY BALANCE AND SURFACE
* WATER BALANCE (EXCLUDING INPUT ATMOSPHERIC FORCINGS OF DOWNWARD
* RADIATION AND PRECIP)
* --------------------------------------------------------------------*/
    int            *FRZGRA, *SNOWNG;

    int            *NSOIL, *VEGTYP;
    int            *ISURBAN;
    int            *NROOT;
    int             KZ, K, iout;

    int             RDLAI2D;
    int             USEMONALB;

    double         *SHDMIN, *SHDMAX, *DT, *DQSDT2, *LWDN, *PRCP, *PRCPRAIN,
       *Q2, *Q2SAT, *SFCPRS, *SFCSPD, *SFCTMP, *SNOALB, *SOLDN, *SOLNET,
       *TBOT, *TH2, *ZLVL, *FFROZP;
    double         *EMBRD;
    double         *ALBEDO;
    double         *COSZ, *SOLARDIRECT, *CH, *CM, *CMC, *SNEQV, *SNCOVR,
       *SNOWH, *T1, *XLAI, *SHDFAC, *Z0BRD, *EMISSI, *ALB;
    double         *SNOTIME1;
    double         *RIBB;
    double         *SLDPTH;
    double         *ET;
    double         *SMAV;
    double         *SH2O, *SMC, *STC;
    double         *RTDIS;
    double          ZSOIL[NOAH->NSOIL];

    double         *ETA_KINEMATIC, *BETA, *DEW, *DRIP, *EC, *EDIR, *ESNOW, *ETA, *ETP, *FLX1, *FLX2, *FLX3, *SHEAT, *PC, *RUNOFF2, *RUNOFF3, *RC, *RSMIN, *RCQ, *RCS, *RCSOIL, *RCT, *SSOIL, *SMCDRY, *SMCMAX, *SMCREF, *SMCWLT, *SNOMLT, *SOILM, *SOILW, *FDOWN, *Q1;
    double         *CFACTR, *CMCMAX, *CSOIL, *CZIL, *DF1, *DKSAT, *ETT, *EPSCA, *F1, *FXEXP, *FRZX, *HS, *QUARTZ, *RCH, *RR, *RGL, *RSMAX, *SNDENS, *SNCOND, *SBETA, *SN_NEW, *SNUP, *SALP, *T24, *T2V, *TOPT, *ZBOT, *Z0, *PRCPF, *ETNS, *PTU;
    double          DF1H, DF1A, DSOIL, DTOT, FRCSNO, FRCSOI, SOILWM, SOILWW;
    double         *LVCOEF;
    double          INTERP_FRACTION;
    double         *LAIMIN, *LAIMAX;
    double         *ALBEDOMIN, *ALBEDOMAX;
    double         *EMISSMIN, *EMISSMAX;
    double         *Z0MIN, *Z0MAX;

#ifdef _FLUX_PIHM_
    double         *VGALPHA, *VGBETA, *SMCMIN, *MACKSAT, *AREAF, *INF;
    double         *PCPDRP;
    int            *NMACD, *NWTBL;
    int            *MAC_STATUS;
#else
    double         *RUNOFF1, *BEXP, *DWSAT, *KDT, *PSISAT, *REFKDT, *SLOPE;
#endif

/*----------------------------------------------------------------------
*   INITIALIZATION
* --------------------------------------------------------------------*/

#ifndef _FLUX_PIHM_
    NOAH->RUNOFF1 = 0.0;
    NOAH->RUNOFF2 = 0.0;
    NOAH->RUNOFF3 = 0.0;
#endif
    NOAH->SNOMLT = 0.0;

/*----------------------------------------------------------------------
* CALCULATE DEPTH (NEGATIVE) BELOW GROUND FROM TOP SKIN SFC TO BOTTOM OF
* EACH SOIL LAYER.  NOTE:  SIGN OF ZSOIL IS NEGATIVE (DENOTING BELOW
* GROUND)
* --------------------------------------------------------------------*/

    ZSOIL[0] = -NOAH->SLDPTH[0];
    for (KZ = 1; KZ < NOAH->NSOIL; KZ++)
        ZSOIL[KZ] = -NOAH->SLDPTH[KZ] + ZSOIL[KZ - 1];

/*----------------------------------------------------------------------
* NEXT IS CRUCIAL CALL TO SET THE LAND-SURFACE PARAMETERS, INCLUDING
* SOIL-TYPE AND VEG-TYPE DEPENDENT PARAMETERS.
* --------------------------------------------------------------------*/

    //  REDPRM(NOAH, LSM, ZSOIL);           /* YS: REDPRM is now called in driver */

    FRZGRA = (int *)malloc (sizeof (int));
    SNOWNG = (int *)malloc (sizeof (int));

    NSOIL = &(NOAH->NSOIL);
#ifndef _FLUX_PIHM_
    SLOPETYP = &(NOAH->SLOPETYP);
    SOILTYP = &(NOAH->SOILTYP);
#endif
    VEGTYP = &(NOAH->VEGTYP);
    ISURBAN = &(NOAH->ISURBAN);
    NROOT = &(NOAH->NROOT);

    RDLAI2D = NOAH->RDLAI2D;
    USEMONALB = NOAH->USEMONALB;

    SHDMIN = &(NOAH->SHDMIN);
    SHDMAX = &(NOAH->SHDMAX);
    DT = &(NOAH->DT);
    DQSDT2 = &(NOAH->DQSDT2);
    LWDN = &(NOAH->LWDN);
    PRCP = &(NOAH->PRCP);
    PRCPRAIN = &(NOAH->PRCPRAIN);
    Q2 = &(NOAH->Q2);
    Q2SAT = &(NOAH->Q2SAT);
    SFCPRS = &(NOAH->SFCPRS);
    SFCSPD = &(NOAH->SFCSPD);
    SFCTMP = &(NOAH->SFCTMP);
    SNOALB = &(NOAH->SNOALB);
    SOLDN = &(NOAH->SOLDN);
    SOLNET = &(NOAH->SOLNET);
    TBOT = &(NOAH->TBOT);
    TH2 = &(NOAH->TH2);
    ZLVL = &(NOAH->ZLVL);
    FFROZP = &(NOAH->FFROZP);

    EMBRD = &(NOAH->EMBRD);
    ALBEDO = &(NOAH->ALBEDO);
    COSZ = &(NOAH->COSZ);
    SOLARDIRECT = &(NOAH->SOLARDIRECT);
    CH = &(NOAH->CH);
    CM = &(NOAH->CM);
    CMC = &(NOAH->CMC);
    SNEQV = &(NOAH->SNEQV);
    SNCOVR = &(NOAH->SNCOVR);
    SNOWH = &(NOAH->SNOWH);
    T1 = &(NOAH->T1);
    XLAI = &(NOAH->XLAI);
    SHDFAC = &(NOAH->SHDFAC);
    Z0BRD = &(NOAH->Z0BRD);
    EMISSI = &(NOAH->EMISSI);
    ALB = &(NOAH->ALB);

    SNOTIME1 = &(NOAH->SNOTIME1);

    RIBB = &(NOAH->RIBB);

    SLDPTH = &(NOAH->SLDPTH[0]);
    ET = &(NOAH->ET[0]);
    SMAV = &(NOAH->SMAV[0]);
    SH2O = &(NOAH->SH2O[0]);
    SMC = &(NOAH->SMC[0]);
    STC = &(NOAH->STC[0]);

    RTDIS = &(NOAH->RTDIS[0]);

    ETA_KINEMATIC = &(NOAH->ETA_KINEMATIC);
    BETA = &(NOAH->BETA);
    DEW = &(NOAH->DEW);
    DRIP = &(NOAH->DRIP);
    EC = &(NOAH->EC);
    EDIR = &(NOAH->EDIR);
    ESNOW = &(NOAH->ESNOW);
    ETA = &(NOAH->ETA);
    ETP = &(NOAH->ETP);
    FLX1 = &(NOAH->FLX1);
    FLX2 = &(NOAH->FLX2);
    FLX3 = &(NOAH->FLX3);
    SHEAT = &(NOAH->SHEAT);
    PC = &(NOAH->PC);
    RUNOFF2 = &(NOAH->RUNOFF2);
    RUNOFF3 = &(NOAH->RUNOFF3);
    RC = &(NOAH->RC);
    RSMIN = &(NOAH->RSMIN);
    RCQ = &(NOAH->RCQ);
    RCS = &(NOAH->RCS);
    RCSOIL = &(NOAH->RCSOIL);
    RCT = &(NOAH->RCT);
    SSOIL = &(NOAH->SSOIL);
    SMCDRY = &(NOAH->SMCDRY);
    SMCMAX = &(NOAH->SMCMAX);
    SMCREF = &(NOAH->SMCREF);
    SMCWLT = &(NOAH->SMCWLT);
    SNOMLT = &(NOAH->SNOMLT);
    SOILM = &(NOAH->SOILM);
    SOILW = &(NOAH->SOILW);
    FDOWN = &(NOAH->FDOWN);
    Q1 = &(NOAH->Q1);
#ifdef _FLUX_PIHM_
    SMCMIN = &(NOAH->SMCMIN);
    VGALPHA = &(NOAH->VGALPHA);
    VGBETA = &(NOAH->VGBETA);
    MACKSAT = &(NOAH->MACKSAT);
    AREAF = &(NOAH->AREAF);
    INF = &(NOAH->INF);
    NMACD = &(NOAH->NMACD);
    MAC_STATUS = &(NOAH->MAC_STATUS);
    NWTBL = &(NOAH->NWTBL);
    PCPDRP = &(NOAH->PCPDRP);
#else
    RUNOFF1 = &(NOAH->RUNOFF1);
    BEXP = &(NOAH->BEXP);
    DWSAT = &(NOAH->DWSAT);
    KDT = &(NOAH->KDT);
    PSISAT = &(NOAH->PSISAT);
    SLOPE = &(NOAH->SLOPE);
    REFKDT = &(NOAH->REFKDT);
#endif
    CFACTR = &(NOAH->CFACTR);
    CMCMAX = &(NOAH->CMCMAX);
    CSOIL = &(NOAH->CSOIL);
    CZIL = &(NOAH->CZIL);
    DF1 = (double *)malloc (sizeof (double));
    DKSAT = &(NOAH->DKSAT);
    ETT = &(NOAH->ETT);
    EPSCA = (double *)malloc (sizeof (double));
    F1 = &(NOAH->F1);
    FXEXP = &(NOAH->FXEXP);
    FRZX = &(NOAH->FRZX);
    HS = &(NOAH->HS);
    QUARTZ = &(NOAH->QUARTZ);
    RCH = (double *)malloc (sizeof (double));
    RR = (double *)malloc (sizeof (double));
    RGL = &(NOAH->RGL);
    RSMAX = &(NOAH->RSMAX);
    SNDENS = (double *)malloc (sizeof (double));
    SNCOND = (double *)malloc (sizeof (double));
    SBETA = &(NOAH->SBETA);
    SN_NEW = (double *)malloc (sizeof (double));
    SNUP = &(NOAH->SNUP);
    SALP = &(NOAH->SALP);
    T24 = (double *)malloc (sizeof (double));
    T2V = (double *)malloc (sizeof (double));
    TOPT = &(NOAH->TOPT);
    ZBOT = &(NOAH->ZBOT);
    Z0 = &(NOAH->Z0);
    PRCPF = (double *)malloc (sizeof (double));
    ETNS = (double *)malloc (sizeof (double));
    PTU = &(NOAH->PTU);

    LVCOEF = &(NOAH->LVCOEF);
    LAIMIN = &(NOAH->LAIMIN);
    LAIMAX = &(NOAH->LAIMAX);
    ALBEDOMIN = &(NOAH->ALBEDOMIN);
    ALBEDOMAX = &(NOAH->ALBEDOMAX);
    EMISSMIN = &(NOAH->EMISSMIN);
    EMISSMAX = &(NOAH->EMISSMAX);
    Z0MIN = &(NOAH->Z0MIN);
    Z0MAX = &(NOAH->Z0MAX);

#ifdef _FLUX_PIHM_
    *PCPDRP = 0.;
#endif


    /*
     * urban 
     */
    if (*VEGTYP == *ISURBAN)
    {
        *SHDFAC = 0.05;
        *RSMIN = 400.0;
        *SMCMAX = 0.45;
#ifdef _FLUX_PIHM_
        *SMCMIN = 0.0;
#endif
        *SMCREF = 0.42;
        *SMCWLT = 0.40;
        *SMCDRY = 0.40;
    }

#ifdef _FLUX_PIHM_

/*----------------------------------------------------------------------
* YS: FLUX-PIHM USES LAI AS A FORCING VARIABLE
* VEGETATION FRACTION IS CALCULATED FROM LAI FOLLOWING NOAH-MP
* --------------------------------------------------------------------*/

    if (*XLAI >= *LAIMAX)
    {
        *EMBRD = *EMISSMAX;
        *ALB = *ALBEDOMIN;
        *Z0BRD = *Z0MAX;
    }
    else if (*XLAI <= *LAIMIN)
    {
        *EMBRD = *EMISSMIN;
        *ALB = *ALBEDOMAX;
        *Z0BRD = *Z0MIN;
    }
    else
    {
        if (*LAIMAX > *LAIMIN)
        {
            INTERP_FRACTION = (*XLAI - *LAIMIN) / (*LAIMAX - *LAIMIN);
            /*
             * Bound INTERP_FRACTION between 0 and 1 
             */
            INTERP_FRACTION = INTERP_FRACTION < 1.0 ? INTERP_FRACTION : 1.0;
            INTERP_FRACTION = INTERP_FRACTION > 0.0 ? INTERP_FRACTION : 0.0;
            /*
             * Scale Emissivity and LAI between EMISSMIN and EMISSMAX by INTERP_FRACTION 
             */
            *EMBRD = ((1.0 - INTERP_FRACTION) * *EMISSMIN) + (INTERP_FRACTION * *EMISSMAX);
            *ALB = ((1.0 - INTERP_FRACTION) * *ALBEDOMAX) + (INTERP_FRACTION * *ALBEDOMIN);
            *Z0BRD = ((1.0 - INTERP_FRACTION) * *Z0MIN) + (INTERP_FRACTION * *Z0MAX);
        }
        else
        {
            *EMBRD = 0.5 * *EMISSMIN + 0.5 * *EMISSMAX;
            *ALB = 0.5 * *ALBEDOMIN + 0.5 * *ALBEDOMAX;
            *Z0BRD = 0.5 * *Z0MIN + 0.5 * *Z0MAX;
        }
    }

//    *SHDFAC = 1. - exp (-0.52 * (*XLAI));
      *SHDFAC = 1. - exp (-0.75 * (*XLAI));
#else
    if (*SHDFAC >= *SHDMAX)
    {
        *EMBRD = *EMISSMAX;
        if (!RDLAI2D)
            *XLAI = *LAIMAX;
        if (!USEMONALB)
            *ALB = *ALBEDOMIN;
        *Z0BRD = *Z0MAX;
    }
    else if (*SHDFAC <= *SHDMIN)
    {
        *EMBRD = *EMISSMIN;
        if (!RDLAI2D)
            *XLAI = *LAIMIN;
        if (!USEMONALB)
            *ALB = *ALBEDOMAX;
        *Z0BRD = *Z0MIN;
    }
    else
    {
        if (*SHDMAX > *SHDMIN)
        {
            INTERP_FRACTION = (*SHDFAC - *SHDMIN) / (*SHDMAX - *SHDMIN);
            /*
             * Bound INTERP_FRACTION between 0 and 1 
             */
            INTERP_FRACTION = INTERP_FRACTION < 1.0 ? INTERP_FRACTION : 1.0;
            INTERP_FRACTION = INTERP_FRACTION > 0.0 ? INTERP_FRACTION : 0.0;
            /*
             * Scale Emissivity and LAI between EMISSMIN and EMISSMAX by INTERP_FRACTION 
             */
            *EMBRD =
               ((1.0 - INTERP_FRACTION) * *EMISSMIN) +
               (INTERP_FRACTION * *EMISSMAX);
            if (!RDLAI2D)
                *XLAI =
                   ((1.0 - INTERP_FRACTION) * *LAIMIN) +
                   (INTERP_FRACTION * *LAIMAX);
            if (!USEMONALB)
                *ALB =
                   ((1.0 - INTERP_FRACTION) * *ALBEDOMAX) +
                   (INTERP_FRACTION * *ALBEDOMIN);
            *Z0BRD =
               ((1.0 - INTERP_FRACTION) * *Z0MIN) +
               (INTERP_FRACTION * *Z0MAX);
        }
        else
        {
            *EMBRD = 0.5 * *EMISSMIN + 0.5 * *EMISSMAX;
            if (!RDLAI2D)
                *XLAI = 0.5 * *LAIMIN + 0.5 * *LAIMAX;
            if (!USEMONALB)
                *ALB = 0.5 * *ALBEDOMIN + 0.5 * *ALBEDOMAX;
            *Z0BRD = 0.5 * *Z0MIN + 0.5 * *Z0MAX;
        }
    }
#endif

    /*
     * INITIALIZE PRECIPITATION LOGICALS. 
     */
    *SNOWNG = 0;
    *FRZGRA = 0;

/*----------------------------------------------------------------------
* IF INPUT SNOWPACK IS NONZERO, THEN COMPUTE SNOW DENSITY "SNDENS" AND
* SNOW THERMAL CONDUCTIVITY "SNCOND" (NOTE THAT CSNOW IS A FUNCTION
* SUBROUTINE)
* --------------------------------------------------------------------*/
    if (*SNEQV <= 1.0e-7)       /* safer IF kmh (2008/03/25) */
    {
        *SNEQV = 0.0;
        *SNDENS = 0.0;
        *SNOWH = 0.0;
        *SNCOND = 1.0;
    }
    else
    {
        *SNDENS = *SNEQV / *SNOWH;
        if (*SNDENS > 1.0)
        {
            printf ("Physical snow depth is less than snow water equiv.\n");
            exit (0);
        }
        CSNOW (SNCOND, SNDENS);
    }

#ifdef _DEBUG_
    //printf ("NSOIL = %d, ISURBAN = %d, NROOT = %d\n", *NSOIL, *ISURBAN, *NROOT);
    //printf ("RDLAI2D = %d, USEMONALB = %d\n", RDLAI2D, USEMONALB);
    //printf ("SHDMIN = %f, SHDMAX = %f, DT = %f, DQSDT2 = %f, LWDN = %f, PRCP = %f, PRCPRAIN = %f, Q2 = %f, Q2SAT = %f, SFCPRS = %f, SFCSPD = %f, SFCTMP = %f, SNOALB = %f, SOLDN = %f, SOLNET = %f, TBOT = %f, TH2, = %f, ZLVL = %f, FFROZP = %f\n", *SHDMIN, *SHDMAX, *DT, *DQSDT2, *LWDN, *PRCP, *PRCPRAIN, *Q2, *Q2SAT, *SFCPRS, *SFCSPD, *SFCTMP, *SNOALB, *SOLDN, *SOLNET, *TBOT, *TH2, *ZLVL, *FFROZP);
    //printf ("CH = %f, CM = %f, CMC = %f, SNEQV = %f, SNCOVR = %f, SNOWH = %f, T1 = %f, XLAI = %f, SHDFAC = %f, Z0BRD = %f, EMISSI = %f, ALB = %f\n", *CH, *CM, *CMC, *SNEQV, *SNCOVR, *SNOWH, *T1, *XLAI, *SHDFAC, *Z0BRD, *EMISSI, *ALB);
    //printf ("SNOTIME1 = %f\n", *SNOTIME1);
    //printf ("RIBB = %f\n", *RIBB);
    //for (KZ = 0; KZ < *NSOIL; KZ++) printf ("SLDPTH[%d] = %f ", KZ, SLDPTH[KZ]);
    //printf ("\n");
    //for (KZ = 0; KZ < *NSOIL; KZ++)
    //    printf ("SH2O[%d] = %f ", KZ, SH2O[KZ]);
    //printf ("\n");
    //for (KZ = 0; KZ < *NSOIL; KZ++)
    //    printf ("SMC[%d] = %f ", KZ, SMC[KZ]);
    //printf ("\n");
    //for (KZ = 0; KZ < *NSOIL; KZ++)
    //    printf ("STC[%d] = %f ", KZ, STC[KZ]);
    //printf ("\n");
    //for (KZ = 0; KZ < *NSOIL; KZ++)
    //    printf ("RTDIS[%d] = %f ", KZ, RTDIS[KZ]);
    //printf ("\n");
    //for (KZ = 0; KZ < *NSOIL; KZ++)
    //    printf ("ZSOIL[%d] = %f ", KZ, ZSOIL[KZ]);
    //printf ("\n");
    //printf ("VGALPHA = %f, VGBETA = %f, CFACTR = %f, CMCMAX = %f, CSOIL = %f, CZIL = %f, DF1 = %f, DKSAT = %f, FXEXP = %f, FRZX = %f, LVH2O = %f, RGL = %f, RSMAX = %f, SBETA = %f, TOPT = %f, HS = %f, ZBOT = %f, SNUP = %f, SALP = %f, SMCMAX = %f, SMCWLT = %f, SMCREF = %f, SMCDRY = %f, QUARTZ = %f, LAIMIN = %f, LAIMAX = %f, EMISSMIN = %f, EMISSIMAX = %f, ALBEDOMIN = %f, ALBEDOMAX = %f, Z0MIN = %f, Z0MAX = %f\n", *VGALPHA, *VGBETA, *CFACTR, *CMCMAX, *CSOIL, *CZIL, *DF1, *DKSAT, *FXEXP, *FRZX, LVH2O, *RGL, *RSMAX, *SBETA, *TOPT, *HS, *ZBOT, *SNUP, *SALP, *SMCMAX, *SMCWLT, *SMCREF, *SMCDRY, *QUARTZ, *LAIMIN, *LAIMAX, *EMISSMIN, *EMISSMAX, *ALBEDOMIN, *ALBEDOMAX, *Z0MIN, *Z0MAX);
#endif

/*----------------------------------------------------------------------
* DETERMINE IF IT'S PRECIPITATING AND WHAT KIND OF PRECIP IT IS.
* IF IT'S PRCPING AND THE AIR TEMP IS COLDER THAN 0 C, IT'S SNOWING!
* IF IT'S PRCPING AND THE AIR TEMP IS WARMER THAN 0 C, BUT THE GRND
* TEMP IS COLDER THAN 0 C, FREEZING RAIN IS PRESUMED TO BE FALLING.
* --------------------------------------------------------------------*/

    if (*PRCP > 0.0)

        /*
         * snow defined when fraction of frozen precip (FFROZP) > 0.5,
         * passed in from model microphysics. 
         */
    {
        if (*FFROZP > 0.5)
            *SNOWNG = 1;
        else
        {
            if (*T1 <= TFREEZ)
                *FRZGRA = 1;
        }
    }

/*----------------------------------------------------------------------
* IF EITHER PRCP FLAG IS SET, DETERMINE NEW SNOWFALL (CONVERTING PRCP
* RATE FROM KG M-2 S-1 TO A LIQUID EQUIV SNOW DEPTH IN METERS) AND ADD
* IT TO THE EXISTING SNOWPACK.
* NOTE THAT SINCE ALL PRECIP IS ADDED TO SNOWPACK, NO PRECIP INFILTRATES
* INTO THE SOIL SO THAT PRCP1 IS SET TO ZERO.
* --------------------------------------------------------------------*/
    if (*SNOWNG || *FRZGRA)
    {
        *SN_NEW = *PRCP * *DT * 0.001;
        *SNEQV = *SNEQV + *SN_NEW;
        *PRCPF = 0.0;

/*----------------------------------------------------------------------
* UPDATE SNOW DENSITY BASED ON NEW SNOWFALL, USING OLD AND NEW SNOW.
* UPDATE SNOW THERMAL CONDUCTIVITY
* --------------------------------------------------------------------*/
        SNOW_NEW (SFCTMP, SN_NEW, SNOWH, SNDENS);
        CSNOW (SNCOND, SNDENS);
    }

/*----------------------------------------------------------------------
* PRECIP IS LIQUID (RAIN), HENCE SAVE IN THE PRECIP VARIABLE THAT
* LATER CAN WHOLELY OR PARTIALLY INFILTRATE THE SOIL (ALONG WITH
* ANY CANOPY "DRIP" ADDED TO THIS LATER)
* --------------------------------------------------------------------*/
    else
        *PRCPF = *PRCP;

/*----------------------------------------------------------------------
* DETERMINE SNOWCOVER AND ALBEDO OVER LAND.
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* IF SNOW DEPTH=0, SET SNOW FRACTION=0, ALBEDO=SNOW FREE ALBEDO.
* --------------------------------------------------------------------*/
    if (*SNEQV == 0.0)
    {
        *SNCOVR = 0.0;
        *ALBEDO = *ALB;
        *EMISSI = *EMBRD;
    }
    else
    {

/*----------------------------------------------------------------------
* DETERMINE SNOW FRACTIONAL COVERAGE.
* DETERMINE SURFACE ALBEDO MODIFICATION DUE TO SNOWDEPTH STATE.
* --------------------------------------------------------------------*/
        SNFRAC (SNEQV, SNUP, SALP, SNOWH, SNCOVR);

        *SNCOVR = *SNCOVR < 0.98 ? *SNCOVR : 0.98;

        ALCALC (ALB, SNOALB, EMBRD, SHDFAC, SHDMIN, SNCOVR, T1, ALBEDO,
           EMISSI, DT, SNOWNG, SNOTIME1, LVCOEF);
    }

/*----------------------------------------------------------------------
* NEXT CALCULATE THE SUBSURFACE HEAT FLUX, WHICH FIRST REQUIRES
* CALCULATION OF THE THERMAL DIFFUSIVITY.  TREATMENT OF THE
* LATTER FOLLOWS THAT ON PAGES 148-149 FROM "HEAT TRANSFER IN
* COLD CLIMATES", BY V. J. LUNARDINI (PUBLISHED IN 1981
* BY VAN NOSTRAND REINHOLD CO.) I.E. TREATMENT OF TWO CONTIGUOUS
* "PLANE PARALLEL" MEDIUMS (NAMELY HERE THE FIRST SOIL LAYER
* AND THE SNOWPACK LAYER, IF ANY). THIS DIFFUSIVITY TREATMENT
* BEHAVES WELL FOR BOTH ZERO AND NONZERO SNOWPACK, INCLUDING THE
* LIMIT OF VERY THIN SNOWPACK.  THIS TREATMENT ALSO ELIMINATES
* THE NEED TO IMPOSE AN ARBITRARY UPPER BOUND ON SUBSURFACE
* HEAT FLUX WHEN THE SNOWPACK BECOMES EXTREMELY THIN.
* ----------------------------------------------------------------------
* FIRST CALCULATE THERMAL DIFFUSIVITY OF TOP SOIL LAYER, USING
* BOTH THE FROZEN AND LIQUID SOIL MOISTURE, FOLLOWING THE
* SOIL THERMAL DIFFUSIVITY FUNCTION OF PETERS-LIDARD ET AL.
* (1998,JAS, VOL 55, 1209-1224), WHICH REQUIRES THE SPECIFYING
* THE QUARTZ CONTENT OF THE GIVEN SOIL CLASS (SEE ROUTINE REDPRM)
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* NEXT ADD SUBSURFACE HEAT FLUX REDUCTION EFFECT FROM THE
* OVERLYING GREEN CANOPY, ADAPTED FROM SECTION 2.1.2 OF
* PETERS-LIDARD ET AL. (1997, JGR, VOL 102(D4))
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_
    TDFCND (DF1, SMC, QUARTZ, SMCMAX, SMCMIN, SH2O);
#else
    TDFCND (DF1, SMC, QUARTZ, SMCMAX, SH2O);
#endif

    /*
     * urban 
     */
    if (*VEGTYP == *ISURBAN)
        *DF1 = 3.24;

    *DF1 = *DF1 * exp (*SBETA * *SHDFAC);

    /*
     * kmh 09/03/2006
     * kmh 03/25/2008  change SNCOVR threshold to 0.97
     */
    if (*SNCOVR > 0.97)
        *DF1 = *SNCOND;

/*----------------------------------------------------------------------
* FINALLY "PLANE PARALLEL" SNOWPACK EFFECT FOLLOWING
* V.J. LINARDINI REFERENCE CITED ABOVE. NOTE THAT DTOT IS
* COMBINED DEPTH OF SNOWDEPTH AND THICKNESS OF FIRST SOIL LAYER
* --------------------------------------------------------------------*/

    DSOIL = -(0.5 * ZSOIL[0]);
    if (*SNEQV == 0.)
        *SSOIL = *DF1 * (*T1 - STC[0]) / DSOIL;
    else
    {
        DTOT = *SNOWH + DSOIL;
        FRCSNO = *SNOWH / DTOT;

        /*
         * 1. HARMONIC MEAN (SERIES FLOW) 
         */
        //      DF1 = (SNCOND*DF1)/(FRCSOI*SNCOND+FRCSNO*DF1)
        FRCSOI = DSOIL / DTOT;

        /*
         * 2. ARITHMETIC MEAN (PARALLEL FLOW) 
         */
        //      DF1 = FRCSNO*SNCOND + FRCSOI*DF1
        DF1H = (*SNCOND * *DF1) / (FRCSOI * *SNCOND + FRCSNO * *DF1);

        /*
         * 3. GEOMETRIC MEAN (INTERMEDIATE BETWEEN HARMONIC AND ARITHMETIC MEAN) 
         */
        //      DF1 = (SNCOND**FRCSNO)*(DF1**FRCSOI)
        /*
         * weigh DF by snow fraction 
         */
        //      DF1 = DF1H*SNCOVR + DF1A*(1.0-SNCOVR)
        //      DF1 = DF1H*SNCOVR + DF1*(1.0-SNCOVR)
        DF1A = FRCSNO * *SNCOND + FRCSOI * *DF1;

/*----------------------------------------------------------------------
* CALCULATE SUBSURFACE HEAT FLUX, SSOIL, FROM FINAL THERMAL DIFFUSIVITY
* OF SURFACE MEDIUMS, DF1 ABOVE, AND SKIN TEMPERATURE AND TOP
* MID-LAYER SOIL TEMPERATURE
* --------------------------------------------------------------------*/
        *DF1 = DF1A * *SNCOVR + *DF1 * (1.0 - *SNCOVR);
        *SSOIL = *DF1 * (*T1 - STC[0]) / DTOT;
    }

/*----------------------------------------------------------------------
* DETERMINE SURFACE ROUGHNESS OVER SNOWPACK USING SNOW CONDITION FROM
* THE PREVIOUS TIMESTEP.
* --------------------------------------------------------------------*/
    if (*SNCOVR > 0.)
        SNOWZ0 (SNCOVR, Z0, Z0BRD, SNOWH);
    else
        *Z0 = *Z0BRD;

/*----------------------------------------------------------------------
* NEXT CALL ROUTINE SFCDIF TO CALCULATE THE SFC EXCHANGE COEF (CH) FOR
* HEAT AND MOISTURE.

* NOTE !!!
* DO NOT CALL SFCDIF UNTIL AFTER ABOVE CALL TO REDPRM, IN CASE
* ALTERNATIVE VALUES OF ROUGHNESS LENGTH (Z0) AND ZILINTINKEVICH COEF
* (CZIL) ARE SET THERE VIA NAMELIST I/O.

* NOTE !!!
* ROUTINE SFCDIF RETURNS A CH THAT REPRESENTS THE WIND SPD TIMES THE
* "ORIGINAL" NONDIMENSIONAL "Ch" TYPICAL IN LITERATURE.  HENCE THE CH
* RETURNED FROM SFCDIF HAS UNITS OF M/S.  THE IMPORTANT COMPANION
* COEFFICIENT OF CH, CARRIED HERE AS "RCH", IS THE CH FROM SFCDIF TIMES
* AIR DENSITY AND PARAMETER "CP".  "RCH" IS COMPUTED IN "CALL PENMAN".
* RCH RATHER THAN CH IS THE COEFF USUALLY INVOKED LATER IN EQNS.

* NOTE !!!
*-----------------------------------------------------------------------
* SFCDIF ALSO RETURNS THE SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM, CM,
* ALSO KNOWN AS THE SURFACE DRAGE COEFFICIENT. Needed as a state variable
* for iterative/implicit solution of CH in SFCDIF
* --------------------------------------------------------------------*/

    /*
     * if (!LCH)
     * {
     * *T1V = *T1 * (1.0 + 0.61 * *Q2);
     * *TH2V = *TH2 * (1.0 + 0.61 * *Q2);
     * SFCDIF_off(ZLVL, Z0, T1V, TH2V, SFCSPD, CZIL, CM, CH);
     * }
     */

/*----------------------------------------------------------------------
* CALL PENMAN SUBROUTINE TO CALCULATE POTENTIAL EVAPORATION (ETP), AND
* OTHER PARTIAL PRODUCTS AND SUMS SAVE IN COMMON/RITE FOR LATER
* CALCULATIONS.
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* CALCULATE TOTAL DOWNWARD RADIATION (SOLAR PLUS LONGWAVE) NEEDED IN
* PENMAN EP SUBROUTINE THAT FOLLOWS
* --------------------------------------------------------------------*/
    //  *FDOWN = *SOLDN * (1.0- *ALBEDO) + *LWDN;
    *FDOWN = *SOLNET + *LWDN;

/*----------------------------------------------------------------------
* CALC VIRTUAL TEMPS AND VIRTUAL POTENTIAL TEMPS NEEDED BY SUBROUTINES
* PENMAN.
* --------------------------------------------------------------------*/
    *T2V = *SFCTMP * (1.0 + 0.61 * *Q2);

#ifdef _DEBUG_
    iout = 1;
#else
    iout = 0;
#endif
    //if (iout == 1)
    //{
    //    printf ("before penman\n");
    //    printf ("SFCTMP = %lf SFCPRS = %lf CH = %lf T2V = %lf TH2 = %lf PRCP = %lf FDOWN = %lf T24 = %lf SSOIL = %lf Q2 = %f Q2SAT = %lf ETP = %lf RCH = %lf EPSCA = %lf RR = %lf SNOWNG = %d FRZGRA = %d DQSDT2 = %lf FLX2 = %lf SNOWH = %lf SNEQV = %lf DSOIL = %lf FRCSNO = %lf SNCOVR = %lf DTOT = %lf ZSOIL(1) = %lf DF1 = %lf T1 = %lf STC1 = %lf ALBEDO = %lf SMC = %lf STC = %lf SH2O = %lf\n", *SFCTMP, *SFCPRS, *CH, *T2V, *TH2, *PRCP, *FDOWN, *T24, *SSOIL, *Q2, *Q2SAT, *ETP, *RCH, *EPSCA, *RR, *SNOWNG, *FRZGRA, *DQSDT2, *FLX2, *SNOWH, *SNEQV, DSOIL, FRCSNO, *SNCOVR, DTOT, ZSOIL[0], *DF1, *T1, STC[0], *ALBEDO, SMC[0], STC[1], SH2O[0]);
    //    //      for (K = 0; K < *NSOIL; K++)
    //    //          printf("SH2O = %f, SMC = %f\t", SH2O[K], SMC[K]);
    //    //      printf("\n");
    //}

    PENMAN (SFCTMP, SFCPRS, CH, T2V, TH2, PRCP, FDOWN, T24, SSOIL, Q2, Q2SAT, ETP, RCH, EPSCA, RR, SNOWNG, FRZGRA, DQSDT2, FLX2, EMISSI, SNEQV, T1, SNCOVR);

/*----------------------------------------------------------------------
* CALL CANRES TO CALCULATE THE CANOPY RESISTANCE AND CONVERT IT INTO PC
* IF NONZERO GREENNESS FRACTION
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
*  FROZEN GROUND EXTENSION: TOTAL SOIL WATER "SMC" WAS REPLACED
*  BY UNFROZEN SOIL WATER "SH2O" IN CALL TO CANRES BELOW
* --------------------------------------------------------------------*/
    if (*SHDFAC > 0.)
        CANRES (SOLDN, CH, SFCTMP, Q2, SFCPRS, SH2O, ZSOIL, NSOIL, SMCWLT,
           SMCREF, RSMIN, RC, PC, NROOT, Q2SAT, DQSDT2, TOPT, RSMAX, RGL, HS,
           XLAI, RCS, RCT, RCQ, RCSOIL, EMISSI);
    else
        *RC = 0.0;

/*----------------------------------------------------------------------
* NOW DECIDE MAJOR PATHWAY BRANCH TO TAKE DEPENDING ON WHETHER SNOWPACK
* EXISTS OR NOT:
* --------------------------------------------------------------------*/
    *ESNOW = 0.0;
    if (*SNEQV == 0.0)
    {
#ifdef _FLUX_PIHM_
        NOPAC (ETP, ETA, PRCP, PCPDRP, SMC, SMCMAX, SMCMIN, SMCWLT, SMCREF,
           SMCDRY, CMC, CMCMAX, NSOIL, DT, SHDFAC, SBETA, Q2, T1, SFCTMP, T24,
           TH2, FDOWN, F1, EMISSI, SSOIL, STC, EPSCA, VGALPHA, VGBETA,
           MACKSAT, AREAF, NMACD, MAC_STATUS, NWTBL, PC, RCH, RR, CFACTR, SH2O, FRZX,
           ZSOIL, DKSAT, TBOT, ZBOT, INF, RUNOFF2, RUNOFF3, EDIR, EC, ET, ETT,
           NROOT, RTDIS, QUARTZ, FXEXP, CSOIL, BETA, DRIP, DEW, FLX1, FLX3,
           VEGTYP, ISURBAN);
#else
        NOPAC (ETP, ETA, PRCP, SMC, SMCMAX, SMCWLT, SMCREF, SMCDRY, CMC,
           CMCMAX, NSOIL, DT, SHDFAC, SBETA, Q2, T1, SFCTMP, T24, TH2, FDOWN,
           F1, EMISSI, SSOIL, STC, EPSCA, BEXP, PC, RCH, RR, CFACTR, SH2O,
           SLOPE, KDT, FRZX, PSISAT, ZSOIL, DKSAT, DWSAT, TBOT, ZBOT, RUNOFF1,
           RUNOFF2, RUNOFF3, EDIR, EC, ET, ETT, NROOT, RTDIS, QUARTZ, FXEXP,
           CSOIL, BETA, DRIP, DEW, FLX1, FLX3, VEGTYP, ISURBAN);
#endif
        *ETA_KINEMATIC = *ETA;
    }
    else
    {
#ifdef _FLUX_PIHM_
        SNOPAC (ETP, ETA, PRCP, PRCPF, PCPDRP, SNOWNG, SMC, SMCMAX, SMCMIN,
           SMCWLT, SMCREF, SMCDRY, CMC, CMCMAX, NSOIL, DT, SBETA, DF1, Q2, T1,
           SFCTMP, T24, TH2, FDOWN, F1, SSOIL, STC, EPSCA, SFCPRS, VGALPHA,
           VGBETA, MACKSAT, AREAF, NMACD, MAC_STATUS, NWTBL, PC, RCH, RR, CFACTR, SNCOVR,
           SNEQV, SNDENS, SNOWH, SH2O, FRZX, ZSOIL, DKSAT, TBOT, ZBOT, SHDFAC,
           INF, RUNOFF2, RUNOFF3, EDIR, EC, ET, ETT, NROOT, SNOMLT, RTDIS,
           QUARTZ, FXEXP, CSOIL, BETA, DRIP, DEW, FLX1, FLX2, FLX3, ESNOW,
           ETNS, EMISSI, RIBB, SOLDN, ISURBAN, VEGTYP);
#else
        SNOPAC (ETP, ETA, PRCP, PRCPF, SNOWNG, SMC, SMCMAX, SMCWLT, SMCREF,
           SMCDRY, CMC, CMCMAX, NSOIL, DT, SBETA, DF1, Q2, T1, SFCTMP, T24,
           TH2, FDOWN, F1, SSOIL, STC, EPSCA, SFCPRS, BEXP, PC, RCH, RR,
           CFACTR, SNCOVR, SNEQV, SNDENS, SNOWH, SH2O, SLOPE, KDT, FRZX,
           PSISAT, ZSOIL, DWSAT, DKSAT, TBOT, ZBOT, SHDFAC, RUNOFF1, RUNOFF2,
           RUNOFF3, EDIR, EC, ET, ETT, NROOT, SNOMLT, RTDIS, QUARTZ, FXEXP,
           CSOIL, BETA, DRIP, DEW, FLX1, FLX2, FLX3, ESNOW, ETNS, EMISSI,
           RIBB, SOLDN, ISURBAN, VEGTYP);
#endif
        *ETA_KINEMATIC = *ESNOW + *ETNS;
    }

    /*
     * Calculate effective mixing ratio at grnd level (skin) 
     */

    //  *Q1 = *Q2 + *ETA * CP / *RCH;
    *Q1 = *Q2 + *ETA_KINEMATIC * CP / *RCH;

/*----------------------------------------------------------------------
* DETERMINE SENSIBLE HEAT (H) IN ENERGY UNITS (W M-2)
* --------------------------------------------------------------------*/

    *SHEAT = -(*CH * CP * *SFCPRS) / (R * *T2V) * (*TH2 - *T1);

/*----------------------------------------------------------------------
* CONVERT EVAP TERMS FROM KINEMATIC (KG M-2 S-1) TO ENERGY UNITS (W M-2)
* --------------------------------------------------------------------*/
    *EDIR = *EDIR * LVH2O;
    *EC = *EC * LVH2O;
    for (K = 0; K < *NSOIL; K++)
        ET[K] = ET[K] * LVH2O;
    *ETT = *ETT * LVH2O;
    *ESNOW = *ESNOW * LSUBS;
    *ETP = *ETP * ((1. - *SNCOVR) * LVH2O + *SNCOVR * LSUBS);
    if (*ETP > 0.)
        *ETA = *EDIR + *EC + *ETT + *ESNOW;
    else
        *ETA = *ETP;

/*----------------------------------------------------------------------
* DETERMINE BETA (RATIO OF ACTUAL TO POTENTIAL EVAP)
* --------------------------------------------------------------------*/
    if (*ETP == 0.0)
        *BETA = 0.0;
    else
        *BETA = *ETA / *ETP;

/*----------------------------------------------------------------------
* CONVERT THE SIGN OF SOIL HEAT FLUX SO THAT:
*   SSOIL>0: WARM THE SURFACE  (NIGHT TIME)
*   SSOIL<0: COOL THE SURFACE  (DAY TIME)
* --------------------------------------------------------------------*/
    *SSOIL = -1.0 * *SSOIL;

/*----------------------------------------------------------------------
*  FOR THE CASE OF LAND:
*  CONVERT RUNOFF3 (INTERNAL LAYER RUNOFF FROM SUPERSAT) FROM M TO M S-1
*  AND ADD TO SUBSURFACE RUNOFF/DRAINAGE/BASEFLOW.  RUNOFF2 IS ALREADY
*  A RATE AT THIS POINT
* --------------------------------------------------------------------*/

    *RUNOFF3 = *RUNOFF3 / *DT;
#ifndef _FLUX_PIHM_
    *RUNOFF2 = *RUNOFF2 + *RUNOFF3;
#endif
    *SOILM = -1.0 * SMC[0] * ZSOIL[0];
    for (K = 1; K < *NSOIL; K++)
        *SOILM = *SOILM + SMC[K] * (ZSOIL[K - 1] - ZSOIL[K]);
    SOILWM = -1.0 * (*SMCMAX - *SMCWLT) * ZSOIL[0];
    SOILWW = -1.0 * (SMC[0] - *SMCWLT) * ZSOIL[0];

    for (K = 0; K < *NSOIL; K++)
        SMAV[K] = (SMC[K] - *SMCWLT) / (*SMCMAX - *SMCWLT);

    if (*NROOT > 0)
    {
        for (K = 1; K < *NROOT; K++)
        {
            SOILWM = SOILWM + (*SMCMAX - *SMCWLT) * (ZSOIL[K - 1] - ZSOIL[K]);
            SOILWW = SOILWW + (SMC[K] - *SMCWLT) * (ZSOIL[K - 1] - ZSOIL[K]);
        }
    }
    if (SOILWM < 1.e-6)
    {
        SOILWM = 0.0;
        *SOILW = 0.0;
        *SOILM = 0.0;
    }
    else
        *SOILW = SOILWW / SOILWM;

    /*
     * #ifdef _DEBUG_
     * printf("PRCP = %f, Q2 = %f, Q2SAT = %f, SFCPRS = %f, SFCSPD = %f, SFCTMP = %f, SOLNET = %f, ZLVL = %f\n", *PRCP, *Q2, *Q2SAT, *SFCPRS, *SFCSPD, *SFCTMP, *SOLNET, *ZLVL);
     * printf("CZIL = %f, CH = %f, CM = %f\n", *CZIL, *CH, *CM);
     * printf("XLAI = %f, SHDFAC = %f, Z0 = %f\n", *XLAI, *SHDFAC, *Z0);
     * printf("EC = %f, EDIR = %f, ETT = %f, ETP = %f\n", *EC, *EDIR, *ETT, *ETP);
     * printf("PC = %f\n", *PC);
     * #endif
     */
    free (FRZGRA);
    free (SNOWNG);
    free (DF1);
    free (EPSCA);
    free (RCH);
    free (RR);
    free (SNDENS);
    free (SNCOND);
    free (SN_NEW);
    free (T24);
    free (T2V);
    free (PRCPF);
    free (ETNS);

/*----------------------------------------------------------------------
  END SUBROUTINE SFLX
* --------------------------------------------------------------------*/
}

void
ALCALC (double *ALB, double *SNOALB, double *EMBRD, double *SHDFAC,
   double *SHDMIN, double *SNCOVR, double *TSNOW, double *ALBEDO,
   double *EMISSI, double *DT, int *SNOWNG, double *SNOTIME1, double *LVCOEF)
{

/*----------------------------------------------------------------------
* CALCULATE ALBEDO INCLUDING SNOW EFFECT (0 -> 1)
*   ALB     SNOWFREE ALBEDO
*   SNOALB  MAXIMUM (DEEP) SNOW ALBEDO
*   SHDFAC    AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
*   SHDMIN    MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION
*   SNCOVR  FRACTIONAL SNOW COVER
*   ALBEDO  SURFACE ALBEDO INCLUDING SNOW EFFECT
*   TSNOW   SNOW SURFACE TEMPERATURE (K)
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* SNOALB IS ARGUMENT REPRESENTING MAXIMUM ALBEDO OVER DEEP SNOW,
* AS PASSED INTO SFLX, AND ADAPTED FROM THE SATELLITE-BASED MAXIMUM
* SNOW ALBEDO FIELDS PROVIDED BY D. ROBINSON AND G. KUKLA
* (1985, JCAM, VOL 24, 402-411)
* --------------------------------------------------------------------*/
    double          SNOALB2;
    double          SNOALB1;
    double          SNACCA = 0.94, SNACCB = 0.58;

    /*
     * turn of vegetation effect 
     */
    //      ALBEDO = ALB + (1.0- (SHDFAC - SHDMIN))* SNCOVR * (SNOALB - ALB)
    //      ALBEDO = (1.0-SNCOVR)*ALB + SNCOVR*SNOALB !this is equivalent to below

    *ALBEDO = *ALB + *SNCOVR * (*SNOALB - *ALB);
    *EMISSI = *EMBRD + *SNCOVR * (EMISSI_S - *EMBRD);

    /*
     * BASE FORMULATION (DICKINSON ET AL., 1986, COGLEY ET AL., 1990) 
     */

    /*
     * IF (TSNOW.LE.263.16) THEN
     * ALBEDO=SNOALB
     * ELSE
     * IF (TSNOW.LT.273.16) THEN
     * TM=0.1*(TSNOW-263.16)
     * SNOALB1=0.5*((0.9-0.2*(TM**3))+(0.8-0.16*(TM**3)))
     * ELSE
     * SNOALB1=0.67
     * IF(SNCOVR.GT.0.95) SNOALB1= 0.6
     * SNOALB1 = ALB + SNCOVR*(SNOALB-ALB)
     * ENDIF
     * ENDIF
     * ALBEDO = ALB + SNCOVR*(SNOALB1-ALB)
     * 
     * ISBA FORMULATION (VERSEGHY, 1991; BAKER ET AL., 1990)
     * SNOALB1 = SNOALB+COEF*(0.85-SNOALB)
     * SNOALB2=SNOALB1
     * m          LSTSNW=LSTSNW+1
     * SNOTIME1 = SNOTIME1 + DT
     * IF (SNOWNG) THEN
     * SNOALB2=SNOALB
     * *m             LSTSNW=0
     * SNOTIME1 = 0.0
     * ELSE
     * IF (TSNOW.LT.273.16) THEN
     * *              SNOALB2=SNOALB-0.008*LSTSNW*DT/86400
     * *m              SNOALB2=SNOALB-0.008*SNOTIME1/86400
     * SNOALB2=(SNOALB2-0.65)*EXP(-0.05*DT/3600)+0.65
     * *              SNOALB2=(ALBEDO-0.65)*EXP(-0.01*DT/3600)+0.65
     * ELSE
     * SNOALB2=(SNOALB2-0.5)*EXP(-0.0005*DT/3600)+0.5
     * *              SNOALB2=(SNOALB-0.5)*EXP(-0.24*LSTSNW*DT/86400)+0.5
     * *m              SNOALB2=(SNOALB-0.5)*EXP(-0.24*SNOTIME1/86400)+0.5
     * ENDIF
     * ENDIF
     * 
     * *               print*,'SNOALB2',SNOALB2,'ALBEDO',ALBEDO,'DT',DT
     * ALBEDO = ALB + SNCOVR*(SNOALB2-ALB)
     * IF (ALBEDO .GT. SNOALB2) ALBEDO=SNOALB2
     * *m          LSTSNW1=LSTSNW
     * *          SNOTIME = SNOTIME1
     */

    /*
     * formulation by Livneh
     * * ----------------------------------------------------------------------
     * * SNOALB IS CONSIDERED AS THE MAXIMUM SNOW ALBEDO FOR NEW SNOW, AT
     * * A VALUE OF 85%. SNOW ALBEDO CURVE DEFAULTS ARE FROM BRAS P.263. SHOULD
     * * NOT BE CHANGED EXCEPT FOR SERIOUS PROBLEMS WITH SNOW MELT.
     * * TO IMPLEMENT ACCUMULATIN PARAMETERS, SNACCA AND SNACCB, ASSERT THAT IT
     * * IS INDEED ACCUMULATION SEASON. I.E. THAT SNOW SURFACE TEMP IS BELOW
     * * ZERO AND THE DATE FALLS BETWEEN OCTOBER AND FEBRUARY
     * * --------------------------------------------------------------------
     */
    SNOALB1 = *SNOALB + *LVCOEF * (0.85 - *SNOALB);
    SNOALB2 = SNOALB1;

/*---------------- Initial LSTSNW ------------------------------------*/
    if (*SNOWNG)
        *SNOTIME1 = 0.0;
    else
    {
        *SNOTIME1 = *SNOTIME1 + *DT;
        //      IF (TSNOW.LT.273.16) THEN
        SNOALB2 = SNOALB1 * pow (SNACCA, pow (*SNOTIME1 / 86400.0, SNACCB));
        //      ELSE
        //          SNOALB2 =SNOALB1*(SNTHWA**((SNOTIME1/86400.0)**SNTHWB))
        //               ENDIF
    }

    SNOALB2 = SNOALB2 > *ALB ? SNOALB2 : *ALB;
    *ALBEDO = *ALB + *SNCOVR * (SNOALB2 - *ALB);
    if (*ALBEDO > SNOALB2)
        *ALBEDO = SNOALB2;

    //          IF (TSNOW.LT.273.16) THEN
    //            ALBEDO=SNOALB-0.008*DT/86400
    //          ELSE
    //            ALBEDO=(SNOALB-0.5)*EXP(-0.24*DT/86400)+0.5
    //          ENDIF

    //      IF (ALBEDO > SNOALB) ALBEDO = SNOALB

/*----------------------------------------------------------------------
  END SUBROUTINE ALCALC
* --------------------------------------------------------------------*/
}

void
CANRES (double *SOLAR, double *CH, double *SFCTMP, double *Q2, double *SFCPRS,
   double *SMC, double *ZSOIL, int *NSOIL, double *SMCWLT, double *SMCREF,
   double *RSMIN, double *RC, double *PC, int *NROOT, double *Q2SAT,
   double *DQSDT2, double *TOPT, double *RSMAX, double *RGL, double *HS,
   double *XLAI, double *RCS, double *RCT, double *RCQ, double *RCSOIL,
   double *EMISSI)
{

/*----------------------------------------------------------------------
* SUBROUTINE CANRES
* ----------------------------------------------------------------------
* CALCULATE CANOPY RESISTANCE WHICH DEPENDS ON INCOMING SOLAR RADIATION,
* AIR TEMPERATURE, ATMOSPHERIC WATER VAPOR PRESSURE DEFICIT AT THE
* LOWEST MODEL LEVEL, AND SOIL MOISTURE (PREFERABLY UNFROZEN SOIL
* MOISTURE RATHER THAN TOTAL)
* ----------------------------------------------------------------------
* SOURCE:  JARVIS (1976), NOILHAN AND PLANTON (1989, MWR), JACQUEMIN AND
* NOILHAN (1990, BLM)
* SEE ALSO:  CHEN ET AL (1996, JGR, VOL 101(D3), 7251-7268), EQNS 12-14
* AND TABLE 2 OF SEC. 3.1.2
* ----------------------------------------------------------------------
* INPUT:
*   SOLAR   INCOMING SOLAR RADIATION
*   CH      SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE
*   SFCTMP  AIR TEMPERATURE AT 1ST LEVEL ABOVE GROUND
*   Q2      AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
*   Q2SAT   SATURATION AIR HUMIDITY AT 1ST LEVEL ABOVE GROUND
*   DQSDT2  SLOPE OF SATURATION HUMIDITY FUNCTION WRT TEMP
*   SFCPRS  SURFACE PRESSURE
*   SMC     VOLUMETRIC SOIL MOISTURE
*   ZSOIL   SOIL DEPTH (NEGATIVE SIGN, AS IT IS BELOW GROUND)
*   NSOIL   NO. OF SOIL LAYERS
*   NROOT   NO. OF SOIL LAYERS IN ROOT ZONE (1.LE.NROOT.LE.NSOIL)
*   XLAI    LEAF AREA INDEX
*   SMCWLT  WILTING POINT
*   SMCREF  REFERENCE SOIL MOISTURE (WHERE SOIL WATER DEFICIT STRESS
*             SETS IN)
* RSMIN, RSMAX, TOPT, RGL, HS ARE CANOPY STRESS PARAMETERS SET IN
*   SURBOUTINE REDPRM
* OUTPUT:
*   PC  PLANT COEFFICIENT
*   RC  CANOPY RESISTANCE
* --------------------------------------------------------------------*/

    int             K;
    double          DELTA, FF, GX, RR;
    double          PART[*NSOIL];
    double          SLV = 2.501000e6;

/*----------------------------------------------------------------------
* INITIALIZE CANOPY RESISTANCE MULTIPLIER TERMS.
* --------------------------------------------------------------------*/
    *RCS = 0.0;
    *RCT = 0.0;
    *RCQ = 0.0;
    *RCSOIL = 0.0;

/*----------------------------------------------------------------------
* CONTRIBUTION DUE TO INCOMING SOLAR RADIATION
* --------------------------------------------------------------------*/
    *RC = 0.0;
    FF = 0.55 * 2.0 * *SOLAR / (*RGL * *XLAI);
    *RCS = (FF + *RSMIN / *RSMAX) / (1.0 + FF);

//    printf("SOLAR = %lf, RGL = %lf, XLAI = %lf, FF = %lf, RSMIN = %lf, RSMAX = %lf, RCS = %lf\n", *SOLAR, *RGL, *XLAI, FF, *RSMIN, *RSMAX, *RCS);

/*----------------------------------------------------------------------
* CONTRIBUTION DUE TO AIR TEMPERATURE AT FIRST MODEL LEVEL ABOVE GROUND
* RCT EXPRESSION FROM NOILHAN AND PLANTON (1989, MWR).
* --------------------------------------------------------------------*/
    *RCS = *RCS > 0.0001 ? *RCS : 0.0001;
    *RCT = 1.0 - 0.0016 * pow (*TOPT - *SFCTMP, 2.0);

/*----------------------------------------------------------------------
* CONTRIBUTION DUE TO VAPOR PRESSURE DEFICIT AT FIRST MODEL LEVEL.
* RCQ EXPRESSION FROM SSIB
* --------------------------------------------------------------------*/
    *RCT = *RCT > 0.0001 ? *RCT : 0.0001;
    *RCQ = 1.0 / (1.0 + *HS * (*Q2SAT - *Q2));

/*----------------------------------------------------------------------
* CONTRIBUTION DUE TO SOIL MOISTURE AVAILABILITY.
* DETERMINE CONTRIBUTION FROM EACH SOIL LAYER, THEN ADD THEM UP.
* --------------------------------------------------------------------*/
    *RCQ = *RCQ > 0.01 ? *RCQ : 0.01;
    GX = (SMC[0] - *SMCWLT) / (*SMCREF - *SMCWLT);
    if (GX > 1.)
        GX = 1.;
    if (GX < 0.)
        GX = 0.;

/*----------------------------------------------------------------------
* USE SOIL DEPTH AS WEIGHTING FACTOR
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* USE ROOT DISTRIBUTION AS WEIGHTING FACTOR
*      PART(1) = RTDIS(1) * GX
* --------------------------------------------------------------------*/
    PART[0] = (ZSOIL[0] / ZSOIL[*NROOT - 1]) * GX;
    for (K = 1; K < *NROOT; K++)
    {
        GX = (SMC[K] - *SMCWLT) / (*SMCREF - *SMCWLT);
        if (GX > 1.)
            GX = 1.;
        if (GX < 0.)
            GX = 0.;

/*----------------------------------------------------------------------
* USE SOIL DEPTH AS WEIGHTING FACTOR
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* USE ROOT DISTRIBUTION AS WEIGHTING FACTOR
*        PART(K) = RTDIS(K) * GX
* --------------------------------------------------------------------*/
        PART[K] = ((ZSOIL[K] - ZSOIL[K - 1]) / ZSOIL[*NROOT - 1]) * GX;
    }
    for (K = 0; K < *NROOT; K++)
        *RCSOIL = *RCSOIL + PART[K];

/*----------------------------------------------------------------------
* DETERMINE CANOPY RESISTANCE DUE TO ALL FACTORS.  CONVERT CANOPY
* RESISTANCE (RC) TO PLANT COEFFICIENT (PC) TO BE USED WITH POTENTIAL
* EVAP IN DETERMINING ACTUAL EVAP.  PC IS DETERMINED BY:
*   PC * LINERIZED PENMAN POTENTIAL EVAP =
*   PENMAN-MONTEITH ACTUAL EVAPORATION (CONTAINING RC TERM).
* --------------------------------------------------------------------*/
    *RCSOIL = *RCSOIL > 0.0001 ? *RCSOIL : 0.0001;

    *RC = *RSMIN / (*XLAI * *RCS * *RCT * *RCQ * *RCSOIL);
    //  RR = (4.* SIGMA * RD / CP)* (SFCTMP **4.)/ (SFCPRS * CH) + 1.0;
    RR =
       (4. * *EMISSI * SIGMA * RD / CP) * pow (*SFCTMP,
       4.) / (*SFCPRS * *CH) + 1.0;

    DELTA = (SLV / CP) * *DQSDT2;

    *PC = (RR + DELTA) / (RR * (1. + *RC * *CH) + DELTA);

/*----------------------------------------------------------------------
  END SUBROUTINE CANRES
* --------------------------------------------------------------------*/
}

void CSNOW (double *SNCOND, double *DSNOW)
{

/*----------------------------------------------------------------------
* SUBROUTINE CSNOW
* FUNCTION CSNOW
* ----------------------------------------------------------------------
* CALCULATE SNOW THERMAL CONDUCTIVITY
* --------------------------------------------------------------------*/
    double          C;
    double          UNIT = 0.11631;

/*----------------------------------------------------------------------
* SNCOND IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C)
* CSNOW IN UNITS OF CAL/(CM*HR*C), RETURNED IN W/(M*C)
* BASIC VERSION IS DYACHKOVA EQUATION (1960), FOR RANGE 0.1-0.4
* --------------------------------------------------------------------*/
    C = 0.328 * pow (10, 2.25 * *DSNOW);
    //  CSNOW=UNIT*C

/*----------------------------------------------------------------------
* DE VAUX EQUATION (1933), IN RANGE 0.1-0.6
* ----------------------------------------------------------------------
*      SNCOND=0.0293*(1.+100.*DSNOW**2)
*      CSNOW=0.0293*(1.+100.*DSNOW**2)

* ----------------------------------------------------------------------
* E. ANDERSEN FROM FLERCHINGER
* ----------------------------------------------------------------------
*      SNCOND=0.021+2.51*DSNOW**2
*      CSNOW=0.021+2.51*DSNOW**2

*      SNCOND = UNIT * C
* double snow thermal conductivity
*/
    *SNCOND = 2.0 * UNIT * C;

/*----------------------------------------------------------------------
  END SUBROUTINE CSNOW
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
DEVAP (double *EDIR, double *ETP1, double *SMC, double *ZSOIL, double *SHDFAC,
   double *SMCMAX, double *DKSAT, double *SMCDRY, double *SMCREF,
   double *SMCWLT, double *FXEXP)
#else
void
DEVAP (double *EDIR, double *ETP1, double *SMC, double *ZSOIL, double *SHDFAC,
   double *SMCMAX, double *BEXP, double *DKSAT, double *DWSAT, double *SMCDRY,
   double *SMCREF, double *SMCWLT, double *FXEXP)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE DEVAP
* FUNCTION DEVAP
* ----------------------------------------------------------------------
* CALCULATE DIRECT SOIL EVAPORATION
* --------------------------------------------------------------------*/
    double          FX, SRATIO;

/*----------------------------------------------------------------------
* DIRECT EVAP A FUNCTION OF RELATIVE SOIL MOISTURE AVAILABILITY, LINEAR
* WHEN FXEXP=1.
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* FX > 1 REPRESENTS DEMAND CONTROL
* FX < 1 REPRESENTS FLUX CONTROL
* --------------------------------------------------------------------*/

    SRATIO = (*SMC - *SMCDRY) / (*SMCMAX - *SMCDRY);
    if (SRATIO > 0.)
    {
        FX = pow (SRATIO, *FXEXP);
        FX = FX > 1. ? 1. : FX;
        FX = FX < 0. ? 0. : FX;
    }
    else
        FX = 0.;

/*----------------------------------------------------------------------
* ALLOW FOR THE DIRECT-EVAP-REDUCING EFFECT OF SHADE
* --------------------------------------------------------------------*/
    *EDIR = FX * (1.0 - *SHDFAC) * *ETP1;

/*----------------------------------------------------------------------
  END SUBROUTINE DEVAP
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
EVAPO (double *ETA1, double *SMC, int *NSOIL, double *CMC, double *ETP1,
   double *DT, double *ZSOIL, double *SH2O, double *SMCMAX, double *PC,
   double *SMCWLT, double *DKSAT, double *SMCREF, double *SHDFAC,
   double *CMCMAX, double *SMCDRY, double *CFACTR, double *EDIR, double *EC,
   double *ET, double *ETT, double *SFCTMP, double *Q2, int *NROOT,
   double *RTDIS, double *FXEXP)
#else
void
EVAPO (double *ETA1, double *SMC, int *NSOIL, double *CMC, double *ETP1,
   double *DT, double *ZSOIL, double *SH2O, double *SMCMAX, double *BEXP,
   double *PC, double *SMCWLT, double *DKSAT, double *DWSAT, double *SMCREF,
   double *SHDFAC, double *CMCMAX, double *SMCDRY, double *CFACTR,
   double *EDIR, double *EC, double *ET, double *ETT, double *SFCTMP,
   double *Q2, int *NROOT, double *RTDIS, double *FXEXP)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE EVAPO
* ----------------------------------------------------------------------
* CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER
* UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH
* PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
* FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND
* CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE.
* --------------------------------------------------------------------*/
    int             K;
    double          CMC2MS;

/*----------------------------------------------------------------------
* EXECUTABLE CODE BEGINS HERE IF THE POTENTIAL EVAPOTRANSPIRATION IS
* GREATER THAN ZERO.
* --------------------------------------------------------------------*/
    *EDIR = 0.;
    *EC = 0.;
    *ETT = 0.;
    for (K = 0; K < *NSOIL; K++)
        ET[K] = 0.;

/*----------------------------------------------------------------------
* RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE.  CALL THIS FUNCTION
* ONLY IF VEG COVER NOT COMPLETE.
* FROZEN GROUND VERSION:  SH2O STATES REPLACE SMC STATES.
* --------------------------------------------------------------------*/
    if (*ETP1 > 0.0)
    {
        if (*SHDFAC < 1.)
#ifdef _FLUX_PIHM_
            DEVAP (EDIR, ETP1, SMC, ZSOIL, SHDFAC, SMCMAX, DKSAT, SMCDRY,
               SMCREF, SMCWLT, FXEXP);
#else
            DEVAP (EDIR, ETP1, SMC, ZSOIL, SHDFAC, SMCMAX, BEXP, DKSAT, DWSAT,
               SMCDRY, SMCREF, SMCWLT, FXEXP);
#endif

/*----------------------------------------------------------------------
* INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION,
* AND ACCUMULATE IT FOR ALL SOIL LAYERS.
* --------------------------------------------------------------------*/

        if (*SHDFAC > 0.0)
        {
            TRANSP (ET, NSOIL, ETP1, SH2O, CMC, ZSOIL, SHDFAC, SMCWLT, CMCMAX,
               PC, CFACTR, SMCREF, SFCTMP, Q2, NROOT, RTDIS);
            for (K = 0; K < *NSOIL; K++)
                *ETT = *ETT + ET[K];

/*----------------------------------------------------------------------
* CALCULATE CANOPY EVAPORATION.
* IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR CMC=0.0.
* --------------------------------------------------------------------*/
            if (*CMC > 0.0)
                *EC =
                   *SHDFAC * pow ((*CMC / *CMCMAX > 1. ? 1. : *CMC / *CMCMAX),
                   *CFACTR) * *ETP1;
            else
                *EC = 0.0;

/*----------------------------------------------------------------------
* EC SHOULD BE LIMITED BY THE TOTAL AMOUNT OF AVAILABLE WATER ON THE
* CANOPY.  -F.CHEN, 18-OCT-1994
* --------------------------------------------------------------------*/
            CMC2MS = *CMC / *DT;
            *EC = CMC2MS < *EC ? CMC2MS : *EC;
        }
    }

/*----------------------------------------------------------------------
* TOTAL UP EVAP AND TRANSP TYPES TO OBTAIN ACTUAL EVAPOTRANSP
* --------------------------------------------------------------------*/
    *ETA1 = *EDIR + *ETT + *EC;

/*----------------------------------------------------------------------
  END SUBROUTINE EVAPO
* --------------------------------------------------------------------*/
}

void FAC2MIT (double *SMCMAX, double *FLIMIT)
{
    *FLIMIT = 0.90;

    if (*SMCMAX == 0.395)
        *FLIMIT = 0.59;
    else if ((*SMCMAX == 0.434) || (*SMCMAX == 0.404))
        *FLIMIT = 0.85;
    else if ((*SMCMAX == 0.465) || (*SMCMAX == 0.406))
        *FLIMIT = 0.86;
    else if ((*SMCMAX == 0.476) || (*SMCMAX == 0.439))
        *FLIMIT = 0.74;
    else if ((*SMCMAX == 0.200) || (*SMCMAX == 0.464))
        *FLIMIT = 0.80;

/*----------------------------------------------------------------------
  END SUBROUTINE FAC2MIT
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
FRH2O (double *FREE, double *TKELV, double *SMC, double *SH2O, double *SMCMAX,
   double *SMCMIN, double *VGALPHA, double *VGBETA)
#else
void
FRH2O (double *FREE, double *TKELV, double *SMC, double *SH2O, double *SMCMAX,
   double *BEXP, double *PSIS)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE FRH2O
* ----------------------------------------------------------------------
* CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT IF
* TEMPERATURE IS BELOW 273.15K (T0).  REQUIRES NEWTON-TYPE ITERATION TO
* SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF KOREN ET AL
* (1999, JGR, VOL 104(D16), 19569-19585).
* ----------------------------------------------------------------------
* NEW VERSION (JUNE 2001): MUCH FASTER AND MORE ACCURATE NEWTON
* ITERATION ACHIEVED BY FIRST TAKING LOG OF EQN CITED ABOVE -- LESS THAN
* 4 (TYPICALLY 1 OR 2) ITERATIONS ACHIEVES CONVERGENCE.  ALSO, EXPLICIT
* 1-STEP SOLUTION OPTION FOR SPECIAL CASE OF PARAMETER CK=0, WHICH
* REDUCES THE ORIGINAL IMPLICIT EQUATION TO A SIMPLER EXPLICIT FORM,
* KNOWN AS THE "FLERCHINGER EQN". IMPROVED HANDLING OF SOLUTION IN THE
* LIMIT OF FREEZING POINT TEMPERATURE T0.
* ----------------------------------------------------------------------
* INPUT:

*   TKELV.........TEMPERATURE (Kelvin)
*   SMC...........TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC)
*   SH2O..........LIQUID SOIL MOISTURE CONTENT (VOLUMETRIC)
*   SMCMAX........SATURATION SOIL MOISTURE CONTENT (FROM REDPRM)
*   B.............SOIL TYPE "B" PARAMETER (FROM REDPRM)
*   PSIS..........SATURATED SOIL MATRIC POTENTIAL (FROM REDPRM)

* OUTPUT:
*   FREE..........SUPERCOOLED LIQUID WATER CONTENT
* --------------------------------------------------------------------*/
    double          DENOM, DF, DSWL, FK, SWL, SWLK;
    int             NLOG, KCOUNT;
    //      PARAMETER(CK = 0.0)
    double          CK = 8.0, ERROR = 0.005, HLICE = 3.335e5, GS = 9.81, T0 = 273.15;

#ifdef _FLUX_PIHM_
    double          MX;
    MX = *VGBETA / (1 - *VGBETA);
#else

/*----------------------------------------------------------------------
* LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)
* SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT IS
* NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES.
* --------------------------------------------------------------------*/
    BX = *BEXP;

/*----------------------------------------------------------------------
* INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
* --------------------------------------------------------------------*/
    if (*BEXP > BLIM)
        BX = BLIM;
#endif
    NLOG = 0;

/*----------------------------------------------------------------------
*  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (T0), SH2O = SMC
* --------------------------------------------------------------------*/
    KCOUNT = 0;
    //      FRH2O = SMC
    if (*TKELV > (T0 - 1.e-3))
        *FREE = *SMC;
    else
    {

/*----------------------------------------------------------------------
* OPTION 1: ITERATED SOLUTION FOR NONZERO CK
* IN KOREN ET AL, JGR, 1999, EQN 17
* ----------------------------------------------------------------------
* INITIAL GUESS FOR SWL (frozen content)
* --------------------------------------------------------------------*/
        if (CK != 0.0)
        {
            SWL = *SMC - *SH2O;

/*----------------------------------------------------------------------
* KEEP WITHIN BOUNDS.
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_
            if (SWL > (*SMC - *SMCMIN - 0.02))
                SWL = *SMC - *SMCMIN - 0.02;
#else
            if (SWL > (*SMC - 0.02))
                SWL = *SMC - 0.02;
#endif

/*----------------------------------------------------------------------
*  START OF ITERATIONS
* --------------------------------------------------------------------*/
            if (SWL < 0.)
                SWL = 0.;
          C1001:
            if (!((NLOG < 10) && (KCOUNT == 0)))
                goto C1002;
            NLOG = NLOG + 1;
#ifdef _FLUX_PIHM_
            DF =
               log ((GS / *VGALPHA / HLICE) * pow (1. + CK * SWL,
                  2.) * pow (pow ((*SMC - SWL - *SMCMIN) / (*SMCMAX -
                        *SMCMIN), MX) - 1.,
                  1. / *VGBETA)) - log (-(*TKELV - T0) / *TKELV);
            DENOM =
               2. * CK / (1. + CK * SWL) - 1. / (1 - *VGBETA) / (*SMCMAX -
               *SMCMIN) * pow ((*SMC - SWL - *SMCMIN) / (*SMCMAX - *SMCMIN),
               MX - 1.) / (pow ((*SMC - SWL - *SMCMIN) / (*SMCMAX - *SMCMIN),
                  MX) - 1.);
#else
            DF =
               log ((*PSIS * GS / HLICE) * pow (1. + CK * SWL,
                  2.) * pow (*SMCMAX / (*SMC - SWL),
                  BX)) - log (-(*TKELV - T0) / *TKELV);
            DENOM = 2. * CK / (1. + CK * SWL) + BX / (*SMC - SWL);
#endif
            SWLK = SWL - DF / DENOM;

/*----------------------------------------------------------------------
* BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_
            if (SWLK > (*SMC - *SMCMIN - 0.02))
                SWLK = *SMC - *SMCMIN - 0.02;
#else
            if (SWLK > (*SMC - 0.02))
                SWLK = *SMC - 0.02;
#endif
            if (SWLK < 0.)
                SWLK = 0.;

/*----------------------------------------------------------------------
* MATHEMATICAL SOLUTION BOUNDS APPLIED.
* --------------------------------------------------------------------*/
            DSWL = fabs (SWLK - SWL);

/*----------------------------------------------------------------------
* IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
* WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
* --------------------------------------------------------------------*/
            SWL = SWLK;
            if (DSWL <= ERROR)
                KCOUNT = KCOUNT + 1;

/*----------------------------------------------------------------------
*  END OF ITERATIONS
* ----------------------------------------------------------------------
* BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
* --------------------------------------------------------------------*/
            //          FRH2O = SMC - SWL
            goto C1001;
          C1002:
            *FREE = *SMC - SWL;
        }

/*----------------------------------------------------------------------
* END OPTION 1
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
* IN KOREN ET AL., JGR, 1999, EQN 17
* APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
* --------------------------------------------------------------------*/
        if (KCOUNT == 0)
        {
            //          printf ("Flerchinger USEd in NEW version. Iterations= %d\n",
            //             NLOG);
#ifdef _FLUX_PIHM_
            FK =
               pow (pow (-(*TKELV - T0) / *TKELV * *VGALPHA * HLICE / GS,
                  *VGBETA) + 1., 1. / MX) * (*SMCMAX - *SMCMIN);
#else
            FK =
               pow ((HLICE / (GS * (-*PSIS))) * ((*TKELV - T0) / *TKELV),
               -1 / BX) * *SMCMAX;
#endif
            //          FRH2O = MIN (FK, SMC)
            if (FK < 0.02)
                FK = 0.02;
            *FREE = FK < *SMC ? FK : *SMC;

/*----------------------------------------------------------------------
* END OPTION 2
* --------------------------------------------------------------------*/
        }
    }

/*----------------------------------------------------------------------
  END SUBROUTINE FRH2O
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
HRT (double *RHSTS, double *STC, double *SMC, double *SMCMAX, double *SMCMIN,
   int *NSOIL, double *ZSOIL, double *YY, double *ZZ1, double *TBOT,
   double *ZBOT, double *SH2O, double *DT, double *VGALPHA, double *VGBETA,
   double *F1, double *DF1, double *QUARTZ, double *CSOIL, double *AI,
   double *BI, double *CI, int *VEGTYP, int *ISURBAN)
#else
void
HRT (double *RHSTS, double *STC, double *SMC, double *SMCMAX, int *NSOIL,
   double *ZSOIL, double *YY, double *ZZ1, double *TBOT, double *ZBOT,
   double *PSISAT, double *SH2O, double *DT, double *BEXP, double *F1,
   double *DF1, double *QUARTZ, double *CSOIL, double *AI, double *BI,
   double *CI, int *VEGTYP, int *ISURBAN)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE HRT
* ----------------------------------------------------------------------
* CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
* THERMAL DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
* COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
* --------------------------------------------------------------------*/
    int             ITAVG;
    int             K;

    double          DDZ, DDZ2, DENOM, DF1K, DTSDZ, DTSDZ2, HCPCT, SSOIL, SICE,
       CSOIL_LOC;
    double         *DF1N, *QTOT, *TAVG, *TBK, *TBK1, *TSNSR, *TSURF;
    double          T0 = 273.15, CAIR = 1004.0, CICE = 2.106e6, CH2O = 4.2e6;

    DF1N = (double *)malloc (sizeof (double));
    QTOT = (double *)malloc (sizeof (double));
    TAVG = (double *)malloc (sizeof (double));
    TBK = (double *)malloc (sizeof (double));
    TBK1 = (double *)malloc (sizeof (double));
    TSNSR = (double *)malloc (sizeof (double));
    TSURF = (double *)malloc (sizeof (double));

    /*
     * urban 
     */
    if (*VEGTYP == *ISURBAN)
        CSOIL_LOC = 3.0e6;
    else
        CSOIL_LOC = *CSOIL;

/*----------------------------------------------------------------------
* INITIALIZE LOGICAL FOR SOIL LAYER TEMPERATURE AVERAGING.
* --------------------------------------------------------------------*/
    ITAVG = 1;

/*----------------------------------------------------------------------
* BEGIN SECTION FOR TOP SOIL LAYER
* ----------------------------------------------------------------------
* CALC THE HEAT CAPACITY OF THE TOP SOIL LAYER
* --------------------------------------------------------------------*/
    HCPCT =
       SH2O[0] * CH2O + (1.0 - *SMCMAX) * CSOIL_LOC + (*SMCMAX -
       SMC[0]) * CAIR + (SMC[0] - SH2O[0]) * CICE;

/*----------------------------------------------------------------------
* CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
* --------------------------------------------------------------------*/
    DDZ = 1.0 / (-0.5 * ZSOIL[1]);
    AI[0] = 0.0;
    CI[0] = (*DF1 * DDZ) / (ZSOIL[0] * HCPCT);

/*----------------------------------------------------------------------
* CALCULATE THE VERTICAL SOIL TEMP GRADIENT BTWN THE 1ST AND 2ND SOIL
* LAYERS.  THEN CALCULATE THE SUBSURFACE HEAT FLUX. USE THE TEMP
* GRADIENT AND SUBSFC HEAT FLUX TO CALC "RIGHT-HAND SIDE TENDENCY
* TERMS", OR "RHSTS", FOR TOP SOIL LAYER.
* --------------------------------------------------------------------*/
    BI[0] = -CI[0] + *DF1 / (0.5 * ZSOIL[0] * ZSOIL[0] * HCPCT * *ZZ1);
    DTSDZ = (STC[0] - STC[1]) / (-0.5 * ZSOIL[1]);
    SSOIL = *DF1 * (STC[0] - *YY) / (0.5 * ZSOIL[0] * *ZZ1);
    //      RHSTS[0] = (*DF1 * DTSDZ - SSOIL) / (ZSOIL[0] * HCPCT);
    DENOM = (ZSOIL[0] * HCPCT);

/*----------------------------------------------------------------------
* NEXT CAPTURE THE VERTICAL DIFFERENCE OF THE HEAT FLUX AT TOP AND
* BOTTOM OF FIRST SOIL LAYER FOR USE IN HEAT FLUX CONSTRAINT APPLIED TO
* POTENTIAL SOIL FREEZING/THAWING IN ROUTINE SNKSRC.
* --------------------------------------------------------------------*/
    //  *QTOT = SSOIL - *DF1*DTSDZ;
    RHSTS[0] = (*DF1 * DTSDZ - SSOIL) / DENOM;

/*----------------------------------------------------------------------
* CALCULATE FROZEN WATER CONTENT IN 1ST SOIL LAYER.
* --------------------------------------------------------------------*/
    *QTOT = -1.0 * RHSTS[0] * DENOM;

/*----------------------------------------------------------------------
* IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):
* SET TEMP "TSURF" AT TOP OF SOIL COLUMN (FOR USE IN FREEZING SOIL
* PHYSICS LATER IN FUNCTION SUBROUTINE SNKSRC).  IF SNOWPACK CONTENT IS
* ZERO, THEN TSURF EXPRESSION BELOW GIVES TSURF = SKIN TEMP.  IF
* SNOWPACK IS NONZERO (HENCE ARGUMENT ZZ1=1), THEN TSURF EXPRESSION
* BELOW YIELDS SOIL COLUMN TOP TEMPERATURE UNDER SNOWPACK.  THEN
* CALCULATE TEMPERATURE AT BOTTOM INTERFACE OF 1ST SOIL LAYER FOR USE
* LATER IN FUNCTION SUBROUTINE SNKSRC
* --------------------------------------------------------------------*/
    SICE = SMC[0] - SH2O[0];
    if (ITAVG)
    {
        *TSURF = (*YY + (*ZZ1 - 1) * STC[0]) / *ZZ1;

/*----------------------------------------------------------------------
* IF FROZEN WATER PRESENT OR ANY OF LAYER-1 MID-POINT OR BOUNDING
* INTERFACE TEMPERATURES BELOW FREEZING, THEN CALL SNKSRC TO
* COMPUTE HEAT SOURCE/SINK (AND CHANGE IN FROZEN WATER CONTENT)
* DUE TO POSSIBLE SOIL WATER PHASE CHANGE
* --------------------------------------------------------------------*/
        TBND (STC, STC + 1, ZSOIL, ZBOT, 0, NSOIL, TBK);
        if ((SICE > 0.) || (STC[0] < T0) || (*TSURF < T0) || (*TBK < T0))
        {
            //          *TSNSR = SNKSRC (TAVG,SMC(1),SH2O(1),
            TMPAVG (TAVG, TSURF, STC, TBK, ZSOIL, NSOIL, 0);
#ifdef _FLUX_PIHM_
            SNKSRC (TSNSR, TAVG, SMC, SH2O, ZSOIL, NSOIL, SMCMAX, SMCMIN,
               VGALPHA, VGBETA, DT, 0, QTOT);
#else
            SNKSRC (TSNSR, TAVG, SMC, SH2O, ZSOIL, NSOIL, SMCMAX, PSISAT,
               BEXP, DT, 0, QTOT);
#endif
            //          RHSTS(1) = RHSTS(1) - *TSNSR / ( ZSOIL(1) * HCPCT )
            RHSTS[0] = RHSTS[0] - *TSNSR / DENOM;
        }
    }
    else
    {
        //          *TSNSR = SNKSRC (STC(1),SMC(1),SH2O(1),
        if ((SICE > 0.) || (STC[0] < T0))
        {
#ifdef _FLUX_PIHM_
            SNKSRC (TSNSR, STC, SMC, SH2O, ZSOIL, NSOIL, SMCMAX, SMCMIN,
               VGALPHA, VGBETA, DT, 0, QTOT);
#else
            SNKSRC (TSNSR, STC, SMC, SH2O, ZSOIL, NSOIL, SMCMAX, PSISAT, BEXP,
               DT, 0, QTOT);
#endif
            //          RHSTS(1) = RHSTS(1) - *TSNSR / ( ZSOIL(1) * HCPCT )
            RHSTS[0] = RHSTS[0] - *TSNSR / DENOM;
        }

/*----------------------------------------------------------------------
* THIS ENDS SECTION FOR TOP SOIL LAYER.
* --------------------------------------------------------------------*/
    }

    /*
     * INITIALIZE DDZ2 
     */

/*--------------------------------------------------------------------*/

    DDZ2 = 0.0;
    DF1K = *DF1;

/*----------------------------------------------------------------------
* LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABOVE PROCESS
* (EXCEPT SUBSFC OR "GROUND" HEAT FLUX NOT REPEATED IN LOWER LAYERS)
* ----------------------------------------------------------------------
* CALCULATE HEAT CAPACITY FOR THIS SOIL LAYER.
* --------------------------------------------------------------------*/
    for (K = 1; K < *NSOIL; K++)
    {
        HCPCT =
           SH2O[K] * CH2O + (1.0 - *SMCMAX) * CSOIL_LOC + (*SMCMAX -
           SMC[K]) * CAIR + (SMC[K] - SH2O[K]) * CICE;

/*----------------------------------------------------------------------
* THIS SECTION FOR LAYER 2 OR GREATER, BUT NOT LAST LAYER.
* ----------------------------------------------------------------------
* CALCULATE THERMAL DIFFUSIVITY FOR THIS LAYER.
* --------------------------------------------------------------------*/
        if (K != *NSOIL - 1)
        {

/*----------------------------------------------------------------------
* CALC THE VERTICAL SOIL TEMP GRADIENT THRU THIS LAYER
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_
            TDFCND (DF1N, SMC + K, QUARTZ, SMCMAX, SMCMIN, SH2O + K);
#else
            TDFCND (DF1N, SMC + K, QUARTZ, SMCMAX, SH2O + K);
#endif

            /*
             * urban 
             */
            if (*VEGTYP == *ISURBAN)
                *DF1N = 3.24;

            DENOM = 0.5 * (ZSOIL[K - 1] - ZSOIL[K + 1]);

/*----------------------------------------------------------------------
* CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
* --------------------------------------------------------------------*/
            DTSDZ2 = (STC[K] - STC[K + 1]) / DENOM;
            DDZ2 = 2. / (ZSOIL[K - 1] - ZSOIL[K + 1]);

/*----------------------------------------------------------------------
* IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE
* TEMP AT BOTTOM OF LAYER.
* --------------------------------------------------------------------*/
            CI[K] = -*DF1N * DDZ2 / ((ZSOIL[K - 1] - ZSOIL[K]) * HCPCT);
            if (ITAVG)
                TBND (STC + K, STC + K + 1, ZSOIL, ZBOT, K, NSOIL, TBK1);
        }
        else
        {

/*----------------------------------------------------------------------
* SPECIAL CASE OF BOTTOM SOIL LAYER:  CALCULATE THERMAL DIFFUSIVITY FOR
* BOTTOM LAYER.
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* CALC THE VERTICAL SOIL TEMP GRADIENT THRU BOTTOM LAYER.
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_
            TDFCND (DF1N, SMC + K, QUARTZ, SMCMAX, SMCMIN, SH2O + K);
#else
            TDFCND (DF1N, SMC + K, QUARTZ, SMCMAX, SH2O + K);
#endif

            /*
             * urban 
             */
            if (*VEGTYP == *ISURBAN)
                *DF1N = 3.24;

            DENOM = 0.5 * (ZSOIL[K - 1] + ZSOIL[K]) - *ZBOT;

/*----------------------------------------------------------------------
* SET MATRIX COEF, CI TO ZERO IF BOTTOM LAYER.
* --------------------------------------------------------------------*/
            DTSDZ2 = (STC[K] - *TBOT) / DENOM;

/*----------------------------------------------------------------------
* IF TEMPERATURE AVERAGING INVOKED (ITAVG=TRUE; ELSE SKIP):  CALCULATE
* TEMP AT BOTTOM OF LAST LAYER.
* --------------------------------------------------------------------*/
            CI[K] = 0.;
            if (ITAVG)
                TBND (STC + K, TBOT, ZSOIL, ZBOT, K, NSOIL, TBK1);

/*----------------------------------------------------------------------
* THIS ENDS SPECIAL LOOP FOR BOTTOM LAYER.
* --------------------------------------------------------------------*/
        }

/*----------------------------------------------------------------------
* CALCULATE RHSTS FOR THIS LAYER AFTER CALC'NG A PARTIAL PRODUCT.
* --------------------------------------------------------------------*/
        DENOM = (ZSOIL[K] - ZSOIL[K - 1]) * HCPCT;
        RHSTS[K] = (*DF1N * DTSDZ2 - DF1K * DTSDZ) / DENOM;
        *QTOT = -1.0 * DENOM * RHSTS[K];

        SICE = SMC[K] - SH2O[K];

        if (ITAVG)
        {
            TMPAVG (TAVG, TBK, STC + K, TBK1, ZSOIL, NSOIL, K);
            //                  *TSNSR = SNKSRC(TAVG,SMC(K),SH2O(K),ZSOIL,NSOIL,
            if ((SICE > 0.) || (STC[K] < T0) || (*TBK < T0) || (*TBK1 < T0))
            {
#ifdef _FLUX_PIHM_
                SNKSRC (TSNSR, TAVG, SMC + K, SH2O + K, ZSOIL, NSOIL, SMCMAX,
                   SMCMIN, VGALPHA, VGBETA, DT, K, QTOT);
#else
                SNKSRC (TSNSR, TAVG, SMC + K, SH2O + K, ZSOIL, NSOIL, SMCMAX,
                   PSISAT, BEXP, DT, K, QTOT);
#endif
                RHSTS[K] = RHSTS[K] - *TSNSR / DENOM;
            }
        }
        else
        {
            //            TSNSR = SNKSRC(STC(K),SMC(K),SH2O(K),ZSOIL,NSOIL,
            if ((SICE > 0.) || (STC[K] < T0))
            {
#ifdef _FLUX_PIHM_
                SNKSRC (TSNSR, STC + K, SMC + K, SH2O + K, ZSOIL, NSOIL,
                   SMCMAX, SMCMIN, VGALPHA, VGBETA, DT, K, QTOT);
#else
                SNKSRC (TSNSR, STC + K, SMC + K, SH2O + K, ZSOIL, NSOIL,
                   SMCMAX, PSISAT, BEXP, DT, K, QTOT);
#endif
                RHSTS[K] = RHSTS[K] - *TSNSR / DENOM;
            }
        }

/*----------------------------------------------------------------------
* CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER.
* --------------------------------------------------------------------*/
        AI[K] = -DF1K * DDZ / ((ZSOIL[K - 1] - ZSOIL[K]) * HCPCT);

/*----------------------------------------------------------------------
* RESET VALUES OF DF1, DTSDZ, DDZ, AND TBK FOR LOOP TO NEXT SOIL LAYER.
* --------------------------------------------------------------------*/
        BI[K] = -(AI[K] + CI[K]);
        *TBK = *TBK1;
        DF1K = *DF1N;
        DTSDZ = DTSDZ2;
        DDZ = DDZ2;
    }

    free (DF1N);
    free (QTOT);
    free (TAVG);
    free (TBK);
    free (TBK1);
    free (TSNSR);
    free (TSURF);

/*----------------------------------------------------------------------
  END SUBROUTINE HRT
* --------------------------------------------------------------------*/
}


void
HSTEP (double *STCOUT, double *STCIN, double *RHSTS, double *DT, int *NSOIL,
   double *AI, double *BI, double *CI)
{

/*----------------------------------------------------------------------
* SUBROUTINE HSTEP
* ----------------------------------------------------------------------
* CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
* --------------------------------------------------------------------*/
    int             K;

    double          RHSTSin[*NSOIL];
    double          CIin[*NSOIL];

/*----------------------------------------------------------------------
* CREATE FINITE DIFFERENCE VALUES FOR USE IN ROSR12 ROUTINE
* --------------------------------------------------------------------*/
    for (K = 0; K < *NSOIL; K++)
    {
        RHSTS[K] = RHSTS[K] * *DT;
        AI[K] = AI[K] * *DT;
        BI[K] = 1. + BI[K] * *DT;
        CI[K] = CI[K] * *DT;
    }

/*----------------------------------------------------------------------
* COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
* --------------------------------------------------------------------*/
    for (K = 0; K < *NSOIL; K++)
        RHSTSin[K] = RHSTS[K];
    for (K = 0; K < *NSOIL; K++)
        CIin[K] = CI[K];

/*----------------------------------------------------------------------
* SOLVE THE TRI-DIAGONAL MATRIX EQUATION
* --------------------------------------------------------------------*/
    ROSR12 (CI, AI, BI, CIin, RHSTSin, RHSTS, NSOIL);

/*----------------------------------------------------------------------
* CALC/UPDATE THE SOIL TEMPS USING MATRIX SOLUTION
* --------------------------------------------------------------------*/
    for (K = 0; K < *NSOIL; K++)
        STCOUT[K] = STCIN[K] + CI[K];

/*----------------------------------------------------------------------
* END SUBROUTINE HSTEP
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
NOPAC (double *ETP, double *ETA, double *PRCP, double *PCPDRP, double *SMC,
   double *SMCMAX, double *SMCMIN, double *SMCWLT, double *SMCREF,
   double *SMCDRY, double *CMC, double *CMCMAX, int *NSOIL, double *DT,
   double *SHDFAC, double *SBETA, double *Q2, double *T1, double *SFCTMP,
   double *T24, double *TH2, double *FDOWN, double *F1, double *EMISSI,
   double *SSOIL, double *STC, double *EPSCA, double *VGALPHA, double *VGBETA,
   double *MACKSAT, double *AREAF, int *NMACD, int *MAC_STATUS, int *NWTBL, double *PC,
   double *RCH, double *RR, double *CFACTR, double *SH2O, double *FRZFACT,
   double *ZSOIL, double *DKSAT, double *TBOT, double *ZBOT, double *INF,
   double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC, double *ET,
   double *ETT, int *NROOT, double *RTDIS, double *QUARTZ, double *FXEXP,
   double *CSOIL, double *BETA, double *DRIP, double *DEW, double *FLX1,
   double *FLX3, int *VEGTYP, int *ISURBAN)
#else
void
NOPAC (double *ETP, double *ETA, double *PRCP, double *SMC, double *SMCMAX,
   double *SMCWLT, double *SMCREF, double *SMCDRY, double *CMC,
   double *CMCMAX, int *NSOIL, double *DT, double *SHDFAC, double *SBETA,
   double *Q2, double *T1, double *SFCTMP, double *T24, double *TH2,
   double *FDOWN, double *F1, double *EMISSI, double *SSOIL, double *STC,
   double *EPSCA, double *BEXP, double *PC, double *RCH, double *RR,
   double *CFACTR, double *SH2O, double *SLOPE, double *KDT, double *FRZFACT,
   double *PSISAT, double *ZSOIL, double *DKSAT, double *DWSAT, double *TBOT,
   double *ZBOT, double *RUNOFF1, double *RUNOFF2, double *RUNOFF3,
   double *EDIR, double *EC, double *ET, double *ETT, int *NROOT,
   double *RTDIS, double *QUARTZ, double *FXEXP, double *CSOIL, double *BETA,
   double *DRIP, double *DEW, double *FLX1, double *FLX3, int *VEGTYP,
   int *ISURBAN)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE NOPAC
* ----------------------------------------------------------------------
* CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES AND UPDATE SOIL MOISTURE
* CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN NO SNOW PACK IS
* PRESENT.
* --------------------------------------------------------------------*/

    int             K;

    double          ET1[*NSOIL];
    double         *EC1, *EDIR1, *ETT1, *DF1, *ETA1, *ETP1, *PRCP1, *YY, *ZZ1;
    double          YYNUM;

    EC1 = (double *)malloc (sizeof (double));
    EDIR1 = (double *)malloc (sizeof (double));
    ETT1 = (double *)malloc (sizeof (double));
    DF1 = (double *)malloc (sizeof (double));
    ETA1 = (double *)malloc (sizeof (double));
    ETP1 = (double *)malloc (sizeof (double));
    PRCP1 = (double *)malloc (sizeof (double));
    YY = (double *)malloc (sizeof (double));
    ZZ1 = (double *)malloc (sizeof (double));

/*----------------------------------------------------------------------
* EXECUTABLE CODE BEGINS HERE:
* CONVERT ETP Fnd PRCP FROM KG M-2 S-1 TO M S-1 AND INITIALIZE DEW.
* --------------------------------------------------------------------*/
    *PRCP1 = *PRCP * 0.001;
    *ETP1 = *ETP * 0.001;
    *DEW = 0.0;

/*----------------------------------------------------------------------
* INITIALIZE EVAP TERMS.
* --------------------------------------------------------------------*/
    *EDIR = 0.;
    *EDIR1 = 0.;
    *EC1 = 0.;
    *EC = 0.;
    for (K = 0; K < *NSOIL; K++)
    {
        ET[K] = 0.;
        ET1[K] = 0.;
    }
    *ETT = 0.;
    *ETT1 = 0.;

    if (*ETP > 0.0)
    {
#ifdef _FLUX_PIHM_
        EVAPO (ETA1, SMC, NSOIL, CMC, ETP1, DT, ZSOIL, SH2O, SMCMAX, PC,
           SMCWLT, DKSAT, SMCREF, SHDFAC, CMCMAX, SMCDRY, CFACTR, EDIR1, EC1,
           ET1, ETT1, SFCTMP, Q2, NROOT, RTDIS, FXEXP);
        SMFLX (SMC, NSOIL, CMC, DT, PRCP1, PCPDRP, ZSOIL, SH2O, FRZFACT,
           SMCMAX, SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, MAC_STATUS, NWTBL,
           SMCWLT, DKSAT, SHDFAC, CMCMAX, INF, RUNOFF2, RUNOFF3, EDIR1, EC1,
           ET1, DRIP);
#else
        EVAPO (ETA1, SMC, NSOIL, CMC, ETP1, DT, ZSOIL, SH2O, SMCMAX, BEXP, PC,
           SMCWLT, DKSAT, DWSAT, SMCREF, SHDFAC, CMCMAX, SMCDRY, CFACTR,
           EDIR1, EC1, ET1, ETT1, SFCTMP, Q2, NROOT, RTDIS, FXEXP);
        SMFLX (SMC, NSOIL, CMC, DT, PRCP1, ZSOIL, SH2O, SLOPE, KDT, FRZFACT,
           SMCMAX, BEXP, SMCWLT, DKSAT, DWSAT, SHDFAC, CMCMAX, RUNOFF1,
           RUNOFF2, RUNOFF3, EDIR1, EC1, ET1, DRIP);
#endif

/*----------------------------------------------------------------------
* CONVERT MODELED EVAPOTRANSPIRATION FROM  M S-1  TO  KG M-2 S-1.
* --------------------------------------------------------------------*/

        *ETA = *ETA1 * 1000.0;
    }

/*----------------------------------------------------------------------
* IF ETP < 0, ASSUME DEW FORMS (TRANSFORM ETP1 INTO DEW AND REINITIALIZE
* ETP1 TO ZERO).
* --------------------------------------------------------------------*/
    else
    {
        *DEW = -*ETP1;

/*----------------------------------------------------------------------
* CONVERT PRCP FROM 'KG M-2 S-1' TO 'M S-1' AND ADD DEW AMOUNT.
* --------------------------------------------------------------------*/

        *PRCP1 = *PRCP1 + *DEW;
#ifdef _FLUX_PIHM_
        SMFLX (SMC, NSOIL, CMC, DT, PRCP1, PCPDRP, ZSOIL, SH2O, FRZFACT,
           SMCMAX, SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, MAC_STATUS, NWTBL,
           SMCWLT, DKSAT, SHDFAC, CMCMAX, INF, RUNOFF2, RUNOFF3, EDIR1, EC1,
           ET1, DRIP);
#else
        SMFLX (SMC, NSOIL, CMC, DT, PRCP1, ZSOIL, SH2O, SLOPE, KDT, FRZFACT,
           SMCMAX, BEXP, SMCWLT, DKSAT, DWSAT, SHDFAC, CMCMAX, RUNOFF1,
           RUNOFF2, RUNOFF3, EDIR1, EC1, ET1, DRIP);
#endif

/*----------------------------------------------------------------------
* CONVERT MODELED EVAPOTRANSPIRATION FROM 'M S-1' TO 'KG M-2 S-1'.
* --------------------------------------------------------------------*/
        //      *ETA = *ETA1 * 1000.0
    }

/*----------------------------------------------------------------------
* BASED ON ETP AND E VALUES, DETERMINE BETA
* --------------------------------------------------------------------*/

    if (*ETP <= 0.0)
    {
        *BETA = 0.0;
        *ETA = *ETP;
        if (*ETP < 0.0)
            *BETA = 1.0;
    }
    else
        *BETA = *ETA / *ETP;

/*----------------------------------------------------------------------
* CONVERT MODELED EVAPOTRANSPIRATION COMPONENTS 'M S-1' TO 'KG M-2 S-1'.
* --------------------------------------------------------------------*/
    *EDIR = *EDIR1 * 1000.;
    *EC = *EC1 * 1000.;
    for (K = 0; K < *NSOIL; K++)
        ET[K] = ET1[K] * 1000.;
    *ETT = *ETT1 * 1000.;

/*----------------------------------------------------------------------
* GET SOIL THERMAL DIFFUXIVITY/CONDUCTIVITY FOR TOP SOIL LYR,
* CALC. ADJUSTED TOP LYR SOIL TEMP AND ADJUSTED SOIL FLUX, THEN
* CALL SHFLX TO COMPUTE/UPDATE SOIL HEAT FLUX AND SOIL TEMPS.
* --------------------------------------------------------------------*/

#ifdef _FLUX_PIHM_
    TDFCND (DF1, SMC, QUARTZ, SMCMAX, SMCMIN, SH2O);
#else
    TDFCND (DF1, SMC, QUARTZ, SMCMAX, SH2O);
#endif

    /*
     * urban 
     */
    if (*VEGTYP == *ISURBAN)
        *DF1 = 3.24;

/*----------------------------------------------------------------------
* VEGETATION GREENNESS FRACTION REDUCTION IN SUBSURFACE HEAT FLUX
* VIA REDUCTION FACTOR, WHICH IS CONVENIENT TO APPLY HERE TO THERMAL
* DIFFUSIVITY THAT IS LATER USED IN HRT TO COMPUTE SUB SFC HEAT FLUX
* (SEE ADDITIONAL COMMENTS ON VEG EFFECT SUB-SFC HEAT FLX IN
* ROUTINE SFLX)
* --------------------------------------------------------------------*/
    *DF1 = *DF1 * exp (*SBETA * *SHDFAC);

/*----------------------------------------------------------------------
* COMPUTE INTERMEDIATE TERMS PASSED TO ROUTINE HRT (VIA ROUTINE
* SHFLX BELOW) FOR USE IN COMPUTING SUBSURFACE HEAT FLUX IN HRT
* --------------------------------------------------------------------*/
    YYNUM = *FDOWN - *EMISSI * SIGMA * *T24;
    *YY = *SFCTMP + (YYNUM / *RCH + *TH2 - *SFCTMP - *BETA * *EPSCA) / *RR;

    *ZZ1 = *DF1 / (-0.5 * ZSOIL[0] * *RCH * *RR) + 1.0;

    /*
     * urban 
     */
#ifdef _FLUX_PIHM_
    SHFLX (SSOIL, STC, SMC, SMCMAX, SMCMIN, NSOIL, T1, DT, YY, ZZ1, ZSOIL,
       TBOT, ZBOT, SMCWLT, SH2O, VGALPHA, VGBETA, F1, DF1, QUARTZ, CSOIL,
       VEGTYP, ISURBAN);
#else
    SHFLX (SSOIL, STC, SMC, SMCMAX, NSOIL, T1, DT, YY, ZZ1, ZSOIL, TBOT, ZBOT,
       SMCWLT, PSISAT, SH2O, BEXP, F1, DF1, QUARTZ, CSOIL, VEGTYP, ISURBAN);
#endif

/*----------------------------------------------------------------------
* SET FLX1 AND FLX3 (SNOPACK PHASE CHANGE HEAT FLUXES) TO ZERO SINCE
* THEY ARE NOT USED HERE IN SNOPAC.  FLX2 (FREEZING RAIN HEAT FLUX) WAS
* SIMILARLY INITIALIZED IN THE PENMAN ROUTINE.
* --------------------------------------------------------------------*/
    *FLX1 = CPH2O * *PRCP * (*T1 - *SFCTMP);
    *FLX3 = 0.0;

    free (EC1);
    free (EDIR1);
    free (ETT1);
    free (DF1);
    free (ETA1);
    free (ETP1);
    free (PRCP1);
    free (YY);
    free (ZZ1);

/*----------------------------------------------------------------------
  END SUBROUTINE NOPAC
* --------------------------------------------------------------------*/
}

void PENMAN (double *SFCTMP, double *SFCPRS, double *CH, double *T2V, double *TH2, double *PRCP, double *FDOWN, double *T24, double *SSOIL, double *Q2, double *Q2SAT, double *ETP, double *RCH, double *EPSCA, double *RR, int *SNOWNG, int *FRZGRA, double *DQSDT2, double *FLX2, double *EMISSI_IN, double *SNEQV, double *T1, double *SNCOVR)
{

/*----------------------------------------------------------------------
* SUBROUTINE PENMAN
* ----------------------------------------------------------------------
* CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT.  VARIOUS
* PARTIAL SUMS/PRODUCTS ARE ALSO CALCULATED AND PASSED BACK TO THE
* CALLING ROUTINE FOR LATER USE.
* --------------------------------------------------------------------*/

    double          A, DELTA, FNET, RAD, RHO, EMISSI, ELCP1, LVS;

    double          ELCP = 2.4888e+3, LSUBC = 2.501000e+6;

/*----------------------------------------------------------------------
* EXECUTABLE CODE BEGINS HERE:
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* PREPARE PARTIAL QUANTITIES FOR PENMAN EQUATION.
* --------------------------------------------------------------------*/
    EMISSI = *EMISSI_IN;
    ELCP1 = (1.0 - *SNCOVR) * ELCP + *SNCOVR * ELCP * LSUBS / LSUBC;
    LVS = (1.0 - *SNCOVR) * LSUBC + *SNCOVR * LSUBS;

    *FLX2 = 0.0;
    //DELTA = ELCP * DQSDT2
    DELTA = ELCP1 * *DQSDT2;
    *T24 = *SFCTMP * *SFCTMP * *SFCTMP * *SFCTMP;
    //RR = T24 * 6.48E-8 / (SFCPRS * CH) + 1.0
    *RR = EMISSI * *T24 * 6.48e-8 / (*SFCPRS * *CH) + 1.0;
    RHO = *SFCPRS / (RD * *T2V);

/*----------------------------------------------------------------------
* ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT
* EFFECTS CAUSED BY FALLING PRECIPITATION.
* --------------------------------------------------------------------*/
    *RCH = RHO * CP * *CH;
    if (!*SNOWNG)
    {
        if (*PRCP > 0.0)
            *RR = *RR + CPH2O * *PRCP / *RCH;
    }
    else
        *RR = *RR + CPICE * *PRCP / *RCH;

/*----------------------------------------------------------------------
* INCLUDE THE LATENT HEAT EFFECTS OF FRZNG RAIN CONVERTING TO ICE ON
* IMPACT IN THE CALCULATION OF FLX2 AND FNET.
* --------------------------------------------------------------------*/
    //      FNET = FDOWN - SIGMA * T24- SSOIL
    FNET = *FDOWN - EMISSI * SIGMA * *T24 - *SSOIL;
    if (*FRZGRA)
    {
        *FLX2 = -LSUBF * (*PRCP);
        FNET = FNET - *FLX2;

/*----------------------------------------------------------------------
* FINISH PENMAN EQUATION CALCULATIONS.
* --------------------------------------------------------------------*/
    }
    RAD = FNET / *RCH + *TH2 - *SFCTMP;
    //  A = ELCP * (Q2SAT - Q2)
    A = ELCP1 * (*Q2SAT - *Q2);
    *EPSCA = (A * *RR + RAD * DELTA) / (DELTA + *RR);
    //  ETP = EPSCA * RCH / LSUBC;
    *ETP = *EPSCA * *RCH / LVS;
#ifdef _DEBUG_
    printf("RR = %lf, RAD = %lf, DELTA = %lf, CH = %lf, EPSCA = %f, ETP = %lg, DQSDT2 = %lf\n", *RR, RAD, DELTA, *CH, *EPSCA, *ETP, *DQSDT2);
#endif

/*----------------------------------------------------------------------
  END SUBROUTINE PENMAN
* --------------------------------------------------------------------*/
}

void REDPRM (GRID_TYPE * NOAH, LSM_STRUCT LSM, double *ZSOIL)
{

/*----------------------------------------------------------------------
* Internally set (default valuess)
* all soil and vegetation parameters required for the execusion oF
* the Noah lsm are defined in VEGPARM.TBL, SOILPARM.TB, and GENPARM.TBL.
* ----------------------------------------------------------------------
*     Vegetation parameters:
*             ALBBRD: SFC background snow-free albedo
*             CMXTBL: MAX CNPY Capacity
*              Z0BRD: Background roughness length
*             SHDFAC: Green vegetation fraction
*              NROOT: Rooting depth
*              RSMIN: Mimimum stomatal resistance
*              RSMAX: Max. stomatal resistance
*                RGL: Parameters used in radiation stress function
*                 HS: Parameter used in vapor pressure deficit functio
*               TOPT: Optimum transpiration air temperature.
*             CMCMAX: Maximum canopy water capacity
*             CFACTR: Parameter used in the canopy inteception calculation
*               SNUP: Threshold snow depth (in water equivalent m) that
*                     implies 100 percent snow cover
*                LAI: Leaf area index
*
* ----------------------------------------------------------------------
*      Soil parameters:
*        SMCMAX: MAX soil moisture content (porosity)
*        SMCREF: Reference soil moisture  (field capacity)
*        SMCWLT: Wilting point soil moisture
*        SMCWLT: Air dry soil moist content limits
*       SSATPSI: SAT (saturation) soil potential
*         DKSAT: SAT soil conductivity
*          BEXP: B parameter
*        SSATDW: SAT soil diffusivity
*           F1: Soil thermal diffusivity/conductivity coef.
*        QUARTZ: Soil quartz content
*  Modified by F. Chen (12/22/97)  to use the STATSGO soil map
*  Modified By F. Chen (01/22/00)  to include PLaya, Lava, and White San
*  Modified By F. Chen (08/05/02)  to include additional parameters for the Noah
* NOTE: SATDW = BB*SATDK*(SATPSI/MAXSMC)
*         F11 = ALOG10(SATPSI) + BB*ALOG10(MAXSMC) + 2.0
*       REFSMC1=MAXSMC*(5.79E-9/SATDK)**(1/(2*BB+3)) 5.79E-9 m/s= 0.5 mm
*       REFSMC=REFSMC1+1./3.(MAXSMC-REFSMC1)
*       WLTSMC1=MAXSMC*(200./SATPSI)**(-1./BB)    (Wetzel and Chang, 198
*       WLTSMC=WLTSMC1-0.5*WLTSMC1
* Note: the values for playa is set for it to have a thermal conductivit
* as sand and to have a hydrulic conductivity as clay
*
* ----------------------------------------------------------------------
* Class parameter 'SLOPETYP' was included to estimate linear reservoir
* coefficient 'SLOPE' to the baseflow runoff out of the bottom layer.
* lowest class (slopetyp=0) means highest slope parameter = 1.
* definition of slopetyp from 'zobler' slope type:
* slope class  percent slope
* 1            0-8
* 2            8-30
* 3            > 30
* 4            0-30
* 5            0-8 & > 30
* 6            8-30 & > 30
* 7            0-8, 8-30, > 30
* 9            GLACIAL ICE
* BLANK        OCEAN/SEA
*       SLOPE_DATA: linear reservoir coefficient
*       SBETA_DATA: parameter used to caluculate vegetation effect on soil heat
*       FXEXP_DAT:  soil evaporation exponent used in DEVAP
*       CSOIL_DATA: soil heat capacity [J M-3 K-1]
*       SALP_DATA: shape parameter of  distribution function of snow cover
*       REFDK_DATA and REFKDT_DATA: parameters in the surface runoff parameteriz
*       FRZK_DATA: frozen ground parameter
*       ZBOT_DATA: depth[M] of lower boundary soil temperature
*       CZIL_DATA: calculate roughness length of heat
*       SMLOW_DATA and MHIGH_DATA: two soil moisture wilt, soil moisture referen
*                 parameters
* Set maximum number of soil-, veg-, and slopetyp in data statement.
* --------------------------------------------------------------------*/

    int             I;

    double          FRZFACT;

    /*
     * SAVE
     * * ----------------------------------------------------------------------
     */
    if (NOAH->SOILTYP > LSM->SOILTBL.SLCATS)
    {
        printf ("Warning: too many input soil types\n");
        exit (0);
    }
    if (NOAH->VEGTYP > LSM->VEGTBL.LUCATS)
    {
        printf ("Warning: too many input landuse types\n");
        exit (0);
    }
    if (NOAH->SLOPETYP > LSM->GENPRMT.SLPCATS)
    {
        printf ("Warning: too many input slope types\n");
        exit (0);
    }

/*----------------------------------------------------------------------
*  SET-UP SOIL PARAMETERS
* --------------------------------------------------------------------*/
    NOAH->CSOIL = LSM->GENPRMT.CSOIL_DATA;
#ifdef _FLUX_PIHM_
    NOAH->VGALPHA = LSM->SOILTBL.VGA[NOAH->SOILTYP - 1];
    NOAH->VGBETA = LSM->SOILTBL.VGB[NOAH->SOILTYP - 1];
    NOAH->SMCMIN = LSM->SOILTBL.MINSMC[NOAH->SOILTYP - 1];
    NOAH->MACKSAT = LSM->SOILTBL.MACKSAT[NOAH->SOILTYP - 1];
    NOAH->AREAF = LSM->SOILTBL.AREAF[NOAH->SOILTYP - 1];
    NOAH->NMACD = LSM->SOILTBL.NMACD[NOAH->SOILTYP - 1];
#else
    NOAH->BEXP = LSM->SOILTBL.BB[NOAH->SOILTYP - 1];
    NOAH->PSISAT = LSM->SOILTBL.SATPSI[NOAH->SOILTYP - 1];
    NOAH->DWSAT = LSM->SOILTBL.SATDW[NOAH->SOILTYP - 1];
#endif
    NOAH->DKSAT = LSM->SOILTBL.SATDK[NOAH->SOILTYP - 1];
    NOAH->F1 = LSM->SOILTBL.F11[NOAH->SOILTYP - 1];
    NOAH->QUARTZ = LSM->SOILTBL.QTZ[NOAH->SOILTYP - 1];
    NOAH->SMCDRY = LSM->SOILTBL.DRYSMC[NOAH->SOILTYP - 1];
    NOAH->SMCMAX = LSM->SOILTBL.MAXSMC[NOAH->SOILTYP - 1];
    NOAH->SMCREF = LSM->SOILTBL.REFSMC[NOAH->SOILTYP - 1];
    NOAH->SMCWLT = LSM->SOILTBL.WLTSMC[NOAH->SOILTYP - 1];

/*----------------------------------------------------------------------
* Set-up universal parameters (not dependent on SOILTYP, VEGTYP or
* SLOPETYP)
* --------------------------------------------------------------------*/
    NOAH->ZBOT = LSM->GENPRMT.ZBOT_DATA;
    NOAH->SALP = LSM->GENPRMT.SALP_DATA;
    NOAH->SBETA = LSM->GENPRMT.SBETA_DATA;
    NOAH->FRZK = LSM->GENPRMT.FRZK_DATA;
    NOAH->FXEXP = LSM->GENPRMT.FXEXP_DATA;
    NOAH->PTU = 0.;             /* (not used yet) to satisify intent(out) */
    NOAH->CZIL = LSM->GENPRMT.CZIL_DATA;
    NOAH->LVCOEF = LSM->GENPRMT.LVCOEF_DATA;
#ifndef _FLUX_PIHM_
    NOAH->SLOPE = LSM->GENPRMT.SLOPE_DATA[NOAH->SLOPETYP - 1];
    NOAH->REFKDT = LSM->GENPRMT.REFKDT_DATA;
    NOAH->REFDK = LSM->GENPRMT.REFDK_DATA;
    NOAH->KDT = NOAH->REFKDT * NOAH->DKSAT / NOAH->REFDK;
#endif

/*----------------------------------------------------------------------
* TO ADJUST FRZK PARAMETER TO ACTUAL SOIL TYPE: FRZK * FRZFACT
* --------------------------------------------------------------------*/
    FRZFACT = (NOAH->SMCMAX / NOAH->SMCREF) * (0.412 / 0.468);
    NOAH->FRZX = NOAH->FRZK * FRZFACT;

/*----------------------------------------------------------------------
* SET-UP VEGETATION PARAMETERS
* --------------------------------------------------------------------*/
    NOAH->TOPT = LSM->VEGTBL.TOPT_DATA;
    NOAH->CFACTR = LSM->VEGTBL.CFACTR_DATA;
    NOAH->RSMAX = LSM->VEGTBL.RSMAX_DATA;
#ifdef _FLUX_PIHM_
    NOAH->CMCMAX = LSM->VEGTBL.CMCFACTRTBL[NOAH->VEGTYP - 1] * NOAH->XLAI;
#else
    NOAH->CMCMAX = LSM->VEGTBL.CMCMAX_DATA;
#endif
    NOAH->NROOT = LSM->VEGTBL.NROTBL[NOAH->VEGTYP - 1];
    NOAH->SNUP = LSM->VEGTBL.SNUPTBL[NOAH->VEGTYP - 1];
    NOAH->RSMIN = LSM->VEGTBL.RSTBL[NOAH->VEGTYP - 1];
    NOAH->RGL = LSM->VEGTBL.RGLTBL[NOAH->VEGTYP - 1];
    NOAH->HS = LSM->VEGTBL.HSTBL[NOAH->VEGTYP - 1];
    NOAH->EMISSMIN = LSM->VEGTBL.EMISSMINTBL[NOAH->VEGTYP - 1];
    NOAH->EMISSMAX = LSM->VEGTBL.EMISSMAXTBL[NOAH->VEGTYP - 1];
    NOAH->LAIMIN = LSM->VEGTBL.LAIMINTBL[NOAH->VEGTYP - 1];
    NOAH->LAIMAX = LSM->VEGTBL.LAIMAXTBL[NOAH->VEGTYP - 1];
    NOAH->Z0MIN = LSM->VEGTBL.Z0MINTBL[NOAH->VEGTYP - 1];
    NOAH->Z0MAX = LSM->VEGTBL.Z0MAXTBL[NOAH->VEGTYP - 1];
    NOAH->ALBEDOMIN = LSM->VEGTBL.ALBEDOMINTBL[NOAH->VEGTYP - 1];
    NOAH->ALBEDOMAX = LSM->VEGTBL.ALBEDOMAXTBL[NOAH->VEGTYP - 1];

    NOAH->ISURBAN = LSM->VEGTBL.ISURBAN;

    if (NOAH->VEGTYP == LSM->VEGTBL.BARE)
        NOAH->SHDFAC = 0.0;

    if (NOAH->NROOT > NOAH->NSOIL)
    {
        printf ("Error: too many root layers %d, %d\n", NOAH->NSOIL,
           NOAH->NROOT);
        exit (0);
    }

/*----------------------------------------------------------------------
* CALCULATE ROOT DISTRIBUTION.  PRESENT VERSION ASSUMES UNIFORM
* DISTRIBUTION BASED ON SOIL LAYER DEPTHS.
* --------------------------------------------------------------------*/
    for (I = 0; I < NOAH->NROOT; I++)
        NOAH->RTDIS[I] = -NOAH->SLDPTH[I] / ZSOIL[NOAH->NROOT - 1];

/*----------------------------------------------------------------------
*  SET-UP SLOPE PARAMETER
* --------------------------------------------------------------------*/

    /*
     * print*,'end of PRMRED'
     * *       print*,'VEGTYP',VEGTYP,'SOILTYP',SOILTYP,'SLOPETYP',SLOPETYP,    &
     * *    & 'CFACTR',CFACTR,'CMCMAX',CMCMAX,'RSMAX',RSMAX,'TOPT',TOPT,        &
     * *    & 'REFKDT',REFKDT,'KDT',KDT,'SBETA',SBETA, 'SHDFAC',SHDFAC,         &
     * *    &  'RSMIN',RSMIN,'RGL',RGL,'HS',HS,'ZBOT',ZBOT,'FRZX',FRZX,         &
     * *    &  'PSISAT',PSISAT,'SLOPE',SLOPE,'SNUP',SNUP,'SALP',SALP,'BEXP',    &
     * *    &   BEXP,                                                           &
     * *    &  'DKSAT',DKSAT,'DWSAT',DWSAT,                                     &
     * *    &  'SMCMAX',SMCMAX,'SMCWLT',SMCWLT,'SMCREF',SMCREF,'SMCDRY',SMCDRY, &
     * *    &  'F1',F1,'QUARTZ',QUARTZ,'FXEXP',FXEXP,                           &
     * *    &  'RTDIS',RTDIS,'SLDPTH',SLDPTH,'ZSOIL',ZSOIL, 'NROOT',NROOT,      &
     * *    &  'NSOIL',NSOIL,'Z0',Z0,'CZIL',CZIL,'LAI',LAI,                     &
     * *    &  'CSOIL',CSOIL,'PTU',PTU,                                         &
     * *    &  'LOCAL', LOCAL
     */

}

void
ROSR12 (double *P, double *A, double *B, double *C, double *D, double *DELTA,
   int *NSOIL)
{

/*----------------------------------------------------------------------
* SUBROUTINE ROSR12
* ----------------------------------------------------------------------
* INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
* ###                                            ### ###  ###   ###  ###
* #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
* #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
* # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
* # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
* # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
* # .                                          .   # #  .   # = #   .  #
* # .                                          .   # #  .   #   #   .  #
* # .                                          .   # #  .   #   #   .  #
* # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
* # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
* # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
* ###                                            ### ###  ###   ###  ###
* --------------------------------------------------------------------*/

    int             K, KK;

/*----------------------------------------------------------------------
* INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
* --------------------------------------------------------------------*/
    C[*NSOIL - 1] = 0.0;
    P[0] = -C[0] / B[0];

/*----------------------------------------------------------------------
* SOLVE THE COEFS FOR THE 1ST SOIL LAYER
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
* --------------------------------------------------------------------*/
    DELTA[0] = D[0] / B[0];
    for (K = 1; K < *NSOIL; K++)
    {
        P[K] = -C[K] * (1.0 / (B[K] + A[K] * P[K - 1]));
        DELTA[K] =
           (D[K] - A[K] * DELTA[K - 1]) * (1.0 / (B[K] + A[K] * P[K - 1]));
    }

/*----------------------------------------------------------------------
* SET P TO DELTA FOR LOWEST SOIL LAYER
* --------------------------------------------------------------------*/
    P[*NSOIL - 1] = DELTA[*NSOIL - 1];

/*----------------------------------------------------------------------
* ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
* --------------------------------------------------------------------*/
    for (K = 1; K < *NSOIL; K++)
    {
        KK = *NSOIL - K - 1;
        P[KK] = P[KK] * P[KK + 1] + DELTA[KK];
    }

/*----------------------------------------------------------------------
  END SUBROUTINE ROSR12
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
SHFLX (double *SSOIL, double *STC, double *SMC, double *SMCMAX,
   double *SMCMIN, int *NSOIL, double *T1, double *DT, double *YY,
   double *ZZ1, double *ZSOIL, double *TBOT, double *ZBOT, double *SMCWLT,
   double *SH2O, double *VGALPHA, double *VGBETA, double *F1, double *DF1,
   double *QUARTZ, double *CSOIL, int *VEGTYP, int *ISURBAN)
#else
void
SHFLX (double *SSOIL, double *STC, double *SMC, double *SMCMAX, int *NSOIL,
   double *T1, double *DT, double *YY, double *ZZ1, double *ZSOIL,
   double *TBOT, double *ZBOT, double *SMCWLT, double *PSISAT, double *SH2O,
   double *BEXP, double *F1, double *DF1, double *QUARTZ, double *CSOIL,
   int *VEGTYP, int *ISURBAN)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE SHFLX
* ----------------------------------------------------------------------
* UPDATE THE TEMPERATURE STATE OF THE SOIL COLUMN BASED ON THE THERMAL
* DIFFUSION EQUATION AND UPDATE THE FROZEN SOIL MOISTURE CONTENT BASED
* ON THE TEMPERATURE.
* --------------------------------------------------------------------*/
    int             I;

    double          AI[*NSOIL], BI[*NSOIL], CI[*NSOIL], STCF[*NSOIL],
       RHSTS[*NSOIL];

/*----------------------------------------------------------------------
* HRT ROUTINE CALCS THE RIGHT HAND SIDE OF THE SOIL TEMP DIF EQN
* --------------------------------------------------------------------*/

    /*
     * Land case 
     */

#ifdef _FLUX_PIHM_
    HRT (RHSTS, STC, SMC, SMCMAX, SMCMIN, NSOIL, ZSOIL, YY, ZZ1, TBOT, ZBOT,
       SH2O, DT, VGALPHA, VGBETA, F1, DF1, QUARTZ, CSOIL, AI, BI, CI, VEGTYP,
       ISURBAN);
#else
    HRT (RHSTS, STC, SMC, SMCMAX, NSOIL, ZSOIL, YY, ZZ1, TBOT, ZBOT, PSISAT,
       SH2O, DT, BEXP, F1, DF1, QUARTZ, CSOIL, AI, BI, CI, VEGTYP, ISURBAN);
#endif

    HSTEP (STCF, STC, RHSTS, DT, NSOIL, AI, BI, CI);

    for (I = 0; I < *NSOIL; I++)
        STC[I] = STCF[I];

/*----------------------------------------------------------------------
* IN THE NO SNOWPACK CASE (VIA ROUTINE NOPAC BRANCH,) UPDATE THE GRND
* (SKIN) TEMPERATURE HERE IN RESPONSE TO THE UPDATED SOIL TEMPERATURE
* PROFILE ABOVE.  (NOTE: INSPECTION OF ROUTINE SNOPAC SHOWS THAT T1
* BELOW IS A DUMMY VARIABLE ONLY, AS SKIN TEMPERATURE IS UPDATED
* DIFFERENTLY IN ROUTINE SNOPAC)
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* CALCULATE SURFACE SOIL HEAT FLUX
* --------------------------------------------------------------------*/
    *T1 = (*YY + (*ZZ1 - 1.0) * STC[0]) / *ZZ1;
    *SSOIL = *DF1 * (STC[0] - *T1) / (0.5 * ZSOIL[0]);

/*----------------------------------------------------------------------
  END SUBROUTINE SHFLX
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
SMFLX (double *SMC, int *NSOIL, double *CMC, double *DT, double *PRCP1,
   double *PCPDRP, double *ZSOIL, double *SH2O, double *FRZFACT,
   double *SMCMAX, double *SMCMIN, double *VGALPHA, double *VGBETA,
   double *MACKSAT, double *AREAF, int *NMACD, int *MAC_STATUS, int *NWTBL, double *SMCWLT,
   double *DKSAT, double *SHDFAC, double *CMCMAX, double *INF,
   double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC, double *ET,
   double *DRIP)
#else
void
SMFLX (double *SMC, int *NSOIL, double *CMC, double *DT, double *PRCP1,
   double *ZSOIL, double *SH2O, double *SLOPE, double *KDT, double *FRZFACT,
   double *SMCMAX, double *BEXP, double *SMCWLT, double *DKSAT, double *DWSAT,
   double *SHDFAC, double *CMCMAX, double *RUNOFF1, double *RUNOFF2,
   double *RUNOFF3, double *EDIR, double *EC, double *ET, double *DRIP)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE SMFLX
* ----------------------------------------------------------------------
* CALCULATE SOIL MOISTURE FLUX.  THE SOIL MOISTURE CONTENT (SMC - A PER
* UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED WITH
* PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
* FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND
* CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE.
* --------------------------------------------------------------------*/

    int             I;
    double          AI[*NSOIL], BI[*NSOIL], CI[*NSOIL];
    double          RHSTT[*NSOIL];
    double          SICE[*NSOIL];
    double         *DUMMY;
    double          EXCESS;
    double         *RHSCT;
    double          TRHSCT;
#ifndef _FLUX_PIHM_
    double         *PCPDRP;
#endif
    //double          FAC2;
    double         *FLIMIT;
#ifdef _FLUX_PIHM_
    double          KD = 6.54e-7, BFACTR = 3.89;
#endif
    DUMMY = (double *)malloc (sizeof (double));
#ifndef _FLUX_PIHM_
    PCPDRP = (double *)malloc (sizeof (double));
#endif
    RHSCT = (double *)malloc (sizeof (double));
    FLIMIT = (double *)malloc (sizeof (double));

/*----------------------------------------------------------------------
* EXECUTABLE CODE BEGINS HERE.
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* COMPUTE THE RIGHT HAND SIDE OF THE CANOPY EQN TERM ( RHSCT )
* --------------------------------------------------------------------*/
    *DUMMY = 0.;

/*----------------------------------------------------------------------
* CONVERT RHSCT (A RATE) TO TRHSCT (AN AMOUNT) AND ADD IT TO EXISTING
* CMC.  IF RESULTING AMT EXCEEDS MAX CAPACITY, IT BECOMES DRIP AND WILL
* FALL TO THE GRND.
* --------------------------------------------------------------------*/
    *RHSCT = *SHDFAC * *PRCP1 - *EC;
    *DRIP = 0.;
    TRHSCT = *DT * *RHSCT;
    EXCESS = *CMC + TRHSCT;

/*----------------------------------------------------------------------
* PCPDRP IS THE COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE
* SOIL
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_

/*----------------------------------------------------------------------
* PIHM DRIP CALCULATION FOLLOWING RUTTER AND MORTAN (1977 JAE)
* --------------------------------------------------------------------*/
    if (EXCESS > 0)
    {
        if (EXCESS >= *CMCMAX)
        {
            *DRIP = (KD * *CMCMAX * exp (BFACTR)) * *DT + EXCESS - *CMCMAX;
            *RHSCT = *RHSCT - KD * *CMCMAX * exp (BFACTR);
        }
        else
        {
            *DRIP = (KD * *CMCMAX * exp (BFACTR * EXCESS / *CMCMAX)) * *DT;
            *RHSCT = *RHSCT - KD * *CMCMAX * exp (BFACTR * EXCESS / *CMCMAX);
        }
    }
#else
    if (EXCESS > *CMCMAX)
        *DRIP = EXCESS - *CMCMAX;
#endif
    *PCPDRP = (1. - *SHDFAC) * *PRCP1 + *DRIP / *DT;

/*----------------------------------------------------------------------
* STORE ICE CONTENT AT EACH SOIL LAYER BEFORE CALLING SRT and SSTEP
* --------------------------------------------------------------------*/
    for (I = 0; I < *NSOIL; I++)
        SICE[I] = SMC[I] - SH2O[I];

/*----------------------------------------------------------------------
* CALL SUBROUTINES SRT AND SSTEP TO SOLVE THE SOIL MOISTURE
* TENDENCY EQUATIONS.
* IF THE INFILTRATING PRECIP RATE IS NONTRIVIAL,
*   (WE CONSIDER NONTRIVIAL TO BE A PRECIP TOTAL OVER THE TIME STEP
*    EXCEEDING ONE ONE-THOUSANDTH OF THE WATER HOLDING CAPACITY OF
*    THE FIRST SOIL LAYER)
* THEN CALL THE SRT/SSTEP SUBROUTINE PAIR TWICE IN THE MANNER OF
*   TIME SCHEME "F" (IMPLICIT STATE, AVERAGED COEFFICIENT)
*   OF SECTION 2 OF KALNAY AND KANAMITSU (1988, MWR, VOL 116,
*   PAGES 1945-1958)TO MINIMIZE 2-DELTA-T OSCILLATIONS IN THE
*   SOIL MOISTURE VALUE OF THE TOP SOIL LAYER THAT CAN ARISE BECAUSE
*   OF THE EXTREME NONLINEAR DEPENDENCE OF THE SOIL HYDRAULIC
*   DIFFUSIVITY COEFFICIENT AND THE HYDRAULIC CONDUCTIVITY ON THE
*   SOIL MOISTURE STATE
* OTHERWISE CALL THE SRT/SSTEP SUBROUTINE PAIR ONCE IN THE MANNER OF
*   TIME SCHEME "D" (IMPLICIT STATE, EXPLICIT COEFFICIENT)
*   OF SECTION 2 OF KALNAY AND KANAMITSU
* PCPDRP IS UNITS OF KG/M**2/S OR MM/S, ZSOIL IS NEGATIVE DEPTH IN M
* ----------------------------------------------------------------------
*  According to Dr. Ken Mitchell's suggestion, add the second contraint
*  to remove numerical instability of runoff and soil moisture
*  FLIMIT is a limit value for FAC2
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_

/*----------------------------------------------------------------------
* FROZEN GROUND VERSION:
* SMC STATES REPLACED BY SH2O STATES IN SRT SUBR.  SH2O & SICE STATES
* INC&UDED IN SSTEP SUBR.  FROZEN GROUND CORRECTION FACTOR, FRZFACT
* ADDED.  ALL WATER BALANCE CALCULATIONS USING UNFROZEN WATER
* --------------------------------------------------------------------*/
    if (*NWTBL == 0)
    {
        for (I = 0; I < *NSOIL; I++)
        {
            SMC[I] = *SMCMAX;
            SH2O[I] = SMC[I] - SICE[I];
        }
    }
    else
    {

        //      if ((*PCPDRP * *DT) > (0.0001 * 1000.0 * (-ZSOIL[0]) * *SMCMAX))
        //      {
        //          SRT(RHSTT, EDIR, ET, SH2O, SH2O, NWTBL, PCPDRP, ZSOIL, DKSAT, SMCMAX, SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, INF, RUNOFF2, DT, SMCWLT, FRZFACT, SICE, AI, BI, CI);
        //          SSTEP(SH2OFG, SH2O, DUMMY, RHSTT, RHSCT, DT, NWTBL, SMCMAX, SMCMIN, CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
        //          for (K = 0; K < *NSOIL; K++)
        //              SH2OA[K] = (SH2O[K] + SH2OFG[K]) * 0.5;
        //          SRT (RHSTT, EDIR, ET, SH2O, SH2OA, NWTBL, PCPDRP, ZSOIL, DKSAT, SMCMAX, SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, INF, RUNOFF2, DT, SMCWLT, FRZFACT, SICE, AI, BI, CI);
        //          SSTEP (SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NWTBL, SMCMAX, CMCMAX, SMCMIN, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
        //      }
        //      else
        //      {
        //          SRT(RHSTT, EDIR, ET, SH2O, SH2O, NWTBL, PCPDRP, ZSOIL, DKSAT, SMCMAX, SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, INF, RUNOFF2, DT, SMCWLT, FRZFACT, SICE, AI, BI, CI);
        //          SSTEP (SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NWTBL, SMCMAX, SMCMIN, CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
        //      }
        SRT (RHSTT, EDIR, ET, SH2O, SH2O, NSOIL, NWTBL, PCPDRP, ZSOIL, DKSAT, SMCMAX,
           SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, MAC_STATUS, INF, RUNOFF2, DT,
           SMCWLT, FRZFACT, SICE, AI, BI, CI);
        SSTEP (SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NSOIL, SMCMAX, SMCMIN,
           CMCMAX, RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
    }

#else
    FAC2 = 0.0;
    for (I = 0; I < *NSOIL; I++)
    {
        FAC2 = FAC2 > (SH2O[I] / *SMCMAX) ? FAC2 : (SH2O[I] / *SMCMAX);
    }
    FAC2MIT (SMCMAX, FLIMIT);

/*----------------------------------------------------------------------
* FROZEN GROUND VERSION:
* SMC STATES REPLACED BY SH2O STATES IN SRT SUBR.  SH2O & SICE STATES
* INC&UDED IN SSTEP SUBR.  FROZEN GROUND CORRECTION FACTOR, FRZFACT
* ADDED.  ALL WATER BALANCE CALCULATIONS USING UNFROZEN WATER
* --------------------------------------------------------------------*/

    if (((*PCPDRP * *DT) > (0.0001 * 1000.0 * (-ZSOIL[0]) * *SMCMAX))
       || (FAC2 > *FLIMIT))
    {
        SRT (RHSTT, EDIR, ET, SH2O, SH2O, NSOIL, PCPDRP, ZSOIL, DWSAT, DKSAT,
           SMCMAX, BEXP, RUNOFF1, RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT,
           SICE, AI, BI, CI);
        SSTEP (SH2OFG, SH2O, DUMMY, RHSTT, RHSCT, DT, NSOIL, SMCMAX, CMCMAX,
           RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
        for (K = 0; K < *NSOIL; K++)
            SH2OA[K] = (SH2O[K] + SH2OFG[K]) * 0.5;
        SRT (RHSTT, EDIR, ET, SH2O, SH2OA, NSOIL, PCPDRP, ZSOIL, DWSAT, DKSAT,
           SMCMAX, BEXP, RUNOFF1, RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT,
           SICE, AI, BI, CI);
        SSTEP (SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NSOIL, SMCMAX, CMCMAX,
           RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
    }
    else
    {
        SRT (RHSTT, EDIR, ET, SH2O, SH2O, NSOIL, PCPDRP, ZSOIL, DWSAT, DKSAT,
           SMCMAX, BEXP, RUNOFF1, RUNOFF2, DT, SMCWLT, SLOPE, KDT, FRZFACT,
           SICE, AI, BI, CI);
        SSTEP (SH2O, SH2O, CMC, RHSTT, RHSCT, DT, NSOIL, SMCMAX, CMCMAX,
           RUNOFF3, ZSOIL, SMC, SICE, AI, BI, CI);
        //      RUNOF = RUNOFF

    }
#endif
    free (DUMMY);
#ifndef _FLUX_PIHM_
    free (PCPDRP);
#endif
    free (RHSCT);
    free (FLIMIT);

/*----------------------------------------------------------------------
  END SUBROUTINE SMFLX
* --------------------------------------------------------------------*/
}

void
SNFRAC (double *SNEQV, double *SNUP, double *SALP, double *SNOWH,
   double *SNCOVR)
{

/*----------------------------------------------------------------------
* SUBROUTINE SNFRAC
* ----------------------------------------------------------------------
* CALCULATE SNOW FRACTION (0 -> 1)
* SNEQV   SNOW WATER EQUIVALENT (M)
* SNUP    THRESHOLD SNEQV DEPTH ABOVE WHICH SNCOVR=1
* SALP    TUNING PARAMETER
* SNCOVR  FRACTIONAL SNOW COVER
* --------------------------------------------------------------------*/

    double          RSNOW;

/*----------------------------------------------------------------------
* SNUP IS VEG-CLASS DEPENDENT SNOWDEPTH THRESHHOLD (SET IN ROUTINE
* REDPRM) ABOVE WHICH SNOCVR=1.
* --------------------------------------------------------------------*/
    if (*SNEQV < *SNUP)
    {
        RSNOW = *SNEQV / *SNUP;
        *SNCOVR = 1. - (exp (-*SALP * RSNOW) - RSNOW * exp (-*SALP));
    }
    else
        *SNCOVR = 1.0;

    /*
     * FORMULATION OF DICKINSON ET AL. 1986 
     */
    //  Z0N = 0.035;

    //  *SNCOVR = *SNOWH / (*SNOWH + 5. * Z0N);

    /*
     * FORMULATION OF MARSHALL ET AL. 1994 
     */
    //  *SNCOVR = *SNEQV / (*SNEQV + 2. * Z0N);

/*----------------------------------------------------------------------
  END SUBROUTINE SNFRAC
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
SNKSRC (double *TSNSR, double *TAVG, double *SMC, double *SH2O, double *ZSOIL,
   int *NSOIL, double *SMCMAX, double *SMCMIN, double *VGALPHA,
   double *VGBETA, double *DT, int K, double *QTOT)
#else
void
SNKSRC (double *TSNSR, double *TAVG, double *SMC, double *SH2O, double *ZSOIL,
   int *NSOIL, double *SMCMAX, double *PSISAT, double *BEXP, double *DT,
   int K, double *QTOT)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE SNKSRC
* ----------------------------------------------------------------------
* CALCULATE SINK/SOURCE TERM OF THE THERMAL DIFFUSION EQUATION. (SH2O)
* IS AVAILABLE LIQUED WATER.
* --------------------------------------------------------------------*/

    double          DZ, *FREE, XH2O;

    double          DH2O = 1.0000e3, HLICE = 3.3350e5;

    FREE = (double *)malloc (sizeof (double));

    if (K == 0)
        DZ = -ZSOIL[0];
    else
        DZ = ZSOIL[K - 1] - ZSOIL[K];

/*----------------------------------------------------------------------
* VIA FUNCTION FRH2O, COMPUTE POTENTIAL OR 'EQUILIBRIUM' UNFROZEN
* SUPERCOOLED FREE WATER FOR GIVEN SOIL TYPE AND SOIL LAYER TEMPERATURE.
* FUNCTION FRH20 INVOKES EQN (17) FROM V. KOREN ET AL (1999, JGR, VOL.
* 104, PG 19573).  (ASIDE:  LATTER EQN IN JOURNAL IN CENTIGRADE UNITS.
* ROUTINE FRH2O USE FORM OF EQN IN KELVIN UNITS.)
* --------------------------------------------------------------------*/

    /*
     * FREE = FRH2O(TAVG,SMC,SH2O,SMCMAX,BEXP,PSISAT)   
     */

/*----------------------------------------------------------------------
* IN NEXT BLOCK OF CODE, INVOKE EQN 18 OF V. KOREN ET AL (1999, JGR,
* VOL. 104, PG 19573.)  THAT IS, FIRST ESTIMATE THE NEW AMOUNTOF LIQUID
* WATER, 'XH2O', IMPLIED BY THE SUM OF (1) THE LIQUID WATER AT THE BEGIN
* OF CURRENT TIME STEP, AND (2) THE FREEZE OF THAW CHANGE IN LIQUID
* WATER IMPLIED BY THE HEAT FLUX 'QTOT' PASSED IN FROM ROUTINE HRT.
* SECOND, DETERMINE IF XH2O NEEDS TO BE BOUNDED BY 'FREE' (EQUIL AMT) OR
* IF 'FREE' NEEDS TO BE BOUNDED BY XH2O.
* --------------------------------------------------------------------*/
#ifdef _FLUX_PIHM_
    FRH2O (FREE, TAVG, SMC, SH2O, SMCMAX, SMCMIN, VGALPHA, VGBETA);
#else
    FRH2O (FREE, TAVG, SMC, SH2O, SMCMAX, BEXP, PSISAT);
#endif

/*----------------------------------------------------------------------
* FIRST, IF FREEZING AND REMAINING LIQUID LESS THAN LOWER BOUND, THEN
* REDUCE EXTENT OF FREEZING, THEREBY LETTING SOME OR ALL OF HEAT FLUX
* QTOT COOL THE SOIL TEMP LATER IN ROUTINE HRT.
* --------------------------------------------------------------------*/
    XH2O = *SH2O + *QTOT * *DT / (DH2O * HLICE * DZ);

    if (XH2O < *SH2O && XH2O < *FREE)
    {
        if (*FREE > *SH2O)
            XH2O = *SH2O;
        else
            XH2O = *FREE;
    }

/*----------------------------------------------------------------------
* SECOND, IF THAWING AND THE INCREASE IN LIQUID WATER GREATER THAN UPPER
* BOUND, THEN REDUCE EXTENT OF THAW, THEREBY LETTING SOME OR ALL OF HEAT
* FLUX QTOT WARM THE SOIL TEMP LATER IN ROUTINE HRT.
* --------------------------------------------------------------------*/

    if (XH2O > *SH2O && XH2O > *FREE)
    {
        if (*FREE < *SH2O)
            XH2O = *SH2O;
        else
            XH2O = *FREE;
    }

/*----------------------------------------------------------------------
* CALCULATE PHASE-CHANGE HEAT SOURCE/SINK TERM FOR USE IN ROUTINE HRT
* AND UPDATE LIQUID WATER TO REFLCET FINAL FREEZE/THAW INCREMENT.
* --------------------------------------------------------------------*/
    //      SNKSRC = -DH2O*HLICE*DZ*(XH2O-SH2O)/ *DT
    if (XH2O < 0.)
        XH2O = 0.;
    if (XH2O > *SMC)
        XH2O = *SMC;
    *TSNSR = -DH2O * HLICE * DZ * (XH2O - *SH2O) / *DT;

    *SH2O = XH2O;

    free (FREE);

/*----------------------------------------------------------------------
  END SUBROUTINE SNKSRC
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
SNOPAC (double *ETP, double *ETA, double *PRCP, double *PRCPF, double *PCPDRP,
   int *SNOWNG, double *SMC, double *SMCMAX, double *SMCMIN, double *SMCWLT,
   double *SMCREF, double *SMCDRY, double *CMC, double *CMCMAX, int *NSOIL,
   double *DT, double *SBETA, double *DF1, double *Q2, double *T1,
   double *SFCTMP, double *T24, double *TH2, double *FDOWN, double *F1,
   double *SSOIL, double *STC, double *EPSCA, double *SFCPRS, double *VGALPHA,
   double *VGBETA, double *MACKSAT, double *AREAF, int *NMACD, int *MAC_STATUS, int *NWTBL,
   double *PC, double *RCH, double *RR, double *CFACTR, double *SNCOVR,
   double *ESD, double *SNDENS, double *SNOWH, double *SH2O, double *FRZFACT,
   double *ZSOIL, double *DKSAT, double *TBOT, double *ZBOT, double *SHDFAC,
   double *INF, double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC,
   double *ET, double *ETT, int *NROOT, double *SNOMLT, double *RTDIS,
   double *QUARTZ, double *FXEXP, double *CSOIL, double *BETA, double *DRIP,
   double *DEW, double *FLX1, double *FLX2, double *FLX3, double *ESNOW,
   double *ETNS, double *EMISSI, double *RIBB, double *SOLDN, int *ISURBAN,
   int *VEGTYP)
#else
void
SNOPAC (double *ETP, double *ETA, double *PRCP, double *PRCPF, int *SNOWNG,
   double *SMC, double *SMCMAX, double *SMCWLT, double *SMCREF,
   double *SMCDRY, double *CMC, double *CMCMAX, int *NSOIL, double *DT,
   double *SBETA, double *DF1, double *Q2, double *T1, double *SFCTMP,
   double *T24, double *TH2, double *FDOWN, double *F1, double *SSOIL,
   double *STC, double *EPSCA, double *SFCPRS, double *BEXP, double *PC,
   double *RCH, double *RR, double *CFACTR, double *SNCOVR, double *ESD,
   double *SNDENS, double *SNOWH, double *SH2O, double *SLOPE, double *KDT,
   double *FRZFACT, double *PSISAT, double *ZSOIL, double *DWSAT,
   double *DKSAT, double *TBOT, double *ZBOT, double *SHDFAC, double *RUNOFF1,
   double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC, double *ET,
   double *ETT, int *NROOT, double *SNOMLT, double *RTDIS, double *QUARTZ,
   double *FXEXP, double *CSOIL, double *BETA, double *DRIP, double *DEW,
   double *FLX1, double *FLX2, double *FLX3, double *ESNOW, double *ETNS,
   double *EMISSI, double *RIBB, double *SOLDN, int *ISURBAN, int *VEGTYP)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE SNOPAC
* ----------------------------------------------------------------------
* CALCULATE SOIL MOISTURE AND HEAT FLUX VALUES & UPDATE SOIL MOISTURE
* CONTENT AND SOIL HEAT CONTENT VALUES FOR THE CASE WHEN A SNOW PACK IS
* PRESENT.
* --------------------------------------------------------------------*/

    int             K;

    double          ET1[*NSOIL];
    double          DENOM;
    double          DSOIL;
    double          DTOT;
    double          ESNOW1;
    double          ESNOW2;
    //double          ETP2;
    double          ETP3;
    double          ETANRG;
    double          EX;
    double          SEH;
    double          SNCOND;
    double          T12;
    double          T12A;
    double          T12B;
    double          T14;
    double         *EC1;
    double         *EDIR1;
    double         *ETT1;
    double         *ETP1;
    double         *ETNS1;
    double         *PRCP1;
    double         *SSOIL1;
    double         *T11;
    double         *YY;
    double         *ZZ1;
    double          ESDMIN = 1.0e-6;
    double          LSUBC = 2.501e6;
    double          SNOEXP = 2.0;

    EC1 = (double *)malloc (sizeof (double));
    EDIR1 = (double *)malloc (sizeof (double));
    ETT1 = (double *)malloc (sizeof (double));
    ETP1 = (double *)malloc (sizeof (double));
    ETNS1 = (double *)malloc (sizeof (double));
    PRCP1 = (double *)malloc (sizeof (double));
    SSOIL1 = (double *)malloc (sizeof (double));
    T11 = (double *)malloc (sizeof (double));
    YY = (double *)malloc (sizeof (double));
    ZZ1 = (double *)malloc (sizeof (double));

/*----------------------------------------------------------------------
* EXECUTABLE CODE BEGINS HERE:
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* INITIALIZE EVAP TERMS.
* ----------------------------------------------------------------------
* conversions:
* ESNOW [KG M-2 S-1]
* ESDFLX [KG M-2 S-1] .le. ESNOW
* ESNOW1 [M S-1]
* ESNOW2 [M]
* ETP [KG M-2 S-1]
* ETP1 [M S-1]
* ETP2 [M]
* --------------------------------------------------------------------*/
    *DEW = 0.;
    *EDIR = 0.;
    *EDIR1 = 0.;
    *EC1 = 0.;
    *EC = 0.;
    //      EMISSI_S=0.95    ! For snow

    for (K = 0; K < *NSOIL; K++)
    {
        ET[K] = 0.;
        ET1[K] = 0.;
    }
    *ETT = 0.;
    *ETT1 = 0.;
    *ETNS = 0.;
    *ETNS1 = 0.;
    *ESNOW = 0.;
    ESNOW1 = 0.;
    ESNOW2 = 0.;

/*----------------------------------------------------------------------
* CONVERT POTENTIAL EVAP (ETP) FROM KG M-2 S-1 TO ETP1 IN M S-1
* --------------------------------------------------------------------*/
    *PRCP1 = *PRCPF * 0.001;

/*----------------------------------------------------------------------
* IF ETP<0 (DOWNWARD) THEN DEWFALL (=FROSTFALL IN THIS CASE).
* --------------------------------------------------------------------*/
    *BETA = 1.0;
    if (*ETP <= 0.0)
    {
        if ((*RIBB >= 0.1) && (*FDOWN > 150.0))
        {
            *ETP =
               ((*ETP * (1.0 - *RIBB) <
                  0. ? *ETP * (1.0 - *RIBB) : 0.0) * *SNCOVR / 0.980 +
               *ETP * (0.980 - *SNCOVR)) / 0.980;
        }
        if (*ETP == 0.)
            *BETA = 0.0;
        *ETP1 = *ETP * 0.001;
        *DEW = -*ETP1;
        ESNOW2 = *ETP1 * *DT;
        ETANRG = *ETP * ((1. - *SNCOVR) * LSUBC + *SNCOVR * LSUBS);
    }
    else
    {
        *ETP1 = *ETP * 0.001;
        /*
         * LAND CASE 
         */
        if (*SNCOVR < 1.)
        {
#ifdef _FLUX_PIHM_
            EVAPO (ETNS1, SMC, NSOIL, CMC, ETP1, DT, ZSOIL, SH2O, SMCMAX, PC,
               SMCWLT, DKSAT, SMCREF, SHDFAC, CMCMAX, SMCDRY, CFACTR, EDIR1,
               EC1, ET1, ETT1, SFCTMP, Q2, NROOT, RTDIS, FXEXP);
#else
            EVAPO (ETNS1, SMC, NSOIL, CMC, ETP1, DT, ZSOIL, SH2O, SMCMAX,
               BEXP, PC, SMCWLT, DKSAT, DWSAT, SMCREF, SHDFAC, CMCMAX, SMCDRY,
               CFACTR, EDIR1, EC1, ET1, ETT1, SFCTMP, Q2, NROOT, RTDIS,
               FXEXP);
#endif

            /*
             * --------------------------------------------------------------------------
             */
            *EDIR1 = *EDIR1 * (1. - *SNCOVR);
            *EC1 = *EC1 * (1. - *SNCOVR);
            for (K = 0; K < *NSOIL; K++)
                ET1[K] = ET1[K] * (1. - *SNCOVR);
            *ETT1 = *ETT1 * (1. - *SNCOVR);
            //          *ETNS1 = *EDIR1+ *EC1+ *ETT1;
            *ETNS1 = *ETNS1 * (1. - *SNCOVR);

/*--------------------------------------------------------------------------*/
            *EDIR = *EDIR1 * 1000.;
            *EC = *EC1 * 1000.;
            for (K = 0; K < *NSOIL; K++)
                ET[K] = ET1[K] * 1000.;
            *ETT = *ETT1 * 1000.;
            *ETNS = *ETNS1 * 1000.;

/*--------------------------------------------------------------------*/

        }
        *ESNOW = *ETP * *SNCOVR;
        ESNOW1 = *ESNOW * 0.001;
        ESNOW2 = ESNOW1 * *DT;
        ETANRG = *ESNOW * LSUBS + *ETNS * LSUBC;
    }

/*----------------------------------------------------------------------
* IF PRECIP IS FALLING, CALCULATE HEAT FLUX FROM SNOW SFC TO NEWLY
* ACCUMULATING PRECIP.  NOTE THAT THIS REFLECTS THE FLUX APPROPRIATE FOR
* THE NOT-YET-UPDATED SKIN TEMPERATURE (T1).  ASSUMES TEMPERATURE OF THE
* SNOWFALL STRIKING THE GROUND IS =SFCTMP (LOWEST MODEL LEVEL AIR TEMP).
* --------------------------------------------------------------------*/
    *FLX1 = 0.0;
    if (*SNOWNG)
        *FLX1 = CPICE * *PRCP * (*T1 - *SFCTMP);
    else
    {
        if (*PRCP > 0.0)
            *FLX1 = CPH2O * *PRCP * (*T1 - *SFCTMP);

/*----------------------------------------------------------------------
* CALCULATE AN 'EFFECTIVE SNOW-GRND SFC TEMP' (T12) BASED ON HEAT FLUXES
* BETWEEN THE SNOW PACK AND THE SOIL AND ON NET RADIATION.
* INCLUDE FLX1 (PRECIP-SNOW SFC) AND FLX2 (FREEZING RAIN LATENT HEAT)
* FLUXES.  FLX1 FROM ABOVE, FLX2 BROUGHT IN VIA COMMOM BLOCK RITE.
* FLX2 REFLECTS FREEZING RAIN LATENT HEAT FLUX USING T1 CALCULATED IN
* PENMAN.
* --------------------------------------------------------------------*/
    }
    DSOIL = -(0.5 * ZSOIL[0]);
    DTOT = *SNOWH + DSOIL;
    DENOM = 1.0 + *DF1 / (DTOT * *RR * *RCH);

    /*
     * surface emissivity weighted by snow cover fraction
     * *      T12A = ( (*FDOWN - *FLX1 - *FLX2 -                                   &
     * *     &       ((*SNCOVR*EMISSI_S)+*EMISSI*(1.0-*SNCOVR))*SIGMA **T24)/ *RCH    &
     * *     &       + *TH2 - *SFCTMP - ETANRG / *RCH ) / *RR
     */
    T12A =
       ((*FDOWN - *FLX1 - *FLX2 - *EMISSI * SIGMA * *T24) / *RCH + *TH2 -
       *SFCTMP - ETANRG / *RCH) / *RR;
    T12B = *DF1 * STC[0] / (DTOT * *RR * *RCH);

/*----------------------------------------------------------------------
* IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS AT OR BELOW FREEZING, NO SNOW
* MELT WILL OCCUR.  SET THE SKIN TEMP TO THIS EFFECTIVE TEMP.  REDUCE
* (BY SUBLIMINATION ) OR INCREASE (BY FROST) THE DEPTH OF THE SNOWPACK,
* DEPENDING ON SIGN OF ETP.
* UPDATE SOIL HEAT FLUX (SSOIL) USING NEW SKIN TEMPERATURE (T1)
* SINCE NO SNOWMELT, SET ACCUMULATED SNOWMELT TO ZERO, SET 'EFFECTIVE'
* PRECIP FROM SNOWMELT TO ZERO, SET PHASE-CHANGE HEAT FLUX FROM SNOWMELT
* TO ZERO.
* ----------------------------------------------------------------------
* SUB-FREEZING BLOCK
* --------------------------------------------------------------------*/
    T12 = (*SFCTMP + T12A + T12B) / DENOM;

    if (T12 <= TFREEZ)
    {
        *T1 = T12;
        *SSOIL = *DF1 * (*T1 - STC[0]) / DTOT;

        //      *ESD = MAX (0.0, *ESD- ETP2)
        *ESD = *ESD - ESNOW2 > 0. ? (*ESD - ESNOW2) : 0.;
        *FLX3 = 0.0;
        EX = 0.0;

        *SNOMLT = 0.0;
    }

/*----------------------------------------------------------------------
* IF THE 'EFFECTIVE SNOW-GRND SFC TEMP' IS ABOVE FREEZING, SNOW MELT
* WILL OCCUR.  CALL THE SNOW MELT RATE,EX AND AMT, SNOMLT.  REVISE THE
* EFFECTIVE SNOW DEPTH.  REVISE THE SKIN TEMP BECAUSE IT WOULD HAVE CHGD
* DUE TO THE LATENT HEAT RELEASED BY THE MELTING. CALC THE LATENT HEAT
* RELEASED, FLX3. SET THE EFFECTIVE PRECIP, PRCP1 TO THE SNOW MELT RATE,
* EX FOR USE IN SMFLX.  ADJUSTMENT TO T1 TO ACCOUNT FOR SNOW PATCHES.
* CALCULATE QSAT VALID AT FREEZING POINT.  NOTE THAT ESAT (SATURATION
* VAPOR PRESSURE) VALUE OF 6.11E+2 USED HERE IS THAT VALID AT FRZZING
* POINT.  NOTE THAT ETP FROM CALL PENMAN IN SFLX IS IGNORED HERE IN
* FAVOR OF BULK ETP OVER 'OPEN WATER' AT FREEZING TEMP.
* UPDATE SOIL HEAT FLUX (S) USING NEW SKIN TEMPERATURE (T1)
* ----------------------------------------------------------------------
* ABOVE FREEZING BLOCK
* --------------------------------------------------------------------*/
    else
    {
        *T1 =
           TFREEZ * pow (*SNCOVR, SNOEXP) + T12 * (1.0 - pow (*SNCOVR,
              SNOEXP));
        *BETA = 1.0;

/*----------------------------------------------------------------------
* IF POTENTIAL EVAP (SUBLIMATION) GREATER THAN DEPTH OF SNOWPACK.
* BETA<1
* SNOWPACK HAS SUBLIMATED AWAY, SET DEPTH TO ZERO.
* --------------------------------------------------------------------*/
        *SSOIL = *DF1 * (*T1 - STC[0]) / DTOT;

        if (*ESD - ESNOW2 <= ESDMIN)
        {
            *ESD = 0.0;
            EX = 0.0;
            *SNOMLT = 0.0;
            *FLX3 = 0.0;

/*----------------------------------------------------------------------
* SUBLIMATION LESS THAN DEPTH OF SNOWPACK
* SNOWPACK (ESD) REDUCED BY ESNOW2 (DEPTH OF SUBLIMATED SNOW)
* --------------------------------------------------------------------*/
        }
        else
        {
            *ESD = *ESD - ESNOW2;
            ETP3 = *ETP * LSUBC;
            SEH = *RCH * (*T1 - *TH2);
            T14 = *T1 * *T1;
            T14 = T14 * T14;
            //          *FLX3 = *FDOWN - *FLX1 - *FLX2 - ((*SNCOVR*EMISSI_S)+*EMISSI*(1-*SNCOVR))*SIGMA*T14 - *SSOIL - SEH - ETANRG;
            *FLX3 =
               *FDOWN - *FLX1 - *FLX2 - *EMISSI * SIGMA * T14 - *SSOIL - SEH -
               ETANRG;
            if (*FLX3 <= 0.0)
                *FLX3 = 0.0;

/*----------------------------------------------------------------------
* SNOWMELT REDUCTION DEPENDING ON SNOW COVER
* --------------------------------------------------------------------*/
            EX = *FLX3 * 0.001 / LSUBF;

/*----------------------------------------------------------------------
* ESDMIN REPRESENTS A SNOWPACK DEPTH THRESHOLD VALUE BELOW WHICH WE
* CHOOSE NOT TO RETAIN ANY SNOWPACK, AND INSTEAD INCLUDE IT IN SNOWMELT.
* --------------------------------------------------------------------*/
            *SNOMLT = EX * *DT;
            if (*ESD - *SNOMLT >= ESDMIN)
            {
                *ESD = *ESD - *SNOMLT;

/*----------------------------------------------------------------------
* SNOWMELT EXCEEDS SNOW DEPTH
* --------------------------------------------------------------------*/
            }
            else
            {
                EX = *ESD / *DT;
                *FLX3 = EX * 1000.0 * LSUBF;
                *SNOMLT = *ESD;

                *ESD = 0.0;

/*----------------------------------------------------------------------
* END OF 'ESD .LE. ETP2' IF-BLOCK
----------------------------------------------------------------------*/
            }
        }

/*----------------------------------------------------------------------
* END OF 'T12 .LE. TFREEZ' IF-BLOCK
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* IF NON-GLACIAL LAND, ADD SNOWMELT RATE (EX) TO PRECIP RATE TO BE USED
* IN SUBROUTINE SMFLX (SOIL MOISTURE EVOLUTION) VIA INFILTRATION.
*
* RUNOFF/BASEFLOW LATER NEAR THE END OF SFLX (AFTER RETURN FROM CALL TO
* SUBROUTINE SNOPAC)
* --------------------------------------------------------------------*/
        *PRCP1 = *PRCP1 + EX;

/*----------------------------------------------------------------------
* SET THE EFFECTIVE POTNL EVAPOTRANSP (ETP1) TO ZERO SINCE THIS IS SNOW
* CASE, SO SURFACE EVAP NOT CALCULATED FROM EDIR, EC, OR ETT IN SMFLX
* (BELOW).
* SMFLX RETURNS UPDATED SOIL MOISTURE VALUES FOR NON-GLACIAL LAND.
* --------------------------------------------------------------------*/
    }
#ifdef _FLUX_PIHM_
    SMFLX (SMC, NSOIL, CMC, DT, PRCP1, PCPDRP, ZSOIL, SH2O, FRZFACT, SMCMAX,
       SMCMIN, VGALPHA, VGBETA, MACKSAT, AREAF, NMACD, MAC_STATUS, NWTBL, SMCWLT, DKSAT,
       SHDFAC, CMCMAX, INF, RUNOFF2, RUNOFF3, EDIR1, EC1, ET1, DRIP);
#else
    SMFLX (SMC, NSOIL, CMC, DT, PRCP1, ZSOIL, SH2O, SLOPE, KDT, FRZFACT,
       SMCMAX, BEXP, SMCWLT, DKSAT, DWSAT, SHDFAC, CMCMAX, RUNOFF1, RUNOFF2,
       RUNOFF3, EDIR1, EC1, ET1, DRIP);
#endif

/*----------------------------------------------------------------------
* BEFORE CALL SHFLX IN THIS SNOWPACK CASE, SET ZZ1 AND YY ARGUMENTS TO
* SPECIAL VALUES THAT ENSURE THAT GROUND HEAT FLUX CALCULATED IN SHFLX
* MATCHES THAT ALREADY COMPUTER FOR BELOW THE SNOWPACK, THUS THE SFC
* HEAT FLUX TO BE COMPUTED IN SHFLX WILL EFFECTIVELY BE THE FLUX AT THE
* SNOW TOP SURFACE.  T11 IS A DUMMY ARGUEMENT SO WE WILL NOT USE THE
* SKIN TEMP VALUE AS REVISED BY SHFLX.
* --------------------------------------------------------------------*/
    *ZZ1 = 1.0;
    *YY = STC[0] - 0.5 * *SSOIL * ZSOIL[0] * *ZZ1 / *DF1;

/*----------------------------------------------------------------------
* SHFLX WILL CALC/UPDATE THE SOIL TEMPS.  NOTE:  THE SUB-SFC HEAT FLUX
* (SSOIL1) AND THE SKIN TEMP (T11) OUTPUT FROM THIS SHFLX CALL ARE NOT
* USED  IN ANY SUBSEQUENT CALCULATIONS. RATHER, THEY ARE DUMMY VARIABLES
* HERE IN THE SNOPAC CASE, SINCE THE SKIN TEMP AND SUB-SFC HEAT FLUX ARE
* UPDATED INSTEAD NEAR THE BEGINNING OF THE CALL TO SNOPAC.
* --------------------------------------------------------------------*/
    *T11 = *T1;
#ifdef _FLUX_PIHM_
    SHFLX (SSOIL1, STC, SMC, SMCMAX, SMCMIN, NSOIL, T11, DT, YY, ZZ1, ZSOIL,
       TBOT, ZBOT, SMCWLT, SH2O, VGALPHA, VGBETA, F1, DF1, QUARTZ, CSOIL,
       VEGTYP, ISURBAN);
#else
    SHFLX (SSOIL1, STC, SMC, SMCMAX, NSOIL, T11, DT, YY, ZZ1, ZSOIL, TBOT,
       ZBOT, SMCWLT, PSISAT, SH2O, BEXP, F1, DF1, QUARTZ, CSOIL, VEGTYP,
       ISURBAN);
#endif

/*----------------------------------------------------------------------
* SNOW DEPTH AND DENSITY ADJUSTMENT BASED ON SNOW COMPACTION.  YY IS
* ASSUMED TO BE THE SOIL TEMPERTURE AT THE TOP OF THE SOIL COLUMN.
* --------------------------------------------------------------------*/
    /*
     * LAND 
     */
    if (*ESD > 0.)
        SNOWPACK (ESD, DT, SNOWH, SNDENS, T1, YY);
    else
    {
        *ESD = 0.;
        *SNOWH = 0.;
        *SNDENS = 0.;
        SNCOND = 1.;
        *SNCOVR = 0.;
    }

    free (EC1);
    free (EDIR1);
    free (ETT1);
    free (ETP1);
    free (ETNS1);
    free (PRCP1);
    free (SSOIL1);
    free (T11);
    free (YY);
    free (ZZ1);

/*----------------------------------------------------------------------
  END SUBROUTINE SNOPAC
* --------------------------------------------------------------------*/
}

void SNOWPACK (double *ESD, double *DTSEC, double *SNOWH, double *SNDENS, double *TSNOW, double *TSOIL)
{
/*----------------------------------------------------------------------
* SUBROUTINE SNOWPACK
* ----------------------------------------------------------------------
* CALCULATE COMPACTION OF SNOWPACK UNDER CONDITIONS OF INCREASING SNOW
* DENSITY, AS OBTAINED FROM AN APPROXIMATE SOLUTION OF E. ANDERSON'S
* DIFFERENTIAL EQUATION (3.29), NOAA TECHNICAL REPORT NWS 19, BY VICTOR
* KOREN, 03/25/95.
* ----------------------------------------------------------------------
* ESD     WATER EQUIVALENT OF SNOW (M)
* DTSEC   TIME STEP (SEC)
* SNOWH   SNOW DEPTH (M)
* SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
* TSNOW   SNOW SURFACE TEMPERATURE (K)
* TSOIL   SOIL SURFACE TEMPERATURE (K)

* SUBROUTINE WILL RETURN NEW VALUES OF SNOWH AND SNDENS
* --------------------------------------------------------------------*/

    int             IPOL;
    int             J;
    double          BFAC;
    double          DSX;
    double          DTHR;
    double          DW;
    double          SNOWHC;
    double          PEXP;
    double          TAVGC;
    double          TSNOWC;
    double          TSOILC;
    double          ESDC;
    double          ESDCX;
    double          C1 = 0.01;
    double          C2 = 21.0;

/*----------------------------------------------------------------------
* CONVERSION INTO SIMULATION UNITS
* --------------------------------------------------------------------*/
    SNOWHC = *SNOWH * 100.0;
    ESDC = *ESD * 100.0;
    DTHR = *DTSEC / 3600.0;
    TSNOWC = *TSNOW - 273.15;
    TSOILC = *TSOIL - 273.15;

/*----------------------------------------------------------------------
* CALCULATING OF AVERAGE TEMPERATURE OF SNOW PACK
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* CALCULATING OF SNOW DEPTH AND DENSITY AS A RESULT OF COMPACTION
*  SNDENS=DS0*(EXP(BFAC*ESD)-1.)/(BFAC*ESD)
*  BFAC=DTHR*C1*EXP(0.08*TAVGC-C2*DS0)
* NOTE: BFAC*ESD IN SNDENS EQN ABOVE HAS TO BE CAREFULLY TREATED
* NUMERICALLY BELOW:
*   C1 IS THE FRACTIONAL INCREASE IN DENSITY (1/(CM*HR))
*   C2 IS A CONSTANT (CM3/G) KOJIMA ESTIMATED AS 21 CMS/G
* --------------------------------------------------------------------*/
    TAVGC = 0.5 * (TSNOWC + TSOILC);
    if (ESDC > 1.e-2)
        ESDCX = ESDC;
    else
        ESDCX = 1.e-2;

    //      DSX = *SNDENS*((DEXP(BFAC*ESDC)-1.)/(BFAC*ESDC))

/*----------------------------------------------------------------------
* THE FUNCTION OF THE FORM (e**x-1)/x EMBEDDED IN ABOVE EXPRESSION
* FOR DSX WAS CAUSING NUMERICAL DIFFICULTIES WHEN THE DENOMINATOR "x"
* (I.E. BFAC*ESDC) BECAME ZERO OR APPROACHED ZERO (DESPITE THE FACT THAT
* THE ANALYTICAL FUNCTION (e**x-1)/x HAS A WELL DEFINED LIMIT AS
* "x" APPROACHES ZERO), HENCE BELOW WE REPLACE THE (e**x-1)/x
* EXPRESSION WITH AN EQUIVALENT, NUMERICALLY WELL-BEHAVED
* POLYNOMIAL EXPANSION.

* NUMBER OF TERMS OF POLYNOMIAL EXPANSION, AND HENCE ITS ACCURACY,
* IS GOVERNED BY ITERATION LIMIT "IPOL".
*      IPOL GREATER THAN 9 ONLY MAKES A DIFFERENCE ON DOUBLE
*            PRECISION (RELATIVE ERRORS GIVEN IN PERCENT %).
*       IPOL=9, FOR REL.ERROR <~ 1.6 E-6 % (8 SIGNIFICANT DIGITS)
*       IPOL=8, FOR REL.ERROR <~ 1.8 E-5 % (7 SIGNIFICANT DIGITS)
*       IPOL=7, FOR REL.ERROR <~ 1.8 E-4 % ...
* --------------------------------------------------------------------*/
    BFAC = DTHR * C1 * exp (0.08 * TAVGC - C2 * *SNDENS);
    IPOL = 4;
    PEXP = 0.;
    //  PEXP = (1. + PEXP)*BFAC*ESDC/REAL(J+1);
    for (J = IPOL; J > 0; J--)
    {
        PEXP = (1. + PEXP) * BFAC * ESDCX / (double)(J + 1);
    }

    PEXP = PEXP + 1.;

/*----------------------------------------------------------------------
* ABOVE LINE ENDS POLYNOMIAL SUBSTITUTION
* ----------------------------------------------------------------------
*     END OF KOREAN FORMULATION
* --------------------------------------------------------------------*/

    /*
     *     BASE FORMULATION (COGLEY ET AL., 1990)
     *     CONVERT DENSITY FROM G/CM3 TO KG/M3
     *       DSM=*SNDENS*1000.0
     
     *       DSX=DSM+*DTSEC*0.5*DSM*G* *ESD/
     *    &      (1E7*EXP(-0.02*DSM+KN/(TAVGC+273.16)-14.643))
     
     *  &   CONVERT DENSITY FROM KG/M3 TO G/CM3
     *       DSX=DSX/1000.0
     
     *     END OF COGLEY ET AL. FORMULATION
     */

/*----------------------------------------------------------------------
* SET UPPER/LOWER LIMIT ON SNOW DENSITY
* --------------------------------------------------------------------*/
    DSX = *SNDENS * (PEXP);
    if (DSX > 0.40)
        DSX = 0.40;
    if (DSX < 0.05)
        DSX = 0.05;

/*----------------------------------------------------------------------
* UPDATE OF SNOW DEPTH AND DENSITY DEPENDING ON LIQUID WATER DURING
* SNOWMELT.  ASSUMED THAT 13% OF LIQUID WATER CAN BE STORED IN SNOW PER
* DAY DURING SNOWMELT TILL SNOW DENSITY 0.40.
* --------------------------------------------------------------------*/
    *SNDENS = DSX;
    if (TSNOWC >= 0.)
    {
        DW = 0.13 * DTHR / 24.;
        *SNDENS = *SNDENS * (1. - DW) + DW;
        if (*SNDENS >= 0.40)
            *SNDENS = 0.40;

/*----------------------------------------------------------------------
* CALCULATE SNOW DEPTH (CM) FROM SNOW WATER EQUIVALENT AND SNOW DENSITY.
* CHANGE SNOW DEPTH UNITS TO METERS
* --------------------------------------------------------------------*/
    }
    SNOWHC = ESDC / *SNDENS;
    *SNOWH = SNOWHC * 0.01;

/*----------------------------------------------------------------------
  END SUBROUTINE SNOWPACK
* --------------------------------------------------------------------*/
}

void SNOWZ0 (double *SNCOVR, double *Z0, double *Z0BRD, double *SNOWH)
{

/*----------------------------------------------------------------------
* SUBROUTINE SNOWZ0
* ----------------------------------------------------------------------
* CALCULATE TOTAL ROUGHNESS LENGTH OVER SNOW
* SNCOVR  FRACTIONAL SNOW COVER
* Z0      ROUGHNESS LENGTH (m)
* Z0S     SNOW ROUGHNESS LENGTH:=0.001 (m)
* --------------------------------------------------------------------*/
    double          Z0S = 0.001;
    double          BURIAL;
    double          Z0EFF;

    //m Z0 = (1.- SNCOVR)* Z0BRD + SNCOVR * Z0S
    BURIAL = 7.0 * *Z0BRD - *SNOWH;
    if (BURIAL < 0.0007)
        Z0EFF = Z0S;
    else
        Z0EFF = BURIAL / 7.0;

    *Z0 = (1. - *SNCOVR) * *Z0BRD + *SNCOVR * Z0EFF;

/*----------------------------------------------------------------------
  END SUBROUTINE SNOWZ0
* --------------------------------------------------------------------*/
}

void SNOW_NEW (double *TEMP, double *NEWSN, double *SNOWH, double *SNDENS)
{

/*----------------------------------------------------------------------
* SUBROUTINE SNOW_NEW
* ----------------------------------------------------------------------
* CALCULATE SNOW DEPTH AND DENSITY TO ACCOUNT FOR THE NEW SNOWFALL.
* NEW VALUES OF SNOW DEPTH & DENSITY RETURNED.

* TEMP    AIR TEMPERATURE (K)
* NEWSN   NEW SNOWFALL (M)
* SNOWH   SNOW DEPTH (M)
* SNDENS  SNOW DENSITY (G/CM3=DIMENSIONLESS FRACTION OF H2O DENSITY)
* --------------------------------------------------------------------*/

    double          DSNEW, HNEWC, SNOWHC, NEWSNC, TEMPC;

/*----------------------------------------------------------------------
* CONVERSION INTO SIMULATION UNITS
* --------------------------------------------------------------------*/
    SNOWHC = *SNOWH * 100.0;
    NEWSNC = *NEWSN * 100.0;

/*----------------------------------------------------------------------
* CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE
* EQUATION FROM GOTTLIB L. 'A GENERAL RUNOFF MODEL FOR SNOWCOVERED
* AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE,
* VEMADOLEN, SWEDEN, 1980, 172-177PP.
*---------------------------------------------------------------------*/
    TEMPC = *TEMP - 273.15;
    if (TEMPC <= -15.)
        DSNEW = 0.05;
    else
        DSNEW = 0.05 + 0.0017 * pow (TEMPC + 15., 1.5);

/*----------------------------------------------------------------------
* ADJUSTMENT OF SNOW DENSITY DEPENDING ON NEW SNOWFALL
* --------------------------------------------------------------------*/
    HNEWC = NEWSNC / DSNEW;
    if (SNOWHC + HNEWC < 0.001)
        *SNDENS = DSNEW > *SNDENS ? DSNEW : *SNDENS;
    else
        *SNDENS = (SNOWHC * *SNDENS + HNEWC * DSNEW) / (SNOWHC + HNEWC);
    SNOWHC = SNOWHC + HNEWC;
    *SNOWH = SNOWHC * 0.01;

/*----------------------------------------------------------------------
  END SUBROUTINE SNOW_NEW
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
SRT (double *RHSTT, double *EDIR, double *ET, double *SH2O, double *SH2OA,
   int *NSOIL, int *NWTBL, double *PCPDRP, double *ZSOIL, double *DKSAT, double *SMCMAX,
   double *SMCMIN, double *VGALPHA, double *VGBETA, double *MACKSAT,
   double *AREAF, int *NMACD, int *MAC_STATUS, double *INF, double *RUNOFF2, double *DT,
   double *SMCWLT, double *FRZX, double *SICE, double *AI, double *BI,
   double *CI)
#else
void
SRT (double *RHSTT, double *EDIR, double *ET, double *SH2O, double *SH2OA,
   int *NSOIL, double *PCPDRP, double *ZSOIL, double *DWSAT, double *DKSAT,
   double *SMCMAX, double *BEXP, double *RUNOFF1, double *RUNOFF2, double *DT,
   double *SMCWLT, double *SLOPE, double *KDT, double *FRZX, double *SICE,
   double *AI, double *BI, double *CI)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE SRT
* ----------------------------------------------------------------------
* CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
* WATER DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
* COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
* --------------------------------------------------------------------*/
    int             IOHINF;
    int             K, KS;
    double          DDZ;
    double          DDZ2;
    double          DENOM;
    double          DENOM2;
    double          NUMER;
    double          PDDUM;
    double          SSTT;
    double         *MXSMC, *MXSMC2;
    double         *SICEMAX;
    double         *WCND;
    double         *WCND2;
    double         *WDF;
    double         *WDF2;
#ifdef _FLUX_PIHM_
    double         *DSMDZ, *DSMDZ2;
    int             MACPORE[*NSOIL];
#else
    double          DSMDZ, DSMDZ2;
#endif

/*----------------------------------------------------------------------
* FROZEN GROUND VERSION:
* REFERENCE FROZEN GROUND PARAMETER, CVFRZ, IS A SHAPE PARAMETER OF
* AREAL DISTRIBUTION FUNCTION OF SOIL ICE CONTENT WHICH EQUALS 1/CV.
* CV IS A COEFFICIENT OF SPATIAL VARIATION OF SOIL ICE CONTENT.  BASED
* ON FIELD DATA CV DEPENDS ON AREAL MEAN OF FROZEN DEPTH, AND IT CLOSE
* TO CONSTANT = 0.6 IF AREAL MEAN FROZEN DEPTH IS ABOVE 20 CM.  THAT IS
* WHY PARAMETER CVFRZ = 3 (INT{1/0.6*0.6}).
* CURRENT LOGIC DOESN'T ALLOW CVFRZ BE BIGGER THAN 3
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF.  INCLUDE THE
* INFILTRATION FORMULE FROM SCHAAKE AND KOREN MODEL.
* MODIFIED BY Q DUAN
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* LET SICEMAX BE THE GREATEST, IF ANY, FROZEN WATER CONTENT WITHIN SOIL
* LAYERS.
* --------------------------------------------------------------------*/
    SICEMAX = (double *)malloc (sizeof (double));
    WCND = (double *)malloc (sizeof (double));
    WCND2 = (double *)malloc (sizeof (double));
    WDF = (double *)malloc (sizeof (double));
    WDF2 = (double *)malloc (sizeof (double));
    IOHINF = 1;
    *SICEMAX = 0.0;

    for (KS = 0; KS < *NSOIL; KS++)
    {
        if (SICE[KS] > *SICEMAX)
            *SICEMAX = SICE[KS];
/*----------------------------------------------------------------------
* DETERMINE RAINFALL INFILTRATION RATE AND RUNOFF
* --------------------------------------------------------------------*/
    }

#ifdef _FLUX_PIHM_
    DSMDZ = (double *)malloc (sizeof (double));
    DSMDZ2 = (double *)malloc (sizeof (double));

    PDDUM = *INF;

    for (K = 0; K < *NSOIL; K++)
        MACPORE[K] = 0;
    for (K = 0; K < *NMACD - 1; K++)
        MACPORE[K] = 1;

///*----------------------------------------------------------------------
//* YS: IF LATERAL RUNOFF (RUNOFF2) IS NEGATIVE (GRID IS A SINK) AND THE
//* INTERFACE LAYER IS CLOSE TO SATURATION, LATERAL RUNOFF IS ADDED TO THE
//* LAYER ABOVE
//* --------------------------------------------------------------------*/
//    DENOM2 = ZSOIL[*NSOIL - 2] - ZSOIL[*NSOIL - 1];
//    if (*NSOIL == 2)
//        DENOM = -ZSOIL[*NSOIL - 1];
//    else
//        DENOM = ZSOIL[*NSOIL - 3] - ZSOIL[*NSOIL - 1];
//
//    *DSMDZ = (SH2O[*NSOIL - 2] - SH2O[*NSOIL - 1]) / (DENOM * 0.5);
//
//    WDFCND (WDF, WCND, SH2OA + *NSOIL - 2, SMCMAX, SMCMIN, VGALPHA, VGBETA,
//       DKSAT, MACKSAT, AREAF, MAC_STATUS, SICEMAX, DSMDZ, MACPORE + *NSOIL - 2);
//
//    if (*RUNOFF2 < 0
//       && (*SMCMAX - SH2O[*NSOIL - 1]) / *DT * DENOM2 <
//       *WDF * *DSMDZ + *WCND - ET[*NSOIL - 1] - *RUNOFF2)
//        RFLAG = 1;
//    else
//        RFLAG = 0;

    MXSMC = SH2OA;

    *DSMDZ = (SH2O[0] - SH2O[1]) / (-0.5 * ZSOIL[1]);
    WDFCND (WDF, WCND, MXSMC, SMCMAX, SMCMIN, VGALPHA, VGBETA, DKSAT, MACKSAT,
       AREAF, MAC_STATUS, SICEMAX, DSMDZ, MACPORE);

/*----------------------------------------------------------------------
* CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
* --------------------------------------------------------------------*/
    DDZ = 1. / (-.5 * ZSOIL[1]);
    AI[0] = 0.0;
    BI[0] = *WDF * DDZ / (-ZSOIL[0]);

/*----------------------------------------------------------------------
* CALC RHSTT FOR THE TOP LAYER AFTER CALC'NG THE VERTICAL SOIL MOISTURE
* GRADIENT BTWN THE TOP AND NEXT TO TOP LAYERS.
* --------------------------------------------------------------------*/
    CI[0] = -BI[0];
    RHSTT[0] = (*WDF * *DSMDZ + *WCND - PDDUM + *EDIR + ET[0]) / ZSOIL[0];

    if (*NWTBL == 1)
        RHSTT[0] += *RUNOFF2 / ZSOIL[0];

/*----------------------------------------------------------------------
* INITIALIZE DDZ2
* --------------------------------------------------------------------*/
    SSTT = *WDF * *DSMDZ + *WCND + *EDIR + ET[0];

/*----------------------------------------------------------------------
* LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
* --------------------------------------------------------------------*/
    DDZ2 = 0.0;
    for (K = 1; K < *NSOIL; K++)
    {
        DENOM2 = (ZSOIL[K - 1] - ZSOIL[K]);
        if (K < *NSOIL - 1)
        {

/*----------------------------------------------------------------------
* AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN
* LINE BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
* 'MXSMC2 = MAX (SH2OA(K), SH2OA(K+1))'
* --------------------------------------------------------------------*/
            MXSMC2 = SH2OA + K;
            DENOM = (ZSOIL[K - 1] - ZSOIL[K + 1]);
            *DSMDZ2 = (SH2O[K] - SH2O[K + 1]) / (DENOM * 0.5);
            WDFCND (WDF2, WCND2, MXSMC2, SMCMAX, SMCMIN, VGALPHA, VGBETA,
               DKSAT, MACKSAT, AREAF, MAC_STATUS, SICEMAX, DSMDZ2, MACPORE + K);

/*-----------------------------------------------------------------------
* CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALC'NG RHSTT
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
* --------------------------------------------------------------------*/
            DDZ2 = 2.0 / DENOM;
            CI[K] = -*WDF2 * DDZ2 / DENOM2;
        }
        else
        {

/*----------------------------------------------------------------------
* SLOPE OF BOTTOM LAYER IS INTRODUCED
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR
* THIS LAYER
* --------------------------------------------------------------------*/
            *WDF2 = 0;
            *WCND2 = 0;

/*----------------------------------------------------------------------
* CALC A PARTIAL PRODUCT FOR LATER USE IN CALC'NG RHSTT
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* SET MATRIX COEF CI TO ZERO
* --------------------------------------------------------------------*/
            *DSMDZ2 = 0.0;
            CI[K] = 0.0;

/*----------------------------------------------------------------------
* CALC RHSTT FOR THIS LAYER AFTER CALC'NG ITS NUMERATOR
* --------------------------------------------------------------------*/
        }

        NUMER = (*WDF2 * *DSMDZ2) + *WCND2 - (*WDF * *DSMDZ) - *WCND + ET[K];
        if (K == *NWTBL - 1)
                NUMER = NUMER + *RUNOFF2;

/*----------------------------------------------------------------------
* CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
* --------------------------------------------------------------------*/
        RHSTT[K] = NUMER / (-DENOM2);
        AI[K] = -*WDF * DDZ / DENOM2;

/*----------------------------------------------------------------------
* RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
* RUNOFF2:  SUB-SURFACE OR BASEFLOW RUNOFF
* --------------------------------------------------------------------*/
        BI[K] = -(AI[K] + CI[K]);

        if (K != *NSOIL - 1)
        {
            *WDF = *WDF2;
            *WCND = *WCND2;
            *DSMDZ = *DSMDZ2;
            DDZ = DDZ2;
        }
    }

    free (DSMDZ);
    free (DSMDZ2);
#else
    PDDUM = *PCPDRP;
    *RUNOFF1 = 0.0;

/*----------------------------------------------------------------------
* MODIFIED BY Q. DUAN, 5/16/94
* --------------------------------------------------------------------*/
    //        IF (IOHINF == 1) THEN

    if (*PCPDRP != 0.0)
    {
        DT1 = *DT / 86400.;
        SMCAV = *SMCMAX - *SMCWLT;

/*----------------------------------------------------------------------
* FROZEN GROUND VERSION:
* --------------------------------------------------------------------*/
        DMAX[0] = -ZSOIL[0] * SMCAV;

        DICE = -ZSOIL[0] * SICE[0];
        DMAX[0] = DMAX[0] * (1.0 - (SH2OA[0] + SICE[0] - *SMCWLT) / SMCAV);

        DD = DMAX[0];

/*----------------------------------------------------------------------
* FROZEN GROUND VERSION:
* --------------------------------------------------------------------*/
        for (KS = 1; KS < *NSOIL; KS++)
        {
            DICE = DICE + (ZSOIL[KS - 1] - ZSOIL[KS]) * SICE[KS];
            DMAX[KS] = (ZSOIL[KS - 1] - ZSOIL[KS]) * SMCAV;
            DMAX[KS] =
               DMAX[KS] * (1.0 - (SH2OA[KS] + SICE[KS] - *SMCWLT) / SMCAV);
            DD = DD + DMAX[KS];

/*----------------------------------------------------------------------
* VAL = (1.-EXP(-KDT*SQRT(DT1)))
* IN BELOW, REMOVE THE SQRT IN ABOVE
* --------------------------------------------------------------------*/
        }
        VAL = (1. - exp (-*KDT * DT1));
        DDT = DD * VAL;
        PX = *PCPDRP * *DT;
        if (PX < 0.0)
            PX = 0.0;

/*----------------------------------------------------------------------
* FROZEN GROUND VERSION:
* REDUCTION OF INFILTRATION BASED ON FROZEN GROUND PARAMETERS
* --------------------------------------------------------------------*/
        INFMAX = (PX * (DDT / (PX + DDT))) / *DT;
        FCR = 1.;
        if (DICE > 1.e-2)
        {
            ACRT = (double)CVFRZ **FRZX / DICE;
            SUM = 1.;
            IALP1 = CVFRZ - 1;
            for (J = 1; J < IALP1 + 1; J++)
            {
                K = 1;
                for (JJ = J + 1; JJ < IALP1; JJ++)
                {
                    K = K * JJ;
                }
                SUM = SUM + pow (ACRT, (double)(CVFRZ - J)) / (double)K;
            }
            FCR = 1. - exp (-ACRT) * SUM;
        }

/*----------------------------------------------------------------------
* CORRECTION OF INFILTRATION LIMITATION:
* IF INFMAX .LE. HYDROLIC CONDUCTIVITY ASSIGN INFMAX THE VALUE OF
* HYDROLIC CONDUCTIVITY
* ---------------------------------------------------------------------*/
        //         MXSMC = MAX ( SH2OA(1), SH2OA(2) )
        INFMAX = INFMAX * FCR;
        MXSMC = SH2OA;
        WDFCND (WDF, WCND, MXSMC, SMCMAX, BEXP, DKSAT, DWSAT, SICEMAX);
        INFMAX = INFMAX > *WCND ? INFMAX : *WCND;

        INFMAX = INFMAX < PX / *DT ? INFMAX : PX / *DT;
        if (*PCPDRP > INFMAX)
        {
            *RUNOFF1 = *PCPDRP - INFMAX;
            PDDUM = INFMAX;
        }

/*----------------------------------------------------------------------
* TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN LINE
* BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
* 'MXSMC = MAX(SH2OA(1), SH2OA(2))'
* --------------------------------------------------------------------*/
    }

    MXSMC = SH2OA;
    WDFCND (WDF, WCND, MXSMC, SMCMAX, BEXP, DKSAT, DWSAT, SICEMAX);

/*----------------------------------------------------------------------
* CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
* --------------------------------------------------------------------*/
    DDZ = 1. / (-.5 * ZSOIL[1]);
    AI[0] = 0.0;
    BI[0] = *WDF * DDZ / (-ZSOIL[0]);

/*----------------------------------------------------------------------
* CALC RHSTT FOR THE TOP LAYER AFTER CALC'NG THE VERTICAL SOIL MOISTURE
* GRADIENT BTWN THE TOP AND NEXT TO TOP LAYERS.
* --------------------------------------------------------------------*/
    CI[0] = -BI[0];
    DSMDZ = (SH2O[0] - SH2O[1]) / (-0.5 * ZSOIL[1]);
    RHSTT[0] = (*WDF * DSMDZ + *WCND - PDDUM + *EDIR + ET[0]) / ZSOIL[0];

/*----------------------------------------------------------------------
* INITIALIZE DDZ2
* --------------------------------------------------------------------*/
    SSTT = *WDF * DSMDZ + *WCND + *EDIR + ET[0];

/*----------------------------------------------------------------------
* LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
* --------------------------------------------------------------------*/
    DDZ2 = 0.0;
    for (K = 1; K < *NSOIL; K++)
    {
        DENOM2 = (ZSOIL[K - 1] - ZSOIL[K]);
        if (K != *NSOIL - 1)
        {

/*----------------------------------------------------------------------
* AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN
* LINE BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
* 'MXSMC2 = MAX (SH2OA(K), SH2OA(K+1))'
* --------------------------------------------------------------------*/
            SLOPX = 1.;

            MXSMC2 = SH2OA + K;
            WDFCND (WDF2, WCND2, MXSMC2, SMCMAX, BEXP, DKSAT, DWSAT, SICEMAX);

/*-----------------------------------------------------------------------
* CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALC'NG RHSTT
* --------------------------------------------------------------------*/
            DENOM = (ZSOIL[K - 1] - ZSOIL[K + 1]);

/*----------------------------------------------------------------------
* CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
* --------------------------------------------------------------------*/
            DSMDZ2 = (SH2O[K] - SH2O[K + 1]) / (DENOM * 0.5);
            DDZ2 = 2.0 / DENOM;
            CI[K] = -*WDF2 * DDZ2 / DENOM2;
        }
        else
        {

/*----------------------------------------------------------------------
* SLOPE OF BOTTOM LAYER IS INTRODUCED
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR
* THIS LAYER
* --------------------------------------------------------------------*/
            SLOPX = *SLOPE;
            WDFCND (WDF2, WCND2, SH2OA + *NSOIL - 1, SMCMAX, BEXP, DKSAT,
               DWSAT, SICEMAX);

/*----------------------------------------------------------------------
* CALC A PARTIAL PRODUCT FOR LATER USE IN CALC'NG RHSTT
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* SET MATRIX COEF CI TO ZERO
* --------------------------------------------------------------------*/
            DSMDZ2 = 0.0;
            CI[K] = 0.0;

/*----------------------------------------------------------------------
* CALC RHSTT FOR THIS LAYER AFTER CALC'NG ITS NUMERATOR
* --------------------------------------------------------------------*/
        }
        NUMER =
           (*WDF2 * DSMDZ2) + SLOPX * *WCND2 - (*WDF * DSMDZ) - *WCND + ET[K];

/*----------------------------------------------------------------------
* CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
* --------------------------------------------------------------------*/
        RHSTT[K] = NUMER / (-DENOM2);
        AI[K] = -*WDF * DDZ / DENOM2;

/*----------------------------------------------------------------------
* RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
* RUNOFF2:  SUB-SURFACE OR BASEFLOW RUNOFF
* --------------------------------------------------------------------*/
        BI[K] = -(AI[K] + CI[K]);
        if (K == *NSOIL - 1)
            *RUNOFF2 = SLOPX * *WCND2;
        if (K != *NSOIL - 1)
        {
            *WDF = *WDF2;
            *WCND = *WCND2;
            DSMDZ = DSMDZ2;
            DDZ = DDZ2;
        }
    }
#endif

    free (SICEMAX);
    free (WCND);
    free (WCND2);
    free (WDF);
    free (WDF2);

/*----------------------------------------------------------------------
  END SUBROUTINE SRT
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
SSTEP (double *SH2OOUT, double *SH2OIN, double *CMC, double *RHSTT,
   double *RHSCT, double *DT, int *NSOIL, double *SMCMAX, double *SMCMIN,
   double *CMCMAX, double *RUNOFF3, double *ZSOIL, double *SMC, double *SICE,
   double *AI, double *BI, double *CI)
#else
void
SSTEP (double *SH2OOUT, double *SH2OIN, double *CMC, double *RHSTT,
   double *RHSCT, double *DT, int *NSOIL, double *SMCMAX, double *CMCMAX,
   double *RUNOFF3, double *ZSOIL, double *SMC, double *SICE, double *AI,
   double *BI, double *CI)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE SSTEP
* ----------------------------------------------------------------------
* CALCULATE/UPDATE SOIL MOISTURE CONTENT VALUES AND CANOPY MOISTURE
* CONTENT VALUES.
* --------------------------------------------------------------------*/
    int             K, KK11;

    double          RHSTTin[*NSOIL], CIin[*NSOIL];
    double          SH2Omid[*NSOIL];
    double          DDZ, STOT, WPLUS;

/*----------------------------------------------------------------------
* CREATE 'AMOUNT' VALUES OF VARIABLES TO BE INPUT TO THE
* TRI-DIAGONAL MATRIX ROUTINE.
* --------------------------------------------------------------------*/
    for (K = 0; K < *NSOIL; K++)
    {
        RHSTT[K] = RHSTT[K] * *DT;
        AI[K] = AI[K] * *DT;
        BI[K] = 1. + BI[K] * *DT;
        CI[K] = CI[K] * *DT;
    }

/*----------------------------------------------------------------------
* COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
* --------------------------------------------------------------------*/
    for (K = 0; K < *NSOIL; K++)
        RHSTTin[K] = RHSTT[K];
    for (K = 0; K < *NSOIL; K++)
        CIin[K] = CI[K];

/*----------------------------------------------------------------------
* CALL ROSR12 TO SOLVE THE TRI-DIAGONAL MATRIX
* --------------------------------------------------------------------*/
    ROSR12 (CI, AI, BI, CIin, RHSTTin, RHSTT, NSOIL);

/*----------------------------------------------------------------------
* SUM THE PREVIOUS SMC VALUE AND THE MATRIX SOLUTION TO GET A
* NEW VALUE.  MIN ALLOWABLE VALUE OF SMC WILL BE 0.02.
* RUNOFF3: RUNOFF WITHIN SOIL LAYERS
* --------------------------------------------------------------------*/
    WPLUS = 0.0;
    *RUNOFF3 = 0.;

    for (K = *NSOIL - 1; K >= 0; K--)
    {
        if (K != 0)
            DDZ = ZSOIL[K - 1] - ZSOIL[K];
        else
            DDZ = -ZSOIL[0];

        SH2Omid[K] = SH2OIN[K] + CI[K] + WPLUS / DDZ;
        STOT = SH2Omid[K] + SICE[K];

        if (STOT > *SMCMAX)
        {
            if (K == 0)
                DDZ = -ZSOIL[0];
            else
            {
                KK11 = K - 1;
                DDZ = -ZSOIL[K] + ZSOIL[KK11];
            }
            WPLUS = (STOT - *SMCMAX) * DDZ;
        }
        else
            WPLUS = 0.0;

        if (STOT < *SMCMAX)
            SMC[K] = STOT;
        else
            SMC[K] = *SMCMAX;

        SH2Omid[K] = SMC[K] - SICE[K];
    }

    DDZ = -ZSOIL[0];
    for (K = 0; K < *NSOIL; K++)
    {
        if (K != 0)
            DDZ = ZSOIL[K - 1] - ZSOIL[K];
        SH2OOUT[K] = SH2Omid[K] + WPLUS / DDZ;
        STOT = SH2OOUT[K] + SICE[K];
        if (STOT > *SMCMAX)
        {
            if (K == 0)
                DDZ = -ZSOIL[0];
            else
            {
                KK11 = K - 1;
                DDZ = -ZSOIL[K] + ZSOIL[KK11];
            }
            WPLUS = (STOT - *SMCMAX) * DDZ;
        }
        else
            WPLUS = 0.;

        SMC[K] = STOT < *SMCMAX ? STOT : *SMCMAX;
#ifdef _FLUX_PIHM_
        SMC[K] = SMC[K] > *SMCMIN + 0.02 ? SMC[K] : *SMCMIN + 0.02;
#else
        SMC[K] = SMC[K] > 0.02 ? SMC[K] : 0.02;
#endif
        SH2OOUT[K] = SMC[K] - SICE[K];
        SH2OOUT[K] = SH2OOUT[K] > 0 ? SH2OOUT[K] : 0;
    }

/*----------------------------------------------------------------------
* UPDATE CANOPY WATER CONTENT/INTERCEPTION (CMC).  CONVERT RHSCT TO
* AN 'AMOUNT' VALUE AND ADD TO PREVIOUS CMC VALUE TO GET NEW CMC.
* --------------------------------------------------------------------*/
    *RUNOFF3 = WPLUS;
    *CMC = *CMC + *DT * *RHSCT;
    if (*CMC < 1.e-20)
        *CMC = 0.0;
    *CMC = *CMC < *CMCMAX ? *CMC : *CMCMAX;

/*----------------------------------------------------------------------
  END SUBROUTINE SSTEP
* --------------------------------------------------------------------*/
}

void
TBND (double *TU, double *TB, double *ZSOIL, double *ZBOT, int K, int *NSOIL,
   double *TBND1)
{

/*----------------------------------------------------------------------
* SUBROUTINE TBND
* ----------------------------------------------------------------------
* CALCULATE TEMPERATURE ON THE BOUNDARY OF THE LAYER BY INTERPOLATION OF
* THE MIDDLE LAYER TEMPERATURES
* --------------------------------------------------------------------*/
    double          ZB, ZUP;

/*----------------------------------------------------------------------
* USE SURFACE TEMPERATURE ON THE TOP OF THE FIRST LAYER
* --------------------------------------------------------------------*/
    if (K == 0)
        ZUP = 0.;
    else
        ZUP = ZSOIL[K - 1];

/*----------------------------------------------------------------------
* USE DEPTH OF THE CONSTANT BOTTOM TEMPERATURE WHEN INTERPOLATE
* TEMPERATURE INTO THE LAST LAYER BOUNDARY
* --------------------------------------------------------------------*/
    if (K == *NSOIL - 1)
        ZB = 2. * *ZBOT - ZSOIL[K];
    else
        ZB = ZSOIL[K + 1];

/*----------------------------------------------------------------------
* LINEAR INTERPOLATION BETWEEN THE AVERAGE LAYER TEMPERATURES
* --------------------------------------------------------------------*/

    *TBND1 = *TU + (*TB - *TU) * (ZUP - ZSOIL[K]) / (ZUP - ZB);

/*----------------------------------------------------------------------
  END SUBROUTINE TBND
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
TDFCND (double *DF, double *SMC, double *QZ, double *SMCMAX, double *SMCMIN,
   double *SH2O)
#else
void
TDFCND (double *DF, double *SMC, double *QZ, double *SMCMAX, double *SH2O)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE TDFCND
* ----------------------------------------------------------------------
* CALCULATE THERMAL DIFFUSIVITY AND CONDUCTIVITY OF THE SOIL FOR A GIVEN
* POINT AND TIME.
* ----------------------------------------------------------------------
* PETERS-LIDARD APPROACH (PETERS-LIDARD et al., 1998)
* June 2001 CHANGES: FROZEN SOIL CONDITION.
* --------------------------------------------------------------------*/
    double          AKE, GAMMD, THKDRY, THKICE, THKO, THKQTZ, THKSAT, THKS,
       THKW, SATRATIO, XU, XUNFROZ;

/*----------------------------------------------------------------------
* WE NOW GET QUARTZ AS AN INPUT ARGUMENT (SET IN ROUTINE REDPRM):
*      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52,
*     &             0.35, 0.60, 0.40, 0.82/
* ----------------------------------------------------------------------
* IF THE SOIL HAS ANY MOISTURE CONTENT COMPUTE A PARTIAL SUM/PRODUCT
* OTHERWISE USE A CONSTANT VALUE WHICH WORKS WELL WITH MOST SOILS
* ----------------------------------------------------------------------
*  THKW ......WATER THERMAL CONDUCTIVITY
*  THKQTZ ....THERMAL CONDUCTIVITY FOR QUARTZ
*  THKO ......THERMAL CONDUCTIVITY FOR OTHER SOIL COMPONENTS
*  THKS ......THERMAL CONDUCTIVITY FOR THE SOLIDS COMBINED(QUARTZ+OTHER)
*  THKICE ....ICE THERMAL CONDUCTIVITY
*  SMCMAX ....POROSITY (= SMCMAX)
*  QZ .........QUARTZ CONTENT (SOIL TYPE DEPENDENT)
* ----------------------------------------------------------------------
* USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).

*                                  PABLO GRUNMANN, 08/17/98
* REFS.:
*      FAROUKI, O.T.,1986: THERMAL PROPERTIES OF SOILS. SERIES ON ROCK
*              AND SOIL MECHANICS, VOL. 11, TRANS TECH, 136 PP.
*      JOHANSEN, O., 1975: THERMAL CONDUCTIVITY OF SOILS. PH.D. THESIS,
*              UNIVERSITY OF TRONDHEIM,
*      PETERS-LIDARD, C. D., ET AL., 1998: THE EFFECT OF SOIL THERMAL
*              CONDUCTIVITY PARAMETERIZATION ON SURFACE ENERGY FLUXES
*              AND TEMPERATURES. JOURNAL OF THE ATMOSPHERIC SCIENCES,
*              VOL. 55, PP. 1209-1224.
* --------------------------------------------------------------------*/
    // NEEDS PARAMETERS
    // POROSITY(SOIL TYPE):
    //      POROS = SMCMAX
    // SATURATION RATIO:
    // PARAMETERS  W/(M.K)
#ifdef _FLUX_PIHM_
    SATRATIO = (*SMC - *SMCMIN) / (*SMCMAX - *SMCMIN);
#else
    SATRATIO = *SMC / *SMCMAX;
#endif

    /*
     * ICE CONDUCTIVITY: 
     */
    THKICE = 2.2;

    /*
     * WATER CONDUCTIVITY: 
     */
    THKW = 0.57;

    /*
     * THERMAL CONDUCTIVITY OF "OTHER" SOIL COMPONENTS 
     */
    //      IF (QZ .LE. 0.2) THKO = 3.0
    THKO = 2.0;

    /*
     * QUARTZ' CONDUCTIVITY 
     */
    THKQTZ = 7.7;

    /*
     * SOLIDS' CONDUCTIVITY 
     */
    THKS = pow (THKQTZ, *QZ) * pow (THKO, 1. - *QZ);

    /*
     * UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN)) 
     */
    XUNFROZ = *SH2O / *SMC;

    /*
     * UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ) 
     */
    XU = XUNFROZ * (*SMCMAX);

    /*
     * SATURATED THERMAL CONDUCTIVITY 
     */
    THKSAT =
       pow (THKS, 1. - *SMCMAX) * pow (THKICE, *SMCMAX - XU) * pow (THKW, XU);

    /*
     * DRY DENSITY IN KG/M3 
     */
    GAMMD = (1. - *SMCMAX) * 2700.;

    /*
     * DRY THERMAL CONDUCTIVITY IN W.M-1.K-1 
     */
    THKDRY = (0.135 * GAMMD + 64.7) / (2700. - 0.947 * GAMMD);

    /*
     * FROZEN 
     */
    if ((*SH2O + 0.0005) < *SMC)
        AKE = SATRATIO;

    /*
     * UNFROZEN
     * RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
     */
    else
    {

        /*
         * KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT
         * LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
         * (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).
         */

        if (SATRATIO > 0.1)
            AKE = log10 (SATRATIO) + 1.0;

        /*
         * USE K = KDRY 
         */
        else
            AKE = 0.0;
    }

    /*
     * THERMAL CONDUCTIVITY 
     */

    *DF = AKE * (THKSAT - THKDRY) + THKDRY;

/*----------------------------------------------------------------------
  END SUBROUTINE TDFCND
* --------------------------------------------------------------------*/
}

void
TMPAVG (double *TAVG, double *TUP, double *TM, double *TDN, double *ZSOIL,
   int *NSOIL, int K)
{

/*----------------------------------------------------------------------
* SUBROUTINE TMPAVG
* ----------------------------------------------------------------------
* CALCULATE SOIL LAYER AVERAGE TEMPERATURE (TAVG) IN FREEZING/THAWING
* LAYER USING UP, DOWN, AND MIDDLE LAYER TEMPERATURES (TUP, TDN, TM),
* WHERE TUP IS AT TOP BOUNDARY OF LAYER, TDN IS AT BOTTOM BOUNDARY OF
* LAYER.  TM IS LAYER PROGNOSTIC STATE TEMPERATURE.
* --------------------------------------------------------------------*/
    double          DZ, DZH, X0, XDN, XUP;
    double          T0 = 2.7315e2;

/*--------------------------------------------------------------------*/

    if (K == 0)
        DZ = -ZSOIL[0];
    else
        DZ = ZSOIL[K - 1] - ZSOIL[K];

    DZH = DZ * 0.5;
    if (*TUP < T0)
    {
        if (*TM < T0)
        {

/*----------------------------------------------------------------------
* TUP, TM, TDN < T0
* --------------------------------------------------------------------*/
            if (*TDN < T0)
                *TAVG = (*TUP + 2.0 * *TM + *TDN) / 4.0;

/*----------------------------------------------------------------------
* TUP & TM < T0,  TDN .ge. T0
* --------------------------------------------------------------------*/
            else
            {
                X0 = (T0 - *TM) * DZH / (*TDN - *TM);
                *TAVG =
                   0.5 * (*TUP * DZH + *TM * (DZH + X0) + T0 * (2. * DZH -
                      X0)) / DZ;
            }
        }
        else
        {

/*----------------------------------------------------------------------
* TUP < T0, TM .ge. T0, TDN < T0
* --------------------------------------------------------------------*/
            if (*TDN < T0)
            {
                XUP = (T0 - *TUP) * DZH / (*TM - *TUP);
                XDN = DZH - (T0 - *TM) * DZH / (*TDN - *TM);
                *TAVG =
                   0.5 * (*TUP * XUP + T0 * (2. * DZ - XUP - XDN) +
                   *TDN * XDN) / DZ;
            }

/*----------------------------------------------------------------------
* TUP < T0, TM .ge. T0, TDN .ge. T0
* --------------------------------------------------------------------*/
            else
            {
                XUP = (T0 - *TUP) * DZH / (*TM - *TUP);
                *TAVG = 0.5 * (*TUP * XUP + T0 * (2. * DZ - XUP)) / DZ;
            }
        }
    }
    else
    {
        if (*TM < T0)
        {

/*----------------------------------------------------------------------
* TUP .ge. T0, TM < T0, TDN < T0
* --------------------------------------------------------------------*/
            if (*TDN < T0)
            {
                XUP = DZH - (T0 - *TUP) * DZH / (*TM - *TUP);
                *TAVG =
                   0.5 * (T0 * (DZ - XUP) + *TM * (DZH + XUP) +
                   *TDN * DZH) / DZ;
            }

/*----------------------------------------------------------------------
* TUP .ge. T0, TM < T0, TDN .ge. T0
* --------------------------------------------------------------------*/
            else
            {
                XUP = DZH - (T0 - *TUP) * DZH / (*TM - *TUP);
                XDN = (T0 - *TM) * DZH / (*TDN - *TM);
                *TAVG =
                   0.5 * (T0 * (2. * DZ - XUP - XDN) + *TM * (XUP +
                      XDN)) / DZ;
            }
        }
        else
        {

/*----------------------------------------------------------------------
* TUP .ge. T0, TM .ge. T0, TDN < T0
* --------------------------------------------------------------------*/
            if (*TDN < T0)
            {
                XDN = DZH - (T0 - *TM) * DZH / (*TDN - *TM);
                *TAVG = (T0 * (DZ - XDN) + 0.5 * (T0 + *TDN) * XDN) / DZ;
            }

/*----------------------------------------------------------------------
* TUP .ge. T0, TM .ge. T0, TDN .ge. T0
* --------------------------------------------------------------------*/
            else
                *TAVG = (*TUP + 2.0 * *TM + *TDN) / 4.0;
        }
    }

/*----------------------------------------------------------------------
  END SUBROUTINE TMPAVG
* --------------------------------------------------------------------*/
}

void
TRANSP (double *ET, int *NSOIL, double *ETP1, double *SMC, double *CMC,
   double *ZSOIL, double *SHDFAC, double *SMCWLT, double *CMCMAX, double *PC,
   double *CFACTR, double *SMCREF, double *SFCTMP, double *Q2, int *NROOT,
   double *RTDIS)
{

/*----------------------------------------------------------------------
* SUBROUTINE TRANSP
* ----------------------------------------------------------------------
* CALCULATE TRANSPIRATION FOR THE VEG CLASS.
* --------------------------------------------------------------------*/
    int             I, K;
    double          DENOM;
    double          ETP1A;
    //.....REAL PART(NSOIL)
    double          GX[*NROOT];
    double          RTX, SGX;

/*----------------------------------------------------------------------
* INITIALIZE PLANT TRANSP TO ZERO FOR ALL SOIL LAYERS.
* --------------------------------------------------------------------*/
    for (K = 0; K < *NSOIL; K++)
        ET[K] = 0.;

/*----------------------------------------------------------------------
* CALCULATE AN 'ADJUSTED' POTENTIAL TRANSPIRATION
* IF STATEMENT BELOW TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
* NOTE: GX AND OTHER TERMS BELOW REDISTRIBUTE TRANSPIRATION BY LAYER,
* ET(K), AS A FUNCTION OF SOIL MOISTURE AVAILABILITY, WHILE PRESERVING
* TOTAL ETP1A.
* --------------------------------------------------------------------*/
    if (*CMC != 0.0)
        ETP1A = *SHDFAC * *PC * *ETP1 * (1.0 - pow (*CMC / *CMCMAX, *CFACTR));
    else
        ETP1A = *SHDFAC * *PC * *ETP1;
    SGX = 0.0;
    for (I = 0; I < *NROOT; I++)
    {
        GX[I] = (SMC[I] - *SMCWLT) / (*SMCREF - *SMCWLT);
        if (GX[I] < 0)
            GX[I] = 0;
        if (GX[I] > 1)
            GX[I] = 1.0;
        SGX = SGX + GX[I];
    }

    SGX = SGX / (double)*NROOT;
    DENOM = 0.;
    for (I = 0; I < *NROOT; I++)
    {
        RTX = RTDIS[I] + GX[I] - SGX;
        GX[I] = GX[I] * (RTX > 0. ? RTX : 0.);
        DENOM = DENOM + GX[I];
    }

    if (DENOM <= 0.0)
        DENOM = 1.;
    for (I = 0; I < *NROOT; I++)
    {
        ET[I] = ETP1A * GX[I] / DENOM;

/*----------------------------------------------------------------------
* ABOVE CODE ASSUMES A VERTICALLY UNIFORM ROOT DISTRIBUTION
* CODE BELOW TESTS A VARIABLE ROOT DISTRIBUTION
* ----------------------------------------------------------------------
*		ET[0] = (ZSOIL[0] / ZSOIL(*NROOT) ) * GX * ETP1A;
*		ET[0] = (ZSOIL[0] / ZSOIL(*NROOT) ) * ETP1A;
* ----------------------------------------------------------------------
* USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
* ----------------------------------------------------------------------
*		ET[0] = RTDIS[0] * ETP1A;
*		ET[0] = ETP1A * PART[0];
* ----------------------------------------------------------------------
* LOOP DOWN THRU THE SOIL LAYERS REPEATING THE OPERATION ABOVE,
* BUT USING THE THICKNESS OF THE SOIL LAYER (RATHER THAN THE
* ABSOLUTE DEPTH OF EACH LAYER) IN THE FINAL CALCULATION.
* ----------------------------------------------------------------------
*		for (K = 1; K < *NROOT; K++)
		{
*			GX = (SMC[K] - *SMCWLT ) / (*SMCREF - *SMCWLT);
*			GX = GX < 1.0 ? GX : 1.0;
*			GX = GX > 0.0 ? GX : 0.0;
* TEST CANOPY RESISTANCE
*			GX = 1.0;
*			ET[K] = ((ZSOIL[K] - ZSOIL[K-1]) / ZSOIL[*NROOT - 1]) * GX * ETP1A;
*			ET[K] = ((ZSOIL[K] - ZSOIL[K-1]) / ZSOIL[*NROOT - 1]) * ETP1A;
* ----------------------------------------------------------------------
* USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
* ----------------------------------------------------------------------
*			ET[K] = RTDIS[K] * ETP1A;
*			ET[K] = ETP1A * PART[K];
*		}
*/
    }

/*----------------------------------------------------------------------
  END SUBROUTINE TRANSP
* --------------------------------------------------------------------*/
}

#ifdef _FLUX_PIHM_
void
WDFCND (double *WDF, double *WCND, double *SMC, double *SMCMAX, double *SMCMIN, double *VGALPHA, double *VGBETA, double *DKSAT, double *MACKSAT, double *AREAF, int *MAC_STATUS, double *SICEMAX, double *DSMDZ, int *MACPORE)
#else
void
WDFCND (double *WDF, double *WCND, double *SMC, double *SMCMAX, double *BEXP, double *DKSAT, double *DWSAT, double *SICEMAX)
#endif
{

/*----------------------------------------------------------------------
* SUBROUTINE WDFCND
* ----------------------------------------------------------------------
* CALCULATE SOIL WATER DIFFUSIVITY AND SOIL HYDRAULIC CONDUCTIVITY.
* --------------------------------------------------------------------*/
    double          EXPON, FACTR1, FACTR2, VKWGT;
#ifdef _FLUX_PIHM_
    double          SATKFUNC, DPSIDSM;
#endif

/*----------------------------------------------------------------------
*     CALC THE RATIO OF THE ACTUAL TO THE MAX PSBL SOIL H2O CONTENT
* --------------------------------------------------------------------*/

#ifdef _FLUX_PIHM_
    FACTR1 = 0.05 / (*SMCMAX - *SMCMIN);
    FACTR2 = (*SMC - *SMCMIN) / (*SMCMAX - *SMCMIN);

/*----------------------------------------------------------------------
* FACTR2 should avoid to be 0 or 1
* --------------------------------------------------------------------*/
    if (FACTR2 > 1. - .0005)
        FACTR2 = 1. - .0005;
    if (FACTR2 < 0. + .0005)
        FACTR2 = .0005;
    FACTR1 = FACTR1 < FACTR2 ? FACTR1 : FACTR2;
    EXPON = 1.0 - 1. / *VGBETA;

    SATKFUNC =
       pow (FACTR2, 0.5) * pow (1. - pow (1. - pow (FACTR2, 1. / EXPON),
          EXPON), 2.);
    DPSIDSM =
       (1. - EXPON) / *VGALPHA / EXPON / (*SMCMAX -
       *SMCMIN) * pow (pow (FACTR2, -1. / EXPON) - 1.,
       0. - EXPON) * pow (FACTR2, -(1. / EXPON + 1.));

    if (*MACPORE == 1)
        *WCND = EFFKV (SATKFUNC, FACTR2, *MAC_STATUS, *MACKSAT, *DKSAT, *AREAF);
    else
        *WCND = *DKSAT * SATKFUNC;

    *WDF = *WCND * DPSIDSM;

    //  *WDF = (1. - EXPON) * (*DKSAT * (1. - *AREAF) + *MACKSAT * *AREAF) / *VGALPHA / EXPON / (*SMCMAX - *SMCMIN) * pow(FACTR2, 0.5 - 1. / EXPON) * (pow(1. - pow(FACTR2, 1. / EXPON) ,-EXPON) + pow(1. - pow(FACTR2, 1. / EXPON) ,EXPON) - 2.);
    //  *WCND = sqrt(FACTR2) * pow(1. - pow(1. - pow(FACTR2, 1. / EXPON), EXPON), 2.) * (*DKSAT * (1. - *AREAF) + *MACKSAT * *AREAF);

    if (*SICEMAX > 0.0)
    {
        VKWGT = 1. / (1. + pow (500. * *SICEMAX, 3.));
        SATKFUNC =
           pow (FACTR1, 0.5) * pow (1. - pow (1. - pow (FACTR1, 1. / EXPON),
              EXPON), 2.);
        DPSIDSM =
           (1. - EXPON) / *VGALPHA / EXPON / (*SMCMAX -
           *SMCMIN) * pow (pow (FACTR1, -1. / EXPON) - 1.,
           0. - EXPON) * pow (FACTR1, -(1. / EXPON + 1.));
        if (*MACPORE == 1)
            *WDF =
               VKWGT * *WDF + (1. - VKWGT) * DPSIDSM * EFFKV (SATKFUNC, FACTR1, *MAC_STATUS, *MACKSAT, *DKSAT, *AREAF);
        else
            *WDF = VKWGT * *WDF + (1. - VKWGT) * DPSIDSM * SATKFUNC * *DKSAT;
    }


#else

/*----------------------------------------------------------------------
* PREP AN EXPNTL COEF AND CALC THE SOIL WATER DIFFUSIVITY
* --------------------------------------------------------------------*/
    FACTR1 = 0.05 / *SMCMAX;
    FACTR2 = *SMC / *SMCMAX;
    FACTR1 = FACTR1 < FACTR2 ? FACTR1 : FACTR2;
    EXPON = *BEXP + 2.0;

/*----------------------------------------------------------------------
* FROZEN SOIL HYDRAULIC DIFFUSIVITY.  VERY SENSITIVE TO THE VERTICAL
* GRADIENT OF UNFROZEN WATER. THE LATTER GRADIENT CAN BECOME VERY
* EXTREME IN FREEZING/THAWING SITUATIONS, AND GIVEN THE RELATIVELY
* FEW AND THICK SOIL LAYERS, THIS GRADIENT SUFFERES SERIOUS
* TRUNCTION ERRORS YIELDING ERRONEOUSLY HIGH VERTICAL TRANSPORTS OF
* UNFROZEN WATER IN BOTH DIRECTIONS FROM HUGE HYDRAULIC DIFFUSIVITY.
* THEREFORE, WE FOUND WE HAD TO ARBITRARILY CONSTRAIN WDF
* --
* VERSION D_10CM: ........  FACTR1 = 0.2/SMCMAX
* WEIGHTED APPROACH...................... PABLO GRUNMANN, 28_SEP_1999.
* --------------------------------------------------------------------*/
    *WDF = *DWSAT * pow (FACTR2, EXPON);
    if (*SICEMAX > 0.0)
    {
        VKWGT = 1. / (1. + pow (500. * *SICEMAX, 3.));
        *WDF = VKWGT * *WDF + (1. - VKWGT) * *DWSAT * pow (FACTR1, EXPON);

/*----------------------------------------------------------------------
* RESET THE EXPNTL COEF AND CALC THE HYDRAULIC CONDUCTIVITY
* --------------------------------------------------------------------*/
    }
    EXPON = (2.0 * *BEXP) + 3.0;
    *WCND = *DKSAT * pow (FACTR2, EXPON);
#endif

/*----------------------------------------------------------------------
  END SUBROUTINE WDFCND
* --------------------------------------------------------------------*/
}


void SFCDIF_off (double *ZLM, double *ZLM_WIND, double *Z0, double *THZ0, double *THLM, double *SFCSPD, double *CZIL, double *AKMS, double *AKHS, int *VEGTYP, int *ISURBAN, int *IZ0TLND)
{

/*----------------------------------------------------------------------
* SUBROUTINE SFCDIF (renamed SFCDIF_off to avoid clash with Eta PBL)
* ----------------------------------------------------------------------
* CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
* SEE CHEN ET AL (1997, BLM)
* --------------------------------------------------------------------*/

    double          ZILFC, ZU, ZT, RDZ, CXCH;
    double          DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT;
    double          RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4;
    //!CC   ......REAL ZTFC

    double          XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN,
       RLMA;

    int             ILECH, ITR;

    double          WWST = 1.2;
    double          WWST2;
    double          G = 9.8, VKRM = 0.40, EXCM = 0.001, BETA = 1. / 270.;
    double          BTG, ELFC;
    double          WOLD = .15;
    double          WNEW;
    int             ITRMX = 5;

    double          EPSU2 = 1.e-4;
    double          EPSUST = 0.07;
    //  double EPSA = 1.e-8;
    double          ZTMIN = -5.;
    double          ZTMAX = 1.;
    double          HPBL = 1000.0;
    double          SQVISC = 258.2;

    WWST2 = WWST * WWST;
    BTG = BETA * G;
    ELFC = VKRM * BTG;
    WNEW = 1. - WOLD;

/*----------------------------------------------------------------------
*     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
*     C......ZTFC=0.1
*     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
* --------------------------------------------------------------------*/
    ILECH = 0;

/*--------------------------------------------------------------------*/
    if ((*IZ0TLND == 0) || (*VEGTYP == *ISURBAN))
    {
        /* Just use the original CZIL value. */
        ZILFC = -*CZIL * VKRM * SQVISC;
    }
    else
    {
        /* Modify CZIL according to Chen & Zhang, 2009
         * CZIL = 10 ** -0.40 H, ( where H = 10*Zo ) */
        *CZIL = pow (10.0, -0.4 * (*Z0 / 0.07));
        ZILFC = -*CZIL * VKRM * SQVISC;
    }
    //     C.......ZT=Z0*ZTFC
    ZU = *Z0;
    RDZ = 1. / *ZLM_WIND;
    CXCH = EXCM * RDZ;
    DTHV = *THLM - *THZ0;

/*----------------------------------------------------------------------
* BELJARS CORRECTION OF USTAR
* --------------------------------------------------------------------*/
    DU2 = (*SFCSPD * *SFCSPD) > EPSU2 ? (*SFCSPD * *SFCSPD) : EPSU2;

    /*
     * cc   If statements to avoid TANGENT LINEAR problems near zero 
     */
    BTGH = BTG * HPBL;
    if (BTGH * *AKHS * DTHV != 0.0)
        WSTAR2 = WWST2 * pow (fabs (BTGH * *AKHS * DTHV), 2. / 3.);
    else
        WSTAR2 = 0.0;

/*----------------------------------------------------------------------
* ZILITINKEVITCH APPROACH FOR ZT
* --------------------------------------------------------------------*/
    USTAR = sqrt (*AKMS * sqrt (DU2 + WSTAR2));
    USTAR = USTAR > EPSUST ? USTAR : EPSUST;

/*--------------------------------------------------------------------*/
    ZT = exp (ZILFC * sqrt (USTAR * *Z0)) * *Z0;
    ZSLU = *ZLM_WIND + ZU;
    //     PRINT*,'ZSLT=',ZSLT
    //     PRINT*,'ZLM=',ZLM
    //     PRINT*,'ZT=',ZT

    ZSLT = *ZLM + ZT;
//    ZSLT = *ZLM_WIND + ZT;
    RLOGU = log (ZSLU / ZU);

    RLOGT = log (ZSLT / ZT);
    //     PRINT*,'RLMO=',RLMO
    //     PRINT*,'ELFC=',ELFC
    //     PRINT*,'AKHS=',AKHS
    //     PRINT*,'DTHV=',DTHV
    //     PRINT*,'USTAR=',USTAR

    RLMO = ELFC * *AKHS * DTHV / pow (USTAR, 3);

/*----------------------------------------------------------------------
* 1./MONIN-OBUKKHOV LENGTH-SCALE
* --------------------------------------------------------------------*/
    for (ITR = 0; ITR < ITRMX; ITR++)
    {
        ZETALT = ZSLT * RLMO;
        ZETALT = ZETALT > ZTMIN ? ZETALT : ZTMIN;
        RLMO = ZETALT / ZSLT;
        ZETALU = ZSLU * RLMO;
        ZETAU = ZU * RLMO;

        ZETAT = ZT * RLMO;
        if (ILECH == 0)
        {
            if (RLMO < 0.)
            {
                XLU4 = 1. - 16. * ZETALU;
                XLT4 = 1. - 16. * ZETALT;
                XU4 = 1. - 16. * ZETAU;
                XT4 = 1. - 16. * ZETAT;
                XLU = sqrt (sqrt (XLU4));
                XLT = sqrt (sqrt (XLT4));
                XU = sqrt (sqrt (XU4));

                XT = sqrt (sqrt (XT4));
                //     PRINT*,'-----------1------------'
                //     PRINT*,'PSMZ=',PSMZ
                //     PRINT*,'PSPMU(ZETAU)=',PSPMU(ZETAU)
                //     PRINT*,'XU=',XU
                //     PRINT*,'------------------------'
                PSMZ = PSPMU (XU);
                SIMM = PSPMU (XLU) - PSMZ + RLOGU;
                PSHZ = PSPHU (XT);
                SIMH = PSPHU (XLT) - PSHZ + RLOGT;
            }
            else
            {
                ZETALU = ZETALU < ZTMAX ? ZETALU : ZTMAX;
                ZETALT = ZETALT < ZTMAX ? ZETALT : ZTMAX;
                //     PRINT*,'-----------2------------'
                //     PRINT*,'PSMZ=',PSMZ
                //     PRINT*,'PSPMS(ZETAU)=',PSPMS(ZETAU)
                //     PRINT*,'ZETAU=',ZETAU
                //     PRINT*,'------------------------'
                PSMZ = PSPMS (ZETAU);
                SIMM = PSPMS (ZETALU) - PSMZ + RLOGU;
                PSHZ = PSPHS (ZETAT);
                SIMH = PSPHS (ZETALT) - PSHZ + RLOGT;
            }
        }

/*----------------------------------------------------------------------
* LECH'S FUNCTIONS
* --------------------------------------------------------------------*/
        else
        {
            if (RLMO < 0.)
            {
                //     PRINT*,'-----------3------------'
                //     PRINT*,'PSMZ=',PSMZ
                //     PRINT*,'PSLMU(ZETAU)=',PSLMU(ZETAU)
                //     PRINT*,'ZETAU=',ZETAU
                //     PRINT*,'------------------------'
                PSMZ = PSLMU (ZETAU);
                SIMM = PSLMU (ZETALU) - PSMZ + RLOGU;
                PSHZ = PSLHU (ZETAT);
                SIMH = PSLHU (ZETALT) - PSHZ + RLOGT;
            }
            else
            {
                ZETALU = ZETALU < ZTMAX ? ZETALU : ZTMAX;
                ZETALT = ZETALT < ZTMAX ? ZETALT : ZTMAX;
                //     PRINT*,'-----------4------------'
                //     PRINT*,'PSMZ=',PSMZ
                //     PRINT*,'PSLMS(ZETAU)=',PSLMS(ZETAU)
                //     PRINT*,'ZETAU=',ZETAU
                //     PRINT*,'------------------------'
                PSMZ = PSLMS (ZETAU);
                SIMM = PSLMS (ZETALU) - PSMZ + RLOGU;
                PSHZ = PSLHS (ZETAT);
                SIMH = PSLHS (ZETALT) - PSHZ + RLOGT;
            }
        }

/*----------------------------------------------------------------------
* BELJAARS CORRECTION FOR USTAR
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* ZILITINKEVITCH FIX FOR ZT
* --------------------------------------------------------------------*/
        USTAR = sqrt (*AKMS * sqrt (DU2 + WSTAR2));
        USTAR = USTAR > EPSUST ? USTAR : EPSUST;

        ZT = exp (ZILFC * sqrt (USTAR * *Z0)) * *Z0;
        ZSLT = *ZLM + ZT;

/*--------------------------------------------------------------------*/
        RLOGT = log (ZSLT / ZT);
        USTARK = USTAR * VKRM;
        *AKMS = USTARK / SIMM > CXCH ? USTARK / SIMM : CXCH;

/*----------------------------------------------------------------------
* IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
*---------------------------------------------------------------------*/
        *AKHS = USTARK / SIMH > CXCH ? USTARK / SIMH : CXCH;
        if (BTGH * *AKHS * DTHV != 0.0)
            WSTAR2 = WWST2 * pow (fabs (BTGH * *AKHS * DTHV), 2. / 3.);
        else
            WSTAR2 = 0.0;

/*--------------------------------------------------------------------*/
        RLMN = ELFC * *AKHS * DTHV / pow (USTAR, 3.0);

/*----------------------------------------------------------------------
*     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
*---------------------------------------------------------------------*/
        RLMA = RLMO * WOLD + RLMN * WNEW;

/*--------------------------------------------------------------------*/
        RLMO = RLMA;
        //     PRINT*,'----------------------------'
        //     PRINT*,'SFCDIF OUTPUT !  ! ! ! ! ! ! ! !  !   !    !'

        //     PRINT*,'ZLM=',ZLM
        //     PRINT*,'Z0=',Z0
        //     PRINT*,'THZ0=',THZ0
        //     PRINT*,'THLM=',THLM
        //     PRINT*,'SFCSPD=',SFCSPD
        //     PRINT*,'CZIL=',CZIL
        //     PRINT*,'AKMS=',AKMS
        //     PRINT*,'AKHS=',AKHS
        //     PRINT*,'----------------------------'
    }

/*----------------------------------------------------------------------
  END SUBROUTINE SFCDIF_off
* --------------------------------------------------------------------*/
}

/*----------------------------------------------------------------------
* NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
* ----------------------------------------------------------------------
* LECH'S SURFACE FUNCTIONS
* --------------------------------------------------------------------*/
double PSLMU (double ZZ)
{
    double          x;
    x = -0.96 * log (1.0 - 4.5 * ZZ);
    return x;
}

double PSLMS (double ZZ)
{
    double          RIC = 0.183, RRIC;
    double          x;
    RRIC = 1.0 / RIC;
    x = ZZ * RRIC - 2.076 * (1. - 1. / (ZZ + 1.));
    return x;
}

double PSLHU (double ZZ)
{
    double          x;
    x = -0.96 * log (1.0 - 4.5 * ZZ);
    return x;
}

double PSLHS (double ZZ)
{
    double          RIC = 0.183;
    double          FHNEU = 0.8, RFC = 0.191;
    double          RFAC;
    double          x;

    RFAC = RIC / (FHNEU * RFC * RFC);
    //  x = ZZ * RFAC -2.076* (1. -1./ (ZZ +1.));
    x = ZZ * RFAC - 2.076 * (1. - exp (-1.2 * ZZ));
    printf("now: %lf, before: %lf\n", x, ZZ * RFAC -2.076* (1. -1./ (ZZ +1.)));
    return x;
}

/*----------------------------------------------------------------------
* PAULSON'S SURFACE FUNCTIONS
* --------------------------------------------------------------------*/

double PSPMU (double XX)
{
    double          PIHF = 3.14159265 / 2.0;
    double          x;
    x = -2. * log ((XX + 1.) * 0.5) - log ((XX * XX + 1.) * 0.5) +
       2. * atan (XX) - PIHF;
    return x;
}

double PSPMS (double YY)
{
    double          x;
    x = 5. * YY;
    return x;
}

double PSPHU (double XX)
{
    double          x;
    x = -2. * log ((XX * XX + 1.) * 0.5);
    return x;
}

/*----------------------------------------------------------------------
* THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
* OVER SOLID SURFACE (LAND, SEA-ICE).
* --------------------------------------------------------------------*/
double PSPHS (double YY)
{
    double          x;
    x = 5. * YY;
    return x;
}

double EFFKV (double KSATFUNC, double ELEMSATN, int STATUS, double MACKV, double KV, double AREAF)
{
    //return (KV * KSATFUNC);
    if (STATUS == SAT_CTRL)
        return (MACKV * AREAF + KV * (1. - AREAF) * KSATFUNC);
    else
    {
        if (STATUS == MAT_CTRL)
            return KV * KSATFUNC;
        else
        {
            if (STATUS == APP_CTRL)
                return (MACKV * AREAF * KSATFUNC + KV * (1. - AREAF) * KSATFUNC);
            else
                return (MACKV * AREAF + KV * (1. - AREAF) * KSATFUNC);
        }
    }
}
