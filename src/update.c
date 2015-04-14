/*****************************************************************************
 * File		: update.c
 * Function	: Update variables during simulation
 ****************************************************************************/

#include "pihm.h"

void update (realtype t, void *DS)
{
    int             k;
    Model_Data      MD;

    MD = (Model_Data) DS;

    for (k = 0; k < MD->NumTS; k++)
    {
        while (MD->TSD_meteo[k].iCounter < MD->TSD_meteo[k].length && t > MD->TSD_meteo[k].TS[MD->TSD_meteo[k].iCounter + 1][0])
            MD->TSD_meteo[k].iCounter++;
    }
}

void summary (Model_Data DS, N_Vector CV_Y, realtype t, realtype stepsize)
{
    realtype       *Y, Yc[3 * DS->NumEle + 2 * DS->NumRiv];
    realtype        WTD0, WTD1, elemSatn0, elemSatn1, RealUnsat0, RealUnsat1, RealGW0, RealGW1, Recharge, Runoff;
    realtype        h, AquiferDepth;
    int             i, j;

    Y = NV_DATA_S (CV_Y);

    for (i = 0; i < DS->NumEle; i++)
    {
        Yc[i] = 0.5 * (Y[i] + DS->EleSurf[i]);
        Yc[i + DS->NumEle] = 0.5 * (Y[i + DS->NumEle] + DS->EleUnsat[i]);
        Yc[i + 2 * DS->NumEle] = 0.5 * (Y[i + 2 * DS->NumEle] + DS->EleGW[i]);
    }

    for (i = 0; i < DS->NumRiv; i++)
    {
        Yc[i + 3 * DS->NumEle] = 0.5 * (Y[i + 3 * DS->NumEle] + DS->RivStg[i]);
        Yc[i + 3 * DS->NumEle + DS->NumRiv] = 0.5 * (Y[i + 3 * DS->NumEle + DS->NumRiv] + DS->EleGW[i + DS->NumEle]);
    }

    for (i = 0; i < 3 * DS->NumEle + 2 * DS->NumRiv; i++)
        DS->DummyY[i] = Y[i] >= 0.0 ? Y[i] : 0.0;

    for (i = 0; i < DS->NumEle; i++)
    {
        h = DS->DummyY[i + 2 * DS->NumEle];
        /* Calculate infiltration based on mass conservation */
        AquiferDepth = DS->Ele[i].zmax - DS->Ele[i].zmin;
        WTD0 = AquiferDepth - (DS->EleGW[i] > 0.0 ? DS->EleGW[i] : 0.0);
        WTD0 = WTD0 < 0.0 ? 0.0 : WTD0;
        elemSatn0 = (WTD0 <= 0.0) ? 1.0 : (DS->EleUnsat[i] < 0.0 ? 0.0 : DS->EleUnsat[i] / WTD0);
        elemSatn0 = elemSatn0 > 1.0 ? 1.0 : (elemSatn0 < 0.0 ? 0.0 : elemSatn0);
        RealUnsat0 = elemSatn0 * WTD0;
        RealGW0 = DS->EleGW[i] > AquiferDepth ? AquiferDepth : (DS->EleGW[i] < 0.0 ? 0.0 : DS->EleGW[i]);
        WTD1 = AquiferDepth - DS->DummyY[i + 2 * DS->NumEle];
        WTD1 = WTD1 < 0.0 ? 0.0 : WTD1;
        elemSatn1 = (WTD1 <= 0.0) ? 1.0 : (DS->DummyY[i + DS->NumEle] < 0.0 ? 0.0 : DS->DummyY[i + DS->NumEle] / WTD1);
        elemSatn1 = elemSatn1 > 1.0 ? 1.0 : (elemSatn1 < 0.0 ? 0.0 : elemSatn1);
        RealUnsat1 = elemSatn1 * WTD1;
        RealGW1 = DS->DummyY[i + 2 * DS->NumEle] > AquiferDepth ? AquiferDepth : DS->DummyY[i + 2 * DS->NumEle];

        /* Subsurface runoff rate */
        Runoff = 0.0;
        for (j = 0; j < 3; j++)
            Runoff = Runoff + DS->FluxSub[i][j] / DS->Ele[i].area;
#ifdef _FLUX_PIHM_
        Recharge = (RealGW1 - RealGW0) * DS->Ele[i].Porosity / stepsize + Runoff + DS->EleETsat[i] * DS->EleET[i][1];
        DS->EleViR[i] = (RealUnsat1 - RealUnsat0) * DS->Ele[i].Porosity / stepsize + Recharge + (1. - DS->EleETsat[i]) * DS->EleET[i][1] + DS->EleET[i][2];
#else
        Recharge = (RealGW1 - RealGW0) * DS->Ele[i].Porosity / stepsize + Runoff + ((DS->EleGW[i] > AquiferDepth - DS->Ele[i].RzD) ? DS->EleET[i][1] : 0.0);
        DS->EleViR[i] = (RealUnsat1 - RealUnsat0) * DS->Ele[i].Porosity / stepsize + Recharge + (DS->EleSurf[i] < EPS / 100.0 ? DS->EleET[i][2] : 0.0) + ((DS->EleGW[i] <= AquiferDepth - DS->Ele[i].RzD) ? DS->EleET[i][1] : 0.0);
#endif
        DS->EleViR[i] = DS->EleViR[i] > 0.0 ? DS->EleViR[i] : 0.0;
    }
    for (i = 0; i < DS->NumEle; i++)
    {
        DS->EleSurf[i] = Y[i];
        DS->EleUnsat[i] = Y[i + DS->NumEle];
        DS->EleGW[i] = Y[i + 2 * DS->NumEle];
    }
    for (i = 0; i < DS->NumRiv; i++)
    {
        DS->RivStg[i] = Y[i + 3 * DS->NumEle];
        DS->EleGW[i + DS->NumEle] = Y[i + 3 * DS->NumEle + DS->NumRiv];
    }
}
