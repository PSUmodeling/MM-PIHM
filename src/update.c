/*****************************************************************************
 * File		: update.c
 * Function	: Update variables during simulation
 * Version	: June, 2014
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pihm.h"

void update (realtype t, void *DS)
{
    int             k;
    //int             NumForc;

    Model_Data      MD;

    MD = (Model_Data) DS;

    //printf ("iCounter before :%d", MD->TSD_meteo[0].iCounter);
    for (k = 0; k < MD->NumTS; k++)
    {
        while (MD->TSD_meteo[k].iCounter < MD->TSD_meteo[k].length && t > MD->TSD_meteo[k].TS[MD->TSD_meteo[k].iCounter + 1][0])
            MD->TSD_meteo[k].iCounter++;
    }
    //printf ("iCounter affter :%d", MD->TSD_meteo[0].iCounter);
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

    //f_update (t, Yc, DS);
    for (i = 0; i < 3 * DS->NumEle + 2 * DS->NumRiv; i++)
        DS->DummyY[i] = Y[i] >= 0. ? Y[i] : 0.;

    for (i = 0; i < DS->NumEle; i++)
    {
        h = DS->DummyY[i + 2 * DS->NumEle];
        /* Calculate infiltration based on mass conservation */
        AquiferDepth = DS->Ele[i].zmax - DS->Ele[i].zmin;
        WTD0 = AquiferDepth - (DS->EleGW[i] > 0 ? DS->EleGW[i] : 0);
        WTD0 = WTD0 < 0 ? 0 : WTD0;
        elemSatn0 = (WTD0 <= 0) ? 1.0 : (DS->EleUnsat[i] < 0 ? 0 : DS->EleUnsat[i] / WTD0);
        elemSatn0 = elemSatn0 > 1.0 ? 1.0 : (elemSatn0 < 0 ? 0 : elemSatn0);
        RealUnsat0 = elemSatn0 * WTD0;
        RealGW0 = DS->EleGW[i] > AquiferDepth ? AquiferDepth : (DS->EleGW[i] < 0 ? 0 : DS->EleGW[i]);
        WTD1 = AquiferDepth - DS->DummyY[i + 2 * DS->NumEle];
        WTD1 = WTD1 < 0 ? 0 : WTD1;
        elemSatn1 = (WTD1 <= 0) ? 1.0 : (DS->DummyY[i + DS->NumEle] < 0 ? 0 : DS->DummyY[i + DS->NumEle] / WTD1);
        elemSatn1 = elemSatn1 > 1.0 ? 1.0 : (elemSatn1 < 0 ? 0 : elemSatn1);
        RealUnsat1 = elemSatn1 * WTD1;
        RealGW1 = DS->DummyY[i + 2 * DS->NumEle] > AquiferDepth ? AquiferDepth : DS->DummyY[i + 2 * DS->NumEle];

        /* Subsurface runoff rate */
        Runoff = 0;
        for (j = 0; j < 3; j++)
            Runoff = Runoff + DS->FluxSub[i][j] / DS->Ele[i].area;
#ifdef _FLUX_PIHM_
        Recharge = (RealGW1 - RealGW0) * DS->Ele[i].Porosity / stepsize + Runoff + DS->EleETsat[i] * DS->EleET[i][1];
        DS->EleViR[i] = (RealUnsat1 - RealUnsat0) * DS->Ele[i].Porosity / stepsize + Recharge + (1. - DS->EleETsat[i]) * DS->EleET[i][1] + DS->EleET[i][2];
#else
        Recharge = (RealGW1 - RealGW0) * DS->Ele[i].Porosity / stepsize + Runoff + ((DS->EleGW[i] > AquiferDepth - DS->Ele[i].RzD) ? DS->EleET[i][1] : 0);
        DS->EleViR[i] = (RealUnsat1 - RealUnsat0) * DS->Ele[i].Porosity / stepsize + Recharge + (DS->EleSurf[i] < EPS / 100 ? DS->EleET[i][2] : 0) + ((DS->EleGW[i] <= AquiferDepth - DS->Ele[i].RzD) ? DS->EleET[i][1] : 0);
#endif
        DS->EleViR[i] = DS->EleViR[i] > 0 ? DS->EleViR[i] : 0;
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

void f_update (realtype t, realtype * Y, void *DS)
{
    int             i, j, inabr;
    //realtype        RivNetPrep;
    realtype        Avg_Y_Surf, Dif_Y_Surf, Grad_Y_Surf, Avg_Sf, Distance;
    //realtype        Cwr, 
    realtype	    TotalY_Riv, TotalY_Riv_down, CrossA, CrossAdown, AvgCrossA, Perem, Perem_down, Avg_Rough, Avg_Perem, Avg_Y_Riv, Dif_Y_Riv, Grad_Y_Riv, Wid, Wid_down, Avg_Wid;
    realtype        Avg_Y_Sub, Dif_Y_Sub, Avg_Ksat, Grad_Y_Sub, nabrAqDepth, AquiferDepth, Deficit, elemSatn, satKfunc, effK, effKnabr, TotalY_Ele, TotalY_Ele_down;
    //realtype        temp;
    realtype        dt;
    Model_Data      MD;
    MD = (Model_Data) DS;

    dt = MD->dt;
    /*
     * Initialization of temporary state variables 
     */
    for (i = 0; i < 3 * MD->NumEle + 2 * MD->NumRiv; i++)
    {
        MD->DummyY[i] = (Y[i] >= 0) ? Y[i] : 0;
        if (i < MD->NumRiv)
        {
            MD->FluxRiv[i][0] = 0;
            MD->FluxRiv[i][10] = 0;
        }
        if ((MD->SurfMode == 2) && (i < MD->NumEle))
        {
            for (j = 0; j < 3; j++)
            {
                MD->Ele[i].surfH[j] = (MD->Ele[i].nabr[j] > 0) ? ((MD->Ele[i].BC[j] > -4) ? (MD->Ele[MD->Ele[i].nabr[j] - 1].zmax + MD->DummyY[MD->Ele[i].nabr[j] - 1]) : ((MD->DummyY[-(MD->Ele[i].BC[j] / 4) - 1 + 3 * MD->NumEle] > MD->Riv[-(MD->Ele[i].BC[j] / 4) - 1].depth) ? MD->Riv[-(MD->Ele[i].BC[j] / 4) - 1].zmin + MD->DummyY[-(MD->Ele[i].BC[j] / 4) - 1 + 3 * MD->NumEle] : MD->Riv[-(MD->Ele[i].BC[j] / 4) - 1].zmax)) : ((MD->Ele[i].BC[j] != 1) ? (MD->Ele[i].zmax + MD->DummyY[i]) : Interpolation (&MD-> TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t));
            }
            MD->Ele[i].dhBYdx = -1 * (MD->Ele[i].surfY[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfX[2] * (MD->Ele[i].surfY[1] - MD->Ele[i].surfY[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfY[0] - MD->Ele[i].surfY[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfY[2] - MD->Ele[i].surfY[1]));
            MD->Ele[i].dhBYdy = -1 * (MD->Ele[i].surfX[2] * (MD->Ele[i].surfH[1] - MD->Ele[i].surfH[0]) + MD->Ele[i].surfX[1] * (MD->Ele[i].surfH[0] - MD->Ele[i].surfH[2]) + MD->Ele[i].surfX[0] * (MD->Ele[i].surfH[2] - MD->Ele[i].surfH[1])) / (MD->Ele[i].surfY[2] * (MD->Ele[i].surfX[1] - MD->Ele[i].surfX[0]) + MD->Ele[i].surfY[1] * (MD->Ele[i].surfX[0] - MD->Ele[i].surfX[2]) + MD->Ele[i].surfY[0] * (MD->Ele[i].surfX[2] - MD->Ele[i].surfX[1]));
        }
    }
    /*
     * Lateral Flux Calculation between Triangular elements Follows  
     */
    for (i = 0; i < MD->NumEle; i++)
    {
        AquiferDepth = (MD->Ele[i].zmax - MD->Ele[i].zmin);
        if (AquiferDepth < MD->Ele[i].macD)
            MD->Ele[i].macD = AquiferDepth;
        for (j = 0; j < 3; j++)
        {
            if (MD->Ele[i].nabr[j] > 0)
            {
                /*
                 * Subsurface Lateral Flux Calculation between Triangular elements Follows 
                 */
                Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - (MD->DummyY[MD->Ele[i].nabr[j] - 1 + 2 * MD->NumEle] + MD->Ele[MD->Ele[i].nabr[j] - 1].zmin);
                //              Avg_Y_Sub=avgY(MD->Ele[i].zmin,MD->Ele[MD->Ele[i].nabr[j]-1].zmin,MD->DummyY[i+2*MD->NumEle],MD->DummyY[MD->Ele[i].nabr[j]-1 + 2*MD->NumEle]);
                Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], MD->DummyY[MD->Ele[i].nabr[j] - 1 + 2 * MD->NumEle]);
                Distance = sqrt (pow ((MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x), 2) + pow ((MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y), 2));
                Grad_Y_Sub = Dif_Y_Sub / Distance;
                /*
                 * take care of macropore effect 
                 */
                effK = effKH (MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
                inabr = MD->Ele[i].nabr[j] - 1;
                nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
                effKnabr = effKH (MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
                /*
                 * It should be weighted average. However, there is an ambiguity about distance used 
                 */
                Avg_Ksat = 0.5 * (effK + effKnabr);
                /*
                 * groundwater flow modeled by Darcy's law 
                 */
                MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
                /*
                 * Surface Lateral Flux Calculation between Triangular elements Follows    
                 */
                Dif_Y_Surf = (MD->SurfMode == 1) ? (MD->Ele[i].zmax - MD->Ele[MD->Ele[i].nabr[j] - 1].zmax) : ((MD->DummyY[i] + MD->Ele[i].zmax) - (MD->DummyY[MD->Ele[i].nabr[j] - 1] + MD->Ele[MD->Ele[i].nabr[j] - 1].zmax));
                //              Avg_Y_Surf=avgY(MD->Ele[i].zmax,MD->Ele[MD->Ele[i].nabr[j] - 1].zmax,MD->DummyY[i],MD->DummyY[MD->Ele[i].nabr[j]-1]);
                Avg_Y_Surf = avgY (Dif_Y_Surf, MD->DummyY[i], MD->DummyY[MD->Ele[i].nabr[j] - 1]);
                Grad_Y_Surf = Dif_Y_Surf / Distance;
                Avg_Sf = 0.5 * (sqrt (pow (MD->Ele[i].dhBYdx, 2) + pow (MD->Ele[i].dhBYdy, 2)) + sqrt (pow (MD->Ele[MD->Ele[i].nabr[j] - 1].dhBYdx, 2) + pow (MD->Ele[MD->Ele[i].nabr[j] - 1].dhBYdy, 2))); //?? Xuan Weighting needed
                //              Avg_Sf=sqrt(pow(MD->Ele[i].dhBYdx,2)+pow(MD->Ele[i].dhBYdy,2));
                Avg_Sf = (MD->SurfMode == 1) ? (Grad_Y_Surf > 0 ? Grad_Y_Surf : EPS / pow (10.0, 6)) : (Avg_Sf > EPS / pow (10.0, 6) ? Avg_Sf : EPS / pow (10.0, 6));
                /*
                 * Weighting needed 
                 */
                Avg_Rough = 0.5 * (MD->Ele[i].Rough + MD->Ele[MD->Ele[i].nabr[j] - 1].Rough);
                CrossA = Avg_Y_Surf * MD->Ele[i].edge[j];
                OverlandFlow (MD->FluxSurf, i, j, Avg_Y_Surf, Grad_Y_Surf, Avg_Sf, CrossA, Avg_Rough);
            }

            /*
             * Boundary condition Flux Calculations Follows 
             */
            else
            {
                /*
                 * No flow (natural) boundary condition is default 
                 */
                if (MD->Ele[i].BC[j] == 0)
                {
                    MD->FluxSurf[i][j] = 0;
                    MD->FluxSub[i][j] = 0;
                }
                else if (MD->Ele[i].BC[j] == 1) /* Note: ideally different boundary conditions need to be incorporated  for surf and subsurf respectively */
                    /*
                     * Note: the formulation assumes only dirichlet TS right now 
                     */
                {
                    MD->FluxSurf[i][j] = 0; /* Note the assumption here is no flow for surface */
                    Dif_Y_Sub = (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin) - Interpolation (&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
                    Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 2 * MD->NumEle], (Interpolation (&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t) - MD->Ele[i].zmin));
                    //                          Avg_Y_Sub = (MD->DummyY[i+2*MD->NumEle] + (Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t) - MD->Ele[i].zmin))/2;
                    /*
                     * Minimum Distance from circumcenter to the edge of the triangle on which BDD. condition is defined
                     */
                    Distance = sqrt (pow (MD->Ele[i].edge[0] * MD->Ele[i].edge[1] * MD->Ele[i].edge[2] / (4 * MD->Ele[i].area), 2) - pow (MD->Ele[i].edge[j] / 2, 2));
                    effK = effKH (MD->Ele[i].Macropore, MD->DummyY[i + 2 * MD->NumEle], AquiferDepth, MD->Ele[i].macD, MD->Ele[i].macKsatH, MD->Ele[i].vAreaF, MD->Ele[i].KsatH);
                    Avg_Ksat = effK;
                    Grad_Y_Sub = Dif_Y_Sub / Distance;
                    MD->FluxSub[i][j] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * MD->Ele[i].edge[j];
                }
                else            /* Neumann BC (Note: MD->Ele[i].BC[j] value have to be = 2+(index of neumann boundary TS) */
                {
                    MD->FluxSurf[i][j] = Interpolation (&MD->TSD_EleBC[(MD->Ele[i].BC[j]) - 1], t);
                    MD->FluxSub[i][j] = Interpolation (&MD->TSD_EleBC[(-MD->Ele[i].BC[j]) - 1], t);
                }
            }
        }

        /*
         * Note: Assumption is OVL flow depth less than EPS/100 is immobile water 
         */
        if (MD->DummyY[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].infD)
        {
            /*
             * Assumption: infD<macD 
             */
            Grad_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmax - (MD->DummyY[i + 2 * MD->NumEle] + MD->Ele[i].zmin)) / MD->Ele[i].infD;
            Grad_Y_Sub = ((MD->DummyY[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
            elemSatn = 1.0;
            satKfunc = pow (elemSatn, 0.5) * pow (-1 + pow (1 - pow (elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
            effK = (MD->Ele[i].Macropore == 1) ? effKV (satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].infKsatV, MD->Ele[i].hAreaF) : MD->Ele[i].infKsatV;
#ifdef _FLUX_PIHM_
            MD->EleViR[i] = MD->EleFCR[i] * effK * Grad_Y_Sub;
#else
            MD->EleViR[i] = effK * Grad_Y_Sub;
#endif

#ifdef _FLUX_PIHM_
            if (MD->DummyY[i] + (MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area - MD->EleViR[i]) *dt  < 0)
#else
            if (MD->DummyY[i] + (MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area - MD->EleViR[i] - (MD->DummyY[i] < (EPS / 100) ? 0 : MD->EleET[i][2])) * dt < 0)
#endif
            {
#ifdef _FLUX_PIHM_
                MD->EleViR[i] = MD->DummyY[i] /dt + MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area;
#else
                MD->EleViR[i] = MD->DummyY[i] / dt + MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area - (MD->DummyY[i] < (EPS / 100) ? 0 : MD->EleET[i][2]);
#endif
                MD->EleViR[i] = MD->EleViR[i] < 0 ? 0 : MD->EleViR[i];
            }
            MD->Recharge[i] = MD->EleViR[i];
        }
        else
        {
            Deficit = AquiferDepth - MD->DummyY[i + 2 * MD->NumEle];
            //          elemSatn = elemSatn>1?1:elemSatn;
            elemSatn = 1.0;
            satKfunc = pow (elemSatn, 0.5) * pow (-1 + pow (1 - pow (elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
            /*
             * Note: for psi calculation using van genuchten relation, cutting the psi-sat tail at small saturation can be performed for computational advantage. If you dont' want to perform this, comment the statement that follows 
             */
#ifdef _FLUX_PIHM_
            elemSatn = MD->SfcSat[i];   //(MD->EleSW[i][0]-MD->Ele[i].ThetaR)/(MD->Ele[i].ThetaS - MD->Ele[i].ThetaR);
#else
            elemSatn = ((MD->DummyY[i + MD->NumEle] / Deficit) > 1) ? 1 : ((MD->DummyY[i + MD->NumEle] <= 0) ? EPS / 1000.0 : MD->DummyY[i + MD->NumEle] / Deficit);
#endif
            //          printf("elemSatn = %f, SW = %f, ThetaS = %f, ThetaR = %f, soiltype = %d\n", elemSatn, MD->EleSW[i][0], MD->Ele[i].ThetaS, MD->Ele[i].ThetaR, MD->Ele[i].geol);
            elemSatn = elemSatn > 1. ? 1. : elemSatn;
            elemSatn = (elemSatn < multF * EPS) ? (multF * EPS) : elemSatn;
            Avg_Y_Sub = (-(pow (pow (1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha) < MINpsi) ? MINpsi : (-(pow (pow (1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha));
            TotalY_Ele = Avg_Y_Sub + MD->Ele[i].zmin + AquiferDepth - MD->Ele[i].infD;
            Grad_Y_Sub = (MD->DummyY[i] + MD->Ele[i].zmax - TotalY_Ele) / MD->Ele[i].infD;
            Grad_Y_Sub = ((MD->DummyY[i] < EPS / 100) && (Grad_Y_Sub > 0)) ? 0 : Grad_Y_Sub;
            //          satKfunc = satKfunc<0.13?0.13:satKfunc;
            effK = (MD->Ele[i].Macropore == 1) ? effKV (satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].infKsatV, MD->Ele[i].hAreaF) : MD->Ele[i].infKsatV;
            //          MD->Ele[i].effKV = effK;
            //              MD->EleViR[i] = 0.5*(effK+MD->Ele[i].infKsatV)*Grad_Y_Sub;
#ifdef _FLUX_PIHM_
            MD->EleViR[i] = MD->EleFCR[i] * 0.5 * (effK) * Grad_Y_Sub;
#else
            MD->EleViR[i] = 0.5 * (effK) * Grad_Y_Sub;
#endif

#ifdef _FLUX_PIHM_
            if (MD->DummyY[i] + (MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area - MD->EleViR[i]) * dt  < 0)
#else
            if (MD->DummyY[i] + (MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area - MD->EleViR[i] - (MD->DummyY[i] < EPS / 100 ? 0 : MD->EleET[i][2])) * dt < 0)
#endif
            {
#ifdef _FLUX_PIHM_
                MD->EleViR[i] = MD->DummyY[i] / dt + MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area;
#else
                MD->EleViR[i] = MD->DummyY[i] / dt + MD->EleNetPrep[i] + (MD->FluxSurf[i][0] + MD->FluxSurf[i][1] + MD->FluxSurf[i][2]) / MD->Ele[i].area - (MD->DummyY[i] < EPS / 100 ? 0 : MD->EleET[i][2]);
#endif
                MD->EleViR[i] = MD->EleViR[i] < 0 ? 0 : MD->EleViR[i];
            }
            /*
             * Harmonic mean formulation. Note that if unsaturated zone has low saturation, satKfunc becomes very small. Use arithmetic mean instead
             */
            //                  MD->Recharge[i] = (elemSatn==0.0)?0:(Deficit<=0)?0:(MD->Ele[i].KsatV*satKfunc*(MD->Ele[i].Alpha*Deficit-2*pow(-1+pow(elemSatn,MD->Ele[i].Beta/(-MD->Ele[i].Beta+1)),1/MD->Ele[i].Beta))/(MD->Ele[i].Alpha*((Deficit+MD->DummyY[i+2*MD->NumEle]*satKfunc))));
            /*
             * Arithmetic Mean Formulation 
             */
            elemSatn = ((MD->DummyY[i + MD->NumEle] / Deficit) > 1) ? 1 : ((MD->DummyY[i + MD->NumEle] <= 0) ? (EPS / 100.0) : (MD->DummyY[i + MD->NumEle] / Deficit));
            elemSatn = (elemSatn < multF * EPS) ? (multF * EPS) : elemSatn;
            satKfunc = pow (elemSatn, 0.5) * pow (-1 + pow (1 - pow (elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)), (MD->Ele[i].Beta - 1) / MD->Ele[i].Beta), 2);
            satKfunc = satKfunc < 0.13 ? 0.13 : satKfunc;
            //          effK=(MD->Ele[i].Macropore==1)?((MD->DummyY[i+2*MD->NumEle]>AquiferDepth-MD->Ele[i].macD)?effK:(MD->Ele[i].KsatV*satKfunc)):(MD->Ele[i].KsatV*satKfunc);
            Avg_Y_Sub = (-(pow (pow (1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha) < MINpsi) ? MINpsi : (-(pow (pow (1 / elemSatn, MD->Ele[i].Beta / (MD->Ele[i].Beta - 1)) - 1, 1 / MD->Ele[i].Beta) / MD->Ele[i].Alpha));
            TotalY_Ele = Avg_Y_Sub + MD->Ele[i].zmax - 0.5 * Deficit;
            Grad_Y_Sub = (TotalY_Ele - (MD->Ele[i].zmax - Deficit)) / (0.5 * AquiferDepth);
            //(MD->DummyY[i]+MD->Ele[i].zmax-TotalY_Ele)/MD->Ele[i].infD;
            //          Grad_Y_Sub=((MD->DummyY[i]<EPS/100)&&(Grad_Y_Sub>0))?0:Grad_Y_Sub;
            effK = (MD->Ele[i].Macropore == 1) ? ((MD->DummyY[i + 2 * MD->NumEle] > AquiferDepth - MD->Ele[i].macD) ? effKV (satKfunc, Grad_Y_Sub, MD->Ele[i].macKsatV, MD->Ele[i].KsatV, MD->Ele[i].hAreaF) : (MD->Ele[i].KsatV * satKfunc)) : (MD->Ele[i].KsatV * satKfunc);

            MD->Recharge[i] = (elemSatn == 0.0) ? 0 : ((Deficit <= 0) ? 0 : (MD->Ele[i].KsatV * MD->DummyY[i + 2 * MD->NumEle] + effK * Deficit) * (MD->Ele[i].Alpha * Deficit - 2 * pow (-1 + pow (elemSatn, MD->Ele[i].Beta / (-MD->Ele[i].Beta + 1)), 1 / MD->Ele[i].Beta)) / (MD->Ele[i].Alpha * pow (Deficit + MD->DummyY[i + 2 * MD->NumEle], 2)));
            MD->Recharge[i] = (MD->Recharge[i] > 0 && MD->DummyY[i + MD->NumEle] <= 0) ? 0 : MD->Recharge[i];   //??BHATT
            MD->Recharge[i] = (MD->Recharge[i] < 0 && MD->DummyY[i + 2 * MD->NumEle] <= 0) ? 0 : MD->Recharge[i];   //??BHATT

            //          MD->EleET[i][2]=(MD->DummyY[i]<EPS/100)?elemSatn*MD->EleET[i][2]:MD->EleET[i][2];
        }
    }
    /*
     * Lateral Flux Calculation between River-River and River-Triangular elements Follows 
     */
    for (i = 0; i < MD->NumRiv; i++)
    {
        TotalY_Riv = MD->DummyY[i + 3 * MD->NumEle] + MD->Riv[i].zmin;
        Perem = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->DummyY[i + 3 * MD->NumEle], MD->Riv[i].coeff, 2);
        if (MD->Riv[i].down > 0)
        {
            /*
             * Lateral Flux Calculation between River-River element Follows 
             */
            TotalY_Riv_down = MD->DummyY[MD->Riv[i].down - 1 + 3 * MD->NumEle] + MD->Riv[MD->Riv[i].down - 1].zmin;
            Perem_down = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd, MD->DummyY[MD->Riv[i].down - 1 + 3 * MD->NumEle], MD->Riv[MD->Riv[i].down - 1].coeff, 2);
            Avg_Perem = (Perem + Perem_down) / 2.0;
            Avg_Rough = (MD->Riv_Mat[MD->Riv[i].material - 1].Rough + MD->Riv_Mat[MD->Riv[MD->Riv[i].down - 1].material - 1].Rough) / 2.0;
            Distance = 0.5 * (MD->Riv[i].Length + MD->Riv[MD->Riv[i].down - 1].Length);
            Dif_Y_Riv = (MD->RivMode == 1) ? (MD->Riv[i].zmin - MD->Riv[MD->Riv[i].down - 1].zmin) : (TotalY_Riv - TotalY_Riv_down);
            Grad_Y_Riv = Dif_Y_Riv / Distance;
            Avg_Sf = (Grad_Y_Riv > 0) ? Grad_Y_Riv : EPS;
            CrossA = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->DummyY[i + 3 * MD->NumEle], MD->Riv[i].coeff, 1);
            CrossAdown = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd, MD->DummyY[MD->Riv[i].down - 1 + 3 * MD->NumEle], MD->Riv[MD->Riv[i].down - 1].coeff, 1);
            AvgCrossA = 0.5 * (CrossA + CrossAdown);
            Avg_Y_Riv = (Avg_Perem == 0) ? 0 : (AvgCrossA / Avg_Perem);
            OverlandFlow (MD->FluxRiv, i, 1, Avg_Y_Riv, Grad_Y_Riv, Avg_Sf, CrossA, Avg_Rough);
            /*
             * accumulate to get in-flow for down segments: [0] for inflow, [1] for outflow 
             */
            MD->FluxRiv[MD->Riv[i].down - 1][0] = MD->FluxRiv[MD->Riv[i].down - 1][0] - MD->FluxRiv[i][1];

            /*
             * Lateral Flux Calculation between Element Beneath River (EBR) and EBR 
             */
            TotalY_Ele = MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv] + MD->Ele[i + MD->NumEle].zmin;
            TotalY_Ele_down = MD->DummyY[MD->Riv[i].down - 1 + 3 * MD->NumEle + MD->NumRiv] + MD->Ele[MD->Riv[i].down - 1 + MD->NumEle].zmin;
            Wid = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->Riv[i].depth, MD->Riv[i].coeff, 3);
            Wid_down = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[MD->Riv[i].down - 1].shape - 1].interpOrd, MD->Riv[MD->Riv[i].down - 1].depth, MD->Riv[MD->Riv[i].down - 1].coeff, 3);
            Avg_Wid = (Wid + Wid_down) / 2.0;
            Distance = 0.5 * (MD->Riv[i].Length + MD->Riv[MD->Riv[i].down - 1].Length);
            Dif_Y_Sub = TotalY_Ele - TotalY_Ele_down;
            //              Avg_Y_Sub=avgY(MD->Ele[i+MD->NumEle].zmin,MD->Ele[MD->Riv[i].down - 1+MD->NumEle].zmin,MD->DummyY[i + 3*MD->NumEle+MD->NumRiv],MD->DummyY[MD->Riv[i].down - 1 + 3*MD->NumEle+MD->NumRiv]);
            Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv], MD->DummyY[MD->Riv[i].down - 1 + 3 * MD->NumEle + MD->NumRiv]);
            Grad_Y_Sub = Dif_Y_Sub / Distance;
            /*
             * take care of macropore effect 
             */
            AquiferDepth = MD->Ele[i + MD->NumEle].zmax - MD->Ele[i + MD->NumEle].zmin;
            //                      effK=MD->Ele[i+MD->NumEle].KsatH;
            effK = 0.5 * (effKH (MD->Ele[MD->Riv[i].LeftEle - 1].Macropore, MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[i].LeftEle - 1].zmax - MD->Ele[MD->Riv[i].LeftEle - 1].zmin, MD->Ele[MD->Riv[i].LeftEle - 1].macD, MD->Ele[MD->Riv[i].LeftEle - 1].macKsatH, MD->Ele[MD->Riv[i].LeftEle - 1].vAreaF, MD->Ele[MD->Riv[i].LeftEle - 1].KsatH) + effKH (MD->Ele[MD->Riv[i].RightEle - 1].Macropore, MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[i].RightEle - 1].zmax - MD->Ele[MD->Riv[i].RightEle - 1].zmin, MD->Ele[MD->Riv[i].RightEle - 1].macD, MD->Ele[MD->Riv[i].RightEle - 1].macKsatH, MD->Ele[MD->Riv[i].RightEle - 1].vAreaF, MD->Ele[MD->Riv[i].RightEle - 1].KsatH));
            inabr = MD->Riv[i].down - 1;
            nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
            //                      effKnabr=MD->Ele[inabr+MD->NumEle].KsatH;
            effKnabr =
               0.5 * (effKH (MD->Ele[MD->Riv[inabr].LeftEle - 1].Macropore, MD->DummyY[MD->Riv[inabr].LeftEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[inabr].LeftEle - 1].zmax - MD->Ele[MD->Riv[inabr].LeftEle - 1].zmin, MD->Ele[MD->Riv[inabr].LeftEle - 1].macD, MD->Ele[MD->Riv[inabr].LeftEle - 1].macKsatH, MD->Ele[MD->Riv[inabr].LeftEle - 1].vAreaF, MD->Ele[MD->Riv[inabr].LeftEle - 1].KsatH) + effKH (MD->Ele[MD->Riv[inabr].RightEle - 1].Macropore, MD->DummyY[MD->Riv[inabr].RightEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[inabr].RightEle - 1].zmax - MD->Ele[MD->Riv[inabr].RightEle - 1].zmin, MD->Ele[MD->Riv[inabr].RightEle - 1].macD, MD->Ele[MD->Riv[inabr].RightEle - 1].macKsatH, MD->Ele[MD->Riv[inabr].RightEle - 1].vAreaF, MD->Ele[MD->Riv[inabr].RightEle - 1].KsatH));
            Avg_Ksat = 0.5 * (effK + effKnabr);
            /*
             * groundwater flow modeled by Darcy's law 
             */
            MD->FluxRiv[i][9] = Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub * Avg_Wid;
            /*
             * accumulate to get in-flow for down segments: [10] for inflow, [9] for outflow 
             */
            MD->FluxRiv[MD->Riv[i].down - 1][10] = MD->FluxRiv[MD->Riv[i].down - 1][10] - MD->FluxRiv[i][9];
        }
        else
        {
            switch (MD->Riv[i].down)
            {
                case -1:
                    /*
                     * Dirichlet boundary condition 
                     */
                    TotalY_Riv_down = Interpolation (&MD->TSD_Riv[(MD->Riv[i].BC) - 1], t) + (MD->Node[MD->Riv[i].ToNode - 1].zmax - MD->Riv[i].depth);
                    Distance = sqrt (pow (MD->Riv[i].x - MD->Node[MD->Riv[i].ToNode - 1].x, 2) + pow (MD->Riv[i].y - MD->Node[MD->Riv[i].ToNode - 1].y, 2));
                    Grad_Y_Riv = (TotalY_Riv - TotalY_Riv_down) / Distance;
                    /*
                     * Note: do i need to change else part here for diff wave 
                     */
                    Avg_Sf = (MD->RivMode == 1) ? Grad_Y_Riv : Grad_Y_Riv;
                    Avg_Rough = MD->Riv_Mat[MD->Riv[i].material - 1].Rough;
                    Avg_Y_Riv = avgY (Grad_Y_Riv, MD->DummyY[i + 3 * MD->NumEle], Interpolation (&MD->TSD_Riv[(MD->Riv[i].BC) - 1], t));
                    Avg_Perem = Perem;
                    CrossA = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->DummyY[i + 3 * MD->NumEle], MD->Riv[i].coeff, 1);
                    Avg_Y_Riv = (Perem == 0) ? 0 : (CrossA / Avg_Perem);
                    OverlandFlow (MD->FluxRiv, i, 1, Avg_Y_Riv, Grad_Y_Riv, Avg_Sf, CrossA, Avg_Rough);
                    break;
                case -2:
                    /*
                     * Neumann boundary condition 
                     */
                    MD->FluxRiv[i][1] = Interpolation (&MD->TSD_Riv[MD->Riv[i].BC - 1], t);
                    break;
                case -3:
                    /*
                     * zero-depth-gradient boundary conditions 
                     */
                    Distance = sqrt (pow (MD->Riv[i].x - MD->Node[MD->Riv[i].ToNode - 1].x, 2) + pow (MD->Riv[i].y - MD->Node[MD->Riv[i].ToNode - 1].y, 2));
                    Grad_Y_Riv = (MD->Riv[i].zmin - (MD->Node[MD->Riv[i].ToNode - 1].zmax - MD->Riv[i].depth)) / Distance;
                    Avg_Rough = MD->Riv_Mat[MD->Riv[i].material - 1].Rough;
                    Avg_Y_Riv = MD->DummyY[i + 3 * MD->NumEle];
                    Avg_Perem = Perem;
                    CrossA = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->DummyY[i + 3 * MD->NumEle], MD->Riv[i].coeff, 1);
                    MD->FluxRiv[i][1] = sqrt (Grad_Y_Riv) * CrossA * ((Avg_Perem > 0) ? pow (CrossA / Avg_Perem, 2.0 / 3.0) : 0) / Avg_Rough;
                    break;
                case -4:
                    /*
                     * Critical Depth boundary conditions 
                     */
                    CrossA = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->DummyY[i + 3 * MD->NumEle], MD->Riv[i].coeff, 1);
                    MD->FluxRiv[i][1] = CrossA * sqrt (GRAV * MD->DummyY[i + 3 * MD->NumEle]);    /* Note the dependence on physical units */
                    break;
                default:
                    printf ("Fatal Error: River Routing Boundary Condition Type Is Wrong!");
                    exit (1);
            }
            /*
             * Note: bdd condition for subsurface element can be changed. Assumption: No flow condition 
             */
            MD->FluxRiv[i][9] = 0;
        }
        if (MD->Riv[i].LeftEle > 0)
        {

            /*
             * Lateral Surface Flux Calculation between River-Triangular element Follows 
             */
            OLFeleToriv (MD->DummyY[MD->Riv[i].LeftEle - 1] + MD->Ele[MD->Riv[i].LeftEle - 1].zmax, MD->Ele[MD->Riv[i].LeftEle - 1].zmax, MD->Riv_Mat[MD->Riv[i].material - 1].Cwr, MD->Riv[i].zmax, TotalY_Riv, MD->FluxRiv, i, 2, MD->Riv[i].Length);

            /*
             * Lateral Sub-surface Flux Calculation between River-Triangular element Follows 
             */
            Dif_Y_Sub = (MD->DummyY[i + 3 * MD->NumEle] + MD->Riv[i].zmin) - (MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle - 1].zmin);
            //          Avg_Y_Sub=(MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle-1].zmin-MD->Riv[i].zmin)>0?MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle-1].zmin-MD->Riv[i].zmin:0;
            /*
             * This is head at river edge representation 
             */
            //          Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])-MD->Riv[i].zmin):0;
            /*
             * This is head in neighboring cell represention 
             */
            Avg_Y_Sub = MD->Ele[MD->Riv[i].LeftEle - 1].zmin > MD->Riv[i].zmin ? MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle] : ((MD->Ele[MD->Riv[i].LeftEle - 1].zmin + MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle]) > MD->Riv[i].zmin ? (MD->Ele[MD->Riv[i].LeftEle - 1].zmin + MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle] - MD->Riv[i].zmin) : 0);
            //          Avg_Y_Sub=avgY(MD->Riv[i].zmin,MD->Riv[i].zmin,MD->DummyY[i+3*MD->NumEle],Avg_Y_Sub);
            Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 3 * MD->NumEle], Avg_Y_Sub);
            effK = MD->Riv[i].KsatH;
            Distance = sqrt (pow ((MD->Riv[i].x - MD->Ele[MD->Riv[i].LeftEle - 1].x), 2) + pow ((MD->Riv[i].y - MD->Ele[MD->Riv[i].LeftEle - 1].y), 2));
            Grad_Y_Sub = Dif_Y_Sub / Distance;
            /*
             * take care of macropore effect 
             */
            inabr = MD->Riv[i].LeftEle - 1;
            AquiferDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
            effKnabr = effKH (MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], AquiferDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
            Avg_Ksat = 0.5 * (effK + effKnabr);
            MD->FluxRiv[i][4] = MD->Riv[i].Length * Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub;

            /*
             * Lateral Flux between rectangular element (beneath river) and triangular element 
             */
            Dif_Y_Sub = (MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv] + MD->Ele[i + MD->NumEle].zmin) - (MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle] + MD->Ele[MD->Riv[i].LeftEle - 1].zmin);
            //          Avg_Y_Sub=((MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].LeftEle-1].zmin-MD->Riv[i].zmin)>0)?MD->Riv[i].zmin-MD->Ele[MD->Riv[i].LeftEle-1].zmin:MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle];
            /*
             * This is head at river edge representation 
             */
            //          Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?MD->Riv[i].zmin-(MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)):MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle];
            /*
             * This is head in neighboring cell represention 
             */
            Avg_Y_Sub = MD->Ele[MD->Riv[i].LeftEle - 1].zmin > MD->Riv[i].zmin ? 0 : ((MD->Ele[MD->Riv[i].LeftEle - 1].zmin + MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle]) > MD->Riv[i].zmin ? (MD->Riv[i].zmin - MD->Ele[MD->Riv[i].LeftEle - 1].zmin) : MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle]);
            //          Avg_Y_Sub=avgY(MD->Ele[i+MD->NumEle].zmin,MD->Ele[MD->Riv[i].LeftEle-1].zmin,MD->DummyY[i+3*MD->NumEle+MD->NumRiv],Avg_Y_Sub); 
            Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv], Avg_Y_Sub);
            AquiferDepth = (MD->Ele[i + MD->NumEle].zmax - MD->Ele[i + MD->NumEle].zmin);
            //          effK=MD->Ele[i+MD->NumEle].KsatH;
            effK = 0.5 * (effKH (MD->Ele[MD->Riv[i].LeftEle - 1].Macropore, MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[i].LeftEle - 1].zmax - MD->Ele[MD->Riv[i].LeftEle - 1].zmin, MD->Ele[MD->Riv[i].LeftEle - 1].macD, MD->Ele[MD->Riv[i].LeftEle - 1].macKsatH, MD->Ele[MD->Riv[i].LeftEle - 1].vAreaF, MD->Ele[MD->Riv[i].LeftEle - 1].KsatH) + effKH (MD->Ele[MD->Riv[i].RightEle - 1].Macropore, MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[i].RightEle - 1].zmax - MD->Ele[MD->Riv[i].RightEle - 1].zmin, MD->Ele[MD->Riv[i].RightEle - 1].macD, MD->Ele[MD->Riv[i].RightEle - 1].macKsatH, MD->Ele[MD->Riv[i].RightEle - 1].vAreaF, MD->Ele[MD->Riv[i].RightEle - 1].KsatH));
            inabr = MD->Riv[i].LeftEle - 1;
            nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
            effKnabr = effKH (MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
            Avg_Ksat = 0.5 * (effK + effKnabr);
            Grad_Y_Sub = Dif_Y_Sub / Distance;  /* take care of macropore effect */
            MD->FluxRiv[i][7] = MD->Riv[i].Length * Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub;

            /*
             * replace flux term 
             */
            for (j = 0; j < 3; j++)
            {
                if (MD->Ele[MD->Riv[i].LeftEle - 1].nabr[j] == MD->Riv[i].RightEle)
                {
                    if (-MD->FluxRiv[i][2] > 0 && -MD->FluxRiv[i][2] > MD->FluxSurf[MD->Riv[i].LeftEle - 1][j])
                    {
                        MD->FluxRiv[i][2] = -MD->DummyY[MD->Riv[i].LeftEle - 1] / dt;
                    }
                    MD->FluxSurf[MD->Riv[i].LeftEle - 1][j] = -MD->FluxRiv[i][2];
                    MD->FluxSub[MD->Riv[i].LeftEle - 1][j] = -MD->FluxRiv[i][4];
                    MD->FluxSub[MD->Riv[i].LeftEle - 1][j] = MD->FluxSub[MD->Riv[i].LeftEle - 1][j] - MD->FluxRiv[i][7];
                    break;
                }
            }
        }
        if (MD->Riv[i].RightEle > 0)
        {

            /*
             * Lateral Surface Flux Calculation between River-Triangular element Follows 
             */
            OLFeleToriv (MD->DummyY[MD->Riv[i].RightEle - 1] + MD->Ele[MD->Riv[i].RightEle - 1].zmax, MD->Ele[MD->Riv[i].RightEle - 1].zmax, MD->Riv_Mat[MD->Riv[i].material - 1].Cwr, MD->Riv[i].zmax, TotalY_Riv, MD->FluxRiv, i, 3, MD->Riv[i].Length);

            /*
             * Lateral Sub-surface Flux Calculation between River-Triangular element Follows 
             */
            Dif_Y_Sub = (MD->DummyY[i + 3 * MD->NumEle] + MD->Riv[i].zmin) - (MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle] + MD->Ele[MD->Riv[i].RightEle - 1].zmin);
            //          Avg_Y_Sub=(MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle-1].zmin-MD->Riv[i].zmin>0)?MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle-1].zmin-MD->Riv[i].zmin:0;
            /*
             * This is head at river edge representation 
             */
            //          Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])-MD->Riv[i].zmin):0;
            /*
             * This is head in neighboring cell represention 
             */
            Avg_Y_Sub = MD->Ele[MD->Riv[i].RightEle - 1].zmin > MD->Riv[i].zmin ? MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle] : ((MD->Ele[MD->Riv[i].RightEle - 1].zmin + MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle]) > MD->Riv[i].zmin ? (MD->Ele[MD->Riv[i].RightEle - 1].zmin + MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle] - MD->Riv[i].zmin) : 0);
            //          Avg_Y_Sub=avgY(MD->Riv[i].zmin,MD->Riv[i].zmin,MD->DummyY[i+3*MD->NumEle],Avg_Y_Sub);
            Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 3 * MD->NumEle], Avg_Y_Sub);
            effK = MD->Riv[i].KsatH;
            Distance = sqrt (pow ((MD->Riv[i].x - MD->Ele[MD->Riv[i].RightEle - 1].x), 2) + pow ((MD->Riv[i].y - MD->Ele[MD->Riv[i].RightEle - 1].y), 2));
            Grad_Y_Sub = Dif_Y_Sub / Distance;
            /*
             * take care of macropore effect 
             */
            inabr = MD->Riv[i].RightEle - 1;
            AquiferDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
            effKnabr = effKH (MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], AquiferDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
            Avg_Ksat = 0.5 * (effK + effKnabr);
            MD->FluxRiv[i][5] = MD->Riv[i].Length * Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub;

            /*
             * Lateral Flux between rectangular element (beneath river) and triangular element 
             */
            Dif_Y_Sub = (MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv] + MD->Ele[i + MD->NumEle].zmin) - (MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle] + MD->Ele[MD->Riv[i].RightEle - 1].zmin);
            //          Avg_Y_Sub=((MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle]+MD->Ele[MD->Riv[i].RightEle-1].zmin-MD->Riv[i].zmin)>0)?MD->Riv[i].zmin-MD->Ele[MD->Riv[i].RightEle-1].zmin:MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle];
            /*
             * This is head at river edge representation 
             */
            //          Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)+MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?MD->Riv[i].zmin-(MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].RightEle-1].zmax-MD->Ele[MD->Riv[i].RightEle-1].zmin)):MD->DummyY[MD->Riv[i].RightEle-1 + 2*MD->NumEle];
            /*
             * This is head in neighboring cell represention 
             */
            Avg_Y_Sub = MD->Ele[MD->Riv[i].RightEle - 1].zmin > MD->Riv[i].zmin ? 0 : ((MD->Ele[MD->Riv[i].RightEle - 1].zmin + MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle]) > MD->Riv[i].zmin ? (MD->Riv[i].zmin - MD->Ele[MD->Riv[i].RightEle - 1].zmin) : MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle]);
            //          Avg_Y_Sub=avgY(MD->Ele[i+MD->NumEle].zmin,MD->Ele[MD->Riv[i].RightEle-1].zmin,MD->DummyY[i+3*MD->NumEle+MD->NumRiv],Avg_Y_Sub); 
            Avg_Y_Sub = avgY (Dif_Y_Sub, MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv], Avg_Y_Sub);
            AquiferDepth = (MD->Ele[i + MD->NumEle].zmax - MD->Ele[i + MD->NumEle].zmin);
            //          effK=MD->Ele[i+MD->NumEle].KsatH;
            effK = 0.5 * (effKH (MD->Ele[MD->Riv[i].LeftEle - 1].Macropore, MD->DummyY[MD->Riv[i].LeftEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[i].LeftEle - 1].zmax - MD->Ele[MD->Riv[i].LeftEle - 1].zmin, MD->Ele[MD->Riv[i].LeftEle - 1].macD, MD->Ele[MD->Riv[i].LeftEle - 1].macKsatH, MD->Ele[MD->Riv[i].LeftEle - 1].vAreaF, MD->Ele[MD->Riv[i].LeftEle - 1].KsatH) + effKH (MD->Ele[MD->Riv[i].RightEle - 1].Macropore, MD->DummyY[MD->Riv[i].RightEle - 1 + 2 * MD->NumEle], MD->Ele[MD->Riv[i].RightEle - 1].zmax - MD->Ele[MD->Riv[i].RightEle - 1].zmin, MD->Ele[MD->Riv[i].RightEle - 1].macD, MD->Ele[MD->Riv[i].RightEle - 1].macKsatH, MD->Ele[MD->Riv[i].RightEle - 1].vAreaF, MD->Ele[MD->Riv[i].RightEle - 1].KsatH));
            inabr = MD->Riv[i].RightEle - 1;
            nabrAqDepth = (MD->Ele[inabr].zmax - MD->Ele[inabr].zmin);
            effKnabr = effKH (MD->Ele[inabr].Macropore, MD->DummyY[inabr + 2 * MD->NumEle], nabrAqDepth, MD->Ele[inabr].macD, MD->Ele[inabr].macKsatH, MD->Ele[inabr].vAreaF, MD->Ele[inabr].KsatH);
            Avg_Ksat = 0.5 * (effK + effKnabr);
            Grad_Y_Sub = Dif_Y_Sub / Distance;  /* take care of macropore effect */
            MD->FluxRiv[i][8] = MD->Riv[i].Length * Avg_Ksat * Grad_Y_Sub * Avg_Y_Sub;
            /*
             * replace flux item 
             */
            for (j = 0; j < 3; j++)
            {
                if (MD->Ele[MD->Riv[i].RightEle - 1].nabr[j] == MD->Riv[i].LeftEle)
                {
                    if (-MD->FluxRiv[i][3] > 0 && -MD->FluxRiv[i][3] > MD->FluxSurf[MD->Riv[i].RightEle - 1][j])
                    {
                        MD->FluxRiv[i][3] = -MD->DummyY[MD->Riv[i].RightEle - 1] / dt;
                    }
                    MD->FluxSurf[MD->Riv[i].RightEle - 1][j] = -MD->FluxRiv[i][3];
                    MD->FluxSub[MD->Riv[i].RightEle - 1][j] = -MD->FluxRiv[i][5];
                    MD->FluxSub[MD->Riv[i].RightEle - 1][j] = MD->FluxSub[MD->Riv[i].RightEle - 1][j] - MD->FluxRiv[i][8];
                    break;
                }
            }
        }
        Avg_Wid = CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd, MD->DummyY[i + 3 * MD->NumEle], MD->Riv[i].coeff, 3);
        Dif_Y_Riv = (MD->Riv[i].zmin - (MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv] + MD->Ele[i + MD->NumEle].zmin)) > 0 ? MD->DummyY[i + 3 * MD->NumEle] : MD->DummyY[i + 3 * MD->NumEle] + MD->Riv[i].zmin - (MD->DummyY[i + 3 * MD->NumEle + MD->NumRiv] + MD->Ele[i + MD->NumEle].zmin);
        Grad_Y_Riv = Dif_Y_Riv / MD->Riv[i].bedThick;
        MD->FluxRiv[i][6] = MD->Riv[i].KsatV * Avg_Wid * MD->Riv[i].Length * Grad_Y_Riv;
    }
}
