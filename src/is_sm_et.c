/********************************************************************************
 * File        : is_sm_et.c                                                     *
 * Function    : for calculation of canopy ET, IS and snow melt       	        *
 * Version     : Nov, 2007 (2.0)                                                *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                   *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)              *
 *------------------------------------------------------------------------------*
 *                                                                              *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0.............................*
 * a) Modification of ET components from canopy, ground and transpiration       *
 * b) Addition of snow melt and throughflow processes                           *
 * c) Fully Coupled Inclusion of ET from transpiration and OVL 		        *
 * d) Rain/Snow Fraction calculation according to USACE (1956)		        *
 * e) Change in maximum interception storage due to snow accretion has been     *
 *	accounted for								*
 * f) Incorporation of Interception storage for rainfall as well as snow	*
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "pihm.h"

#define multF1	1
#define multF2	5//??
#define multF3	1
#define multF4	1.0
#define C_air 1004.0
#define Lv (2.503*pow(10,6))
#define SIGMA (5.67*pow(10,-8)*60)
#define R_dry 287.04
#define R_v 461.5


void is_sm_et(realtype t, realtype stepsize, void *DS, N_Vector VY)
{
	int             i;
	realtype        totEvap;
	realtype	elemSatn;
	realtype	AquiferDepth;
	realtype	ThetaRef, ThetaW;
	realtype	beta_s, Rmax, f_r, alpha_r, eta_s, gamma_s, r_s, P_c;
	realtype        Delta, Gamma, Lambda;
	realtype        Rn, G, T, Vel, RH, VP, P, LAI, zero_dh, cnpy_h, rl,
	                r_a, qv_sat, qv, ETp;
	realtype        isval = 0, etval = 0;
	realtype        fracSnow, snowRate, MeltRateGrnd, MeltRateCanopy, eltRate,
	                MF, Ts = -3.0, Tr = 1.0, To = 0.0, ret;

	Model_Data      MD;

	MD = (Model_Data) DS;

//	stepsize = stepsize / UNIT_C;
	for (i = 0; i < MD->NumEle; i++)
	{
		/* Note the dependence on physical units */
		MD->ElePrep[i] = Interpolation(&MD->Forcing[0][MD->Ele[i].prep - 1], t) / 1000.;
//		Rn = Interpolation(&MD->Forcing[4][MD->Ele[i].Rn - 1], t);
		Rn = Interpolation(&MD->Forcing[4][MD->Ele[i].Sdown - 1], t);
		//G = Interpolation(&MD->TSD_G[MD->Ele[i].G - 1], t);
		//G = 0.1 * Rn;
		T = Interpolation(&MD->Forcing[1][MD->Ele[i].temp - 1], t) - 273.15;
		Vel = Interpolation(&MD->Forcing[3][MD->Ele[i].WindVel - 1], t);
		RH = Interpolation(&MD->Forcing[2][MD->Ele[i].humidity - 1], t) / 100.;
		//VP = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure - 1], t);
		VP = 611.2 * exp(17.67 * T / (T + 243.5)) * RH;
		P = 101.325 * pow(10, 3) * pow((293 - 0.0065 * MD->Ele[i].zmax) / 293, 5.26);
		qv = 0.622 * VP / P;
		qv_sat = 0.622 * (VP / RH) / P;
		LAI = Interpolation(&MD->Forcing[7][MD->Ele[i].LC - 1], t);

		MF = multF2 * Interpolation(&MD->Forcing[9][MD->Ele[i].meltF - 1], t);
		/******************************************************************************************/
		/* Snow Accumulation/Melt Calculation				  */
		/******************************************************************************************/
		fracSnow = T < Ts ? 1.0 : T > Tr ? 0 : (Tr - T) / (Tr - Ts);
		snowRate = fracSnow * MD->ElePrep[i];
		/*
		 * EleSnowGrnd, EleSnowCanopy, EleISsnowmax,
		 * MeltRateGrnd,MeltRateCanopy are the average value prorated
		 * over the whole elemental area
		 */
		MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] + (1 - MD->Ele[i].VegFrac) * snowRate * stepsize;
		MD->EleSnowCanopy[i] = MD->EleSnowCanopy[i] + MD->Ele[i].VegFrac * snowRate * stepsize;
		MD->EleISsnowmax[i] = MD->EleSnowCanopy[i] > 0 ? 0.003 * LAI * MD->Ele[i].VegFrac : 0;
		MD->EleISsnowmax[i] = multF1 * MD->EleISsnowmax[i];
		if (MD->EleSnowCanopy[i] > MD->EleISsnowmax[i])
		{
			MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] + MD->EleSnowCanopy[i] - MD->EleISsnowmax[i];
			MD->EleSnowCanopy[i] = MD->EleISsnowmax[i]; /* update in 2.2  */
		}
		MeltRateGrnd = MeltRateCanopy = (T > To ? (T - To) * MF : 0);	/* Note the units for
										 * MF. */
		if (MD->EleSnowGrnd[i] > MeltRateGrnd * stepsize)
		{
			MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] - MeltRateGrnd * stepsize;
		}
		else
		{
			MeltRateGrnd = MD->EleSnowGrnd[i] / stepsize;
			MD->EleSnowGrnd[i] = 0;
		}
		if (MD->EleSnowCanopy[i] > MeltRateCanopy * stepsize)
		{
			MD->EleSnowCanopy[i] = MD->EleSnowCanopy[i] - MeltRateCanopy * stepsize;
		}
		else
		{
			MeltRateCanopy = MD->EleSnowCanopy[i] / stepsize;
			MD->EleSnowCanopy[i] = 0;
		}
		/************************************************************************/
		/* ThroughFall and Evaporation from canopy			 */
		/************************************************************************/
		/*
		 * EleIS, EleET[0] and ret are prorated for the whole
		 * element. Logistics are simpler if assumed in volumetric
		 * form by multiplication of Area on either side of equation
		 */
		MD->EleISmax[i] = multF1 * MD->ISFactor[MD->Ele[i].LC - 1] * LAI * MD->Ele[i].VegFrac;

		rl = Interpolation(&MD->Forcing[8][MD->Ele[i].LC - 1], t);
		r_a = log(MD->Ele[i].windH / rl) * log(10 * MD->Ele[i].windH / rl) / (Vel * 0.16);
		//r_a = 12 * 4.72 * log(MD->Ele[i].windH / rl) / (0.54 * Vel / UNIT_C / 60 + 1) / UNIT_C / 60;

		Gamma = 4 * 0.7 * SIGMA * R_dry / C_air * pow(T + 273.15, 4) / (P / r_a) + 1;
		Delta = Lv * Lv * 0.622 / R_v / C_air / pow(T + 273.15, 2) * qv_sat;

		ETp = (Rn * Delta + Gamma * (1.2 * Lv * (qv_sat - qv) / r_a)) / (1000.0 * Lv * (Delta + Gamma));

		if ((MD->Ele[i].zmax - MD->Ele[i].zmin) - MD->EleGW[i] < MD->Ele[i].RzD)
			elemSatn = 1.0;
		else
			elemSatn = ((MD->EleUnsat[i] / (AquiferDepth - MD->EleGW[i])) > 1.) ? 1. : ((MD->EleUnsat[i] / (AquiferDepth - MD->EleGW[i])) < 0) ? 0 : 0.5 * (1 - cos(3.14 * (MD->EleUnsat[i] / (AquiferDepth - MD->EleGW[i]))));

		beta_s = (elemSatn * MD->Ele[i].Porosity + MD->Ele[i].ThetaR - MD->Ele[i].ThetaW) / (MD->Ele[i].ThetaRef - MD->Ele[i].ThetaW);
		beta_s = (beta_s < 0.0001) ? 0.0001 : (beta_s > 1 ? 1 : beta_s);
		MD->EleET[i][2] = MD->pcCal.Et2 * (1 - MD->Ele[i].VegFrac) * pow(beta_s, 2.) * ETp;
		MD->EleET[i][2] = MD->EleET[i][2] < 0 ? 0 : MD->EleET[i][2];

		/* Note the dependence on physical units */
		if (LAI > 0.0)
		{

			MD->EleET[i][0] = MD->pcCal.Et0 * MD->Ele[i].VegFrac * (pow((MD->EleIS[i] < 0 ? 0 : (MD->EleIS[i] > MD->EleISmax[i] ? MD->EleISmax[i] : MD->EleIS[i])) / MD->EleISmax[i], 1.0 / 2.0)) * ETp;
			MD->EleET[i][0] = MD->EleET[i][0] < 0 ? 0 : MD->EleET[i][0];

			Rmax = MD->Rmax;	/* Unit day_per_m */
			f_r = 1.1 * Rn / (MD->Ele[i].Rs_ref * LAI);
			f_r = f_r < 0 ? 0 : f_r;
			alpha_r = (1 + f_r) / (f_r + (MD->Ele[i].Rmin / Rmax));
			alpha_r = alpha_r > 10000 ? 10000 : alpha_r;
			eta_s = 1 - 0.0016 * (pow((MD->Tref - 273.15 - T), 2));
			eta_s = eta_s < 0.0001 ? 0.0001 : eta_s;
			gamma_s = 1 / (1 + 0.00025 * (VP / RH - VP));
			gamma_s = (gamma_s < 0.01) ? 0.01 : gamma_s;
			r_s = ((MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s)) > Rmax) ? Rmax : (MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s));
			P_c = (1 + Delta / Gamma) / (1 + r_s / r_a + Delta / Gamma);
			MD->EleET[i][1] = MD->pcCal.Et1 * MD->Ele[i].VegFrac * P_c * (1 - pow(((MD->EleIS[i] + MD->EleSnowCanopy[i] < 0) ? 0 : (MD->EleIS[i] + MD->EleSnowCanopy[i])) / (MD->EleISmax[i] + MD->EleISsnowmax[i]), MD->fx_canopy)) * ETp;
			MD->EleET[i][1] = MD->EleET[i][1] < 0 ? 0 : MD->EleET[i][1];
			AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
			MD->EleET[i][1] = ((MD->DummyY[i + 2 * MD->NumEle] < (AquiferDepth - MD->Ele[i].RzD)) && MD->DummyY[i + MD->NumEle] <= 0) ? 0 : MD->EleET[i][1];

			MD->EleTF[i] = MD->EleIS[i] <= 0 ? 0 : 5.65 * pow(10, -2) * MD->EleISmax[i] * exp(3.89 * (MD->EleIS[i] < 0 ? 0 : MD->EleIS[i]) / MD->EleISmax[i]);	/* Note the dependece on  physical units */
			MD->EleTF[i] = multF3 * MD->EleTF[i];
		}
		else
		{
			MD->EleET[i][1] = 0.0;
			MD->EleET[i][0] = 0.0;
			MD->EleTF[i] = 0.0;
		}

		if(MD->EleTF[i]<0)MD->EleTF[i] = 0.0;
		if(MD->EleTF[i]*stepsize>MD->EleIS[i])MD->EleTF[i]=MD->EleIS[i]/stepsize;
		if (MD->EleIS[i] >= MD->EleISmax[i]) {
			if (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) >= MD->EleET[i][0] + MD->EleTF[i]) {
//				MD->EleETloss[i] = MD->EleET[i][0];
				ret = MD->EleTF[i] + (MD->EleIS[i] - MD->EleISmax[i])/stepsize + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - (MD->EleET[i][0] + MD->EleTF[i]));
				isval = MD->EleISmax[i];
				//MD->EleIS[i] = MD->EleISmax[i];
			} else if ((((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) < MD->EleET[i][0] + MD->EleTF[i]) && (MD->EleIS[i] + stepsize * ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy - MD->EleET[i][0] - MD->EleTF[i]) <= 0)) {
				MD->EleET[i][0] = (MD->EleET[i][0] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
				ret = (MD->EleTF[i] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
//				MD->EleET[i][0]=MD->EleETloss[i];
				//MD->EleIS[i] = 0;
				isval = 0;
//				MD->EleETloss[i] = MD->EleET[i][0];
			} else {
				isval = MD->EleIS[i] + stepsize * (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]);
				//MD->EleIS[i] = MD->EleIS[i] + stepsize * (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]);
				ret = MD->EleTF[i];
//				MD->EleETloss[i] = MD->EleET[i][0];
			}
		} else if ((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) >= MD->EleISmax[i])) {
//			MD->EleETloss[i] = MD->EleET[i][0];
			isval = MD->EleISmax[i];
			ret = MD->EleTF[i] + (((MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) - MD->EleISmax[i]))/stepsize;
		} else if ((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) <= 0)) {
			if ((MD->EleET[i][0] > 0) || (MD->EleTF[i] > 0)) {
				MD->EleET[i][0] = (MD->EleET[i][0] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
				ret = (MD->EleTF[i] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
//				MD->EleET[i][0]=MD->EleETloss[i];
			} else {
				MD->EleET[i][0] = 0;
				ret = 0;
			}
//			MD->EleETloss[i] = MD->EleET[i][0];
			isval = 0;
		} else {
			isval = MD->EleIS[i] + (((1 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize;
//			MD->EleETloss[i] = MD->EleET[i][0];
			ret = MD->EleTF[i];
		}
		MD->EleNetPrep[i] = (1 - MD->Ele[i].VegFrac) * (1 - fracSnow) * MD->ElePrep[i] + ret + MeltRateGrnd;
		MD->EleTF[i] = ret;
		MD->EleIS[i] = isval;
		//MD->EleNetPrep[i] = MD->ElePrep[i];
	}
}
