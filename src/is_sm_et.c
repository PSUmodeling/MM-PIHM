/********************************************************************************
 * File        : is_sm_et.c                                                     *
 * Function    : for calculation of canopy ET, IS and snow melt       	        *
 *------------------------------------------------------------------------------*
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0.............................*
 * a) Modification of ET components from canopy, ground and transpiration       *
 * b) Addition of snow melt and throughflow processes                           *
 * c) Fully Coupled Inclusion of ET from transpiration and OVL 		        *
 * d) Rain/Snow Fraction calculation according to USACE (1956)		        *
 * e) Change in maximum interception storage due to snow accretion has been     *
 *	accounted for								*
 * f) Incorporation of Interception storage for rainfall as well as snow	*
 ********************************************************************************/

#include "pihm.h"

void is_sm_et(realtype t, realtype stepsize, void *DS, N_Vector VY)
{
	int             i;
	realtype	elemSatn;
	realtype	AquiferDepth;
	realtype	beta_s;
        realtype        Rmax;
        realtype        f_r;
        realtype        alpha_r;
        realtype        eta_s;
        realtype        gamma_s;
        realtype        r_s;
        realtype        P_c;
	realtype        Delta, Gamma;
	realtype        Rn;
        realtype        T;
        realtype        Vel;
        realtype        RH;
        realtype        VP;
        realtype        P;
        realtype        LAI;
	realtype	rl;
        realtype        r_a;
        realtype        qv_sat;
        realtype        qv;
        realtype        ETp;
	realtype        isval = 0;
	realtype        fracSnow;
        realtype        snowRate;
        realtype        MeltRateGrnd;
        realtype        MeltRateCanopy;
	realtype	MF;
        realtype        Ts = -3.0;
        realtype        Tr = 1.0;
        realtype        To = 0.0;
        realtype        ret;
        realtype        *metarr;

	Model_Data      MD;

	MD = (Model_Data) DS;

        metarr = (realtype *) malloc (7 * sizeof (realtype));
	for (i = 0; i < MD->NumEle; i++)
        {
            MultiInterpolation (&MD->TSD_meteo[MD->Ele[i].meteo - 1], t, metarr, 7);

            /* Note the dependence on physical units */
            MD->ElePrep[i] = metarr[PRCP_TS] / 1000.0;
            Rn = metarr[SOLAR_TS];
            T = metarr[SFCTMP_TS] - 273.15;
            Vel = metarr[SFCSPD_TS];
            RH = metarr[RH_TS] / 100.0;

            VP = 611.2 * exp (17.67 * T / (T + 243.5)) * RH;
            P = 101.325 * 1.0e3 * pow((293.0 - 0.0065 * MD->Ele[i].zmax) / 293.0, 5.26);
            qv = 0.622 * VP / P;
            qv_sat = 0.622 * (VP / RH) / P;
            if (MD->Ele[i].LAI > 0)
                LAI = Interpolation(&MD->TSD_lai[MD->Ele[i].LAI - 1], t);
            else
                LAI = monthly_lai (t, MD->Ele[i].LC);

            MF = monthly_mf (t);

            /* Snow Accumulation/Melt Calculation */
            fracSnow = T < Ts ? 1.0 : T > Tr ? 0 : (Tr - T) / (Tr - Ts);
            snowRate = fracSnow * MD->ElePrep[i];
            /* EleSnowGrnd, EleSnowCanopy, EleISsnowmax, MeltRateGrnd,
             * MeltRateCanopy are the average value prorated over the whole
             * elemental area */
            MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] + (1.0 - MD->Ele[i].VegFrac) * snowRate * stepsize;
            MD->EleSnowCanopy[i] = MD->EleSnowCanopy[i] + MD->Ele[i].VegFrac * snowRate * stepsize;
            MD->EleISsnowmax[i] = MD->EleSnowCanopy[i] > 0.0 ? 0.003 * LAI * MD->Ele[i].VegFrac : 0.0;
            MD->EleISsnowmax[i] = MD->EleISsnowmax[i];
            if (MD->EleSnowCanopy[i] > MD->EleISsnowmax[i])
            {
                MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] + MD->EleSnowCanopy[i] - MD->EleISsnowmax[i];
                MD->EleSnowCanopy[i] = MD->EleISsnowmax[i];
            }
            MeltRateGrnd = MeltRateCanopy = (T > To ? (T - To) * MF : 0.0);

            if (MD->EleSnowGrnd[i] > MeltRateGrnd * stepsize)
                MD->EleSnowGrnd[i] = MD->EleSnowGrnd[i] - MeltRateGrnd * stepsize;
            else
            {
                MeltRateGrnd = MD->EleSnowGrnd[i] / stepsize;
                MD->EleSnowGrnd[i] = 0.0;
            }
            if (MD->EleSnowCanopy[i] > MeltRateCanopy * stepsize)
                MD->EleSnowCanopy[i] = MD->EleSnowCanopy[i] - MeltRateCanopy * stepsize;
            else
            {
                MeltRateCanopy = MD->EleSnowCanopy[i] / stepsize;
                MD->EleSnowCanopy[i] = 0.0;
            }
            /* ThroughFall and Evaporation from canopy */
            /*
             * EleIS, EleET[0] and ret are prorated for the whole element.
             * Logistics are simpler if assumed in volumetric form by
             * multiplication of Area on either side of equation */
            MD->EleISmax[i] = MD->ISFactor[MD->Ele[i].LC - 1] * LAI * MD->Ele[i].VegFrac;

            rl = monthly_rl (t, MD->Ele[i].LC);

            r_a = log (MD->Ele[i].windH / rl) * log (10.0 * MD->Ele[i].windH / rl) / (Vel * 0.16);

            Gamma = 4.0 * 0.7 * SIGMA * R_dry / C_air * pow (T + 273.15, 4) / (P / r_a) + 1.0;
            Delta = Lv * Lv * 0.622 / R_v / C_air / pow (T + 273.15, 2) * qv_sat;

            ETp = (Rn * Delta + Gamma * (1.2 * Lv * (qv_sat - qv) / r_a)) / (1000.0 * Lv * (Delta + Gamma));

            if ((MD->Ele[i].zmax - MD->Ele[i].zmin) - MD->EleGW[i] < MD->Ele[i].RzD)
                elemSatn = 1.0;
            else
                elemSatn = ((MD->EleUnsat[i] / (AquiferDepth - MD->EleGW[i])) > 1.0) ? 1.0 : ((MD->EleUnsat[i] / (AquiferDepth - MD->EleGW[i])) < 0.0) ? 0.0 : 0.5 * (1.0 - cos (3.14 * (MD->EleUnsat[i] / (AquiferDepth - MD->EleGW[i]))));

            beta_s = (elemSatn * MD->Ele[i].Porosity + MD->Ele[i].ThetaR - MD->Ele[i].ThetaW) / (MD->Ele[i].ThetaRef - MD->Ele[i].ThetaW);
            beta_s = (beta_s < 0.0001) ? 0.0001 : (beta_s > 1.0 ? 1.0 : beta_s);
            MD->EleET[i][2] = MD->pcCal.Et2 * (1.0 - MD->Ele[i].VegFrac) * pow(beta_s, 2) * ETp;
            MD->EleET[i][2] = MD->EleET[i][2] < 0.0 ? 0.0 : MD->EleET[i][2];

            /* Note the dependence on physical units */
            if (LAI > 0.0)
            {
                MD->EleET[i][0] = MD->pcCal.Et0 * MD->Ele[i].VegFrac * (pow ((MD->EleIS[i] < 0.0 ? 0.0 : (MD->EleIS[i] > MD->EleISmax[i] ? MD->EleISmax[i] : MD->EleIS[i])) / MD->EleISmax[i], 1.0 / 2.0)) * ETp;
                MD->EleET[i][0] = MD->EleET[i][0] < 0.0 ? 0.0 : MD->EleET[i][0];

                Rmax = MD->Rmax;
                f_r = 1.1 * Rn / (MD->Ele[i].Rs_ref * LAI);
                f_r = f_r < 0.0 ? 0.0 : f_r;
                alpha_r = (1.0 + f_r) / (f_r + (MD->Ele[i].Rmin / Rmax));
                alpha_r = alpha_r > 10000.0 ? 10000.0 : alpha_r;
                eta_s = 1.0 - 0.0016 * (pow ((MD->Tref - 273.15 - T), 2));
                eta_s = eta_s < 0.0001 ? 0.0001 : eta_s;
                gamma_s = 1.0 / (1.0 + 0.00025 * (VP / RH - VP));
                gamma_s = (gamma_s < 0.01) ? 0.01 : gamma_s;
                r_s = ((MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s)) > Rmax) ? Rmax : (MD->Ele[i].Rmin * alpha_r / (beta_s * LAI * eta_s * gamma_s));
                P_c = (1.0 + Delta / Gamma) / (1.0 + r_s / r_a + Delta / Gamma);
                MD->EleET[i][1] = MD->pcCal.Et1 * MD->Ele[i].VegFrac * P_c * (1.0 - pow(((MD->EleIS[i] + MD->EleSnowCanopy[i] < 0.0) ? 0.0 : (MD->EleIS[i] + MD->EleSnowCanopy[i])) / (MD->EleISmax[i] + MD->EleISsnowmax[i]), MD->fx_canopy)) * ETp;
                MD->EleET[i][1] = MD->EleET[i][1] < 0.0 ? 0.0 : MD->EleET[i][1];
                AquiferDepth = MD->Ele[i].zmax - MD->Ele[i].zmin;
                MD->EleET[i][1] = ((MD->DummyY[i + 2 * MD->NumEle] < (AquiferDepth - MD->Ele[i].RzD)) && MD->DummyY[i + MD->NumEle] <= 0.0) ? 0.0 : MD->EleET[i][1];

                MD->EleTF[i] = MD->EleIS[i] <= 0.0 ? 0.0 : 5.65e-2 * MD->EleISmax[i] * exp (3.89 * (MD->EleIS[i] < 0.0 ? 0.0 : MD->EleIS[i]) / MD->EleISmax[i]);
            }
            else
            {
                MD->EleET[i][1] = 0.0;
                MD->EleET[i][0] = 0.0;
                MD->EleTF[i] = 0.0;
            }

            if (MD->EleTF[i] < 0.0)
                MD->EleTF[i] = 0.0;
            if (MD->EleTF[i] * stepsize > MD->EleIS[i])
                MD->EleTF[i] = MD->EleIS[i] / stepsize;
            if (MD->EleIS[i] >= MD->EleISmax[i])
            {
                if (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) >= MD->EleET[i][0] + MD->EleTF[i])
                {
                    ret = MD->EleTF[i] + (MD->EleIS[i] - MD->EleISmax[i])/stepsize + (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - (MD->EleET[i][0] + MD->EleTF[i]));
                    isval = MD->EleISmax[i];
                }
                else if ((((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) < MD->EleET[i][0] + MD->EleTF[i]) && (MD->EleIS[i] + stepsize * ((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy - MD->EleET[i][0] - MD->EleTF[i]) <= 0.0))
                {
                    MD->EleET[i][0] = (MD->EleET[i][0] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                    ret = (MD->EleTF[i] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                    isval = 0.0;
                }
                else
                {
                    isval = MD->EleIS[i] + stepsize * (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]);
                    ret = MD->EleTF[i];
                }
            }
            else if ((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) >= MD->EleISmax[i]))
            {
                isval = MD->EleISmax[i];
                ret = MD->EleTF[i] + (((MD->EleIS[i] + (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) - MD->EleISmax[i]))/stepsize;
            }
            else if ((MD->EleIS[i] < MD->EleISmax[i]) && ((MD->EleIS[i] + (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize) <= 0.0))
            {
                if ((MD->EleET[i][0] > 0.0) || (MD->EleTF[i] > 0.0))
                {
                    MD->EleET[i][0] = (MD->EleET[i][0] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                    ret = (MD->EleTF[i] / (MD->EleET[i][0] + MD->EleTF[i])) * (MD->EleIS[i] / stepsize + ((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy));
                }
                else
                {
                    MD->EleET[i][0] = 0.0;
                    ret = 0.0;
                }
                isval = 0.0;
            }
            else
            {
                isval = MD->EleIS[i] + (((1.0 - fracSnow) * MD->ElePrep[i] * MD->Ele[i].VegFrac + MeltRateCanopy) - MD->EleET[i][0] - MD->EleTF[i]) * stepsize;
                ret = MD->EleTF[i];
            }
            MD->EleNetPrep[i] = (1.0 - MD->Ele[i].VegFrac) * (1.0 - fracSnow) * MD->ElePrep[i] + ret + MeltRateGrnd;
            MD->EleTF[i] = ret;
            MD->EleIS[i] = isval;
        }

    free (metarr);
}
