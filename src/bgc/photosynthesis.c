#include "pihm.h"

void TotalPhotosynthesis(const epconst_struct *epc, epvar_struct *epv,
    const pstate_struct *ps, cflux_struct *cf, psn_struct *psn_sun,
    psn_struct *psn_shade, const daily_struct *daily)
{
    /*
     * This function is a wrapper and replacement for the photosynthesis code
     * which used to be in the central bgc.c code. At Mott Jolly's request, all
     * of the science code is being moved into functions.
     */
    double          tday;

    tday = daily->tday - TFREEZ;

    /* SUNLIT canopy fraction photosynthesis */
    /* Set the input variables */
    psn_sun->c3 = epc->c3_flag;
    psn_sun->co2 = ps->co2;
    psn_sun->pa = daily->avg_sfcprs;
    psn_sun->t = tday;
    psn_sun->lnc = 1.0 / (epv->sun_proj_sla * epc->leaf_cn);
    psn_sun->flnr = epc->flnr;
    psn_sun->ppfd = ps->ppfd_per_plaisun;
    /* Convert conductance from m/s --> umol/m2/s/Pa, and correct for CO2 vs.
     * water vapor */
    psn_sun->g = epv->gl_t_wv_sun * 1.0e6 / (1.6 * 8.3143 * (tday + TFREEZ));
    psn_sun->dlmr = epv->dlmr_area_sun;

    /* Calculate photosynthesis for sunlit leaves */
    Photosynthesis(psn_sun);

    epv->assim_sun = psn_sun->A;

    /* For the final flux assignment, the assimilation output needs to have the
     * maintenance respiration rate added, this sum multiplied by the projected
     * leaf area in the relevant canopy fraction, and this total converted from
     * umol/m2/s -> kgC/m2/d */
    cf->psnsun_to_cpool = (epv->assim_sun + epv->dlmr_area_sun) * ps->plaisun *
        epv->dayl * 12.011e-9;

    /* SHADED canopy fraction photosynthesis */
    psn_shade->c3 = epc->c3_flag;
    psn_shade->co2 = ps->co2;
    psn_shade->pa = daily->avg_sfcprs;
    psn_shade->t = tday;
    psn_shade->lnc = 1.0 / (epv->shade_proj_sla * epc->leaf_cn);
    psn_shade->flnr = epc->flnr;
    psn_shade->ppfd = ps->ppfd_per_plaishade;
    /* Convert conductance from m/s --> umol/m2/s/Pa, and correct for CO2 vs.
     * water vapor */
    psn_shade->g =
        epv->gl_t_wv_shade * 1.0e6 / (1.6 * 8.3143 * (tday + TFREEZ));
    psn_shade->dlmr = epv->dlmr_area_shade;
    /* Calculate photosynthesis for shaded leaves */
    Photosynthesis(psn_shade);

    epv->assim_shade = psn_shade->A;

    /* For the final flux assignment, the assimilation output needs to have the
     * maintenance respiration rate added, this sum multiplied by the projected
     * leaf area in the relevant canopy fraction, and this total converted from
     * umol/m2/s -> kgC/m2/d */
    cf->psnshade_to_cpool = (epv->assim_shade + epv->dlmr_area_shade) *
        ps->plaishade * epv->dayl * 12.011e-9;
}

void Photosynthesis(psn_struct *psn)
{
    /*
     * The following variables are assumed to be defined in the psn struct
     * at the time of the function call:
     * c3         (flag) set to 1 for C3 model, 0 for C4 model
     * pa         (Pa) atmospheric pressure
     * co2        (ppm) atmospheric [CO2]
     * t          (deg C) air temperature
     * lnc        (kg Nleaf/m2) leaf N concentration, per unit projected LAI
     * flnr       (kg NRub/kg Nleaf) fraction of leaf N in Rubisco
     * ppfd       (umol photons/m2/s) PAR flux density, per unit projected LAI
     * g          (umol CO2/m2/s/Pa) leaf-scale conductance to CO2, proj area
     *              basis
     * dlmr       (umol CO2/m2/s) day leaf maint resp, on projected leaf area
     *              basis
     *
     * The following variables in psn struct are defined upon function return:
     * Ci         (Pa) intercellular [CO2]
     * Ca         (Pa) atmospheric [CO2]
     * O2         (Pa) atmospheric [O2]
     * gamma      (Pa) CO2 compensation point, in the absence of maint resp.
     * Kc         (Pa) MM constant for carboxylation
     * Ko         (Pa) MM constant for oxygenation
     * Vmax       (umol CO2/m2/s) max rate of carboxylation
     * Jmax       (umol electrons/m2/s) max rate electron transport
     * J          (umol RuBP/m2/s) rate of RuBP regeneration
     * Av         (umol CO2/m2/s) carboxylation limited assimilation
     * Aj         (umol CO2/m2/s) RuBP regen limited assimilation
     * A          (umol CO2/m2/s) final assimilation rate
     */

    /* The weight proportion of Rubisco to its nitrogen content, fnr, is
     * calculated from the relative proportions of the basic amino acids that
     * make up the enzyme, as listed in the Handbook of Biochemistry, Proteins,
     * Vol III, p. 510, which references: Kuehn and McFadden, Biochemistry,
     * 8:2403, 1969 */
    static double   fnr = 7.16;    /* kg Rub/kg NRub */

    /* The following enzyme kinetic constants are from:
     * Woodrow, I.E., and J.A. Berry, 1980. Enzymatic regulation of
     * photosynthetic CO2 fixation in C3 plants. Ann. Rev. Plant Physiol.
     * Plant Mol. Biol., 39:533-594.
     * Note that these values are given in the units used in the paper, and
     * that they are converted to units appropriate to the rest of this function
     * before they are used. */
    /* I've changed the values for Kc and Ko from the Woodrow and Berry
     * reference, and am now using the values from De Pury and Farquhar, 1997.
     * Simple scaling of photosynthesis from leaves to canopies without the
     * errors of big-leaf models. Plant, Cell and Env. 20: 537-557.
     * All other parameters, including the q10's for Kc and Ko are the same as
     * in Woodrow and Berry. */
    static double   Kc25 = 404.0;    /* (ubar) MM const carboxylase, 25 deg C */
    static double   q10Kc = 2.1;     /* (DIM) Q_10 for Kc */
    static double   Ko25 = 248.0;    /* (mbar) MM const oxygenase, 25 deg C */
    static double   q10Ko = 1.2;     /* (DIM) Q_10 for Ko */
    static double   act25 = 3.6;     /* (umol/mgRubisco/min) Rubisco activity */
    static double   q10act = 2.4;    /* (DIM) Q_10 for Rubisco activity */
    static double   pabs = 0.85;     /* (DIM) fPAR effectively absorbed by PSII
                                      */

    /* local variables */
    double          t;       /* (deg C) temperature */
    double          tk;      /* (K) absolute temperature */
    double          Kc;      /* (Pa) MM constant for carboxylase reaction */
    double          Ko;      /* (Pa) MM constant for oxygenase reaction */
    double          act;     /* (umol CO2/kgRubisco/s) Rubisco activity */
    double          Jmax;    /* (umol/m2/s) max rate electron transport */
    double          ppe;     /* (mol/mol) photons absorbed by PSII per e-
                              * transported */
    double          Vmax, J, gamma, Ca, Rd, O2, g;
    double          a, b, c, det;
    double          Av, Aj, A;

    /* Begin by assigning local variables */
    g = psn->g;
    t = psn->t;
    tk = t + 273.15;
    Rd = psn->dlmr;

    /* Convert atmospheric CO2 from ppm --> Pa */
    Ca = psn->co2 * psn->pa / 1e6;

    /* Set parameters for C3 vs C4 model */
    if (psn->c3)
    {
        ppe = 2.6;
    }
    else    /* C4 */
    {
        ppe = 3.5;
        Ca *= 10.0;
    }
    psn->Ca = Ca;

    /* Calculate atmospheric O2 in Pa, assumes 21% O2 by volume */
    psn->O2 = O2 = 0.21 * psn->pa;

    /* Correct kinetic constants for temperature, and do unit conversions */
    Ko = Ko25 * pow(q10Ko, (t - 25.0) / 10.0);
    psn->Ko = Ko = Ko * 100.0;    /* mbar --> Pa */
    if (t > 15.0)
    {
        Kc = Kc25 * pow(q10Kc, (t - 25.0) / 10.0);
        act = act25 * pow(q10act, (t - 25.0) / 10.0);
    }
    else
    {
        Kc = Kc25 * pow(1.8 * q10Kc, (t - 15.0) / 10.0) / q10Kc;
        act = act25 * pow(1.8 * q10act, (t - 15.0) / 10.0) / q10act;
    }
    psn->Kc = Kc = Kc * 0.10;    /* ubar --> Pa */
    act = act * 1e6 / 60.0;      /* umol/mg/min --> umol/kg/s */

    /* Calculate gamma (Pa), assumes Vomax/Vcmax = 0.21 */
    psn->gamma = gamma = 0.5 * 0.21 * Kc * psn->O2 / Ko;

    /* Calculate Vmax from leaf nitrogen data and Rubisco activity */

    /* kg Nleaf   kg NRub    kg Rub      umol            umol
     * -------- X -------  X ------- X ---------   =   --------
     *  m2         kg Nleaf   kg NRub   kg RUB * s       m2 * s
     *
     * (lnc)  X  (flnr)  X  (fnr)  X   (act)     =    (Vmax) */
    psn->Vmax = Vmax = psn->lnc * psn->flnr * fnr * act;

    /* Calculate Jmax = f(Vmax), reference:
     * Wullschleger, S.D., 1993.  Biochemical limitations to carbon assimilation
     * in C3 plants - A retrospective analysis of the A/Ci curves from 109
     * species. Journal of Experimental Botany, 44:907-920. */
    psn->Jmax = Jmax = 2.1 * Vmax;

    /* Calculate J = f(Jmax, ppfd), reference:
     * de Pury and Farquhar 1997 Plant Cell and Env. */
    a = 0.7;
    b = -Jmax - (psn->ppfd * pabs / ppe);
    c = Jmax * psn->ppfd * pabs / ppe;
    psn->J = J = (-b - sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

    /* Solve for Av and Aj using the quadratic equation, substitution for Ci
     * from A = g(Ca-Ci) into the equations from Farquhar and von Caemmerer:
     *
     *        Vmax (Ci - gamma)
     * Av =  -------------------   -   Rd
     *        Ci + Kc (1 + O2/Ko)
     *
     *          J (Ci - gamma)
     * Aj  =  -------------------  -   Rd
     *        4.5 Ci + 10.5 gamma
     */

    /* Quadratic solution for Av */
    a = -1.0 / g;
    b = Ca + (Vmax - Rd) / g + Kc * (1.0 + O2 / Ko);
    c = Vmax * (gamma - Ca) + Rd * (Ca + Kc * (1.0 + O2 / Ko));

    if ((det = b * b - 4.0 * a * c) < 0.0)
    {
        PIHMprintf(VL_ERROR, "ERROR: negative root error in psn routine\n");
        PIHMexit(EXIT_FAILURE);
    }

    psn->Av = Av = (-b + sqrt(det)) / (2.0 * a);

    /* Quadratic solution for Aj */
    a = -4.5 / g;
    b = 4.5 * Ca + 10.5 * gamma + J / g - 4.5 * Rd / g;
    c = J * (gamma - Ca) + Rd * (4.5 * Ca + 10.5 * gamma);

    if ((det = b * b - 4.0 * a * c) < 0.0)
    {
        PIHMprintf(VL_ERROR, "ERROR: negative root error in psn routine\n");
        PIHMexit(EXIT_FAILURE);
    }

    psn->Aj = Aj = (-b + sqrt(det)) / (2.0 * a);

    /* Estimate A as the minimum of (Av,Aj) */
    if (Av < Aj)
    {
        A = Av;
    }
    else
    {
        A = Aj;
    }
    psn->A = A;
    psn->Ci = Ca - (A / g);
}
