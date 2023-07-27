#include "pihm.h"

void ReadCalib(const char fn[], calib_struct *calib)
{
    char            cmdstr[MAXSTRING];
    FILE           *fp;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", 'd', fn, lno, &calib->ksath);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", 'd', fn, lno, &calib->ksatv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KINF", 'd', fn, lno, &calib->kinfv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", 'd', fn, lno, &calib->kmach);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", 'd', fn, lno, &calib->kmacv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DROOT", 'd', fn, lno, &calib->rzd);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", 'd', fn, lno, &calib->dmac);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", 'd', fn, lno, &calib->porosity);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", 'd', fn, lno, &calib->alpha);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", 'd', fn, lno, &calib->beta);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", 'd', fn, lno, &calib->areafv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", 'd', fn, lno, &calib->areafh);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "VEGFRAC", 'd', fn, lno, &calib->vegfrac);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALBEDO", 'd', fn, lno, &calib->albedo);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH", 'd', fn, lno, &calib->rough);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH_RIV", 'd', fn, lno, &calib->rivrough);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KRIVH", 'd', fn, lno, &calib->rivksath);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_DPTH", 'd', fn, lno, &calib->rivdepth);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_WDTH", 'd', fn, lno, &calib->rivshpcoeff);

#if defined(_NOAH_)
    FindLine(fp, "LSM_CALIBRATION", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DRIP", 'd', fn, lno, &calib->drip);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CMCMAX", 'd', fn, lno, &calib->cmcmax);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RS", 'd', fn, lno, &calib->rsmin);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CZIL", 'd', fn, lno, &calib->czil);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "FXEXP", 'd', fn, lno, &calib->fxexp);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CFACTR", 'd', fn, lno, &calib->cfactr);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RGL", 'd', fn, lno, &calib->rgl);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "HS", 'd', fn, lno, &calib->hs);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "REFSMC", 'd', fn, lno, &calib->smcref);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WLTSMC", 'd', fn, lno, &calib->smcwlt);
#endif

#if defined(_CYCLES_)
    FindLine(fp, "AG_CALIBRATION", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "FERTILIZATION", 'd', fn, lno, &calib->fert);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SOC_DECOMP_RATE", 'd', fn, lno, &calib->soc_decomp_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RESIDUE_DECOMP_RATE", 'd', fn, lno, &calib->residue_decomp_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROOT_DECOMP_RATE", 'd', fn, lno, &calib->root_decomp_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RHIZO_DECOMP_RATE", 'd', fn, lno, &calib->rhizo_decomp_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MANURE_DECOMP_RATE", 'd', fn, lno, &calib->manure_decomp_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MICROB_DECOMP_RATE", 'd', fn, lno, &calib->microb_decomp_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SOC_HUMIF_POWER", 'd', fn, lno, &calib->soc_humif_power);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NITRIF_RATE", 'd', fn, lno, &calib->nitrif_const);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "POT_DENITRIF_RATE", 'd', fn, lno, &calib->pot_denitrif);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DENITRIF_HALF_RATE", 'd', fn, lno, &calib->denitrif_half_rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DECOMP_HALF_RESP", 'd', fn, lno, &calib->decomp_half_resp);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DECOMP_RESP_POWER", 'd', fn, lno, &calib->decomp_resp_power);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KD_NO3", 'd', fn, lno, &calib->kd_no3);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KD_NH4", 'd', fn, lno, &calib->kd_nh4);
#endif

#if defined(_DGW_)
    FindLine(fp, "DGW_CALIBRATION", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", 'd', fn, lno, &calib->geol_ksath);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", 'd', fn, lno, &calib->geol_ksatv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", 'd', fn, lno, &calib->geol_porosity);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", 'd', fn, lno, &calib->geol_alpha);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", 'd', fn, lno, &calib->geol_beta);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", 'd', fn, lno, &calib->geol_kmach);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", 'd', fn, lno, &calib->geol_kmacv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", 'd', fn, lno, &calib->geol_dmac);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", 'd', fn, lno, &calib->geol_areafv);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", 'd', fn, lno, &calib->geol_areafh);
#endif

#if defined(_BGC_)
    FindLine(fp, "BGC_CALIBRATION", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MORTALITY", 'd', fn, lno, &calib->mortality);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SLA", 'd', fn, lno, &calib->sla);
#endif

#if defined(_RT_)
    FindLine(fp, "RT_CALIBRATION", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "rate", 'd', fn, lno, &calib->rate);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ssa", 'd', fn, lno, &calib->ssa);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "Xsorption", 'd', fn, lno, &calib->Xsorption);
#endif

    // Scenarios
    FindLine(fp, "SCENARIO", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRCP", 'd', fn, lno, &calib->prcp);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SFCTMP", 'd', fn, lno, &calib->sfctmp);

    fclose(fp);
}
