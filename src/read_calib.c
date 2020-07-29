#include "pihm.h"

void ReadCalib(const char filename[], calib_struct *calib)
{
    char            cmdstr[MAXSTRING];
    FILE           *global_calib;
    int             lno = 0;

    global_calib = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", &calib->ksath, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", &calib->ksatv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KINF", &calib->kinfv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", &calib->kmach, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", &calib->kmacv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DINF", &calib->dinf, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DROOT", &calib->rzd, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", &calib->dmac, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", &calib->porosity, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", &calib->alpha, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", &calib->beta, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", &calib->areafv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", &calib->areafh, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "VEGFRAC", &calib->vegfrac, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALBEDO", &calib->albedo, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH", &calib->rough, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "EC", &calib->ec, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ETT", &calib->ett, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "EDIR", &calib->edir, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH_RIV", &calib->rivrough, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KRIVH", &calib->rivksath, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_DPTH", &calib->rivdepth, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_WDTH", &calib->rivshpcoeff, 'd', filename, lno);

#if defined(_NOAH_)
    FindLine(global_calib, "LSM_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DRIP", &calib->drip, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CMCMAX", &calib->cmcmax, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RS", &calib->rsmin, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CZIL", &calib->czil, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "FXEXP", &calib->fxexp, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CFACTR", &calib->cfactr, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RGL", &calib->rgl, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "HS", &calib->hs, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "REFSMC", &calib->smcref, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "WLTSMC", &calib->smcwlt, 'd', filename, lno);
#endif

#if defined(_FBR_)
    FindLine(global_calib, "DGW_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", &calib->geol_ksath, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", &calib->geol_ksatv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", &calib->geol_porosity, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", &calib->geol_alpha, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", &calib->geol_beta, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", &calib->geol_kmach, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", &calib->geol_kmacv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", &calib->geol_dmac, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", &calib->geol_areafv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", &calib->geol_areafh, 'd', filename, lno);
#endif

#if defined(_BGC_)
    FindLine(global_calib, "BGC_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MORTALITY", &cal->mortality, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "SLA", &cal->sla, 'd', filename, lno);
#endif

#if defined(_RT_)
    FindLine(global_calib, "RT_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "rate", &calib->rate, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ssa", &calib->ssa, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "prcpconc", &calib->prcpconc, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "initconc", &calib->initconc, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "Xsorption", &calib->Xsorption, 'd', filename, lno);
#endif

    /*
     * Scenarios
     */
    FindLine(global_calib, "SCENARIO", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRCP", &calib->prcp, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "SFCTMP", &calib->sfctmp, 'd', filename, lno);

    fclose(global_calib);
}
