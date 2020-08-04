#include "pihm.h"

void ReadCalib(const char filename[], calib_struct *calib)
{
    char            cmdstr[MAXSTRING];
    FILE           *global_calib;
    int             lno = 0;

    global_calib = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", 'd', filename, lno, &calib->ksath);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", 'd', filename, lno, &calib->ksatv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KINF", 'd', filename, lno, &calib->kinfv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", 'd', filename, lno, &calib->kmach);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", 'd', filename, lno, &calib->kmacv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DINF", 'd', filename, lno, &calib->dinf);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DROOT", 'd', filename, lno, &calib->rzd);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", 'd', filename, lno, &calib->dmac);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", 'd', filename, lno, &calib->porosity);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", 'd', filename, lno, &calib->alpha);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", 'd', filename, lno, &calib->beta);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", 'd', filename, lno, &calib->areafv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", 'd', filename, lno, &calib->areafh);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "VEGFRAC", 'd', filename, lno, &calib->vegfrac);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALBEDO", 'd', filename, lno, &calib->albedo);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH", 'd', filename, lno, &calib->rough);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH_RIV", 'd', filename, lno, &calib->rivrough);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KRIVH", 'd', filename, lno, &calib->rivksath);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_DPTH", 'd', filename, lno, &calib->rivdepth);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_WDTH", 'd', filename, lno, &calib->rivshpcoeff);

#if defined(_NOAH_)
    FindLine(global_calib, "LSM_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DRIP", 'd', filename, lno, &calib->drip);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CMCMAX", 'd', filename, lno, &calib->cmcmax);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RS", 'd', filename, lno, &calib->rsmin);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CZIL", 'd', filename, lno, &calib->czil);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "FXEXP", 'd', filename, lno, &calib->fxexp);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CFACTR", 'd', filename, lno, &calib->cfactr);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RGL", 'd', filename, lno, &calib->rgl);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "HS", 'd', filename, lno, &calib->hs);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "REFSMC", 'd', filename, lno, &calib->smcref);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "WLTSMC", 'd', filename, lno, &calib->smcwlt);
#endif

#if defined(_DGW_)
    FindLine(global_calib, "DGW_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", 'd', filename, lno, &calib->geol_ksath);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", 'd', filename, lno, &calib->geol_ksatv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", 'd', filename, lno, &calib->geol_porosity);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", 'd', filename, lno, &calib->geol_alpha);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", 'd', filename, lno, &calib->geol_beta);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", 'd', filename, lno, &calib->geol_kmach);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", 'd', filename, lno, &calib->geol_kmacv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", 'd', filename, lno, &calib->geol_dmac);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", 'd', filename, lno, &calib->geol_areafv);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", 'd', filename, lno, &calib->geol_areafh);
#endif

#if defined(_BGC_)
    FindLine(global_calib, "BGC_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MORTALITY", 'd', filename, lno, &calib->mortality);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "SLA", 'd', filename, lno, &calib->sla);
#endif

#if defined(_RT_)
    FindLine(global_calib, "RT_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "rate", 'd', filename, lno, &calib->rate);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ssa", 'd', filename, lno, &calib->ssa);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "Xsorption", 'd', filename, lno, &calib->Xsorption);
#endif

    /*
     * Scenarios
     */
    FindLine(global_calib, "SCENARIO", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRCP", 'd', filename, lno, &calib->prcp);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "SFCTMP", 'd', filename, lno, &calib->sfctmp);

    fclose(global_calib);
}
