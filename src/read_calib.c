#include "pihm.h"

void ReadCalib(const char *filename, calib_struct *cal)
{
    char            cmdstr[MAXSTRING];
    FILE           *global_calib;   /* Pointer to .calib file */
    int             lno = 0;

    global_calib = PIHMfopen(filename, "r");
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", &cal->ksath, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", &cal->ksatv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KINF", &cal->kinfv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATH", &cal->kmach, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KMACSATV", &cal->kmacv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DINF", &cal->dinf, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DROOT", &cal->rzd, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DMAC", &cal->dmac, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", &cal->porosity, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", &cal->alpha, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", &cal->beta, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACVF", &cal->areafv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "MACHF", &cal->areafh, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "VEGFRAC", &cal->vegfrac, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALBEDO", &cal->albedo, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH", &cal->rough, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "EC", &cal->ec, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ETT", &cal->ett, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "EDIR", &cal->edir, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ROUGH_RIV", &cal->rivrough, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KRIVH", &cal->rivksath, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KRIVV", &cal->rivksatv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BEDTHCK", &cal->rivbedthick, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_DPTH", &cal->rivdepth, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_WDTH", &cal->rivshpcoeff, 'd', filename, lno);

#if defined(_NOAH_)
    FindLine(global_calib, "LSM_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "DRIP", &cal->drip, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CMCMAX", &cal->cmcmax, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RS", &cal->rsmin, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CZIL", &cal->czil, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "FXEXP", &cal->fxexp, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "CFACTR", &cal->cfactr, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "RGL", &cal->rgl, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "HS", &cal->hs, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "REFSMC", &cal->smcref, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "WLTSMC", &cal->smcwlt, 'd', filename, lno);
#endif

#if defined(_FBR_)
    FindLine(global_calib, "DGW_CALIBRATION", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATH", &cal->geol_ksath, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "KSATV", &cal->geol_ksatv, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "POROSITY", &cal->geol_porosity, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ALPHA", &cal->geol_alpha, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "BETA", &cal->geol_beta, 'd', filename, lno);
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
    ReadKeyword(cmdstr, "rate", &cal->rate, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "ssa", &cal->ssa, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "gwinflux", &cal->gwinflux, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "prcpconc", &cal->prcpconc, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "initconc", &cal->initconc, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "Xsorption", &cal->Xsorption, 'd', filename, lno);
#endif

    /*
     * Scenarios
     */
    FindLine(global_calib, "SCENARIO", &lno, filename);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "PRCP", &cal->prcp, 'd', filename, lno);

    NextLine(global_calib, cmdstr, &lno);
    ReadKeyword(cmdstr, "SFCTMP", &cal->sfctmp, 'd', filename, lno);

    fclose(global_calib);
}
