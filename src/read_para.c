#include "pihm.h"

void ReadPara(const char filename[], ctrl_struct *ctrl)
{
    FILE           *para_file;    /* Pointer to .para file */
    char            cmdstr[MAXSTRING];
    int             i;
    int             lno = 0;
    pihm_t_struct   start_time, end_time;

    for (i = 0; i < MAXPRINT; i++)
    {
        ctrl->prtvrbl[i] = 0;
    }

    para_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    /* Start reading para_file */
    /* Read through parameter file to find parameters */
    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SIMULATION_MODE", 'i', filename, lno, &spinup_mode);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_MODE", 'i', filename, lno, &ctrl->init_type);
    ctrl->init_type = (ctrl->init_type > RELAX) ? RST_FILE : RELAX;

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "ASCII_OUTPUT", 'i', filename, lno, &ctrl->ascii);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "WATBAL_OUTPUT", 'i', filename, lno, &ctrl->waterbal);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", 'i', filename, lno, &ctrl->write_ic);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "UNSAT_MODE", 'i', filename, lno, &ctrl->unsat_mode);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SURF_MODE", 'i', filename, lno, &ctrl->surf_mode);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_MODE", 'i', filename, lno, &ctrl->riv_mode);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "START", 't', filename, lno, &ctrl->starttime);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "END", 't', filename, lno, &ctrl->endtime);

    /* In spinup mode, simulation time should be full years */
    start_time = PIHMTime(ctrl->starttime);
    end_time = PIHMTime(ctrl->endtime);

    if (end_time.t <= start_time.t)
    {
        pihm_printf(VL_ERROR,
            "Error: simulation end time should be after start time.\n"
            "Please check your .para input file.\n");
        pihm_exit(EXIT_FAILURE);
    }
    else if (spinup_mode)
    {
        if (start_time.month != end_time.month ||
            start_time.day != end_time.day ||
            start_time.hour != end_time.hour ||
            start_time.minute != end_time.minute)
        {
            pihm_printf(VL_ERROR,
                "Error: In BGC spinup mode, simulation period should be full "
                "years.\nPlease check your .para input file.\n");
            pihm_exit(EXIT_FAILURE);
        }
    }

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MAX_SPINUP_YEAR", 'i', filename, lno, &ctrl->maxspinyears);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MODEL_STEPSIZE", 'i', filename, lno, &ctrl->stepsize);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "LSM_STEP", 'i', filename, lno, &ctrl->etstep);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "ABSTOL", 'd', filename, lno, &ctrl->abstol);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "RELTOL", 'd', filename, lno, &ctrl->reltol);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_SOLVER_STEP", 'd', filename, lno, &ctrl->initstep);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUM_NONCOV_FAIL", 'd', filename, lno, &ctrl->nncfn);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MAX_NONLIN_ITER", 'd', filename, lno, &ctrl->nnimax);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_NONLIN_ITER", 'd', filename, lno, &ctrl->nnimin);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "DECR_FACTOR", 'd', filename, lno, &ctrl->decr);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "INCR_FACTOR", 'd', filename, lno, &ctrl->incr);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_MAXSTEP", 'd', filename, lno, &ctrl->stmin);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SURF_CTRL] = ReadPrintCtrl(cmdstr, "SURF", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[UNSAT_CTRL] = ReadPrintCtrl(cmdstr, "UNSAT", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[GW_CTRL] = ReadPrintCtrl(cmdstr, "GW", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[STAGE_CTRL] = ReadPrintCtrl(cmdstr, "RIVSTG", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVGW_CTRL] = ReadPrintCtrl(cmdstr, "RIVGW", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SNOW_CTRL] = ReadPrintCtrl(cmdstr, "SNOW", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[CMC_CTRL] = ReadPrintCtrl(cmdstr, "CMC", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[INFIL_CTRL] = ReadPrintCtrl(cmdstr, "INFIL", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RECHARGE_CTRL] = ReadPrintCtrl(cmdstr, "RECHARGE", filename,
        lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[EC_CTRL] = ReadPrintCtrl(cmdstr, "EC", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[ETT_CTRL] = ReadPrintCtrl(cmdstr, "ETT", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[EDIR_CTRL] = ReadPrintCtrl(cmdstr, "EDIR", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX0_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX0", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX1_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX1", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX2_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX2", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX3_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX3", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX4_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX4", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX5_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX5", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SUBFLX_CTRL] = ReadPrintCtrl(cmdstr, "SUBFLX", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SURFFLX_CTRL] = ReadPrintCtrl(cmdstr, "SURFFLX", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[IC_CTRL] = ReadPrintCtrl(cmdstr, "IC", filename, lno);

    fclose(para_file);

#if defined(_CYCLES_)
    if (ctrl->etstep != DAYINSEC)
    {
        pihm_printf(VL_ERROR, "Warning: When coupled to Cycles, daily land "
            "surface step is required. Thus land\nsurface model step is "
            "changed to daily (86400 seconds).\n\n");
        ctrl->etstep = DAYINSEC;
    }
#endif

    if (ctrl->etstep < ctrl->stepsize || ctrl->etstep % ctrl->stepsize > 0)
    {
        pihm_printf(VL_ERROR,
            "Error: Land surface model (ET) step size "
            "should be an integral multiple of model step size.\n");
        pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
        pihm_exit(EXIT_FAILURE);
    }
}
