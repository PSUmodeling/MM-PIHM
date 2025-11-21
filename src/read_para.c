#include "pihm.h"

void ReadPara(const char fn[], ctrl_struct *ctrl)
{
    FILE           *fp;                     // Pointer to .para file
    char            cmdstr[MAXSTRING];
    int             i;
    int             lno = 0;
    pihm_t_struct   start_time, end_time;

    for (i = 0; i < MAXPRINT; i++)
    {
        ctrl->prtvrbl[i] = 0;
    }

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    // Start reading para_file
    // Read through parameter file to find parameters
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SIMULATION_MODE", 'i', fn, lno, &spinup_mode);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_MODE", 'i', fn, lno, &ctrl->init_type);
    ctrl->init_type = (ctrl->init_type > RELAX) ? RST_FILE : RELAX;

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ASCII_OUTPUT", 'i', fn, lno, &ctrl->ascii);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WATBAL_OUTPUT", 'i', fn, lno, &ctrl->waterbal);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", 'i', fn, lno, &ctrl->write_ic);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "START", 't', fn, lno, &ctrl->starttime);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "END", 't', fn, lno, &ctrl->endtime);

    // In spinup mode, simulation time should be full years
    start_time = PIHMTime(ctrl->starttime);
    end_time = PIHMTime(ctrl->endtime);

    if (end_time.t <= start_time.t)
    {
        pihm_printf(VL_ERROR, "Error: simulation end time should be after start time.\n"
            "Please check your .para input file.\n");
        pihm_exit(EXIT_FAILURE);
    }
    else if (spinup_mode)
    {
        if (start_time.month != end_time.month || start_time.day != end_time.day || start_time.hour != end_time.hour || start_time.minute != end_time.minute)
        {
            pihm_printf(VL_ERROR, "Error: In spinup mode, simulation period should be full years. Please check your\n.para input file.\n");
            pihm_exit(EXIT_FAILURE);
        }
    }

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MAX_SPINUP_YEAR", 'i', fn, lno, &ctrl->maxspinyears);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MODEL_STEPSIZE", 'i', fn, lno, &ctrl->stepsize);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "LSM_STEP", 'i', fn, lno, &ctrl->etstep);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "ABSTOL", 'd', fn, lno, &ctrl->abstol);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RELTOL", 'd', fn, lno, &ctrl->reltol);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_SOLVER_STEP", 'd', fn, lno, &ctrl->initstep);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUM_NONCOV_FAIL", 'd', fn, lno, &ctrl->nncfn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MAX_NONLIN_ITER", 'd', fn, lno, &ctrl->nnimax);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_NONLIN_ITER", 'd', fn, lno, &ctrl->nnimin);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "DECR_FACTOR", 'd', fn, lno, &ctrl->decr);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INCR_FACTOR", 'd', fn, lno, &ctrl->incr);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_MAXSTEP", 'd', fn, lno, &ctrl->stmin);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SURF_CTRL] = ReadPrintCtrl(cmdstr, "SURF", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[UNSAT_CTRL] = ReadPrintCtrl(cmdstr, "UNSAT", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[GW_CTRL] = ReadPrintCtrl(cmdstr, "GW", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[STAGE_CTRL] = ReadPrintCtrl(cmdstr, "RIVSTG", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SNOW_CTRL] = ReadPrintCtrl(cmdstr, "SNOW", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[CMC_CTRL] = ReadPrintCtrl(cmdstr, "CMC", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[INFIL_CTRL] = ReadPrintCtrl(cmdstr, "INFIL", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RECHARGE_CTRL] = ReadPrintCtrl(cmdstr, "RECHARGE", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[EC_CTRL] = ReadPrintCtrl(cmdstr, "EC", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[ETT_CTRL] = ReadPrintCtrl(cmdstr, "ETT", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[EDIR_CTRL] = ReadPrintCtrl(cmdstr, "EDIR", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX0_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX0", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX1_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX1", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX2_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX2", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX3_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX3", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX4_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX4", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX5_CTRL] = ReadPrintCtrl(cmdstr, "RIVFLX5", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SUBFLX_CTRL] = ReadPrintCtrl(cmdstr, "SUBFLX", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SURFFLX_CTRL] = ReadPrintCtrl(cmdstr, "SURFFLX", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[IC_CTRL] = ReadPrintCtrl(cmdstr, "IC", fn, lno);

    fclose(fp);

#if defined(_CYCLES_)
    if (ctrl->etstep != DAYINSEC)
    {
        pihm_printf(VL_ERROR, "Warning: When coupled to Cycles, daily land surface step is required. Thus land\n"
            "surface model step is changed to daily (86400 seconds).\n\n");
        ctrl->etstep = DAYINSEC;
    }
#endif

    if (ctrl->etstep < ctrl->stepsize || ctrl->etstep % ctrl->stepsize > 0)
    {
        pihm_printf(VL_ERROR, "Error: Land surface model (ET) step size should be an integral multiple of model step size.\n");
        pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
        pihm_exit(EXIT_FAILURE);
    }
}
