#include "pihm.h"

void ReadPara(const char *filename, ctrl_struct *ctrl)
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

    para_file = fopen(filename, "r");
    CheckFile(para_file, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    /* Start reading para_file */
    /* Read through parameter file to find parameters */
    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SIMULATION_MODE", &spinup_mode, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_MODE", &ctrl->init_type, 'i', filename, lno);
    ctrl->init_type = (ctrl->init_type > RELAX) ? RST_FILE : RELAX;

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "ASCII_OUTPUT", &ctrl->ascii, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "WATBAL_OUTPUT", &ctrl->waterbal, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", &ctrl->write_ic, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "UNSAT_MODE", &ctrl->unsat_mode, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SURF_MODE", &ctrl->surf_mode, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "RIV_MODE", &ctrl->riv_mode, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "START", &ctrl->starttime, 't', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "END", &ctrl->endtime, 't', filename, lno);

    /* In spinup mode, simulation time should be full years */
    start_time = PIHMTime(ctrl->starttime);
    end_time = PIHMTime(ctrl->endtime);

    if (end_time.t <= start_time.t)
    {
        PIHMprintf(VL_ERROR,
            "Error: simulation end time should be after start time.\n");
        PIHMprintf(VL_ERROR, "Please check your .para input file.\n");
        PIHMexit(EXIT_FAILURE);
    }
    else if (spinup_mode)
    {
        if (start_time.month != end_time.month ||
            start_time.day != end_time.day ||
            start_time.hour != end_time.hour ||
            start_time.minute != end_time.minute)
        {
            PIHMprintf(VL_ERROR,
                "Error: In BGC spinup mode, "
                "simulation period should be full years.\n");
            PIHMprintf(VL_ERROR, "Please check your .para input file.\n");
            PIHMexit(EXIT_FAILURE);
        }
    }

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MAX_SPINUP_YEAR", &ctrl->maxspinyears, 'i', filename,
        lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MODEL_STEPSIZE", &ctrl->stepsize, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "LSM_STEP", &ctrl->etstep, 'i', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "ABSTOL", &ctrl->abstol, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "RELTOL", &ctrl->reltol, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "INIT_SOLVER_STEP", &ctrl->initstep, 'd', filename,
        lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUM_NONCOV_FAIL", &ctrl->nncfn, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MAX_NONLIN_ITER", &ctrl->nnimax, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_NONLIN_ITER", &ctrl->nnimin, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "DECR_FACTOR", &ctrl->decr, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "INCR_FACTOR", &ctrl->incr, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MIN_MAXSTEP", &ctrl->stmin, 'd', filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SURF_CTRL] = ReadPrtCtrl(cmdstr, "SURF", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[UNSAT_CTRL] = ReadPrtCtrl(cmdstr, "UNSAT", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[GW_CTRL] = ReadPrtCtrl(cmdstr, "GW", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVSTG_CTRL] = ReadPrtCtrl(cmdstr, "RIVSTG", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVGW_CTRL] = ReadPrtCtrl(cmdstr, "RIVGW", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SNOW_CTRL] = ReadPrtCtrl(cmdstr, "SNOW", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[CMC_CTRL] = ReadPrtCtrl(cmdstr, "CMC", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[INFIL_CTRL] = ReadPrtCtrl(cmdstr, "INFIL", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RECHARGE_CTRL] = ReadPrtCtrl(cmdstr, "RECHARGE", filename,
        lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[EC_CTRL] = ReadPrtCtrl(cmdstr, "EC", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[ETT_CTRL] = ReadPrtCtrl(cmdstr, "ETT", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[EDIR_CTRL] = ReadPrtCtrl(cmdstr, "EDIR", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX0_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX0", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX1_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX1", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX2_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX2", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX3_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX3", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX4_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX4", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX5_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX5", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX6_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX6", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX7_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX7", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX8_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX8", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX9_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX9", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[RIVFLX10_CTRL] = ReadPrtCtrl(cmdstr, "RIVFLX10", filename,
        lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SUBFLX_CTRL] = ReadPrtCtrl(cmdstr, "SUBFLX", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[SURFFLX_CTRL] = ReadPrtCtrl(cmdstr, "SURFFLX", filename, lno);

    NextLine(para_file, cmdstr, &lno);
    ctrl->prtvrbl[IC_CTRL] = ReadPrtCtrl(cmdstr, "IC", filename, lno);

    fclose(para_file);

    if (ctrl->etstep < ctrl->stepsize || ctrl->etstep % ctrl->stepsize > 0)
    {
        PIHMprintf(VL_ERROR,
            "Error: Land surface model (ET) step size "
            "should be an integral multiple of model step size.\n");
        PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
        PIHMexit(EXIT_FAILURE);
    }
}
