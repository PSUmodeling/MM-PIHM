#include "pihm.h"

void ReadTecplot(const char *filename, ctrl_struct *ctrl)
{
    FILE           *tecplot_file;
    char            cmdstr[MAXSTRING];
    int             i;
    int             lno = 0;

    for (i = 0; i < MAXPRINT; i++)
    {
        ctrl->tpprtvrbl[i] = 0;
    }

    tecplot_file = fopen(filename, "r");
    CheckFile(tecplot_file, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    NextLine(tecplot_file, cmdstr, &lno);
    ctrl->tpprtvrbl[SURF_CTRL] = ReadPrtCtrl(cmdstr, "SURF", filename, lno);

    NextLine(tecplot_file, cmdstr, &lno);
    ctrl->tpprtvrbl[UNSAT_CTRL] = ReadPrtCtrl(cmdstr, "UNSAT", filename, lno);

    NextLine(tecplot_file, cmdstr, &lno);
    ctrl->tpprtvrbl[GW_CTRL] = ReadPrtCtrl(cmdstr, "GW", filename, lno);

    NextLine(tecplot_file, cmdstr, &lno);
    ctrl->tpprtvrbl[RIVSTG_CTRL] = ReadPrtCtrl(cmdstr, "RIVSTG", filename, lno);

    NextLine(tecplot_file, cmdstr, &lno);
    ctrl->tpprtvrbl[RIVGW_CTRL] = ReadPrtCtrl(cmdstr, "RIVGW", filename, lno);

    fclose(tecplot_file);
}
