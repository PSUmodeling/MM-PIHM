#include "pihm.h"

void WriteRtIc(const char *outputdir, const chemtbl_struct chemtbl[],
    const rttbl_struct *rttbl, elem_struct elem[], river_struct river[])
{
    int             i;
    FILE           *fp;
    char            restart_fn[MAXSTRING];

    sprintf(restart_fn, "%s/restart/%s.rtic", outputdir, project);

    fp = fopen(restart_fn, "wb");
    CheckFile(fp, restart_fn);
    PIHMprintf(VL_VERBOSE, "Writing RT initial conditions.\n");

    for (i = 0; i < nelem; i++)
    {
        int             j, k;

        for (k = 0; k < MAXSPS; k++)
        {
            if (k < rttbl->NumStc && chemtbl[k].itype == MINERAL)
            {
                elem[i].chms_unsat.t_conc[k] /= (rttbl->RelMin == 0) ?
                    1000.0 / chemtbl[k].MolarVolume / elem[i].soil.smcmax :
                    (1.0 - elem[i].soil.smcmax) * 1000.0 /
                    chemtbl[k].MolarVolume / elem[i].soil.smcmax;
                elem[i].chms_gw.t_conc[k] /= (rttbl->RelMin == 0) ?
                    1000.0 / chemtbl[k].MolarVolume / elem[i].soil.smcmax :
                    (1.0 - elem[i].soil.smcmax) * 1000.0 /
                    chemtbl[k].MolarVolume / elem[i].soil.smcmax;
#if defined(_FBR_)
                elem[i].chms_fbrunsat.t_conc[k] /= (rttbl->RelMin == 0) ?
                    1000.0 / chemtbl[k].MolarVolume / elem[i].geol.smcmax :
                    (1.0 - elem[i].geol.smcmax) * 1000.0 /
                    chemtbl[k].MolarVolume / elem[i].geol.smcmax;
                elem[i].chms_fbrgw.t_conc[k] /= (rttbl->RelMin == 0) ?
                    1000.0 / chemtbl[k].MolarVolume / elem[i].geol.smcmax :
                    (1.0 - elem[i].geol.smcmax) * 1000.0 /
                    chemtbl[k].MolarVolume / elem[i].geol.smcmax;
#endif
            }
            else if (k < rttbl->NumStc && (chemtbl[k].itype == CATION_ECHG ||
                chemtbl[k].itype == ADSORPTION))
            {
                elem[i].chms_unsat.t_conc[k] /=
                    (1.0 - elem[i].soil.smcmax) * 2650.0;
                elem[i].chms_gw.t_conc[k] /=
                    (1.0 - elem[i].soil.smcmax) * 2650.0;
#if defined(_FBR_)
                elem[i].chms_fbrunsat.t_conc[k] /=
                    (1.0 - elem[i].geol.smcmax) * 2650.0;
                elem[i].chms_fbrgw.t_conc[k] /=
                    (1.0 - elem[i].geol.smcmax) * 2650.0;
#endif
            }

            elem[i].restart_output[UNSAT_CHMVOL].t_conc[k] =
                elem[i].chms_unsat.t_conc[k];
            elem[i].restart_output[UNSAT_CHMVOL].ssa[k] =
                elem[i].chms_unsat.ssa[k];

            elem[i].restart_output[GW_CHMVOL].t_conc[k] =
                elem[i].chms_gw.t_conc[k];
            elem[i].restart_output[GW_CHMVOL].ssa[k] =
                elem[i].chms_gw.ssa[k];
#if defined(_FBR_)
            elem[i].restart_output[FBRUNSAT_CHMVOL].t_conc[k] =
                elem[i].chms_fbrunsat.t_conc[k];
            elem[i].restart_output[FBRUNSAT_CHMVOL].ssa[k] =
                elem[i].chms_fbrunsat.ssa[k];

            elem[i].restart_output[FBRGW_CHMVOL].t_conc[k] =
                elem[i].chms_fbrgw.t_conc[k];
            elem[i].restart_output[FBRGW_CHMVOL].ssa[k] =
                elem[i].chms_fbrgw.ssa[k];
#endif
        }

        for (j = 0; j < NCHMVOL; j++)
        {
            fwrite(&(elem[i].restart_output[j]), sizeof(rtic_struct), 1, fp);
        }
    }

    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < MAXSPS; k++)
        {
            river[i].restart_output.tconc_stream[k] =
                river[i].chms_stream.t_conc[k];
            river[i].restart_output.tconc_rivbed[k] =
                river[i].chms_rivbed.t_conc[k];
        }

        fwrite(&(river[i].restart_output), sizeof(river_rtic_struct), 1, fp);
    }

    fclose(fp);
}

void ReadRtIc(const char *fn, elem_struct elem[], river_struct river[])
{
    FILE           *fp;
    int             i;

    fp = fopen(fn, "rb");
    CheckFile(fp, fn);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", fn);

    for (i = 0; i < nelem; i++)
    {
        int             j;

        for (j = 0; j < NCHMVOL; j++)
        {
            fread(&elem[i].restart_input[j], sizeof(rtic_struct), 1, fp);
        }
    }

    for (i = 0; i < nriver; i++)
    {
        fread(&river[i].restart_input, sizeof(river_rtic_struct), 1, fp);
    }

    fclose(fp);
}
