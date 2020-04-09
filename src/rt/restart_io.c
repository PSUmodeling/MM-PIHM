#include "pihm.h"

void WriteRtIc(const char *outputdir, const chemtbl_struct chemtbl[],
    const rttbl_struct *rttbl, elem_struct elem[])
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
                elem[i].chms.t_conc[k] /= (rttbl->RelMin == 0) ?
                    1000.0 / chemtbl[k].MolarVolume / elem[i].soil.smcmax :
                    (1.0 - elem[i].soil.smcmax) * 1000.0 /
                    chemtbl[k].MolarVolume / elem[i].soil.smcmax;
#if defined(_FBR_)
                elem[i].chms_geol.t_conc[k] /= (rttbl->RelMin == 0) ?
                    1000.0 / chemtbl[k].MolarVolume / elem[i].geol.smcmax :
                    (1.0 - elem[i].geol.smcmax) * 1000.0 /
                    chemtbl[k].MolarVolume / elem[i].geol.smcmax;
#endif
            }
            else if (k < rttbl->NumStc && (chemtbl[k].itype == CATION_ECHG ||
                chemtbl[k].itype == ADSORPTION))
            {
                elem[i].chms.t_conc[k] /=
                    (1.0 - elem[i].soil.smcmax) * 2650.0;
#if defined(_FBR_)
                elem[i].chms_geol.t_conc[k] /=
                    (1.0 - elem[i].geol.smcmax) * 2650.0;
#endif
            }

            elem[i].restart_output[SOIL_CHMVOL].t_conc[k] =
                elem[i].chms.t_conc[k];
            elem[i].restart_output[SOIL_CHMVOL].ssa[k] =
                elem[i].chms.ssa[k];

#if defined(_FBR_)
            elem[i].restart_output[GEOL_CHMVOL].t_conc[k] =
                elem[i].chms_geol.t_conc[k];
            elem[i].restart_output[GEOL_CHMVOL].ssa[k] =
                elem[i].chms_geol.ssa[k];
#endif
        }

        for (j = 0; j < NCHMVOL; j++)
        {
            fwrite(&(elem[i].restart_output[j]), sizeof(rtic_struct), 1, fp);
        }
    }

    fclose(fp);
}

void ReadRtIc(const char *fn, elem_struct elem[])
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

    fclose(fp);
}
