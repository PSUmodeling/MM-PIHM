#include "pihm.h"

void WriteRtIc(const char *outputdir, const chemtbl_struct chemtbl[], const rttbl_struct *rttbl, elem_struct elem[])
{
    int             i;
    FILE           *fp;
    char            fn[MAXSTRING];

    sprintf(fn, "%s/restart/%s.rtic", outputdir, project);

    fp = pihm_fopen(fn, "wb");
    pihm_printf(VL_VERBOSE, "Writing RT initial conditions.\n");

    for (i = 0; i < nelem; i++)
    {
        int             j, k;

        for (k = 0; k < MAXSPS; k++)
        {
            if (k < rttbl->num_stc && chemtbl[k].itype == MINERAL)
            {
                elem[i].chms.tot_conc[k] /= (rttbl->rel_min == 0) ?
                    1000.0 / chemtbl[k].molar_vol / elem[i].soil.smcmax :
                    (1.0 - elem[i].soil.smcmax) * 1000.0 / chemtbl[k].molar_vol / elem[i].soil.smcmax;
#if defined(_DGW_)
                elem[i].chms_geol.tot_conc[k] /= (rttbl->rel_min == 0) ?
                    1000.0 / chemtbl[k].molar_vol / elem[i].geol.smcmax :
                    (1.0 - elem[i].geol.smcmax) * 1000.0 / chemtbl[k].molar_vol / elem[i].geol.smcmax;
#endif
            }
            else if (k < rttbl->num_stc && (chemtbl[k].itype == CATION_ECHG || chemtbl[k].itype == ADSORPTION))
            {
                elem[i].chms.tot_conc[k] /= (1.0 - elem[i].soil.smcmax) * 2650.0;
#if defined(_DGW_)
                elem[i].chms_geol.tot_conc[k] /= (1.0 - elem[i].geol.smcmax) * 2650.0;
#endif
            }

            elem[i].restart_output[SOIL_CHMVOL].tot_conc[k] = elem[i].chms.tot_conc[k];
            elem[i].restart_output[SOIL_CHMVOL].ssa[k] = elem[i].chms.ssa[k];

#if defined(_DGW_)
            elem[i].restart_output[GEOL_CHMVOL].tot_conc[k] = elem[i].chms_geol.tot_conc[k];
            elem[i].restart_output[GEOL_CHMVOL].ssa[k] = elem[i].chms_geol.ssa[k];
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

    fp = pihm_fopen(fn, "rb");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

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
