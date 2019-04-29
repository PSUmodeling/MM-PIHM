#include "pihm.h"

void WriteRtIc(const char *outputdir, elem_struct elem[], river_struct river[])
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
        int             k;

        for (k = 0; k < MAXSPS; k++)
        {
            elem[i].restart_output.tconc_unsat[k] =
                elem[i].chms_unsat.t_conc[k];
            elem[i].restart_output.ssa_unsat[k] =
                elem[i].chms_unsat.ssa[k];

            elem[i].restart_output.tconc_gw[k] =
                elem[i].chms_gw.t_conc[k];
            elem[i].restart_output.ssa_gw[k] =
                elem[i].chms_gw.ssa[k];
#if defined(_FBR_)
            elem[i].restart_output.tconc_fbrunsat[k] =
                elem[i].chms_fbrunsat.t_fbrconc[k];
            elem[i].restart_output.ssa_fbrunsat[k] =
                elem[i].chms_fbrunsat.ssa[k];

            elem[i].restart_output.tconc_fbrgw[k] =
                elem[i].chms_fbrgw.t_fbrconc[k];
            elem[i].restart_output.ssa_fbrgw[k] =
                elem[i].chms_fbrgw.ssa[k];
#endif
        }

        fwrite(&(elem[i].restart_output), sizeof(rtic_struct), 1, fp);
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
        fread(&elem[i].restart_input, sizeof(rtic_struct), 1, fp);
    }

    for (i = 0; i < nriver; i++)
    {
        fread(&river[i].restart_input, sizeof(river_rtic_struct), 1, fp);
    }

    fclose(fp);
}
