#include "pihm.h"

void ReadIc(const char fn[], elem_struct elem[], river_struct river[])
{
    FILE           *fp;
    int             i;
    int             size;

    fp = pihm_fopen(fn, "rb");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    fseek(fp, 0L, SEEK_END);
    size = ftell(fp);

    if (size !=
        (int)(sizeof(ic_struct) * nelem + sizeof(river_ic_struct) * nriver))
    {
        pihm_printf(VL_ERROR,
            "Error in initial condition file %s.\n"
            "The file size does not match requirement.\n", fn);
        pihm_printf(VL_ERROR, "Please use a correct initial condition file.\n");
        pihm_exit(EXIT_FAILURE);
    }

    fseek(fp, 0L, SEEK_SET);

    for (i = 0; i < nelem; i++)
    {
        fread(&elem[i].ic, sizeof(ic_struct), 1, fp);
    }

    for (i = 0; i < nriver; i++)
    {
        fread(&river[i].ic, sizeof(river_ic_struct), 1, fp);
    }

    fclose(fp);
}
