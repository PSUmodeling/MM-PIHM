#include "pihm.h"

void ReadIc(const char filename[], elem_struct elem[], river_struct river[])
{
    FILE           *ic_file;
    int             i;
    int             size;

    ic_file = pihm_fopen(filename, "rb");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    fseek(ic_file, 0L, SEEK_END);
    size = ftell(ic_file);

    if (size !=
        (int)(sizeof(ic_struct) * nelem + sizeof(river_ic_struct) * nriver))
    {
        pihm_printf(VL_ERROR,
            "Error in initial condition file %s.\n"
            "The file size does not match requirement.\n", filename);
        pihm_printf(VL_ERROR, "Please use a correct initial condition file.\n");
        pihm_exit(EXIT_FAILURE);
    }

    fseek(ic_file, 0L, SEEK_SET);

    for (i = 0; i < nelem; i++)
    {
        fread(&elem[i].ic, sizeof(ic_struct), 1, ic_file);
    }

    for (i = 0; i < nriver; i++)
    {
        fread(&river[i].ic, sizeof(river_ic_struct), 1, ic_file);
    }

    fclose(ic_file);
}
