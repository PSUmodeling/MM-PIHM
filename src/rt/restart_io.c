#include "pihm.h"

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
