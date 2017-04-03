#include "pihm.h"

int             verbose_mode;
int             debug_mode;
int             corr_mode;
char            project[MAXSTRING];
#ifdef _OPENMP
int             nthreads;
#endif

int main (int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    int             spec_output_mode = 0;

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    /*
     * Read command line arguments
     */
    ParseCmdLineParam (argc, argv, &spec_output_mode, outputdir);

    /*
     * Print AscII art
     */
    AsciiArt ();

    /*
     * Create output directory
     */
    CreateOutputDir (outputdir, spec_output_mode);

    /*
     * Run PIHM (or EnKF system)
     */
    PIHM (project, outputdir, 1);

    PIHMprintf (VL_NORMAL, "\nSimulation completed.\n");

    return (EXIT_SUCCESS);
}
