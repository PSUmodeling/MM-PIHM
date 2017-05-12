#include "pihm.h"

int             verbose_mode;
int             debug_mode;
char            project[MAXSTRING];

int main (int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    int             spec_output_mode = 0;
#ifdef _ENKF_
    int             id;
    int             p;
    int             ierr;

    ierr = MPI_Init (&argc, &argv);
    ierr = MPI_Comm_rank (MPI_COMM_WORLD, &id);
    ierr = MPI_Comm_size (MPI_COMM_WORLD, &p);

    #ifdef _DEBUG_
    PauseParal (id);
    #endif
#endif


#ifdef _ENKF_
    if (0 == id)
    {
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
#ifdef _ENKF_
    }
#endif

    /*
     * Run PIHM (or EnKF system)
     */
#ifdef _ENKF_
    PIHMParal (id, p, outputdir);
#else
    PIHM (project, outputdir, 1);
#endif

#ifdef _ENKF_
    ierr = MPI_Finalize ();
#endif

#ifdef _ENKF_
    if (0 == id)
    {
#endif
    PIHMprintf (VL_NORMAL, "\nSimulation completed.\n");
#ifdef _ENKF_
    }
#endif

    return (EXIT_SUCCESS);
}
