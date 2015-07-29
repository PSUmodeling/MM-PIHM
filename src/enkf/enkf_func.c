#include "pihm.h"
#include "enkf.h"

void JobHandout (int msg, int total_jobs)
{
    int             dest;

    for (dest = 1; dest < total_jobs + 1; dest++)
    {
        MPI_Send (&msg, 1, MPI_INT, dest, NE_TAG, MPI_COMM_WORLD);
        printf ("Job sent to %d\n", dest);
        fflush (stdout);
    }
}

void PrintEnKFStatus (int starttime, int endtime)
{
    time_t          rawtime;
    struct tm      *timestamp;

    printf ("Running ensemble members from ");

    rawtime = (int) starttime;
    timestamp = gmtime (&rawtime);
    printf ("%4.4d-%2.2d-%2.2d %2.2d:%2.2d to ",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);

    rawtime = (int) endtime;
    timestamp = gmtime (&rawtime);
    printf ("%4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);
}

void JobHandIn (int total_jobs)
{
    int             ierr;
    int             success;
    int             received = 0;
    int             source;

    MPI_Status      status;

    while (received < total_jobs)
    {
        ierr = MPI_Recv (&success, 1, MPI_INT, MPI_ANY_SOURCE, SUCCESS_TAG, MPI_COMM_WORLD, &status);
        received++;
        source = status.MPI_SOURCE;
        printf("PIHM job handed in from Node: %d\n", source);
    }
}
