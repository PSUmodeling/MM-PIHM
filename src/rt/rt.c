/******************************************************************************
* RT-Flux-PIHM is a finite volume based, reactive transport module that
* operate on top of the hydrological land surface processes described by
* Flux-PIHM. RT-Flux-PIHM track the transportation and reaction in a given
* watershed. PIHM-RT uses operator splitting technique to couple transport
* and reaction.
*
* RT-Flux-PIHM requires two additional input files:
*     a. chemical condition file:     projectname.chem
*     b. index of initial conditions: projectname.cini
*
*
*
* If you have any questions, concerns, suggestions, please contact me at
* the following address:
*
*     Developer: Chen Bao <baochen.d.s@gmail.com>
*     Version  : 0.2
*     Date     : Feb, 2014
*****************************************************************************/

#include "pihm.h"               // 09.23

// 09.23
// /* Begin system library calls */
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//#include <time.h>
//#include <assert.h>
//#include <sys/time.h>
// 09.23

// 09.23
// /* Begin other library calls */
// //#include "sundialstypes.h" /* realtype, integertype, booleantype defination */
//#include "sundials_types.h"   // 09.16
//#include "pihm.h"     /* Data Model and Variable Declarations     */
//#include "rt.h"           /* Data Model and Variable Declarations for chemical processes */
// 09.23


/* Begin global variable definition (MACRO) */
#define UNIT_C 1440
#define ZERO   1E-20
#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80
#define TIME_MIN  1E-5
#define EPS       0.05
#define INFTYSMALL  1E-6
#define RTdepth 5.0
#define MIN(a,b) (((a)<(b))? (a):(b))
#define MAX(a,b) (((a)>(b))? (a):(b))

/* Fucntion declarations finished   */

// Timer
static double timer()
{
    struct timeval  tp;
    gettimeofday(&tp, NULL);
    return ((double)(tp.tv_sec) + 1e-6 * tp.tv_usec);
}


void Monitor(realtype t, realtype stepsize, const pihm_struct pihm, Chem_Data CD)    // 09.30
{
    /* unit of t and stepsize: min */
    /* DS: model data              */
    //    FILE*        logfile = fopen("logfile/pihmlog.out","w");

    int             i = 0;
    struct tm      *timestamp;
    time_t         *rawtime;
    //double timelps, sumflux1, sumflux2, correction, swi, inv_swi, partratio, sumlateral, depth, hu, hg, hn, ht, imrecharge, flux, A;
    double          timelps;

    //Model_Data MD;
    //MD = (Model_Data)DS;

    pihm_struct     pihmMonitor;
    pihmMonitor = pihm;

    rawtime = (time_t *)malloc(sizeof(time_t));
    *rawtime = (int)(t * 60);
    timestamp = gmtime(rawtime);
    timelps = t - CD->StartTime;

    double          unit_c = stepsize / UNIT_C;

    //    fprintf(logfile,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\n",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
    //    fprintf(logfile," Time step is %6.4f\n", stepsize);

    double         *resflux = (double *)malloc(CD->NumOsv * sizeof(double));


    for (i = 0; i < CD->NumOsv; i++)
    {
        resflux[i] = 0.0;
    }

    for (i = 0; i < CD->NumFac; i++)
    {
        if (!CD->Flux[i].flux_type)
        {
            resflux[CD->Flux[i].nodeup - 1] -= CD->Flux[i].flux * unit_c;
            // sum lateral fluxes
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < CD->NumEle; i++)
    {
        // 11.08
        double          sumflux1, sumflux2, correction, A;

        // correct recharge in the saturated zone
        // The Err Dumper here is the recharge for saturated zone
        sumflux1 =
            (CD->Vcele[i].height_t -
            CD->Vcele[i].height_o) * pihmMonitor->elem[i].topo.area *
            CD->Vcele[i].porosity;
        // 09.30
        //fprintf(stderr, "  (Monitor) pihmMonitor->elem[i, %d].topo.area = %6.4f \n", i, pihmMonitor->elem[i].topo.area);
        sumflux2 = sumflux1 - resflux[i];
        correction =
            -sumflux2 * UNIT_C / stepsize /
            CD->Flux[CD->Vcele[i].ErrDumper].flux;
        A = -CD->Flux[CD->Vcele[i].ErrDumper].flux;
        CD->Flux[CD->Vcele[i].ErrDumper].flux = -sumflux2 * UNIT_C / stepsize;
        CD->Flux[CD->Vcele[i].ErrDumper].velocity =
            CD->Flux[CD->Vcele[i].ErrDumper].flux /
            CD->Flux[CD->Vcele[i].ErrDumper].s_area;
        CD->Flux[CD->Vcele[i].ErrDumper - CD->NumEle].flux =
            -CD->Flux[CD->Vcele[i].ErrDumper].flux;
        CD->Flux[CD->Vcele[i].ErrDumper - CD->NumEle].velocity =
            -CD->Flux[CD->Vcele[i].ErrDumper].velocity;
        // fprintf(logfile, "%d\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n", i+1, CD->Vcele[i].height_t, CD->Vcele[i].height_o, sumflux1, tmpflux[i], resflux[i], sumflux2, A, sumflux2 * UNIT_C/ stepsize , (sumflux1-resflux[i] + CD->Flux[CD->Vcele[i].ErrDumper].flux * unit_c )/sumflux1 * 100, correction);
    }

    // Since flux are corrected for saturated zones, resflux need re-calculate

    for (i = 0; i < CD->NumOsv; i++)
    {
        resflux[i] = 0.0;
    }

    for (i = 0; i < CD->NumFac; i++)
    {
        if (!CD->Flux[i].flux_type)
        {
            resflux[CD->Flux[i].nodeup - 1] -= CD->Flux[i].flux * unit_c;
            // sum lateral fluxes
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = CD->NumEle; i < CD->NumEle * 2; i++)
    {
        // 11.08
        double          sumflux1, sumflux2, correction, A;

        // correct infiltration in the unsaturated zone
        // The Err Dumper here is the infiltration for the saturated zone.

        //sumflux1 = (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * MD->Ele[i - CD->NumEle].area * CD->Vcele[i].porosity;  // 09.30
        sumflux1 =
            (CD->Vcele[i].height_t -
            CD->Vcele[i].height_o) * pihmMonitor->elem[i -
            CD->NumEle].topo.area * CD->Vcele[i].porosity;
        // 09.30
        //fprintf(stderr, "  (Monitor) pihmMonitor->elem[i, %d].topo.area = %6.4f \n", i, pihmMonitor->elem[i - CD->NumEle].topo.area);
        sumflux2 = sumflux1 - resflux[i];
        correction = -sumflux2 * UNIT_C / stepsize / CD->Vcele[i].q;
        A = CD->Flux[CD->Vcele[i].ErrDumper].flux;
        CD->Vcele[i].q = sumflux2 * UNIT_C / stepsize;
        CD->Vcele[i].q = MAX(CD->Vcele[i].q, 0.0);

        // input of rain water chemistry can not be negative;
        //CD->Vcele[i].q += fabs(MD->EleET[i - MD->NumEle][2]) * MD->Ele[i - MD->NumEle].area;  // 09.30
        CD->Vcele[i].q +=
            fabs(pihmMonitor->elem[i - CD->NumEle].wf.edir_unsat +
            pihmMonitor->elem[i -
                CD->NumEle].wf.edir_gw) * 86400 * pihmMonitor->elem[i -
            CD->NumEle].topo.area;
        // 09.30
        //fprintf(stderr, "  (Monitor) elem[%d].wf.edir_unsat = %g \n", (i - CD->NumEle), (pihmMonitor->elem[i - CD->NumEle].wf.edir_unsat + pihmMonitor->elem[i - CD->NumEle].wf.edir_gw) * 86400);

        // in addition, the soil evaporation leaves chemicals inside
        // The above code is , ensure the q term, which is the net input of water resulted from precipitation, should be net precipitation plus soil evaporation. Note
        // that soil evaporation itself might not be accurate in flux-PIHM. If flux-PIHM underestimates soil evaporation, RT need overestimate the incoming concentration
        // to compensate.


        //      fprintf(stderr, " %f = %f / %f, %f;\n",  MD->EleET[i-MD->NumEle][2]/ MD->EleNetPrep[i-MD->NumEle], MD->EleET[i-MD->NumEle][2],MD->EleNetPrep[i-MD->NumEle], MD->ElePrep[i-MD->NumEle]);
        //      CD->Condensation += fabs(MD->EleET[i-MD->NumEle][2]/ MD->EleNetPrep[i-MD->NumEle]);

        // q: in positier, out negative
        // flux: in negative, out positive
        //      if ( CD->Vcele[i].q < A) fprintf(logfile, " !!! A:%f, Q:%f, %f %f\n",A, CD->Vcele[i].q, sumflux2 * UNIT_C / stepsize, fabs(MD->EleET[i-MD->NumEle][2] * MD->Ele[i-MD->NumEle].area));
        //  fprintf(logfile, "%d\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n", i+1, CD->Vcele[i].height_t, CD->Vcele[i].height_o, MD->Ele[i-MD->NumEle].area, MD->Ele[i-MD->NumEle].Porosity, resflux[i], sumflux2, A, CD->Vcele[i].q , (sumflux1 - resflux[i] - CD->Vcele[i].q * unit_c ) , MD->EleET[i- MD->NumEle][2] * MD->Ele[i-MD->NumEle].area);
    }

    //    fclose(logfile);

    free(resflux);
    free(rawtime);
}

int
//upstream(element up, element lo, const void *DS)  // 09.26 update
upstream(elem_struct up, elem_struct lo, const pihm_struct pihm)
{
    /* Locate the upstream grid of up -> lo flow */
    /* Require verification                      */
    /* only determines points in triangular elements */

    double          x_, y_, x_a, x_b, x_c, y_a, y_b, y_c,
        dot00, dot01, dot02, dot11, dot12, u, v, invDenom;

    int             i, upstreamfound = 0;
    //x_ = 2 * up.x - lo.x;  // 09.26 update
    //y_ = 2 * up.y - lo.y;  // 09.26 update
    x_ = 2 * up.topo.x - lo.topo.x;
    y_ = 2 * up.topo.y - lo.topo.y;

    //Model_Data MD;         // 09.26 update
    //MD = (Model_Data)DS;   // 09.26 update
    pihm_struct     pihmStreamcopy;
    pihmStreamcopy = pihm;

    //  fprintf(stderr,"%d\t%d\t%d\t%d\t",MD->NumEle,MD->Ele[up.nabr[0]].index,MD->Ele[up.nabr[1]].index,MD->Ele[up.nabr[2]].index);

    for (i = 0; i < nelem; i++)
    {
        /* find point lies in which triangular element, a very interesting method */
        //if ((i != (up.index - 1)) && (i != (lo.index - 1)))  // 09.26 update
        if ((i != (up.ind - 1)) && (i != (lo.ind - 1)))
        {
            //x_a = MD->Node[MD->Ele[i].node[0] - 1].x;  // 09.26 update
            //x_b = MD->Node[MD->Ele[i].node[1] - 1].x;
            //x_c = MD->Node[MD->Ele[i].node[2] - 1].x;
            //y_a = MD->Node[MD->Ele[i].node[0] - 1].y;
            //y_b = MD->Node[MD->Ele[i].node[1] - 1].y;
            //y_c = MD->Node[MD->Ele[i].node[2] - 1].y;
            x_a =
                pihmStreamcopy->meshtbl.x[pihmStreamcopy->elem[i].node[0] - 1];
            x_b =
                pihmStreamcopy->meshtbl.x[pihmStreamcopy->elem[i].node[1] - 1];
            x_c =
                pihmStreamcopy->meshtbl.x[pihmStreamcopy->elem[i].node[2] - 1];
            y_a =
                pihmStreamcopy->meshtbl.y[pihmStreamcopy->elem[i].node[0] - 1];
            y_b =
                pihmStreamcopy->meshtbl.y[pihmStreamcopy->elem[i].node[1] - 1];
            y_c =
                pihmStreamcopy->meshtbl.y[pihmStreamcopy->elem[i].node[2] - 1];
            dot00 = (x_c - x_a) * (x_c - x_a) + (y_c - y_a) * (y_c - y_a);
            dot01 = (x_c - x_a) * (x_b - x_a) + (y_c - y_a) * (y_b - y_a);
            dot02 = (x_c - x_a) * (x_ - x_a) + (y_c - y_a) * (y_ - y_a);
            dot11 = (x_b - x_a) * (x_b - x_a) + (y_b - y_a) * (y_b - y_a);
            dot12 = (x_b - x_a) * (x_ - x_a) + (y_b - y_a) * (y_ - y_a);
            invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
            u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            v = (dot00 * dot12 - dot01 * dot02) * invDenom;
            if ((u > 0.0) && (v > 0.0) && (u + v < 1))
            {
                upstreamfound = 1;
                //      fprintf(stderr, "The upstream of %d and %d is %d!\n",up.index, lo.index, MD->Ele[i].index);
                //return (MD->Ele[i].index);  // 09.26 update
                return pihmStreamcopy->elem[i].ind;
            }
        }
    }

    return 0;
}

int realcheck(const char *words)
{

    int             flg = 1, i;
    if (((words[0] < 58) && (words[0] > 47)) || (words[0] == 46)
        || (words[0] == 45) || (words[0] == 43))
    {
        for (i = 0; i < (int)strlen(words); i++)
            if ((words[i] > 57 || words[i] < 43) && (words[i] != 69)
                && (words[i] != 101) && (words[i] != 10) && (words[i] != 13))
                flg = 0;
    }
    else
        flg = 0;
    return (flg);
}


int
keymatch(const char *line, const char *keyword, double *value, char **strval)
{
    /* A very general and convinient way of reading datafile and input file */
    /* find keyword in line, assign the value after keyword to value array if there is any */
    /* store both numbers and strings in order for later use, buffer required */
    /* if is keyword not found return 0. If comments, return 2. Otherwise return 1 */
    int             i;

    for (i = 0; i < WORDS_LINE; i++)
        value[i] = 0.0;

    if ((line[0] == '!') || (line[0] == '#'))
    {

        return (2);
        /* assign a special flag for comments */
    }

    int             j, k;
    int             words_line = WORDS_LINE;
    int             keyfoundflag = 0;

    char          **words;
    words = (char **)malloc(WORDS_LINE * sizeof(char *));

    for (i = 0; i < WORDS_LINE; i++)
    {
        words[i] = (char *)malloc(WORD_WIDTH * sizeof(char));
        memset(words[i], 0, WORD_WIDTH);
    }
    i = j = k = 0;

    /* Partition the line into words */
    while (i < (int)strlen(line))
    {
        if (line[i] != 39)
        {
            while (line[i] != 9 && line[i] != 0 && line[i] != 10
                && line[i] != 32 && line[i] != 13)
            {
                words[k][j++] = line[i++];
                if (line[i] == 9 || line[i] == 32 || line[i] == 13)
                {
                    k++;
                    j = 0;
                }
            }
        }
        else
        {
            words[k][j++] = line[i++];
            while (line[i] != 39)
            {
                words[k][j++] = line[i++];
            }
            words[k++][j] = line[i++];
            j = 0;
        }
        i++;
    }

    words_line = k + 1;

    for (i = 0; i < words_line; i++)
        if (strcmp(words[i], keyword) == 0)
            keyfoundflag = 1;

    j = k = 0;
    for (i = 0; i < words_line; i++)
    {
        //    fprintf(stderr, "word#%d=%s, length=%d\n" , i, words[i], strlen(words[i]));
        strcpy(strval[k++], words[i]);
        //    if ((( words[i][0] < 58)&&(words[i][0] > 47))||(words[i][0] == 46)||(words[i][0]==45)||(words[i][0]==43))
        if (realcheck(words[i]) == 1)
            value[j++] = atof(words[i]);
    }

    for (i = 0; i < WORDS_LINE; i++)
        free(words[i]);
    free(words);
    return (keyfoundflag);

}

void ConditionAssign(int condition, char *str, int *index)
{
    /* This subroutine takes in input strings and output an index array that record the conditions each blocks assigned to */
    /* input strings could use separators like - and , */

    int             i, j, k, l, length = strlen(str);
    char          **words = (char **)malloc(length * sizeof(char *));
    for (i = 0; i < length; i++)
        words[i] = (char *)malloc(WORD_WIDTH * sizeof(char));

    char           *tmpstr = (char *)malloc(WORD_WIDTH * sizeof(char));

    char           *separator = (char *)malloc(length * sizeof(char));
    int            *value = (int *)malloc(length * sizeof(char));

    i = j = k = l = 0;
    while (i < length)
    {
        while (str[i] != 0 && str[i] != 45 && str[i] != 44)
        {
            words[k][j++] = str[i++];
            if (str[i] == 45 || str[i] == 44)
            {
                k++;
                j = 0;
                separator[l++] = str[i];
            }
        }
        i++;
    }

    for (i = 0; i < length; i++)
    {
        strcpy(tmpstr, words[i]);
        value[i] = atoi(tmpstr);
    }
    /*
     * for ( i = 0; i < length; i ++)
     * fprintf(stderr, "%s\t%c\n",words[i],separator[i]);
     *
     * for ( i = 0; i < length; i ++)
     * fprintf(stderr, "%d\n",value[i]);
     *
     * fprintf(stderr, "condition = %d\n", condition);
     */
    for (i = 0; i < length; i++)
    {
        if (separator[i] == ',' || separator[i] == 0)
            index[value[i]] = condition;
        if (separator[i] == '-')
            for (j = value[i]; j <= value[i + 1]; j++)
                index[j] = condition;
    }

    for (i = 0; i < length; i++)
        free(words[i]);
    free(words);
    free(separator);
    free(value);
    free(tmpstr);

}




void
//chem_alloc(char *filename, const void *DS, const Control_Data * CS, Chem_Data CD, realtype t)       // original RT
//chem_alloc(char *filename, const Model_Data DS, const Control_Data * CS, Chem_Data CD, realtype t)  // 09.25 try change the second term
chem_alloc(char *filename, const pihm_struct pihm, N_Vector CV_Y, Chem_Data CD, realtype t) // 09.26 new MMPIHM
{

    int             i, j, k, num_face =
        0, num_species, num_mineral, num_ads, num_cex, num_other,
        num_conditions = 0;
    int             id;         // 01.21 for read-in '*.maxwater'
    int             line_width = LINE_WIDTH, words_line =
        WORDS_LINE, word_width = WORD_WIDTH;
    int             Global_diff = 0, Global_disp = 0;
    int             error_flag = 0, speciation_flg = 0, specflg;
    double          total_area = 0.0, tmpval[WORDS_LINE];
    time_t          rawtime;
    struct tm      *timeinfo;

    pihm_struct     pihmRTcopy;
    pihmRTcopy = pihm;

    // 09.25 new trick
    assert(pihm != NULL);
    assert(pihmRTcopy != NULL);

    char            line[256];
    char          **tmpstr = (char **)malloc(WORDS_LINE * sizeof(char *));

    timeinfo = (struct tm *)malloc(sizeof(struct tm));

    for (i = 0; i < words_line; i++)
        tmpstr[i] = (char *)malloc(WORD_WIDTH * sizeof(char));

    char           *chemfn =
        (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
    sprintf(chemfn, "input/%s/%s.chem", filename, filename);
    FILE           *chemfile = fopen(chemfn, "r");

    char           *datafn =
        (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
    sprintf(datafn, "input/%s/%s.cdbs", filename, filename);
    FILE           *database = fopen(datafn, "r");

    char           *forcfn =
        (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
    sprintf(forcfn, "input/%s/%s.prep", filename, filename);
    FILE           *prepconc = fopen(forcfn, "r");

    // 01.21
    char           *maxwaterfn =
        (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
    sprintf(maxwaterfn, "input/%s/%s.maxwater", filename, filename);
    FILE           *maxwater = fopen(maxwaterfn, "r");

    // 10.01 check input files
    //assert(chemfile != NULL);
    //assert(database != NULL);
    //assert(prepconc != NULL);

    if (chemfile == NULL)
    {
        fprintf(stderr, "\n  Fatal Error: %s.chem does not exist! \n",
            filename);
        exit(1);
    }

    if (database == NULL)
    {
        fprintf(stderr, "\n  Fatal Error: %s.cdbs does not exist! \n",
            filename);
        exit(1);
    }

    if (prepconc == NULL)
    {
        fprintf(stderr, "\n  Fatal Error: %s.prep does not exist! \n",
            filename);
        exit(1);
    }

    // 01.21
    if (maxwater == NULL)
    {
        fprintf(stderr, "\n  Fatal Error: %s.maxwater does not exist! \n",
            filename);
        exit(1);
    }

    // 09.25 test
    printf("\n RTtest: nelem = %d \n", nelem);
    printf(" RTtest: nriver = %d \n", nriver);
    printf(" RTtest: pihm->ctrl.starttime = %d [s] \n", pihm->ctrl.starttime);

    // *******************************************************************************************************
    // 09.25 beginning updating variables
    CD->NumVol = 2 * (nelem + nriver) + 2;
    CD->NumOsv = CD->NumVol - 2;
    CD->NumEle = nelem;
    CD->NumRiv = nriver;
    printf(" RTtest: CD->NumEle = %d \n", nelem);
    printf(" RTtest: CD->NumRiv = %d \n", nriver);

    /* default control variable if not found in input file */
    CD->StartTime = t;
    printf(" RTtest: CD->StartTime = %d [min] \n\n", (int)t);
    CD->TVDFlg = 1;
    CD->OutItv = 1;
    CD->Cementation = 1.0;
    CD->ACTmod = 0;
    CD->DHEdel = 0;
    CD->TEMcpl = 0;
    CD->EffAds = 0;
    CD->RelMin = 0;
    CD->AvgScl = 1;
    CD->CptFlg = 1;
    CD->RivOff = 0;
    CD->TimRiv = 1.0;
    CD->React_delay = 10;
    CD->Condensation = 1.0;
    CD->NumBTC = 0;
    CD->NumPUMP = 0;
    CD->SUFEFF = 1;
    CD->CnntVelo = 0.01;
    CD->TimLst = 0.0;           // 10.24 fix uninitialised value(s)

    // 09.25 reading "*.chem"
    // RUNTIME block
    fprintf(stderr, "\n Reading '%s.chem' RUNTIME: \n", filename);
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "RUNTIME", tmpval, tmpstr) != 1)
        fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        fgets(line, line_width, chemfile);
        if (keymatch(line, "tvd", tmpval, tmpstr) == 1)
        {
            if (strcmp(tmpstr[1], "false") == 0)
                CD->TVDFlg = 0;
            if (strcmp(tmpstr[1], "true") == 0)
                CD->TVDFlg = 1;
            if (strcmp(tmpstr[1], "false") && strcmp(tmpstr[1], "true"))
                fprintf(stderr, "  TVD FLAG INPUT ERROR! \n");
            fprintf(stderr, "  Total variation diminishing set to %d %s. \n",
                CD->TVDFlg, tmpstr[1]);
        }
        if (keymatch(line, "output", tmpval, tmpstr) == 1)
        {
            CD->OutItv = (int)tmpval[0];
            fprintf(stderr, "  Output interval set to %d hours. \n",
                CD->OutItv);
        }
        if (keymatch(line, "activity", tmpval, tmpstr) == 1)
        {
            CD->ACTmod = (int)tmpval[0];
            fprintf(stderr, "  Activity correction is set to %d. \n",
                CD->ACTmod);
            // 0 for unity activity coefficient and 1 for DH equation update
        }
        if (keymatch(line, "act_coe_delay", tmpval, tmpstr) == 1)
        {
            CD->DHEdel = (int)tmpval[0];
            fprintf(stderr,
                "  Activity coefficient update delay is set to %d. \n",
                CD->DHEdel);
            // 0 for delay and 1 for no delay (solving together )
        }
        if (keymatch(line, "thermo", tmpval, tmpstr) == 1)
        {
            CD->TEMcpl = (int)tmpval[0];
            fprintf(stderr, "  Coupling of thermo modelling is set to %d. \n",
                CD->DHEdel);
            // 0 for delay and 1 for no delay (solving together )
        }
        if (keymatch(line, "relmin", tmpval, tmpstr) == 1)
        {
            CD->RelMin = (int)tmpval[0];
            switch (CD->RelMin)
            {
                case 0:
                    fprintf(stderr,
                        "  Using absolute mineral volume fraction. \n");
                    break;
                case 1:
                    fprintf(stderr,
                        "  Using relative mineral volume fraction. \n");
                    break;
            }
        }
        if (keymatch(line, "effads", tmpval, tmpstr) == 1)
        {
            CD->EffAds = (int)tmpval[0];
            switch (CD->EffAds)
            {
                case 0:
                    fprintf(stderr, "  Using the normal adsorption model. \n");
                    break;
                case 1:
                    fprintf(stderr,
                        "  Using the coupled MIM and adsorption model. \n");
                    break;
                    // under construction.
            }
        }
        if (keymatch(line, "transport_only", tmpval, tmpstr) == 1)
        {
            CD->RecFlg = (int)tmpval[0];
            switch (CD->RecFlg)
            {
                case 0:
                    fprintf(stderr, "  Transport only mode disabled.\n");
                    break;
                case 1:
                    fprintf(stderr, "  Transport only mode enabled. \n");
                    break;
                    // under construction.
            }
        }
        if (keymatch(line, "precipitation", tmpval, tmpstr) == 1)
        {
            CD->PrpFlg = (int)tmpval[0];
            switch (CD->PrpFlg)
            {
                case 0:
                    fprintf(stderr, "  No precipitation condition. \n");
                    break;
                case 1:
                    fprintf(stderr,
                        "  Precipitation condition is to be specified. \n");
                    break;
                case 2:
                    fprintf(stderr,
                        "  Precipitation condition is specified via file *.prep. \n");
                    break;
                    // under construction.
            }
        }
        if (keymatch(line, "RT_delay", tmpval, tmpstr) == 1)
        {
            CD->Delay = (int)tmpval[0];
            fprintf(stderr,
                "  Flux-PIHM-RT will start after running PIHM for %d days. \n",
                CD->Delay);
            CD->Delay *= UNIT_C;
            // under construction.
        }
        if (keymatch(line, "Condensation", tmpval, tmpstr) == 1)
        {
            CD->Condensation = tmpval[0];
            fprintf(stderr,
                "  The concentrations of infiltrating rainfall is set to be %f times of concentrations in precipitation. \n",
                CD->Condensation);
            // under construction.
            //CD->Condensation *= CS->Cal.Prep_conc;  // 09.25 temporal comment-out
            fprintf(stderr,
                "  The concentrations of infiltrating rainfall is set to be %f times of concentrations in precipitation. \n",
                CD->Condensation);
        }
        if (keymatch(line, "AvgScl", tmpval, tmpstr) == 1)
        {
            CD->React_delay = tmpval[0];
            fprintf(stderr,
                "  Averaging window for asynchronous reaction %d. \n",
                CD->React_delay);
            // under construction.
        }
        if (keymatch(line, "SUFEFF", tmpval, tmpstr) == 1)
        {
            CD->SUFEFF = tmpval[0];
            fprintf(stderr, "  Effective surface area mode set to %d. \n\n",
                CD->SUFEFF);
            // under construction.
        }
        if (keymatch(line, "Mobile_exchange", tmpval, tmpstr) == 1)
        {
            CD->TimRiv = tmpval[0];
            fprintf(stderr, "  Ratio of immobile ion exchange site %f. \n",
                CD->TimRiv);
            // under construction.
        }

        if (keymatch(line, "Connectivity_threshold", tmpval, tmpstr) == 1)
        {
            CD->CnntVelo = tmpval[0];
            fprintf(stderr,
                "  Minimum velocity to be deemed as connected is %f m/d. \n",
                CD->CnntVelo);
            // under construction.
        }
    }

    // OUTPUT block
    fprintf(stderr, "\n Reading '%s.chem' OUTPUT: \n", filename);
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while ((keymatch(line, "OUTPUT", tmpval, tmpstr) != 1) && (!feof(chemfile)))
        fgets(line, line_width, chemfile);
    CD->NumBTC = tmpval[0];
    fprintf(stderr, "  %d breakthrough points specified. \n", CD->NumBTC);
    CD->BTC_loc = (int *)malloc(CD->NumBTC * sizeof(int));
    i = 0;
    fprintf(stderr, "  --");
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        fgets(line, line_width, chemfile);
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            CD->BTC_loc[i] = (int)tmpval[0] - 1;
            fprintf(stderr, " Grid %d ", CD->BTC_loc[i] + 1);
            i++;
        }
        if (i >= CD->NumBTC)
            break;
    }
    fprintf(stderr, "are breakthrough points.\n\n");

    // GLOBAL block
    fprintf(stderr, " Reading '%s.chem' GLOBAL: \n", filename);
    species         Global_type;
    Global_type.ChemName = (char *)malloc(WORD_WIDTH * sizeof(char));
    strcpy(Global_type.ChemName, "GLOBAL");

    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "GLOBAL", tmpval, tmpstr) != 1)
        fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        fgets(line, line_width, chemfile);
        if (keymatch(line, "t_species", tmpval, tmpstr) == 1)
        {
            CD->NumStc = (int)tmpval[0];
            fprintf(stderr, "  %d chemical species specified. \n", CD->NumStc);
            /* H2O is always a primary species */
        }
        if (keymatch(line, "s_species", tmpval, tmpstr) == 1)
        {
            CD->NumSsc = (int)tmpval[0];
            fprintf(stderr, "  %d secondary species specified. \n",
                (int)tmpval[0]);
        }
        if (keymatch(line, "minerals", tmpval, tmpstr) == 1)
        {
            CD->NumMin = (int)tmpval[0];
            fprintf(stderr, "  %d minerals specified. \n", CD->NumMin);
        }
        if (keymatch(line, "adsorption", tmpval, tmpstr) == 1)
        {
            CD->NumAds = (int)tmpval[0];
            fprintf(stderr, "  %d surface complexation specified. \n",
                CD->NumAds);
        }
        if (keymatch(line, "cation_exchange", tmpval, tmpstr) == 1)
        {
            CD->NumCex = (int)tmpval[0];
            fprintf(stderr, "  %d cation exchange specified. \n", CD->NumCex);
        }
        if (keymatch(line, "mineral_kinetic", tmpval, tmpstr) == 1)
        {
            CD->NumMkr = (int)tmpval[0];
            fprintf(stderr, "  %d mineral kinetic reaction(s) specified. \n",
                CD->NumMkr);
        }
        if (keymatch(line, "aqueous_kinetic", tmpval, tmpstr) == 1)
        {
            CD->NumAkr = (int)tmpval[0];
            fprintf(stderr, "  %d aqueous kinetic reaction(s) specified. \n",
                CD->NumAkr);
        }
        if (keymatch(line, "diffusion", tmpval, tmpstr) == 1)
        {
            fprintf(stderr, "  Diffusion coefficient = %g [cm2/s] \n",
                tmpval[0]);
            Global_type.DiffCoe = tmpval[0] * 60.0 * 60.0 * 24.0 / 10000.0;
            Global_diff = 1;
            /* Require unit conversion ! */
        }
        if (keymatch(line, "dispersion", tmpval, tmpstr) == 1)
        {
            fprintf(stderr, "  Dispersion coefficient = %2.2f [m] \n",
                tmpval[0]);
            Global_type.DispCoe = tmpval[0];
            Global_disp = 1;
            /* Set global flags to indicate the global values are present */
        }
        if (keymatch(line, "cementation", tmpval, tmpstr) == 1)
        {
            fprintf(stderr, "  Cementation factor = %2.1f \n", tmpval[0]);
            CD->Cementation = tmpval[0];
        }
        if (keymatch(line, "temperature", tmpval, tmpstr) == 1)
        {
            CD->Temperature = tmpval[0];
            fprintf(stderr, "  Temperature = %3.1f \n\n", CD->Temperature);
        }
    }

    CD->NumSpc = CD->NumStc - (CD->NumMin + CD->NumAds + CD->NumCex);
    /* the number of species that are mobile, later used in the OS3D subroutine */
    CD->NumSdc = CD->NumStc - (CD->NumMin);
    /* the number of species that others depend on */

    CD->Dependency = (double **)malloc(CD->NumSsc * sizeof(double *));
    for (i = 0; i < CD->NumSsc; i++)
    {
        CD->Dependency[i] = (double *)malloc(CD->NumSdc * sizeof(double));
        /* convert secondary species as an expression of primary species */
        for (j = 0; j < CD->NumSdc; j++)
            CD->Dependency[i][j] = 0.0;
    }

    CD->Dep_kinetic =
        (double **)malloc((CD->NumMkr + CD->NumAkr) * sizeof(double *));
    for (i = 0; i < CD->NumMkr + CD->NumAkr; i++)
    {
        CD->Dep_kinetic[i] = (double *)malloc(CD->NumStc * sizeof(double));
        /* express kinetic species as function of primary species */
        for (j = 0; j < CD->NumStc; j++)
            CD->Dep_kinetic[i][j] = 0.0;
    }

    CD->Dep_kinetic_all = (double **)malloc((CD->NumMin) * sizeof(double *));
    for (i = 0; i < CD->NumMin; i++)
    {
        CD->Dep_kinetic_all[i] = (double *)malloc(CD->NumStc * sizeof(double));
        /* Dependencies of minearls, all */
        for (j = 0; j < CD->NumStc; j++)
            CD->Dep_kinetic_all[i][j] = 0.0;
    }

    CD->Keq = (double *)malloc(CD->NumSsc * sizeof(double));
    CD->KeqKinect =
        (double *)malloc((CD->NumMkr + CD->NumAkr) * sizeof(double));
    CD->KeqKinect_all = (double *)malloc(CD->NumMin * sizeof(double));
    /* Keqs of equilibrium/ kinetic and kinetic all */

    CD->Totalconc = (double **)malloc(CD->NumStc * sizeof(double *));
    for (i = 0; i < CD->NumStc; i++)
        CD->Totalconc[i] =
            (double *)malloc((CD->NumStc + CD->NumSsc) * sizeof(double));
    /* convert total concentration as an expression of all species */

    CD->Totalconck = (double **)malloc(CD->NumStc * sizeof(double *));
    for (i = 0; i < CD->NumStc; i++)
        CD->Totalconck[i] =
            (double *)malloc((CD->NumStc + CD->NumSsc) * sizeof(double));
    /* convert total concentration as an expression of all species */

    for (i = 0; i < CD->NumStc; i++)
        for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
        {
            CD->Totalconc[i][j] = 0.0;
            CD->Totalconck[i][j] = 0.0;
        }
    num_species = CD->NumSpc;

    // INITIAL_CONDITIONS block
    fprintf(stderr, " Reading '%s.chem' INITIAL_CONDITIONS: \n", filename);
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "INITIAL_CONDITIONS", tmpval, tmpstr) != 1)
        fgets(line, line_width, chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            num_conditions++;
        }
        fgets(line, line_width, chemfile);
    }
    fprintf(stderr, "  %d conditions assigned. \n", num_conditions);

    char          **chemcon = (char **)malloc(num_conditions * sizeof(char *));
    for (i = 0; i < num_conditions; i++)
        chemcon[i] = (char *)malloc(word_width * sizeof(char));
    char         ***con_chem_name =
        (char ***)malloc((num_conditions + 1) * sizeof(char **));
    for (i = 0; i < num_conditions + 1; i++)
    {                           // all conditions + precipitation
        con_chem_name[i] = (char **)malloc(CD->NumStc * sizeof(char *));
        for (j = 0; j < CD->NumStc; j++)
            con_chem_name[i][j] = (char *)malloc(WORD_WIDTH * sizeof(char));
    }

    int            *condition_index =
        (int *)malloc((CD->NumVol + 1) * sizeof(int));
    /* when user assign conditions to blocks, they start from 1 */

    for (i = 0; i < CD->NumVol; i++)
        condition_index[i] = 0;

    vol_conc       *Condition_vcele =
        (vol_conc *) malloc(num_conditions * sizeof(vol_conc));
    for (i = 0; i < num_conditions; i++)
    {
        Condition_vcele[i].index = i + 1;
        Condition_vcele[i].t_conc =
            (double *)malloc(CD->NumStc * sizeof(double));
        Condition_vcele[i].p_conc =
            (double *)malloc(CD->NumStc * sizeof(double));
        Condition_vcele[i].p_para =
            (double *)malloc(CD->NumStc * sizeof(double));
        Condition_vcele[i].p_type = (int *)malloc(CD->NumStc * sizeof(int));
        Condition_vcele[i].s_conc = NULL;
        /* we do not input cocentration for secondary speices in rt */
        for (j = 0; j < CD->NumStc; j++)
        {
            Condition_vcele[i].t_conc[j] = ZERO;
            Condition_vcele[i].p_conc[j] = ZERO;
        }
    }

    if (CD->PrpFlg)
    {
        CD->Precipitation.t_conc =
            (double *)malloc(CD->NumStc * sizeof(double));
        CD->Precipitation.p_conc =
            (double *)malloc(CD->NumStc * sizeof(double));
        CD->Precipitation.p_para =
            (double *)malloc(CD->NumStc * sizeof(double));
        CD->Precipitation.p_type = (int *)malloc(CD->NumStc * sizeof(int));
        CD->Precipitation.s_conc = NULL;
        for (i = 0; i < CD->NumStc; i++)
        {
            CD->Precipitation.t_conc[i] = ZERO;
            CD->Precipitation.p_conc[i] = ZERO;
        }

    }

    CD->chemtype =
        (species *) malloc((CD->NumStc + CD->NumSsc) * sizeof(species));
    if (CD->chemtype == NULL)
        fprintf(stderr, " Memory allocation error\n");

    for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {
        if (Global_diff == 1)
            CD->chemtype[i].DiffCoe = Global_type.DiffCoe;
        /*
         * else
         * CD->chemtype[i].DiffCoe = ZERO;
         */
        /* in squre m per day */
        if (Global_disp == 1)
            CD->chemtype[i].DispCoe = Global_type.DispCoe;
        else
            CD->chemtype[i].DispCoe = ZERO;

        CD->chemtype[i].ChemName = (char *)malloc(WORD_WIDTH * sizeof(char));
        assert(CD->chemtype[i].ChemName != NULL);
        memset(CD->chemtype[i].ChemName, 0, WORD_WIDTH);
        CD->chemtype[i].Charge = 0.0;
        CD->chemtype[i].SizeF = 1.0;
        CD->chemtype[i].itype = 0;
    }

    k = 0;
    int             initfile = 0;
    FILE           *cheminitfile = NULL;
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "INITIAL_CONDITIONS", tmpval, tmpstr) != 1)
        fgets(line, line_width, chemfile);
    if (strcmp(tmpstr[1], "FILE") == 0)
    {                           // Initialize chemical distribution from file evoked. This will nullify all the condition assignment given in the next lines.
        // But for now, please keep those lines to let the code work.

        initfile = 1;
        //fprintf(stderr, "  Initializing the initial chemical distribution from file '%s'. \n\n", tmpstr[2]);      // from '*.chem'
        fprintf(stderr, "  Specifiying the initial chemical distribution from file '%s.cini'. \n", filename);   // 10.02 CORRECT

        //char cheminit[20];
        //strcpy(cheminit, "input/example/");
        //strcat(cheminit, tmpstr[2]);
        //cheminitfile = fopen(cheminit, "r");
        char           *cheminit =
            (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
        sprintf(cheminit, "input/%s/%s.cini", filename, filename);
        cheminitfile = fopen(cheminit, "r");

        if (cheminitfile == NULL)
        {
            fprintf(stderr, "  Fatal Error: %s.cini does not exist! \n",
                filename);
            exit(1);
        }
        else
        {
            fprintf(stderr, "  Reading the '%s.cini'!! \n", filename);
        }

        free(cheminit);         // 10.02
    }

    fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            strcpy(chemcon[k++], tmpstr[0]);
            if (initfile == 0)
            {
                fprintf(stderr, " Condition %s assigned to cells %s.\n",
                    chemcon[k - 1], tmpstr[1]);
                ConditionAssign(k, tmpstr[1], condition_index);
            }
        }
        fgets(line, line_width, chemfile);
    }
    if (initfile == 1)
    {
        for (i = 0; i < CD->NumVol; i++)
        {
            fscanf(cheminitfile, "%d %d", &k, condition_index + i + 1);
            // fprintf(stderr, "%6d %6d %6s\n", i+1, condition_index[i+1], chemcon[condition_index[i+1]-1]);
        }
    }

    if (cheminitfile != NULL)
        fclose(cheminitfile);

    // CONDITIONS block
    fprintf(stderr, "\n Reading '%s.chem' CONDITIONS: ", filename);
    for (i = 0; i < num_conditions; i++)
    {
        rewind(chemfile);
        num_species = 0;
        num_mineral = 0;
        num_ads = 0;
        num_cex = 0;
        num_other = 0;
        fgets(line, line_width, chemfile);
        while ((keymatch(line, "Condition", tmpval, tmpstr) != 1) ||
            (keymatch(line, chemcon[i], tmpval, tmpstr) != 1))
            fgets(line, line_width, chemfile);
        if (strcmp(tmpstr[1], chemcon[i]) == 0)
            fprintf(stderr, "\n  %s", line);
        fgets(line, line_width, chemfile);
        while (keymatch(line, "END", tmpval, tmpstr) != 1)
        {
            if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
            {
                specflg = SpeciationType(database, tmpstr[0]);

                if (specflg == 1)
                {
                    num_other = num_mineral + num_ads + num_cex;
                    Condition_vcele[i].t_conc[num_species - num_other] =
                        tmpval[0];
                    strcpy(con_chem_name[i][num_species - num_other],
                        tmpstr[0]);
                    //fprintf(stderr, "  %s\t\t %g \n", con_chem_name[i][num_species - num_other], tmpval[0]);
                    fprintf(stderr, "  %-28s %g \n",
                        con_chem_name[i][num_species - num_other], tmpval[0]);
                    Condition_vcele[i].p_type[num_species - num_other] = 1;
                }
                // arrange the concentration of the primary species in such a way that all the mobile species are at the beginning.
                if (specflg == 4)
                {
                    Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds +
                        CD->NumCex + num_mineral] = tmpval[0];
                    if (strcmp(tmpstr[2], "-ssa") == 0)
                        //Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds + CD->NumCex + num_mineral] = tmpval[1] * CS->Cal.SSA;  // 09.25 temporal comment-out
                        Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds +
                            CD->NumCex + num_mineral] = tmpval[1] * 1.0;
                    strcpy(con_chem_name[i][CD->NumSpc + CD->NumAds +
                            CD->NumCex + num_mineral], tmpstr[0]);
                    fprintf(stderr,
                        "  mineral %-20s %6.4f \t specific surface area \t%6.4f \n",
                        con_chem_name[i][CD->NumSpc + CD->NumAds + CD->NumCex +
                            num_mineral], tmpval[0], tmpval[1]);
                    Condition_vcele[i].p_type[CD->NumSpc + CD->NumAds +
                        CD->NumCex + num_mineral] = 4;
                    num_mineral++;
                }
                if ((tmpstr[0][0] == '>') || (specflg == 2))
                {               // adsorptive sites and species start with >
                    //Condition_vcele[i].t_conc[CD->NumSpc + num_ads] = tmpval[0] * CS->Cal.Site_den;    // 09.25 temporal comment-out
                    Condition_vcele[i].t_conc[CD->NumSpc + num_ads] =
                        tmpval[0] * 1.0;
                    Condition_vcele[i].p_type[CD->NumSpc + num_ads] = 2;
                    Condition_vcele[i].p_para[CD->NumSpc + num_ads] = 0;
                    // update when fill in the parameters for adsorption.
                    strcpy(con_chem_name[i][CD->NumSpc + num_ads], tmpstr[0]);
                    fprintf(stderr, "  surface complex %s\t\t%6.4f \n",
                        con_chem_name[i][CD->NumSpc + num_ads], tmpval[0]);
                    num_ads++;
                    // under construction
                }
                if (specflg == 3)
                {
                    Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds +
                        num_cex] = tmpval[0];
                    Condition_vcele[i].p_type[CD->NumSpc + CD->NumAds +
                        num_cex] = 3;
                    Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds +
                        num_cex] = 0;
                    // update when fill in the parameters for cation exchange.
                    strcpy(con_chem_name[i][CD->NumSpc + CD->NumAds + num_cex],
                        tmpstr[0]);
                    fprintf(stderr, "  cation exchange %s\t\t%6.4f \n",
                        con_chem_name[i][CD->NumSpc + CD->NumAds + num_cex],
                        tmpval[0]);
                    num_cex++;
                    // under construction
                }
                num_species++;
            }
            fgets(line, line_width, chemfile);
        }
    }

    // PRECIPITATION block
    fprintf(stderr, "\n Reading '%s.chem' PRECIPITATION: ", filename);
    if (CD->PrpFlg)
    {
        rewind(chemfile);
        fgets(line, line_width, chemfile);
        while (keymatch(line, "PRECIPITATION", tmpval, tmpstr) != 1)
            fgets(line, line_width, chemfile);
        fgets(line, line_width, chemfile);
        fprintf(stderr, " \n");
        fprintf(stderr, "  ---------------------------------\n");
        fprintf(stderr, "  The condition of precipitation is \n");
        fprintf(stderr, "  ---------------------------------\n");
        num_species = 0;
        num_mineral = 0;
        num_ads = 0;
        num_cex = 0;
        num_other = 0;
        while (keymatch(line, "END", tmpval, tmpstr) != 1)
        {
            if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
            {
                specflg = SpeciationType(database, tmpstr[0]);

                if (specflg == 1)
                {
                    num_other = num_mineral + num_ads + num_cex;
                    CD->Precipitation.t_conc[num_species - num_other] =
                        tmpval[0];
                    strcpy(con_chem_name[num_conditions][num_species -
                            num_other], tmpstr[0]);
                    fprintf(stderr, "  %-28s %g \n",
                        con_chem_name[num_conditions][num_species - num_other],
                        tmpval[0]);
                    CD->Precipitation.p_type[num_species - num_other] = 1;
                }
                // arrange the concentration of the primary species in such a way that all the mobile species are at the beginning.
                if (specflg == 4)
                {
                    CD->Precipitation.t_conc[CD->NumSpc + CD->NumAds +
                        CD->NumCex + num_mineral] = tmpval[0];
                    if (strcmp(tmpstr[2], "-ssa") == 0)
                        CD->Precipitation.p_para[CD->NumSpc + CD->NumAds +
                            CD->NumCex + num_mineral] = tmpval[1];
                    strcpy(con_chem_name[num_conditions][CD->NumSpc +
                            CD->NumAds + CD->NumCex + num_mineral], tmpstr[0]);
                    fprintf(stderr,
                        "  mineral %-20s %6.4f \t specific surface area %6.4f\n",
                        con_chem_name[num_conditions][CD->NumSpc + CD->NumAds +
                            CD->NumCex + num_mineral], tmpval[0], tmpval[1]);
                    CD->Precipitation.p_type[CD->NumSpc + CD->NumAds +
                        CD->NumCex + num_mineral] = 4;
                    num_mineral++;
                }
                if ((tmpstr[0][0] == '>') || (specflg == 2))
                {               // adsorptive sites and species start with >
                    CD->Precipitation.t_conc[CD->NumSpc + num_ads] = tmpval[0]; // this is the site density of the adsorptive species.
                    CD->Precipitation.p_type[CD->NumSpc + num_ads] = 2;
                    CD->Precipitation.p_para[CD->NumSpc + num_ads] = 0;
                    // update when fill in the parameters for adsorption.
                    strcpy(con_chem_name[num_conditions][CD->NumSpc + num_ads],
                        tmpstr[0]);
                    fprintf(stderr, " surface complex %s\t %6.4f\n",
                        con_chem_name[num_conditions][CD->NumSpc + num_ads],
                        tmpval[0]);
                    num_ads++;
                    // under construction
                }
                if (specflg == 3)
                {
                    CD->Precipitation.t_conc[CD->NumSpc + CD->NumAds +
                        num_cex] = tmpval[0];
                    CD->Precipitation.p_type[CD->NumSpc + CD->NumAds +
                        num_cex] = 3;
                    CD->Precipitation.p_para[CD->NumSpc + CD->NumAds +
                        num_cex] = 0;
                    // update when fill in the parameters for cation exchange.
                    strcpy(con_chem_name[num_conditions][CD->NumSpc +
                            CD->NumAds + num_cex], tmpstr[0]);
                    fprintf(stderr, " cation exchange %s\t %6.4f\n",
                        con_chem_name[num_conditions][CD->NumSpc + CD->NumAds +
                            num_cex], tmpval[0]);
                    num_cex++;
                    // under construction
                }
                num_species++;
            }
            fgets(line, line_width, chemfile);
        }

    }

    int             check_conditions_num;

    if (CD->PrpFlg)
        check_conditions_num = num_conditions + 1;
    else
        check_conditions_num = num_conditions;

    if (num_species != CD->NumStc)
        fprintf(stderr, " Number of species does not match indicated value!\n");

    for (i = 1; i < check_conditions_num; i++)
    {
        for (j = 0; j < num_species; j++)
        {
            if (strcmp(con_chem_name[i][j], con_chem_name[i - 1][j]) != 0)
            {
                error_flag = 1;
            }
        }
        if (error_flag == 1)
            fprintf(stderr,
                " The order of the chemicals in condition <%s> is incorrect!\n",
                chemcon[i - 1]);
    }

    // Primary species table
    fprintf(stderr,
        "\n Primary species and their types: [1], aqueous; [2], adsorption; [3], cation exchange; [4], mineral. \n");
    for (i = 0; i < CD->NumStc; i++)    // 08.19 Number of total species in the rt simulator
    {
        strcpy(CD->chemtype[i].ChemName, con_chem_name[0][i]);
        CD->chemtype[i].itype = Condition_vcele[0].p_type[i];
        fprintf(stderr, "  %-20s %10d\n", CD->chemtype[i].ChemName,
            CD->chemtype[i].itype);
    }

    // Precipitation conc table
    if (CD->PrpFlg)
    {
        fprintf(stderr, "\n Total concentraions in precipitataion: \n");
        for (i = 0; i < CD->NumSpc; i++)
        {
            if (!strcmp(con_chem_name[num_conditions][i], "pH"))
            {
                if (CD->Precipitation.t_conc[i] < 7)
                {
                    CD->Precipitation.t_conc[i] =
                        pow(10, -CD->Precipitation.t_conc[i]);
                }
                else
                {
                    CD->Precipitation.t_conc[i] =
                        -pow(10, CD->Precipitation.t_conc[i] - 14);
                }
            }
            // change the pH of precipitation into total concentraion of H
            // We skip the speciation for rain and assume it is OK to calculate this way.
            fprintf(stderr, "  %-20s %-10.3g [M] \n",
                con_chem_name[num_conditions][i], CD->Precipitation.t_conc[i]);
        }
    }

    // SECONDARY_SPECIES block
    fprintf(stderr, "\n Reading 'shp.chem' SECONDARY_SPECIES: \n");
    fprintf(stderr, "  Secondary species specified in the input file: \n");
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "SECONDARY_SPECIES", tmpval, tmpstr) != 1)
        fgets(line, line_width, chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
        {
            strcpy(CD->chemtype[num_species++].ChemName, tmpstr[0]);
            fprintf(stderr, "  %s \n", CD->chemtype[num_species - 1].ChemName);
        }
        fgets(line, line_width, chemfile);
    }

    // MINERALS block
    fprintf(stderr, "\n Reading 'shp.chem' MINERALS: \n");

    // int num_dep = 2;
    int             num_dep = 4;    // 08.19 maximum  # of dependence, monod, and inhibition term

    CD->kinetics =
        (Kinetic_Reaction *) malloc(CD->NumMkr * sizeof(Kinetic_Reaction));
    for (i = 0; i < CD->NumMkr; i++)
    {
        CD->kinetics[i].species = (char *)malloc(WORD_WIDTH * sizeof(char));
        CD->kinetics[i].Label = (char *)malloc(WORD_WIDTH * sizeof(char));
        CD->kinetics[i].biomass_species = (char *)malloc(WORD_WIDTH * sizeof(char));    // 08.19
        CD->kinetics[i].dep_species = (char **)malloc(num_dep * sizeof(char *));
        CD->kinetics[i].dep_power = (double *)malloc(num_dep * sizeof(double));
        CD->kinetics[i].monod_species = (char **)malloc(num_dep * sizeof(char *));  // 08.19
        CD->kinetics[i].monod_para = (double *)malloc(num_dep * sizeof(double));
        CD->kinetics[i].inhib_species = (char **)malloc(num_dep * sizeof(char *));  // 08.19
        CD->kinetics[i].inhib_para = (double *)malloc(num_dep * sizeof(double));
        CD->kinetics[i].dep_position = (int *)malloc(num_dep * sizeof(int));
        CD->kinetics[i].monod_position = (int *)malloc(num_dep * sizeof(int));  // 08.19
        CD->kinetics[i].inhib_position = (int *)malloc(num_dep * sizeof(int));  // 08.19
        for (j = 0; j < num_dep; j++)
        {
            CD->kinetics[i].dep_species[j] =
                (char *)malloc(WORD_WIDTH * sizeof(char));
            CD->kinetics[i].monod_species[j] = (char *)malloc(WORD_WIDTH * sizeof(char));   // 08.19
            CD->kinetics[i].inhib_species[j] = (char *)malloc(WORD_WIDTH * sizeof(char));   // 08.19
            CD->kinetics[i].dep_position[j] = 0;
            CD->kinetics[i].monod_position[j] = 0;  // 08.19
            CD->kinetics[i].inhib_position[j] = 0;  // 08.19
        }
    }

    k = 0;
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "MINERALS", tmpval, tmpstr) != 1)
        fgets(line, line_width, chemfile);
    fgets(line, line_width, chemfile);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            strcpy(CD->kinetics[k].species, tmpstr[0]);
            if (strcmp(tmpstr[1], "-label") == 0)
                strcpy(CD->kinetics[k].Label, tmpstr[2]);
            k++;
        }
        fgets(line, line_width, chemfile);
    }

    for (i = 0; i < k; i++)
        fprintf(stderr,
            "  Kinetic reaction on '%s' is specified, label '%s'. \n",
            CD->kinetics[i].species, CD->kinetics[i].Label);

    // Precipitation conc read in
    fprintf(stderr, "\n Reading 'shp.prep': \n");
    if (CD->PrpFlg == 2)
    {
        CD->TSD_prepconc = (TSD *) malloc(sizeof(TSD));
        fscanf(prepconc, "%s %d %d", &(CD->TSD_prepconc->name),
            &(CD->TSD_prepconc->index), &(CD->TSD_prepconc->length));

        CD->prepconcindex =
            (int *)malloc(CD->TSD_prepconc->index * sizeof(int));
        // here prepconc.index is used to save the number of primary species. Must be equal to the number of primary species specified before.
        for (i = 0; i < CD->TSD_prepconc->index; i++)
        {
            fscanf(prepconc, "%d", &(CD->prepconcindex[i]));
            if (CD->prepconcindex[i] > 0)
            {
                assert(CD->prepconcindex[i] <= CD->NumSpc);
                fprintf(stderr,
                    "  Precipitation conc of '%s' is a time series. \n",
                    CD->chemtype[CD->prepconcindex[i] - 1].ChemName);
            }
        }

        CD->TSD_prepconc->TS =
            (realtype **)malloc((CD->TSD_prepconc->length) *
            sizeof(realtype *));
        for (i = 0; i < CD->TSD_prepconc->length; i++)
        {
            CD->TSD_prepconc->TS[i] =
                (realtype *)malloc((1 +
                    CD->TSD_prepconc->index) * sizeof(realtype));
            fscanf(prepconc, "%d-%d-%d %d:%d:%d", &timeinfo->tm_year,
                &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour,
                &timeinfo->tm_min, &timeinfo->tm_sec);
            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            rawtime = timegm(timeinfo);
            CD->TSD_prepconc->TS[i][0] = (realtype)rawtime / 60 / 60 / 24;
            for (j = 0; j < CD->TSD_prepconc->index; j++)
            {
                fscanf(prepconc, "%lf", &CD->TSD_prepconc->TS[i][j + 1]);   // [i][0] stores the time
            }
            /*
             * fprintf(stderr, " %8.4f\t", CD->TSD_prepconc->TS[i][0] );
             * for ( j = 0 ; j < CD->NumSpc; j ++)
             * fprintf(stderr, " %d %d %6.4f\t", i, j , CD->TSD_prepconc->TS[i][j]);
             * fprintf(stderr, "\n");
             */
        }
        CD->TSD_prepconc->iCounter = 0;
    }

    // PUMP block
    // 02.12
    CD->CalGwinflux = pihmRTcopy->cal.gwinflux;



    // PUMP block
    fprintf(stderr, "\n Reading 'shp.chem' PUMP: \n");
    rewind(chemfile);
    fgets(line, line_width, chemfile);
    while ((keymatch(line, "PUMP", tmpval, tmpstr) != 1) && (!feof(chemfile)))
        fgets(line, line_width, chemfile);
    CD->NumPUMP = tmpval[0];
    fprintf(stderr, "  %d pumps specified. \n", CD->NumPUMP);
    CD->pumps = (Pump *) malloc(CD->NumPUMP * sizeof(Pump));
    i = 0;

    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        fgets(line, line_width, chemfile);
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            CD->pumps[i].Pump_Location = (int)tmpval[0];
            CD->pumps[i].Injection_rate = (double)tmpval[1];
            CD->pumps[i].Injection_conc = (double)tmpval[2];
            CD->pumps[i].flow_rate =
                CD->pumps[i].Injection_rate / CD->pumps[i].Injection_conc /
                365 * 1E-3;
            CD->pumps[i].Name_Species = (char *)malloc(20 * sizeof(char));
            strcpy(CD->pumps[i].Name_Species, tmpstr[1]);
            //      wrap(CD->pumps[i].Name_Species);
            CD->pumps[i].Position_Species = -1;
            for (j = 0; j < CD->NumStc; j++)
            {
                if (!strcmp(CD->pumps[i].Name_Species,
                        CD->chemtype[j].ChemName))
                {
                    CD->pumps[i].Position_Species = j;
                }
            }

            fprintf(stderr,
                "  -- Rate %g [moles/year] of '%s' (pos: %d) at Grid '%d' with a concentration of %g [M/L]. \n",
                CD->pumps[i].Injection_rate, CD->pumps[i].Name_Species,
                (CD->pumps[i].Position_Species + 1), CD->pumps[i].Pump_Location,
                CD->pumps[i].Injection_conc);
            fprintf(stderr, "  -- Flow rate is then %g [m3/d]. \n",
                CD->pumps[i].flow_rate);
            //      CD->pumps[i].Injection_rate *= 1E-3 / 365;

            // 02.12 calibration
            CD->pumps[i].Injection_rate =
                CD->pumps[i].Injection_rate * CD->CalGwinflux;
            CD->pumps[i].flow_rate = CD->pumps[i].flow_rate * CD->CalGwinflux;
            fprintf(stderr,
                "  -- after calibration: injection_rate %g [moles/year], flow _rate %g [m3/d], CD->CalGwinflux = %f. \n",
                CD->pumps[i].Injection_rate, CD->pumps[i].flow_rate,
                CD->CalGwinflux);
            i++;
        }
        if (i >= CD->NumPUMP)
            break;
    }
    // Ending of reading input files

    // 01.21 reading '*.maxwater' input file
    fprintf(stderr, "\n Reading 'coalcreek_952.maxwater': \n");
    CD->Vcele = (vol_conc *) malloc(CD->NumVol * sizeof(vol_conc));
    for (i = 0; i < CD->NumVol; i++)
    {
        CD->Vcele[i].maxwater = 0;  // 01.21 initialize, including ghost cells
    }

    fscanf(maxwater, "%*[^\n]%*c"); // jump over the first header line

    for (i = 0; i < nelem; i++) // 01.21, GW cells
    {
        fscanf(maxwater, "%d %lf", &id, &(CD->Vcele[i].maxwater));  // read-in
        //printf(" EleID = %d, maxwater = %3.2f \n", id, CD->Vcele[i].maxwater);
    }

    for (i = nelem; i < 2 * nelem; i++) // 01.21, unsat cells
    {
        CD->Vcele[i].maxwater = CD->Vcele[i - nelem].maxwater;
        //printf(" EleID = %d, maxwater = %3.2f \n", i+1, CD->Vcele[i].maxwater);
    }



    /* Initializing volumetric parameters, inherit from pihm
     * That is, if pihm is started from a hot start, rt is also
     * initialized with the hot data */

    //printf("\n CD->NumVol = %d \n", CD->NumVol);
    //CD->Vcele = (vol_conc *)malloc(CD->NumVol * sizeof(vol_conc));  02.06 debug maxwater

    // 10.24 initializing BC, fix uninitialised value(s)
    for (i = 0; i < CD->NumVol; i++)
    {
        CD->Vcele[i].BC = 0;
    }

    /* Initializing volumetrics for groundwater (GW) cells */
    fprintf(stderr, "\n Initializing 'GW' cells, Vcele [i, 0 ~ nelem]... \n");
    for (i = 0; i < nelem; i++)
    {
        CD->Vcele[i].height_v =
            pihmRTcopy->elem[i].topo.zmax - pihmRTcopy->elem[i].topo.zmin;
        //fprintf(stderr, "  pihmRTcopy->elem[i, %d]: height_v = %f, topo.zmax & zmin = %6.4f, %6.4f \n", i+1, CD->Vcele[i].height_v, pihmRTcopy->elem[i].topo.zmax, pihmRTcopy->elem[i].topo.zmin);
        CD->Vcele[i].height_o = NV_Ith(CV_Y, GW(i));
        CD->Vcele[i].height_t = NV_Ith(CV_Y, GW(i));
        //fprintf(stderr, "  CD->Vcele[i, %d].height_o = %6.4f \n", i, NV_Ith(CV_Y, GW(i)));
        CD->Vcele[i].area = pihmRTcopy->elem[i].topo.area;
        //fprintf(stderr, "  CD->Vcele[i, %d].area = %6.4f \n", i, pihmRTcopy->elem[i].topo.area);
        CD->Vcele[i].porosity = pihmRTcopy->elem[i].soil.smcmax;
        //fprintf(stderr, "  CD->Vcele[i, %d].porosity = %6.4f \n", i, pihmRTcopy->elem[i].soil.smcmax);
        // fprintf(stderr, " PIHM Pe: %6.4f, RT P: %6.4f\n", MD->Ele[i].Porosity, CD->Vcele[i].porosity);
        // Porosity in PIHM is Effective Porosity = Porosity - Residue Water Porosity
        // Porosity in RT   is total Porosity, therefore, the water height in the unsaturated zone needs be converted as well
        CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
        CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        CD->Vcele[i].sat = 1.0;
        CD->Vcele[i].sat_o = 1.0;
        CD->Vcele[i].temperature = pihmRTcopy->elem[i].attrib.meteo_type;
        //fprintf(stderr, "  CD->Vcele[i, %d].temperature = %d \n", i, pihmRTcopy->elem[i].attrib.meteo_type);

        CD->Vcele[i].reset_ref = 0;
        for (j = 0; j < 3; j++)
        {
            //fprintf(stderr, "  pihmRTcopy->meshtbl.nabr[i, %d][j, %d] = %d \n", i, j, pihmRTcopy->meshtbl.nabr[i][j]);
            if (condition_index[pihmRTcopy->meshtbl.nabr[i][j]] ==
                condition_index[i + 1])
                CD->Vcele[i].reset_ref = pihmRTcopy->meshtbl.nabr[i][j];
        }
    }

    /* Initializing volumetrics for unsaturated cells */
    fprintf(stderr,
        "\n Initializing 'UNSAT' cells, Vcele [i, nelem ~ 2*nelem]... \n");
    for (i = nelem; i < 2 * nelem; i++)
    {
        j = i - nelem;
        CD->Vcele[i].height_v = CD->Vcele[i - nelem].height_v;
        //fprintf(stderr, "  pihmRTcopy->elem[i, %d]: height_v = %f \n", i+1, CD->Vcele[i].height_v);
        CD->Vcele[i].height_o =
            (NV_Ith(CV_Y,
                UNSAT(j)) * (pihmRTcopy->elem[j].soil.smcmax -
                pihmRTcopy->elem[j].soil.smcmin) + (CD->Vcele[i].height_v -
                CD->Vcele[i -
                    nelem].height_o) * pihmRTcopy->elem[j].soil.smcmin) /
            (pihmRTcopy->elem[j].soil.smcmax);
        //fprintf(stderr, "  (unsat) CD->Vcele[i, %d].height_v = %6.4f \n", i, CD->Vcele[i - nelem].height_v);
        //fprintf(stderr, "  (unsat) pihmRTcopy->elem[j, %d].soil.smcmax = %6.4f \n", j, pihmRTcopy->elem[j].soil.smcmax);
        //fprintf(stderr, "  (unsat) pihmRTcopy->elem[j, %d].soil.smcmin = %6.4f \n", j, pihmRTcopy->elem[j].soil.smcmin);
        //fprintf(stderr, "  (unsat) MD->DummyY[i, %d] = NV_Ith(CV_Y, UNSAT(j, %d)), %6.4f \n", i, j, NV_Ith(CV_Y, UNSAT(j)));
        //fprintf(stderr, "  (unsat) CD->Vcele[i, %d].height_o =  %6.4f \n", i, CD->Vcele[i].height_o);

        CD->Vcele[i].height_t = CD->Vcele[i].height_o;
        CD->Vcele[i].area = CD->Vcele[i - nelem].area;
        //fprintf(stderr, "  CD->Vcele[i, %d].area = %6.4f \n", i, CD->Vcele[i - nelem].area);
        CD->Vcele[i].porosity = CD->Vcele[i - nelem].porosity;
        // Unsaturated zone has the same porosity as saturated zone
        CD->Vcele[i].sat =
            CD->Vcele[i].height_o / (CD->Vcele[i].height_v - CD->Vcele[i -
                nelem].height_o);
        //fprintf(stderr, "  (unsat) CD->Vcele[i, %d].sat = %6.4f \n", i, CD->Vcele[i].sat);
        CD->Vcele[i].sat_o = CD->Vcele[i].sat;
        CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
        CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        //fprintf(stderr, "  (unsat) CD->Vcele[i, %d].vol = %6.4f \n", i, CD->Vcele[i].vol);
        CD->Vcele[i].temperature = CD->Vcele[i - nelem].temperature;
        //fprintf(stderr, "  (unsat) CD->Vcele[i, %d].temperature = %6.4f \n", i, CD->Vcele[i].temperature);

        /* The saturation of unsaturated zone is the Hu divided by height of this cell */
        if (CD->Vcele[i].sat > 1.0)
            fprintf(stderr,
                "Fatal Error, Unsaturated Zone Initialization For RT Failed!\n");
        CD->Vcele[i].reset_ref = i - nelem + 1;
        // default reset reference of unsaturated cells are the groundwater cells underneath
    }

    CD->CalPorosity = pihmRTcopy->cal.porosity;
    //fprintf(stderr, "  CD->CalPorosity = %6.4f \n", pihmRTcopy->cal.porosity);

    // 02.12
    CD->CalPorosity = pihmRTcopy->cal.porosity;
    CD->CalRate = pihmRTcopy->cal.rate;
    CD->CalSSA = pihmRTcopy->cal.ssa;
    //CD->CalGwinflux = pihmRTcopy->cal.gwinflux;
    CD->CalPrcpconc = pihmRTcopy->cal.prcpconc;
    CD->CalInitconc = pihmRTcopy->cal.initconc;
    CD->CalXsorption = pihmRTcopy->cal.Xsorption;   // 03.06

    //fprintf(stderr, "  CD->CalXsorption = %6.4f \n", CD->CalXsorption);

    /* Initializing volumetrics for river cells */
    fprintf(stderr,
        "\n Initializing 'RIV' cells, Vcele [i, 2*nelem ~ 2*nelem + nriver]... \n");
    for (i = 2 * nelem; i < 2 * nelem + nriver; i++)
    {
        j = i - 2 * nelem;
        CD->Vcele[i].height_v = NV_Ith(CV_Y, RIVSTG(j));
        //fprintf(stderr, "  (riv) CD->Vcele[i, %d].height_v = %6.4f \n", i, CD->Vcele[i].height_v);
        CD->Vcele[i].height_o = CD->Vcele[i].height_v;
        CD->Vcele[i].height_t = CD->Vcele[i].height_o;
        CD->Vcele[i].area = pihmRTcopy->river[i - 2 * nelem].topo.area;
        //fprintf(stderr, "  (riv) CD->Vcele[i, %d].area = %6.4f \n", i, CD->Vcele[i].area);
        CD->Vcele[i].porosity = 1.0;
        CD->Vcele[i].sat = 1.0;
        CD->Vcele[i].sat_o = CD->Vcele[i].sat;
        CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
        CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        //fprintf(stderr, "  (riv) CD->Vcele[i, %d].vol = %6.4f \n", i, CD->Vcele[i].vol);
        CD->Vcele[i].reset_ref = i + nriver + 1;
        // default reset reference of river segments are the EBR underneath
    }

    /* Initializing volumetrics for river EBR cells */
    fprintf(stderr,
        "\n Initializing 'RIV EBR' cells, Vcele [i, 2*nelem + nriver ~ 2*nelem + 2*nriver]... \n");
    for (i = 2 * nelem + nriver; i < 2 * nelem + 2 * nriver; i++)
    {
        j = i - 2 * nelem - nriver;
        CD->Vcele[i].height_v = NV_Ith(CV_Y, RIVGW(j));;
        //fprintf(stderr, "  (riv EBR) CD->Vcele[i, %d].height_v = %6.4f \n", i, CD->Vcele[i].height_v);
        CD->Vcele[i].height_o = CD->Vcele[i].height_v;
        CD->Vcele[i].height_t = CD->Vcele[i].height_o;
        CD->Vcele[i].area = pihmRTcopy->river[i - 2 * nelem - nriver].topo.area;
        //fprintf(stderr, "  (riv EBR) CD->Vcele[i, %d].area = %6.4f \n", i, pihmRTcopy->riv[i - 2 * nelem - nriver].topo.area);
        CD->Vcele[i].porosity = 1.0;
        CD->Vcele[i].sat = 1.0;
        CD->Vcele[i].sat_o = CD->Vcele[i].sat;
        CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
        CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        //fprintf(stderr, "  (riv EBR) CD->Vcele[i, %d].vol = %6.4f \n", i, CD->Vcele[i].vol);
        CD->Vcele[i].reset_ref = 0;
        if (condition_index[pihmRTcopy->river[j].leftele] ==
            condition_index[i + 1])
            CD->Vcele[i].reset_ref = pihmRTcopy->river[j].leftele;
        if (condition_index[pihmRTcopy->river[j].rightele] ==
            condition_index[i + 1])
            CD->Vcele[i].reset_ref = pihmRTcopy->river[j].rightele;
        //fprintf(stderr, "  (riv EBR) pihmRTcopy->riv[j, %d].leftele & rightele = %d, %d \n", j, pihmRTcopy->riv[j].leftele, pihmRTcopy->riv[j].rightele);
    }

    tmpval[0] = 0.0;
    for (i = 0; i < nelem; i++)
    {
        tmpval[0] += CD->Vcele[i].height_v;
    }
    tmpval[0] = tmpval[0] / nelem;
    fprintf(stderr, "  Average bedrock depth is %f [m]. \n", tmpval[0]);

    for (i = 0; i < CD->NumSpc; i++)
        if (strcmp(CD->chemtype[i].ChemName, "pH") == 0)
        {
            strcpy(CD->chemtype[i].ChemName, "H+");
            speciation_flg = 1;
        }

    for (i = CD->NumOsv; i < CD->NumVol; i++)
    {
        //j = i - 2 * nelem - nriver;
        CD->Vcele[i].height_v = 1.0;
        CD->Vcele[i].height_o = 1.0;
        CD->Vcele[i].height_t = 1.0;
        CD->Vcele[i].area = 1.0;
        CD->Vcele[i].porosity = 1.0;
        CD->Vcele[i].sat = 1.0;
        CD->Vcele[i].sat_o = 1.0;
        CD->Vcele[i].vol_o = 1.0;
        CD->Vcele[i].vol = 1.0;
        CD->Vcele[i].reset_ref = 2 * nelem + nriver + 1;
    }

    /* Initializing concentration distributions */
    fprintf(stderr,
        "\n Initializing concentration, Vcele [i, 0 ~ NumVol]... \n");

    for (i = 0; i < CD->NumVol; i++)
    {
        CD->Vcele[i].index = i + 1;
        CD->Vcele[i].NumStc = CD->NumStc;
        CD->Vcele[i].NumSsc = CD->NumSsc;
        CD->Vcele[i].t_conc = (double *)malloc(CD->NumStc * sizeof(double));
        CD->Vcele[i].t_rate = (double *)malloc(CD->NumStc * sizeof(double));
        CD->Vcele[i].t_tol = (double *)malloc(CD->NumStc * sizeof(double));
        CD->Vcele[i].p_conc = (double *)malloc(CD->NumStc * sizeof(double));
        CD->Vcele[i].s_conc = (double *)malloc(CD->NumSsc * sizeof(double));
        CD->Vcele[i].p_actv = (double *)malloc(CD->NumStc * sizeof(double));
        CD->Vcele[i].s_actv = (double *)malloc(CD->NumSsc * sizeof(double));
        CD->Vcele[i].p_para = (double *)malloc(CD->NumStc * sizeof(double));
        CD->Vcele[i].p_type = (int *)malloc(CD->NumStc * sizeof(int));

        //NewCell(&(CD->Vcele[i]), CD->NumStc, CD->NumSsc);

        CD->Vcele[i].q = 0.0;
        CD->Vcele[i].illness = 0;
        //  fprintf(stderr, " Reference of this cell is %d \n", CD->Vcele[i].reset_ref);

        // q could be used for precipitataion volume
        // net precipitation into the unsaturated zone is equal to the infiltration - discharge - evaportransporation.
        // net precipitation defined here are not the same concept in the pihm itselt
        // under construction.

        for (j = 0; j < CD->NumStc; j++)
        {
            if ((speciation_flg == 1) &&
                (strcmp(CD->chemtype[j].ChemName, "H+") == 0))
            {
                CD->Vcele[i].p_conc[j] =
                    pow(10,
                    -(Condition_vcele[condition_index[i + 1] - 1].t_conc[j]));
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_type[j] = 1;
            }
            else if (CD->chemtype[j].itype == 4)
            {
                CD->Vcele[i].t_conc[j] =
                    Condition_vcele[condition_index[i + 1] - 1].t_conc[j];
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                CD->Vcele[i].p_actv[j] = 1.0;
                CD->Vcele[i].p_para[j] =
                    Condition_vcele[condition_index[i + 1] - 1].p_para[j];
                CD->Vcele[i].p_type[j] =
                    Condition_vcele[condition_index[i + 1] - 1].p_type[j];
            }
            else
            {
                // 02.12 calibration
                if (strcmp(CD->chemtype[j].ChemName, "DOC") == 0)
                {
                    CD->Vcele[i].t_conc[j] =
                        CD->CalInitconc * Condition_vcele[condition_index[i +
                            1] - 1].t_conc[j];
                    //printf(" cell = %d, spec = %s, t_conc = %g, CD->CalInitconc = %f \n", i+1, CD->chemtype[j].ChemName, CD->Vcele[i].t_conc[j], CD->CalInitconc);
                }
                else
                {
                    CD->Vcele[i].t_conc[j] =
                        Condition_vcele[condition_index[i + 1] - 1].t_conc[j];
                    //printf(" cell = %d, spec = %s, t_conc = %g \n", i+1, CD->chemtype[j].ChemName, CD->Vcele[i].t_conc[j]);
                }
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j] * 0.5;
                // 01.24 fix IAP[i] = nan bug
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_para[j] =
                    Condition_vcele[condition_index[i + 1] - 1].p_para[j];
                CD->Vcele[i].p_type[j] =
                    Condition_vcele[condition_index[i + 1] - 1].p_type[j];
            }
        }
        for (j = 0; j < CD->NumSsc; j++)
        {
            CD->Vcele[i].s_conc[j] = ZERO;
        }
    }

    // 09.28 Beginning configuring the connectivity for flux
    /*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  */
    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < 3; j++)
            if (pihmRTcopy->meshtbl.nabr[i][j] > 0) // 09.26 "0" is on boundary
            {
                num_face++;     // 09.28 GW cell faces
            }
        total_area += pihmRTcopy->elem[i].topo.area;
    }
    CD->PIHMFac = num_face;
    num_face *= 2;              // 09.28 plus infilitration and rechage
    num_face += 2 * nelem + 6 * nriver; // A river + EBR taks 6 faces
    CD->NumFac = num_face;
    //fprintf(stderr, "\n Num_face = %d. \n", num_face);
    fprintf(stderr, "\n Total area of the watershed is %f [m^2]. \n",
        total_area);

    for (i = 0; i < CD->NumPUMP; i++)
    {
        CD->pumps[i].flow_rate = CD->pumps[i].flow_rate;
        fprintf(stderr, "\n PUMP rate is specified %g [m^3/d]. \n",
            CD->pumps[i].flow_rate);
    }


    /* Configuring the lateral connectivity of GW grid blocks */
    fprintf(stderr,
        "\n Configuring the lateral connectivity of GW grid blocks... \n");

    CD->Flux = (face *) malloc(CD->NumFac * sizeof(face));
    k = 0;

    double          dist1, dist2, para_a, para_b, para_c, x_0, x_1, y_0, y_1;
    int             index_0, index_1, rivi;

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < 3; j++)
            if (pihmRTcopy->meshtbl.nabr[i][j] > 0)
            {
                /* node indicates the index of grid blocks, not nodes at corners */
                CD->Flux[k].nodeup = i + 1;
                CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
                CD->Flux[k].nodelo = pihmRTcopy->meshtbl.nabr[i][j];    // 09.26 good so far
                CD->Flux[k].nodeuu = upstream(pihmRTcopy->elem[i], pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] - 1], pihm); // 09.26 good so far
                CD->Flux[k].nodell = upstream(pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] - 1], pihmRTcopy->elem[i], pihm); // 09.26 good so far
                CD->Flux[k].flux_type = 0;
                CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
                CD->Flux[k].BC = 0;
                //fprintf(stderr, " [i=%d] CD->Flux[k, %d].nodeup, nodelo, nodeuu, nodell = %d, %d, %d, %d \n", i,k,CD->Flux[k].nodeup,CD->Flux[k].nodelo,CD->Flux[k].nodeuu,CD->Flux[k].nodell);

                if (CD->Flux[k].nodeuu > 0) // 09.26 good so far
                    CD->Flux[k].distuu =
                        sqrt(pow(pihmRTcopy->elem[i].topo.x -
                            pihmRTcopy->elem[CD->Flux[k].nodeuu - 1].topo.x,
                            2) + pow(pihmRTcopy->elem[i].topo.y -
                            pihmRTcopy->elem[CD->Flux[k].nodeuu - 1].topo.y,
                            2));
                else
                    CD->Flux[k].distuu = 0.0;

                if (CD->Flux[k].nodell > 0) // 09.26 good so far
                    CD->Flux[k].distll =
                        sqrt(pow(pihmRTcopy->elem[pihmRTcopy->
                                meshtbl.nabr[i][j] - 1].topo.x -
                            pihmRTcopy->elem[CD->Flux[k].nodell - 1].topo.x,
                            2) +
                        pow(pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                                1].topo.y -
                            pihmRTcopy->elem[CD->Flux[k].nodell - 1].topo.y,
                            2));
                else
                    CD->Flux[k].distll = 0.0;

                CD->Flux[k].distance =
                    sqrt(pow(pihmRTcopy->elem[i].topo.x -
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.x,
                        2) + pow(pihmRTcopy->elem[i].topo.y -
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.y, 2));
                //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].distuu, distll, distance  = %f, %f, %f \n", i, k, CD->Flux[k].distuu, CD->Flux[k].distll, CD->Flux[k].distance);

                CD->Flux[k].distuu =
                    (pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1].topo.x -
                    pihmRTcopy->elem[i].topo.x) / CD->Flux[k].distance;
                // x component of the normal vector flux direction
                CD->Flux[k].distll =
                    (pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1].topo.y -
                    pihmRTcopy->elem[i].topo.y) / CD->Flux[k].distance;
                // y component of the normal vector flux direction
                //fprintf(stderr, " CD->Flux[k, %d].distuu & distll = %f \n", k, CD->Flux[k].distuu, CD->Flux[k].distll);   // 09.26 good so far

                index_0 = pihmRTcopy->elem[i].node[(j + 1) % 3] - 1;    // 09.27 good so far
                index_1 = pihmRTcopy->elem[i].node[(j + 2) % 3] - 1;
                x_0 = pihmRTcopy->meshtbl.x[index_0];   // 09.27 good so far
                y_0 = pihmRTcopy->meshtbl.y[index_0];
                x_1 = pihmRTcopy->meshtbl.x[index_1];
                y_1 = pihmRTcopy->meshtbl.y[index_1];
                para_a = y_1 - y_0;
                para_b = x_0 - x_1;
                para_c = (x_1 - x_0) * y_0 - (y_1 - y_0) * x_0;
                //fprintf(stderr, " [i, %d; j, %d] para_a, para_b, para_c = %f, %f, %f \n", i, j, para_a, para_b, para_c);  // 09.27 good so far

                if ((para_a * pihmRTcopy->elem[i].topo.x +
                        para_b * pihmRTcopy->elem[i].topo.y +
                        para_c) * (para_a *
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.x +
                        para_b *
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.y + para_c) > 0.0)
                    fprintf(stderr, " Two points at the same side of edge!\n");
                dist1 =
                    fabs(para_a * pihmRTcopy->elem[i].topo.x +
                    para_b * pihmRTcopy->elem[i].topo.y +
                    para_c) / sqrt(para_a * para_a + para_b * para_b);
                dist2 =
                    fabs(para_a *
                    pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1].topo.x +
                    para_b * pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1].topo.y + para_c) / sqrt(para_a * para_a +
                    para_b * para_b);
                //fprintf(stderr, " [i, %d; j, %d] dist1, dist2 = %f, %f \n", i, j, dist1, dist2);   // 09.27 good so far

                if ((CD->Flux[k].distance < dist1) ||
                    (CD->Flux[k].distance < dist2))
                {
                    fprintf(stderr,
                        "\n Checking the distance calculation wrong, total is %f, from ele to edge is %f, from edge to neighbor is %f\n",
                        CD->Flux[k].distance, dist1, dist2);
                    fprintf(stderr,
                        " Ele coordinates, x: %f, y: %f\n Node_1 coordinates, x: %f, y: %f\n Node_2 coordinates, x: %f, y: %f\n Nabr x:%f, y:%f\n",
                        pihmRTcopy->elem[i].topo.x, pihmRTcopy->elem[i].topo.y,
                        x_0, y_0, x_1, y_1,
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.x,
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.y);
                    fprintf(stderr,
                        " node_1 x:%.3f y:%.3f\t node_2 x:%.3f y:%.3f node_3 x:%.3f y:%.3f nabr x:%.3f y:%.3f, no:%d\n",
                        pihmRTcopy->meshtbl.x[pihmRTcopy->elem[i].node[0] - 1],
                        pihmRTcopy->meshtbl.y[pihmRTcopy->elem[i].node[0] - 1],
                        pihmRTcopy->meshtbl.x[pihmRTcopy->elem[i].node[1] - 1],
                        pihmRTcopy->meshtbl.y[pihmRTcopy->elem[i].node[1] - 1],
                        pihmRTcopy->meshtbl.x[pihmRTcopy->elem[i].node[2] - 1],
                        pihmRTcopy->meshtbl.y[pihmRTcopy->elem[i].node[2] - 1],
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.x,
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.y, j);
                }
                else
                {
                    CD->Flux[k].distance = dist1;
                }

                for (rivi = 0; rivi < nriver; rivi++)
                {
                    if ((CD->Flux[k].nodeup == pihmRTcopy->river[rivi].leftele)
                        && (CD->Flux[k].nodelo ==
                            pihmRTcopy->river[rivi].rightele))
                    {
                        CD->Flux[k].nodelo = 2 * nelem + nriver + rivi + 1;
                        CD->Flux[k].nodeuu = 0;
                        CD->Flux[k].nodell = 0;
                        CD->Flux[k].distance = dist1;
                        //fprintf(stderr, " Riv - Watershed connection corrected, now from cell %d to cell %d\n", CD->Flux[k].nodeup, CD->Flux[k].nodelo);
                    }
                    if ((CD->Flux[k].nodeup == pihmRTcopy->river[rivi].rightele)
                        && (CD->Flux[k].nodelo ==
                            pihmRTcopy->river[rivi].leftele))
                    {
                        CD->Flux[k].nodelo = 2 * nelem + nriver + rivi + 1;
                        CD->Flux[k].nodeuu = 0;
                        CD->Flux[k].nodell = 0;
                        CD->Flux[k].distance = dist1;
                        //fprintf(stderr, " Riv - Watershed connection corrected, now from cell %d to cell %d\n", CD->Flux[k].nodeup, CD->Flux[k].nodelo);
                    }
                }

                //fprintf(stderr, " CD->Flux[k, %d].nodeup, nodelo, nodeuu, nodell = %d, %d, %d, %d \n", k, CD->Flux[k].nodeup, CD->Flux[k].nodelo, CD->Flux[k].nodeuu, CD->Flux[k].nodell);
                //fprintf(stderr, " CD->Flux[k, %d].distance = %f \n", k, CD->Flux[k].distance);   // 09.27 good so far
                k++;
            }
    }

    // 09.28 11:52 am, the following is for UNSAT face/flux
    CD->NumLateral = k;

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < 3; j++)
            if (pihmRTcopy->meshtbl.nabr[i][j] > 0)
            {
                /* node indicates the index of grid blocks, not nodes at corners */
                CD->Flux[k].nodeup = i + 1 + nelem;
                CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
                CD->Flux[k].nodelo = pihmRTcopy->meshtbl.nabr[i][j] + nelem;
                CD->Flux[k].nodeuu = upstream(pihmRTcopy->elem[i], pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] - 1], pihm) + nelem; // 09.28 note the "+nelem"
                CD->Flux[k].nodell =
                    upstream(pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1], pihmRTcopy->elem[i], pihm) + nelem;
                CD->Flux[k].flux_type = 0;  // 02.21 unsat lateral flux, flux = 0.0
                CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
                CD->Flux[k].BC = 0;

                //if (CD->Flux[k].nodeuu  > 0)  // 09.28 11:53 fix Chen's bug
                if (CD->Flux[k].nodeuu - nelem > 0)
                {
                    CD->Flux[k].distuu =
                        sqrt(pow(pihmRTcopy->elem[i].topo.x -
                            pihmRTcopy->elem[CD->Flux[k].nodeuu - nelem -
                                1].topo.x,
                            2) + pow(pihmRTcopy->elem[i].topo.y -
                            pihmRTcopy->elem[CD->Flux[k].nodeuu - nelem -
                                1].topo.y, 2));
                }
                else
                    CD->Flux[k].distuu = 0.0;

                //if (CD->Flux[k].nodell > 0)  // 09.28 11:53 fix Chen's bug
                if (CD->Flux[k].nodell - nelem > 0)
                    CD->Flux[k].distll =
                        sqrt(pow(pihmRTcopy->elem[pihmRTcopy->
                                meshtbl.nabr[i][j] - 1].topo.x -
                            pihmRTcopy->elem[CD->Flux[k].nodell - nelem -
                                1].topo.x,
                            2) +
                        pow(pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                                1].topo.y -
                            pihmRTcopy->elem[CD->Flux[k].nodell - nelem -
                                1].topo.y, 2));
                else
                    CD->Flux[k].distll = 0.0;

                CD->Flux[k].distance =
                    sqrt(pow(pihmRTcopy->elem[i].topo.x -
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.x,
                        2) + pow(pihmRTcopy->elem[i].topo.y -
                        pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                            1].topo.y, 2));
                //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].distuu, distll, distance  = %f, %f, %f \n", i, k, CD->Flux[k].distuu, CD->Flux[k].distll, CD->Flux[k].distance); // 09.28 good

                index_0 = pihmRTcopy->elem[i].node[(j + 1) % 3] - 1;
                index_1 = pihmRTcopy->elem[i].node[(j + 2) % 3] - 1;
                x_0 = pihmRTcopy->meshtbl.x[index_0];
                y_0 = pihmRTcopy->meshtbl.y[index_0];
                x_1 = pihmRTcopy->meshtbl.x[index_1];
                y_1 = pihmRTcopy->meshtbl.y[index_1];
                para_a = y_1 - y_0;
                para_b = x_0 - x_1;
                para_c = (x_1 - x_0) * y_0 - (y_1 - y_0) * x_0;
                dist1 =
                    fabs(para_a * pihmRTcopy->elem[i].topo.x +
                    para_b * pihmRTcopy->elem[i].topo.y +
                    para_c) / sqrt(para_a * para_a + para_b * para_b);
                dist2 =
                    fabs(para_a *
                    pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1].topo.x +
                    para_b * pihmRTcopy->elem[pihmRTcopy->meshtbl.nabr[i][j] -
                        1].topo.y + para_c) / sqrt(para_a * para_a +
                    para_b * para_b);

                CD->Flux[k].distance = dist1;
                for (rivi = 0; rivi < nriver; rivi++)
                {
                    if ((CD->Flux[k].nodeup - nelem ==
                            pihmRTcopy->river[rivi].leftele) &&
                        (CD->Flux[k].nodelo - nelem ==
                            pihmRTcopy->river[rivi].rightele))
                    {
                        CD->Flux[k].nodelo = 2 * nelem + nriver + rivi + 1;
                        CD->Flux[k].nodeuu = 0;
                        CD->Flux[k].nodell = 0;
                        CD->Flux[k].distance = dist1;
                    }
                    if ((CD->Flux[k].nodeup - nelem ==
                            pihmRTcopy->river[rivi].rightele) &&
                        (CD->Flux[k].nodelo - nelem ==
                            pihmRTcopy->river[rivi].leftele))
                    {
                        CD->Flux[k].nodelo = 2 * nelem + nriver + rivi + 1;
                        CD->Flux[k].nodeuu = 0;
                        CD->Flux[k].nodell = 0;
                        CD->Flux[k].distance = dist1;
                    }
                }
                //fprintf(stderr, " CD->Flux[k, %d].nodeup, nodelo, nodeuu, nodell = %d, %d, %d, %d \n", k, CD->Flux[k].nodeup, CD->Flux[k].nodelo, CD->Flux[k].nodeuu, CD->Flux[k].nodell);
                //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].distance  = %f \n", i, k, CD->Flux[k].distance);  // 09.28 good so far
                k++;
            }
    }

    /* Configuring the vertical connectivity of UNSAT - GW blocks */
    fprintf(stderr,
        "\n Configuring the vertical connectivity of UNSAT - GW grid blocks... \n");

    /* centered at unsat blocks */
    for (i = nelem; i < 2 * nelem; i++)
    {
        CD->Vcele[i].ErrDumper = k;
        CD->Flux[k].nodeup = i + 1;
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodelo = i + 1 - nelem;

        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;  // 02.21 although vertical flow, but intended for infiltration correction in Monitor ( ); 0 for lateral, 1 for vertical
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].BC = 0;
        CD->Flux[k].distance = 0.1 * CD->Vcele[i].height_v;
        //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].nodeup, nodelo, nodeuu, nodell, distance  = %d, %d, %d, %d, %f \n", i, k, CD->Flux[k].nodeup, CD->Flux[k].nodelo, CD->Flux[k].nodeuu, CD->Flux[k].nodell, CD->Flux[k].distance);  // 09.28 good so far
        k++;
    }

    /* centered at gw blocks */
    for (i = 0; i < nelem; i++)
    {
        CD->Vcele[i].ErrDumper = k;
        CD->Flux[k].nodeup = i + 1;
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodelo = i + 1 + nelem;
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 1;  // 02.21 vertical flow, but intended for recharge correction in Monitor ( ); 0 for lateral, 1 for vertical
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].BC = 0;
        CD->Flux[k].distance = 0.1 * CD->Vcele[i].height_v;
        //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].nodeup, nodelo, nodeuu, nodell, distance  = %d, %d, %d, %d, %f \n", i, k, CD->Flux[k].nodeup, CD->Flux[k].nodelo, CD->Flux[k].nodeuu, CD->Flux[k].nodell, CD->Flux[k].distance);  // 09.28 good so far
        k++;
    }
    CD->NumDis = k;

    //printf(" CD->NumDis = %d \n", k);

    /* Configuring the connectivity of RIVER and EBR blocks */
    fprintf(stderr,
        "\n Configuring the connectivity of RIVER & EBR grid blocks... \n");

    /* Between River and Left */
    // River to left OFL 2
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[k].nodeup = i + 2 * nelem + 1 + nriver;
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodelo = CD->NumVol;
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;
        CD->Flux[k].BC = 1;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].distance = 1.0;
        CD->Flux[k].s_area = 0.0;
        k++;
    }
    //printf(" river2 flux = %d, CD->NumVol = %d \n", k-1, CD->NumVol);

    /* Between River and Right */
    // River to right OFL 3
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[k].nodeup = i + 2 * nelem + 1 + nriver;
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodelo = CD->NumVol;
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;
        CD->Flux[k].BC = 1;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].distance = 1.0;
        CD->Flux[k].s_area = 0.0;
        k++;
    }
    //printf(" river3 flux = %d \n", k-1);

    /* Between Left and EBR */
    // EBR to left  7 + 4
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[k].nodeup = i + 2 * nelem + nriver + 1;
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodelo = pihmRTcopy->river[i].leftele;
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;
        CD->Flux[k].BC = 0;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].distance = 1.0;
        CD->Flux[k].s_area =
            pihmRTcopy->river[i].shp.length *
            CD->Vcele[pihmRTcopy->river[i].leftele - 1].height_v;
        //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].nodeup, nodelo, riv_length, s_area = %d, %d, %f, %f \n", i, k, CD->Flux[k].nodeup, CD->Flux[k].nodelo, pihmRTcopy->river[i].shp.length, CD->Flux[k].s_area);  // 09.28 good so far
        k++;
    }
    //printf(" river4+7 flux = %d \n", k-1);

    /* Between Right and EBR */
    // EBR to right 8 + 5
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[k].nodeup = i + 2 * nelem + nriver + 1;
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodelo = pihmRTcopy->river[i].rightele;
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;
        CD->Flux[k].BC = 0;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].distance = 1.0;
        CD->Flux[k].s_area =
            pihmRTcopy->river[i].shp.length *
            CD->Vcele[pihmRTcopy->river[i].rightele - 1].height_v;
        //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].nodeup, nodelo, riv_length, s_area = %d, %d, %f, %f \n", i, k, CD->Flux[k].nodeup, CD->Flux[k].nodelo, pihmRTcopy->river[i].shp.length, CD->Flux[k].s_area);  // 09.28 good so far
        k++;
    }
    //printf(" river5+8 flux = %d \n", k-1);

    /* Between EBR */
    // To downstream EBR 9
    for (i = 0; i < nriver; i++)
    {
        // 01.12 fix bug
        //if (i == nriver - 1)
        if (pihmRTcopy->river[i].down < 0)
        {
            CD->Flux[k].nodeup = i + 2 * nelem + nriver + 1;
            CD->Flux[k].nodelo = CD->NumVol;
        }
        else
        {
            CD->Flux[k].nodeup = i + 2 * nelem + nriver + 1;
            // 01.12 fix bug
            //CD->Flux[k].nodelo = i + 2 * nelem + nriver + 2;
            // 02.21 fix critical bug
            //CD->Flux[k].nodelo = pihmRTcopy->river[i].down;
            CD->Flux[k].nodelo = 2 * nelem + nriver + pihmRTcopy->river[i].down;    // 02.21 fix
        }
        CD->Flux[k].node_trib = 0;  // 01.14 by Wei Zhi
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;
        CD->Flux[k].BC = 1;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].distance = 1.0;
        CD->Flux[k].s_area = 0.0;
        //fprintf(stderr, "  [i=%d] CD->Flux[k, %d].nodeup, nodelo = %d, %d \n", i, k, CD->Flux[k].nodeup, CD->Flux[k].nodelo);  // 09.28 good so far  02.21 fix bug, now OK
        k++;
    }
    //printf(" river9 flux = %d \n", k-1);

    // From upstream EBR 10
    for (i = 0; i < nriver; i++)
    {
        // 01.12 fix bug
        //if (i == 0)
        if (pihmRTcopy->river[i].up[0] < 0)     /* No upstream */
        {
            CD->Flux[k].nodelo = CD->NumVol;
            CD->Flux[k].nodeup = i + 2 * nelem + nriver + 1;
        }
        else
        {
            CD->Flux[k].nodeup = i + 2 * nelem + nriver + 1;
            CD->Flux[k].nodelo = 2 * nelem + nriver + pihmRTcopy->river[i].up[0];
        }
        if (pihmRTcopy->river[i].up[1] < 0)     /* Only one upstream segment */
        {
            CD->Flux[k].node_trib = pihmRTcopy->river[i].up[1];
        }
        else
        {
            CD->Flux[k].node_trib = 2 * nelem + nriver + pihmRTcopy->river[i].up[1];
        }
        CD->Flux[k].nodeuu = 0;
        CD->Flux[k].nodell = 0;
        CD->Flux[k].flux_type = 0;
        CD->Flux[k].BC = 1;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].distance = 1.0;
        CD->Flux[k].s_area = 0.0;
        k++;
    }
    //printf(" river10 flux = %d \n", k-1);

    for (k = 0; k < CD->NumFac; k++)
    {
        CD->Flux[k].velocity = 0.0;
        CD->Flux[k].flux = 0.0; // 01.14 initialize 0.0 for above sections of GW-GW, UNSAT-UNSAT, GW-UNSAT, UNSAT-GW
        CD->Flux[k].flux_trib = 0.0;    // 01.14 by Wei Zhi
        CD->Flux[k].s_area = 0.0;
    }

    CD->SPCFlg = speciation_flg;
    Lookup(database, CD, 0);
    // update the concentration of mineral after get the molar volume of minera

    // 09.28 temporal comment-out
    double          Cal_PCO2 = 1.0;
    double          Cal_Keq = 1.0;
    for (i = 0; i < CD->NumAkr + CD->NumMkr; i++)
    {
        if (!strcmp(CD->chemtype[i + CD->NumSpc + CD->NumAds +
                    CD->NumCex].ChemName, "'CO2(*g)'"))
        {
            CD->KeqKinect[i] += log10(Cal_PCO2);
        }
        else
        {
            CD->KeqKinect[i] += log10(Cal_Keq);
        }
    }

    // 08.19 comment out by Wei, why print cell 200 concentration
    // for (i = 0; i < CD->NumStc; i++)
    //{
    // fprintf (stderr, " \n Conc and SSA of each species: %6.4g %6.4g\n",
    //CD->Vcele[200].t_conc[i], CD->Vcele[200].p_para[i]);
    // }

    fprintf(stderr, "\n Kinetic Mass Matrx (calibrated Keq)! \n\t\t");  // 08.19 Wei
    for (i = 0; i < CD->NumStc; i++)
        fprintf(stderr, " %6s\t", CD->chemtype[i].ChemName);
    fprintf(stderr, "\n");
    for (j = 0; j < CD->NumMkr + CD->NumAkr; j++)
    {
        fprintf(stderr, " %6s\t",
            CD->chemtype[j + CD->NumSpc + CD->NumAds + CD->NumCex].ChemName);
        for (i = 0; i < CD->NumStc; i++)
        {
            fprintf(stderr, " %6.2f\t", CD->Dep_kinetic[j][i]);
        }
        fprintf(stderr, "\tKeq = %6.2f\n", CD->KeqKinect[j]);
    }
    fprintf(stderr, "\n");
    // use calibration coefficient to produce new Keq values for 1) CO2, 2) other kinetic reaction

    fprintf(stderr,
        " \n Mass action species type determination (0: immobile, 1: mobile, 2: Mixed) \n");
    for (i = 0; i < CD->NumSpc; i++)
    {
        if (CD->chemtype[i].itype == 1)
            CD->chemtype[i].mtype = 1;
        else
            CD->chemtype[i].mtype = 0;
        for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
        {
            if ((CD->Totalconc[i][j] != 0) &&
                (CD->chemtype[j].itype != CD->chemtype[i].mtype))
                CD->chemtype[i].mtype = 2;
        }
        /*
         * if (strcmp( CD->chemtype[i].ChemName, "'H+'") == 0)
         * CD->chemtype[i].mtype = 1;
         */
        fprintf(stderr, " %12s\t%10d\n", CD->chemtype[i].ChemName,
            CD->chemtype[i].mtype);
    }

    fprintf(stderr,
        " \n Individual species type determination (1: aqueous, 2: adsorption, 3: ion exchange, 4: solid) \n");
    for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {
        fprintf(stderr, " %12s\t%10d\n", CD->chemtype[i].ChemName,
            CD->chemtype[i].itype);
    }

    for (i = 0; i < CD->NumVol; i++)
    {
        for (j = 0; j < CD->NumStc; j++)
        {
            if (CD->chemtype[j].itype == 4)
            {
                if (CD->RelMin == 0)
                {
                    CD->Vcele[i].t_conc[j] =
                        CD->Vcele[i].t_conc[j] * 1000 /
                        CD->chemtype[j].MolarVolume / CD->Vcele[i].porosity;
                    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                    // absolute mineral volume fraction
                }
                if (CD->RelMin == 1)
                {
                    CD->Vcele[i].t_conc[j] =
                        CD->Vcele[i].t_conc[j] * (1 - CD->Vcele[i].porosity +
                        INFTYSMALL) * 1000 / CD->chemtype[j].MolarVolume /
                        CD->Vcele[i].porosity;
                    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                    // relative mineral volume fraction
                    // porosity can be 1.0 so the relative fraction option need a small modification
                }
            }
            if ((CD->chemtype[j].itype == 3) || (CD->chemtype[j].itype == 2))
            {
                CD->Vcele[i].t_conc[j] =
                    CD->Vcele[i].t_conc[j] * (1 - CD->Vcele[i].porosity) * 2650;
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                //     fprintf (stderr, "site concentration is %f mole per porous space L\n", CD->Vcele[i].t_conc[j]);
                // change the unit of CEC (eq/g) into C(ion site) (eq/L porous space), assuming density of solid is always 2650 g/L solid
            }
        }
    }

    CD->SPCFlg = 1;
    if (!CD->RecFlg)
    {
        for (i = 0; i < nelem; i++)
            Speciation(CD, i);
    }
    CD->SPCFlg = 0;

    for (i = nelem; i < nelem * 2; i++)
        for (k = 0; k < CD->NumStc; k++)
        {
            CD->Vcele[i].t_conc[k] = CD->Vcele[i - CD->NumEle].t_conc[k];
            CD->Vcele[i].p_conc[k] = CD->Vcele[i - CD->NumEle].p_conc[k];
            CD->Vcele[i].p_actv[k] = CD->Vcele[i - CD->NumEle].p_actv[k];
            //      fprintf(stderr, " %d %s: %g\n", i, CD->chemtype[k].ChemName, CD->Vcele[i].t_conc[k]);
        }

    for (i = nelem * 2; i < CD->NumVol; i++)
    {
        for (k = 0; k < CD->NumStc; k++)
        {
            if (CD->chemtype[k].itype != 1)
            {
                CD->Vcele[i].t_conc[k] = 1.0E-20;
                CD->Vcele[i].p_conc[k] = 1.0E-20;
                CD->Vcele[i].p_actv[k] = 1.0E-20;
            }
        }
    }

    /*
     * #pragma omp parallel num_threads(2)
     * {
     * int tid, nthreads;
     *
     * tid = omp_get_thread_num();
     * nthreads = omp_get_num_threads();
     *
     * printf("Hello world from thread %3d of %3d\n", tid, nthreads);
     * }
     */

    for (i = 0; i < CD->PIHMFac * 2; i++)
        for (j = CD->PIHMFac * 2; j < CD->NumFac; j++)
        {
            if ((CD->Flux[i].nodeup == CD->Flux[j].nodelo) &&
                (CD->Flux[i].nodelo == CD->Flux[j].nodeup))
            {
                //      fprintf (stderr, " Flux between %d and %d identified \n",
                //CD->Flux[i].nodelo, CD->Flux[i].nodeup);
                CD->Flux[j].BC = CD->Flux[i].BC = 0;
                CD->Flux[j].distance = CD->Flux[i].distance;
                //   fprintf(stderr, " Flux from river to ele corrected as BC %d, distance %f.\n",CD->Flux[j].BC, CD->Flux[j].distance);
            }
        }

    for (i = 0; i < num_conditions; i++)
    {
        free(Condition_vcele[i].t_conc);
        free(Condition_vcele[i].p_conc);
        free(Condition_vcele[i].p_para);
        free(Condition_vcele[i].p_type);
    }
    free(Condition_vcele);

    for (i = 0; i < num_conditions; i++)
        free(chemcon[i]);
    free(chemcon);

    for (i = 0; i < num_conditions + 1; i++)
    {
        for (j = 0; j < CD->NumStc; j++)
            free(con_chem_name[i][j]);
        free(con_chem_name[i]);
    }
    free(con_chem_name);

    free(chemfn);
    free(datafn);
    free(forcfn);
    free(condition_index);

    // 10.23 memory leak
    free(timeinfo);
    free(Global_type.ChemName);
    for (i = 0; i < words_line; i++)
    {
        free(tmpstr[i]);
    }
    free(tmpstr);


    fclose(chemfile);
    fclose(database);
    fclose(prepconc);
}


void
//fluxtrans(realtype t, realtype stepsize, const Model_Data DS, Chem_Data CD, N_Vector VY)  // 09.25 try change the third type
//fluxtrans(int t, int stepsize, const pihm_struct pihm, Chem_Data CD, N_Vector VY)  10.05
fluxtrans(int t, int stepsize, const pihm_struct pihm, Chem_Data CD, N_Vector VY, double *t_duration_transp, double *t_duration_react)  // 10.05 add two timers
{
    /* unit of t and stepsize: min */

    /*
     * swi irreducible water saturation
     * hn  non mobile water height
     * ht  transient zone height
     */

    int             i, j, k = 0;
    struct tm      *timestamp;
    time_t         *rawtime;
    double          timelps, rt_step, temp_rt_step, peclet, invavg, unit_c;
    realtype       *Y;

    // 10.05

    rt_step = stepsize * (double)CD->AvgScl;    // by default, the largest averaging period is per 10 mins. Longer default averaging value will fail

    unit_c = stepsize / UNIT_C;

    Y = NV_DATA_S(VY);

    rawtime = (time_t *)malloc(sizeof(time_t));
    *rawtime = (int)(t * 60);
    timestamp = gmtime(rawtime);
    timelps = t - CD->StartTime;
    // 09.29 test time variables
    //fprintf(stderr, " CD->StartTime = %f \n", CD->StartTime);
    //fprintf(stderr, " int t [min] = %d \n", t);
    //fprintf(stderr, " rawtime [s] = %d \n", *rawtime);
    //fprintf(stderr, " timestamp [date] = %4.4d-%2.2d-%2.2d %2.2d:%2.2d \n", timestamp->tm_year+1900, timestamp->tm_mon+1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
    //fprintf(stderr, " timelps [min] = %f \n", timelps);

    pihm_struct     pihmFlux;
    pihmFlux = pihm;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = nelem; i < 2 * nelem; i++)
    {
        // 11.08
        int             j;
        j = i - nelem;
        CD->Vcele[i].q +=
            MAX(pihmFlux->elem[j].wf.prcp * 86400, 0.0) * CD->Vcele[i].area;
        // 2018.01.08
        //if ((int)timelps % (15) == 0 )
        //{
        //if (i == nelem)
        //{
        //fprintf(stderr, " timelps = %f, elem[i = %d].wf.prcp, infil [m/d] = %6.4g, %6.4g \n", timelps, i+1, pihmFlux->elem[j].wf.prcp * 86400, pihmFlux->elem[j].wf.infil * 86400);  // 09.29 good so far
        //}
        //}
        // 2018.01.08
    }

    /* flux for GW lateral flow */
    for (i = 0; i < nelem; i++)
    {
        CD->Vcele[i].height_tl = Y[i + 2 * nelem];
        CD->Vcele[i + nelem].height_tl = Y[i + nelem];
        //fprintf(stderr, " [timelps = %f, cell = %d ] GW.height_tl, UNSAT.height_tl = %f, %f \n", timelps, (i+1), CD->Vcele[i].height_tl, CD->Vcele[i + nelem].height_tl);  // 09.29 good

        for (j = 0; j < 3; j++)
            if (pihmFlux->meshtbl.nabr[i][j] > 0)
            {
                /* node indicates the index of grid blocks, not nodes at corners */
                //CD->Flux[k].flux += pihmFlux->elem[i].wf.subsurf[j] * 86400;          // 09.29 from [m3/s] to [m3/d]
                CD->Flux[k].flux += 1.0 * pihmFlux->elem[i].wf.subsurf[j] * 86400;  // 02.14 test lateral dilution
                CD->Flux[k].s_area += pihmFlux->elem[i].wf.subareaRT[j];    // 09.29 RT new add, in lat_flow.c
                CD->Flux[k].velocity += pihmFlux->elem[i].wf.subveloRT[j] * 86400;  // 09.29 RT new add, [m/s] to [m/d]
                // 10.27 test
                //fprintf(stderr, " [timelps = %.2f, cell = %d, j = %d, k = %d ] MD->FluxSub, AreaSub, VeloSub = %f, %f, %f \n", timelps, (i+1), (j+1), (k+1), pihmFlux->elem[i].wf.subsurf[j] * 86400, pihmFlux->elem[i].wf.subareaRT[j], pihmFlux->elem[i].wf.subveloRT[j] * 86400);
                //fprintf(stderr, "  -- CD->Flux[k, %d].flux, s_area, velocity = %g, %f, %g \n", (k+1), CD->Flux[k].flux, CD->Flux[k].s_area, CD->Flux[k].velocity);

                /*      if ( fabs(MD->VeloSub[i][j]) < 0.001)
                 * fprintf(stderr, "%d\t %d\t %f\t %f\t %f\t %f\n", CD->Flux[k].nodeup, CD->Flux[k].nodelo, MD->VeloSub[i][j], CD->Vcele[CD->Flux[k].nodeup-1].height_tl,
                 * CD->Vcele[CD->Flux[k].nodelo-1].height_tl, CD->Vcele[CD->Flux[k].nodeup-1].height_tl - CD->Vcele[CD->Flux[k].nodelo-1].height_tl);
                 */
                k++;
            }
    }

    //printf(" GW lateral flux: 0 to k = %d \n", k-1);

    /* flux for UNSAT lateral flow */
    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < 3; j++)
            if (pihmFlux->meshtbl.nabr[i][j] > 0)
            {
                CD->Flux[k].flux = 0.0;
                CD->Flux[k].s_area = 1.0;
                CD->Flux[k].velocity = 0.0;
                k++;
            }
    }

    //printf(" Unsat lateral flux: k = %d \n", k-1);

    /* flux for UNSAT - GW vertical flow */
    for (i = nelem; i < 2 * nelem; i++)
    {
        CD->Flux[k].velocity += pihmFlux->elem[i - nelem].wf.rechg * 86400; // 09.29 [m/s] to [m/d]
        CD->Flux[k].flux +=
            (pihmFlux->elem[i - nelem].wf.rechg * 86400) * CD->Vcele[i].area;
        CD->Flux[k].s_area += CD->Vcele[i].area;
        //fprintf(stderr, " [timelps = %.2f, cell = %d, k = %d ] CD->Flux[k].velocity, flux, s_area = %f, %f, %f \n", timelps, (i+1), (k+1), CD->Flux[k].velocity, CD->Flux[k].flux, CD->Flux[k].s_area);
        k++;
    }

    //printf(" Unsat-GW vertical flux: k = %d \n", k-1);

    for (i = 0; i < nelem; i++)
    {
        CD->Flux[k].velocity += -pihmFlux->elem[i].wf.rechg * 86400;
        CD->Flux[k].flux +=
            (-pihmFlux->elem[i].wf.rechg * 86400) * CD->Vcele[i].area;
        CD->Flux[k].s_area += CD->Vcele[i].area;
        k++;
    }

    //printf(" GW-unsat vertical flux: k = %d \n", k-1);

    /* flux for RIVER flow */
    //CD->riv += MD->FluxRiv[MD->NumRiv - 1][1];  // works for Shale Hills, outlet at 20; but doesn't work for Coal Creek
    //fprintf(stderr, " pihmFlux->riv_outlet = %d \n", pihmFlux->riv_outlet);

    for (i = 0; i < nriver; i++)
    {
        if (pihmFlux->river[i].down < 0)
        {
            CD->riv += pihmFlux->river[i].wf.rivflow[1] * 86400;
        }
    }

    if ((int)timelps % 1440 == 0)
    {
        CD->rivd = CD->riv / 1440;  // 09.29 averaging the sum of 1440 mins for a daily discharge, rivFlx1
        CD->riv = 0;
        //fprintf(stderr, "   Rivd (Qd) is %.3f [m3/d] this day. \n", CD->rivd);
        //printf("   Rivd (Qd) is %.3f [m3/d] this day. \n", CD->rivd);
    }


    // 01.11
    //printf(" k = %d", k);

    for (i = 0; i < nriver; i++)
    {
        //    CD->Flux[k].flux   += MIN(MD->FluxRiv[i][2],0.0);
        CD->Flux[k].flux += pihmFlux->river[i].wf.rivflow[2] * 86400;
        //fprintf(stderr, "    - CD->Flux[k, %d].flux = %.3f [m3/d]. \n", k, CD->Flux[k].flux);
        //if (i == 37 - 1)
        //printf("timelps = %.2f,  riv = 37: rivflx2 (k = %d) = %.3f [m3/d]; ", timelps, k, CD->Flux[k].flux);
        k++;
    }

    //printf(" k = %d", k);

    for (i = 0; i < nriver; i++)
    {
        //    CD->Flux[k].flux   += MIN(MD->FluxRiv[i][3],0.0);
        CD->Flux[k].flux += pihmFlux->river[i].wf.rivflow[3] * 86400;
        //fprintf(stderr, "  [timelps = %.2f] riv[i, %d].wf.rivflow[3] = %.3f [m3/d]. \n", timelps, i, (pihmFlux->riv[i].wf.rivflow[3] * 86400));
        //if (i == 37 - 1)
        //printf(" rivflx3 (k = %d) = %.3f [m3/d]; ", k, CD->Flux[k].flux);
        k++;
    }

    /*
     * for ( i = 0; i <MD->NumRiv; i++){
     * //    CD->Flux[k].flux   += MD->FluxRiv[i][6];
     * k++;
     * }
     * for ( i = 0; i <MD->NumRiv; i++){
     * //  CD->Flux[k].flux   += -MD->FluxRiv[i][6];
     * k++;
     * }
     */

    for (i = 0; i < nriver; i++)
    {
        //    CD->Flux[k].flux   += MD->FluxRiv[i][7];
        //fprintf(stderr, "  [timelps = %.2f], riv = %d, CD->Flux[k = %d].flux = %6.4g [m3/d]. \n", timelps, i+1, k, CD->Flux[k].flux);
        CD->Flux[k].flux +=
            pihmFlux->river[i].wf.rivflow[7] * 86400 +
            pihmFlux->river[i].wf.rivflow[4] * 86400;
        //fprintf(stderr, "      riv[i, %d].wf.rivflow[7],[4] = %6.4g, %6.4g [m3/d] \n", i+1, pihmFlux->river[i].wf.rivflow[7] * 86400, pihmFlux->river[i].wf.rivflow[4] * 86400);
        //if (i == 37 - 1)
        //printf(" rivflx4+7 (k = %d) = %.3f [m3/d]; ", k, CD->Flux[k].flux);
        k++;
    }
    for (i = 0; i < nriver; i++)
    {
        //    CD->Flux[k].flux   += MD->FluxRiv[i][8];
        CD->Flux[k].flux +=
            pihmFlux->river[i].wf.rivflow[8] * 86400 +
            pihmFlux->river[i].wf.rivflow[5] * 86400;
        //fprintf(stderr, "  [timelps = %.2f] riv[i, %d].wf.rivflow[8] = %.3f [m3/d] \n", timelps, i, (pihmFlux->riv[i].wf.rivflow[8] * 86400));
        //if (i == 37 - 1)
        //printf(" rivflx5+8 (k = %d) = %.3f [m3/d]; ", k, CD->Flux[k].flux);
        k++;
    }

    //fprintf(stderr, "\n");
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[k].flux +=
            pihmFlux->river[i].wf.rivflow[9] * 86400 +
            pihmFlux->river[i].wf.rivflow[1] * 86400;
        //fprintf(stderr, "  [timelps = %.2f] riv[i, %d].wf.rivflow[1] = %.3f [m3/d] \n", timelps, i, (pihmFlux->riv[i].wf.rivflow[1] * 86400));
        // 02.21 check
        //fprintf(stderr, "  [timelps = %.2f] riv[i, %d].down[%d] = %.3f [m3/d] \n", timelps, i+1, CD->Flux[k].nodelo - 2 * nelem - nriver, CD->Flux[k].flux);
        //if (i == 37 - 1)
        //printf(" rivflx9+1 (k = %d) = %.3f [m3/d]; ", k, CD->Flux[k].flux);
        k++;
    }



    //printf(" k = %d \n", k);
    for (i = 0; i < nriver; i++)
    {
        // 01.14 change for tributary flux
        //CD->Flux[k].flux += pihmFlux->river[i].wf.rivflow[10] * 86400 + pihmFlux->river[i].wf.rivflow[0] * 86400;  // 01.14 change for tributary flux
        if (pihmFlux->river[i].up[0] > 0) // not tributary starting point
            CD->Flux[k].flux +=
                -(pihmFlux->river[CD->Flux[k].nodelo - 1 - 2 * nelem -
                    nriver].wf.rivflow[9] * 86400 +
                pihmFlux->river[CD->Flux[k].nodelo - 1 - 2 * nelem -
                    nriver].wf.rivflow[1] * 86400);
        else                    // tributary starting point = 0.0 [m3/d]
            CD->Flux[k].flux +=
                pihmFlux->river[i].wf.rivflow[10] * 86400 +
                pihmFlux->river[i].wf.rivflow[0] * 86400;

        if (CD->Flux[k].node_trib > 0)
            CD->Flux[k].flux_trib +=
                -(pihmFlux->river[CD->Flux[k].node_trib - 1 - 2 * nelem -
                    nriver].wf.rivflow[9] * 86400 +
                pihmFlux->river[CD->Flux[k].node_trib - 1 - 2 * nelem -
                    nriver].wf.rivflow[1] * 86400);
        else
            CD->Flux[k].flux_trib += 0.0;

        k++;
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* update the cell volumetrics every averaging cycle */
    // 09.30 5:08PM

    if (((int)timelps - (int)CD->TimLst) == CD->AvgScl * stepsize)  // 10.24 fix uninitialised value(s) on "CD->TimLst"
    {
        // 09.30
        /*fprintf(stderr,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\", rawtime = %d [s] \n", \
         * timestamp->tm_year+1900, timestamp->tm_mon+1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min, (int) * rawtime); */

        // update the concentration in precipitation here.
        if (CD->PrpFlg == 2)
        {
            while (CD->TSD_prepconc->TS[CD->TSD_prepconc->iCounter][0] <
                (realtype)*rawtime / (24 * 60 * 60))
                CD->TSD_prepconc->iCounter++;
            CD->TSD_prepconc->iCounter--;
            // fprintf(stderr, " %d %6.4f time %6.4f\n", CD->TSD_prepconc->iCounter, CD->TSD_prepconc->TS[CD->TSD_prepconc->iCounter][0], (realtype) * rawtime/(24*60*60));
            for (i = 0; i < CD->TSD_prepconc->index; i++)
            {
                if (CD->prepconcindex[i] > 0)
                {
                    j = CD->prepconcindex[i] - 1;
                    if (CD->Precipitation.t_conc[j] !=
                        CD->TSD_prepconc->TS[CD->TSD_prepconc->iCounter][i + 1])
                    {
                        CD->Precipitation.t_conc[j] =
                            CD->TSD_prepconc->TS[CD->TSD_prepconc->iCounter][i +
                            1];
                        fprintf(stderr,
                            "  %s in precipitation is changed to %6.4g\n",
                            CD->chemtype[j].ChemName,
                            CD->Precipitation.t_conc[j]);
                    }
                }
            }
        }

        // 1.11
        //printf(" NumFac = %d, AvgScl = %d ", CD->NumFac, CD->AvgScl);

        invavg = stepsize / (double)CD->AvgScl;

        for (i = 0; i < 2 * CD->NumEle; i++)
        {
            CD->Vcele[i].height_tl = 0.0;
        }

        for (k = 0; k < CD->NumFac; k++)
        {
            CD->Flux[k].velocity *= invavg;
            CD->Flux[k].flux *= invavg;
            CD->Flux[k].s_area *= invavg;

            if (CD->Flux[k].s_area > 1.0E-4)
            {
                CD->Flux[k].velocity = CD->Flux[k].flux / CD->Flux[k].s_area;
            }
            else
            {
                CD->Flux[k].velocity = 1.0E-10;
            }

            // 1.11
            /*
             * printf (" k = 10012, CD->Flux[k].nodeup, nodelo = %d, %d \n", CD->Flux[10012-1].nodeup, CD->Flux[10012-1].nodelo, CD->Flux[10012-1].nodeuu, CD->Flux[10012-1].nodell);
             * printf (" k = 10013, CD->Flux[k].nodeup, nodelo = %d, %d \n", CD->Flux[10013-1].nodeup, CD->Flux[10013-1].nodelo, CD->Flux[10013-1].nodeuu, CD->Flux[10013-1].nodell);
             * printf ("   = 10122, CD->Flux[k].nodeup, nodelo = %d, %d \n", CD->Flux[10122-1].nodeup, CD->Flux[10122-1].nodelo, CD->Flux[10122-1].nodeuu, CD->Flux[10122-1].nodell);
             * printf ("   = 10232, CD->Flux[k].nodeup, nodelo = %d, %d \n", CD->Flux[10232-1].nodeup, CD->Flux[10232-1].nodelo, CD->Flux[10232-1].nodeuu, CD->Flux[10232-1].nodell);
             * printf ("   = 10233, CD->Flux[k].nodeup, nodelo = %d, %d \n", CD->Flux[10233-1].nodeup, CD->Flux[10233-1].nodelo, CD->Flux[10233-1].nodeuu, CD->Flux[10233-1].nodell);
             */

            CD->Vcele[CD->Flux[k].nodeup - 1].height_tl +=
                CD->Flux[k].velocity * CD->Flux[k].distuu;
            //CD->Vcele[CD->Flux[k].nodeup - 1 + CD->NumEle].height_tl += CD->Flux[k].velocity * CD->Flux[k].distll;   // 10.23 comment-out, fix invalid read/write
            //    CD->Vcele[CD->Flux[k].nodelo-1].height_tl += fabs(CD->Flux[k].velocity);

            /*
             * if (( fabs(CD->Flux[k].velocity) < CD->CnntVelo) && (k < CD->NumLateral))
             * {
             * CD->Vcele[CD->Flux[k].nodeup-1].height_tl = 0.0;
             * CD->Vcele[CD->Flux[k].nodelo-1].height_tl = 0.0;
             */
            /*fprintf(stderr, "%d\t %d\t %f\t %f\t %f\t %f\n", CD->Flux[k].nodeup, CD->Flux[k].nodelo, CD->Flux[k].velocity, CD->Vcele[CD->Flux[k].nodeup-1].height_tl,
             * CD->Vcele[CD->Flux[k].nodelo-1].height_tl, CD->Vcele[CD->Flux[k].nodeup-1].height_tl - CD->Vcele[CD->Flux[k].nodelo-1].height_tl); */
            /*
             * kk ++ ;
             * }
             */
            // velocity is calculated according to the flux_avg and the area_avg.
        }                       // for gw cells, contact area is needed for dispersion;
        //      fprintf(stderr, " Number of faces disconnected %d, ratio %f\n", kk, (double)kk/(double)CD->NumLateral);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < 2 * CD->NumEle; i++)
        {
            CD->Vcele[i].height_tl = CD->Vcele[i].height_tl / 3.0;
        }


#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < CD->PIHMFac * 2; i++)
        {
            int             j;  // 11.09 fix bug
            for (j = CD->PIHMFac * 2; j < CD->NumFac; j++)
            {
                if ((CD->Flux[i].nodeup == CD->Flux[j].nodelo) &&
                    (CD->Flux[i].nodelo == CD->Flux[j].nodeup))
                {
                    CD->Flux[j].s_area = CD->Flux[i].s_area;
                    CD->Flux[j].velocity = -CD->Flux[i].velocity;
                    //   fprintf(stderr, " Flux from river to ele corrected as BC %d, distance %f.\n",CD->Flux[j].BC, CD->Flux[j].distance);
                }
            }
        }


#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            CD->Vcele[i].height_o = CD->Vcele[i].height_t;
            CD->Vcele[i].height_t = MAX(Y[i + 2 * nelem], 1.0E-5);
            CD->Vcele[i].height_int = CD->Vcele[i].height_t;
            CD->Vcele[i].height_sp =
                (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
            CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
            CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        }


#ifdef _OPENMP
#pragma omp parallel for
#endif
        /* update the unsaturated zone (vadoze) */
        for (i = nelem; i < 2 * nelem; i++)
        {
            // 11.08
            int             j;
            j = i - nelem;
            CD->Vcele[i].height_o = CD->Vcele[i].height_t;
            CD->Vcele[i].height_t =
                MAX(((Y[i] * (pihmFlux->elem[j].soil.smcmax -
                            pihmFlux->elem[j].soil.smcmin) +
                        (CD->Vcele[i].height_v - CD->Vcele[i -
                                nelem].height_t) *
                        pihmFlux->elem[j].soil.smcmin) /
                    pihmFlux->elem[j].soil.smcmax), 1.0E-5);
            //    fprintf(stderr, " PIHM Hu: %6.4f, RT Hu*: %6.4f\n", MD->DummyY[i], CD->Vcele[i].height_o);
            // 09.30
            //fprintf(stderr, "  (unsat) CD->Vcele[i, %d].height_o & height_t = %6.4f, %6.4f \n", i, CD->Vcele[i].height_o, CD->Vcele[i].height_t);
            CD->Vcele[i].height_int = CD->Vcele[i].height_t;
            CD->Vcele[i].height_sp =
                (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
            CD->Vcele[i].sat_o = CD->Vcele[i].sat;
            CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
            CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
            //      CD->Vcele[i].sat      = CD->Vcele[i].height_t/ (CD->Vcele[i].height_v - CD->Vcele[i-MD->NumEle].height_t);
            CD->Vcele[i].sat =
                CD->Vcele[i].height_t / (CD->Vcele[i].height_v - CD->Vcele[i -
                    nelem].height_t);
            CD->Vcele[i].q *= invavg;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        /* update river cells */
        for (i = 2 * nelem; i < 2 * nelem + nriver; i++)
        {
            // 11.08
            int             j;
            j = i - 2 * nelem;
            CD->Vcele[i].height_o = CD->Vcele[i].height_t;
            CD->Vcele[i].height_t = MAX(Y[i + nelem], 1.0E-5);
            CD->Vcele[i].height_int = CD->Vcele[i].height_t;
            CD->Vcele[i].height_sp =
                (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
            CD->Vcele[i].area = pihmFlux->river[j].topo.area;
            // 09.30
            //fprintf(stderr, "  (RIV) CD->Vcele[i, %d].area = %6.4f \n", i, CD->Vcele[i].area);
            CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
            CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        }


#ifdef _OPENMP
#pragma omp parallel for
#endif
        /* update EBR cells */
        for (i = 2 * nelem + nriver; i < 2 * (nelem + nriver); i++)
        {
            // 11.08
            int             j;
            j = i - 2 * nelem - nriver;
            CD->Vcele[i].height_o = CD->Vcele[i].height_t;
            CD->Vcele[i].height_t =
                MAX(Y[i + nelem], 1.0E-5) + MAX(Y[i + nelem - nriver],
                1.0E-5) / CD->Vcele[i].porosity;
            // 09.30
            //fprintf(stderr, "  (EBR) CD->Vcele[i, %d].height_o & height_t = %6.4f, %6.4f \n", i, CD->Vcele[i].height_o, CD->Vcele[i].height_t);
            CD->Vcele[i].height_int = CD->Vcele[i].height_t;
            CD->Vcele[i].height_sp =
                (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
            CD->Vcele[i].area = pihmFlux->river[j].topo.area;
            // 09.30
            //fprintf(stderr, "  (EBR) CD->Vcele[i, %d].area = %6.4f \n", i, CD->Vcele[i].area);
            CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
            CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
        }

        // 09.30 8:25PM
        Monitor(t, stepsize * (double)CD->AvgScl, pihm, CD);

        // 09.30 11:19PM
        for (k = 0; k < CD->NumStc; k++)
        {
            CD->Vcele[CD->NumVol - 1].t_conc[k] =
                CD->Precipitation.t_conc[k] * CD->Condensation;
            CD->Vcele[CD->NumVol - 1].p_conc[k] =
                CD->Precipitation.t_conc[k] * CD->Condensation;
        }

        /*
         * for (i = CD->PIHMFac; i < CD->PIHMFac * 2; i++)
         * for (j = CD->PIHMFac; j < CD->PIHMFac * 2; j++)
         * {
         * if ((CD->Flux[i].nodeup == CD->Flux[j].nodelo) && (CD->Flux[i].nodelo == CD->Flux[j].nodeup))
         * {
         * if (fabs(CD->Flux[i].flux + CD->Flux[j].flux) > 1.0E-5)
         * {
         * fprintf(stderr, " Dismatch %d %d %f %f\n", CD->Flux[i].nodeup, CD->Flux[i].nodelo, CD->Flux[i].flux, CD->Flux[j].flux);
         * fprintf(stderr, " Corresponding %d %d %f %f\n", CD->Flux[i - CD->PIHMFac].nodeup, CD->Flux[i - CD->PIHMFac].nodelo, CD->Flux[i - CD->PIHMFac].flux, CD->Flux[j - CD->PIHMFac].flux);
         * }
         * }
         * }
         */

        for (i = 0; i < CD->NumVol; i++)
            CD->Vcele[i].rt_step = 0.0;

        for (i = 0; i < CD->NumDis; i++)
        {                       // 02.03 fix left '{'
            for (j = 0; j < CD->NumSpc; j++)
            {
                peclet =
                    fabs(CD->Flux[i].velocity * CD->Flux[i].distance /
                    (CD->chemtype[j].DiffCoe +
                        CD->chemtype[j].DispCoe * CD->Flux[i].velocity));
                peclet = MAX(peclet, 1.0E-8);
            }                   // 02.03 fix right '}'

            if (i < CD->NumDis)
            {
                temp_rt_step =
                    fabs(CD->Flux[i].flux / CD->Vcele[CD->Flux[i].nodeup -
                        1].vol) * (1 + peclet) / peclet;
            }
            else
                temp_rt_step =
                    fabs(CD->Flux[i].flux / CD->Vcele[CD->Flux[i].nodeup -
                        1].vol);

            CD->Vcele[CD->Flux[i].nodeup - 1].rt_step += temp_rt_step;
        }

        /*
         * for ( i = 0 ; i < CD->NumEle; i ++)
         * if ( MD->ElePrep[i] > 1.0E-8)
         * fprintf(stderr, " Ratio %f from %f/%f or %f %f %f\n", MD->EleET[i][2]/MD->EleNetPrep[i], MD->EleET[i][2], MD->EleNetPrep[i], MD->ElePrep[i], MD->EleET[i][1], MD->EleViR[i]);
         */

        //    FILE * rt_dist = fopen("logfile/rtdist.out", "w");

        k = 0;                  // used to count the number of slow cells.
        invavg = 0.0;
        for (i = 0; i < CD->NumOsv; i++)
        {
            CD->Vcele[i].rt_step = 0.6 * UNIT_C / CD->Vcele[i].rt_step;
            //      fprintf(rt_dist, " %d\t %f\n", i+1, CD->Vcele[i].rt_step);
            if (rt_step > CD->Vcele[i].rt_step)
                rt_step = CD->Vcele[i].rt_step;

            if (CD->Vcele[i].rt_step >= (double)CD->AvgScl)
                CD->Vcele[i].rt_step = (double)CD->AvgScl;
            else
            {
                k++;
            }
            invavg += CD->Vcele[i].rt_step;
        }

        invavg = invavg / (double)CD->NumOsv;
        //    fprintf(rt_dist, " Number of slow cell is %d\n ",k);
        //    fclose(rt_dist);

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        // RT step control begins

        //for (i = 0; i < CD->NumEle; i++)
        //printf(" t = %f, [Cell,%d]: height_v = %f \n", timelps, i+1, CD->Vcele[i].height_v);  // 11.28 check fd

        if (CD->TimLst >= CD->Delay)
        {
            //fprintf(stderr, "\n Calling RT from %.2f to %.2f, mean timestep = %.3f [min]. \n", CD->TimLst, stepsize * (double)CD->AvgScl + CD->TimLst, invavg);
            rt_step = stepsize * (double)CD->AvgScl;

            // AdptTime(CD, CD->TimLst, stepsize * (double)CD->AvgScl , stepsize * (double)CD->AvgScl, &rt_step, MD->NumEle * 2);
            // AdptTime (CD, CD->TimLst, rt_step, stepsize * (double) CD->AvgScl, &rt_step, MD->NumEle * 2);  // 08.22 comment out
            // AdptTime(CD, CD->TimLst, rt_step, stepsize * (double)CD->AvgScl, &rt_step, MD->NumEle * 2, DS); //08.22 mdata (DC is converted to MD)
            // AdptTime(CD, CD->TimLst, rt_step, stepsize * (double)CD->AvgScl, &rt_step, nelem * 2, pihm);  10.05
            AdptTime(CD, CD->TimLst, rt_step, stepsize * (double)CD->AvgScl, &rt_step, nelem * 2, pihm, t_duration_transp, t_duration_react);   // 10.05 add two timers

            for (i = 0; i < CD->NumEle; i++)
                for (j = 0; j < CD->NumStc; j++)
                {
                    if (CD->chemtype[j].itype == 4)
                    {
                        CD->Vcele[i].t_conc[j] =
                            (CD->Vcele[i].t_conc[j] * CD->Vcele[i].height_t +
                            CD->Vcele[i +
                                CD->NumEle].t_conc[j] * (CD->Vcele[i].height_v -
                                CD->Vcele[i].height_t)) / CD->Vcele[i].height_v;
                        CD->Vcele[i + CD->NumEle].t_conc[j] =
                            CD->Vcele[i].t_conc[j];
                        CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                        CD->Vcele[i + CD->NumEle].p_conc[j] =
                            CD->Vcele[i].t_conc[j];
                        // Averaging mineral concentration to ensure mass conservation !!
                    }
                }

            for (i = 0; i < CD->NumOsv; i++)
            {
                if (fabs(CD->Vcele[i].height_t - CD->Vcele[i].height_int) >
                    1.0E-6)
                    fprintf(stderr, "%d %6.4f\t%6.4f\n", i,
                        CD->Vcele[i].height_t, CD->Vcele[i].height_int);
                assert(fabs(CD->Vcele[i].height_t - CD->Vcele[i].height_int) <
                    1.0E-6);
                if (CD->Vcele[i].illness >= 20)
                {
                    for (j = 0; j < CD->NumStc; j++)
                        CD->Vcele[i].t_conc[j] = 1.0E-10;
                    fprintf(stderr,
                        " Cell %d isolated due to proneness to err!\n",
                        CD->Vcele[i].index);
                }
            }
            // to make sure intrapolation worked well.
        }
        /* RT step control ends */

        //chem_updater(CD, pihm); // 10.01
        CD->TimLst = timelps;

        for (i = nelem; i < 2 * nelem; i++)
        {
            j = i - nelem;
            CD->Vcele[i].q = 0.0;
        }

        // Reset fluxes for next averaging stage
        for (k = 0; k < CD->NumDis; k++)
        {
            CD->Flux[k].velocity = 0.0;
            CD->Flux[k].flux = 0.0;
            CD->Flux[k].flux_trib = 0.0;    // 01.14 do not forget
            CD->Flux[k].s_area = 0.0;
        }                       // for gw cells, contact area is needed for dispersion;

        for (k = 0; k < CD->NumFac; k++)
        {
            CD->Flux[k].flux = 0.0;
            CD->Flux[k].flux_trib = 0.0;    // 01.14 do not forget
            CD->Flux[k].velocity = 0.0;
        }                       // for riv cells, contact area is not needed;


        free(rawtime);
    }
}

void InitialChemFile(char *outputdir, char *filename, int NumBTC, int *BTC_loc)
{
    FILE           *Cfile[20];
    char           *cfn[20];
    int             i;

    // 09.28
    //cfn[0] = (char *)malloc((strlen(filename) + 20) * sizeof(char));
    cfn[0] =
        (char *)malloc((strlen(outputdir) + strlen(filename) +
            100) * sizeof(char));
    //sprintf(cfn[0], "output/%s.conc", filename);
    sprintf(cfn[0], "%sRT_%s.conc", outputdir, filename);

    //cfn[1] = (char *)malloc((strlen(filename) + 20) * sizeof(char));
    cfn[1] =
        (char *)malloc((strlen(outputdir) + strlen(filename) +
            100) * sizeof(char));
    //sprintf(cfn[1], "output/%s.gwt", filename);
    sprintf(cfn[1], "%sRT_%s.gwt", outputdir, filename);

    //cfn[2] = (char *)malloc((strlen(filename) + 20) * sizeof(char));
    cfn[2] =
        (char *)malloc((strlen(outputdir) + strlen(filename) +
            100) * sizeof(char));
    //sprintf(cfn[2], "output/%s.btcv", filename);
    sprintf(cfn[2], "%sRT_%s.btcv", outputdir, filename);

    for (i = 3; i < 3 + NumBTC; i++)
    {
        cfn[i] =
            (char *)malloc((strlen(outputdir) + strlen(filename) +
                20) * sizeof(char));
        //sprintf(cfn[i], "output/%s%d.btcv", filename, BTC_loc[i - 3] + 1);
        sprintf(cfn[i], "%sRT_%s%d.btcv", outputdir, filename,
            BTC_loc[i - 3] + 1);
    }

    cfn[3 + NumBTC] =
        (char *)malloc((strlen(outputdir) + strlen(filename) +
            20) * sizeof(char));
    sprintf(cfn[3 + NumBTC], "%sRT_%s.vol", outputdir, filename);

    for (i = 0; i <= 3 + NumBTC; i++)
    {
        Cfile[i] = fopen(cfn[i], "w");
        free(cfn[i]);
        fclose(Cfile[i]);
    }
}

void PrintChem(char *outputdir, char *filename, Chem_Data CD, int t)    // 10.01
{
    FILE           *Cfile[20];
    char           *cfn[20];
    int             i, j, k;
    struct tm      *timestamp;
    time_t         *rawtime;
    double          timelps, tmpconc;

    rawtime = (time_t *)malloc(sizeof(time_t));
    *rawtime = (int)(t * 60);
    timestamp = gmtime(rawtime);
    timelps = t - CD->StartTime;
    tmpconc = 0.0;

    //if ((int)timelps % (720) == 0)
    if ((int)timelps % (60) == 0)   // 11.07 update every 1 hour instead every 12 hour
    {
        CD->SPCFlg = 0;

        if (!CD->RecFlg)
        {
            for (i = 0; i < CD->NumStc; i++)
            {
                for (j = CD->NumEle * 2; j < CD->NumOsv; j++)
                    if (CD->chemtype[i].itype == 4)
                        CD->Vcele[j].p_conc[i] = CD->Vcele[j].t_conc[i];
                    else
                        CD->Vcele[j].p_conc[i] =
                            fabs(CD->Vcele[j].t_conc[i] * 0.1);
            }
        }

        if (!CD->RecFlg)
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = CD->NumEle * 2 + CD->NumRiv; i < CD->NumOsv; i++)
                Speciation(CD, i);
        else
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = 0; i < CD->NumOsv; i++)
                Speciation(CD, i);
    }

    if ((int)timelps % (CD->OutItv * 60) == 0)
    {
        cfn[0] =
            (char *)malloc((strlen(outputdir) + strlen(filename) +
                100) * sizeof(char));
        sprintf(cfn[0], "%sRT_%s.conc", outputdir, filename);
        Cfile[0] = fopen(cfn[0], "a+");

        cfn[1] =
            (char *)malloc((strlen(outputdir) + strlen(filename) +
                100) * sizeof(char));
        sprintf(cfn[1], "%sRT_%s.btcv", outputdir, filename);
        Cfile[1] = fopen(cfn[1], "a+");

        cfn[2] =
            (char *)malloc((strlen(outputdir) + strlen(filename) +
                100) * sizeof(char));
        sprintf(cfn[2], "%sRT_%s.gwt", outputdir, filename);
        Cfile[2] = fopen(cfn[2], "a+");

        for (i = 3; i < 3 + CD->NumBTC; i++)
        {
            cfn[i] =
                (char *)malloc((strlen(outputdir) + strlen(filename) +
                    100) * sizeof(char));
            sprintf(cfn[i], "%sRT_%s%d.btcv", outputdir, filename,
                CD->BTC_loc[i - 3] + 1);
            Cfile[i] = fopen(cfn[i], "a+");
            if (Cfile[i] == NULL)
                fprintf(stderr, " Output BTC not found \n");
        }
        cfn[3 + CD->NumBTC] =
            (char *)malloc((strlen(outputdir) + strlen(filename) +
                100) * sizeof(char));
        sprintf(cfn[3 + CD->NumBTC], "%sRT_%s.vol", outputdir, filename);
        Cfile[3 + CD->NumBTC] = fopen(cfn[3 + CD->NumBTC], "a+");

        fprintf(Cfile[0], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\n",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
        fprintf(Cfile[0], "Cell\t");
        for (i = 0; i < CD->NumStc; i++)
            fprintf(Cfile[0], "%6s\t", CD->chemtype[i].ChemName);
        for (i = 0; i < CD->NumSsc; i++)
            fprintf(Cfile[0], "%6s\t", CD->chemtype[i + CD->NumStc].ChemName);
        fprintf(Cfile[0], "\n");

        for (i = 0; i < CD->NumEle * 2; i++)
        {
            fprintf(Cfile[0], "%d\t", i + 1);
            for (j = 0; j < CD->NumStc; j++)
                fprintf(Cfile[0], "%12.8f\t", log10(CD->Vcele[i].p_conc[j]));
            for (j = 0; j < CD->NumSsc; j++)
                fprintf(Cfile[0], "%12.8f\t", log10(CD->Vcele[i].s_conc[j]));

            fprintf(Cfile[0], "\n");
        }

        for (i = CD->NumEle * 2; i < CD->NumOsv; i++)
        {
            fprintf(Cfile[0], "%d\t", i + 1);
            for (j = 0; j < CD->NumStc; j++)
                fprintf(Cfile[0], "%12.8f\t", log10(CD->Vcele[i].p_conc[j]));
            for (j = 0; j < CD->NumSsc; j++)
                fprintf(Cfile[0], "%12.8f\t", log10(CD->Vcele[i].s_conc[j]));

            fprintf(Cfile[0], "\n");
        }

        // Output the breakthrough curves "stream chemistry"
        // River is not updated in reaction stage, so a speciation is required before output

        if (t == CD->StartTime + CD->OutItv * 60)
        {
            fprintf(Cfile[1], "\t\t\t");
            for (i = 0; i < CD->NumStc; i++)
                fprintf(Cfile[1], "%6s\t", CD->chemtype[i].ChemName);
            for (i = 0; i < CD->NumSsc; i++)
                fprintf(Cfile[1], "%6s\t",
                    CD->chemtype[i + CD->NumStc].ChemName);
            fprintf(Cfile[1], "\n");
        }
        /* print the header of file if first time entered */
        fprintf(Cfile[1], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

        for (j = 0; j < CD->NumStc; j++)
        {
            tmpconc = 0.0;
            for (i = 0; i < 2 * CD->NumEle; i++)
                tmpconc += CD->Vcele[i].p_conc[j];
            tmpconc = tmpconc / (double)(CD->NumEle) * 0.5;
            fprintf(Cfile[1], "%12.8f\t", log10(tmpconc));
        }

        for (j = 0; j < CD->NumSsc; j++)
        {
            tmpconc = 0.0;
            for (i = 0; i < 2 * CD->NumEle; i++)
                tmpconc += CD->Vcele[i].s_conc[j];
            tmpconc = tmpconc / (double)(CD->NumEle) * 0.5;
            fprintf(Cfile[1], "%12.8f\t", log10(tmpconc));
        }
        fprintf(Cfile[1], "\n");

        /*
         * if ( t == CD->StartTime + CD->OutItv * 60){
         * fprintf(Cfile[2], "\t\t\t");
         * for ( i = 0 ; i < CD->NumStc; i++)
         * fprintf(Cfile[2], "%6s\t", CD->chemtype[i].ChemName);
         * for ( i = 0 ; i < CD->NumSsc; i++)
         * fprintf(Cfile[2], "%6s\t", CD->chemtype[i+CD->NumStc].ChemName);
         * fprintf(Cfile[2],"\n");
         * }
         * // print the header of file if first time entered
         * fprintf(Cfile[2],"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
         * for (j = 0; j < CD->NumStc; j++)
         * fprintf(Cfile[2], "%12.8f\t",log10(CD->Vcele[15].p_conc[j]));
         * for (j = 0; j < CD->NumSsc; j++)
         * fprintf(Cfile[2], "%12.8f\t",log10(CD->Vcele[15].s_conc[j]));
         * fprintf(Cfile[2], "\n");
         */

        fprintf(Cfile[2], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

        for (j = 0; j < CD->NumEle * 2; j++)
        {
            fprintf(Cfile[2], "%6.4f\t", CD->Vcele[j].height_t);
        }
        fprintf(Cfile[2], "\n");


        fprintf(Cfile[3 + CD->NumBTC], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
        for (j = 0; j < 2 * CD->NumEle; j++)
        {
            fprintf(Cfile[3 + CD->NumBTC], "%6.4f\t", CD->Vcele[j].height_tl);
        }
        fprintf(Cfile[3 + CD->NumBTC], "\n");

        for (k = 3; k < 3 + CD->NumBTC; k++)
        {
            if (t == CD->StartTime + CD->OutItv * 60)
            {
                fprintf(Cfile[k], "\t\t\t");
                for (i = 0; i < CD->NumStc; i++)
                    fprintf(Cfile[k], "%6s\t", CD->chemtype[i].ChemName);
                for (i = 0; i < CD->NumSsc; i++)
                    fprintf(Cfile[k], "%6s\t",
                        CD->chemtype[i + CD->NumStc].ChemName);
                fprintf(Cfile[k], "\n");
            }
            /* print the header of file if first time entered */
            fprintf(Cfile[k], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
            for (j = 0; j < CD->NumStc; j++)
            {

                if ((CD->BTC_loc[k - 3] >= CD->pumps[0].Pump_Location - 1) &&
                    (j == CD->pumps[0].Position_Species))
                {
                    fprintf(Cfile[k], "%12.8f\t",
                        log10((CD->Vcele[CD->BTC_loc[k -
                                        3]].p_conc[j] * CD->rivd +
                                CD->pumps[0].Injection_conc *
                                CD->pumps[0].flow_rate) / (CD->rivd +
                                CD->pumps[0].flow_rate)));
                    // 10.02 comment out
                    /*
                     * fprintf(stderr, "i %d infd %6.4g, rivd %6.4g, tconc %6.4g, nconc %6.4g\n", CD->BTC_loc[k - 3] + 1, CD->pumps[0].flow_rate / 84170, \
                     * CD->rivd / 84170, CD->Vcele[CD->BTC_loc[k - 3]].p_conc[j], ((CD->Vcele[CD->BTC_loc[k - 3]].p_conc[j] * CD->rivd + \
                     * CD->pumps[0].Injection_conc * CD->pumps[0].flow_rate) / (CD->rivd + CD->pumps[0].flow_rate)));
                     */
                }
                else if (CD->BTC_loc[k - 3] > 2 * CD->NumEle)
                    fprintf(Cfile[k], "%12.8f\t",
                        log10(CD->Vcele[CD->BTC_loc[k - 3]].p_conc[j]));
                else
                    fprintf(Cfile[k], "%12.8f\t",
                        log10(CD->Vcele[CD->BTC_loc[k - 3]].p_conc[j]));
            }
            for (j = 0; j < CD->NumSsc; j++)
                fprintf(Cfile[k], "%12.8f\t",
                    log10(CD->Vcele[CD->BTC_loc[k - 3]].s_conc[j]));
            fprintf(Cfile[k], "\n");
        }

        for (i = 0; i <= 3 + CD->NumBTC; i++)
        {
            free(cfn[i]);
            fclose(Cfile[i]);
        }

    }
    free(rawtime);
}





void
//AdptTime(Chem_Data CD, realtype timelps, double rt_step, double hydro_step, double *start_step, int num_blocks, const void * Model_Data) // 08.22 add model data (from fluxtrans, MD)
//AdptTime(Chem_Data CD, realtype timelps, double rt_step, double hydro_step, double *start_step, int num_blocks, const pihm_struct pihm)  10.05
AdptTime(Chem_Data CD, realtype timelps, double rt_step, double hydro_step, double *start_step, int num_blocks, const pihm_struct pihm, double *t_duration_transp, double *t_duration_react)    // 10.05 add two timer
{
    //double stepsize, org_time, step_rst, end_time, substep;
    double          stepsize, org_time, step_rst, end_time; // 10.05 'substep' location changed
    double          timer1, timer2;
    //int i, j, k, m, nr_tmp, nr_max, int_flg, tot_nr;
    int             i, j, k, m, nr_max, int_flg, tot_nr;    // 10.05 'nr_temp' location changed
    time_t          t_start_transp, t_end_transp, t_start_react, t_end_react;

    stepsize = *start_step;
    org_time = timelps;
    step_rst = *start_step;
    tot_nr = 0;

    // 10.01 test purpose
    //printf("  Doing AdptTime() ... \n");

    // 10.05 timing
    t_start_transp = time(NULL);

    if (rt_step >= hydro_step)
        int_flg = 0;
    else
    {
        int_flg = 1;
        fprintf(stderr, " Sub time step intrapolation performed. \n");
    }

    end_time = org_time + hydro_step;

    // simple check to determine whether or not to intrapolate the gw height;

    while (timelps < end_time)
    {
        nr_max = 5;
        //nr_tmp = 1;  // 10.05 location change

        if (stepsize > end_time - timelps)
        {
            // before adjusting, write the current timestep to file
            step_rst = stepsize;
            stepsize = end_time - timelps;
        }

        //    fprintf(stderr, "  Transport From %10.2f to %10.2f", timelps,  timelps +stepsize);

        if (int_flg)
        {
            // do intrapolation. note that height_int always store the end time height.

            for (i = 0; i < CD->NumOsv; i++)
            {
                CD->Vcele[i].height_t =
                    CD->Vcele[i].height_o + CD->Vcele[i].height_sp * stepsize;
                CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
            }
        }

        timer1 = timer();

        for (i = 0; i < num_blocks; i++)
            for (j = 0; j < CD->NumSpc; j++)
            {
                if (CD->chemtype[j].mtype == 2)
                {
                    for (k = 0; k < CD->NumSsc; k++)
                        if ((CD->Totalconc[j][k + CD->NumStc] != 0) &&
                            (CD->chemtype[k + CD->NumStc].itype != 1))
                        {
                            //   if ( i == 58) fprintf(stderr, " ##%s Tconc %f\t %s Sconc %f\t Mconc %f\n", CD->chemtype[j].ChemName,  CD->Vcele[i].t_conc[j], CD->chemtype[k+ CD->NumStc].ChemName, CD->Vcele[i].s_conc[k],CD->Vcele[i].t_conc[j] - CD->Totalconc[j][k + CD->NumStc] * CD->Vcele[i].s_conc[k]);
                            CD->Vcele[i].t_conc[j] =
                                CD->Vcele[i].t_conc[j] - CD->Totalconc[j][k +
                                CD->NumStc] * CD->Vcele[i].s_conc[k] *
                                CD->TimRiv;
                            //   if ( CD->Vcele[i].t_conc[j] < 0.0) fprintf( stderr, " Mixed mobility error! %s %f\n", CD->chemtype[j].ChemName, CD->Vcele[i].t_conc[j]);
                        }
                }
            }

        /*
         * for (i = CD->NumEle * 2; i < CD->NumEle * 2 + CD->NumRiv; i++)
         * {
         * for (j = 0; j < CD->NumStc; j++)
         * {
         * //       CD->Vcele[i].t_conc[j] = CD->Precipitation.t_conc[j] * CD->Condensation ;
         * CD->Vcele[i].t_conc[j] = 1.0E-20;
         * }
         * }
         */

        OS3D(timelps, stepsize, CD);

        // 10.05 test
        //printf(" !! Pssing OS3D? \n\n");

        for (i = 0; i < num_blocks; i++)
            for (j = 0; j < CD->NumSpc; j++)
            {
                if (CD->chemtype[j].mtype == 2)
                {
                    for (k = 0; k < CD->NumSsc; k++)
                        if ((CD->Totalconc[j][k + CD->NumStc] != 0) &&
                            (CD->chemtype[k + CD->NumStc].itype != 1))
                        {
                            //    if ( i == 1) fprintf(stderr, " ##%s Tconc %f\t %s Sconc %f\t Mconc %f\n", CD->chemtype[j].ChemName,  CD->Vcele[i].t_conc[j], CD->chemtype[k+ CD->NumStc].ChemName, CD->Vcele[i].s_conc[k],CD->Vcele[i].t_conc[j] + CD->Totalconc[j][k + CD->NumStc] * CD->Vcele[i].s_conc[k]);
                            CD->Vcele[i].t_conc[j] =
                                CD->Vcele[i].t_conc[j] + CD->Totalconc[j][k +
                                CD->NumStc] * CD->Vcele[i].s_conc[k] *
                                CD->TimRiv;
                            //    if ( CD->Vcele[i].t_conc[j] < 0.0) fprintf( stderr, " Mixed mobility error! %s %f\n", CD->chemtype[j].ChemName, CD->Vcele[i].t_conc[j]);
                            // total concentration except for adsorptions have been transported and adjusted by the volume. For example, if no transport but volume increased by rain, the concentration need be decreased. However, the adsorption part has not been treated yet, so they need be adjusted by the volume change.
                            // the porosity is not changed during the period, so the ratio between pore space before and after OS3D is the same ratio between volume of porous media before and after OS3D.
                        }
                }
            }

        // 10.05
        t_end_transp = time(NULL);
        *t_duration_transp += (t_end_transp - t_start_transp);

        /*
         * for (i = 0; i < CD->NumPUMP; i++)
         * {
         * //j = CD->pumps[i].Pump_Location;
         * for ( j = 0 ; j < CD->NumEle; j ++)
         * //      fprintf(stderr, " T_conc binj %6.4g", CD->Vcele[j].t_conc[CD->pumps[i].Position_Species]);
         * {
         * CD->Vcele[j].t_conc[CD->pumps[i].Position_Species] +=
         * CD->pumps[i].flow_rate * CD->pumps[i].Injection_conc * stepsize * CD->Vcele[i].area / (CD->Vcele[j].porosity *
         * CD->Vcele[j].vol);
         * fprintf (stderr, " m_inverse %6.4f, ainj %6.4g, add %6.4g, time %6.4f\n",
         * CD->Vcele[i].area/(CD->Vcele[j].porosity * CD->Vcele[j].vol),
         * CD->Vcele[j].t_conc[CD->pumps[i].Position_Species],
         * CD->pumps[i].flow_rate * CD->pumps[i].Injection_conc * stepsize * CD->Vcele[i].area
         * / (CD->Vcele[j].porosity * CD->Vcele[j].vol), stepsize);
         *
         * }
         * }
         */

        t_start_react = time(NULL);
        timer2 = timer();
        //    fprintf(stderr, " takes %6.4f seconds\n", timer2- timer1);

        if (int_flg)
        {
            for (i = 0; i < CD->NumVol; i++)
            {
                CD->Vcele[i].height_o = CD->Vcele[i].height_t;
                CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
            }
        }

        timer1 = timer();
        m = 0;

        if ((!CD->RecFlg) &&
            ((int)(timelps + stepsize) % (int)(CD->React_delay * stepsize) ==
                0))
        {

// 10.05
#ifdef _OPENMP
#pragma omp parallel for        // 11.08 exactly the same
#endif
            for (i = 0; i < num_blocks; i++)
            {
                // 10.01 test purpose
                //fprintf(stderr, "  CD->Vcele[i, %d].illness = %f \n", (i+1), CD->Vcele[i].illness);

                // 10.05 test purpose
                //int ID = omp_get_thread_num();
                //printf(" (timelps = %f) i = %d, by # of thread =  %d \n", timelps, i, ID);

                // 10.05 need privating!!
                int             k, j, nr_tmp = 1;
                double          substep;

                if (CD->Vcele[i].illness < 20)
                    //if (React(timelps - (CD->React_delay * stepsize), stepsize * CD->React_delay, CD, i, &nr_tmp, Model_Data)) // 08.22 add model data
                    if (React(timelps - (CD->React_delay * stepsize), stepsize * CD->React_delay, CD, i, &nr_tmp, pihm))    // 10.01
                    {
                        fprintf(stderr, "  ---> React failed at cell %12d.\t",
                            CD->Vcele[i].index);

                        substep = 0.5 * stepsize;
                        k = 2;

                        //while (j = React(timelps, substep, CD, i, &nr_tmp, Model_Data)) // 08.22 add model data
                        while ((j = React(timelps, substep, CD, i, &nr_tmp, pihm)))   // 10.01
                        {
                            substep = 0.5 * substep;
                            k = 2 * k;
                            if (substep < 0.5)
                                break;
                        }

                        if (j == 0)
                        {
                            tot_nr += nr_tmp;
                            fprintf(stderr,
                                " Reaction passed with step equals to %f (1/%d)\n",
                                substep, k);
                            for (j = 1; j < k; j++)
                            {
                                //React(timelps + j * substep, substep, CD, i, &nr_tmp, Model_Data);  // 08.22 add model data
                                React(timelps + j * substep, substep, CD, i, &nr_tmp, pihm);    // 10.01
                                tot_nr += nr_tmp;
                            }
                        }

                        //              else
                        //        {
                        //          fprintf(stderr, " Taking smaller step does not work, altering initial guess!\n");
                        /*
                         * for ( j = 0; j < CD->NumStc; j ++){
                         * if (CD->chemtype[j].itype == 4)
                         * CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                         * else
                         * CD->Vcele[i].p_conc[j] = fabs(CD->Vcele[i].t_conc[j]*0.1);
                         * }
                         * //           CD->Vcele[i].illness ++ ;
                         * for (j = 0; j < k; j ++){
                         * m = MAX(React(timelps + j * substep, substep, CD, i , &nr_tmp),m);
                         * }
                         */
                        //              if ( m == 1) ReportError(CD->Vcele[i],CD);
                        //        }
                    }
                tot_nr += nr_tmp;
            }
// 10.05
// ending of OpenMP looping
        }

        timer2 = timer();
        timelps += stepsize;
        if (timelps >= org_time + hydro_step)
            break;
    }

    if ((!CD->RecFlg) &&
        ((int)(timelps) % (int)(CD->React_delay * stepsize) == 0))
    {
        // 10.04
        /*
         * fprintf(stderr, "    React from %.1f to %.1f [min]. Average newton iteration = %.3f. Elapsed time = %6.4f seconds. \n", timelps - (CD->React_delay * stepsize), \
         * timelps, (double)tot_nr / (double)CD->NumEle / 2.0, timer2 - timer1);
         */
    }

    // 10.05
    t_end_react = time(NULL);
    *t_duration_react += (t_end_react - t_start_react);
}

void FreeChem(Chem_Data CD)     // 10.01
{
    int             i;
    int             num_dep = 4;

    free(CD->BTC_loc);
    free(CD->prepconcindex);

    for (i = 0; i < CD->NumSsc; i++)
    {
        free(CD->Dependency[i]);
    }
    free(CD->Dependency);

    for (i = 0; i < CD->NumMkr + CD->NumAkr; i++)
    {
        free(CD->Dep_kinetic[i]);
    }
    free(CD->Dep_kinetic);

    for (i = 0; i < CD->NumMin; i++)
    {
        free(CD->Dep_kinetic_all[i]);
    }
    free(CD->Dep_kinetic_all);

    for (i = 0; i < CD->NumStc; i++)
    {
        free(CD->Totalconc[i]);
        free(CD->Totalconck[i]);
    }
    free(CD->Totalconc);
    free(CD->Totalconck);

    free(CD->Keq);
    free(CD->KeqKinect);
    free(CD->KeqKinect_all);

    // CD->Vcele
    free(CD->Vcele->t_conc);
    free(CD->Vcele->t_rate);
    free(CD->Vcele->t_tol);
    free(CD->Vcele->p_conc);
    free(CD->Vcele->s_conc);
    free(CD->Vcele->p_actv);
    free(CD->Vcele->s_actv);
    free(CD->Vcele->p_para);
    free(CD->Vcele->p_type);
    free(CD->Vcele);

    free(CD->Flux);

    free(CD->chemtype->ChemName);
    free(CD->chemtype);

    // CD->kinetics
    if (CD->NumMkr > 0)
    {
        free(CD->kinetics->species);
        free(CD->kinetics->Label);
        for (i = 0; i < num_dep; i++)
        {
            free(CD->kinetics->dep_species[i]);
            free(CD->kinetics->monod_species[i]);
            free(CD->kinetics->inhib_species[i]);
        }
        free(CD->kinetics->dep_species);
        free(CD->kinetics->monod_species);
        free(CD->kinetics->inhib_species);
        free(CD->kinetics->dep_position);
        free(CD->kinetics->dep_power);
        free(CD->kinetics->biomass_species);
        free(CD->kinetics->monod_position);
        free(CD->kinetics->monod_para);
        free(CD->kinetics->inhib_position);
        free(CD->kinetics->inhib_para);
        free(CD->kinetics);
    }

    free(CD->pumps->Name_Species);
    free(CD->pumps);

    // CD->TSD_prepconc
    for (i = 0; i < CD->TSD_prepconc->length; i++)
    {
        free(CD->TSD_prepconc->TS[i]);
    }
    free(CD->TSD_prepconc->TS);
    free(CD->TSD_prepconc);


}
