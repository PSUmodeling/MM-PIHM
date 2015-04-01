/******************************************************************************
 * PIHM-RT is a finite volume based, reactive transport module that operate
 * on top of the hydrological processes described by PIHM. PIHM-RT track the 
 * transportation and reaction in a given watershed. PIHM-RT uses operator
 * splitting technique to couple transport and reaction. 
 *
 * PIHM-RT requires two additional input files: 
 *     a. chemical condition file:     projectname.chem
 *     b. index of initial conditions: projectname.cini
 *
 *
 *
 * If you have any questions, concerns, suggestions, please contact me at 
 * the following address:
 *
 *     Developer: Chen Bao <baochen.d.s@gmail.com>
 *     Version  : 0.4
 *     Date     : Feb, 2015   
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>

#include "../pihm.h"		/* Data Model and Variable Declarations     */
#include "rt.h"			/* Data Model and Variable Declarations for chemical processes */

#define UNIT_C 1440
#define ZERO   1E-20
#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80
#define TIME_MIN  1E-5
#define EPS       0.05
#define INFTYSMALL  1E-6
#define RTdepth 5.0

/* Functions declarations and usage */
realtype rivArea (int, realtype, realtype);
realtype returnVal (realtype rArea, realtype rPerem, realtype eqWid,
		    realtype ap_Bool);
realtype CS_AreaOrPerem (int rivOrder, realtype rivDepth,
			 realtype rivCoeff, realtype a_pBool);
static double timer ();
// timer, system function called to time subroutines
void Monitor (realtype, realtype, void *, Chem_Data);
// adjust unphysical PIHM flux outputs by mass balance
int upstream (element, element, const void *);
// locate upstream nodes for TVD calculation
int realcheck (const char *);
// check real number or real number range
int keymatch (const char *, const char *, double *, char **);
// keyword matching and data reading
void ConditionAssign (int, char *, int *);
// Assign conditions to different cells
void chem_alloc (char *, const void *, const Control_Data, Chem_Data,
		 realtype);
// chemical component initialization
void fluxtrans (realtype, realtype, const void *, Chem_Data, N_Vector);
// translational function to take PIHM outputs to RT calculation
void chem_updater (Chem_Data, void *);
// unused subroutine to update field properties from chemical reactions
void OS3D (realtype, realtype, Chem_Data);
// operator splitting 3D (finite volume) for transport
int React (realtype, realtype, Chem_Data, int, int *);
// kinetic reaction component
void Lookup (FILE *, Chem_Data, int);
// database fetching
int Speciation (Chem_Data, int);
// chemical condition initialization from 1) total concentration. 2) total concentration with pH
int SpeciationType (FILE *, char *);
// determine the type of solutes from database
void AdptTime (Chem_Data, realtype, double, double, double *, int);
// adaptive time stepping
void Reset (Chem_Data, int);
// unused subroutine to reset chemical conditions for failed cell. Not used anymore given successful time control

/* Fucntion declarations finished   */
// Timer
     static double timer ()
{
  struct timeval tp;
  gettimeofday (&tp, NULL);
  return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}

void
Monitor (realtype t, realtype stepsize, void *DS, Chem_Data CD)
{
  /* unit of t and stepsize: min */
  /* DS: model data              */

  // this is to obtain the infiltration rate that can not be obtained from reading f.c 
  // f.c outputs the last trial value, rather the best values of infiltration

  int i, j, k = 0, num_face = 0;
  struct tm *timestamp;
  time_t *rawtime;
  double timelps, sumflux1, sumflux2, correction, swi,
    inv_swi, partratio, sumlateral, depth, hu, hg, hn, ht,
    imrecharge, flux, A;
  realtype MF_CONVERT;

  MF_CONVERT = (realtype) (24 * 60 * 60);

  Model_Data MD;
  MD = (Model_Data) DS;
  rawtime = (time_t *) malloc (sizeof (time_t));
  *rawtime = (int) (t * 60);
  timestamp = gmtime (rawtime);
  timelps = t - CD->StartTime;

  double unit_c = stepsize / UNIT_C;

  //    fprintf(logfile,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\n",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
  //    fprintf(logfile," Time step is %6.4f\n", stepsize);


  double *tmpflux = (double *) malloc (CD->NumOsv * sizeof (double));
  double *resflux = (double *) malloc (CD->NumOsv * sizeof (double));
  swi = 0.2;
  inv_swi = 2.0 / (1.0 - swi);
  for (i = 0; i < CD->NumOsv; i++)
    {
      tmpflux[i] = 0.0;
      resflux[i] = 0.0;
    }

  for (j = 0; j < CD->NumEle; j++)
    {
      hu = CD->Vcele[j + MD->NumEle].height_t;
      hg = CD->Vcele[j].height_t;
      depth = CD->Vcele[j].height_v;
      hn = ((swi + 1.0) * 0.5 * (depth - hg) - hu) * inv_swi;
      ht = depth - hg - hn;

      if ((hg * ht > 0))
	partratio = ht * 0.30 / hg;	// One forth indicate third order relationship between S and Kr, and so on
      if (hg <= 0)
	partratio = 1.00E3;	// no groundwater and flow essentially
      if (ht <= 0)
	partratio = 1.00E-3;	// no transient zone flow essentially

      A = partratio / (1 + partratio);
      tmpflux[j] = A;
    }


  for (i = 0; i < CD->NumFac; i++)
    {
      if (!CD->Flux[i].flux_type)
	{
	  resflux[CD->Flux[i].nodeup - 1] -= CD->Flux[i].flux * unit_c;	// sum lateral fluxes
	  //      if ( CD->Flux[i].nodeup == 1)
	  //      fprintf(logfile, " 1 flux: %f\t", CD->Flux[i].flux);
	}
    }

  for (i = 0; i < CD->PIHMFac; i++)
    {

      if (CD->Flux[i].nodelo <= CD->NumEle)	// averaging the partition ratio between neighbors. 
	A = (tmpflux[CD->Flux[i].nodeup - 1]
	     + tmpflux[CD->Flux[i].nodelo - 1]) * 0.5;
      else
	{
	  A = tmpflux[CD->Flux[i].nodeup - 1];
	  //      fprintf(stderr, " %d to %d with %f\n", CD->Flux[i].nodeup, CD->Flux[i].nodelo, A);
	}
      //      if ( A > 1) fprintf(stderr, " partition is wrong\n");

      /*      
         flux   = CD->Flux[i].flux;

         CD->Flux[ i + CD->PIHMFac].flux     = A * CD->Flux[i].flux;
         CD->Flux[ i + CD->PIHMFac].velocity = A * CD->Flux[i].velocity ; // this is not accurate, but check if we really used velocity in CFL and OS3D

         CD->Flux[ i ].flux                  = ( 1.0 - A ) * CD->Flux[i].flux;
         CD->Flux[ i ].velocity              = ( 1.0 - A ) * CD->Flux[i].velocity;
       */
      //      fprintf(stderr, " %f=%f\t", 1.0-A, CD->Flux[i].flux/flux);
    }

  for (i = 0; i < CD->NumFac; i++)
    {
      if (!CD->Flux[i].flux_type)
	{
	  tmpflux[CD->Flux[i].nodeup - 1] = resflux[CD->Flux[i].nodeup - 1];
	}
    }

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
	  //      if ( CD->Flux[i].nodeup== 1)
	  //        fprintf(logfile, " 1 flux: %f\t", CD->Flux[i].flux);
	}
    }

  //    fprintf(logfile, "h_n+1, h_n, dV, SumLateral_o, SumLateral_n, dV-SL, Recharge_o, Recharge_n, RelDiff, Correction \n hu, hg, depth, hn, ht, partratio, outSL, ImRecharge, TotRecharge\n");
  for (i = 0; i < CD->NumEle; i++)
    {
      sumflux1 = (CD->Vcele[i].height_t - CD->Vcele[i].height_o)
	* MD->Ele[i].area * MD->Ele[i].Porosity;
      sumflux2 = sumflux1 - resflux[i];
      correction = -sumflux2 * UNIT_C / stepsize
	/ CD->Flux[CD->Vcele[i].ErrDumper].flux;
      A = -CD->Flux[CD->Vcele[i].ErrDumper].flux;
      CD->Flux[CD->Vcele[i].ErrDumper].flux = -sumflux2 * UNIT_C / stepsize;
      CD->Flux[CD->Vcele[i].ErrDumper].velocity =
	CD->Flux[CD->Vcele[i].ErrDumper].flux
	/ CD->Flux[CD->Vcele[i].ErrDumper].s_area;
      CD->Flux[CD->Vcele[i].ErrDumper - MD->NumEle].flux =
	-CD->Flux[CD->Vcele[i].ErrDumper].flux;
      CD->Flux[CD->Vcele[i].ErrDumper - MD->NumEle].velocity =
	-CD->Flux[CD->Vcele[i].ErrDumper].velocity;
      // fprintf(logfile, "%d\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n", i+1, CD->Vcele[i].height_t, CD->Vcele[i].height_o, sumflux1, tmpflux[i], resflux[i], sumflux2, A, sumflux2 * UNIT_C/ stepsize , (sumflux1-resflux[i] + CD->Flux[CD->Vcele[i].ErrDumper].flux * unit_c )/sumflux1 * 100, correction);
    }

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
	  //      if ( CD->Flux[i].nodeup== 1)
	  //        fprintf(logfile, " 1 flux: %f\t", CD->Flux[i].flux);
	}
    }


  for (i = CD->NumEle; i < CD->NumEle * 2; i++)
    {
      sumflux1 = (CD->Vcele[i].height_t - CD->Vcele[i].height_o)
	* MD->Ele[i - MD->NumEle].area * MD->Ele[i - MD->NumEle].Porosity;
      sumflux2 = sumflux1 - resflux[i];
      correction = -sumflux2 * UNIT_C / stepsize / CD->Vcele[i].q;
      A = CD->Flux[CD->Vcele[i].ErrDumper].flux;
      CD->Vcele[i].q = sumflux2 * UNIT_C / stepsize;
      CD->Vcele[i].q = MAX (CD->Vcele[i].q, 0.0);
      // input of rain water chemistry can not be negative;
      CD->Vcele[i].q += fabs (MD->EleET[i - MD->NumEle][2] * MF_CONVERT)
	* MD->Ele[i - MD->NumEle].area;
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

  //    CD->Condensation = CD->Condensation / CD->NumEle;
  //    fprintf(stderr, " ## Condensation factor %f ##\n", CD->Condensation);

  //    fclose(logfile);

  free (tmpflux);
  free (resflux);
  free (rawtime);
}

int
upstream (element up, element lo, const void *DS)
{
  /* Locate the upstream grid of up -> lo flow */
  /* Require verification                      */
  /* only determines points in triangular elements */

  double x_, y_, x_a, x_b, x_c, y_a, y_b, y_c,
    dot00, dot01, dot02, dot11, dot12, u, v, invDenom;

  int i, j, upstreamfound = 0;
  x_ = 2 * up.x - lo.x;
  y_ = 2 * up.y - lo.y;

  Model_Data MD;
  MD = (Model_Data) DS;

  //  fprintf(stderr,"%d\t%d\t%d\t%d\t",MD->NumEle,MD->Ele[up.nabr[0]].index,MD->Ele[up.nabr[1]].index,MD->Ele[up.nabr[2]].index);

  for (i = 0; i < MD->NumEle; i++)
    {
      /* find point lies in which triangular element, a very interesting method */
      if ((i != (up.index - 1)) && (i != (lo.index - 1)))
	{
	  x_a = MD->Node[MD->Ele[i].node[0] - 1].x;
	  x_b = MD->Node[MD->Ele[i].node[1] - 1].x;
	  x_c = MD->Node[MD->Ele[i].node[2] - 1].x;
	  y_a = MD->Node[MD->Ele[i].node[0] - 1].y;
	  y_b = MD->Node[MD->Ele[i].node[1] - 1].y;
	  y_c = MD->Node[MD->Ele[i].node[2] - 1].y;
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
	      return (MD->Ele[i].index);
	    }
	}
    }
  if (upstreamfound == 0)
    {
      //    fprintf(stderr, "The upstream of %d and %d is beyond boundaries!\n",up.index, lo.index);
      return (0);
    }
  return (0);
}

int
realcheck (const char *words)
{

  int flg = 1, i;
  if (((words[0] < 58) && (words[0] > 47)) || (words[0] == 46)
      || (words[0] == 45) || (words[0] == 43))
    {
      for (i = 0; i < strlen (words); i++)
	if ((words[i] > 57 || words[i] < 43) && (words[i] != 69)
	    && (words[i] != 101) && (words[i] != 10) && (words[i] != 13))
	  flg = 0;
    }
  else
    flg = 0;
  return (flg);
}


int
keymatch (const char *line, const char *keyword, double *value, char **strval)
{
  /* A very general and convinient way of reading datafile and input file */
  /* find keyword in line, assign the value after keyword to value array if there is any */
  /* store both numbers and strings in order for later use, buffer required */
  /* if is keyword not found return 0. If comments, return 2. Otherwise return 1 */
  int i;

  for (i = 0; i < WORDS_LINE; i++)
    value[i] = 0.0;

  if ((line[0] == '!') || (line[0] == '#'))
    {

      return (2);
      /* assign a special flag for comments */
    }

  int j, k, line_width, word_width = WORD_WIDTH, quoteflg = 0;
  int words_line = WORDS_LINE;
  int keyfoundflag = 0;

  char **words;
  words = (char **) malloc (WORDS_LINE * sizeof (char *));

  for (i = 0; i < WORDS_LINE; i++)
    {
      words[i] = (char *) malloc (WORD_WIDTH * sizeof (char));
      memset (words[i], 0, WORD_WIDTH);
    }
  i = j = k = 0;

  /* Partition the line into words */
  while (i < strlen (line))
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
    if (strcmp (words[i], keyword) == 0)
      keyfoundflag = 1;

  j = k = 0;
  for (i = 0; i < words_line; i++)
    {
      //    fprintf(stderr, "word#%d=%s, length=%d\n" , i, words[i], strlen(words[i]));
      strcpy (strval[k++], words[i]);
      //    if ((( words[i][0] < 58)&&(words[i][0] > 47))||(words[i][0] == 46)||(words[i][0]==45)||(words[i][0]==43))
      if (realcheck (words[i]) == 1)
	value[j++] = atof (words[i]);
    }

  for (i = 0; i < WORDS_LINE; i++)
    free (words[i]);
  free (words);
  return (keyfoundflag);

}

void
ConditionAssign (int condition, char *str, int *index)
{
  /* This subroutine takes in input strings and output an index array that record the conditions each blocks assigned to */
  /* input strings could use separators like - and , */

  int i, j, k, l, length = strlen (str);
  char **words = (char **) malloc (length * sizeof (char *));
  for (i = 0; i < length; i++)
    words[i] = (char *) malloc (WORD_WIDTH * sizeof (char));

  char *tmpstr = (char *) malloc (WORD_WIDTH * sizeof (char));

  char *separator = (char *) malloc (length * sizeof (char));
  int *value = (int *) malloc (length * sizeof (char));

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
      strcpy (tmpstr, words[i]);
      value[i] = atoi (tmpstr);
    }
  /* 
     for ( i = 0; i < length; i ++)
     fprintf(stderr, "%s\t%c\n",words[i],separator[i]);

     for ( i = 0; i < length; i ++)
     fprintf(stderr, "%d\n",value[i]);

     fprintf(stderr, "condition = %d\n", condition);
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
    free (words[i]);
  free (words);
  free (separator);
  free (value);
  free (tmpstr);

}




void
chem_alloc (char *filename, const void *DS, const Control_Data CS,
	    Chem_Data CD, realtype t)
{

  int i, j, k, num_face =
    0, num_species, num_mineral, num_ads, num_cex, num_other, num_conditions =
    0, line_width = LINE_WIDTH, words_line = WORDS_LINE, word_width =
    WORD_WIDTH, Global_diff = 0, Global_disp = 0, error_flag =
    0, speciation_flg = 0, specflg;
  double total_flux = 0.0, total_area = 0.0, tmpval[WORDS_LINE];
  time_t rawtime;
  struct tm *timeinfo;

  Model_Data MD;
  MD = (Model_Data) DS;

  char keyword[WORD_WIDTH], line[256], word[WORD_WIDTH];
  char **tmpstr = (char **) malloc (WORDS_LINE * sizeof (char *));

  timeinfo = (struct tm *) malloc (sizeof (struct tm));

  for (i = 0; i < words_line; i++)
    tmpstr[i] = (char *) malloc (WORD_WIDTH * sizeof (char));

  char *chemfn = (char *) malloc ((strlen (filename) + 30) * sizeof (char));
  sprintf (chemfn, "input/%s/%s.chem", filename, filename);
  FILE *chemfile = fopen (chemfn, "r");
  char *datafn = (char *) malloc ((strlen (filename) + 30) * sizeof (char));
  sprintf (datafn, "input/%s/%s.cdbs", filename, filename);
  FILE *database = fopen (datafn, "r");
  char *forcfn = (char *) malloc ((strlen (filename) + 30) * sizeof (char));
  sprintf (forcfn, "input/%s/%s.prep", filename, filename);
  FILE *prepconc = fopen (forcfn, "r");


  assert (chemfile != NULL);
  assert (database != NULL);
  assert (prepconc != NULL);

  /* get rid of the following paragraph after testing real input */


  CD->NumVol = 2 * (MD->NumEle + MD->NumRiv) + 2;
  CD->NumOsv = CD->NumVol - 2;
  CD->NumEle = MD->NumEle;
  CD->NumRiv = MD->NumRiv;
  // sat, unsat gw elements, river elements, ERB elements and their outlets

  CD->StartTime = t;
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
  /* default control variable if not found in input file */


  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "RUNTIME", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      fgets (line, line_width, chemfile);
      if (keymatch (line, "tvd", tmpval, tmpstr) == 1)
	{
	  if (strcmp (tmpstr[1], "false") == 0)
	    CD->TVDFlg = 0;
	  if (strcmp (tmpstr[1], "true") == 0)
	    CD->TVDFlg = 1;
	  if (strcmp (tmpstr[1], "false") && strcmp (tmpstr[1], "true"))
	    fprintf (stderr, "TVD FLAG INPUT ERROR!\n");
	  fprintf (stderr, " Total variation diminishing set to %d %s.\n",
		   CD->TVDFlg, tmpstr[1]);
	}
      if (keymatch (line, "output", tmpval, tmpstr) == 1)
	{
	  CD->OutItv = (int) tmpval[0];
	  fprintf (stderr, " Output interval set to %d hours.\n", CD->OutItv);
	}
      if (keymatch (line, "activity", tmpval, tmpstr) == 1)
	{
	  CD->ACTmod = (int) tmpval[0];
	  fprintf (stderr, " Activity correction is set to %d.\n",
		   CD->ACTmod);
	  // 0 for unity activity coefficient and 1 for DH equation update
	}
      if (keymatch (line, "act_coe_delay", tmpval, tmpstr) == 1)
	{
	  CD->DHEdel = (int) tmpval[0];
	  fprintf (stderr,
		   " Activity coefficient update delay is set to %d.\n",
		   CD->DHEdel);
	  // 0 for delay and 1 for no delay (solving together )
	}
      if (keymatch (line, "thermo", tmpval, tmpstr) == 1)
	{
	  CD->TEMcpl = (int) tmpval[0];
	  fprintf (stderr, " Coupling of thermo modelling is set to %d.\n",
		   CD->DHEdel);
	  // 0 for delay and 1 for no delay (solving together )                            
	}
      if (keymatch (line, "relmin", tmpval, tmpstr) == 1)
	{
	  CD->RelMin = (int) tmpval[0];
	  switch (CD->RelMin)
	    {
	    case 0:
	      fprintf (stderr, " Using absolute mineral volume fraction.\n");
	      break;
	    case 1:
	      fprintf (stderr, " Using relative mineral volume fraction.\n");
	      break;
	    }
	}
      if (keymatch (line, "effads", tmpval, tmpstr) == 1)
	{
	  CD->EffAds = (int) tmpval[0];
	  switch (CD->EffAds)
	    {
	    case 0:
	      fprintf (stderr, " Using the normal adsorption model.\n");
	      break;
	    case 1:
	      fprintf (stderr,
		       " Using the coupled MIM and adsorption model. \n");
	      break;
	      // under construction.
	    }
	}
      if (keymatch (line, "transport_only", tmpval, tmpstr) == 1)
	{
	  CD->RecFlg = (int) tmpval[0];
	  switch (CD->RecFlg)
	    {
	    case 0:
	      fprintf (stderr, " Transport only mode disabled.\n");
	      break;
	    case 1:
	      fprintf (stderr, " Transport only mode enabled. \n");
	      break;
	      // under construction.
	    }
	}
      if (keymatch (line, "precipitation", tmpval, tmpstr) == 1)
	{
	  CD->PrpFlg = (int) tmpval[0];
	  switch (CD->PrpFlg)
	    {
	    case 0:
	      fprintf (stderr, " No precipitation condition.\n");
	      break;
	    case 1:
	      fprintf (stderr,
		       " Precipitation condition is to be specified. \n");
	      break;
	    case 2:
	      fprintf (stderr,
		       " Precipitation condition is specified via file .prep. \n");
	      break;
	      // under construction.
	    }
	}
      if (keymatch (line, "RT_delay", tmpval, tmpstr) == 1)
	{
	  CD->Delay = (int) tmpval[0];
	  fprintf (stderr,
		   " Flux-PIHM-RT will start after running PIHM for %d days.\n",
		   CD->Delay);
	  CD->Delay *= UNIT_C;
	  // under construction.
	}
      if (keymatch (line, "Condensation", tmpval, tmpstr) == 1)
	{
	  CD->Condensation = tmpval[0];
	  fprintf (stderr,
		   " The concentrations of infiltrating rainfall is set to be %f times of concentrations in precipitation.\n",
		   CD->Condensation);
	  // under construction.
	  CD->Condensation *= CS->Cal.Prep_conc;
	  fprintf (stderr,
		   " The concentrations of infiltrating rainfall is set to be %f times of concentrations in precipitation.\n",
		   CD->Condensation);
	}
      if (keymatch (line, "SUFEFF", tmpval, tmpstr) == 1)
	{
	  CD->SUFEFF = tmpval[0];
	  fprintf (stderr, " Effective surface area mode set to %d.\n",
		   CD->SUFEFF);
	  // under construction.
	}
      if (keymatch (line, "AvgScl", tmpval, tmpstr) == 1)
	{
	  CD->React_delay = tmpval[0];
	  fprintf (stderr, " Averaging window for asynchronous reaction%d.\n",
		   CD->React_delay);
	  // under construction.
	}
      if (keymatch (line, "Mobile_exchange", tmpval, tmpstr) == 1)
	{
	  CD->TimRiv = tmpval[0];
	  fprintf (stderr, " Ratio of immobile ion exchange site %f.\n",
		   CD->TimRiv);
	  // under construction.
	}
    }


  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "OUTPUT", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  CD->NumBTC = tmpval[0];
  fprintf (stderr, " %d breakthrough points specified\n", CD->NumBTC);
  CD->BTC_loc = (int *) malloc (CD->NumBTC * sizeof (int));
  i = 0;
  fprintf (stderr, " --");
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      fgets (line, line_width, chemfile);
      if (keymatch (line, " ", tmpval, tmpstr) != 2)
	{
	  CD->BTC_loc[i] = (int) tmpval[0];
	  fprintf (stderr, " Grid %d ", CD->BTC_loc[i]);
	  i++;
	}
      if (i >= CD->NumBTC)
	break;
    }
  fprintf (stderr, " are breakthrough points\n");

  species Global_type;
  Global_type.ChemName = (char *) malloc (WORD_WIDTH * sizeof (char));
  strcpy (Global_type.ChemName, "GLOBAL");

  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "GLOBAL", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      fgets (line, line_width, chemfile);
      if (keymatch (line, "t_species", tmpval, tmpstr) == 1)
	{
	  CD->NumStc = (int) tmpval[0];
	  fprintf (stderr, " %d chemical species specified.\n", CD->NumStc);
	  /* H2O is always a primary species */
	}
      if (keymatch (line, "s_species", tmpval, tmpstr) == 1)
	{
	  CD->NumSsc = (int) tmpval[0];
	  fprintf (stderr, " %d secondary species specified.\n",
		   (int) tmpval[0]);
	}
      if (keymatch (line, "minerals", tmpval, tmpstr) == 1)
	{
	  CD->NumMin = (int) tmpval[0];
	  fprintf (stderr, " %d minerals specified.\n", CD->NumMin);
	}
      if (keymatch (line, "adsorption", tmpval, tmpstr) == 1)
	{
	  CD->NumAds = (int) tmpval[0];
	  fprintf (stderr, " %d surface complexation specified.\n",
		   CD->NumAds);
	}
      if (keymatch (line, "cation_exchange", tmpval, tmpstr) == 1)
	{
	  CD->NumCex = (int) tmpval[0];
	  fprintf (stderr, " %d cation exchange specified.\n", CD->NumCex);
	}
      if (keymatch (line, "mineral_kinetic", tmpval, tmpstr) == 1)
	{
	  CD->NumMkr = (int) tmpval[0];
	  fprintf (stderr, " %d mineral kinetic reaction(s) specified.\n",
		   CD->NumMkr);
	}
      if (keymatch (line, "aqueous_kinetic", tmpval, tmpstr) == 1)
	{
	  CD->NumAkr = (int) tmpval[0];
	  fprintf (stderr, " %d aqueous kinetic reaction(s) specified.\n",
		   CD->NumAkr);
	}
      if (keymatch (line, "diffusion", tmpval, tmpstr) == 1)
	{
	  fprintf (stderr, " Diffusion coefficient =%6.4f cm2/s.\n",
		   tmpval[0]);
	  Global_type.DiffCoe = tmpval[0] * 60.0 * 60.0 * 24.0 / 10000.0;
	  Global_diff = 1;
	  /* Require unit conversion ! */
	}
      if (keymatch (line, "dispersion", tmpval, tmpstr) == 1)
	{
	  fprintf (stderr, " Dispersion coefficient =%6.4f m.\n", tmpval[0]);
	  Global_type.DispCoe = tmpval[0];
	  Global_disp = 1;
	  /* Set global flags to indicate the global values are present */
	}
      if (keymatch (line, "cementation", tmpval, tmpstr) == 1)
	{
	  fprintf (stderr, " Cementation factor set to %6.4f. \n", tmpval[0]);
	  CD->Cementation = tmpval[0];
	}
      if (keymatch (line, "temperature", tmpval, tmpstr) == 1)
	{
	  CD->Temperature = tmpval[0];
	  fprintf (stderr, " Temperature set to %6.4f. \n", CD->Temperature);
	}
    }

  CD->NumSpc = CD->NumStc - (CD->NumMin + CD->NumAds + CD->NumCex);
  /* the number of species that are mobile, later used in the OS3D subroutine */
  CD->NumSdc = CD->NumStc - (CD->NumMin);
  /* the number of species that others depend on */

  CD->Dependency = (double **) malloc (CD->NumSsc * sizeof (double *));
  for (i = 0; i < CD->NumSsc; i++)
    {
      CD->Dependency[i] = (double *) malloc (CD->NumSdc * sizeof (double));
      /* convert secondary species as an expression of primary species */
      for (j = 0; j < CD->NumSdc; j++)
	CD->Dependency[i][j] = 0.0;
    }
  CD->Dep_kinetic =
    (double **) malloc ((CD->NumMkr + CD->NumAkr) * sizeof (double *));
  for (i = 0; i < CD->NumMkr + CD->NumAkr; i++)
    {
      CD->Dep_kinetic[i] = (double *) malloc (CD->NumStc * sizeof (double));
      /* express kinetic species as function of primary species */
      for (j = 0; j < CD->NumStc; j++)
	CD->Dep_kinetic[i][j] = 0.0;
    }
  CD->Dep_kinetic_all = (double **) malloc ((CD->NumMin) * sizeof (double *));
  for (i = 0; i < CD->NumMin; i++)
    {
      CD->Dep_kinetic_all[i] =
	(double *) malloc (CD->NumStc * sizeof (double));
      /* Dependencies of minearls, all */
      for (j = 0; j < CD->NumStc; j++)
	CD->Dep_kinetic_all[i][j] = 0.0;
    }
  CD->Keq = (double *) malloc (CD->NumSsc * sizeof (double));
  CD->KeqKinect =
    (double *) malloc ((CD->NumMkr + CD->NumAkr) * sizeof (double));
  CD->KeqKinect_all = (double *) malloc (CD->NumMin * sizeof (double));
  /* Keqs of equilibrium/ kinetic and kinetic all */


  CD->Totalconc = (double **) malloc (CD->NumStc * sizeof (double *));
  for (i = 0; i < CD->NumStc; i++)
    CD->Totalconc[i] =
      (double *) malloc ((CD->NumStc + CD->NumSsc) * sizeof (double));
  /* convert total concentration as an expression of all species */

  CD->Totalconck = (double **) malloc (CD->NumStc * sizeof (double *));
  for (i = 0; i < CD->NumStc; i++)
    CD->Totalconck[i] =
      (double *) malloc ((CD->NumStc + CD->NumSsc) * sizeof (double));
  /* convert total concentration as an expression of all species */

  for (i = 0; i < CD->NumStc; i++)
    for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
      {
	CD->Totalconc[i][j] = 0.0;
	CD->Totalconck[i][j] = 0.0;
      }


  num_species = CD->NumSpc;

  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "INITIAL_CONDITIONS", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      if (keymatch (line, " ", tmpval, tmpstr) != 2)
	{
	  num_conditions++;
	}
      fgets (line, line_width, chemfile);
    }
  fprintf (stderr, " %d conditions assigned.\n", num_conditions);



  char **chemcon = (char **) malloc (num_conditions * sizeof (char *));
  for (i = 0; i < num_conditions; i++)
    chemcon[i] = (char *) malloc (word_width * sizeof (char));
  char ***con_chem_name =
    (char ***) malloc ((num_conditions + 1) * sizeof (char **));
  for (i = 0; i < num_conditions + 1; i++)
    {				// all conditions + precipitation
      con_chem_name[i] = (char **) malloc (CD->NumStc * sizeof (char *));
      for (j = 0; j < CD->NumStc; j++)
	con_chem_name[i][j] = (char *) malloc (WORD_WIDTH * sizeof (char));
    }





  int *condition_index = (int *) malloc ((CD->NumVol + 1) * sizeof (int));
  /* when user assign conditions to blocks, they start from 1 */

  for (i = 0; i < CD->NumVol; i++)
    condition_index[i] = 0;

  vol_conc *Condition_vcele =
    (vol_conc *) malloc (num_conditions * sizeof (vol_conc));
  for (i = 0; i < num_conditions; i++)
    {
      Condition_vcele[i].index = i + 1;
      Condition_vcele[i].t_conc =
	(double *) malloc (CD->NumStc * sizeof (double));
      Condition_vcele[i].p_conc =
	(double *) malloc (CD->NumStc * sizeof (double));
      Condition_vcele[i].p_para =
	(double *) malloc (CD->NumStc * sizeof (double));
      Condition_vcele[i].p_type = (int *) malloc (CD->NumStc * sizeof (int));
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
	(double *) malloc (CD->NumStc * sizeof (double));
      CD->Precipitation.p_conc =
	(double *) malloc (CD->NumStc * sizeof (double));
      CD->Precipitation.p_para =
	(double *) malloc (CD->NumStc * sizeof (double));
      CD->Precipitation.p_type = (int *) malloc (CD->NumStc * sizeof (int));
      CD->Precipitation.s_conc = NULL;
      for (i = 0; i < CD->NumStc; i++)
	{
	  CD->Precipitation.t_conc[i] = ZERO;
	  CD->Precipitation.p_conc[i] = ZERO;
	}

    }

  CD->chemtype =
    (species *) malloc ((CD->NumStc + CD->NumSsc) * sizeof (species));
  if (CD->chemtype == NULL)
    fprintf (stderr, " Memory allocation error\n");

  for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {

      if (Global_diff == 1)
	CD->chemtype[i].DiffCoe = Global_type.DiffCoe;
      /*
         else
         CD->chemtype[i].DiffCoe = ZERO;
       */
      /* in squre m per day */
      if (Global_disp == 1)
	CD->chemtype[i].DispCoe = Global_type.DispCoe;
      else
	CD->chemtype[i].DispCoe = ZERO;
      CD->chemtype[i].ChemName = (char *) malloc (WORD_WIDTH * sizeof (char));
      assert (CD->chemtype[i].ChemName != NULL);
      memset (CD->chemtype[i].ChemName, 0, WORD_WIDTH);
      CD->chemtype[i].Charge = 0.0;
      CD->chemtype[i].SizeF = 1.0;
      CD->chemtype[i].itype = 0;
    }



  k = 0;
  int initfile = 0;
  FILE *cheminitfile = NULL;
  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "INITIAL_CONDITIONS", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  if (strcmp (tmpstr[1], "FILE") == 0)
    {				// initialize chemical distribution from file evoked. This will nullify all the condition assignment given in the next lines. But for now, please keep those lines to let the code work.
      fprintf (stderr,
	       " Initializing the initial chemical distribution from file %s\n",
	       tmpstr[2]);
      char cheminit[30];
      strcpy (cheminit, "input/");
      strcat (cheminit, filename);
      strcat (cheminit, "/");
      strcat (cheminit, tmpstr[2]);
      cheminitfile = fopen (cheminit, "r");
      if (cheminitfile == NULL)
	fprintf (stderr,
		 " Chemical initial distribution file is missing from input directory!\n");
      initfile = 1;
    }
  fgets (line, line_width, chemfile);
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      if (keymatch (line, " ", tmpval, tmpstr) != 2)
	{
	  strcpy (chemcon[k++], tmpstr[0]);
	  if (initfile == 0)
	    {
	      fprintf (stderr, " Condition %s assigned to cells %s.\n",
		       chemcon[k - 1], tmpstr[1]);
	      ConditionAssign (k, tmpstr[1], condition_index);
	    }
	}
      fgets (line, line_width, chemfile);
    }
  if (initfile == 1)
    {
      for (i = 0; i < CD->NumVol; i++)
	{
	  fscanf (cheminitfile, "%d %d", &k, condition_index + i + 1);
	  // fprintf(stderr, "%6d %6d %6s\n", i+1, condition_index[i+1], chemcon[condition_index[i+1]-1]);
	}
    }

  if (cheminitfile != NULL)
    fclose (cheminitfile);

  for (i = 0; i < num_conditions; i++)
    {
      rewind (chemfile);
      num_species = 0;
      num_mineral = 0;
      num_ads = 0;
      num_cex = 0;
      num_other = 0;
      fgets (line, line_width, chemfile);
      while ((keymatch (line, "Condition", tmpval, tmpstr) != 1)
	     || (keymatch (line, chemcon[i], tmpval, tmpstr) != 1))
	fgets (line, line_width, chemfile);
      if (strcmp (tmpstr[1], chemcon[i]) == 0)
	fprintf (stderr, " %s", line);
      fgets (line, line_width, chemfile);
      while (keymatch (line, "END", tmpval, tmpstr) != 1)
	{
	  if (keymatch (line, "NULL", tmpval, tmpstr) != 2)
	    {

	      specflg = SpeciationType (database, tmpstr[0]);

	      if (specflg == 1)
		{
		  num_other = num_mineral + num_ads + num_cex;
		  Condition_vcele[i].t_conc[num_species - num_other] =
		    tmpval[0];
		  strcpy (con_chem_name[i][num_species - num_other],
			  tmpstr[0]);
		  fprintf (stderr, " %s\t%6.4f\n",
			   con_chem_name[i][num_species - num_other],
			   tmpval[0]);
		  Condition_vcele[i].p_type[num_species - num_other] = 1;
		}
	      // arrange the concentration of the primary species in such a way that all the mobile species are at the beginning.
	      if (specflg == 4)
		{
		  Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds +
					    CD->NumCex + num_mineral] =
		    tmpval[0];
		  if (strcmp (tmpstr[2], "-ssa") == 0)
		    Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds +
					      CD->NumCex + num_mineral] =
		      tmpval[1] * CS->Cal.SSA;
		  strcpy (con_chem_name[i]
			  [CD->NumSpc + CD->NumAds + CD->NumCex +
			   num_mineral], tmpstr[0]);
		  fprintf (stderr,
			   " mineral %s\t%6.4f specific surface area %6.4f\n",
			   con_chem_name[i][CD->NumSpc + CD->NumAds +
					    CD->NumCex + num_mineral],
			   tmpval[0], tmpval[1]);
		  Condition_vcele[i].p_type[CD->NumSpc + CD->NumAds +
					    CD->NumCex + num_mineral] = 4;
		  num_mineral++;
		}
	      if ((tmpstr[0][0] == '>') || (specflg == 2))
		{		// adsorptive sites and species start with >
		  Condition_vcele[i].t_conc[CD->NumSpc + num_ads] = tmpval[0] * CS->Cal.Site_den;	// this is the site density of the adsorptive species.
		  Condition_vcele[i].p_type[CD->NumSpc + num_ads] = 2;
		  Condition_vcele[i].p_para[CD->NumSpc + num_ads] = 0;
		  // update when fill in the parameters for adsorption.
		  strcpy (con_chem_name[i][CD->NumSpc + num_ads], tmpstr[0]);
		  fprintf (stderr, " surface complex %s\t %6.4f\n",
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
		  strcpy (con_chem_name[i][CD->NumSpc + CD->NumAds + num_cex],
			  tmpstr[0]);
		  fprintf (stderr, " cation exchange %s\t %6.4f\n",
			   con_chem_name[i][CD->NumSpc + CD->NumAds +
					    num_cex], tmpval[0]);
		  num_cex++;
		  // under construction
		}
	      num_species++;
	    }
	  fgets (line, line_width, chemfile);
	}
    }
  if (CD->PrpFlg)
    {
      rewind (chemfile);
      fgets (line, line_width, chemfile);
      while (keymatch (line, "PRECIPITATION", tmpval, tmpstr) != 1)
	fgets (line, line_width, chemfile);
      fgets (line, line_width, chemfile);
      fprintf (stderr, " ---------------------------------\n");
      fprintf (stderr, " The condition of precipitation is \n");
      fprintf (stderr, " ---------------------------------\n");
      num_species = 0;
      num_mineral = 0;
      num_ads = 0;
      num_cex = 0;
      num_other = 0;
      while (keymatch (line, "END", tmpval, tmpstr) != 1)
	{
	  if (keymatch (line, "NULL", tmpval, tmpstr) != 2)
	    {

	      specflg = SpeciationType (database, tmpstr[0]);

	      if (specflg == 1)
		{
		  num_other = num_mineral + num_ads + num_cex;
		  CD->Precipitation.t_conc[num_species - num_other] =
		    tmpval[0];
		  strcpy (con_chem_name[num_conditions]
			  [num_species - num_other], tmpstr[0]);
		  fprintf (stderr, " %s\t%6.4f\n",
			   con_chem_name[num_conditions][num_species -
							 num_other],
			   tmpval[0]);
		  CD->Precipitation.p_type[num_species - num_other] = 1;
		}
	      // arrange the concentration of the primary species in such a way that all the mobile species are at the beginning.
	      if (specflg == 4)
		{
		  CD->Precipitation.t_conc[CD->NumSpc + CD->NumAds +
					   CD->NumCex + num_mineral] =
		    tmpval[0];
		  if (strcmp (tmpstr[2], "-ssa") == 0)
		    CD->Precipitation.p_para[CD->NumSpc + CD->NumAds +
					     CD->NumCex + num_mineral] =
		      tmpval[1];
		  strcpy (con_chem_name[num_conditions]
			  [CD->NumSpc + CD->NumAds + CD->NumCex +
			   num_mineral], tmpstr[0]);
		  fprintf (stderr,
			   " mineral %s\t%6.4f specific surface area %6.4f\n",
			   con_chem_name[num_conditions][CD->NumSpc +
							 CD->NumAds +
							 CD->NumCex +
							 num_mineral],
			   tmpval[0], tmpval[1]);
		  CD->Precipitation.p_type[CD->NumSpc + CD->NumAds +
					   CD->NumCex + num_mineral] = 4;
		  num_mineral++;
		}
	      if ((tmpstr[0][0] == '>') || (specflg == 2))
		{		// adsorptive sites and species start with >
		  CD->Precipitation.t_conc[CD->NumSpc + num_ads] = tmpval[0];	// this is the site density of the adsorptive species.
		  CD->Precipitation.p_type[CD->NumSpc + num_ads] = 2;
		  CD->Precipitation.p_para[CD->NumSpc + num_ads] = 0;
		  // update when fill in the parameters for adsorption.
		  strcpy (con_chem_name[num_conditions][CD->NumSpc + num_ads],
			  tmpstr[0]);
		  fprintf (stderr, " surface complex %s\t %6.4f\n",
			   con_chem_name[num_conditions][CD->NumSpc +
							 num_ads], tmpval[0]);
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
		  strcpy (con_chem_name[num_conditions]
			  [CD->NumSpc + CD->NumAds + num_cex], tmpstr[0]);
		  fprintf (stderr, " cation exchange %s\t %6.4f\n",
			   con_chem_name[num_conditions][CD->NumSpc +
							 CD->NumAds +
							 num_cex], tmpval[0]);
		  num_cex++;
		  // under construction
		}
	      num_species++;
	    }
	  fgets (line, line_width, chemfile);
	}

    }

  int check_conditions_num;

  if (CD->PrpFlg)
    check_conditions_num = num_conditions + 1;
  else
    check_conditions_num = num_conditions;

  if (num_species != CD->NumStc)
    fprintf (stderr, " Number of species does not match indicated value!\n");

  for (i = 1; i < check_conditions_num; i++)
    {
      for (j = 0; j < num_species; j++)
	{
	  if (strcmp (con_chem_name[i][j], con_chem_name[i - 1][j]) != 0)
	    {
	      error_flag = 1;
	    }
	}
      if (error_flag == 1)
	fprintf (stderr,
		 " The order of the chemicals in condition <%s> is incorrect!\n",
		 chemcon[i - 1]);
    }


  for (i = 0; i < CD->NumStc; i++)
    {
      strcpy (CD->chemtype[i].ChemName, con_chem_name[0][i]);
      CD->chemtype[i].itype = Condition_vcele[0].p_type[i];
      fprintf (stderr, " %12s\t%10d\n", CD->chemtype[i].ChemName,
	       CD->chemtype[i].itype);
    }

  if (CD->PrpFlg)
    {
      fprintf (stderr, " Total concentraions in precipitataion:\n");
      for (i = 0; i < CD->NumSpc; i++)
	{
	  if (!strcmp (con_chem_name[num_conditions][i], "pH"))
	    if (CD->Precipitation.t_conc[i] < 7)
	      CD->Precipitation.t_conc[i] =
		pow (10, -CD->Precipitation.t_conc[i]);
	    else
	      CD->Precipitation.t_conc[i] =
		-pow (10, CD->Precipitation.t_conc[i] - 14);
	  // change the pH of precipitation into total concentraion of H
	  // We skip the speciation for rain and assume it is OK to calculate this way.
	  fprintf (stderr, " %12s: %10g M\n",
		   con_chem_name[num_conditions][i],
		   CD->Precipitation.t_conc[i]);
	}
    }


  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "SECONDARY_SPECIES", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      if (keymatch (line, "NULL", tmpval, tmpstr) != 2)
	{
	  strcpy (CD->chemtype[num_species++].ChemName, tmpstr[0]);
	  fprintf (stderr, " %s\n", CD->chemtype[num_species - 1].ChemName);
	}
      fgets (line, line_width, chemfile);
    }
  int num_dep = 2;

  CD->kinetics =
    (Kinetic_Reaction *) malloc (CD->NumMkr * sizeof (Kinetic_Reaction));
  for (i = 0; i < CD->NumMkr; i++)
    {
      CD->kinetics[i].species = (char *) malloc (WORD_WIDTH * sizeof (char));
      CD->kinetics[i].Label = (char *) malloc (WORD_WIDTH * sizeof (char));
      CD->kinetics[i].dep_species =
	(char **) malloc (num_dep * sizeof (char *));
      CD->kinetics[i].dep_power =
	(double *) malloc (num_dep * sizeof (double));
      CD->kinetics[i].monod = (char **) malloc (num_dep * sizeof (char *));
      CD->kinetics[i].monod_para =
	(double *) malloc (num_dep * sizeof (double));
      CD->kinetics[i].inhibition =
	(char **) malloc (num_dep * sizeof (char *));
      CD->kinetics[i].inhibition_para =
	(double *) malloc (num_dep * sizeof (double));
      CD->kinetics[i].dep_position = (int *) malloc (num_dep * sizeof (int));
      for (j = 0; j < num_dep; j++)
	{
	  CD->kinetics[i].dep_species[j] =
	    (char *) malloc (WORD_WIDTH * sizeof (char));
	  CD->kinetics[i].monod[j] =
	    (char *) malloc (WORD_WIDTH * sizeof (char));
	  CD->kinetics[i].inhibition[j] =
	    (char *) malloc (WORD_WIDTH * sizeof (char));
	  CD->kinetics[i].dep_position[j] = 0;
	}
    }
  k = 0;
  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "MINERALS", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      if (keymatch (line, " ", tmpval, tmpstr) != 2)
	{
	  strcpy (CD->kinetics[k].species, tmpstr[0]);
	  if (strcmp (tmpstr[1], "-label") == 0)
	    strcpy (CD->kinetics[k].Label, tmpstr[2]);
	  k++;
	}
      fgets (line, line_width, chemfile);
    }
  for (i = 0; i < k; i++)
    fprintf (stderr, " Kinetic reaction on %s is specified, label %s\n",
	     CD->kinetics[i].species, CD->kinetics[i].Label);

  // start of concentration in precipitation read in
  if (CD->PrpFlg == 2)
    {

      CD->TSD_prepconc = (TSD *) malloc (sizeof (TSD));

      fscanf (prepconc, "%s %d %d", &(CD->TSD_prepconc->name),
	      &(CD->TSD_prepconc->index), &(CD->TSD_prepconc->length));

      CD->prepconcindex =
	(int *) malloc (CD->TSD_prepconc->index * sizeof (int));
      // here prepconc.index is used to save the number of primary species. Must be equal to the number of primary species specified before.
      for (i = 0; i < CD->TSD_prepconc->index; i++)
	{
	  fscanf (prepconc, "%d", &(CD->prepconcindex[i]));
	  if (CD->prepconcindex[i] > 0)
	    {
	      assert (CD->prepconcindex[i] <= CD->NumSpc);
	      fprintf (stderr,
		       " Concentration of %s specified in input file is a time series\n",
		       CD->chemtype[CD->prepconcindex[i] - 1].ChemName);
	    }
	}


      CD->TSD_prepconc->TS =
	(realtype **) malloc ((CD->TSD_prepconc->length) *
			      sizeof (realtype *));
      for (i = 0; i < CD->TSD_prepconc->length; i++)
	{
	  CD->TSD_prepconc->TS[i] =
	    (realtype *) malloc ((1 + CD->TSD_prepconc->index) *
				 sizeof (realtype));
	  fscanf (prepconc, "%d-%d-%d %d:%d:%d", &timeinfo->tm_year,
		  &timeinfo->tm_mon, &timeinfo->tm_mday, &timeinfo->tm_hour,
		  &timeinfo->tm_min, &timeinfo->tm_sec);
	  timeinfo->tm_year = timeinfo->tm_year - 1900;
	  timeinfo->tm_mon = timeinfo->tm_mon - 1;
	  rawtime = timegm (timeinfo);
	  CD->TSD_prepconc->TS[i][0] = (realtype) rawtime;
	  for (j = 0; j < CD->TSD_prepconc->index; j++)
	    {
	      fscanf (prepconc, "%lf", &CD->TSD_prepconc->TS[i][j + 1]);	// [i][0] stores the time
	    }
	  /*
	     fprintf(stderr, " %8.4f\t", CD->TSD_prepconc->TS[i][0] );
	     for ( j = 1 ; j <= CD->TSD_prepconc->index; j ++)
	     fprintf(stderr, " %d %6.4g\t", j , CD->TSD_prepconc->TS[i][j]);
	     fprintf(stderr, "\n");
	   */
	}

      CD->TSD_prepconc->iCounter = 0;
    }


  rewind (chemfile);
  fgets (line, line_width, chemfile);
  while (keymatch (line, "PUMP", tmpval, tmpstr) != 1)
    fgets (line, line_width, chemfile);
  CD->NumPUMP = tmpval[0];
  fprintf (stderr, " %d pumps specified\n", CD->NumPUMP);
  CD->pumps = (Pump *) malloc (CD->NumPUMP * sizeof (Pump));
  i = 0;

  while (keymatch (line, "END", tmpval, tmpstr) != 1)
    {
      fgets (line, line_width, chemfile);
      if (keymatch (line, " ", tmpval, tmpstr) != 2)
	{
	  CD->pumps[i].Pump_Location = (int) tmpval[0];
	  CD->pumps[i].Injection_rate = (double) tmpval[1];
	  CD->pumps[i].Injection_conc = (double) tmpval[2];
	  CD->pumps[i].flow_rate =
	    CD->pumps[i].Injection_rate / CD->pumps[i].Injection_conc / 365 *
	    1E-3;
	  CD->pumps[i].Name_Species = (char *) malloc (20 * sizeof (char));
	  strcpy (CD->pumps[i].Name_Species, tmpstr[1]);
	  //      wrap(CD->pumps[i].Name_Species);
	  CD->pumps[i].Position_Species = -1;
	  for (j = 0; j < CD->NumStc; j++)
	    {
	      if (!strcmp
		  (CD->pumps[i].Name_Species, CD->chemtype[j].ChemName))
		{
		  CD->pumps[i].Position_Species = j;
		}
	    }
	  fprintf (stderr,
		   " -- Rate %f moles/year of %s (%d) at Grid %d with a concentration of %f moles/L\n",
		   CD->pumps[i].Injection_rate, CD->pumps[i].Name_Species,
		   CD->pumps[i].Position_Species, CD->pumps[i].Pump_Location,
		   CD->pumps[i].Injection_conc);
	  fprintf (stderr, " -- Flow rate is then %f cubic meter per day\n",
		   CD->pumps[i].flow_rate);
	  //      CD->pumps[i].Injection_rate *= 1E-3 / 365;
	  i++;
	}
      if (i >= CD->NumPUMP)
	break;
    }


  // End of input file read in

  //  if( num_species != CD->NumSpc + CD->NumSsc) fprintf(stderr, "%d\t%d\t%d\n",num_species, CD->NumSpc, CD->NumSsc);

  /*  switch(error_flag){
     case 1:
     fprintf(stderr, " The order of the chemicals in condition <%s> is incorrect!\n",chemcon[i]);
     break;
     }
   */
  /* 0         ~   NumEle-1       : Triangular GW Element; 
   * NumEle    ~ 2*NumEle-1       : Triangular unsat Element;
   * 2*NumEle  ~ 2*NumEle+NumRiv-1: EBR (Element Beneath River);
   * 2*NumEle+NumRiv ~ 2*NumEle+2*NumRiv-1 : River Element;
   * index = subscript + 1;
   */

  CD->Vcele = (vol_conc *) malloc (CD->NumVol * sizeof (vol_conc));


  /* Initializing volumetric parameters, inherit from pihm
   * That is, if pihm is started from a hot start, rt is also 
   * initialized with the hot data */

  /* Initializing volumetrics for groundwater (GW) cells */

  for (i = 0; i < MD->NumEle; i++)
    {

      CD->Vcele[i].height_v = MD->Ele[i].zmax - MD->Ele[i].zmin;
      CD->Vcele[i].height_o = MD->DummyY[i + 2 * MD->NumEle];
      CD->Vcele[i].height_t = MD->DummyY[i + 2 * MD->NumEle];
      CD->Vcele[i].area = MD->Ele[i].area;
      CD->Vcele[i].porosity = MD->Ele[i].Porosity;
      CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
      CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
      CD->Vcele[i].sat = 1.0;
      CD->Vcele[i].sat_o = 1.0;
      //      CD->Vcele[i].temperature = (realtype) MD->Ele[i].temp;
      CD->Vcele[i].reset_ref = 0;
      for (j = 0; j < 3; j++)
	{
	  if (condition_index[MD->Ele[i].nabr[j]] == condition_index[i + 1])
	    CD->Vcele[i].reset_ref = MD->Ele[i].nabr[j];
	}
    }

  /* Initializing unsaturated zone (vadoze) */

  for (i = MD->NumEle; i < 2 * MD->NumEle; i++)
    {
      CD->Vcele[i].height_v = CD->Vcele[i - MD->NumEle].height_v;
      CD->Vcele[i].height_o = MD->DummyY[i];
      CD->Vcele[i].height_t = MD->DummyY[i];
      CD->Vcele[i].area = MD->Ele[i - MD->NumEle].area;
      CD->Vcele[i].porosity = MD->Ele[i - MD->NumEle].Porosity;
      //    CD->Vcele[i].sat      = MD->DummyY[i]/ (CD->Vcele[i].height_v-CD->Vcele[i-MD->NumEle].height_o) ;
      CD->Vcele[i].sat =
	MD->DummyY[i] / (CD->Vcele[i].height_v -
			 CD->Vcele[i - MD->NumEle].height_o);
      CD->Vcele[i].sat_o = CD->Vcele[i].sat;
      CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
      CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
      //      CD->Vcele[i].temperature = MD->Ele[i - MD->NumEle].temp;
      /* The saturation of unsaturated zone is the Hu divided by height of this cell */
      if (CD->Vcele[i].sat > 1.0)
	fprintf (stderr,
		 "Fatal Error, Unsaturated Zone Initialization For RT Failed!\n");
      CD->Vcele[i].reset_ref = i - MD->NumEle + 1;
      // default reset reference of unsaturated cells are the groundwater cells underneath
    }

  /* Initializing River cells */
  for (i = 2 * MD->NumEle; i < 2 * MD->NumEle + MD->NumRiv; i++)
    {
      j = i - 2 * MD->NumEle;
      CD->Vcele[i].height_v = MD->DummyY[i + MD->NumEle];
      CD->Vcele[i].height_o = CD->Vcele[i].height_v;
      CD->Vcele[i].height_t = CD->Vcele[i].height_o;
      CD->Vcele[i].area =
	MD->Riv[j].Length *
	CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,
			MD->DummyY[j + 3 * MD->NumEle], MD->Riv[j].coeff, 3);
      CD->Vcele[i].porosity = 1.0;
      CD->Vcele[i].sat = 1.0;
      CD->Vcele[i].sat_o = CD->Vcele[i].sat;
      CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
      CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
      CD->Vcele[i].reset_ref = i + MD->NumRiv + 1;
      // default reset reference of river segments are the EBR underneath
    }


  /* Initializing EBR cells */
  for (i = 2 * MD->NumEle + MD->NumRiv; i < 2 * MD->NumEle + 2 * MD->NumRiv;
       i++)
    {
      j = i - 2 * MD->NumEle - MD->NumRiv;
      CD->Vcele[i].height_v = MD->DummyY[i + MD->NumEle];
      CD->Vcele[i].height_o = CD->Vcele[i].height_v;
      CD->Vcele[i].height_t = CD->Vcele[i].height_o;
      CD->Vcele[i].area =
	MD->Riv[j].Length *
	CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,
			MD->DummyY[j + 3 * MD->NumEle], MD->Riv[j].coeff, 3);
      CD->Vcele[i].porosity = 1.0;
      CD->Vcele[i].sat = 1.0;
      CD->Vcele[i].sat_o = CD->Vcele[i].sat;
      CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
      CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
      CD->Vcele[i].reset_ref = 0;
      if (condition_index[MD->Riv[j].LeftEle] == condition_index[i + 1])
	CD->Vcele[i].reset_ref = MD->Riv[j].LeftEle;
      if (condition_index[MD->Riv[j].RightEle] == condition_index[i + 1])
	CD->Vcele[i].reset_ref = MD->Riv[j].RightEle;
    }

  tmpval[0] = 0.0;
  for (i = 0; i < CD->NumEle; i++)
    {
      tmpval[0] += CD->Vcele[i].height_v;
    }
  tmpval[0] = tmpval[0] / CD->NumEle;
  fprintf (stderr, " Average bedrock depth is %f\n", tmpval[0]);



  for (i = 0; i < CD->NumSpc; i++)
    if (strcmp (CD->chemtype[i].ChemName, "pH") == 0)
      {
	strcpy (CD->chemtype[i].ChemName, "H+");
	speciation_flg = 1;
      }

  for (i = CD->NumOsv; i < CD->NumVol; i++)
    {
      j = i - 2 * MD->NumEle - MD->NumRiv;
      CD->Vcele[i].height_v = 1.0;
      CD->Vcele[i].height_o = 1.0;
      CD->Vcele[i].height_t = 1.0;
      CD->Vcele[i].area = 1.0;
      CD->Vcele[i].porosity = 1.0;
      CD->Vcele[i].sat = 1.0;
      CD->Vcele[i].sat_o = 1.0;
      CD->Vcele[i].vol_o = 1.0;
      CD->Vcele[i].vol = 1.0;
      CD->Vcele[i].reset_ref = 2 * MD->NumEle + MD->NumRiv + 1;
    }

  /* Initializing concentration distributions */

  for (i = 0; i < CD->NumVol; i++)
    {
      CD->Vcele[i].index = i + 1;
      CD->Vcele[i].NumStc = CD->NumStc;
      CD->Vcele[i].NumSsc = CD->NumSsc;
      CD->Vcele[i].t_conc = (double *) malloc (CD->NumStc * sizeof (double));
      CD->Vcele[i].t_rate = (double *) malloc (CD->NumStc * sizeof (double));
      CD->Vcele[i].t_tol = (double *) malloc (CD->NumStc * sizeof (double));
      CD->Vcele[i].p_conc = (double *) malloc (CD->NumStc * sizeof (double));
      CD->Vcele[i].s_conc = (double *) malloc (CD->NumSsc * sizeof (double));
      CD->Vcele[i].p_actv = (double *) malloc (CD->NumStc * sizeof (double));
      CD->Vcele[i].s_actv = (double *) malloc (CD->NumSsc * sizeof (double));
      CD->Vcele[i].p_para = (double *) malloc (CD->NumStc * sizeof (double));
      CD->Vcele[i].p_type = (int *) malloc (CD->NumStc * sizeof (int));

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
	  if ((speciation_flg == 1)
	      && (strcmp (CD->chemtype[j].ChemName, "H+") == 0))
	    {
	      CD->Vcele[i].p_conc[j] =
		pow (10,
		     -(Condition_vcele[condition_index[i + 1] - 1].t_conc
		       [j]));
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
	      CD->Vcele[i].t_conc[j] =
		Condition_vcele[condition_index[i + 1] - 1].t_conc[j];
	      CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j] * 0.5;
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

  for (i = 0; i < MD->NumEle; i++)
    {
      for (j = 0; j < 3; j++)
	if (MD->Ele[i].nabr[j] > 0)
	  {
	    num_face++;
	  }
      total_area += MD->Ele[i].area;
    }
  CD->PIHMFac = num_face;

  num_face *= 2;
  num_face += 2 * MD->NumEle + 6 * MD->NumRiv;
  // A river + EBR taks 6 faces

  fprintf (stderr, " Total area of the watershed is %f m^2\n", total_area);

  CD->NumFac = num_face;

  /* Configure the lateral connectivity of gw grid blocks */

  CD->Flux = (face *) malloc (CD->NumFac * sizeof (face));
  k = 0;

  double dist1, dist2, para_a, para_b, para_c, x_0, x_1, y_0, y_1;
  int index_0, index_1, rivi, control;

  for (i = 0; i < MD->NumEle; i++)
    {
      for (j = 0; j < 3; j++)
	if (MD->Ele[i].nabr[j] > 0)
	  {
	    /* node indicates the index of grid blocks, not nodes at corners */
	    CD->Flux[k].nodeup = i + 1;
	    CD->Flux[k].nodelo = MD->Ele[i].nabr[j];
	    CD->Flux[k].nodeuu =
	      upstream (MD->Ele[i], MD->Ele[MD->Ele[i].nabr[j] - 1], MD);
	    CD->Flux[k].nodell =
	      upstream (MD->Ele[MD->Ele[i].nabr[j] - 1], MD->Ele[i], MD);
	    CD->Flux[k].flux_type = 0;
	    CD->Flux[k].BC = 0;
	    if (CD->Flux[k].nodeuu > 0)
	      CD->Flux[k].distuu =
		sqrt (pow
		      (MD->Ele[i].x - MD->Ele[CD->Flux[k].nodeuu - 1].x,
		       2) + pow (MD->Ele[i].y - MD->Ele[CD->Flux[k].nodeuu -
							1].y, 2));
	    else
	      CD->Flux[k].distuu = 0.0;
	    if (CD->Flux[k].nodell > 0)
	      CD->Flux[k].distll =
		sqrt (pow
		      (MD->Ele[MD->Ele[i].nabr[j] - 1].x -
		       MD->Ele[CD->Flux[k].nodell - 1].x,
		       2) + pow (MD->Ele[MD->Ele[i].nabr[j] - 1].y -
				 MD->Ele[CD->Flux[k].nodell - 1].y, 2));
	    else
	      CD->Flux[k].distll = 0.0;

	    CD->Flux[k].distance =
	      sqrt (pow (MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x, 2)
		    + pow (MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y,
			   2));


	    index_0 = MD->Ele[i].node[(j + 1) % 3] - 1;
	    index_1 = MD->Ele[i].node[(j + 2) % 3] - 1;

	    x_0 = MD->Node[index_0].x;
	    y_0 = MD->Node[index_0].y;
	    x_1 = MD->Node[index_1].x;
	    y_1 = MD->Node[index_1].y;

	    para_a = y_1 - y_0;
	    para_b = x_0 - x_1;
	    para_c = (x_1 - x_0) * y_0 - (y_1 - y_0) * x_0;

	    if ((para_a * MD->Ele[i].x + para_b * MD->Ele[i].y +
		 para_c) * (para_a * MD->Ele[MD->Ele[i].nabr[j] - 1].x +
			    para_b * MD->Ele[MD->Ele[i].nabr[j] - 1].y +
			    para_c) > 0.0)
	      fprintf (stderr, " two points at the same side of edge!\n");
	    dist1 =
	      fabs (para_a * MD->Ele[i].x + para_b * MD->Ele[i].y +
		    para_c) / sqrt (para_a * para_a + para_b * para_b);
	    dist2 =
	      fabs (para_a * MD->Ele[MD->Ele[i].nabr[j] - 1].x +
		    para_b * MD->Ele[MD->Ele[i].nabr[j] - 1].y +
		    para_c) / sqrt (para_a * para_a + para_b * para_b);
	    if ((CD->Flux[k].distance < dist1)
		|| (CD->Flux[k].distance < dist2))
	      {
		fprintf (stderr,
			 "\n Checking the distance calculation wrong, total is %f, from ele to edge is %f, from edge to neighbor is %f\n",
			 CD->Flux[k].distance, dist1, dist2);
		fprintf (stderr,
			 " Ele coordinates, x: %f, y: %f\n Node_1 coordinates, x: %f, y: %f\n Node_2 coordinates, x: %f, y: %f\n Nabr x:%f, y:%f\n",
			 MD->Ele[i].x, MD->Ele[i].y, x_0, y_0, x_1, y_1,
			 MD->Ele[MD->Ele[i].nabr[j] - 1].x,
			 MD->Ele[MD->Ele[i].nabr[j] - 1].y);
		fprintf (stderr,
			 " node_1 x:%f y:%f\t node_2 x:%f y:%f node_3 x:%f y:%f nabr x:%f y:%f, no:%d\n",
			 MD->Node[MD->Ele[i].node[0] - 1].x,
			 MD->Node[MD->Ele[i].node[0] - 1].y,
			 MD->Node[MD->Ele[i].node[1] - 1].x,
			 MD->Node[MD->Ele[i].node[1] - 1].y,
			 MD->Node[MD->Ele[i].node[2] - 1].x,
			 MD->Node[MD->Ele[i].node[2] - 1].y,
			 MD->Ele[MD->Ele[i].nabr[j] - 1].x,
			 MD->Ele[MD->Ele[i].nabr[j] - 1].y, j);
	      }
	    else
	      {
		CD->Flux[k].distance = dist1;
	      }
	    for (rivi = 0; rivi < CD->NumRiv; rivi++)
	      {
		if ((CD->Flux[k].nodeup == MD->Riv[rivi].LeftEle)
		    && (CD->Flux[k].nodelo == MD->Riv[rivi].RightEle))
		  {
		    CD->Flux[k].nodelo =
		      2 * CD->NumEle + CD->NumRiv + rivi + 1;
		    CD->Flux[k].nodeuu = 0;
		    CD->Flux[k].nodell = 0;
		    CD->Flux[k].distance = dist1;
		    //      fprintf(stderr, " Riv - Watershed connection corrected, now from cell %d to cell %d\n", CD->Flux[k].nodeup, CD->Flux[k].nodelo);
		  }
		if ((CD->Flux[k].nodeup == MD->Riv[rivi].RightEle)
		    && (CD->Flux[k].nodelo == MD->Riv[rivi].LeftEle))
		  {
		    CD->Flux[k].nodelo =
		      2 * CD->NumEle + CD->NumRiv + rivi + 1;
		    CD->Flux[k].nodeuu = 0;
		    CD->Flux[k].nodell = 0;
		    CD->Flux[k].distance = dist1;
		    //            fprintf(stderr, " Riv - Watershed connection corrected, now from cell %d to cell %d\n", CD->Flux[k].nodeup, CD->Flux[k].nodelo);
		  }
	      }
	    k++;
	  }
    }

  for (i = 0; i < MD->NumEle; i++)
    {
      for (j = 0; j < 3; j++)
	if (MD->Ele[i].nabr[j] > 0)
	  {
	    /* node indicates the index of grid blocks, not nodes at corners */
	    CD->Flux[k].nodeup = i + 1 + MD->NumEle;
	    CD->Flux[k].nodelo = MD->Ele[i].nabr[j] + MD->NumEle;
	    CD->Flux[k].nodeuu =
	      upstream (MD->Ele[i], MD->Ele[MD->Ele[i].nabr[j] - 1],
			MD) + MD->NumEle;
	    CD->Flux[k].nodell =
	      upstream (MD->Ele[MD->Ele[i].nabr[j] - 1], MD->Ele[i],
			MD) + MD->NumEle;
	    CD->Flux[k].flux_type = 0;
	    CD->Flux[k].BC = 0;
	    if (CD->Flux[k].nodeuu > 0)
	      CD->Flux[k].distuu =
		sqrt (pow
		      (MD->Ele[i].x - MD->Ele[CD->Flux[k].nodeuu - 1].x,
		       2) + pow (MD->Ele[i].y - MD->Ele[CD->Flux[k].nodeuu -
							1].y, 2));
	    else
	      CD->Flux[k].distuu = 0.0;
	    if (CD->Flux[k].nodell > 0)
	      CD->Flux[k].distll =
		sqrt (pow
		      (MD->Ele[MD->Ele[i].nabr[j] - 1].x -
		       MD->Ele[CD->Flux[k].nodell - 1].x,
		       2) + pow (MD->Ele[MD->Ele[i].nabr[j] - 1].y -
				 MD->Ele[CD->Flux[k].nodell - 1].y, 2));
	    else
	      CD->Flux[k].distll = 0.0;
	    CD->Flux[k].distance =
	      sqrt (pow (MD->Ele[i].x - MD->Ele[MD->Ele[i].nabr[j] - 1].x, 2)
		    + pow (MD->Ele[i].y - MD->Ele[MD->Ele[i].nabr[j] - 1].y,
			   2));
	    index_0 = MD->Ele[i].node[(j + 1) % 3] - 1;
	    index_1 = MD->Ele[i].node[(j + 2) % 3] - 1;
	    x_0 = MD->Node[index_0].x;
	    y_0 = MD->Node[index_0].y;
	    x_1 = MD->Node[index_1].x;
	    y_1 = MD->Node[index_1].y;
	    para_a = y_1 - y_0;
	    para_b = x_0 - x_1;
	    para_c = (x_1 - x_0) * y_0 - (y_1 - y_0) * x_0;
	    dist1 =
	      fabs (para_a * MD->Ele[i].x + para_b * MD->Ele[i].y +
		    para_c) / sqrt (para_a * para_a + para_b * para_b);
	    dist2 =
	      fabs (para_a * MD->Ele[MD->Ele[i].nabr[j] - 1].x +
		    para_b * MD->Ele[MD->Ele[i].nabr[j] - 1].y +
		    para_c) / sqrt (para_a * para_a + para_b * para_b);
	    CD->Flux[k].distance = dist1;
	    for (rivi = 0; rivi < CD->NumRiv; rivi++)
	      {
		if ((CD->Flux[k].nodeup - CD->NumEle == MD->Riv[rivi].LeftEle)
		    && (CD->Flux[k].nodelo - CD->NumEle ==
			MD->Riv[rivi].RightEle))
		  {
		    CD->Flux[k].nodelo =
		      2 * CD->NumEle + CD->NumRiv + rivi + 1;
		    CD->Flux[k].nodeuu = 0;
		    CD->Flux[k].nodell = 0;
		    CD->Flux[k].distance = dist1;
		  }
		if ((CD->Flux[k].nodeup - CD->NumEle ==
		     MD->Riv[rivi].RightEle)
		    && (CD->Flux[k].nodelo - CD->NumEle ==
			MD->Riv[rivi].LeftEle))
		  {
		    CD->Flux[k].nodelo =
		      2 * CD->NumEle + CD->NumRiv + rivi + 1;
		    CD->Flux[k].nodeuu = 0;
		    CD->Flux[k].nodell = 0;
		    CD->Flux[k].distance = dist1;
		  }
	      }
	    k++;
	  }
    }

  /* Configure the connectivity of unsat - gw grid blocks follows */

  /* centered at unsat blocks */

  for (i = MD->NumEle; i < 2 * MD->NumEle; i++)
    {
      CD->Vcele[i].ErrDumper = k;
      CD->Flux[k].nodeup = i + 1;
      CD->Flux[k].nodelo = i + 1 - MD->NumEle;
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 0;
      CD->Flux[k].distance = 0.1 * CD->Vcele[i].height_v;
      k++;
    }

  /* centered at gw blocks */

  for (i = 0; i < MD->NumEle; i++)
    {
      CD->Vcele[i].ErrDumper = k;
      CD->Flux[k].nodeup = i + 1;
      CD->Flux[k].nodelo = i + 1 + MD->NumEle;
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 1;
      CD->Flux[k].BC = 0;
      CD->Flux[k].distance = 0.1 * CD->Vcele[i].height_v;
      k++;
    }
  CD->NumDis = k;

  /* Between River Segment */

  // From upstream 0
  /*
     for ( i = 0; i <MD->NumRiv; i++){
     if ( i == 0){
     CD->Flux[k].nodelo = CD->NumVol - 1;
     CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1;
     }
     else{
     CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1;
     CD->Flux[k].nodelo = i + 2 * MD->NumEle;
     }
     CD->Flux[k].nodeuu = 0;
     CD->Flux[k].nodell = 0;
     CD->Flux[k].flux_type = 0;
     CD->Flux[k].BC     = 1;
     CD->Flux[k].flux   = 0.0;
     CD->Flux[k].distance = MD->Riv[i].Length * 0.5;
     CD->Flux[k].s_area   = 0.0;
     k++;
     }

     // To downstream 1

     for ( i = 0; i <MD->NumRiv; i++){
     if ( i == MD->NumRiv - 1){
     CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1;
     CD->Flux[k].nodelo = CD->NumVol - 1;
     }
     else{
     CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1;
     CD->Flux[k].nodelo = i + 2 * MD->NumEle + 2;
     }
     CD->Flux[k].nodeuu = 0;
     CD->Flux[k].nodell = 0;
     CD->Flux[k].flux_type = 0;
     CD->Flux[k].BC     = 1;
     CD->Flux[k].flux   = 0.0;
     CD->Flux[k].distance = MD->Riv[i].Length * 0.5;
     CD->Flux[k].s_area   = 0.0;
     k++;
     }
   */

  /* Between River and Left */
  // River to left OFL 2
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1 + MD->NumRiv;
      CD->Flux[k].nodelo = CD->NumVol;
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 1;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].distance = 1.0;
      CD->Flux[k].s_area = 0.0;
      k++;
    }

  /* Between River and Right */
  // River to right OFL 3
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1 + MD->NumRiv;
      CD->Flux[k].nodelo = CD->NumVol;
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 1;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].distance = 1.0;
      CD->Flux[k].s_area = 0.0;
      k++;
    }

  /* Between Riv and EBR */
  /*
     // From river to EBR 6
     for ( i = 0; i <MD->NumRiv; i++){
     CD->Flux[k].nodeup = i + 2 * MD->NumEle + 1;
     CD->Flux[k].nodelo = i + 2 * MD->NumEle + MD->NumRiv + 1;
     CD->Flux[k].nodeuu = 0;
     CD->Flux[k].nodell = 0;
     CD->Flux[k].flux_type = 0;
     CD->Flux[k].BC     = 1;
     CD->Flux[k].flux   = 0.0;
     CD->Flux[k].distance = 1.0;
     CD->Flux[k].s_area   = 0.0;
     k++;
     }

     // From EBR to river -6
     for ( i = 0; i <MD->NumRiv; i++){
     CD->Flux[k].nodelo = i + 2 * MD->NumEle + 1;
     CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
     CD->Flux[k].nodeuu = 0;
     CD->Flux[k].nodell = 0;
     CD->Flux[k].flux_type = 0;
     CD->Flux[k].BC     = 1;
     CD->Flux[k].flux   = 0.0;
     CD->Flux[k].distance = 1.0;
     CD->Flux[k].s_area   = 0.0;
     k++;
     }
   */
  /* Between Left and EBR */
  // EBR to left  7
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
      CD->Flux[k].nodelo = MD->Riv[i].LeftEle;
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 0;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].distance = 1.0;
      CD->Flux[k].s_area =
	MD->Riv[i].Length * CD->Vcele[MD->Riv[i].LeftEle - 1].height_v;
      k++;
    }


  /* Between Right and EBR */
  // EBR to right 8
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
      CD->Flux[k].nodelo = MD->Riv[i].RightEle;
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 0;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].distance = 1.0;
      CD->Flux[k].s_area =
	MD->Riv[i].Length * CD->Vcele[MD->Riv[i].RightEle - 1].height_v;
      k++;
    }

  /* Between EBR */
  // To downstream EBR 9
  for (i = 0; i < MD->NumRiv; i++)
    {
      if (i == MD->NumRiv - 1)
	{
	  CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
	  CD->Flux[k].nodelo = CD->NumVol;
	}
      else
	{
	  CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
	  CD->Flux[k].nodelo = i + 2 * MD->NumEle + MD->NumRiv + 2;
	}
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 1;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].distance = 1.0;
      CD->Flux[k].s_area = 0.0;
      k++;
    }

  // From upstream EBR 10
  for (i = 0; i < MD->NumRiv; i++)
    {
      if (i == 0)
	{
	  CD->Flux[k].nodelo = CD->NumVol;
	  CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
	}
      else
	{
	  CD->Flux[k].nodeup = i + 2 * MD->NumEle + MD->NumRiv + 1;
	  CD->Flux[k].nodelo = i + 2 * MD->NumEle + MD->NumRiv;
	}
      CD->Flux[k].nodeuu = 0;
      CD->Flux[k].nodell = 0;
      CD->Flux[k].flux_type = 0;
      CD->Flux[k].BC = 1;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].distance = 1.0;
      CD->Flux[k].s_area = 0.0;
      k++;
    }

  for (k = 0; k < CD->NumFac; k++)
    {
      CD->Flux[k].velocity = 0.0;
      CD->Flux[k].flux = 0.0;
      CD->Flux[k].s_area = 0.0;
    }


  CD->SPCFlg = speciation_flg;
  Lookup (database, CD, 0);
  // update the concentration of mineral after get the molar volume of mineral

  for (i = 0; i < CD->NumAkr + CD->NumMkr; i++)
    {
      if (!strcmp
	  (CD->chemtype[i + CD->NumSpc + CD->NumAds + CD->NumCex].ChemName,
	   "'CO2(*g)'"))
	{
	  CD->KeqKinect[i] += log10 (CS->Cal.PCO2);
	}
      else
	{
	  CD->KeqKinect[i] += log10 (CS->Cal.Keq);
	}
    }



  for (i = 0; i < CD->NumStc; i++)
    {
      fprintf (stderr, " Conc and SSA of each species: %6.4g %6.4g\n",
	       CD->Vcele[200].t_conc[i], CD->Vcele[200].p_para[i]);
    }

  fprintf (stderr, " Kinetic Mass Matrx!\n\t\t");
  for (i = 0; i < CD->NumStc; i++)
    fprintf (stderr, " %6s\t", CD->chemtype[i].ChemName);
  fprintf (stderr, "\n");
  for (j = 0; j < CD->NumMkr + CD->NumAkr; j++)
    {
      fprintf (stderr, " %6s\t",
	       CD->chemtype[j + CD->NumSpc + CD->NumAds +
			    CD->NumCex].ChemName);
      for (i = 0; i < CD->NumStc; i++)
	{
	  fprintf (stderr, " %6.2f\t", CD->Dep_kinetic[j][i]);
	}
      fprintf (stderr, "\tKeq = %6.2f\n", CD->KeqKinect[j]);
    }
  // use calibration coefficient to produce new Keq values for 1) CO2, 2) other kinetic reaction



  fprintf (stderr,
	   " Mass action species type determination (0: immobile, 1: mobile, 2: Mixed)\n");
  for (i = 0; i < CD->NumSpc; i++)
    {
      if (CD->chemtype[i].itype == 1)
	CD->chemtype[i].mtype = 1;
      else
	CD->chemtype[i].mtype = 0;
      for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
	{
	  if ((CD->Totalconc[i][j] != 0)
	      && (CD->chemtype[j].itype != CD->chemtype[i].mtype))
	    CD->chemtype[i].mtype = 2;
	}
      /*
         if (strcmp( CD->chemtype[i].ChemName, "'H+'") == 0)
         CD->chemtype[i].mtype = 1;
       */
      fprintf (stderr, " %12s\t%10d\n", CD->chemtype[i].ChemName,
	       CD->chemtype[i].mtype);
    }


  fprintf (stderr,
	   " Individual species type determination (1: aqueous, 2: adsorption, 3: ion exchange, 4: solid)\n");
  for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {
      fprintf (stderr, " %12s\t%10d\n", CD->chemtype[i].ChemName,
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
					      INFTYSMALL) * 1000 /
		    CD->chemtype[j].MolarVolume / CD->Vcele[i].porosity;
		  CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
		  // relative mineral volume fraction
		  // porosity can be 1.0 so the relative fraction option need a small modification 
		}
	    }
	  if (CD->chemtype[j].itype == 3)
	    {
	      CD->Vcele[i].t_conc[j] =
		CD->Vcele[i].t_conc[j] * (1 / CD->Vcele[i].porosity -
					  1) * 2650;
	      CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
	      //      fprintf(stderr, "site concentration is %f mole per aqueous L\n",CD->Vcele[i].t_conc[j]);
	      // change the unit of CEC (eq/g) into C(ion site) (eq/L water), assuming density of solid is always 2650 g/L solid
	    }
	}

    }





  CD->SPCFlg = 1;
  if (!CD->RecFlg)
    {

      for (i = 0; i < CD->NumEle; i++)
	Speciation (CD, i);
    }
  CD->SPCFlg = 0;

  for (i = CD->NumEle; i < CD->NumEle * 2; i++)
    for (k = 0; k < CD->NumStc; k++)
      {
	CD->Vcele[i].t_conc[k] = CD->Vcele[i - CD->NumEle].t_conc[k];
	CD->Vcele[i].p_conc[k] = CD->Vcele[i - CD->NumEle].p_conc[k];
	CD->Vcele[i].p_actv[k] = CD->Vcele[i - CD->NumEle].p_actv[k];
	//      fprintf(stderr, " %d %s: %g\n", i, CD->chemtype[k].ChemName, CD->Vcele[i].t_conc[k]);
      }

  for (i = CD->NumEle * 2; i < CD->NumVol; i++)
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
     #pragma omp parallel num_threads(2)
     {
     int tid, nthreads;

     tid = omp_get_thread_num();
     nthreads = omp_get_num_threads();

     printf("Hello world from thread %3d of %3d\n", tid, nthreads);
     }
   */

  CD->TimLst = 0;

  //  fprintf(stderr, " Total primary species %d\n", CD->NumSpc);

  for (i = 0; i < CD->PIHMFac * 2; i++)
    for (j = CD->PIHMFac * 2; j < CD->NumFac; j++)
      {
	if ((CD->Flux[i].nodeup == CD->Flux[j].nodelo)
	    && (CD->Flux[i].nodelo == CD->Flux[j].nodeup))
	  {
	    //      fprintf (stderr, " Flux between %d and %d identified\n",
	    //               CD->Flux[i].nodelo, CD->Flux[i].nodeup);
	    CD->Flux[j].BC = CD->Flux[i].BC = 0;
	    CD->Flux[j].distance = CD->Flux[i].distance;
	    //   fprintf(stderr, " Flux from river to ele corrected as BC %d, distance %f.\n",CD->Flux[j].BC, CD->Flux[j].distance);
	  }
      }


  for (i = 0; i < num_conditions; i++)
    {
      free (Condition_vcele[i].t_conc);
      free (Condition_vcele[i].p_conc);
      free (Condition_vcele[i].p_para);
      free (Condition_vcele[i].p_type);
    }

  free (Condition_vcele);


  for (i = 0; i < num_conditions; i++)
    free (chemcon[i]);
  free (chemcon);


  for (i = 0; i < num_conditions + 1; i++)
    {
      for (j = 0; j < CD->NumStc; j++)
	free (con_chem_name[i][j]);
      free (con_chem_name[i]);
    }
  free (con_chem_name);


  free (chemfn);
  free (datafn);
  free (forcfn);
  free (condition_index);


  fclose (chemfile);
  fclose (database);
  fclose (prepconc);

}

void
fluxtrans (realtype t, realtype stepsize, const void *DS, Chem_Data CD,
	   N_Vector VY)
{

  /* unit of t and stepsize: min */
  /* DS: model data              */

  int i, j, k, NumEle, NumVol;
  struct tm *timestamp;
  time_t *rawtime;
  double timelps, rt_step, temp_rt_step, peclet, tmptime, tmpval, invavg,
    unit_c, tmpconc, timer1, timer2;
  realtype *Y;
  realtype MF_CONVERT;


  NumEle = CD->NumEle;
  NumVol = CD->NumVol;
  k = 0;
  MF_CONVERT = (realtype) (24 * 60 * 60);
  rt_step = stepsize * (double) CD->AvgScl;	// by default, the largest averaging period is per 10 mins. Longer default averaging value will fail
  Y = NV_DATA_S (VY);
  rawtime = (time_t *) malloc (sizeof (time_t));
  *rawtime = (int) (t * 60);
  timestamp = gmtime (rawtime);
  timelps = (t - CD->StartTime);
  Model_Data MD;
  MD = (Model_Data) DS;


  for (i = NumEle; i < 2 * NumEle; i++)
    {
      j = i - NumEle;
      CD->Vcele[i].q +=
	MAX (MD->ElePrep[j], 0.0) * CD->Vcele[i].area * MF_CONVERT;
    }


  /* update flux for the subsurface lateral flow */

  /*
     FILE * btcv = fopen("shp.btcv2","a+");

     fprintf(btcv, " 1071 %6.4g\n", CD->Vcele[1070].p_conc[2]);

     fclose(btcv);

   */

  for (i = 0; i < NumEle; i++)
    {
      CD->Vcele[i].height_tl = Y[i + 2 * NumEle];
      CD->Vcele[i + NumEle].height_tl = Y[i + NumEle];

      for (j = 0; j < 3; j++)
	if (MD->Ele[i].nabr[j] > 0)
	  {
	    /* node indicates the index of grid blocks, not nodes at corners */
	    CD->Flux[k].flux += MD->FluxSub[i][j] * MF_CONVERT;
	    CD->Flux[k].s_area += MD->AreaSub[i][j];
	    CD->Flux[k].velocity += 0;
	    k++;
	  }

    }

  /* connection between unsaturated cells */

  for (i = 0; i < NumEle; i++)
    {
      for (j = 0; j < 3; j++)
	if (MD->Ele[i].nabr[j] > 0)
	  {
	    CD->Flux[k].flux = 0.0;
	    CD->Flux[k].s_area = 1.0;
	    CD->Flux[k].velocity = 0.0;
	    k++;
	  }

    }

  /* Configure the connectivity of unsat - gw grid blocks follows */


  for (i = NumEle; i < 2 * NumEle; i++)
    {
      CD->Flux[k].velocity += MD->Recharge[i - NumEle] * MF_CONVERT;
      CD->Flux[k].flux +=
	MD->Recharge[i - NumEle] * MF_CONVERT * CD->Vcele[i].area;
      CD->Flux[k].s_area += CD->Vcele[i].area;
      k++;
    }

  for (i = 0; i < NumEle; i++)
    {
      CD->Flux[k].velocity += -MD->Recharge[i] * MF_CONVERT;
      CD->Flux[k].flux += -MD->Recharge[i] * MF_CONVERT * CD->Vcele[i].area;
      CD->Flux[k].s_area += CD->Vcele[i].area;
      k++;
    }


  /* Update the fluxes for river related block */
  /*  
     for ( i = 0; i <MD->NumRiv; i++){
     //    CD->Flux[k].flux     += MD->FluxRiv[i][0];
     //    CD->Flux[k].velocity += CD->Flux[k].flux / CD->Vcele[CD->Flux[k].nodelo-1].area;
     k++;
     }
     for ( i = 0; i <MD->NumRiv; i++){
     //    CD->Flux[k].flux     += MD->FluxRiv[i][1];
     //    CD->Flux[k].velocity += CD->Flux[k].flux / CD->Vcele[CD->Flux[k].nodeup-1].area;
     k++;
     }
   */

  //  CD->rivd = MD->FluxRiv[MD->NumRiv - 1][1];

  //  fprintf(stderr, " Rivd is %f this day\n", CD->rivd);

  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].flux += MD->FluxRiv[i][2] * MF_CONVERT;
      k++;
    }
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].flux += MD->FluxRiv[i][3] * MF_CONVERT;
      k++;
    }
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].flux +=
	(MD->FluxRiv[i][7] + MD->FluxRiv[i][4]) * MF_CONVERT;
      k++;
    }
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].flux +=
	(MD->FluxRiv[i][8] + MD->FluxRiv[i][5]) * MF_CONVERT;
      k++;
    }
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].flux +=
	(MD->FluxRiv[i][9] + MD->FluxRiv[i][1]) * MF_CONVERT;
      k++;
    }
  for (i = 0; i < MD->NumRiv; i++)
    {
      CD->Flux[k].flux +=
	(MD->FluxRiv[i][10] + MD->FluxRiv[i][0]) * MF_CONVERT;
      k++;
    }
  /* update the cell volumetrics every averaging cycle */
  if (((int) timelps - (int) CD->TimLst) >= CD->AvgScl * stepsize)
    {
      //    fprintf(stderr,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\n",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
      //    fprintf(stderr, "%6.4f %6.4f", (realtype) * rawtime, (realtype) * rawtime/24/60/60);
      // update the concentration in precipitation here.
      if (CD->PrpFlg == 2)
	{
	  while (CD->TSD_prepconc->TS[CD->TSD_prepconc->iCounter][0] <
		 (realtype) * rawtime)
	    CD->TSD_prepconc->iCounter++;
	  CD->TSD_prepconc->iCounter--;
	  //      fprintf(stderr, " %d %6.4f time %6.4f\n", CD->TSD_prepconc->iCounter, CD->TSD_prepconc->TS[CD->TSD_prepconc->iCounter][0], (realtype) * rawtime);
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
		      fprintf (stderr,
			       " %s in precipitation is changed to %6.4g\n",
			       CD->chemtype[j].ChemName,
			       CD->Precipitation.t_conc[j]);
		    }
		}
	    }
	}

      invavg = stepsize / (double) CD->AvgScl;

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
	  // velocity is calculated according to the flux_avg and the area_avg. 
	}			// for gw cells, contact area is needed for dispersion;                                                                                                                     

      for (i = 0; i < CD->PIHMFac * 2; i++)
	for (j = CD->PIHMFac * 2; j < CD->NumFac; j++)
	  {
	    if ((CD->Flux[i].nodeup == CD->Flux[j].nodelo)
		&& (CD->Flux[i].nodelo == CD->Flux[j].nodeup))
	      {
		CD->Flux[j].s_area = CD->Flux[i].s_area;
		CD->Flux[j].velocity = -CD->Flux[i].velocity;
		//   fprintf(stderr, " Flux from river to ele corrected as BC %d, distance %f.\n",CD->Flux[j].BC, CD->Flux[j].distance);
	      }
	  }

      for (i = 0; i < NumEle; i++)
	{
	  CD->Vcele[i].height_o = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_t = MAX (Y[i + 2 * NumEle], 1.0E-5);
	  CD->Vcele[i].height_int = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_sp =
	    (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
	  CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
	  CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
	}

      /* update the unsaturated zone (vadoze) */
      for (i = NumEle; i < 2 * NumEle; i++)
	{
	  j = i - NumEle;
	  CD->Vcele[i].height_o = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_t = MAX (Y[i], 1.0E-5);
	  CD->Vcele[i].height_int = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_sp =
	    (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
	  CD->Vcele[i].sat_o = CD->Vcele[i].sat;
	  CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
	  CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
	  //      CD->Vcele[i].sat      = CD->Vcele[i].height_t/ (CD->Vcele[i].height_v - CD->Vcele[i-NumEle].height_t);
	  CD->Vcele[i].sat =
	    CD->Vcele[i].height_t / (CD->Vcele[i].height_v -
				     CD->Vcele[i - NumEle].height_t);
	  CD->Vcele[i].q *= invavg;
	}
      /* update river cells */
      for (i = 2 * NumEle; i < 2 * NumEle + MD->NumRiv; i++)
	{
	  j = i - 2 * NumEle;
	  CD->Vcele[i].height_o = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_t = MAX (Y[i + NumEle], 1.0E-5);
	  CD->Vcele[i].height_int = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_sp =
	    (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
	  CD->Vcele[i].area =
	    MD->Riv[j].Length *
	    CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,
			    Y[j + 3 * NumEle], MD->Riv[j].coeff, 3);
	  CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
	  CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
	}
      /* update EBR cells */
      for (i = 2 * NumEle + MD->NumRiv; i < 2 * (NumEle + MD->NumRiv); i++)
	{
	  j = i - 2 * NumEle - MD->NumRiv;
	  CD->Vcele[i].height_o = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_t =
	    MAX (Y[i + NumEle],
		 1.0E-5) + MAX (Y[i + NumEle - MD->NumRiv],
				1.0E-5) / CD->Vcele[i].porosity;
	  CD->Vcele[i].height_int = CD->Vcele[i].height_t;
	  CD->Vcele[i].height_sp =
	    (CD->Vcele[i].height_t - CD->Vcele[i].height_o) * invavg;
	  CD->Vcele[i].area =
	    MD->Riv[j].Length *
	    CS_AreaOrPerem (MD->Riv_Shape[MD->Riv[j].shape - 1].interpOrd,
			    Y[j + 3 * NumEle], MD->Riv[j].coeff, 3);
	  CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
	  CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
	}

      Monitor (t, stepsize * (double) CD->AvgScl, MD, CD);

      for (k = 0; k < CD->NumStc; k++)
	{
	  CD->Vcele[NumVol - 1].t_conc[k] =
	    CD->Precipitation.t_conc[k] * CD->Condensation;
	  CD->Vcele[NumVol - 1].p_conc[k] =
	    CD->Precipitation.t_conc[k] * CD->Condensation;
	}


      /*
         for ( i = CD->PIHMFac ; i < CD->PIHMFac * 2; i ++)
         for ( j = CD->PIHMFac ; j < CD->PIHMFac * 2; j ++)
         {
         if( ( CD->Flux[i].nodeup == CD->Flux[j].nodelo ) && ( CD->Flux[i].nodelo == CD->Flux[j].nodeup))
         {
         if (fabs( CD->Flux[i].flux + CD->Flux[j].flux )> 1.0E-5 ) {
         fprintf(stderr, " Dismatch %d %d %f %f\n", CD->Flux[i].nodeup, CD->Flux[i].nodelo, CD->Flux[i].flux, CD->Flux[j].flux);
         fprintf(stderr, " Corresponding %d %d %f %f\n", CD->Flux[i-CD->PIHMFac].nodeup, CD->Flux[i-CD->PIHMFac].nodelo, CD->Flux[i-CD->PIHMFac].flux, CD->Flux[j-CD->PIHMFac].flux);
         }
         }
         }

       */

      for (i = 0; i < NumVol; i++)
	CD->Vcele[i].rt_step = 0.0;

      for (i = 0; i < CD->NumDis; i++)
	for (j = 0; j < CD->NumSpc; j++)
	  {
	    peclet =
	      fabs (CD->Flux[i].velocity * CD->Flux[i].distance /
		    (CD->chemtype[j].DiffCoe +
		     CD->chemtype[j].DispCoe * CD->Flux[i].velocity));
	    peclet = MAX (peclet, 1.0E-8);
	    if (i < CD->NumDis)
	      {
		temp_rt_step =
		  fabs (CD->Flux[i].flux /
			CD->Vcele[CD->Flux[i].nodeup - 1].vol) * (1 +
								  peclet) /
		  peclet;
	      }
	    else
	      temp_rt_step =
		fabs (CD->Flux[i].flux /
		      CD->Vcele[CD->Flux[i].nodeup - 1].vol);
	    CD->Vcele[CD->Flux[i].nodeup - 1].rt_step += temp_rt_step;
	  }
      /*
         for ( i = 0 ; i < NumEle; i ++)
         if ( MD->ElePrep[i] > 1.0E-8)
         fprintf(stderr, " Ratio %f from %f/%f or %f %f %f\n", MD->EleET[i][2]/MD->EleNetPrep[i], MD->EleET[i][2], MD->EleNetPrep[i], MD->ElePrep[i], MD->EleET[i][1], MD->EleViR[i]);
       */

      //    FILE * rt_dist = fopen("logfile/rtdist.out","w");

      k = 0;			// used to count the number of slow cells.
      invavg = 0.0;
      for (i = 0; i < CD->NumOsv; i++)
	{
	  CD->Vcele[i].rt_step = 0.6 * UNIT_C / CD->Vcele[i].rt_step;
	  //      fprintf(rt_dist, " %d\t %f\n", i+1, CD->Vcele[i].rt_step);
	  if (rt_step > CD->Vcele[i].rt_step)
	    rt_step = CD->Vcele[i].rt_step;

	  if (CD->Vcele[i].rt_step >= (double) CD->AvgScl)
	    CD->Vcele[i].rt_step = (double) CD->AvgScl;
	  else
	    {
	      k++;
	    }
	  invavg += CD->Vcele[i].rt_step;
	}

      invavg = invavg / (double) CD->NumOsv;
      //    fprintf(rt_dist, " Number of slow cell is %d\n ",k);
      //    fclose(rt_dist);

      //      fprintf(stderr, "TL: %6.4f, DL: %6.4f\n",CD->TimLst, (double) CD->Delay);

      if (CD->TimLst >= (double) CD->Delay)
	{

	  //      fprintf (stderr," Calling RT from %f to %f, mean timestep = %f\n", CD->TimLst, stepsize * (double)CD->AvgScl + CD->TimLst, invavg);
	  rt_step = stepsize * (double) CD->AvgScl;
	  AdptTime (CD, CD->TimLst, rt_step, stepsize * (double) CD->AvgScl,
		    &rt_step, NumEle * 2);
	  // AdptTime(CD, CD->TimLst, stepsize * (double)CD->AvgScl , stepsize * (double)CD->AvgScl, &rt_step, NumEle * 2);

	  if ((int) CD->TimLst % 60 == 0)
	    {
	      tmpconc = 0.0;
	      for (i = 0; i < NumEle; i++)
		{
		  tmpval = MD->Ele[i].temp;
		  CD->Vcele[i].temperature = tmpval;
		  CD->Vcele[i + NumEle].temperature = tmpval;
		  tmpconc += tmpval;
		}
	      tmpconc *= 1 / ((realtype) NumEle);
	      CD->Temperature_avg = tmpconc;
	    }

	  for (i = 0; i < NumEle; i++)
	    for (j = 0; j < CD->NumStc; j++)
	      {
		if (CD->chemtype[j].itype == 4)
		  {
		    CD->Vcele[i].t_conc[j] =
		      (CD->Vcele[i].t_conc[j] * CD->Vcele[i].height_t +
		       CD->Vcele[i +
				 NumEle].t_conc[j] *
		       (CD->Vcele[i].height_v -
			CD->Vcele[i].height_t)) / CD->Vcele[i].height_v;
		    CD->Vcele[i + NumEle].t_conc[j] = CD->Vcele[i].t_conc[j];
		    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
		    CD->Vcele[i + NumEle].p_conc[j] = CD->Vcele[i].t_conc[j];
		    // Averaging mineral concentration to ensure mass conservation !!
		  }
	      }


	  for (i = 0; i < CD->NumOsv; i++)
	    {
	      if (fabs (CD->Vcele[i].height_t - CD->Vcele[i].height_int) >
		  1.0E-6)
		fprintf (stderr, "%d %6.4f\t%6.4f\n", i,
			 CD->Vcele[i].height_t, CD->Vcele[i].height_int);
	      assert (fabs (CD->Vcele[i].height_t - CD->Vcele[i].height_int) <
		      1.0E-6);
	      if (CD->Vcele[i].illness >= 20)
		{
		  for (j = 0; j < CD->NumStc; j++)
		    CD->Vcele[i].t_conc[j] = 1.0E-10;
		  fprintf (stderr,
			   " Cell %d isolated due to proneness to err!\n",
			   CD->Vcele[i].index);
		}
	    }
	  // to make sure intrapolation worked well.
	}
      /* step control ends */

      chem_updater (CD, MD);	/* update essential chem/physical information for pihm */
      CD->TimLst = timelps;

      for (i = NumEle; i < 2 * NumEle; i++)
	{
	  j = i - NumEle;
	  CD->Vcele[i].q = 0.0;
	}

      // Reset fluxes for next averaging stage

      for (k = 0; k < CD->NumDis; k++)
	{
	  CD->Flux[k].velocity = 0.0;
	  CD->Flux[k].flux = 0.0;
	  CD->Flux[k].s_area = 0.0;
	}			// for gw cells, contact area is needed for dispersion;  

      for (k; k < CD->NumFac; k++)
	{
	  CD->Flux[k].flux = 0.0;
	  CD->Flux[k].velocity = 0.0;
	}			// for riv cells, contact area is not needed;
    }
  free (rawtime);

}

void
InitialChemFile (char *filename, int NumBTC, int *BTC_loc)
{

  FILE *Cfile[10];
  char *cfn[10];
  int i;

  cfn[0] = (char *) malloc ((strlen (filename) + 20) * sizeof (char));
  sprintf (cfn[0], "output/%s.conc", filename);
  cfn[1] = (char *) malloc ((strlen (filename) + 20) * sizeof (char));
  sprintf (cfn[1], "output/%s.gwt", filename);
  cfn[2] = (char *) malloc ((strlen (filename) + 20) * sizeof (char));
  sprintf (cfn[2], "output/%s.btcv", filename);
  for (i = 3; i < 3 + NumBTC; i++)
    {
      cfn[i] = (char *) malloc ((strlen (filename) + 20) * sizeof (char));
      sprintf (cfn[i], "output/%s%d.btcv", filename, BTC_loc[i - 3]);
    }

  for (i = 0; i < 3 + NumBTC; i++)
    {
      Cfile[i] = fopen (cfn[i], "w");
      free (cfn[i]);
      fclose (Cfile[i]);
    }


}


void
PrintChem (char *filename, Chem_Data CD, realtype t)
{

  FILE *Cfile[10];
  char *cfn[10];
  int i, j, k, NumEle, NumVol;
  struct tm *timestamp;
  time_t *rawtime;
  double timelps;

  NumEle = CD->NumEle;
  NumVol = CD->NumVol;
  rawtime = (time_t *) malloc (sizeof (time_t));
  *rawtime = (int) (t * 60);
  timestamp = gmtime (rawtime);
  timelps = t - CD->StartTime;

  if ((int) timelps % (720) == 0)
    {
      CD->SPCFlg = 0;

      if (!CD->RecFlg)
	{
	  for (i = 0; i < CD->NumStc; i++)
	    {
	      for (j = NumEle * 2; j < CD->NumOsv; j++)
		if (CD->chemtype[i].itype == 4)
		  CD->Vcele[j].p_conc[i] = CD->Vcele[j].t_conc[i];
		else
		  CD->Vcele[j].p_conc[i] =
		    fabs (CD->Vcele[j].t_conc[i] * 0.1);
	    }
	}
      if (!CD->RecFlg)
	for (i = NumEle * 2 + CD->NumRiv; i < CD->NumOsv; i++)
	  Speciation (CD, i);
      else
	for (i = 0; i < CD->NumOsv; i++)
	  Speciation (CD, i);
    }

  if ((int) timelps % (CD->OutItv * 60) == 0)
    {

      cfn[0] = (char *) malloc ((strlen (filename) + 12) * sizeof (char));
      sprintf (cfn[0], "output/%s.conc", filename);
      Cfile[0] = fopen (cfn[0], "a+");
      cfn[1] = (char *) malloc ((strlen (filename) + 12) * sizeof (char));
      sprintf (cfn[1], "output/%s.btcv", filename);
      Cfile[1] = fopen (cfn[1], "a+");
      cfn[2] = (char *) malloc ((strlen (filename) + 12) * sizeof (char));
      sprintf (cfn[2], "output/%s.gwt", filename);
      Cfile[2] = fopen (cfn[2], "a+");
      for (i = 3; i < 3 + CD->NumBTC; i++)
	{
	  cfn[i] = (char *) malloc ((strlen (filename) + 20) * sizeof (char));
	  sprintf (cfn[i], "output/%s%d.btcv", filename, CD->BTC_loc[i - 3]);
	  Cfile[i] = fopen (cfn[i], "a+");
	  if (Cfile[i] == NULL)
	    fprintf (stderr, " Output BTC not found \n");
	}


      fprintf (Cfile[0], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\n",
	       timestamp->tm_year + 1900, timestamp->tm_mon + 1,
	       timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
      fprintf (Cfile[0], "Cell\t");
      for (i = 0; i < CD->NumStc; i++)
	fprintf (Cfile[0], "%6s\t", CD->chemtype[i].ChemName);
      for (i = 0; i < CD->NumSsc; i++)
	fprintf (Cfile[0], "%6s\t", CD->chemtype[i + CD->NumStc].ChemName);
      fprintf (Cfile[0], "\n");

      for (i = 0; i < NumEle * 2; i++)
	{
	  fprintf (Cfile[0], "%d\t", i + 1);
	  for (j = 0; j < CD->NumStc; j++)
	    fprintf (Cfile[0], "%12.8f\t", log10 (CD->Vcele[i].p_conc[j]));
	  for (j = 0; j < CD->NumSsc; j++)
	    fprintf (Cfile[0], "%12.8f\t", log10 (CD->Vcele[i].s_conc[j]));

	  fprintf (Cfile[0], "\n");
	}

      for (i = NumEle * 2; i < CD->NumOsv; i++)
	{
	  fprintf (Cfile[0], "%d\t", i + 1);
	  for (j = 0; j < CD->NumStc; j++)
	    fprintf (Cfile[0], "%12.8f\t", log10 (CD->Vcele[i].p_conc[j]));
	  for (j = 0; j < CD->NumSsc; j++)
	    fprintf (Cfile[0], "%12.8f\t", log10 (CD->Vcele[i].s_conc[j]));

	  fprintf (Cfile[0], "\n");
	}

      // Output the breakthrough curves "stream chemistry"
      // River is not updated in reaction stage, so a speciation is required before output 

      if (t == CD->StartTime + CD->OutItv * 60)
	{
	  fprintf (Cfile[1], "\t\t\t");
	  for (i = 0; i < CD->NumStc; i++)
	    fprintf (Cfile[1], "%6s\t", CD->chemtype[i].ChemName);
	  for (i = 0; i < CD->NumSsc; i++)
	    fprintf (Cfile[1], "%6s\t",
		     CD->chemtype[i + CD->NumStc].ChemName);
	  fprintf (Cfile[1], "\n");
	}
      /* print the header of file if first time entered */
      fprintf (Cfile[1], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
	       timestamp->tm_year + 1900, timestamp->tm_mon + 1,
	       timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
      for (j = 0; j < CD->NumStc; j++)
	fprintf (Cfile[1], "%12.8f\t", log10 (CD->Vcele[1109].p_conc[j]));
      for (j = 0; j < CD->NumSsc; j++)
	fprintf (Cfile[1], "%12.8f\t", log10 (CD->Vcele[1109].s_conc[j]));
      fprintf (Cfile[1], "\n");

      /*
         if ( t == CD->StartTime + CD->OutItv * 60){
         fprintf(Cfile[2], "\t\t\t");
         for ( i = 0 ; i < CD->NumStc; i++)
         fprintf(Cfile[2], "%6s\t", CD->chemtype[i].ChemName);
         for ( i = 0 ; i < CD->NumSsc; i++)
         fprintf(Cfile[2], "%6s\t", CD->chemtype[i+CD->NumStc].ChemName);
         fprintf(Cfile[2],"\n");
         }
         // print the header of file if first time entered 
         fprintf(Cfile[2],"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
         for (j = 0; j < CD->NumStc; j++)
         fprintf(Cfile[2], "%12.8f\t",log10(CD->Vcele[15].p_conc[j]));
         for (j = 0; j < CD->NumSsc; j++)
         fprintf(Cfile[2], "%12.8f\t",log10(CD->Vcele[15].s_conc[j]));
         fprintf(Cfile[2], "\n");
       */
      fprintf (Cfile[2], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
	       timestamp->tm_year + 1900, timestamp->tm_mon + 1,
	       timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
      for (j = 0; j < NumEle * 2; j++)
	{
	  fprintf (Cfile[2], "%6.4f\t", CD->Vcele[j].height_t);
	}
      fprintf (Cfile[2], "\n");


      for (k = 3; k < 3 + CD->NumBTC; k++)
	{
	  if (t == CD->StartTime + CD->OutItv * 60)
	    {
	      fprintf (Cfile[k], "\t\t\t");
	      for (i = 0; i < CD->NumStc; i++)
		fprintf (Cfile[k], "%6s\t", CD->chemtype[i].ChemName);
	      for (i = 0; i < CD->NumSsc; i++)
		fprintf (Cfile[k], "%6s\t",
			 CD->chemtype[i + CD->NumStc].ChemName);
	      fprintf (Cfile[k], "\n");
	    }
	  /* print the header of file if first time entered */
	  fprintf (Cfile[k], "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t",
		   timestamp->tm_year + 1900, timestamp->tm_mon + 1,
		   timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
	  for (j = 0; j < CD->NumStc; j++)
	    {
	      /*
	         if ((CD->BTC_loc[k - 3] > CD->pumps[0].Pump_Location)
	         && (j == CD->pumps[0].Position_Species))
	         fprintf (Cfile[k], "%12.8f\t",
	         log10 ((CD->Vcele[CD->BTC_loc[k - 3]].t_conc[j] *
	         CD->rivd +
	         CD->pumps[0].Injection_conc * CD->pumps[0].flow_rate) / (CD->rivd +
	         CD->pumps[i].flow_rate)));
	         else
	       */
	      fprintf (Cfile[k], "%12.8f\t",
		       log10 (CD->Vcele[CD->BTC_loc[k - 3]].t_conc[j]));
	    }
	  for (j = 0; j < CD->NumSsc; j++)
	    fprintf (Cfile[k], "%12.8f\t",
		     log10 (CD->Vcele[CD->BTC_loc[k - 3]].s_conc[j]));
	  fprintf (Cfile[k], "\n");
	}

      for (i = 0; i < 3 + CD->NumBTC; i++)
	{
	  free (cfn[i]);
	  fclose (Cfile[i]);
	}

    }
  free (rawtime);

}

void
AdptTime (Chem_Data CD, realtype timelps, double rt_step, double hydro_step,
	  double *start_step, int num_blocks)
{
  double stepsize, org_time, step_rst, end_time, substep;
  double timer1, timer2;
  int i, j, k, m, nr_tmp, nr_max, int_flg, tot_nr, NumEle, NumVol;

  NumEle = CD->NumEle;
  NumVol = CD->NumVol;
  stepsize = *start_step;
  org_time = timelps;
  step_rst = *start_step;
  tot_nr = 0;

  if (rt_step >= hydro_step)
    int_flg = 0;
  else
    int_flg = 1;

  end_time = org_time + hydro_step;

  // simple check to determine whether or not to intrapolate the gw height;

  while (timelps < end_time)
    {

      nr_max = 5;
      nr_tmp = 1;

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

      timer1 = timer ();

      for (i = 0; i < num_blocks; i++)
	for (j = 0; j < CD->NumSpc; j++)
	  {
	    if (CD->chemtype[j].mtype == 2)
	      {
		for (k = 0; k < CD->NumSsc; k++)
		  if ((CD->Totalconc[j][k + CD->NumStc] != 0)
		      && (CD->chemtype[k + CD->NumStc].itype != 1))
		    {
		      //              if ( i == 58) fprintf(stderr, " ##%s Tconc %f\t %s Sconc %f\t Mconc %f\n", CD->chemtype[j].ChemName,  CD->Vcele[i].t_conc[j], CD->chemtype[k+ CD->NumStc].ChemName, CD->Vcele[i].s_conc[k],CD->Vcele[i].t_conc[j] - CD->Totalconc[j][k + CD->NumStc] * CD->Vcele[i].s_conc[k]);
		      CD->Vcele[i].t_conc[j] =
			CD->Vcele[i].t_conc[j] - CD->Totalconc[j][k +
								  CD->NumStc]
			* CD->Vcele[i].s_conc[k] * CD->TimRiv;
		      //              if ( CD->Vcele[i].t_conc[j] < 0.0) fprintf( stderr, " Mixed mobility error! %s %f\n", CD->chemtype[j].ChemName, CD->Vcele[i].t_conc[j]);
		    }
	      }
	  }
      /*
         for ( i = NumEle * 2; i < NumEle * 2 + CD->NumRiv ; i++){
         for ( j = 0 ; j < CD->NumStc; j ++){
         //       CD->Vcele[i].t_conc[j] = CD->Precipitation.t_conc[j] * CD->Condensation ;
         CD->Vcele[i].t_conc[j] = 1.0E-20;
         }
         }
       */
      OS3D (timelps, stepsize, CD);


      for (i = 0; i < num_blocks; i++)
	for (j = 0; j < CD->NumSpc; j++)
	  {
	    if (CD->chemtype[j].mtype == 2)
	      {
		for (k = 0; k < CD->NumSsc; k++)
		  if ((CD->Totalconc[j][k + CD->NumStc] != 0)
		      && (CD->chemtype[k + CD->NumStc].itype != 1))
		    {
		      //              if ( i == 1) fprintf(stderr, " ##%s Tconc %f\t %s Sconc %f\t Mconc %f\n", CD->chemtype[j].ChemName,  CD->Vcele[i].t_conc[j], CD->chemtype[k+ CD->NumStc].ChemName, CD->Vcele[i].s_conc[k],CD->Vcele[i].t_conc[j] + CD->Totalconc[j][k + CD->NumStc] * CD->Vcele[i].s_conc[k]);
		      CD->Vcele[i].t_conc[j] =
			CD->Vcele[i].t_conc[j] + CD->Totalconc[j][k +
								  CD->NumStc]
			* CD->Vcele[i].s_conc[k] * CD->TimRiv;
		      //              if ( CD->Vcele[i].t_conc[j] < 0.0) fprintf( stderr, " Mixed mobility error! %s %f\n", CD->chemtype[j].ChemName, CD->Vcele[i].t_conc[j]);
		    }
	      }
	  }
      /*
         for (i = 0; i < CD->NumPUMP; i++)
         {
         j = CD->pumps[i].Pump_Location;
         //      fprintf(stderr, " T_conc binj %6.4g", CD->Vcele[j].t_conc[CD->pumps[i].Position_Species]);
         CD->Vcele[j].t_conc[CD->pumps[i].Position_Species] +=
         CD->pumps[i].Injection_rate * stepsize / (CD->Vcele[j].porosity *
         CD->Vcele[j].vol);
         fprintf (stderr, " ainj %6.4g, add %6.4g, time %6.4f\n",
         CD->Vcele[j].t_conc[CD->pumps[i].Position_Species],
         CD->pumps[i].Injection_rate, stepsize);
         }
       */

      timer2 = timer ();
      //    fprintf(stderr, " takes %6.4f seconds\n", timer2- timer1);

      if (int_flg)
	{
	  for (i = 0; i < NumVol; i++)
	    {
	      CD->Vcele[i].height_o = CD->Vcele[i].height_t;
	      CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
	    }
	}

      timer1 = timer ();
      m = 0;

      if ((!CD->RecFlg)
	  && ((int) (timelps + stepsize) %
	      (int) (CD->React_delay * stepsize) == 0))
	{
	  for (i = 0; i < num_blocks; i++)
	    {
	      if (CD->Vcele[i].illness < 20)
		if (React
		    (timelps - (CD->React_delay * stepsize),
		     stepsize * CD->React_delay, CD, i, &nr_tmp))
		  {
		    fprintf (stderr, "  ---> React failed at cell %12d.\t",
			     CD->Vcele[i].index);

		    substep = 0.5 * stepsize;
		    k = 2;

		    while (j = React (timelps, substep, CD, i, &nr_tmp))
		      {
			substep = 0.5 * substep;
			k = 2 * k;
			if (substep < 0.5)
			  break;
		      }

		    if (j == 0)
		      {
			tot_nr += nr_tmp;
			fprintf (stderr,
				 " Reaction passed with step equals to %f (1/%d)\n",
				 substep, k);
			for (j = 1; j < k; j++)
			  {
			    React (timelps + j * substep, substep, CD, i,
				   &nr_tmp);
			    tot_nr += nr_tmp;
			  }

		      }
		    //              else
		    //        {
		    //          fprintf(stderr, " Taking smaller step does not work, altering initial guess!\n");
		    /*
		       for ( j = 0; j < CD->NumStc; j ++){
		       if (CD->chemtype[j].itype == 4)
		       CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
		       else
		       CD->Vcele[i].p_conc[j] = fabs(CD->Vcele[i].t_conc[j]*0.1);
		       }
		       //           CD->Vcele[i].illness ++ ;
		       for (j = 0; j < k; j ++){
		       m = MAX(React(timelps + j * substep, substep, CD, i , &nr_tmp),m);
		       }
		     */
		    //              if ( m == 1) ReportError(CD->Vcele[i],CD);
		    //        }
		  }
	      tot_nr += nr_tmp;
	    }
	}
      timer2 = timer ();
      timelps += stepsize;
      if (timelps >= org_time + hydro_step)
	break;
    }
  //  ReportError(CD->Vcele[1075], CD);
  if ((!CD->RecFlg)
      && ((int) (timelps) % (int) (CD->React_delay * stepsize) == 0))
    {
      fprintf (stderr,
	       "  React from %f to %f, Temp:%f/%f, Average NR taken: %f, time elapsed %6.4f seconds\n",
	       timelps - (CD->React_delay * stepsize), timelps,
	       CD->Temperature_avg, CD->Temperature + 273.15,
	       (double) tot_nr / (double) NumEle / 2.0, timer2 - timer1);
    }
}

void
chem_updater (Chem_Data CD, void *DS)
{
}
