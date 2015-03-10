/******************************************************************************
 * This subroutine is used to calculate the advection diffusion dispersion
 * It uses a similar OS3D scheme as detailed in Crunchflow user's manual.
 *
 * If you have any questions, concerns, suggestions, please contact me at
 * the following address:
 *     Developer: Chen Bao <baochen.d.s@gmail.com>
 *     Version  : pre-alpha
 *     Date     : June, 2013
 *****************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

/* SUNDIAL Header Files */
#include "../pihm.h"
#include "rt.h"			/* Data Model and Variable Declarations for chemical processes */

#define EPSILON 1.0E-20

void
OS3D (realtype t, realtype stepsize, Chem_Data CD)
{

  // input t and stepsize in the unit of minute

  double **dconc = (double **) malloc (CD->NumOsv * sizeof (double *));
  int i, j, k, jj, node_1, node_2, node_3, node_4, abnormalflg, nr_tmp;
  double
    flux_t, diff_flux, disp_flux, distance, temp_dconc, velocity,
    temp_conc, inv_dist, diff_conc, unit_c, area, r_, beta_, var_height,
    total_prep_mass, timelps, invavg, adpstep;
  double *tmpconc = (double *) malloc (CD->NumSpc * sizeof (double));

  abnormalflg = 0;
  unit_c = 1.0 / 1440;
  total_prep_mass = 0.0;

  // Initalize the allocated array

  for (i = 0; i < CD->NumOsv; i++)
    {
      dconc[i] = (double *) malloc (CD->NumSpc * sizeof (double));
      for (j = 0; j < CD->NumSpc; j++)
	dconc[i][j] = 0.0;
    }


  for (i = 0; i < CD->NumFac; i++)
    {
      CD->Flux[i].q = 0.0;
      node_1 = CD->Flux[i].nodeup - 1;
      node_2 = CD->Flux[i].nodelo - 1;
      node_3 = CD->Flux[i].nodeuu - 1;
      node_4 = CD->Flux[i].nodell - 1;
      flux_t = -CD->Flux[i].flux;
      distance = CD->Flux[i].distance;
      velocity = -CD->Flux[i].velocity;
      area = CD->Flux[i].s_area;
      inv_dist = 1.0 / distance;

      for (j = 0; j < CD->NumSpc; j++)
	{

	  diff_conc = CD->Vcele[node_2].t_conc[j] - CD->Vcele[node_1].t_conc[j];	/* difference in concentration, in M/kg w */
	  diff_flux = 0.0;
	  disp_flux = 0.0;
	  if (CD->Flux[i].BC == 0)
	    {
	      diff_flux =
		CD->chemtype[j].DiffCoe * pow (CD->Vcele[node_1].porosity,
					       CD->Cementation);
	      /* diffusion flux, effective diffusion coefficient  */
	      if (velocity < 0.0)
		disp_flux = velocity * CD->chemtype[j].DispCoe;
	      else
		disp_flux = -velocity * CD->chemtype[j].DispCoe;
	      /* longitudinal dispersion */
	      diff_flux = -diff_flux * inv_dist * diff_conc * area;
	      /* diffusion is in the opposite direction of conc gradient */
	      disp_flux = disp_flux * inv_dist * diff_conc * area;

	    }


	  temp_dconc = 0.0;
	  temp_conc = 0.0;

	  /* uses temp_conc to store the concentration at the surfaces */
	  /* uses temp_dconc to store the concentration changes at the cell */

	  if (CD->TVDFlg == 0)
	    {
	      if (flux_t > 0)
		{
		  temp_conc = CD->Vcele[node_2].t_conc[j];
		}
	      else
		{
		  temp_conc = CD->Vcele[node_1].t_conc[j];
		}
	    }
	  if (CD->TVDFlg == 1)
	    {
	      if (flux_t > 0)
		{
		  if (node_4 > 0)
		    {
		      r_ =
			(CD->Vcele[node_2].t_conc[j] -
			 CD->Vcele[node_4].t_conc[j] +
			 EPSILON) / (CD->Vcele[node_1].t_conc[j] -
				     CD->Vcele[node_2].t_conc[j] + EPSILON);
		      beta_ = MAX (0, MIN (MIN (2, 2 * r_), (2 + r_) / 3));
		      temp_conc =
			CD->Vcele[node_2].t_conc[j] +
			beta_ * (CD->Vcele[node_1].t_conc[j] -
				 CD->Vcele[node_2].t_conc[j]) * 0.5;
		    }
		  else
		    temp_conc = CD->Vcele[node_2].t_conc[j];
		}
	      else
		{
		  if (node_3 > 0)
		    {
		      r_ =
			(CD->Vcele[node_1].t_conc[j] -
			 CD->Vcele[node_3].t_conc[j] +
			 EPSILON) / (CD->Vcele[node_2].t_conc[j] -
				     CD->Vcele[node_1].t_conc[j] + EPSILON);
		      beta_ = MAX (0, MIN (MIN (2, 2 * r_), (2 + r_) / 3));
		      temp_conc =
			CD->Vcele[node_1].t_conc[j] +
			beta_ * (CD->Vcele[node_2].t_conc[j] -
				 CD->Vcele[node_1].t_conc[j]) * 0.5;
		    }
		  else
		    temp_conc = CD->Vcele[node_1].t_conc[j];
		}
	    }

	  else
	    {
	      if (flux_t > 0)
		temp_conc = CD->Vcele[node_2].t_conc[j];
	      else
		temp_conc = CD->Vcele[node_1].t_conc[j];

	    }

	  // Flux[i].BC = 0 normal cell face
	  //            = 1 flux boundary
	  //            = 2 noflow boundary

	  if (CD->Flux[i].BC != 2)
	    temp_dconc += temp_conc * flux_t;
	  if (CD->Flux[i].BC == 0)
	    temp_dconc -= diff_flux + disp_flux;
	  temp_dconc *= unit_c;
	  CD->Flux[i].q = temp_dconc;
	  dconc[node_1][j] += temp_dconc;
	  /*
	     if (( CD->Flux[i].nodeup == 1109) && (j == 2) && ((int)t % 1440 == 0)){
	     fprintf(stderr, "%s fromid %d, toid %d, flux: %f, velo: %f, dist: %f, vol:%f, conc: %f, sconc: %f, nconc: %f, adv_trans: %g, diff_trans: %g, s_area: %f\n",CD->chemtype[j].ChemName, CD->Flux[i].nodeup, CD->Flux[i].nodelo,  CD->Flux[i].flux, CD->Flux[i].velocity, CD->Flux[i].distance, CD->Vcele[CD->Flux[i].nodeup-1].vol, CD->Vcele[CD->Flux[i].nodeup-1].t_conc[j],temp_conc, CD->Vcele[CD->Flux[i].nodelo-1].t_conc[j], temp_conc * flux_t, - diff_flux  - disp_flux, CD->Flux[i].s_area);
	     }
	     if (( CD->Flux[i].nodeup == 58) && (j == 2) && ((int)t % 1440 == 0)){
	     fprintf(stderr, "%s fromid %d, toid %d, flux: %f, velo: %f, dist: %f, vol:%f, conc: %f, sconc: %f, nconc: %f, adv_trans: %g, diff_trans: %g, s_area: %f\n",CD->chemtype[j].ChemName, CD->Flux[i].nodeup, CD->Flux[i].nodelo,  CD->Flux[i].flux, CD->Flux[i].velocity, CD->Flux[i].distance, CD->Vcele[CD->Flux[i].nodeup-1].vol, CD->Vcele[CD->Flux[i].nodeup-1].t_conc[j],temp_conc, CD->Vcele[CD->Flux[i].nodelo-1].t_conc[j], temp_conc * flux_t, - diff_flux  - disp_flux, CD->Flux[i].s_area);
	     }
	     if (( CD->Flux[i].nodeup == 59) && (j == 2) && ((int)t % 1440 == 0)){
	     fprintf(stderr, "%s fromid %d, toid %d, flux: %f, velo: %f, dist: %f, vol:%f, conc: %f, sconc: %f, nconc: %f, adv_trans: %g, diff_trans: %g, s_area: %f\n",CD->chemtype[j].ChemName, CD->Flux[i].nodeup, CD->Flux[i].nodelo,  CD->Flux[i].flux, CD->Flux[i].velocity, CD->Flux[i].distance, CD->Vcele[CD->Flux[i].nodeup-1].vol, CD->Vcele[CD->Flux[i].nodeup-1].t_conc[j],temp_conc, CD->Vcele[CD->Flux[i].nodelo-1].t_conc[j], temp_conc * flux_t, - diff_flux  - disp_flux, CD->Flux[i].s_area);
	     } */
	}
    }


  //Local time step part

  for (i = 0; i < CD->NumOsv; i++)
    {
      if (CD->CptFlg == 1)
	{
	  if ((CD->Vcele[i].rt_step < stepsize)
	      && (CD->Vcele[i].height_t > 1.0E-3)
	      && (CD->Vcele[i].height_o > 1.0E-3))
	    {
	      // use its intrinsic smaller step for small cells/ fast flowing cells ~= slow cells (in term of time marching).
	      if (i < 2 * CD->NumEle + CD->NumRiv - CD->RivOff)
		{
		  timelps = t;
		  invavg = 1.0 / stepsize;
		  adpstep = CD->Vcele[i].rt_step;
		  //      fprintf(stderr, " Local time step at cell %d is %f\n", CD->Vcele[i].index, CD->Vcele[i].rt_step);

		  //      fprintf(stderr, " %d cell is flowing at AdptTime\n",i+1);
		  while (timelps < t + stepsize)
		    {
		      if (adpstep > t + stepsize - timelps)
			{
			  adpstep = t + stepsize - timelps;
			}
		      diff_conc = 0.0;
		      for (j = 0; j < CD->NumSpc; j++)
			{
			  tmpconc[j] =
			    dconc[i][j] * adpstep +
			    CD->Vcele[i].t_conc[j] * (CD->Vcele[i].porosity *
						      0.5 *
						      (CD->Vcele[i].vol_o +
						       CD->Vcele[i].vol));
			  if (CD->PrpFlg)
			    {
			      if (CD->Vcele[i].q > 0.0)
				tmpconc[j] +=
				  CD->Precipitation.t_conc[j] *
				  CD->Vcele[i].q * adpstep * unit_c *
				  CD->Condensation;
			      if (CD->Vcele[i].q < 0.0)
				{
				  tmpconc[j] += 0.0;	// n_0 design
				}
			    }
			  if ((tmpconc[j] < 0.0)
			      && (strcmp (CD->chemtype[j].ChemName, "'H+'")))
			    {
			      fprintf (stderr, "Local time stepping check\n");
			      fprintf (stderr,
				       "negative concentration change at species %s !\n",
				       CD->chemtype[j].ChemName);
			      fprintf (stderr, "Change from fluxes: %8.4g\n",
				       dconc[i][j] * adpstep);
			      fprintf (stderr,
				       "Change from precipitation: %8.4g\n",
				       CD->Precipitation.t_conc[j] *
				       CD->Vcele[i].q * adpstep * unit_c *
				       CD->Condensation);
			      fprintf (stderr, "Original mass: %8.4g\n",
				       CD->Vcele[i].t_conc[j] *
				       (CD->Vcele[i].porosity * 0.5 *
					(CD->Vcele[i].vol_o +
					 CD->Vcele[i].vol)));
			      fprintf (stderr,
				       "Old Conc: %8.4g\t New Conc: %8.4g\n",
				       CD->Vcele[i].t_conc[j], tmpconc[j]);
			      ReportError (CD->Vcele[i], CD);
			      abnormalflg = 1;
			      CD->Vcele[i].illness++;
			    }
			  tmpconc[j] =
			    tmpconc[j] / (CD->Vcele[i].porosity * 0.5 *
					  (CD->Vcele[i].vol +
					   CD->Vcele[i].vol_o));
			  if (CD->Vcele[i].illness < 20)
			    {
			      diff_conc =
				MAX (fabs
				     (CD->Vcele[i].t_conc[j] - tmpconc[j]),
				     diff_conc);
			      CD->Vcele[i].t_conc[j] = tmpconc[j];
			    }
			}
		      if (diff_conc > 1.0E-6)	// which lead to the change in the flux of mass between cells
			{
			  for (j = 0; j < CD->NumSpc; j++)
			    tmpconc[j] = 0.0;

			  for (k = 0; k < CD->NumFac; k++)
			    {
			      if (CD->Flux[k].nodeup == CD->Vcele[i].index)
				{

				  node_1 = CD->Flux[k].nodeup - 1;
				  node_2 = CD->Flux[k].nodelo - 1;
				  node_3 = CD->Flux[k].nodeuu - 1;
				  node_4 = CD->Flux[k].nodell - 1;
				  flux_t = -CD->Flux[k].flux;
				  distance = CD->Flux[k].distance;
				  velocity = -CD->Flux[k].velocity;
				  area = CD->Flux[k].s_area;
				  inv_dist = 1.0 / distance;

				  for (j = 0; j < CD->NumSpc; j++)
				    {
				      diff_conc =
					CD->Vcele[node_2].t_conc[j] -
					CD->Vcele[node_1].t_conc[j];
				      diff_flux = 0.0;
				      disp_flux = 0.0;
				      if (CD->Flux[k].BC == 0)
					{
					  diff_flux =
					    CD->chemtype[j].DiffCoe *
					    pow (CD->Vcele[node_1].porosity,
						 CD->Cementation);

					  if (velocity < 0.0)
					    disp_flux =
					      velocity *
					      CD->chemtype[j].DispCoe;
					  else
					    disp_flux =
					      -velocity *
					      CD->chemtype[j].DispCoe;

					  diff_flux =
					    -diff_flux * inv_dist *
					    diff_conc * area;

					  disp_flux =
					    disp_flux * inv_dist * diff_conc *
					    area;
					}

				      temp_dconc = 0.0;
				      temp_conc = 0.0;

				      if (CD->TVDFlg == 0)
					{
					  if (flux_t > 0)
					    {
					      temp_conc =
						CD->Vcele[node_2].t_conc[j];
					    }
					  else
					    {
					      temp_conc =
						CD->Vcele[node_1].t_conc[j];
					    }
					}
				      if (CD->TVDFlg == 1)
					{
					  if (flux_t > 0)
					    {
					      if (node_4 > 0)
						{
						  r_ =
						    (CD->Vcele[node_2].
						     t_conc[j] -
						     CD->Vcele[node_4].
						     t_conc[j] +
						     EPSILON) /
						    (CD->Vcele[node_1].
						     t_conc[j] -
						     CD->Vcele[node_2].
						     t_conc[j] + EPSILON);
						  beta_ =
						    MAX (0,
							 MIN (MIN (2, 2 * r_),
							      (2 + r_) / 3));
						  temp_conc =
						    CD->Vcele[node_2].
						    t_conc[j] +
						    beta_ *
						    (CD->Vcele[node_1].
						     t_conc[j] -
						     CD->Vcele[node_2].
						     t_conc[j]) * 0.5;
						}
					      else
						temp_conc =
						  CD->Vcele[node_2].t_conc[j];
					    }
					  else
					    {
					      if (node_3 > 0)
						{
						  r_ =
						    (CD->Vcele[node_1].
						     t_conc[j] -
						     CD->Vcele[node_3].
						     t_conc[j] +
						     EPSILON) /
						    (CD->Vcele[node_2].
						     t_conc[j] -
						     CD->Vcele[node_1].
						     t_conc[j] + EPSILON);
						  beta_ =
						    MAX (0,
							 MIN (MIN (2, 2 * r_),
							      (2 + r_) / 3));
						  temp_conc =
						    CD->Vcele[node_1].
						    t_conc[j] +
						    beta_ *
						    (CD->Vcele[node_2].
						     t_conc[j] -
						     CD->Vcele[node_1].
						     t_conc[j]) * 0.5;
						}
					      else
						temp_conc =
						  CD->Vcele[node_1].t_conc[j];
					    }
					}

				      else
					{
					  if (flux_t > 0)
					    temp_conc =
					      CD->Vcele[node_2].t_conc[j];
					  else
					    temp_conc =
					      CD->Vcele[node_1].t_conc[j];

					}

				      if (CD->Flux[k].BC != 2)
					temp_dconc += temp_conc * flux_t;
				      if (CD->Flux[k].BC == 0)
					temp_dconc -= diff_flux + disp_flux;

				      temp_dconc *= unit_c;

				      tmpconc[j] += temp_dconc;
				    }

				}
			    }

			  for (j = 0; j < CD->NumSpc; j++)
			    {
			      dconc[i][j] = tmpconc[j];
			    }
			}

		      timelps += adpstep;
		      if (timelps >= t + stepsize)
			break;
		    }
		}
	    }
	}

      if ((CD->Vcele[i].height_t > 1.0E-3)
	  && (CD->Vcele[i].height_o > 1.0E-3))
	{

	  if (CD->CptFlg)
	    {
	      if (CD->Vcele[i].rt_step < stepsize)
		continue;	// treated in the above section
	      if ((i >= 2 * CD->NumEle) && (i < 2 * CD->NumEle + CD->NumRiv))
		continue;	// treated in the above section
	    }
	  // For blocks with very small content, we just skip it.
	  if (CD->Vcele[i].BC != 2)
	    {
	      for (j = 0; j < CD->NumSpc; j++)
		{
		  tmpconc[j] =
		    dconc[i][j] * stepsize +
		    CD->Vcele[i].t_conc[j] * (CD->Vcele[i].porosity *
					      CD->Vcele[i].vol_o);
		  // need consider the change of concentration at the unsat zone from precipitation.
		  if (CD->PrpFlg)
		    {
		      if (CD->Vcele[i].q > 0.0)
			tmpconc[j] +=
			  CD->Precipitation.t_conc[j] * CD->Vcele[i].q *
			  stepsize * unit_c * CD->Condensation;
		      if (CD->Vcele[i].q < 0.0)
			{
			  tmpconc[j] += 0.0;	// n_0 design
			}
		    }
		  tmpconc[j] =
		    tmpconc[j] / (CD->Vcele[i].porosity * CD->Vcele[i].vol);
		}
	      if (CD->Vcele[i].illness < 20)
		for (j = 0; j < CD->NumSpc; j++)
		  {

		    if ((tmpconc[j] < 0.0)
			&& (strcmp (CD->chemtype[j].ChemName, "'H+'"))
			&& (i < CD->NumEle * 2))
		      {
			fprintf (stderr,
				 "negative concentration change at species %s !\n",
				 CD->chemtype[j].ChemName);
			fprintf (stderr, "Change from fluxes: %8.4g\t",
				 dconc[i][j] * stepsize);
			fprintf (stderr, "Original mass: %8.4g\n",
				 CD->Vcele[i].t_conc[j] *
				 (CD->Vcele[i].porosity *
				  CD->Vcele[i].vol_o));
			fprintf (stderr,
				 "New mass: %8.4g\t New Volume: %8.4g\t Old Conc: %8.4g\t New Conc: %8.4g\t Timestep: %8.4g\n",
				 dconc[i][j] * stepsize +
				 CD->Vcele[i].t_conc[j] *
				 (CD->Vcele[i].porosity * CD->Vcele[i].vol_o),
				 CD->Vcele[i].porosity *
				 CD->Vcele[i].height_t * CD->Vcele[i].area,
				 CD->Vcele[i].t_conc[j], tmpconc[j],
				 CD->Vcele[i].rt_step);
			ReportError (CD->Vcele[i], CD);
			CD->Vcele[i].illness++;
			abnormalflg = 1;
		      }
		    CD->Vcele[i].t_conc[j] = tmpconc[j];
		  }
	    }
	}
    }
  for (i = 0; i < CD->NumOsv; i++)
    free (dconc[i]);

  free (dconc);
  free (tmpconc);
}
