/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : f.c   (PIHM v.2.0)                                            *
 * Function    : Model Kernel: Building ODE system for each physical process   *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * Developer of PIHM v.2.0:  Mukesh Kumar (muk139@psu.edu)		       *
 * Developer of PIHM v.1.0:  Yizhong Qu	(quyizhong@gmail.com)		       *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * NOTE: f.c has gone a massive revamp (essentially rewritten) since PIHM v.1.0*
 *                                                                             *
 *...........MODFICATIONS/ADDITIONS incorporated in f.c (PIHM v.2.0)...........*
 * a) Surface Flow: 							       *
 *	--> Correction of diffusion wave approximation (calculation of dh/ds)  *
 *              i. Calculation of dh/ds performed using planar slope connecting*
 *                 neighboring centroids				       *
 *              ii.Reflection of elements at boundaries and rivers for dh/ds   *
 *		   calculation
 *	--> Correction of kinematic wave approximation (dh/ds calculation based*
 *	    on elevation only instead of head				       *
 *	--> Correction of gradient for cases with steep change in topography   *
 * b) Subsurface Flow:							       *
 *	--> Addition of macropore phenomena				       *
 *	--> Addition of rectangular cell beneath a river element	       *
 *	--> Implementation of two layered subsurface model(sat/unsat) based on *
 *	Richard's eqn							       *
 *	--> Incorporation of Vertical and Horizontal Anisotropy                *
 *	--> Use of geologic data					       *
 * c) River Flow:							       *
 *	--> Correction of kinematic and diff. wave approximation of SV eqn     *
 *	--> Incorporation of flexible river shapes			       *
 *	--> Separate incorporation of leakage and lateral flow		       *
 *	--> Correction of bank overland flow for extreme cases		       *
 *	--> Addition of aquifer cells below river elements		       *
 * c) Surface/Subsurface Coupling:					       *
 *	--> Implementation of First order coupling through (in/ex)filtration   *
 *		based on head continuity 				       *
 * d) Evaporation:							       *
 *	--> Incorporation of ET from ground/subsurface/vegetation	       *
 *	--> Incorporation of landcover properties for calculation of each ET   *
 *	    component							       *
 * e) Computational:							       *
 *	--> Use of temporary state variables in calculation. Note: Never change*
 *		core state variables					       *
 * f) Miscellaneous (other advantages realtive to PIHM1.0): No maximum         *
 *    constraint on gw level. Accordingly, no numerical constraints on subsur- *
 *    face flux terms.Faster Implementation. Led to first large scale model    *
 *    application.
 *-----------------------------------------------------------------------------*
 *									       *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact 				       *
 *	--> Mukesh Kumar (muk139@psu.edu)				       *
 *	--> Prof. Chris Duffy (cxd11@psu.edu)				       *
 * This code is free for research purpose only.		       		       *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * REFERENCES:                                                                 *
 * PIHM2.0:                                                                    *
 *      a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *      Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *      b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *      a Mesoscale Watershed", Advances in Water Resources (submitted)        *
 * PIHM1.0:                                                                    *
 *      a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *      simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *      b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *      for multiprocess watershed simulation". Water Resources Research       *
 *******************************************************************************/

#include "pihm.h"               // 09.23

/* 09.23
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
//#include "sundialstypes.h"
#include "sundials_types.h"   // 09.16
#include "pihm.h"
*/

#define multF	2
#define MINpsi	-70
#define EPS 0.05
#define THRESH 0.0
#define UNIT_C 1440             /* Note 60*24 for calculation of yDot in m/min units while forcing is in m/day. */
//#define GRAV 9.8*60*60    /* Note the dependence on physical units */  // 09.23 redefined "GRAV" in pihm_const.h


realtype returnVal(realtype rArea, realtype rPerem, realtype eqWid,
    realtype ap_Bool)
{
    if (ap_Bool == 1)
    {
        return rArea;
    }
    else if (ap_Bool == 2)
    {
        return rPerem;
    }
    else
    {
        return eqWid;
    }
}

realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff,
    realtype a_pBool)
{
    realtype        rivArea, rivPerem, eq_Wid;
    switch (rivOrder)
    {
        case 1:
            rivArea = rivDepth * rivCoeff;
            rivPerem = 2.0 * rivDepth + rivCoeff;
            eq_Wid = rivCoeff;
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        case 2:
            rivArea = pow(rivDepth, 2) / rivCoeff;
            rivPerem =
                2.0 * rivDepth * pow(1 + pow(rivCoeff, 2), 0.5) / rivCoeff;
            eq_Wid =
                2.0 * pow(rivDepth + EPS, 1 / (rivOrder - 1)) / pow(rivCoeff,
                1 / (rivOrder - 1));
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        case 3:
            rivArea = 4 * pow(rivDepth, 1.5) / (3 * pow(rivCoeff, 0.5));
            rivPerem =
                (pow(rivDepth * (1 + 4 * rivCoeff * rivDepth) / rivCoeff,
                    0.5)) + (log(2 * pow(rivCoeff * rivDepth,
                        0.5) + pow(1 + 4 * rivCoeff * rivDepth,
                        0.5)) / (2 * rivCoeff));
            eq_Wid =
                2.0 * pow(rivDepth + EPS, 1 / (rivOrder - 1)) / pow(rivCoeff,
                1 / (rivOrder - 1));
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        case 4:
            rivArea =
                3 * pow(rivDepth, 4.0 / 3.0) / (2 * pow(rivCoeff, 1.0 / 3.0));
            rivPerem =
                2 * ((pow(rivDepth * (1 + 9 * pow(rivCoeff,
                                2.0 / 3.0) * rivDepth),
                        0.5) / 3) + (log(3 * pow(rivCoeff,
                            1.0 / 3.0) * pow(rivDepth,
                            0.5) + pow(1 + 9 * pow(rivCoeff,
                                2.0 / 3.0) * rivDepth,
                            0.5)) / (9 * pow(rivCoeff, 1.0 / 3.0))));
            eq_Wid =
                2.0 * pow(rivDepth + EPS, 1 / (rivOrder - 1)) / pow(rivCoeff,
                1 / (rivOrder - 1));
            return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
        default:
            printf("\n Relevant Values entered are wrong");
            //  printf("\n Depth: %lf\tCoeff: %lf\tOrder: %d\t");
            return 0;
    }
}
