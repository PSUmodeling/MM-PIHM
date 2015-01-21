#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pihm.h"

realtype returnVal (realtype rArea, realtype rPerem, realtype eqWid, realtype ap_Bool)
{
    if (ap_Bool == 1)
        return rArea;
    else if (ap_Bool == 2)
        return rPerem;
    else
        return eqWid;
}

realtype CS_AreaOrPerem (int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool)
{
    realtype        rivArea, rivPerem, eq_Wid;
    switch (rivOrder)
    {
        case 1:
            rivArea = rivDepth * rivCoeff;
            rivPerem = 2.0 * rivDepth + rivCoeff;
            eq_Wid = rivCoeff;
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        case 2:
            rivArea = pow (rivDepth, 2) / rivCoeff;
            rivPerem = 2.0 * rivDepth * pow (1 + pow (rivCoeff, 2), 0.5) / rivCoeff;
            eq_Wid = 2.0 * pow (rivDepth + EPS, 1 / (rivOrder - 1)) / pow (rivCoeff, 1 / (rivOrder - 1));
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        case 3:
            rivArea = 4 * pow (rivDepth, 1.5) / (3 * pow (rivCoeff, 0.5));
            rivPerem = (pow (rivDepth * (1 + 4 * rivCoeff * rivDepth) / rivCoeff, 0.5)) + (log (2 * pow (rivCoeff * rivDepth, 0.5) + pow (1 + 4 * rivCoeff * rivDepth, 0.5)) / (2 * rivCoeff));
            eq_Wid = 2.0 * pow (rivDepth + EPS, 1 / (rivOrder - 1)) / pow (rivCoeff, 1 / (rivOrder - 1));
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        case 4:
            rivArea = 3 * pow (rivDepth, 4.0 / 3.0) / (2 * pow (rivCoeff, 1.0 / 3.0));
            rivPerem = 2 * ((pow (rivDepth * (1 + 9 * pow (rivCoeff, 2.0 / 3.0) * rivDepth), 0.5) / 3) + (log (3 * pow (rivCoeff, 1.0 / 3.0) * pow (rivDepth, 0.5) + pow (1 + 9 * pow (rivCoeff, 2.0 / 3.0) * rivDepth, 0.5)) / (9 * pow (rivCoeff, 1.0 / 3.0))));
            eq_Wid = 2.0 * pow (rivDepth + EPS, 1 / (rivOrder - 1)) / pow (rivCoeff, 1 / (rivOrder - 1));
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        default:
            printf ("\n Relevant Values entered are wrong");
            printf ("\n Depth: %lf\tCoeff: %lf\tOrder: %d\t");
            return 0;
    }
}

void OverlandFlow (realtype ** flux, int loci, int locj, realtype avg_y, realtype grad_y, realtype avg_sf, realtype crossA, realtype avg_rough)
{
    flux[loci][locj] = crossA * pow (avg_y, 2.0 / 3.0) * grad_y / (sqrt (fabs (avg_sf)) * avg_rough);
    //  flux[loci][locj] = (grad_y>0?1:-1)*crossA*pow(avg_y, 2.0/3.0)*sqrt(fabs(grad_y))/(avg_rough);
}

void OLFeleToriv (realtype eleYtot, realtype EleZ, realtype cwr, realtype rivZmax, realtype rivYtot, realtype ** fluxriv, int loci, int locj, realtype length)
{
    realtype        threshEle;
    if (rivZmax < EleZ)
        threshEle = EleZ;
    else
        threshEle = rivZmax;
    if (rivYtot > eleYtot)
    {
        if (eleYtot > threshEle)
            fluxriv[loci][locj] = cwr * 2.0 * sqrt (2 * GRAV ) * length * sqrt (rivYtot - eleYtot) * (rivYtot - threshEle) / 3.0;
        else
        {
            if (threshEle < rivYtot)
                fluxriv[loci][locj] = cwr * 2.0 * sqrt (2 * GRAV) * length * sqrt (rivYtot - threshEle) * (rivYtot - threshEle) / 3.0;
            else
                fluxriv[loci][locj] = 0.0;
        }
    }
    else
    {
        if (rivYtot > threshEle)
            fluxriv[loci][locj] = -cwr * 2.0 * sqrt (2 * GRAV) * length * sqrt (eleYtot - rivYtot) * (eleYtot - threshEle) / 3.0;
        else
        {
            if (threshEle < eleYtot)
                fluxriv[loci][locj] = -cwr * 2.0 * sqrt (2 * GRAV) * length * sqrt (eleYtot - threshEle) * (eleYtot - threshEle) / 3.0;
            else
                fluxriv[loci][locj] = 0.0;
        }
    }
}

/*
 * realtype avgY(realtype zi,realtype zinabr,realtype yi,realtype yinabr)
 * {
 * if(zinabr>zi)
 * {
 * if(zinabr>zi+yi)
 * {
 * return yinabr/2;
 * }
 * else
 * {
 * return (yi+zi-zinabr+yinabr)/2;
 * }
 * }
 * else
 * {
 * if(zi>zinabr+yinabr)
 * {
 * return yi/2;
 * }
 * else
 * {
 * return (yi+yinabr+zinabr-zi)/2;
 * }
 * }    
 * }
 */

/*
 * Note: above formulation doesnt' satisfies upwind/downwind scheme.
 */
realtype avgY (realtype diff, realtype yi, realtype yinabr)
{
    if (diff > 0)
    {
        if (yi > 1 * EPS / 100)
        {
            //          return 0.5*(yi+yinabr);
            //          return ((yinabr>yi)?0:1.0*yi);  /* Note the if-else TRUE case can be possible only for Kinematic case */
            return 1.0 * yi;
        }
        else
            return 0;
    }
    else
    {
        if (yinabr > 1 * EPS / 100)
        {
            //          return 0.5*(yi+yinabr);
            //          return ((yi>yinabr)?0:1.0*yinabr);  /* Note the if-else TRUE case can be possible only for Kinematic case */
            return 1.0 * yinabr;
        }
        else
            return 0;
    }
}

realtype effKV (realtype ksatFunc, realtype gradY, realtype macKV, realtype KV, realtype areaF)
{
    if (ksatFunc >= 0.98)
        return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
    else
    {
        if (fabs (gradY) * ksatFunc * KV <= 1 * KV * ksatFunc)
            return KV * ksatFunc;
        else
        {
            if (fabs (gradY) * ksatFunc * KV < (macKV * areaF + KV * (1 - areaF) * ksatFunc))
                return (macKV * areaF * ksatFunc + KV * (1 - areaF) * ksatFunc);
            else
                return (macKV * areaF + KV * (1 - areaF) * ksatFunc);
        }
    }
}

realtype effKH (int mp, realtype tmpY, realtype aqDepth, realtype MacD, realtype MacKsatH, realtype areaF, realtype ksatH)
{
    if (mp == 1)
    {
        if (tmpY > aqDepth - MacD)
        {
            if (tmpY > aqDepth)
                return (MacKsatH * MacD * areaF + ksatH * (aqDepth - MacD * areaF)) / aqDepth;
            else
                return (MacKsatH * (tmpY - (aqDepth - MacD)) * areaF + ksatH * (aqDepth - MacD + (tmpY - (aqDepth - MacD)) * (1 - areaF))) / tmpY;
        }
        else
            return ksatH;
    }
    else
        return ksatH;
}

