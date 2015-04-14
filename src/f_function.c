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
    realtype        rivArea;
    realtype        rivPerem;
    realtype        eq_Wid;

    switch (rivOrder)
    {
        case 1:
            rivArea = rivDepth * rivCoeff;
            rivPerem = 2.0 * rivDepth + rivCoeff;
            eq_Wid = rivCoeff;
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        case 2:
            rivArea = pow (rivDepth, 2) / rivCoeff;
            rivPerem = 2.0 * rivDepth * pow (1.0 + pow (rivCoeff, 2), 0.5) / rivCoeff;
            eq_Wid = 2.0 * pow (rivDepth + EPS, 1.0 / (rivOrder - 1)) / pow (rivCoeff, 1.0 / (rivOrder - 1));
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        case 3:
            rivArea = 4.0 * pow (rivDepth, 1.5) / (3.0 * pow (rivCoeff, 0.5));
            rivPerem = (pow (rivDepth * (1.0 + 4.0 * rivCoeff * rivDepth) / rivCoeff, 0.5)) + (log (2.0 * pow (rivCoeff * rivDepth, 0.5) + pow (1.0 + 4.0 * rivCoeff * rivDepth, 0.5)) / (2.0 * rivCoeff));
            eq_Wid = 2.0 * pow (rivDepth + EPS, 1.0 / (rivOrder - 1)) / pow (rivCoeff, 1.0 / (rivOrder - 1));
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        case 4:
            rivArea = 3.0 * pow (rivDepth, 4.0 / 3.0) / (2.0 * pow (rivCoeff, 1.0 / 3.0));
            rivPerem = 2.0 * ((pow (rivDepth * (1.0 + 9.0 * pow (rivCoeff, 2.0 / 3.0) * rivDepth), 0.5) / 3.0) + (log (3.0 * pow (rivCoeff, 1.0 / 3.0) * pow (rivDepth, 0.5) + pow (1.0 + 9.0 * pow (rivCoeff, 2.0 / 3.0) * rivDepth, 0.5)) / (9.0 * pow (rivCoeff, 1.0 / 3.0))));
            eq_Wid = 2.0 * pow (rivDepth + EPS, 1.0/ (rivOrder - 1)) / pow (rivCoeff, 1.0 / (rivOrder - 1));
            return returnVal (rivArea, rivPerem, eq_Wid, a_pBool);
        default:
            printf ("\n Relevant Values entered are wrong\n");
            printf (" Depth: %lf\tCoeff: %lf\tOrder: %d\t\n", rivDepth, rivCoeff, rivOrder);
            return 0;
    }
}

void OverlandFlow (realtype **flux, int loci, int locj, realtype avg_y, realtype grad_y, realtype avg_sf, realtype crossA, realtype avg_rough)
{
    flux[loci][locj] = crossA * pow (avg_y, 2.0 / 3.0) * grad_y / (sqrt (fabs (avg_sf)) * avg_rough);
}

void OLFeleToriv (realtype eleYtot, realtype EleZ, realtype cwr, realtype rivZmax, realtype rivYtot, realtype **fluxriv, int loci, int locj, realtype length)
{
    realtype        threshEle;

    if (rivZmax < EleZ)
        threshEle = EleZ;
    else
        threshEle = rivZmax;

    if (rivYtot > eleYtot)
    {
        if (eleYtot > threshEle)
            fluxriv[loci][locj] = cwr * 2.0 * sqrt (2.0 * GRAV ) * length * sqrt (rivYtot - eleYtot) * (rivYtot - threshEle) / 3.0;
        else
        {
            if (threshEle < rivYtot)
                fluxriv[loci][locj] = cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (rivYtot - threshEle) * (rivYtot - threshEle) / 3.0;
            else
                fluxriv[loci][locj] = 0.0;
        }
    }
    else
    {
        if (rivYtot > threshEle)
            fluxriv[loci][locj] = -cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (eleYtot - rivYtot) * (eleYtot - threshEle) / 3.0;
        else
        {
            if (threshEle < eleYtot)
                fluxriv[loci][locj] = -cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (eleYtot - threshEle) * (eleYtot - threshEle) / 3.0;
            else
                fluxriv[loci][locj] = 0.0;
        }
    }
}

realtype avgY (realtype diff, realtype yi, realtype yinabr)
{
    if (diff > 0.0)
    {
        if (yi > 1.0 * EPS / 100.0)
            return 1.0 * yi;
        else
            return 0.0;
    }
    else
    {
        if (yinabr > 1.0 * EPS / 100.0)
            return 1.0 * yinabr;
        else
            return 0.0;
    }
}

realtype effKV (realtype ksatFunc, realtype gradY, realtype macKV, realtype KV, realtype areaF)
{
    if (ksatFunc >= 0.98)
        return (macKV * areaF + KV * (1.0 - areaF) * ksatFunc);
    else
    {
        if (fabs (gradY) * ksatFunc * KV <= 1.0 * KV * ksatFunc)
            return KV * ksatFunc;
        else
        {
            if (fabs (gradY) * ksatFunc * KV < (macKV * areaF + KV * (1.0 - areaF) * ksatFunc))
                return (macKV * areaF * ksatFunc + KV * (1.0 - areaF) * ksatFunc);
            else
                return (macKV * areaF + KV * (1.0 - areaF) * ksatFunc);
        }
    }
}

int macpore_status (realtype ksatFunc, realtype elemSatn, realtype gradY, realtype macKV, realtype KV, realtype areaF)
{
    if (elemSatn >= 0.98)
        return SAT_CTRL;
    else
    {
        if (gradY * ksatFunc * KV <= 1.0 * KV * ksatFunc)
            return MAT_CTRL;
        else
        {
            if (gradY * ksatFunc * KV < (macKV * areaF + KV * (1.0 - areaF) * ksatFunc))
                return APP_CTRL;
            else
                return MAC_CTRL;
        }
    }
}

realtype effKV_new (realtype ksatFunc, realtype elemSatn, int status, realtype macKV, realtype KV, realtype areaF)
{
    if (status == SAT_CTRL)
        return (macKV * areaF + KV * (1.0 - areaF) * ksatFunc);
    else
    {
        if (status == MAT_CTRL)
            return KV * ksatFunc;
        else
        {
            if (status == APP_CTRL)
                return (macKV * areaF * elemSatn + KV * (1.0 - areaF) * ksatFunc);
            else
                return (macKV * areaF + KV * (1.0 - areaF) * ksatFunc);
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
                return (MacKsatH * (tmpY - (aqDepth - MacD)) * areaF + ksatH * (aqDepth - MacD + (tmpY - (aqDepth - MacD)) * (1.0 - areaF))) / tmpY;
        }
        else
            return ksatH;
    }
    else
        return ksatH;
}

