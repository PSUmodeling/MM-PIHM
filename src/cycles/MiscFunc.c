#include "Cycles.h"

int IsLeapYear (int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

int doy (int year, int month, int mday, int leap_year_mode)
{
    static const int days[2][13] = {
        {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
        {0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
    };

    int             leap;

    if (leap_year_mode)
        leap = IsLeapYear (year);
    else
        leap = 0;

    return (days[leap][month] + mday);
}

void doy2date (int year, int jday, int *month, int *mday, int leap_year_mode)
{
    static const int eomday[2][13] = {
        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
    };
    int             i = 1;
    int             leap;

    if (leap_year_mode)
        leap = IsLeapYear (year);
    else
        leap = 0;

    while (i <= 12)
    {
        if (eomday[leap][i] >= jday)
        {
            *month = i;
            *mday = jday - eomday[leap][i - 1];
            break;
        }
        else
	    i++;
    }
}

int t2doy (time_t *rawtime)
{
    struct tm      *timestamp;
    int             jday;

    timestamp = gmtime (rawtime);

    jday = doy (timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, 1);

    return (jday);
}

int LT (double x, double y)
{
    int             lt;

    if (x < y && fabs (x - y) < 1e-7)
        lt = 1;
    else
        lt = 0;

    return (lt);
}

int LE (double x, double y)
{
    int             le;

    if (x <= y)
        le = 1;
    else if (fabs (x - y) < 1e-7)
        le = 1;
    else
        le = 0;

    return (le);
}

int EQ (double x, double y)
{
    int             eq;

    if (fabs (x - y) < 1e-7)
        eq = 1;
    else
        eq = 0;

    return (eq);
}

int GE (double x, double y)
{
    int             ge;

    if (x >= y)
        ge = 1;
    else if (fabs (x - y) < 1e-7)
        ge = 1;
    else
        ge = 0;

    return (ge);
}

int GT (double x, double y)
{
    int             gt;

    if (x > y && fabs (x - y) > 1e-7)
        gt = 1;
    else
        gt = 0;

    return (gt);
}
