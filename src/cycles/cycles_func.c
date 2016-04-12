#include "pihm.h"

int DOY (int t)
{
    struct tm      *timestamp;
    int             year, month, mday;
    time_t          rawtime;
    static const int days[2][13] = {
        {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
        {0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
    };
    int             leap;

    rawtime = (time_t)t;
    timestamp = gmtime (&rawtime);

    year = timestamp->tm_year + 1900;
    month = timestamp->tm_mon + 1;
    mday = timestamp->tm_mday;

    leap = IsLeapYear (year);

    return (days[leap][month] + mday);
}

int IsLeapYear (int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}
