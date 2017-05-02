#include "pihm.h"

pihm_t_struct PIHMTime (int t)
{
    pihm_t_struct   pihm_time;
    struct tm      *timestamp;
    time_t          rawtime;

    rawtime = (time_t)t;
    timestamp = gmtime (&rawtime);

    pihm_time.t = t;
    pihm_time.year = timestamp->tm_year + 1900;
    pihm_time.month = timestamp->tm_mon + 1;
    pihm_time.day = timestamp->tm_mday;
    pihm_time.hour = timestamp->tm_hour;
    pihm_time.minute = timestamp->tm_min;
    pihm_time.second = timestamp->tm_sec;
    strftime (pihm_time.str, 17, "%Y-%m-%d %H:%M", timestamp);

    return (pihm_time);
}

