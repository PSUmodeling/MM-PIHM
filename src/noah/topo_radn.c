#include "pihm.h"

double TopoRadn(double sdir, double sdif, double zenith, double azimuth180, const topo_struct *topo)
{
    double          incidence;              // Sun incidence angle (degree)
    double          tcf;                    // terrain configuration factor (-)
    double          soldown;

    azimuth180 = Mod(360.0 + azimuth180, 360.0);

    // If the Sun is blocked, set direct solar radiation to 0.0
    sdir = (zenith > topo->h_phi[(int)floor(azimuth180 / 10.0)]) ? 0.0 : sdir;

    // Calculate Sun incidence angle
    incidence = acos(cos(zenith * PI / 180.0) * cos(topo->slope * PI / 180.0) +
        sin(zenith * PI / 180.0) * sin(topo->slope * PI / 180.0) * cos((azimuth180 - topo->aspect) * PI / 180.0));
    incidence *= 180.0 / PI;
    incidence = MIN(incidence, 90.0);

    // Calculate terrain configuration factor
    // Dozier and Frew 1990, IEEE Transactions on Geoscience and Remote Sensing, 28(5), 963--969
    tcf = (1.0 + cos(topo->slope * PI / 180.0)) / 2.0 - topo->svf;
    tcf = MAX(tcf, 0.0);

    soldown = sdir * cos(incidence * PI / 180.0) + topo->svf * sdif +
        0.2 * tcf * (sdir * cos(zenith * PI / 180.0) + sdif);
    soldown = MAX(soldown, 0.0);

    return soldown;
}

void SunPos(int t, const siteinfo_struct *siteinfo, spa_data *spa)
{
    int             spa_result;
    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime(t);

    spa->year          = pihm_time.year;
    spa->month         = pihm_time.month;
    spa->day           = pihm_time.day;
    spa->hour          = pihm_time.hour;
    spa->minute        = pihm_time.minute;
    spa->second        = 0;
    spa->timezone      = 0;

    spa->delta_t       = 67;
    spa->delta_ut1     = 0;
    spa->atmos_refract = 0.5667;

    spa->longitude = siteinfo->longitude;
    spa->latitude  = siteinfo->latitude;
    spa->elevation = siteinfo->zmax;

    // Calculate surface pressure based on FAO 1998 method (Narasimhan 2002)
    spa->pressure = 1013.25 * pow((293.0 - 0.0065 * spa->elevation) / 293.0, 5.26);
    spa->temperature = siteinfo->tavg;

    spa->function = SPA_ZA_RTS;
    spa_result = spa_calculate(spa);

    if (spa_result != 0)
    {
        pihm_printf(VL_ERROR, "Error with spa error code: %d.\n", spa_result);
        pihm_exit(EXIT_FAILURE);
    }

}
