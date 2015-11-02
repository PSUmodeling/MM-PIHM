#include "Cycles.h"

void ReadSimControl (char *project, SimControlStruct *SimControl)
{
    /*
     * Read simulation control file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * simctrl_file	    FILE*	File pointer of simulation control
     *					  file
     * filename		    char*	Simulation control file name
     */
    FILE           *simctrl_file;
    char            filename[MAXSTRING];
    char            cmdstr[MAXSTRING];

    /* Open simulation control file */
    sprintf (filename, "input/%s.ctrl", project);
    simctrl_file = fopen (filename, "r");
    printf ("%-30s %s.\n", "Read simulation control file:", filename);

    if (simctrl_file == NULL)
    {
        printf ("ERROR: Cannot find the simulation control file %s!\n", filename);
        exit (1);
    }

    /* Read simulation control file */
    FindLine (simctrl_file, "BOF");

    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "SIMULATION_START_YEAR", &SimControl->simStartYear);

    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "SIMULATION_END_YEAR", &SimControl->simEndYear);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "ROTATION_SIZE", &SimControl->yearsInRotation);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "ADJUSTED_YIELDS", &SimControl->adjustedYields);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "HOURLY_INFILTRATION", &SimControl->hourlyInfiltration);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "AUTOMATIC_NITROGEN", &SimControl->automaticNitrogen);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "AUTOMATIC_PHOSPHORUS", &SimControl->automaticPhosphorus);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "AUTOMATIC_SULFUR", &SimControl->automaticSulfur);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_WEATHER_OUT", &SimControl->weatherDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_CROP_OUT", &SimControl->cropDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_REDISUE", &SimControl->residueDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_WATER_OUT", &SimControl->waterDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_NITROGEN_OUT", &SimControl->nitrogenDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_SOIL_CARBON", &SimControl->soilCarbonDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "DAILY_SOIL_OUT", &SimControl->soilDailyOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "ANNUAL_SOIL_OUT", &SimControl->annualSoilOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "ANNUAL_PROFILE_OUT", &SimControl->profileOutput);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "SEASON_OUT", &SimControl->seasonOutput);

    NextLine (simctrl_file, cmdstr);
    ReadKeywordStr (cmdstr, "CROP_FILE", SimControl->crop_filename);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordStr (cmdstr, "OPERATION_FILE", SimControl->operation_filename);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordStr (cmdstr, "SOIL_FILE", SimControl->soil_filename);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordStr (cmdstr, "WEATHER_FILE", SimControl->weather_filename);

    fclose (simctrl_file);

    SimControl->totalYears = SimControl->simEndYear - SimControl->simStartYear + 1;
}
