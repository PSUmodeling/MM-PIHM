#ifndef PIHM_ERRORS_HEADER
#define PIHM_ERRORS_HEADER

#define ERR_WRONG_FORMAT        "input file format error in %s at Line %d."
#define ERR_FORC_NOT_FOUND      "%s forcing for current time step not found.\nPlease check your %s forcing file.\n"
#define ERR_LC_NOT_FOUND        "land cover type %d is not described in the vegetation file.\nPlease check your .att file."
#define ERR_SHAPE_NOT_FOUND     "river shape %d is not defined.\nPlease check your .riv file."
#define ERR_SOIL_NOT_FOUND      "soil type %d is not described in the soil file.\nPlease check your .att file."
#define ERR_POROSITY_OUT_BOUND  "porosity value for soil type %d out of bounds."

#define ERR_INVALID_CMD_LINE_OPTION     "%s: %s."
#define ERR_INVALID_DOY                 "%s %d is invalid in %s at Line %d."
#define ERR_EXPECTED_HEADER             "Expect header \"%s\"."
#define ERR_FILE_CLOSE_ERROR            "Error closing %s."
#define ERR_FILL_MISSING_VALUE          "%s missing in %s at Line %d. Default value %s will be used."
#define ERR_WRONG_HEADER                "File header error in %s at Line %d."
#define ERR_WRONG_KEYWORD               "Expected keyword \"%s\", detected keyword \"%s\"in %s at Line %d."
#define ERR_INCONSIST_LEAP_YEAR         "Some leap years are missing DOY 366 input in %s."
#define ERR_LOWER_THAN_LOW_BOUND        "%s is lower than %s in %s at Line %d."
#define ERR_PTF_LOWER_THAN_LOW_BOUND    "%s is lower than %s from pedotransfer functions in Layer %d."
#define ERR_N_BALANCE                   "Soil nitrogen balance error.\nInitial N = %lf, final N = %lf."
#define ERR_NO_PROJECT_NAME             "Please specify the name of project!\n" \
                                        "\nUsage: " \
                                        "./Cycles [-v or -b] [-c] [-d] [-l] [-m] [-s] [-e <extension>] [-p <path>] <project>\n" \
                                        "    -b Brief mode\n" \
                                        "    -v Verbose mode\n" \
                                        "    -c Calibration mode\n" \
                                        "    -d Debug mode\n" \
                                        "    -l Baseline simulation\n" \
                                        "    -m Multiple simulation mode\n" \
                                        "    -s Spin-up mode\n" \
                                        "    -e Output file extension\n" \
                                        "    -p Path to input/output directories\n\n" \
                                        "To check Cycles version and platform:\n\n" \
                                        "    ./Cycles --version or ./Cycles -V\n"
#define ERR_NON_OPTION                  "%s option %s in %s at Line %d is not supported."
#define ERR_NOT_FOUND                   "Cannot find %s of %s in %s."
#define ERR_OLD_TILLAGE                 "This version of Cycles does not support GRAIN_HARVEST and\n" \
                                        "FORAGE_HARVEST parameters in KILL_CROP tillage operations. To schedule a\n" \
                                        "grain/forage harvest operation, please use create a tillage operation with tool\n" \
                                        "name of GRAIN_HARVEST or FORAGE_HARVEST.\n\n" \
                                        "Note that the scheduled grain/forage harvest operations will not kill crops,\n" \
                                        "regardless of the crop file settings. To kill the crop after scheduled harvest,\n" \
                                        "please use a KILL_CROP tillage operation.\n\n"
#define ERR_NO_ORGANIC_C                "Fertilizer %s in %s at Line %d has %.1f%% organic N but no organic C.\n" \
                                        "Organic C fraction has been changed to %.1f%% using a default CN ratio of 15."
#define ERR_OUT_OF_RANGE                "%s is out of range (%s -- %s) in %s at Line %d."
#define ERR_INVALID_ROTATION_YEAR       "The %s operation scheduled on Day %d Year %d\n" \
                                        "in %s at Line %d will not be performed\n" \
                                        "because operation year is larger than years in rotation."
#define ERR_SHORT_INPUT_FILE            "Input file %s is too short."
#define ERR_LARGE_CN_RATIO              "CN ratio too large in %s."

#endif
