#ifndef PIHM_ERRORS_HEADER
#define PIHM_ERRORS_HEADER

#define ERR_WRONG_FORMAT        "input file format error in %s at Line %d."
#define ERR_FORC_NOT_FOUND      "%s forcing for current time step not found.\nPlease check your %s forcing file.\n"
#define ERR_LC_NOT_FOUND        "land cover type %d is not described in the vegetation file.\nPlease check your .att file."
#define ERR_SHAPE_NOT_FOUND     "river shape %d is not defined.\nPlease check your .riv file."
#define ERR_SOIL_NOT_FOUND      "soil type %d is not described in the soil file.\nPlease check your .att file."
#define ERR_POROSITY_OUT_BOUND  "porosity value for soil type %d out of bounds."

#endif
