
/******************************************************************************
 * File		: print.h                                                         *
 * Function	: Declaration and Definition of print structures and function     *
 * Developer of PIHM 2.3    :	Yuning Shi	(yshi@psu.edu)                    *
 * Developer of PIHM 2.2    :	Xuan Yu	    (xxy113@psu.edu)                  *
 * Developer of PIHM 2.0    :	Mukesh Kumar	(muk139@psu.edu)              *
 * Developer of PIHM 1.0    :	Yizhong Qu	(quyizhong@gmail.com)             *
 * Version                  :   September, 2014 (PIHM V2.3)                   *
 *----------------------------------------------------------------------------*
 * This code is free for research purpose only.				                  *
 * Please provide relevant references if you use this code in your research   *
 * 	work								                                      *
 *----------------------------------------------------------------------------*
 *****************************************************************************/
#ifndef PRINT_HEADER
#define PRINT_HEADER

typedef struct Print_Ctrl_structure /* YS */
{
    char            name[100];
    int             Interval;
    int             NumVar;
    double        **PrintVar;
    double         *buffer;
} Print_Ctrl;

void            PrintData (Print_Ctrl, double, int);
#endif
