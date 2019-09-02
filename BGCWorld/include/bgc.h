#ifndef BGC_H
#define BGC_H

/*
bgc.h
header file to hold includes needed for all of the bgc library.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
Andrew A Neuschwander, andrew@ntsg.umt.edu

*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <errno.h>

#include "bgc_struct.h"
#include "bgc_func.h"
#include "bgc_constants.h"
#include "ini.h"
#include "bgc_epclist.h"
#include "bgc_io.h"
#include "misc_func.h"
#include "soil_func.h"   /*add 2014-07-30*/

#ifdef __cplusplus
extern "C"
{
#endif

/* non include stuff here */
	//extern soilpar_struct vspar;      
	extern soilvar_struct svar[];
#ifdef __cplusplus
}
#endif

#endif
