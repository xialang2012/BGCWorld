#ifndef POINTBGC_H
#define POINTBGC_H

/*
pointbgc.h
header file to hold includes needed for pointbgc.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
Andrew A Neuschwander, andrew@ntsg.umt.edu

*/

#include <time.h>

#if __GNUC__
#include "getopt.h"
#if __x86_64__ || __ppc64__
#define HOST "GNU 64bit"
#else
#define HOST "GNU 32bit"
#endif
#endif

#if _WIN32 || _WIN64
#include "getopt.h"
#define VERS "4.2"
#define USER "unknow"
#if _WIN64
#define HOST "Windows 64bit"
#else
#define HOST "Windows 32bit"
#endif
#endif

#include "bgc.h"
#include "pointbgc_struct.h"
#include "pointbgc_func.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* non include stuff here */

#ifdef __cplusplus
}
#endif

#endif
