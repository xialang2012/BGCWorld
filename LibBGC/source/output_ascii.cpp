/* ASCII output function for Biome-BGC 		*/
/* Written by: Matt Jolly 			*/
/* Dev. Date : 15 Feb 2005			*/
/* Rev. Date :					*/
/* Purpose: take a set of binary inputs and write them as ascii outputs */

/**
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
**/

#include "bgc.h"

int output_ascii(float arr[],int nvars,  FILE *ptr)
{
	int i = 0;
	
	for(i = 0;i < nvars;i++){ fprintf(ptr,"%13.8f\t",arr[i]);}
	fprintf(ptr,"\n");
	
	return(EXIT_SUCCESS);

}
	

int output_soil(soilvar_struct arr[], FILE *ptr)
{
	int i = 0;
	fprintf(ptr,"depth\ttemp \twater\t ice \t Twater\t head \n");
	for(i = 0;i < N;i++)
	{ 
		fprintf(ptr,"%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",i,arr[i].s_t,arr[i].s_w,arr[i].s_i,arr[i].s_tw,arr[i].s_h);
	}
	fprintf(ptr,"\n");
	
	return(EXIT_SUCCESS);

}