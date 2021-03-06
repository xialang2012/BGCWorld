/* 
get_co2.c
retrieve the appropriate CO2 concentration for the current simulation year

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

double get_co2(co2control_struct * co2,int simyr)
{
	int i;
	for(i = 0;i < co2->co2vals;i++)
	{
		//printf("co2year = %d, simyr=%d \n",co2->co2year_array[i],simyr);
		if(co2->co2year_array[i] == simyr)
		{
			//printf("in get_co2(), simyr=%d ,co2ppm=%lf \n",simyr,co2->co2ppm_array[i]);
			return (co2->co2ppm_array[i]);
		}
	}
	return(-999.9);
}
