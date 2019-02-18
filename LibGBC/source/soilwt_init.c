#include "bgc.h"

/* The variable of soil including temperature,water content,ice content,water head */
int soilwt_init(soilvar_struct svar1[])
{  
	//FILE *fpp;
	
	
	/* typedef	  soilvar_struct
    double s_t;          (*C)soil temperature
	double s_w;            soil liquid water content
	double s_i;            soil ice content
	double s_h;            (m)soil water head
	double s_w_max;        the maxmum of the soil liquid water content
	double s_kt;          (Wm-3K-1)the thermal conductivity of the soil
	double s_ct;          (Jm-3K-1)the capacity of the soil 
	double s_kw;          (ms-1)the hydralic conductivity of the soil 
	*/
	//extern soilpar_struct vspar;      
	//extern soilvar_struct svar[10];    /* N is the number of layers of soil sign constant :define N 10 */
	int i;
	int ok=1;
  // fpp=fopen("d:\\s_init.txt","r");
	
   
	for(i=0;i<N;i++)
	{
		svar1[i].s_t=2.01;
		svar1[i].s_w=0.2188;
		svar1[i].s_tw=0.2188;
		svar1[i].s_i=0;
		svar1[i].s_h=-2.3;//Waterhead1(svar1[i].s_w,svar1[i].s_i);
		svar1[i].s_tt=0.0;
	}
   
   //Uusing the data of run 2years result as the initization
 /* for(i=0;i<N;i++)
   {
	   if (ok && fscanf(fpp,"%*s%lf%lf%lf%lf",&svar1[i].s_h,&svar1[i].s_w,&svar1[i].s_i,&svar1[i].s_t)==EOF)
		{
			printf("Error reading s_init file, svar_init()\n");
			ok=0;
		}
	
   }
  fclose(fpp);*/
	
	return (!ok);
 
}