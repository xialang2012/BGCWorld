/* 
sitec_init.c
Initialize the site physical constants for bgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "pointbgc.h"

int sitec_init(file init, siteconst_struct* sitec)
{
	/* reads the site physical constants from *.init */ 

	int ok=1;
	char key[] = "SITE";
	char keyword[80];
	/*double sand,silt,clay;*/ /* percent sand, silt, and clay */

	/* first scan keyword to ensure proper *.init format */ 
	if (ok && scan_value(init, keyword, 's'))
	{
		bgc_printf(BV_ERROR, "Error reading keyword, sitec_init()\n");
		ok=0;
	}
	if (ok && strcmp(keyword,key))
	{
		bgc_printf(BV_ERROR, "Expecting keyword --> %s in %s\n",key,init.name);
		ok=0;
	}

	/* begin reading constants from *.init */
	if (ok && scan_value(init, &sitec->soil_depth, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading soil depth, sitec_init()\n");
		ok=0;
	}
	//printf("读取soil depth\n");
	/*if (ok && scan_value(init, &sand, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading percent sand, sitec_init()\n");
		ok=0;
	}
	if (ok && scan_value(init, &silt, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading percent clay, sitec_init()\n");
		ok=0;
	}
	if (ok && scan_value(init, &clay, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading percent clay, sitec_init()\n");
		ok=0;
	}
	
	if (ok && scan_value(init, &sitec->elev, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading elevation, sitec_init()\n");
		ok=0;
	}
	if (ok && scan_value(init, &sitec->lat, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading site latitude, sitec_init()\n");
		ok=0;
	}
	if (ok && scan_value(init, &sitec->sw_alb, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading shortwave albedo, sitec_init()\n");
		ok=0;
	}*/
	if (ok && scan_value(init, &sitec->ndep, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading N deposition, sitec_init()\n");
		ok=0;
	}
	//printf("读取soil n deposition\n");
	if (ok && scan_value(init, &sitec->nfix, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading N fixation, sitec_init()\n");
		ok=0;
	}
	//printf("读取soil n fixation\n");
	/* calculate the soil pressure-volume coefficients from texture data */
	/* Uses the multivariate regressions from Cosby et al., 1984 */
	/* first check that the percentages add to 100.0 */
	//if (ok && fabs(sand+silt+clay-100) > FLT_COND_TOL)
	//{
	//	bgc_printf(BV_ERROR, "Error: %%sand + %%silt + %%clay  MUST EQUAL 100%%\n");
	//	bgc_printf(BV_ERROR, "Check values in initialization file.\n");
	//	ok=0;
	//}
	//if (ok)
	//{
	//	sitec->soil_b = -(3.10 + 0.157*clay - 0.003*sand);
	//	sitec->vwc_sat = (50.5 - 0.142*sand - 0.037*clay)/100.0;
	//	sitec->psi_sat = -(exp((1.54 - 0.0095*sand + 0.0063*silt)*log(10.0))*9.8e-5);
	//	sitec->vwc_fc = sitec->vwc_sat*pow((-0.015/sitec->psi_sat),1.0/sitec->soil_b);
	//
	//	/* define maximum soilwater content, for outflow calculation
	//	converts volumetric water content (m3/m3) --> (kg/m2) */
	//	sitec->soilw_fc = sitec->soil_depth * sitec->vwc_fc * 1000.0;
	//	sitec->soilw_sat = sitec->soil_depth * sitec->vwc_sat * 1000.0;
	//}
	
	return (!ok);
}
int sitec_init2(file sitefile,siteconst_struct* sitec, epconst_struct* epc, cstate_struct* cs,
cinit_struct* cinit, nstate_struct* ns,co2control_struct* co2)
{

	int ok=1;
	double x,y,soilc_total=0,litrc_total=0,leafc2010=0,stemc2010=0,litrc2010=0,soilc2010=0,temp=0;
	//int climate_id;
	double sand,silt,clay,clim=0,tavg=0.0,soildepth=0.0,vegc88=0.0,rootc88=0.0,vegc10=0.0,rootc10=0.0; /* percent sand, silt, and clay */
	int soil;
	char junk[1000];
	int i;
	/*read sand percentage*/
	//printf("climate_id_site1=%d\n",&sitec->climate_id);		
	fscanf(sitefile.ptr,"%*lf,%*lf,%d,%*lf,%*lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,",
		&sitec->plotID,&y,&sitec->lat,&sitec->elev,&sand,&clay,&soildepth,&clim);
	fscanf(sitefile.ptr,"%lf,%lf,%lf,%lf,%lf,%*lf,%lf,%lf,%lf,%lf,%lf,%*lf,%d,",
		&cinit->max_leafc,&cinit->max_stemc,&cinit->croot1988,&litrc_total,&soilc_total,
		&cinit->leafc2010,&cinit->stemc2010,&cinit->croot2010,&cinit->litrc2010,&cinit->soilc2010,&cinit->frost);
	//土壤质地初始化，2014-08-04
	sitec->sand = sand;
	sitec->clay = clay;
	sitec->silt = 100.0 - sand - clay;
	//printf("sand =%lf,clay=%lf \n",sand,clay);
//	for(i=0;i<co2->co2vals;i++)	
//		printf("1 co2year = %d,co2ppm= %lf\n",co2->co2year_array[i],co2->co2ppm_array[i]);

//	if(sitec->epctype <16 && sitec->epctype >0)
//		soilc_total = 97.104 + 24.609 * log(vegc88) - 5.782 * (sitec->elev/1000.0) * (sitec->elev/1000.0) + 0.363 * soildepth - 7.924 * tavg;
	//printf("soilc88_new=%lf\n",soilc_total);
	//根据建模均值与取样均值调整;
	//soilc_total = soilc_total * 1.09;
	//cinit->soilc2010 = cinit->soilc2010 * 1.09;
	cinit->max_leafc /=10.0;	
	cinit->max_stemc /=10.0;
	litrc_total /=10.0;			
	soilc_total /=10.0;		
	sitec->climate_id=(int)clim;
	cinit->leafc2010 /=10.0;	
	cinit->stemc2010 /=10.0;
	cinit->litrc2010 /=10.0;	
	cinit->soilc2010 /=10.0;
	cinit->croot1988 /= 10.0;	
	cinit->croot2010 /=10.0;
	
	//printf("rootc2012=%lf\n",cinit->croot2010);
//	if(sitec->epctype>0 && sitec->epctype<16 && cinit->croot1988!=0 && cinit->max_stemc!=0)
//		epc->alloc_crootc_stemc = (cinit->croot1988 / cinit->max_stemc);
	/*soil depth*/
	soildepth /= 100.0;
	if(soildepth >= 1.0)
		soildepth = 1.0;
	if(soildepth < 0.20)
		soildepth = 0.20;
	sitec->soil_depth = soildepth;


	/*灌木、竹林、建筑用地等未处理的部分,2010年类型为0、16、17、18；*/
	switch((int)(sitec->epctype_change))
	{
	case 0:
		cinit->leafc2010 = 0.0; cinit->stemc2010 = 0.0; cinit->litrc2010 = 0.0; cinit->soilc2010 = 0.0;
		break;
	case 16:
		cinit->leafc2010 = 0.062; cinit->stemc2010 = 0.76; cinit->litrc2010 = 0.05; cinit->soilc2010 = 10.77;
		break;
	case 17:
		cinit->leafc2010 = 0.062; cinit->stemc2010 = 0.66; cinit->litrc2010 = 0.05; cinit->soilc2010 = 10.62;
		break;
	case 18:
		cinit->leafc2010 = 0.15; cinit->stemc2010 = 0.0; cinit->litrc2010 = 0.1; cinit->soilc2010 = 11.7;
		break;
	default:
		break;
	}
	switch((int)(sitec->epctype))
	{
	case 0:
		cinit->max_leafc = 0.0; cinit->max_stemc = 0.0; litrc_total = 0.0; soilc_total = 0.0;
		break;
	case 16:
		cinit->max_leafc = 0.062; cinit->max_stemc = 0.76; litrc_total = 0.05; soilc_total = 10.77;
		break;
	case 17:
		cinit->max_leafc = 0.062; cinit->max_stemc = 0.66; litrc_total = 0.05; soilc_total = 10.62;
		break;
	case 18:
		cinit->max_leafc = 0.15; cinit->max_stemc = 0.0; litrc_total = 0.1; soilc_total = 11.7;
		break;
	default:
		break;
	}
	/*处理异常值2010*/
	if(cinit->leafc2010 <0)		cinit->leafc2010 = 0.03;
	if(cinit->stemc2010 <0)		cinit->stemc2010 = 0.5;
	if(cinit->litrc2010 <0)		cinit->litrc2010 = 0.05;
	if(cinit->soilc2010 < 0)		cinit->soilc2010 = 2.0;
	if(cinit->leafc2010 > 2.5)		cinit->leafc2010 = 2.5;
	if(cinit->stemc2010 > 40)		cinit->stemc2010 = 40.0;
	if(cinit->litrc2010 > 3)		cinit->litrc2010 = 3.0;
	if(cinit->soilc2010 > 100)		cinit->soilc2010 = 100.0;
	/*处理异常值1988*/
	if(cinit->max_leafc < 0)		cinit->max_leafc = 0.03;
	if(cinit->max_stemc <0)		cinit->max_stemc = 0.5;
	/*乔木叶子最大1kgc/m2，灌木竹林最大0.5，草地最大0.2*/
	if(cinit->max_leafc >2.5)		cinit->max_leafc = 2.5;
	/*乔木树干最大40kgc/m2，灌木竹林最大20，草地为0*/
	if(cinit->max_stemc >40)		cinit->max_stemc = 40;
	if(litrc_total < 0)			litrc_total = 0.05;
	if(soilc_total < 0)			soilc_total = 2.0;
	if(litrc_total > 3)			litrc_total = 3.0;
	if(soilc_total > 100)			soilc_total = 100;
	//printf("max_stemc=%lf\n",cinit->max_stemc);
	//printf("max_leafc=%lf\n",cinit->max_leafc);
	//printf("epc=%lf\n",sitec->epctype);
	//printf("climate_id_site=%d\n",sitec->climate_id);
	silt = 100 - sand - clay;  //2013-03-26
	//printf("sand=%lf,clay=%lf,silt=%lf",sand,clay,silt);
	//fscanf(sitefile.ptr,"%lf",&ns->litr1n);
	//fscanf(sitefile.ptr,"%lf",&ns->sminn);
	ns->litr1n = 0.0002;
	ns->sminn = 0.0;
	if(sitec->epctype < 11)//wood
	{
		/*cs->cwdc = litrc_total * 0.87;
		ns->cwdn = cs->cwdc/epc->deadwood_cn;
		cs->litr1c=litrc_total * 0.003;
		cs->litr2c=litrc_total * 0.034;
		cs->litr3c=litrc_total * 0.024;
		cs->litr4c=litrc_total * 0.069;*/

		cs->cwdc = litrc_total * 0.5;
		ns->cwdn = cs->cwdc/epc->deadwood_cn;
		cs->litr1c=litrc_total * 0.05;
		cs->litr2c=litrc_total * 0.1;
		cs->litr3c=litrc_total * 0.15;
		cs->litr4c=litrc_total * 0.2;

		ns->litr2n = cs->litr2c / epc->leaflitr_cn;
		ns->litr3n = cs->litr3c / epc->leaflitr_cn;
		ns->litr4n = cs->litr4c / epc->leaflitr_cn;


		cs->soil1c=soilc_total * 0.001;
		cs->soil2c=soilc_total * 0.010;
		cs->soil3c=soilc_total * 0.129;
		cs->soil4c=soilc_total * 0.860;
		/*
		cs->soil1c=soilc_total * 0.005;
		cs->soil2c=soilc_total * 0.015;
		cs->soil3c=soilc_total * 0.18;
		cs->soil4c=soilc_total * 0.80;*/

		ns->soil1n = cs->soil1c/SOIL1_CN;
		ns->soil2n = cs->soil2c/SOIL2_CN;
		ns->soil3n = cs->soil3c/SOIL3_CN;
		ns->soil4n = cs->soil4c/SOIL4_CN;
		//printf("leafc=%lf\n",cs->soil4c);
	}
	else//grass
	{
		cs->cwdc=0;
		ns->cwdn = cs->cwdc/epc->deadwood_cn;
		cs->litr1c=litrc_total * 0.069;
		cs->litr2c=litrc_total * 0.319;
		cs->litr3c=litrc_total * 0.189;
		cs->litr4c=litrc_total * 0.423;
		ns->litr2n = cs->litr2c / epc->leaflitr_cn;
		ns->litr3n = cs->litr3c / epc->leaflitr_cn;
		ns->litr4n = cs->litr4c / epc->leaflitr_cn;
		cs->soil1c=soilc_total * 0.0037;
		cs->soil2c=soilc_total * 0.0133;
		cs->soil3c=soilc_total * 0.1334;
		cs->soil4c=soilc_total * 0.8496;
		ns->soil1n = cs->soil1c/SOIL1_CN;
		ns->soil2n = cs->soil2c/SOIL2_CN;
		ns->soil3n = cs->soil3c/SOIL3_CN;
		ns->soil4n = cs->soil4c/SOIL4_CN;
	}
	//	for(i=0;i<co2->co2vals;i++)	
	//		printf("3 co2year = %d,co2ppm= %lf\n",co2->co2year_array[i],co2->co2ppm_array[i]);



	/* calculate the soil pressure-volume coefficients from texture data */
	/* Uses the multivariate regressions from Cosby et al., 1984，
	Cosby, B. J., Hornberger, G.M., Clapp, R. B., and Ginn, T. R.: 1984, 
	‘A Statistical Exploration of the Relationships of Soil Moisture Characteristics to 
	the Physical Properties of Soils’, Water Resour. Res. 20, 682–690. */
	/* first check that the percentages add to 100.0 */
	if (ok && fabs(sand+silt+clay-100) > FLT_COND_TOL)
	{
		bgc_printf(BV_ERROR, "Error: %%sand + %%silt + %%clay  MUST EQUAL 100%%\n");
		bgc_printf(BV_ERROR, "Check values in initialization file.\n");
		ok=0;
	}
	if (ok)
	{
		sitec->soil_b = -(3.10 + 0.157*clay - 0.003*sand);
		sitec->vwc_sat = (50.5 - 0.142*sand - 0.037*clay)/100.0;
		sitec->psi_sat = -(exp((1.54 - 0.0095*sand + 0.0063*silt)*log(10.0))*9.8e-5);
		sitec->vwc_fc = sitec->vwc_sat*pow((-0.015/sitec->psi_sat),1.0/sitec->soil_b);
	
		/* define maximum soilwater content, for outflow calculation
		converts volumetric water content (m3/m3) --> (kg/m2) */
		sitec->soilw_fc = sitec->soil_depth * sitec->vwc_fc * 1000.0;
		sitec->soilw_sat = sitec->soil_depth * sitec->vwc_sat * 1000.0;
	}
//		for(i=0;i<co2->co2vals;i++)	
//			printf("4 co2year = %d,co2ppm= %lf\n",co2->co2year_array[i],co2->co2ppm_array[i]);

	return (!ok);
}
/*back up of sitec_init2*/
int sitec_init3(file sitefile,siteconst_struct* sitec)
{
	int ok=1;
	double x,y;
	double sand,silt,clay; /* percent sand, silt, and clay */
	int soil;
	char junk[10];
	/*read sand percentage*/
	fscanf(sitefile.ptr,"%s%*[^\n]",junk);
	fscanf(sitefile.ptr,"%lf",&sitec->lat);
	fscanf(sitefile.ptr,"%lf",&y);
	fscanf(sitefile.ptr,"%lf",&sitec->elev);
	fscanf(sitefile.ptr,"%lf",&sitec->soiltype);
    
	
	if (ok && scan_value(sitefile, &sitec->epctype, 'd'))
	{
		bgc_printf(BV_ERROR, "Error reading epctype, sitec_init()\n");
		ok=0;
	}

	soil=floor(sitec->soiltype);
	switch(soil)
	{
	case 1:
		sand=85;
		clay=10;
		silt=5;
		break;
	case 2:
		sand=50;
		clay=35;
		silt=15;
		break;
	case 3:
		sand=15;
		clay=25;
		silt=60;
		break;
	case 4:
		sand=60;
		clay=30;
		silt=10;
		break;
	case 5:
		sand=55;
		clay=5;
		silt=40;
		break;
	case 6:
		sand=10;
		clay=45;
		silt=45;
		break;
	case 7:
		sand=30;
		clay=35;
		silt=35;
		break;
	case 10:
		sand=30;
		clay=20;
		silt=50;
		break;
	default:
		break;
	}

	/* calculate the soil pressure-volume coefficients from texture data */
	/* Uses the multivariate regressions from Cosby et al., 1984 */
	/* first check that the percentages add to 100.0 */
	if (ok && fabs(sand+silt+clay-100) > FLT_COND_TOL)
	{
		bgc_printf(BV_ERROR, "Error: %%sand + %%silt + %%clay  MUST EQUAL 100%%\n");
		bgc_printf(BV_ERROR, "Check values in initialization file.\n");
		ok=0;
	}
	if (ok)
	{
		sitec->soil_b = -(3.10 + 0.157*clay - 0.003*sand);
		sitec->vwc_sat = (50.5 - 0.142*sand - 0.037*clay)/100.0;
		sitec->psi_sat = -(exp((1.54 - 0.0095*sand + 0.0063*silt)*log(10.0))*9.8e-5);
		sitec->vwc_fc = sitec->vwc_sat*pow((-0.015/sitec->psi_sat),1.0/sitec->soil_b);
	
		/* define maximum soilwater content, for outflow calculation
		converts volumetric water content (m3/m3) --> (kg/m2) */
		sitec->soilw_fc = sitec->soil_depth * sitec->vwc_fc * 1000.0;
		sitec->soilw_sat = sitec->soil_depth * sitec->vwc_sat * 1000.0;
	}
	
	return (!ok);
}