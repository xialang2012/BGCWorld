/*
metarr_init.c
Initialize meteorological data arrays for pointbgc simulation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions from version 4.1.1:
The ratio of PAR:shortwave radiation is now being used as a
macro definition from bgc_constants.h (RAD2PAR)
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "pointbgc.h"

/* It is assumed here that the meteorological datafile contains the following 
list of variables, with the indicated units.  All other variables in the metv
arrays are derived from this basic set:

VARIABLE    UNITS
yday        (none) (yearday)
prcp        cm     (daily total precipitation, water equivalent)  
tmax        deg C  (daily maximum temperature) 
tmin        deg C  (daily minimum temperature) 
VPD         Pa     (daylight average VPD)
swavgfd     W/m2   (daylight average shortwave flux density)
daylength   s      (daylight duration)

*/

int metarr_init(file metf, metarr_struct* metarr, const climchange_struct* scc,
int nyears) 
{
	
	int ok = 1;
	int i;
	int ndays;
	int year;
	double tmax,tmin,prcp,vpd,swavgfd,dayl,tair_avg,tavg,ws,rh,lrad;
	char junk_head[100];
	//added 临时
	int station,month,day;
	ndays = 365 * nyears;

	/* allocate space for the metv arrays */
	if (ok && !(metarr->tmax = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for tmax array\n");
		ok=0;
	}
	if (ok && !(metarr->tmin = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for tmin array\n");
		ok=0;
	}
	if (ok && !(metarr->prcp = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for prcp array\n");
		ok=0;
	}
	if (ok && !(metarr->vpd = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for vpd array\n");
		ok=0;
	}
	if (ok && !(metarr->tavg = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for tavg array\n");
		ok=0;
	}
	if (ok && !(metarr->tavg_ra = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for tavg_ra array\n");
		ok=0;
	}
	if (ok && !(metarr->swavgfd = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for swavgfd array\n");
		ok=0;
	}
	if (ok && !(metarr->par = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for par array\n");
		ok=0;
	}
	if (ok && !(metarr->dayl = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for dayl array\n");
		ok=0;
	}
	//added
	if (ok && !(metarr->rhumid = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for rh array\n");
		ok=0;
	}
	if (ok && !(metarr->wspeed = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for wind speed array\n");
		ok=0;
	}
	if (ok && !(metarr->lrad = (double*) malloc(ndays * sizeof(double))))
	{
		bgc_printf(BV_ERROR, "Error allocating for lrad array\n");
		ok=0;
	}
	/*read and discard header of met file, move from met_ini. 10,2010
	
	for (i = 0 ; ok && i<1 ; i++)
	{
		if (scan_value(metf, junk_head, 's'))
		{
			bgc_printf(BV_ERROR, "Error reading met file header line #%d\n",i+1);
			ok=0;
		}
	}*/
	/* begin daily loop: read input file, generate array values */
	for(i = 0 ; ok && i<ndays ; i++)
	{
		/*读取插值的气候数据*/
		if (ok && fscanf(metf.ptr,"%lf%lf%lf%lf%lf%lf",&tmax,&tmin,&tavg,&prcp,&vpd,&swavgfd)==EOF)
		{
			bgc_printf(BV_ERROR, "Error reading met file, metarr_init()\n");
			ok=0;
		}

		if (ok && fscanf(metf.ptr,"%lf%lf%lf",&dayl,&ws,&lrad)==EOF)
		{
			bgc_printf(BV_ERROR, "Error reading met file, metv_init()\n");
			ok=0;
		}
		
				/* Fixed 02/05/04 
		if( swavgfd < 0.0 )
		{
			swavgfd = 0.0;
		}

		if( dayl < 0.0 )
		{
			dayl = 0.0;
		}*/
		/* apply the climate change scenario and store 
		metarr->tmax[i] = tmax + scc->s_tmax;
		metarr->tmin[i] = tmin + scc->s_tmin;
		metarr->prcp[i] = prcp * scc->s_prcp;
		metarr->vpd[i] = vpd * scc->s_vpd;
		metarr->swavgfd[i] = swavgfd * scc->s_swavgfd;
		metarr->par[i] = swavgfd * RAD2PAR * scc->s_swavgfd;
		metarr->dayl[i] = dayl;
		metarr->tavg[i] = (metarr->tmax[i] + metarr->tmin[i]) / 2.0;

		//added 根据VPD计算rh
		metarr->rhumid[i] = 1 - (metarr->vpd[i]/(0.61078 * 1000 * exp(17.269 * metarr->tavg[i]/(237.3 + metarr->tavg[i]))));
*/
		/*added 临时修改气候数据读取,读取气象站点的测试数据
		if (ok && fscanf(metf.ptr,"%d%d%d%d",&station,&year,&month,&day)==EOF)
		{
			bgc_printf(BV_ERROR, "Error reading met file, metarr_init()\n");
			ok=0;
		}
		if (ok && fscanf(metf.ptr,"%lf%lf%lf%lf%lf%*lf%*lf%*lf%*lf%lf%*lf%lf%lf%*lf",&tavg,&tmax,&tmin,&ws,&prcp,&rh,&dayl,&swavgfd)==EOF)
		{
			bgc_printf(BV_ERROR, "Error reading met file, metarr_init()\n");
			ok=0;
		}*/

		
	//	vpd = (1 - rh * 0.01) * 0.61078 * exp(17.269 * 0.1 * tavg/(237.3 + tavg*0.1)) * 1000;

		metarr->tmax[i] = tmax;
		metarr->tmin[i] = tmin;
		metarr->tavg[i] = tavg;
		metarr->prcp[i] = prcp * 0.1;//mm to cm
		metarr->vpd[i] = vpd;
		metarr->swavgfd[i] = swavgfd;
		metarr->par[i] = swavgfd * RAD2PAR;
		//metarr->dayl[i] = dayl * 360; //0.1h=>s
		metarr->dayl[i] = dayl;
		metarr->wspeed[i] = ws; //m/s
		//metarr->rhumid[i] = rh;
		metarr->rhumid[i] = 1 - (metarr->vpd[i]/(0.61078 * 1000 * exp(17.269 * metarr->tavg[i]/(237.3 + metarr->tavg[i]))));
		metarr->lrad[i] = lrad;
		//added 根据VPD计算rh
		//metarr->rhumid[i] = 1 - (metarr->vpd[i]/(0.61078 * 1000 * exp(17.269 * metarr->tavg[i]/(237.3 + metarr->tavg[i]))));

		//printf("rhumid= %lf,ws=%lf,prec=%lf,tmax=%lf, rad=%lf\n",metarr->rhumid[i],metarr->wspeed[i],metarr->prcp[i],tmax,metarr->swavgfd[i]);
		//printf("tmax= %lf,tmin=%lf,tavg=%lf,prcp=%lf, vpd=%lf, rad=%lf, dayl=%lf,ws=%lf,lrad=%lf \n"
		//	,metarr->tmax[i],metarr->tmin[i],metarr->tavg[i],metarr->prcp[i],metarr->vpd[i],metarr->swavgfd[i],metarr->dayl[i],metarr->wspeed[i],metarr->lrad[i]);
	}
	
	/* perform running averages of daily average temperature for 
	use in soil temperature routine. 

	This implementation uses a linearly ramped 11-day running average 
	of daily mean air temperature, with days 1-10 based on a 1-10 day
	running average, respectively. 
	*/
	
	if (ok && run_avg(metarr->tavg, metarr->tavg_ra, ndays, 11, 1))
	{
		bgc_printf(BV_ERROR, "Error: run_avg() in metv_init.c \n");
		ok = 0;
	}
	
	return (!ok);
}
