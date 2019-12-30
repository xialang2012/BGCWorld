/* 
maint_resp.c
daily maintenance respiration

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int maint_resp1(const cstate_struct* cs, const nstate_struct* ns,
const epconst_struct* epc, const metvar_struct* metv, cflux_struct* cf,
epvar_struct* epv)
{
	/*
	maintenance respiration routine
	Uses reference values at 20 deg C and an empirical relationship between
	tissue N content and respiration rate given in:

	Ryan, M.G., 1991. Effects of climate change on plant respiration.
	Ecological Applications, 1(2):157-167.
	
	Uses the same value of Q_10 (2.0) for all compartments, leaf, stem, 
	coarse and fine roots.
	
	From Ryan's figures and regressions equations, the maintenance respiration
	in kgC/day per kg of tissue N is:
	mrpern = 0.218 (kgC/kgN/d)
	
	Leaf maintenance respiration is calculated separately for day and
	night, since the PSN routine needs the daylight value.
	
	Leaf and fine root respiration are dependent on phenology.
	*/
	
	int ok=1;
	double t1;
	double q10 = 2.0;
	double mrpern = 0.218;
	double exponent;
	double n_area_sun, n_area_shade, dlmr_area_sun, dlmr_area_shade;
	
	/*对于云冷杉，调低维持呼吸 2013-07-19*/

	/* leaf day and night maintenance respiration when leaves on */
	if (cs->leafc)
	{
		t1 = ns->leafn * mrpern;
		
		/* leaf, day */
		exponent = (metv->tday - 20.0) / 10.0;
		cf->leaf_day_mr = t1 * pow(q10, exponent) * metv->dayl / 86400.0;

		/* for day respiration, also determine rates of maintenance respiration
		per unit of projected leaf area in the sunlit and shaded portions of
		the canopy, for use in the photosynthesis routine */
		/* first, calculate the mass of N per unit of projected leaf area
		in each canopy fraction (kg N/m2 projected area) */
		n_area_sun   = 1.0/(epv->sun_proj_sla * epc->leaf_cn);
		n_area_shade = 1.0/(epv->shade_proj_sla * epc->leaf_cn);
		/* convert to respiration flux in kg C/m2 projected area/day, and
		correct for temperature */
		dlmr_area_sun   = n_area_sun * mrpern * pow(q10, exponent);
		dlmr_area_shade = n_area_shade * mrpern * pow(q10, exponent);
		/* finally, convert from mass to molar units, and from a daily rate to 
		a rate per second */
		epv->dlmr_area_sun = dlmr_area_sun/(86400.0 * 12.011e-9);
		epv->dlmr_area_shade = dlmr_area_shade/(86400.0 * 12.011e-9);
		
		/* leaf, night */
		exponent = (metv->tnight - 20.0) / 10.0;
		cf->leaf_night_mr = t1 * pow(q10, exponent) * 
			(86400.0 - metv->dayl) / 86400.0;
	}
	else /* no leaves on */
	{
		cf->leaf_day_mr = 0.0;
		epv->dlmr_area_sun = 0.0;
		epv->dlmr_area_shade = 0.0;
		cf->leaf_night_mr = 0.0;
	}

	/* fine root maintenance respiration when fine roots on */
	/* ammended to consider only the specified n concentration,
	to avoid excessive MR with n-loading to fine roots */
	if (cs->frootc)
	{
		exponent = (metv->tsoil - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->froot_mr = ns->frootn * mrpern * t1;
	}
	else /* no fine roots on */
		cf->froot_mr = 0.0;

	/* TREE-specific fluxes */
	if (epc->woody)
	{
		/* live stem maintenance respiration */
		exponent = (metv->tavg - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->livestem_mr = ns->livestemn * mrpern * t1;

		/* live coarse root maintenance respiration */
		exponent = (metv->tsoil - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->livecroot_mr = ns->livecrootn * mrpern * t1;
	}
	
	return (!ok);
}
/* 
maint_resp.c
daily maintenance respiration

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int maint_resp(wstate_struct* ws,const siteconst_struct* sitec,const cstate_struct* cs, const nstate_struct* ns,
const epconst_struct* epc, const metvar_struct* metv, cflux_struct* cf,
epvar_struct* epv, const high_time_resolution& highT, std::vector<StationDataFlux*> &sfData, 
const int simyr, const int yday, int mode)
{
	// for comparison between the daily and hourly, 
	// 不需要使用小时数据计算替换天呼吸为false，默认为true，计算高时间分辨率模型时必须为true
	bool comPasison = true;

	/*
	maintenance respiration routine
	Uses reference values at 20 deg C and an empirical relationship between
	tissue N content and respiration rate given in:

	Ryan, M.G., 1991. Effects of climate change on plant respiration.
	Ecological Applications, 1(2):157-167.
	
	Uses the same value of Q_10 (2.0) for all compartments, leaf, stem, 
	coarse and fine roots.
	
	From Ryan's figures and regressions equations, the maintenance respiration
	in kgC/day per kg of tissue N is:
	mrpern = 0.218 (kgC/kgN/d)
	
	Leaf maintenance respiration is calculated separately for day and
	night, since the PSN routine needs the daylight value.
	
	Leaf and fine root respiration are dependent on phenology.
	*/
	
	int ok=1;
	double t1;
	double q10 = 2.0;
	double mrpern = 0.218;
	double exponent;
	double n_area_sun, n_area_shade, dlmr_area_sun, dlmr_area_shade;
	double t2 = 1.0,t3 = 1.0,w=1;
	double soilw_vwc;
	/************************
	if(epc->evergreen)
		t2 = 0.026 * sitec->lat - 0.375;
	else
		t2 = 1.0;
	if(sitec->epctype==1)
	{
		//t2=0.8;
		//q10=1.9;
	}*/
	int timeRes = 0;
	int dayCurr = 0;

	
	if(comPasison && highT.active)
	{
		dayCurr = simyr * 365 + yday;
		timeRes = sfData[1]->HRMIN - sfData[0]->HRMIN;
	}

	/* leaf day and night maintenance respiration when leaves on */
	if (cs->leafc)
	{
		t1 = ns->leafn * mrpern;
		
		/* leaf, day */
		exponent = (metv->tday - 20.0) / 10.0;
		cf->leaf_day_mr = t2 * t1 * pow(q10, exponent) * metv->dayl / 86400.0;

		/* for day respiration, also determine rates of maintenance respiration
		per unit of projected leaf area in the sunlit and shaded portions of
		the canopy, for use in the photosynthesis routine */
		/* first, calculate the mass of N per unit of projected leaf area
		in each canopy fraction (kg N/m2 projected area) */
		n_area_sun   = 1.0/(epv->sun_proj_sla * epc->leaf_cn);
		n_area_shade = 1.0/(epv->shade_proj_sla * epc->leaf_cn);
		/* convert to respiration flux in kg C/m2 projected area/day, and
		correct for temperature */
		dlmr_area_sun   = n_area_sun * mrpern * pow(q10, exponent);
		dlmr_area_shade = n_area_shade * mrpern * pow(q10, exponent);
		/* finally, convert from mass to molar units, and from a daily rate to 
		a rate per second */
		epv->dlmr_area_sun = dlmr_area_sun/(86400.0 * 12.011e-9);
		epv->dlmr_area_shade = dlmr_area_shade/(86400.0 * 12.011e-9);
		
		/* leaf, night */
		exponent = (metv->tnight - 20.0) / 10.0;
		cf->leaf_night_mr =t2 * t1 * pow(q10, exponent) * (86400.0 - metv->dayl) / 86400.0;

		// for high time MR
		if (comPasison && highT.active && mode == MODE_MODEL)
		{
			cf->leaf_day_mr = 0.0;
			cf->leaf_night_mr = 0.0;
			for (int i = 0; i < 24 * 60 / timeRes; ++i)
			{
				double temp= sfData[dayCurr * 24 * 60 / timeRes + i]->TA;
				exponent = (temp - 20.0) / 10.0;

				if (sfData[dayCurr * 24 * 60 / timeRes + i]->Srad > 0)	// daytime
				{					
					cf->leaf_day_mr += t2 * t1 * pow(q10, exponent) * timeRes *60 / 86400.0;
				}
				else  // night
				{
					cf->leaf_night_mr = t2 * t1 * pow(q10, exponent) * timeRes * 60 / 86400.0;
				}

			}			
		}
	}
	else /* no leaves on */
	{
		cf->leaf_day_mr = 0.0;
		epv->dlmr_area_sun = 0.0;
		epv->dlmr_area_shade = 0.0;
		cf->leaf_night_mr = 0.0;
	}

	/* fine root maintenance respiration when fine roots on */
	/* ammended to consider only the specified n concentration,
	to avoid excessive MR with n-loading to fine roots */
	if (cs->frootc)
	{
		exponent = (metv->tsoil - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->froot_mr = ns->frootn * mrpern * t1 * t2;
	}
	else /* no fine roots on */
		cf->froot_mr = 0.0;

	/* TREE-specific fluxes */
	if (epc->woody)
	{
		/* live stem maintenance respiration */
		exponent = (metv->tavg - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->livestem_mr = ns->livestemn * mrpern * t1 * t2;

		/* live coarse root maintenance respiration */
		exponent = (metv->tsoil - 20.0) / 10.0;
		t1 = pow(q10, exponent);
		cf->livecroot_mr = ns->livecrootn * mrpern * t1 * t2;

		// for high time MR
		if (comPasison && highT.active && mode == MODE_MODEL)
		{
			cf->livestem_mr = 0.0;
			for (int i = 0; i < 24 * 60 / timeRes; ++i)
			{
				double temp = sfData[dayCurr * 24 * 60 / timeRes + i]->TA;
				exponent = (temp - 20.0) / 10.0;
				t1 = pow(q10, exponent);
				cf->livestem_mr += ns->livestemn * mrpern * t1 * t2 * timeRes * 60 / 86400.0;				
			}
		}
	}


	//add 2014-12-08；假如类型为湿地，则调整湿地呼吸的系数；
	if(sitec->epctype ==13)
	{
		soilw_vwc = (ws->soilw/(1000*sitec->soil_depth))/sitec->vwc_sat;
		if(soilw_vwc >= 0.7)
			//w = (1 - (ws->soilw/(1000*sitec->soil_depth))/sitec->vwc_sat) / 0.3; // y = -2.5 * x +2.5
			//w = (13 - 10 * soilw_vwc) /6.0;    //体积含水量的百分比，饱和为1。呼吸为0.5倍，为70%饱和时，系数为1.0；
			w = (1 - soilw_vwc) /0.3;
		if(soilw_vwc < 0.7)
			w = 1.0;	
	//	printf("w=%lf\n",w);
	}

	if(sitec->epctype ==13)
	{
		if (cs->frootc)
		{
			exponent = (metv->tsoil - 20.0) / 10.0;
			t1 = pow(q10, exponent);
			cf->froot_mr = ns->frootn * mrpern * t1 * t2 * w;
		}
		else /* no fine roots on */
		cf->froot_mr = 0.0;

		if (epc->woody)
		{
			 /*live stem maintenance respiration */
			exponent = (metv->tavg - 20.0) / 10.0;
			t1 = pow(q10, exponent);
			cf->livestem_mr = ns->livestemn * mrpern * t1 * t2;

			/* live coarse root maintenance respiration */
			exponent = (metv->tsoil - 20.0) / 10.0;
			t1 = pow(q10, exponent);
			cf->livecroot_mr = ns->livecrootn * mrpern * t1 * t2 * w;
		}
	
	}
	
	return (!ok);
}
