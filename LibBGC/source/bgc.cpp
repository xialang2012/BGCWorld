/*
bgc.c
Core BGC model logic

Includes in-line output handling routines that write to daily and annual
output files. This is the only library module that has external
I/O connections, and so it is the only module that includes
bgc_io.h.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions since 4.1.2
	Merged spinup_bgc.c with bgc.c to eliminate
	code duplication
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

/* These DEBUG defines are now depricated. Please use 
   bgc_printf(BV_DIAG,...) instead. The only place where 
/* #define DEBUG */
/* #define DEBUG_SPINUP set this to see the spinup details on-screen */

/*	SANE = Do 'Pan-Arctic' style summary. INSANE is traditional style 
		summary. See the '-p' cli flag in USAGE.TXT */
signed char summary_sanity = INSANE ;
 
int bgc(bgcin_struct* bgcin, bgcout_struct* bgcout, int mode)
{
	// read high time resolution station data
	if (bgcin->hModel.active && mode == MODE_MODEL)
	{
		bgcin->hModel.tmpHighFile.open("./outputs/" + std::to_string( bgcin->sitec.climate_id ) + "_Psn_High.csv");
		bgcin->hModel.tmpHighFile << "year, " << "daya, " << "temperature, " << "temperature_High, " << "psnsun_cpool, " <<
			"psnsun_cpool_High, " << "psnshade_cpool, " << "psnshade_cpool_High, " << "carbon limited, " << "carbon limited Per" << std::endl;

		if (bgcin->hModel.output_stress)
		{
			bgcin->hModel.highFile_stress.open("./outputs/" + std::to_string(bgcin->sitec.climate_id) + "_Stress_High.csv");
		}

		if (bgcin->hModel.output_carbon)
		{
			bgcin->hModel.highFile_carbon.open("./outputs/" + std::to_string(bgcin->sitec.climate_id) + "_Carbon_High.csv");
		}
	}
	else if(mode == MODE_MODEL)
	{
		bgcin->hModel.tmpHighFile.open("./outputs/" + std::to_string(bgcin->sitec.climate_id) + "_Psn.csv");
		bgcin->hModel.tmpHighFile << "year, " << "daya, " << "temperature, " << "temperature_High(invalid), " << "psnsun_cpool, " <<
			"psnsun_cpool_High(invalid), " << "psnshade_cpool, " << "psnshade_cpool_High(invalid), " << "carbon limited, " << "carbon limited Per(invalid)" << std::endl;
	}

	extern signed char summary_sanity;
	/* variable declarations */
	int ok=1;

	/* iofiles and program control variables */
	control_struct     ctrl;

	/* meteorological variables */
	metarr_struct      metarr;
	metvar_struct      metv;
	co2control_struct  co2;
	ramp_ndep_struct ramp_ndep;
	
	/* state and flux variables for water, carbon, and nitrogen */
	wstate_struct      ws, zero_ws;
	wflux_struct       wf, zero_wf;
	cinit_struct       cinit;
	cstate_struct      cs, zero_cs;
	cflux_struct       cf, zero_cf;
	nstate_struct      ns, zero_ns;
	nflux_struct       nf, zero_nf;

	/* primary ecophysiological variables */
	epvar_struct       epv;

	/* site physical constants */
	siteconst_struct   sitec;

	/* phenological data */
	phenarray_struct   phenarr;
	phenology_struct   phen;

	/* ecophysiological constants */
	epconst_struct     epc;

	/* photosynthesis constructs */
	psn_struct         psn_sun, psn_shade;

	/* temporary nitrogen variables for decomposition and allocation */
	ntemp_struct       nt;
	
	/* summary variable structure */
	summary_struct     summary;
	/*添加气候变量结构体，用于存储半小时的气候*/
	pmet_struct pmetvar; //add
	
	/* output mapping (array of pointers to double) */
	double **output_map;
	
	/* local storage for daily and annual output variables */
	float *dayarr, *monavgarr, *annavgarr, *annarr;

	/* miscelaneous variables for program control in main */
	int simyr, yday, metyr, metday;
	int first_balance;
	int annual_alloc;
	int outv;
	int i, nmetdays;
	double tair_avg, tdiff;
	int dayout;

	/* mode == MODE_MODEL only */
	double daily_ndep, daily_nfix, ndep_scalar, ndep_diff, ndep;
	int ind_simyr;
	
	/* variables used for monthly average output */
	int curmonth;
	int mondays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	int endday[12] = {30,58,89,119,150,180,211,242,272,303,333,364};
	float monmaxlai = 0.0,annmaxlai = 0.0,monmaxsnoww = 0.0;
	float eomsnoww = 0.0,eomsoilw = 0.0;

	int tmpyears;

	/* mode == MODE_SPINUP only */
	/* spinup control */
	int ntimesmet, nblock;
	int steady1, steady2, rising, metcycle = 0, spinyears;
	double tally1 = 0.0, tally1b = 0.0, tally2 = 0.0, tally2b = 0.0, t1 = 0.0;
	double naddfrac;
	
	/* mode == MODE_MODEL only */
	/* simple annual variables for text output */
	double annmaxplai,annet,anntrans,anneva,anncanopy,annoutflow,annnpp,annnbp,annprcp,anntavg,ann_nep,
		   ann_gpp,ann_nee,ann_mr,ann_gr,ann_hr,ann_root_mr,ann_leaf_mr;
	double ann_soilw,ann_snoww,ann_FSM_soilw,ann_FSM_icew,ann_FSM_snoww,value_c,yN=0;
//	int typec_day,typec_ok=0;//tc is type change
	char epctype_char[100];//2013-06-27
	int write_title=1;
	file epc_file;
	/*new added struct 2014-07-30*/
	//soilpar_struct *vspar;           /* add */
	//soilvar_struct svar[N];         /* add */
	soilvar_struct soil_sum[N];	/* add */
	FILE *fp1;						/* add */
	double tsoil,tsoil_sum,epv_psi,epv_psi_sum;/* add */
	char pointID[50],soilT_path[100]="./bgc-data/outputs/soilT_";
	int phenology_lastday=0;

	_itoa(bgcin->sitec.pointcounter,pointID,10);	
	strcat(soilT_path,pointID);
	strcat(soilT_path,".csv");
	//printf(soilT_path);
	if(bgcin->cinit.frost)
	{
		//fp1 = fopen(soilT_path,"w");
		/*fprintf(fp1,"DOY,空气温度,bgc土壤温度,FSM土壤温度,雪层温度,bgc土壤雪,FSM_snowW,FSM土壤水,FSM土壤冰,FSM冰+水,bgc土壤水分,bgc水势,FSM水势,FSM_lrad_up,FSM_lrad_down,bgc_lrad_input,");
		for(i = 0;i < N;i++)
		{
			fprintf(fp1,"%d Tavg,",i);
		}
		for(i = 0;i < N;i++)
		{
			fprintf(fp1,"%d Water,",i);
		}
		fprintf(fp1,"\n");*/
	}
	//fprintf(fp1,"bgc_air_T,bgc_soil_T,snow_T,soil_T1m,bgc_水势,FSM_水势,bgc_饱和水势,FSM_饱和水势,FSM_water,sitec->vwc_sat,");
	
	

	/*if(!typec_ok)
	{
		typec_day = (int)(rand()%8395);
		//printf("转化的日期%d\n",typec_day);
		typec_ok=1;
	}*/	
	
	/*****************输出气孔导度********************/
	if (mode != MODE_SPINUP && mode != MODE_MODEL)
	{
		bgc_printf(BV_ERROR, "Error: Unknown MODE given when calling bgc()\n");
		ok=0;
	}
	
	/* copy the input structures into local structures */
	ws = bgcin->ws;
	cinit = bgcin->cinit;
	cs = bgcin->cs;
	ns = bgcin->ns;
	sitec = bgcin->sitec;
	epc = bgcin->epc;
	/* note that the following three structures have dynamic memory elements,
	and so the notion of copying the input structure to a local structure
	value-by-value is not the same as above. In this case, the array pointers
	are being copied, so the local members use the same memory that was
	allocated in the calling function. Note also that bgc() does not modify
	the contents of these structures. */
	ctrl = bgcin->ctrl;
	metarr = bgcin->metarr;
	co2 = bgcin->co2;

	/*初始化冻土模块*/
	sparameters(&bgcin->vspar);	
	//printf("dddd\n");
//	printf("metsoilT= %lf\n",metarr.tavg_ra[0]);
	svar_init(&metarr);
		
	/*for(i=0;i<co2.co2vals;i++)
	{
		printf("bgcin->co2 = %d,co2= %d\n",bgcin->co2.co2year_array[i],co2.co2year_array[i]);
		printf("bgcin->co2ppm = %lf,co2ppm= %lf\n",bgcin->co2.co2ppm_array[i],co2.co2ppm_array[i]);
	}*/
	if (mode == MODE_MODEL)
	{
		ramp_ndep = bgcin->ramp_ndep;
	}
	
	bgc_printf(BV_DIAG, "done copy input\n");

	/* local variable that signals the need for daily output array */
	dayout = (ctrl.dodaily || ctrl.domonavg || ctrl.doannavg);

	/* allocate memory for local output arrays */
	if (ok && dayout &&	!(dayarr = (float*) malloc(ctrl.ndayout * sizeof(float))))
	{
		bgc_printf(BV_ERROR, "Error allocating for local daily output array in bgc()\n");
		ok=0;
	}
	if (ok && ctrl.domonavg && !(monavgarr = (float*) malloc(ctrl.ndayout * sizeof(float))))
	{
		bgc_printf(BV_ERROR, "Error allocating for monthly average output array in bgc()\n");
		ok=0;
	}
	if (ok && ctrl.doannavg && !(annavgarr = (float*) malloc(ctrl.ndayout * sizeof(float))))
	{
		bgc_printf(BV_ERROR, "Error allocating for annual average output array in bgc()\n");
		ok=0;
	}
	if (ok && ctrl.doannual && !(annarr = (float*) malloc(ctrl.nannout * sizeof(float))))
	{
		bgc_printf(BV_ERROR, "Error allocating for local annual output array in bgc()\n");
		ok=0;
	}
	/* allocate space for the output map pointers */
	if (ok && !(output_map = (double**) malloc(NMAP * sizeof(double*))))
	{
		bgc_printf(BV_ERROR, "Error allocating for output map in output_map_init()\n");
		ok=0;
	}
	
	bgc_printf(BV_DIAG, "done allocate out arrays\n");
	
	/* initialize monavg and annavg to 0.0 */
	if (ctrl.domonavg)
	{
		for (outv=0 ; outv<ctrl.ndayout ; outv++)
		{
			monavgarr[outv] = 0.0;
		}
	}
	if (ctrl.doannavg)
	{
		for (outv=0 ; outv<ctrl.ndayout ; outv++)
		{
			annavgarr[outv] = 0.0;
		}
	}

	/* initialize the output mapping array */
	if (ok && output_map_init(output_map,&metv,&ws,&wf,&cs,&cf,&ns,&nf,&phen,
		&epv,&psn_sun,&psn_shade,&summary))
	{
		bgc_printf(BV_ERROR, "Error in call to output_map_init() from bgc()\n");
		ok=0;
	}
	
	bgc_printf(BV_DIAG, "done initialize outmap\n");
	
	/* make zero-flux structures for use inside annual and daily loops */
	if (ok && make_zero_flux_struct(&zero_wf, &zero_cf, &zero_nf))
	{
		bgc_printf(BV_ERROR, "Error in call to make_zero_flux_struct()0 from bgc()\n");
		ok=0;
	}
	
	bgc_printf(BV_DIAG, "done make_zero_flux\n");
	
	/* atmospheric pressure (Pa) as a function of elevation (m) */
	if (ok && atm_pres(sitec.elev, &metv.pa))
	{
		bgc_printf(BV_ERROR, "Error in atm_pres() from bgc()\n");
		ok=0;
	}
	
	bgc_printf(BV_DIAG, "done atm_pres\n");
	
	// init gsi
	if (bgcin->gsiM.active)
	{
		if (!initgsiData(bgcin->gsiM.gsiFile, &epc)) return -1;
		GSI_calculation(&metarr, &sitec, &epc, &phenarr, &ctrl);
	}

	/* determine phenological signals */
	if (ok && prephenology(bgcin->gsiM, &ctrl, &epc, &sitec, &metarr, &phenarr))
	{
		bgc_printf(BV_ERROR, "Error in call to prephenology(), from bgc()\n");
		ok=0;
	}
	
	bgc_printf(BV_DIAG, "done prephenology\n");
	
	/* calculate the annual average air temperature for use in soil 
	temperature corrections. This code added 9 February 1999, in
	conjunction with soil temperature testing done with Mike White. */
	tair_avg = 0.0;
	nmetdays = ctrl.metyears * 365;
	for (i=0 ; i<nmetdays ; i++)
	{
		tair_avg += metarr.tavg[i];
	}
	tair_avg /= (double)nmetdays;
	
	/* if this simulation is using a restart file for its initial
	conditions, then copy restart info into structures */
	if (ok && ctrl.read_restart)
	{
		if (ok && restart_input(&ctrl, &ws, &cs, &ns, &epv, &metyr,
			&(bgcin->restart_input)))
		{
			bgc_printf(BV_ERROR, "Error in call to restart_input() from bgc()\n");
			ok=0;
		}
		
		bgc_printf(BV_DIAG, "done restart_input\n");
	
	}
	else
	/* no restart file, user supplies initial conditions */
	{
		/* initialize leaf C and N pools depending on phenology signals for
		the first metday */
		if (ok && firstday(&epc, &cinit, &epv, &phenarr, &cs, &ns))
		{
			bgc_printf(BV_ERROR, "Error in call to firstday(), from bgc()\n");
			ok=0;
		}
		//printf("deadcrootc=%lf, deadcrootc_storage=%lf, deadcrootc_transfer=%lf\n",cs.deadcrootc,cs.deadcrootc_storage,cs.deadcrootc_transfer);
		/* initial value for metyr 2013-03-26*/
		metyr = 0;//气候数据从1985年读取，实际模拟年份从1987年开始，因此metyr=2，从第三年开始读取气候资料
		//假如1988年地类为0，则将所有的状态量设置为0
		if(sitec.epctype==0)
		{
			if(ok && make_zero_state_struct(&ws, &cs, &ns))
			{
				bgc_printf(BV_ERROR, "Error in call to make_zero_state_struct(), from bgc()\n");
				ok=0;
			}		
		}
		bgc_printf(BV_DIAG, "done firstday\n");
	}

	/* zero water, carbon, and nitrogen source and sink variables */
	if (ok && zero_srcsnk(&cs,&ns,&ws,&summary))
	{
		bgc_printf(BV_ERROR, "Error in call to zero_srcsnk(), from bgc()\n");
		ok=0;
	}

	bgc_printf(BV_DIAG, "done zero_srcsnk\n");
	
	/* initialize the indicator for first day of current simulation, so
	that the checks for mass balance can have two days for comparison */
	first_balance = 1;

	/* mode == MODE_SPINUP only*/
	if (mode == MODE_SPINUP)
	{
		/* for simulations with fewer than 50 metyears, find the multiple of
		metyears that gets close to 100, use this as the block size in
		spinup control */
		
		if (ctrl.metyears < 50)
		{
			ntimesmet = 100 / ctrl.metyears;
			nblock = ctrl.metyears * ntimesmet;
		}
		else
		{
			nblock = ctrl.metyears;
		}

		/* initialize spinup control variables */
		spinyears = 0;
		metcycle = 0;
		steady1 = 0;
		steady2 = 0;
		rising = 1;
	}

	if (mode == MODE_MODEL)
	{
		tmpyears = ctrl.simyears;
	}
	else if (mode == MODE_SPINUP)
	{
		tmpyears = nblock;
	}

	// read high time resolution data to sfData
	std::vector<StationDataFlux*> sfData;
	if (bgcin->hModel.active)
		if (!readStationFluxData(sfData, bgcin->hModel.stationFile, tmpyears)) return 1;

	// read high temp correction file
	std::vector<float> tempCorrFactor;
	if (bgcin->hModel.tempCorr && bgcin->hModel.active)
		if (!readTempCorrFactor(tempCorrFactor, bgcin->hModel.tempCorrFile)) return 1;

	// read lai data
	std::vector<float> laiData;
	if (bgcin->laiM.active)
		if (!readLaiData(laiData, bgcin->laiM.laiFile, tmpyears)) return 1;

	/* do loop for spinup. will only execute once for MODE_MODEL */
	do
	{
		/* begin the annual model loop */
		for (simyr=0 ; ok && simyr<tmpyears ; simyr++)
		{
			//printf("year= %d, simyear=%d\n",ctrl.simstartyear+simyr,simyr);
			if (mode == MODE_MODEL)
			{
				/* reset the simple annual output variables for text output */
				annmaxlai = 0.0;
				annet = 0.0;
				anntrans = 0.0;
				anneva =0.0;
				anncanopy = 0.0;

				ann_soilw = 0.0;
				ann_snoww = 0.0;
				ann_FSM_soilw = 0.0;
				ann_FSM_icew = 0.0;
				ann_FSM_snoww = 0.0;
				value_c = 0.0;//zlyadd
				yN =0;
				annoutflow = 0.0;
				annnpp = 0.0;
				annnbp = 0.0;
				annprcp = 0.0;
				anntavg = 0.0;
				/*添加年gpp,nep,nee,mr,gr*/
				ann_nep=0.0;
				ann_nee=0.0;
				ann_gpp=0.0;
				ann_mr=0.0;
				ann_gr=0.0;
				ann_hr=0.0;
				ann_root_mr=0.0;
				ann_leaf_mr=0.0;
			}
		
			/* set current month to 0 (january) at the beginning of each year */
			curmonth = 0;

			if (mode == MODE_SPINUP)
			{
				/* calculate scaling for N additions (decreasing with
				time since the beginning of metcycle = 0 block */
				naddfrac = 1.0 - ((double)simyr/(double)nblock);

				if (metcycle == 0)
				{
					tally1 = 0.0;
					tally1b = 0.0;
					tally2 = 0.0;
					tally2b = 0.0;
				}
			}
		
		/* test whether metyr needs to be reset */
		if (metyr == ctrl.metyears)
		{
			if (mode == MODE_MODEL)
			{
				if (ctrl.onscreen) bgc_printf(BV_DETAIL, "Resetting met data for cyclic input\n");
			}
			if (mode == MODE_SPINUP)
			{
				bgc_printf(BV_DIAG, "Resetting met data for cyclic input\n");
			}
			metyr = 0;
		}

		if (mode == MODE_MODEL)
		{
			/* output to screen to indicate start of simulation year */
			//if (ctrl.onscreen) bgc_printf(BV_DETAIL, "Year: %6d\n",ctrl.simstartyear+simyr);
		}
		else if (mode == MODE_SPINUP)
		{
			/* output to screen to indicate start of simulation year */
			//if (ctrl.onscreen) bgc_printf(BV_DETAIL, "Year: %6d\n",spinyears);
		}

		/* set the max lai variable, for annual diagnostic output */
		epv.ytd_maxplai = 0.0;
		
		if (mode == MODE_MODEL)
		{
			/* atmospheric CO2 and Ndep handling */
			if (!(co2.varco2))
			{
				/* constant CO2, constant Ndep */
				metv.co2 = co2.co2ppm;
				daily_ndep = sitec.ndep/365.0;
				daily_nfix = sitec.nfix/365.0;
				//printf("bgc,9 daily_ndep= %lf daily_nfix=%lf\n",daily_ndep,daily_nfix);

			}
			else 
			{
				/* when varco2 = 1, use file for co2 */
				if (co2.varco2 == 1) metv.co2 = get_co2(&co2,(ctrl.simstartyear+simyr));
				//printf("year=%d,metv.co2 = %lf \n",ctrl.simstartyear+simyr,metv.co2);
                bgc_printf(BV_DIAG,"CO2 val: %lf Year: %i\n",metv.co2,(ctrl.simstartyear+simyr));
				if(metv.co2 < -999)
				{
					bgc_printf(BV_ERROR,"Error finding CO2 value for year: %i\n",(ctrl.simstartyear+simyr));
					return(EXIT_FAILURE);
				}

				/* when varco2 = 2, use the constant CO2 value, but vary Ndep */
				if (co2.varco2 == 2) metv.co2 = co2.co2ppm;
				
				if (ramp_ndep.doramp && !bgcin->ndepctrl.varndep)
				{
					/* increasing CO2, ramped Ndep */
					ind_simyr = ramp_ndep.ind_year - ctrl.simstartyear;
			//		printf("ind_simyr= %d\n",ind_simyr);
					ndep_scalar = (ramp_ndep.ind_ndep - ramp_ndep.preind_ndep) / 
						(co2.co2ppm_array[ind_simyr]-co2.co2ppm_array[0]);
			/*		printf("ndep_scalar = %lf\n",ndep_scalar);
					printf("ramp_ndep.ind_ndep - ramp_ndep.preind_ndep = %lf\n",ramp_ndep.ind_ndep - ramp_ndep.preind_ndep);
					printf("ndep_scalar = %lf\n",ndep_scalar);
					printf("co2.co2ppm_array= %lf\n",co2.co2ppm_array[ind_simyr]);
					printf("bgc,10 ndep_scalar = %lf co2=%lf\n",daily_ndep,daily_nfix);
					*/
					ndep_diff = (co2.co2ppm_array[simyr] - co2.co2ppm_array[0]) * 
						ndep_scalar;
					ndep = ramp_ndep.preind_ndep + ndep_diff;
			//		printf("ramp_ndep.preind_ndep= %lf, ndep_diff = %lf\n",ramp_ndep.preind_ndep,ndep_diff);

					/* don't allow the industrial ndep levels to be less than
					the preindustrial levels */
					if (ndep < ramp_ndep.preind_ndep) ndep = ramp_ndep.preind_ndep; 
					daily_ndep = ndep/365.0;
					daily_nfix = sitec.nfix/365.0;
					//printf("bgc_ramp,10 daily_ndep= %lf daily_nfix=%lf\n",daily_ndep,daily_nfix);
	
				}
				else
				{
					/* increasing CO2, constant Ndep */
					daily_ndep = sitec.ndep/365.0;
					daily_nfix = sitec.nfix/365.0;	
					

				}
			}
			if(bgcin->ndepctrl.varndep && mode == MODE_MODEL)
			{
				//printf("bgc,12 daily_ndep= %lf daily_nfix=%lf\n",daily_ndep,daily_nfix);

				daily_ndep = get_ndep(&bgcin->ndepctrl,(ctrl.simstartyear + simyr));
				if(daily_ndep < -999)
				{
					bgc_printf(BV_ERROR, "Error finding NDEP for year: %i\n",(ctrl.simstartyear+simyr));
					return(EXIT_FAILURE);
				}
				else
				{
					bgc_printf(BV_DIAG, "Using annual NDEP value: %lf\n",daily_ndep);				
					daily_ndep /= 365.0;
				}
				//printf("bgc_file,11 daily_ndep= %lf daily_nfix=%lf\n",daily_ndep,daily_nfix);

			}
		}
				//printf("bgc, daily_ndep= %lf daily_nfix=%lf\n",daily_ndep,daily_nfix);
		else if (mode == MODE_SPINUP)
		{
			/* atmospheric concentration of CO2 (ppm) */
			/* Always assign a fixed CO2 value for spinups */
			metv.co2 = co2.co2ppm;

			/*if (!(co2.varco2)) 
			else metv.co2 = co2.co2ppm_array[simyr]; */
		}
		/* begin the daily model loop */
		for (yday=0 ; ok && yday<365 ; yday++)
		{
			bgc_printf(BV_DIAG, "year %d\tyday %d\n",simyr,yday);
			tsoil_sum = 0;	/* add 2014-07-30 */
			epv_psi_sum = 0;	/* add 2014-07-30 */
			ws.soilw_add = 0;
			//农田修改，在每年第一天重新初始化transfer
			if(epc.epctype == 141|| epc.epctype == 142 || epc.epctype == 143)
			{
				cs.leafc_transfer = 0.01;///cinit.max_leafc;
				//cs.frootc_transfer = 0.01;
				cs.livecrootc_transfer = 0;
				cs.deadcrootc_transfer = 0;
				if(yday ==0)
					 phenology_lastday=0;
			}
			/* Test for very low state variable values and force them
			to 0.0 to avoid rounding and floating point overflow errors */
			if (ok && precision_control(&ws, &cs, &ns))
			{
				bgc_printf(BV_ERROR, "Error in call to precision_control() from bgc()\n");
				ok=0;
			}
			
			/* set the day index for meteorological and phenological arrays */
			metday = metyr*365 + yday;

			/* zero all the daily flux variables */
			wf = zero_wf;
			cf = zero_cf;
			nf = zero_nf;

			/* daily meteorological variables from metarrays */
			if (ok && daymet(&metarr, &metv, metday))
			{
				bgc_printf(BV_ERROR, "Error in daymet() from bgc()\n");
				ok=0;
			}
		
			bgc_printf(BV_DIAG, "%d\t%d\tdone daymet\n",simyr,yday);
	
			/* soil temperature correction using difference from annual average tair */
			tdiff = tair_avg - metv.tsoil;
			if (ws.snoww)
				metv.tsoil += 0.83 * tdiff;
			else
				metv.tsoil += 0.2 * tdiff;
			
			/* daily phenological variables from phenarrays */
			if (ok && dayphen(&phenarr, &phen, metday))
			{
				bgc_printf(BV_ERROR, "Error in dayphen() from bgc()\n");
				ok=0;
			}
		//	if(simyr*365+yday == typec_day ||simyr*365+yday == typec_day +1)
		//		printf("物候期%lf,%lf,%lf,%lf,%lf\n",phen.predays_litfall,phen.predays_transfer,phen.remdays_curgrowth,phen.remdays_litfall,phen.remdays_transfer);
			bgc_printf(BV_DIAG, "%d\t%d\tdone dayphen\n",simyr,yday);
	
			/* test for the annual allocation day */
			if (phen.remdays_litfall == 1) annual_alloc = 1;
			else annual_alloc = 0;
			
			if(phenology_lastday != (-999))
			{
				phenology_lastday = phen.remdays_litfall;
				//printf("year=%lf,day=%lf,凋落=%lf\n",simyr,yday,phenology_lastday);
			}

			if(epc.epctype == 141|| epc.epctype == 142 || epc.epctype == 143)
			{				
				if(phenology_lastday > 0)
				{
					phenology_lastday = -999;
					summary.leafc_max = cs.leafc + cs.leafc_storage + cs.leafc_transfer;					
					summary.frootc_max = cs.frootc + cs.frootc_storage + cs.frootc_transfer;			
					//计算最大生物量
					summary.crootc_max = cs.livecrootc + cs.livecrootc_storage + cs.livecrootc_transfer +
																cs.deadcrootc + cs.deadcrootc_storage + cs.deadcrootc_transfer;

					summary.vegc_max = cs.leafc + cs.leafc_storage + cs.leafc_transfer + 
																cs.frootc + cs.frootc_storage + cs.frootc_transfer +
																cs.livestemc + cs.livestemc_storage + cs.livestemc_transfer +
																cs.deadstemc + cs.deadstemc_storage + cs.deadstemc_transfer +
																cs.livecrootc + cs.livecrootc_storage + cs.livecrootc_transfer +
																cs.deadcrootc + cs.deadcrootc_storage + cs.deadcrootc_transfer +
																cs.gresp_storage + cs.gresp_transfer + cs.cpool;
					summary.litrc_max = cs.cwdc + cs.litr1c + cs.litr2c + cs.litr3c + cs.litr4c;
					summary.soilc_max = cs.soil1c + cs.soil2c + cs.soil3c + cs.soil4c;
					summary.totalc_max = summary.vegc_max + summary.litrc_max + summary.soilc_max;
					//printf("year=%d,day=%d,flag=%d,叶子= %lf,细根= %lf,其他= %lf\n",simyr,yday,phenology_lastday,summary.leafc_max,summary.frootc_max,summary.crootc_max);
				}
			}


		//	printf("year=%d,yday = %d, 凋落时间=%lf\n",simyr,yday,phen.remdays_litfall);
			/* phenology fluxes */
			if (ok && phenology(&epc, &phen, &epv, &cs, &cf, &ns, &nf))
			{
				bgc_printf(BV_ERROR, "Error in phenology() from bgc()\n");
				ok=0;
			}
			//printf("year=%d,yday = %d, 凋落时间=%lf,叶子=%lf,stor=%lf,trans=%lf\n"
			//	,simyr,yday,phen.remdays_litfall,cs.leafc,cs.leafc_storage,cs.leafc_transfer);
			bgc_printf(BV_DIAG, "%d\t%d\tdone phenology\n",simyr,yday);
			/* calculate leaf area index, sun and shade fractions, and specific
			leaf area for sun and shade canopy fractions, then calculate
			canopy radiation interception and transmission */
			//std::cout << cs.leafc << yday << std::endl;
			if (ok && radtrans(&cs, &epc, &metv, &epv, sitec.sw_alb, laiData, metday))
			{
				bgc_printf(BV_ERROR, "Error in radtrans() from bgc()\n");
				ok=0;
			}

			/* update the ann max LAI for annual diagnostic output */
			if (epv.proj_lai > epv.ytd_maxplai) epv.ytd_maxplai = epv.proj_lai;
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone radtrans\n",simyr,yday);
			
			/* precip routing (when there is precip) */
			if (ok && metv.prcp && prcp_route(&metv, epc.int_coef, epv.all_lai,&wf))
			{
				bgc_printf(BV_ERROR, "Error in prcp_route() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone prcp_route\n",simyr,yday);

			/* snowmelt (when there is a snowpack) */
			if (ok && ws.snoww && snowmelt(&metv, &wf, ws.snoww))
			{
				bgc_printf(BV_ERROR, "Error in snowmelt() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone snowmelt\n",simyr,yday);

			/* bare-soil evaporation (when there is no snowpack) */
			if (ok && !ws.snoww && baresoil_evap(&metv, &wf, &epv.dsr))
			{
				bgc_printf(BV_ERROR, "Error in baresoil_evap() from bgc()\n");
				ok=0;
			}

			bgc_printf(BV_DIAG, "%d\t%d\tdone bare_soil evap\n",simyr,yday);
			/*2014-11-28 add 湿地，修改土壤水分*/
			//农田修改土壤水分，施肥，14为水浇地，15为旱地，16为耕地按照15处理
			if(bgcin->sitec.epctype == 13 || bgcin->sitec.epctype ==141 || bgcin->sitec.epctype ==142 || bgcin->sitec.epctype ==143)
			{
				if(bgcin->sitec.epctype ==13)
				{
					if((ws.soilw / (1000.0 * sitec.soil_depth)) < sitec.vwc_sat * 0.75)
					{
						//vwc = soilw / (1000.0 * sitec->soil_depth);
						//epv.vwc = sitec.vwc_sat * 0.70;	//设定最小水分含量为饱和水分的70%	
						ws.soilw_add = sitec.vwc_sat * 0.75 * 1000 * sitec.soil_depth - ws.soilw; //需要添加的水分 KgH20
						ws.soilw = sitec.vwc_sat * 0.75 * 1000 * sitec.soil_depth;
					}
				}
				else
				{
					//施肥量150kg/hm2，分七天施肥尿素N含量46%，复合肥15%，按照50kgN/hm2施肥【N：P：K=1:0.5:0.5 即是氮肥0.0075kg/m2，可能有误】
					if(yday== 100)// && yday<=104)
					{
						//daily_ndep += 0.0010; //50kg
						if(bgcin->sitec.epctype ==141 || bgcin->sitec.epctype ==142)	
							daily_ndep += 0.0050;  //25kg,
						if(bgcin->sitec.epctype ==143)	
							daily_ndep += 0.005;  //25kg,
					}
					if(bgcin->sitec.epctype ==141)  //水浇地修改土壤水分
					{
						if((ws.soilw / (1000.0 * sitec.soil_depth)) < sitec.vwc_fc * 0.75)
						{
							//vwc = soilw / (1000.0 * sitec->soil_depth);
							//epv.vwc = sitec.vwc_sat * 0.70;	//设定最小水分含量为饱和水分的70%	
							ws.soilw_add = sitec.vwc_fc * 0.75 * 1000 * sitec.soil_depth - ws.soilw; //需要添加的水分 KgH20
							ws.soilw = sitec.vwc_fc * 0.75 * 1000 * sitec.soil_depth;
							//printf("ws.soilw_add=%lf\n",ws.soilw_add);
						}
					}
				}
			}
			
			/* soil water potential */
			if (ok && soilpsi(&sitec, ws.soilw, &epv.psi, &epv.vwc))
			{
				bgc_printf(BV_ERROR, "Error in soilpsi() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone soilpsi\n",simyr,yday);

			/* daily maintenance respiration */
			if (ok && maint_resp(&ws,&sitec,&cs, &ns, &epc, &metv, &cf, &epv, bgcin->hModel, sfData, simyr, yday, mode))
			{
				bgc_printf(BV_ERROR, "Error in m_resp() from bgc()\n");
				ok=0;
			}

			bgc_printf(BV_DIAG, "%d\t%d\tdone maint resp\n",simyr,yday);

			/* begin canopy bio-physical process simulation */
			/* do canopy ET calculations whenever there is leaf area
			displayed, since there may be intercepted water on the 
			canopy that needs to be dealt with */
			wflux_struct wfT = wf;
			if (ok && cs.leafc && metv.dayl)
			{
				/* conductance and evapo-transpiration */
				if (ok && canopy_et(&metv, &epc, &epv, &wf, 0))
				{
					bgc_printf(BV_ERROR, "Error in canopy_et() from bgc()\n");
					ok=0;
				}
				
				bgc_printf(BV_DIAG, "%d\t%d\tdone canopy_et\n",simyr,yday);

			}
			/* do photosynthesis only when it is part of the current
			growth season, as defined by the remdays_curgrowth flag.  This
			keeps the occurrence of new growth consistent with the treatment
			of litterfall and allocation */

			// add for high time resolution, 30 minutes
			epvar_struct epvT = epv;
			cflux_struct cfT = cf;
			//end

			if (ok && cs.leafc && phen.remdays_curgrowth && metv.dayl)
			{
				if (ok && total_photosynthesis(&metv, &epc, &epv, &cf, &psn_sun, &psn_shade))//,fp_gs,mode))
				{
					bgc_printf(BV_ERROR, "Error in total_photosynthesis() from bgc()\n");
					ok=0;
				}
				
			} /* end of photosynthesis calculations */
			else
			{
				epv.assim_sun = epv.assim_shade = 0.0;
			}

			// add for high time resolution, 30 minutes
			psn_struct psn_sunT = psn_sun;
			psn_struct psn_shadeT = psn_shade;
			// add for high time resolution
			if (bgcin->hModel.active && mode == MODE_MODEL && cs.leafc && phen.remdays_curgrowth && metv.dayl)
			{
				total_photosynthesisTimeRes(tempCorrFactor, &bgcin->hModel, &wfT, sfData, &cs, sitec.sw_alb, &metv,
					&epc, &epvT, &cfT, &psn_sunT, &psn_shadeT, simyr, yday);
			}

			if(cs.leafc && phen.remdays_curgrowth && metv.dayl && mode == MODE_MODEL)
				replacePhotosynthesisResults(&bgcin->hModel, &epv, &cf, &psn_sun, &psn_shade, 
				&epvT, &cfT, &psn_sunT, &psn_shadeT, simyr, yday);
			//end

			if (mode == MODE_MODEL)
			{
				/* nitrogen deposition and fixation */
				nf.ndep_to_sminn = daily_ndep;
				nf.nfix_to_sminn = daily_nfix;
				//printf("开始调用冻土模块 ok=%d \n",ok);

				/*if (ok && soilwh(&sitec,svar,soil_sum,&metv,&wf,fp))   //add
				{
					bgc_printf(BV_ERROR,"Error in soilwh.c \n");
					ok=0;
				}*/
				
				if(bgcin->cinit.frost)
				{
					/*if (ok && FSM_core(&metv,metday,fp1,&pmetvar,&bgcin->vspar,soil_sum,&epv,&sitec,&ws))   //add
					{
						bgc_printf(BV_ERROR,"Error in soilwh.c \n");
						ok=0;
					}*/
					if (ok && FSM_core(&metv,metday,&pmetvar,&bgcin->vspar,soil_sum,&epv,&sitec,&ws))   //add
					{
						bgc_printf(BV_ERROR,"Error in soilwh.c \n");
						ok=0;
					}
				}
				//printf("冻土模块运行结束\n");
			}
			else if (mode == MODE_SPINUP)
			{
				/* nitrogen deposition and fixation */
				nf.ndep_to_sminn = sitec.ndep/365.0;
				nf.nfix_to_sminn = sitec.nfix/365.0;
			}

			/* calculate outflow */
			if (ok && outflow(&sitec, &ws, &wf))
			{
				bgc_printf(BV_ERROR, "Error in outflow() from bgc.c\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone outflow\n",simyr,yday);

			/* daily litter and soil decomp and nitrogen fluxes */
			if(mode == MODE_MODEL)
			{
				//判断是否运行冻土模块
				if(bgcin->cinit.frost)
				{
					if(ok && decomp_soil(soil_sum,metv.tsoil,&epc,&epv,&sitec,&cs,&cf,&ns,&nf,&nt))  //add
					{
						bgc_printf(BV_ERROR,"Error in decomp_soil() from bgc.c\n");
						ok=0;
					}
					bgcin->cs.soilT_error_mark = cs.soilT_error_mark;
					for(i = 0;i < 4;i++)                                                     //add
					{ 
						tsoil_sum += soil_sum[i].s_t;
						epv_psi_sum += soil_sum[i].s_h;
					}
					tsoil = tsoil_sum/3;
					epv_psi = epv_psi_sum/3/102;   //unit transter m->MPa				
				}
				else
				{
					if (ok && (ok && decomp(&ws,metv.tsoil,&epc,&epv,&sitec,&cs,&cf,&ns,&nf,&nt)))					             
					{
						bgc_printf(BV_ERROR, "Error in decomp() from bgc.c\n");
						ok=0;
					}
				}
		
			}
			else if(mode == MODE_SPINUP)
			{
			    if (ok && decomp(&ws,metv.tsoil,&epc,&epv,&sitec,&cs,&cf,&ns,&nf,&nt))
			    {
				    bgc_printf(BV_ERROR, "Error in decomp() from bgc.c\n");
				    ok=0;
			    }
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone decomp\n",simyr,yday);

			/* Daily allocation gets called whether or not this is a
			current growth day, because the competition between decomp
			immobilization fluxes and plant growth N demand is resolved
			here.  On days with no growth, no allocation occurs, but
			immobilization fluxes are updated normally */
			if (mode == MODE_MODEL)
			{
				if (ok && daily_allocation(&cf,&cs,&nf,&ns,&epc,&epv,&nt,1.0,MODE_MODEL))
				{
					bgc_printf(BV_ERROR, "Error in daily_allocation() from bgc.c\n");
					ok=0;
				}
			}
			else if (mode == MODE_SPINUP)
			{
				/* spinup control */
				/* in the rising limb, use the spinup allocation code
				that supplements N supply */
				if (!steady1 && rising && metcycle == 0)
				{
					if (ok && daily_allocation(&cf,&cs,&nf,&ns,&epc,&epv,&nt,naddfrac,MODE_SPINUP))
					{
						bgc_printf(BV_ERROR, "Error in daily_allocation() from bgc.c\n");
						ok=0;
					}
				}
				else
				{
					if (ok && daily_allocation(&cf,&cs,&nf,&ns,&epc,&epv,&nt,1.0,MODE_MODEL))
					{
						bgc_printf(BV_ERROR, "Error in daily_allocation() from bgc.c\n");
						ok=0;
					}
				}
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone daily_allocation\n",simyr,yday);

			/* reassess the annual turnover rates for livewood --> deadwood,
			and for evergreen leaf and fine root litterfall. This happens
			once each year, on the annual_alloc day (the last litterfall day) */
			if (ok && annual_alloc)
			{
				if (ok && annual_rates(&epc,&epv))
				{
					bgc_printf(BV_ERROR, "Error in annual_rates() from bgc()\n");
					ok=0;
				}
				
				bgc_printf(BV_DIAG, "%d\t%d\tdone annual rates\n",simyr,yday);
			} 


			/* daily growth respiration */
			if (ok && growth_resp(&epc, &cf))
			{
				bgc_printf(BV_ERROR, "Error in daily_growth_resp() from bgc.c\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone growth_resp\n",simyr,yday);

			/* daily update of the water state variables */
			if (ok && daily_water_state_update(&wf, &ws))
			{
				bgc_printf(BV_ERROR, "Error in daily_water_state_update() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone water state update\n",simyr,yday);

			/* daily update of carbon state variables */
			if (ok && daily_carbon_state_update(&cf, &cs, annual_alloc,
				epc.woody, epc.evergreen))
			{
				bgc_printf(BV_ERROR, "Error in daily_carbon_state_update() from bgc()\n");
				ok=0;
			}
		
			bgc_printf(BV_DIAG, "%d\t%d\tdone carbon state update\n",simyr,yday);

			/* daily update of nitrogen state variables */
			if (ok && daily_nitrogen_state_update(&nf, &ns, annual_alloc,
				epc.woody, epc.evergreen))
			{
				bgc_printf(BV_ERROR, "Error in daily_nitrogen_state_update() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone nitrogen state update\n",simyr,yday);

			/* calculate N leaching loss.  This is a special state variable
			update routine, done after the other fluxes and states are
			reconciled in order to avoid negative sminn under heavy leaching
			potential */
			if (ok && nleaching(&ns, &nf, &ws, &wf))
			{
				bgc_printf(BV_ERROR, "Error in nleaching() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone nitrogen leaching\n",simyr,yday);

			/* calculate daily mortality fluxes and update state variables */
			/* this is done last, with a special state update procedure, to
			insure that pools don't go negative due to mortality fluxes
			conflicting with other proportional fluxes */
			if (ok && mortality(&epc,&cs,&cf,&ns,&nf))
			{
				bgc_printf(BV_ERROR, "Error in mortality() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone mortality\n",simyr,yday);

			/* test for water balance */
			if (ok && check_water_balance(&ws, first_balance))
			{
				bgc_printf(BV_ERROR, "Error in check_water_balance() from bgc()\n");
				bgc_printf(BV_ERROR, "%d\n",metday);
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone water balance\n",simyr,yday);

			/* test for carbon balance */
			if (ok && check_carbon_balance(&cs, first_balance))
			{
				bgc_printf(BV_ERROR, "Error in check_carbon_balance() from bgc()\n");
				bgc_printf(BV_ERROR, "%d\n",metday);
				ok=0;
			}
	
			bgc_printf(BV_DIAG, "%d\t%d\tdone carbon balance\n",simyr,yday);

			/* test for nitrogen balance */
			if (ok && check_nitrogen_balance(&ns, first_balance))
			{
				bgc_printf(BV_ERROR, "Error in check_nitrogen_balance() from bgc()\n");
				bgc_printf(BV_ERROR, "%d\n",metday);
				ok=0;
			}
		
			bgc_printf(BV_DIAG, "%d\t%d\tdone nitrogen balance\n",simyr,yday);

			/* calculate carbon summary variables */
			if (ok && csummary(&cf, &cs, &summary))
			{
				bgc_printf(BV_ERROR, "Error in csummary() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone carbon summary\n",simyr,yday);

			/* calculate water summary variables */
			if (ok && wsummary(&ws,&wf,&summary))
			{
				printf("Error in wsummary() from bgc()\n");
				ok=0;
			}
			
			bgc_printf(BV_DIAG, "%d\t%d\tdone water summary\n", simyr,yday);

			if (yday == 1)
			{
				/*输出每年土壤碳库*/
					//printf("土壤碳库%lf\n",cs.soil4c+cs.soil3c+cs.soil2c+cs.soil1c);
			}
			/* DAILY OUTPUT HANDLING */
			/* fill the daily output array if daily output is requested,
			or if the monthly or annual average of daily output variables
			have been requested */
			bgc_printf(BV_DIAG, "Number of daily outputs: %d\n", ctrl.ndayout);
			if (ok && dayout)
			{
				/* fill the daily output array */
				for (outv=0 ; outv<ctrl.ndayout ; outv++)
				{
					bgc_printf(BV_DIAG, "Outv: %d, ", outv);
					bgc_printf(BV_DIAG, "DayCode: %d, ", ctrl.daycodes[outv]);
					bgc_printf(BV_DIAG, "Output: %f\n", *output_map[ctrl.daycodes[outv]]);
					dayarr[outv] = (float) *output_map[ctrl.daycodes[outv]];
				}
			}
			/* only write daily outputs if requested */
			if (ok && ctrl.dodaily)
			{
				/* write the daily output array to daily output file 
				if (fwrite(dayarr, sizeof(float), ctrl.ndayout, bgcout->dayout.ptr)
					!= (size_t)ctrl.ndayout)
				{
					bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
						bgcout->dayout.name,simyr,yday);
					ok=0;
				}
				*/
				bgc_printf(BV_DIAG, "%d\t%d\tdone daily output\n",simyr,yday);

				/**/	if(ok && bgcout->bgc_ascii)
				{				
					output_ascii(dayarr,ctrl.ndayout,bgcout->dayoutascii.ptr);
					//printf("write dayfile\n");
					//output_soil(soil_sum,bgcout->dayout_soil.ptr);					
				}			
			}
			/*******************/
			/* MONTHLY OUTPUTS */
			/*******************/
			/* MONTHLY AVERAGE OF DAILY OUTPUT VARIABLES */
			if (ctrl.domonavg)
			{
				/* update the monthly average array */
				for (outv=0 ; outv<ctrl.ndayout ; outv++)
				{
					monavgarr[outv] += dayarr[outv];

					switch (ctrl.daycodes[outv])
					{
						/* Leaf area index */
						case 545:   
							if(dayarr[outv] > monmaxlai) monmaxlai = dayarr[outv]; 
							break;
					}
				}
				
				/* if this is the last day of the current month, output... */

				if (yday == endday[curmonth])
				{
					/* finish the averages */
					for (outv=0 ; outv<ctrl.ndayout ; outv++)
					{
						if (summary_sanity == SANE)
						{
							switch (ctrl.daycodes[outv])
							{
								/* Leaf area index */
								/* Maximum monthly */
								case 545:
									monavgarr[outv] = monmaxlai; 
									break;
								/* Snow water */
								case 21:
									monavgarr[outv] = dayarr[outv] - eomsnoww; 
									eomsnoww = dayarr[outv]; 
									break;
								/* Soil water content */
								case 20:
									monavgarr[outv] = dayarr[outv] - eomsoilw;
									eomsoilw = dayarr[outv];
									//printf("eomsoilw=%lf\n",eomsoilw);
									break;
								default:
									monavgarr[outv] /= (float)mondays[curmonth];
									break;
							}
						}
						else 
						{
							monavgarr[outv] /= (float)mondays[curmonth];
						}
					}

					/* write to file 
					if (fwrite(monavgarr, sizeof(float), ctrl.ndayout, bgcout->monavgout.ptr)
						!= (size_t)ctrl.ndayout)
					{
						bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
							bgcout->monavgout.name,simyr,yday);
						ok=0;
					}
					//2013-08-27*/
					if(ok && bgcout->bgc_ascii)
					{
						output_ascii(monavgarr,ctrl.ndayout, bgcout->monoutascii.ptr);
					//	printf("write monthfile\n");
					}
					
					/* reset monthly average variables for next month */
					for (outv=0 ; outv<ctrl.ndayout ; outv++)
					{
						monavgarr[outv] = 0.0;
						monmaxlai = 0.0;
						monmaxsnoww = 0.0;
					}
					
					/* increment current month counter */
					curmonth++;
					
					bgc_printf(BV_DIAG, "%d\t%d\tdone monavg output\n",simyr,yday);
				
				}
			}
			
			/* ANNUAL AVERAGE OF DAILY OUTPUT VARIABLES */
			if (ctrl.doannavg)
			{
				/* update the annual average array */
				for (outv=0 ; outv<ctrl.ndayout ; outv++)
				{
					annavgarr[outv] += dayarr[outv];
					switch (ctrl.daycodes[outv])
					{
						/* Leaf area index */
						case 545:
							if(dayarr[outv] > annmaxplai) annmaxplai = dayarr[outv];
							break;
					}
				}
				
				/* if this is the last day of the year, output... */
				if (yday == 364)
				{
	
					/* finish averages */
					for (outv=0 ; outv<ctrl.ndayout ; outv++)
					{
						if (summary_sanity == SANE)
						{
							switch (ctrl.daycodes[outv])
							{
								/* Leaf area index*/ 
								case 545:
									annavgarr[outv] = (float)annmaxplai;
									break;
								default: 
									annavgarr[outv] /= 365.0;
									break;
							}
						}
						else
						{
							annavgarr[outv] /= 365.0;
						}
					}
					
					/* write to file
					if (fwrite(annavgarr, sizeof(float), ctrl.ndayout, bgcout->annavgout.ptr)
						!= (size_t)ctrl.ndayout)
					{
						bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
							bgcout->annavgout.name,simyr,yday);
						ok=0;
					} */
					
					/* reset annual average variables for next month */
					for (outv=0 ; outv<ctrl.ndayout ; outv++)
					{
						annavgarr[outv] = 0.0;
						annmaxplai = 0.0;
					}
					
					bgc_printf(BV_DIAG, "%d\t%d\tdone annavg output\n",simyr,yday);
				
					
					
				}
				
			}
			
			if (mode == MODE_MODEL)
			{
				/* very simple annual summary variables for text file output */
				if (epv.proj_lai > (double)annmaxlai) annmaxlai = (float)epv.proj_lai;
				annet += wf.canopyw_evap + wf.snoww_subl + wf.soilw_evap +
					wf.soilw_trans;
				//蒸发量
				anneva += wf.canopyw_evap + wf.snoww_subl + wf.soilw_evap;
				//蒸腾量
				anntrans += wf.soilw_trans;
				anncanopy += wf.prcp_to_canopyw;
				ann_soilw += ws.soilw;
				ann_snoww += ws.snoww;
				ann_FSM_soilw += ws.FSM_soilw;
				ann_FSM_icew += ws.FSM_icew;
				ann_FSM_snoww += ws.FSM_snoww;
				yN ++;

				annoutflow += wf.soilw_outflow;
				annnpp += summary.daily_npp * 1000.0;
				annnbp += summary.daily_nee * 1000.0;
				annprcp += metv.prcp;
				anntavg += metv.tavg/365.0;
				/*添加年gpp,nep,nee,mr,gr,root_mr*/
				ann_nep += summary.daily_nep*1000.0;
				ann_nee += summary.daily_nee*1000.0;
				ann_gpp += summary.daily_gpp*1000.0;
				ann_mr += summary.daily_mr*1000.0;
				ann_gr += summary.daily_gr*1000.0;
				ann_hr += summary.daily_hr*1000.0;
				ann_root_mr += summary.daily_root_mr*1000.0;
				ann_leaf_mr += summary.daily_leaf_mr*1000.0;
			}
			else if (mode == MODE_SPINUP)
			{
				/* spinup control */
				/* keep a tally of total soil C during successive
				met cycles for comparison */
				if (metcycle == 1)
				{
					tally1 += summary.soilc;
					tally1b += summary.totalc;
				}
				if (metcycle == 2)
				{
					tally2 += summary.soilc;
					tally2b += summary.totalc;
				}
			}
			
			/* at the end of first day of simulation, turn off the 
			first_balance switch */
			if (first_balance) first_balance = 0;

			/*从某一天开始进行土地利用类型的转换，将模型从新进行初始化，调整碳氮库。
			修改转换之后的类型*/
		//	if(simyr*365+yday==typec_day+1)
		//		printf("转化之后第%d年%d天的stemc= %lf\n",simyr,yday,summary.stemc_sum);
			if (mode == MODE_MODEL)
			{
				if(cinit.change_date == (double)(simyr*365+yday))
				{
					//make_zero_flux_struct();
					//考虑从林地1转为草地3，草地转为林地，灌木2转为草地，草地转为灌木
					//printf("地类转换前后，%lf,%lf\n",sitec.epctype,sitec.epctype_change);
					
					//&& (sitec.epctype>15 || sitec.epctype==0||sitec.epctype_change>15||sitec.epctype_change==0))
					if((sitec.epctype!=sitec.epctype_change))//判断地类是否转化，不转化则跳过
					{

						//1乔木、灌木、竹林、建筑用地转化为草地，4种
						if(sitec.epctype<=17 && sitec.epctype_change==18)
						{
							//printf("1林地转化为草地\n");
							//1:重新读取epc函数
							_itoa(sitec.epctype_change,epctype_char,10);
							strcat(sitec.epcfilepath,epctype_char);
							strcat(sitec.epcfilepath,".epc");
							fopen_s(&(epc_file.ptr),sitec.epcfilepath,"r");		
							//printf("干与叶子的分配alloc_newstemc_newleafc=%lf\n",epc.alloc_newstemc_newleafc);
							//读取生理生态参数 *.epc
							if(epc_init(epc_file, &epc, &sitec))
							{
								bgc_printf(BV_ERROR, "Error in call to epc_init() from pointbgc.c... Exiting\n");
								exit(EXIT_FAILURE);
							}
							fclose(epc_file.ptr);
							//printf("干与叶子的分配alloc_newstemc_newleafc=%lf\n",epc.alloc_newstemc_newleafc);

							//2:调整各部分的碳库；
							/* make zero-flux structures for use inside annual and daily loops*/
							if (ok && make_zero_flux_struct(&zero_wf, &zero_cf, &zero_nf))
							{
								bgc_printf(BV_ERROR, "Error in call to make_zero_flux_struct()1 from bgc()\n");
								ok=0;
							}
							bgc_printf(BV_DIAG, "done make_zero_flux\n");
										/* zero all the daily flux variables 
							wf = zero_wf;
							cf = zero_cf;
							nf = zero_nf;*/

							//printf("deadstemc=%lf\n",cs.deadstemc);
							
							//林下草的生物量作为森林转化为草地之后的生物量。该比例需要查阅文献。
							cs.leafc=0.1;//summary.leafc_sum*0.02;
							cs.frootc=cs.leafc;
							/*如果storage不为0，则存储的C按照已经经过的天数占全年的比例分配*/
							if(cs.leafc_storage)
								cs.leafc_storage=0;//summary.vegc*0.02*yday/365;
							if(cs.leafc_transfer)
								cs.leafc_transfer=0;
							if(cs.frootc_storage)
								cs.frootc_storage=0;
							if(cs.frootc_transfer)
								cs.frootc_transfer=0;

							/**/
							if(re_init(&sitec, &epc, &cs, &cinit,&ns))
							{
								bgc_printf(BV_ERROR, "Error in call to sitec_init2() from pointbgc.c... Exiting\n");
								exit(EXIT_FAILURE);
							}

							/*3 重新计算物候 2013-07-02*/
							if (ok && free_phenmem(&phenarr))
							{
								bgc_printf(BV_ERROR, "Error in free_phenmem() from bgc()\n");
								ok=0;
							}

							if (ok && prephenology(bgcin->gsiM, &ctrl, &epc, &sitec, &metarr, &phenarr))
							{
								bgc_printf(BV_ERROR, "Error in call to prephenology(), from bgc()\n");
								ok=0;
							}
						}
						/*2 草地转化为乔木、灌木或竹林，3种*/
						if(sitec.epctype ==18 && sitec.epctype_change<=17 && sitec.epctype_change>0)
						{
							//printf("2草地转化为林地\n");
							//1:重新读取epc函数，读取生理生态参数 *.epc
							_itoa(sitec.epctype_change,epctype_char,10);
							strcat(sitec.epcfilepath,epctype_char);
							strcat(sitec.epcfilepath,".epc");
							fopen_s(&(epc_file.ptr),sitec.epcfilepath,"r");
							if(epc_init(epc_file, &epc, &sitec))
							{
								bgc_printf(BV_ERROR, "Error in call to epc_init() from pointbgc.c... Exiting\n");
								exit(EXIT_FAILURE);
							}
							fclose(epc_file.ptr);

							//2:调整各部分的碳库；
							/* make zero-flux structures for use inside annual and daily loops*/
							if (ok && make_zero_flux_struct(&zero_wf, &zero_cf, &zero_nf))
							{
								bgc_printf(BV_ERROR, "Error in call to make_zero_flux_struct()2 from bgc()\n");
								ok=0;
							}
							bgc_printf(BV_DIAG, "done make_zero_flux\n");	

							/*凋落物和土壤保持不变，对最大叶子和最大树干初始化*/
							if(re_init(&sitec, &epc, &cs, &cinit,&ns))
							{
								bgc_printf(BV_ERROR, "Error in call to sitec_init2() from pointbgc.c... Exiting\n");
								exit(EXIT_FAILURE);
							}
							/*3 重新计算物候 2013-07-02*/
							if (ok && free_phenmem(&phenarr))
							{
								bgc_printf(BV_ERROR, "Error in free_phenmem() from bgc()\n");
								ok=0;
							}

							if (ok && prephenology(bgcin->gsiM, &ctrl, &epc, &sitec, &metarr, &phenarr))
							{
								bgc_printf(BV_ERROR, "Error in call to prephenology(), from bgc()\n");
								ok=0;
							}
							/*4 根据物候将max_leafC和max_stemC分配到各个部分*/
							if (ok && firstday(&epc, &cinit, &epv, &phenarr, &cs, &ns))
							{
								bgc_printf(BV_ERROR, "Error in call to firstday(), from bgc()\n");
								ok=0;
							}

						}
						/*3 乔木或竹林转化为灌木，乔木转化为竹林，建筑用地转灌木或竹林，5种 2013-07-26*/
						if((sitec.epctype<=15 || sitec.epctype==17) && (sitec.epctype_change==16 ||sitec.epctype_change==17))
						{
							if(sitec.epctype != sitec.epctype_change)
							{
								//	printf("3乔木或竹林转化为灌木,或乔木转化为竹林\n");
								//1:重新读取epc函数，读取生理生态参数 *.epc
								_itoa(sitec.epctype_change,epctype_char,10);
								strcat(sitec.epcfilepath,epctype_char);
								strcat(sitec.epcfilepath,".epc");
								fopen_s(&(epc_file.ptr),sitec.epcfilepath,"r");
								if(epc_init(epc_file, &epc, &sitec))
								{
									bgc_printf(BV_ERROR, "Error in call to epc_init() from pointbgc.c... Exiting\n");
									exit(EXIT_FAILURE);
								}
								fclose(epc_file.ptr);

								//2:调整各部分的碳库；
								/* make zero-flux structures for use inside annual and daily loops*/
								if (ok && make_zero_flux_struct(&zero_wf, &zero_cf, &zero_nf))
								{
									bgc_printf(BV_ERROR, "Error in call to make_zero_flux_struct()3 from bgc()\n");
									ok=0;
								}
								bgc_printf(BV_DIAG, "done make_zero_flux\n");	

								/*凋落物和土壤保持不变，对最大叶子和最大树干初始化*/
								if(re_init(&sitec, &epc, &cs, &cinit,&ns))
								{
									bgc_printf(BV_ERROR, "Error in call to sitec_init2() from pointbgc.c... Exiting\n");
									exit(EXIT_FAILURE);
								}
								/*3 重新计算物候 2013-07-02*/
								if (ok && free_phenmem(&phenarr))
								{
									bgc_printf(BV_ERROR, "Error in free_phenmem() from bgc()\n");
									ok=0;
								}

								if (ok && prephenology(bgcin->gsiM, &ctrl, &epc, &sitec, &metarr, &phenarr))
								{
									bgc_printf(BV_ERROR, "Error in call to prephenology(), from bgc()\n");
									ok=0;
								}
								/*4 根据物候将max_leafC和max_stemC分配到各个部分*/
								if (ok && firstday(&epc, &cinit, &epv, &phenarr, &cs, &ns))
								{
									bgc_printf(BV_ERROR, "Error in call to firstday(), from bgc()\n");
									ok=0;
								}
							}
						}
						/*4 灌木转化为乔木或竹林,竹林、建筑用地转化为乔木2013-07-26 4种*/
						if((sitec.epctype ==17 || sitec.epctype ==16|| (sitec.epctype ==0 && sitec.epctype_change!=17)) && ((sitec.epctype_change<=15 && sitec.epctype_change>0) || sitec.epctype_change==17))
						{
							if(sitec.epctype != sitec.epctype_change)
							{
							//	printf("4灌木转化为乔木或竹林，竹林转化为乔木，或建筑用地转化为森林\n");
								//1:重新读取epc函数，读取生理生态参数 *.epc
								_itoa(sitec.epctype_change,epctype_char,10);
								strcat(sitec.epcfilepath,epctype_char);
								strcat(sitec.epcfilepath,".epc");
								fopen_s(&(epc_file.ptr),sitec.epcfilepath,"r");
								if(epc_init(epc_file, &epc, &sitec))
								{
									bgc_printf(BV_ERROR, "Error in call to epc_init() from pointbgc.c... Exiting\n");
									exit(EXIT_FAILURE);
								}
								fclose(epc_file.ptr);

								//2:调整各部分的碳库；
								/* make zero-flux structures for use inside annual and daily loops*/
								if (ok && make_zero_flux_struct(&zero_wf, &zero_cf, &zero_nf))
								{
									bgc_printf(BV_ERROR, "Error in call to make_zero_flux_struct()4 from bgc()\n");
									ok=0;
								}
								bgc_printf(BV_DIAG, "done make_zero_flux\n");	

								/*凋落物和土壤保持不变，对最大叶子和最大树干初始化*/
								if(re_init(&sitec, &epc, &cs, &cinit,&ns))
								{
									bgc_printf(BV_ERROR, "Error in call to sitec_init2() from pointbgc.c... Exiting\n");
									exit(EXIT_FAILURE);
								}
								/*3 重新计算物候 2013-07-02*/
									/* free phenology memory */
								if (ok && free_phenmem(&phenarr))
								{
									bgc_printf(BV_ERROR, "Error in free_phenmem() from bgc()\n");
									ok=0;
								}
								if (ok && prephenology(bgcin->gsiM, &ctrl, &epc, &sitec, &metarr, &phenarr))
								{
									bgc_printf(BV_ERROR, "Error in call to prephenology(), from bgc()\n");
									ok=0;
								}
								/*4 根据物候将max_leafC和max_stemC分配到各个部分*/
								if (ok && firstday(&epc, &cinit, &epv, &phenarr, &cs, &ns))
								{
									bgc_printf(BV_ERROR, "Error in call to firstday(), from bgc()\n");
									ok=0;
								}
							}
						}//结束4
						/*5. 森林，灌木，竹林，草地转为建筑用地，4种*/
						if(sitec.epctype_change==0)
						{		
							//所有的量都变为0；
						//	printf("5.森林，灌木，竹林，草地转为建筑用地，4种\n");
							if(ok && make_zero_state_struct(&ws, &cs, &ns))
							{
								bgc_printf(BV_ERROR, "Error in call to make_zero_state_struct(), from bgc()\n");
								ok=0;
							}							
							if (ok && make_zero_flux_struct(&zero_wf, &zero_cf, &zero_nf))
							{
								bgc_printf(BV_ERROR, "Error in call to make_zero_flux_struct()5 from bgc()\n");
								ok=0;
							}
							/* zero water, carbon, and nitrogen source and sink variables */
							if (ok && zero_srcsnk(&cs,&ns,&ws,&summary))
							{
								bgc_printf(BV_ERROR, "Error in call to zero_srcsnk(), from bgc()\n");
								ok=0;
							}

						}

					}//结束判断，不发生类型变化则跳过。
				}//在随机的某一天开始转化
			}

		}   
		/* end of daily model loop */
		


		/* ANNUAL OUTPUT HANDLING */
		/* only write annual outputs if requested */
		if (ok && ctrl.doannual)
		{
			/* fill the annual output array */
			for (outv=0 ; outv<ctrl.nannout ; outv++)
			{
				annarr[outv] = (float) *output_map[ctrl.anncodes[outv]];
			}
			/* write the annual output array to annual output file 
			if (fwrite(annarr, sizeof(float), ctrl.nannout, bgcout->annout.ptr)
				!= (size_t)ctrl.nannout)
			{
				bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
					bgcout->annout.name,simyr,yday);
				ok=0;
			}*/
			
			if(ok && bgcout->bgc_ascii)
			{
				/*if(write_title)
				{
					fprintf(bgcout->annoutascii.ptr,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n","LAI","vegC",
					"plotID","litrC","soilC","totalC","leafC_sum","leafC_max","leafC","stemC","crootC","cwdc","litr1C","litr2C","litr3C",
					"litr4C","soil1C","soil2C","soil3C","soil4C");
					write_title = 0;
				}*/
				fprintf(bgcout->annoutascii.ptr,"%10d  ",bgcin->sitec.plotID);
				output_ascii(annarr,ctrl.nannout,bgcout->annoutascii.ptr);
				
			}
			bgc_printf(BV_DIAG, "%d\t%d\tdone annual output\n",simyr,yday);
		}
		
		/*******************************************************************************/
		if(bgcin->sitec.epctype >= 0)
		{	
			if (mode == MODE_MODEL && bgcout->bgc_ascii)
			{
				ann_soilw = ann_soilw /yN;
				ann_snoww = ann_snoww /yN;
				ann_FSM_soilw = ann_FSM_soilw /yN;
				ann_FSM_icew = ann_FSM_icew /yN;
				ann_FSM_snoww = ann_FSM_snoww /yN;
				value_c= ann_gpp*10;//zlyadd
				/* write the simple annual text output */
				fprintf(bgcout->anntext.ptr,"%10d%6d%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",
					bgcin->sitec.plotID,ctrl.simstartyear+simyr,annprcp,anntavg,annmaxlai,annet,anntrans,anneva,anncanopy,annoutflow,annnpp,annnbp
					,ann_nep,ann_nee,ann_gpp,ann_mr,ann_gr,ann_hr,ann_root_mr,ann_leaf_mr
					,ann_soilw,ann_snoww);
				/*fprintf(bgcout->anntext.ptr,"%10d%6d%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",
					bgcin->sitec.plotID,ctrl.simstartyear+simyr,annprcp,anntavg,annmaxlai,annet,anntrans,anneva,anncanopy,annoutflow,annnpp,annnbp
					,ann_nep,ann_nee,ann_gpp,ann_mr,ann_gr,ann_hr,ann_root_mr,ann_leaf_mr
					,ann_soilw,ann_snoww,ann_FSM_soilw,ann_FSM_icew,ann_FSM_snoww,value_c);*/
			}
		}
		metyr++;

		if (mode == MODE_SPINUP)
		{
			/* spinup control */
			spinyears++;
		}

	}   /* end of annual model loop */

	if (mode == MODE_SPINUP)
	{
		/* spinup control */
		/* if this is the third pass through metcycle, do comparison */
		/* first block is during the rising phase */
		if (!steady1 && metcycle == 2)
		{
			/* convert tally1 and tally2 to average daily soilc */
			tally1 /= (double)nblock * 365.0;
			tally2 /= (double)nblock * 365.0;
			rising = (tally2 > tally1);
			t1 = (tally2-tally1)/(double)nblock;
			steady1 = (fabs(t1) < SPINUP_TOLERANCE);

			//bgc_printf(BV_DIAG, "spinyears = %d rising = %d steady1 = %d\n",spinyears,
			//	rising,steady1);
			//bgc_printf(BV_DIAG, "metcycle = %d tally1 = %lf tally2 = %lf pdif = %lf\n\n",
			//	metcycle,tally1,tally2,t1);
			if (steady1) bgc_printf(BV_DIAG, "SWITCH\n\n");

			metcycle = 0;
		}
		/* second block is after supplemental N turned off */
		else if (steady1 && metcycle == 2)
		{
			/* convert tally1 and tally2 to average daily soilc */
			tally1 /= (double)nblock * 365.0;
			tally2 /= (double)nblock * 365.0;
			t1 = (tally2-tally1)/(double)nblock;
			steady2 = (fabs(t1) < SPINUP_TOLERANCE);

			/* if rising above critical rate, back to steady1=0 */
			if (t1 > SPINUP_TOLERANCE)
			{
				bgc_printf(BV_DIAG, "\nSWITCH BACK\n");

				steady1 = 0;
				rising = 1;
			}

			//bgc_printf(BV_DIAG, "spinyears = %d rising = %d steady2 = %d\n",spinyears,
			//	rising,steady2);
			//bgc_printf(BV_DIAG, "metcycle = %d tally1 = %lf tally2 = %lf pdif = %lf\n\n",
			//	metcycle,tally1,tally2,t1);

			metcycle = 0;
		}
		else
		{
			bgc_printf(BV_DIAG, "spinyears = %d rising = %d steady1 = %d\n",spinyears,
				rising,steady1);
			bgc_printf(BV_DIAG, "metcycle = %d tally1 = %lf tally2 = %lf pdif = %lf\n",
				metcycle,tally1,tally2,t1);

			metcycle++;
		}
	}

	/* end of do block, test for steady state */
	} while (mode == MODE_SPINUP && (!(steady1 && steady2) && (spinyears < ctrl.maxspinyears ||
		metcycle != 0)) );

	/* mode == MODE_SPINUP only */
	if (mode == MODE_SPINUP)
	{
		/* save some information on the end status of spinup */
		tally1b /= (double)nblock * 365.0;
		tally2b /= (double)nblock * 365.0;
		bgcout->spinup_resid_trend = (tally2b-tally1b)/(double)nblock;
		bgcout->spinup_years = spinyears;
	}
	
	/* RESTART OUTPUT HANDLING */
	/* if write_restart flag is set, copy data to the output restart struct */
	/* Removed 'write_restart' restriction to ensure that restart data are */
	/* available for spin and go operation.  WMJ 3/16/2005 */
	if (ok)
	{
		if (restart_output(&ctrl, &ws, &cs, &ns, &epv, metyr, 
			&(bgcout->restart_output)))
		{
			bgc_printf(BV_ERROR, "Error in call to restart_output() from bgc()\n");
			ok=0;
		}
		
		bgc_printf(BV_DIAG, "%d\t%d\tdone restart output\n",simyr,yday);
	}

	/* free phenology memory */
	if (ok && free_phenmem(&phenarr))
	{
		bgc_printf(BV_ERROR, "Error in free_phenmem() from bgc()\n");
		ok=0;
	}

	bgc_printf(BV_DIAG, "%d\t%d\tdone free phenmem\n",simyr,yday);
	
	/* free memory for local output arrays */	
	if (dayout) free(dayarr);
	if (ctrl.domonavg) free(monavgarr);
	if (ctrl.doannavg) free(annavgarr);
	if (ctrl.doannual) free(annarr);
	free(output_map);
	
	/* print timing info if error */
	if (!ok)
	{
		bgc_printf(BV_ERROR, "ERROR at year %d\n",simyr-1);
		bgc_printf(BV_ERROR, "ERROR at yday %d\n",yday-1);
	}
	
	/*Close files
	if(bgcin->cinit.frost)
		fclose(fp1);
*/

	bgcin->hModel.tmpHighFile.close();
	if (bgcin->hModel.active && mode != MODE_SPINUP)
		free(bgcin->hModel.stationFile);
	if (bgcin->laiM.active && mode != MODE_SPINUP)
		free(bgcin->laiM.laiFile);
	if (bgcin->gsiM.active && mode != MODE_SPINUP)
		free(bgcin->gsiM.gsiFile);
	/* return error status */	
	return (!ok);
}
