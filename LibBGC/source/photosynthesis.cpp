/*
photosynthesis.c
C3/C4 photosynthesis model

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/
#include "bgc.h"

int total_photosynthesis(const metvar_struct* metv, const epconst_struct* epc, 
			epvar_struct* epv, cflux_struct* cf, psn_struct *psn_sun, psn_struct *psn_shade)//,FILE *fp_ci, int mode)
{
	/* This function is a wrapper and replacement for the photosynthesis code which used
		to be in the central bgc.c code.  At Mott Jolly's request, all of the science code
		is being moved into funtions. */
	int ok=1;
	/* psn_struct psn_sun, psn_shade; */

	/* SUNLIT canopy fraction photosynthesis */
  /* set the input variables */
	psn_sun->c3 = epc->c3_flag;
	psn_sun->co2 = metv->co2;
	psn_sun->pa = metv->pa;
	psn_sun->t = metv->tday;
	psn_sun->lnc = 1.0 / (epv->sun_proj_sla * epc->leaf_cn);
	psn_sun->flnr = epc->flnr;
	psn_sun->ppfd = metv->ppfd_per_plaisun;
	/* convert conductance from m/s --> umol/m2/s/Pa, and correct
	for CO2 vs. water vapor */
	psn_sun->g = epv->gl_t_wv_sun * 1e6/(1.6*R*(metv->tday+273.15));
	psn_sun->dlmr = epv->dlmr_area_sun;
	if (ok && photosynthesis(psn_sun,metv))//,1,fp_ci,mode))
	{
		bgc_printf(BV_ERROR, "Error in photosynthesis() from bgc()\n");
		ok=0;
	}

	bgc_printf(BV_DIAG, "\t\tdone sun psn\n");

	epv->assim_sun = psn_sun->A;
	
	/* for the final flux assignment, the assimilation output
		needs to have the maintenance respiration rate added, this
		sum multiplied by the projected leaf area in the relevant canopy
		fraction, and this total converted from umol/m2/s -> kgC/m2/d */
	cf->psnsun_to_cpool = (epv->assim_sun + epv->dlmr_area_sun) *
		epv->plaisun * metv->dayl * 12.011e-9;

	/* SHADED canopy fraction photosynthesis */
	psn_shade->c3 = epc->c3_flag;
	psn_shade->co2 = metv->co2;
	psn_shade->pa = metv->pa;
	psn_shade->t = metv->tday;
	psn_shade->lnc = 1.0 / (epv->shade_proj_sla * epc->leaf_cn);
	psn_shade->flnr = epc->flnr;
	psn_shade->ppfd = metv->ppfd_per_plaishade;
	/* convert conductance from m/s --> umol/m2/s/Pa, and correct
	for CO2 vs. water vapor */
	psn_shade->g = epv->gl_t_wv_shade * 1e6/(1.6*R*(metv->tday+273.15));
	psn_shade->dlmr = epv->dlmr_area_shade;
	if (ok && photosynthesis(psn_shade,metv))//,0,fp_ci,mode))
	{
		bgc_printf(BV_ERROR, "Error in photosynthesis() from bgc()\n");
		ok=0;
	}

	bgc_printf(BV_DIAG, "\t\tdone shade_psn\n");

	epv->assim_shade = psn_shade->A;

	/* for the final flux assignment, the assimilation output
		needs to have the maintenance respiration rate added, this
		sum multiplied by the projected leaf area in the relevant canopy
		fraction, and this total converted from umol/m2/s -> kgC/m2/d */
	cf->psnshade_to_cpool = (epv->assim_shade + epv->dlmr_area_shade) *
		epv->plaishade * metv->dayl * 12.011e-9;
	return (!ok);
}


int photosynthesis(psn_struct* psn, const metvar_struct* metv)//,int sun,FILE *fp_ci, int mode)
{
	/*
	The following variables are assumed to be defined in the psn struct
	at the time of the function call:
	c3         (flag) set to 1 for C3 model, 0 for C4 model
	pa         (Pa) atmospheric pressure 
	co2        (ppm) atmospheric [CO2] 
	t          (deg C) air temperature
	lnc        (kg Nleaf/m2) leaf N concentration, per unit projected LAI 
	flnr       (kg NRub/kg Nleaf) fraction of leaf N in Rubisco
	ppfd       (umol photons/m2/s) PAR flux density, per unit projected LAI
	g          (umol CO2/m2/s/Pa) leaf-scale conductance to CO2, proj area basis
	dlmr       (umol CO2/m2/s) day leaf maint resp, on projected leaf area basis
	
	The following variables in psn struct are defined upon function return:
	Ci(Pa) intercellular [CO2]
	Ca         (Pa) atmospheric [CO2]
	O2         (Pa) atmospheric [O2]
	gamma      (Pa) CO2 compensation point, in the absence of maint resp.
	Kc         (Pa) MM constant for carboxylation
	Ko         (Pa) MM constant for oxygenation
	Vmax       (umol CO2/m2/s) max rate of carboxylation
	Jmax       (umol electrons/m2/s) max rate electron transport
	J          (umol RuBP/m2/s) rate of RuBP regeneration
	Av         (umol CO2/m2/s) carboxylation limited assimilation
	Aj         (umol CO2/m2/s) RuBP regen limited assimilation
	A          (umol CO2/m2/s) final assimilation rate
	*/
	
	/* the weight proportion of Rubisco to its nitrogen content, fnr, is 
	calculated from the relative proportions of the basic amino acids 
	that make up the enzyme, as listed in the Handbook of Biochemistry, 
	Proteins, Vol III, p. 510, which references:
	Kuehn and McFadden, Biochemistry, 8:2403, 1969 */
	static double fnr = 7.16;   /* kg Rub/kg NRub */
	
	/* the following enzyme kinetic constants are from: 
	Woodrow, I.E., and J.A. Berry, 1980. Enzymatic regulation of photosynthetic
	CO2 fixation in C3 plants. Ann. Rev. Plant Physiol. Plant Mol. Biol.,
	39:533-594.
	Note that these values are given in the units used in the paper, and that
	they are converted to units appropriate to the rest of this function before
	they are used. */
	/* I've changed the values for Kc and Ko from the Woodrow and Berry
	reference, and am now using the values from De Pury and Farquhar,
	1997. Simple scaling of photosynthesis from leaves to canopies
	without the errors of big-leaf models. Plant, Cell and Env. 20: 537-557. 
	All other parameters, including the q10's for Kc and Ko are the same
	as in Woodrow and Berry. */
	static double Kc25 = 404.0;   /* (ubar) MM const carboxylase, 25 deg C */ 
	static double q10Kc = 2.1;    /* (DIM) Q_10 for Kc */
	static double Ko25 = 248.0;   /* (mbar) MM const oxygenase, 25 deg C */
	static double q10Ko = 1.2;    /* (DIM) Q_10 for Ko */
	static double act25 = 3.6;    /* (umol/mgRubisco/min) Rubisco activity */
	static double q10act = 2.4;   /* (DIM) Q_10 for Rubisco activity */
	static double pabs = 0.85;    /* (DIM) fPAR effectively absorbed by PSII */
	
	/* local variables */
	int ok=1;	
	double t;      /* (deg C) temperature */
	double tk;     /* (K) absolute temperature */
	double Kc;     /* (Pa) MM constant for carboxylase reaction */
	double Ko;     /* (Pa) MM constant for oxygenase reaction */
	double act;    /* (umol CO2/kgRubisco/s) Rubisco activity */
	double Jmax;   /* (umol/m2/s) max rate electron transport */
	double ppe;    /* (mol/mol) photons absorbed by PSII per e- transported */
	double Vmax, J, gamma, Ca, Rd, O2, g;
	double a,b,c,det;
	double Av,Aj, A;

	double adp=1;
	/* begin by assigning local variables */
	g = psn->g;
	t = psn->t;
	tk = t + 273.15;
	Rd = psn->dlmr;
	
	/* convert atmospheric CO2 from ppm --> Pa */
	Ca = psn->co2 * psn->pa / 1e6;
	
	/* set parameters for C3 vs C4 model */
	if (psn->c3)
	{
		ppe = 2.6;
	}
	else /* C4 */
	{
		ppe = 3.5;
		Ca *= 10.0;
	}
	psn->Ca = Ca;		
	
	/* calculate atmospheric O2 in Pa, assumes 21% O2 by volume */
	psn->O2 = O2 = 0.21 * psn->pa;
	
	/* correct kinetic constants for temperature, and do unit conversions */
	Ko = Ko25 * pow(q10Ko, (t-25.0)/10.0);
	psn->Ko = Ko = Ko * 100.0;   /* mbar --> Pa */
	if (t > 15.0)
	{
		Kc = Kc25 * pow(q10Kc, (t-25.0)/10.0);
		act = act25 * pow(q10act, (t-25.0)/10.0);
	}
	else
	{
		Kc = Kc25 * pow(1.8*q10Kc, (t-15.0)/10.0) / q10Kc;
		act = act25 * pow(1.8*q10act, (t-15.0)/10.0) / q10act;
	}
	psn->Kc = Kc = Kc * 0.10;   /* ubar --> Pa */
	act = act * 1e6 / 60.0;     /* umol/mg/min --> umol/kg/s */
	
	/* calculate gamma (Pa), assumes Vomax/Vcmax = 0.21 */
	psn->gamma = gamma = 0.5 * 0.21 * Kc * psn->O2 / Ko;
	 
	/* calculate Vmax from leaf nitrogen data and Rubisco activity */
	
	/* kg Nleaf   kg NRub    kg Rub      umol            umol 
	   -------- X -------  X ------- X ---------   =   --------
	      m2      kg Nleaf   kg NRub   kg RUB * s       m2 * s       
	   
	     (lnc)  X  (flnr)  X  (fnr)  X   (act)     =    (Vmax)
	*/
		/* CO2 multiplier 
	if(metv->co2 < 350)
		adp = 1.0;
	else
		adp = metv->co2 * (-0.000286)+1.1;
	if(metv->co2 < 350)
		adp = 1.0;
	else
		adp = metv->co2 * (-0.000429)+1.15;*/

	adp=1.0;
	psn->Vmax = Vmax = psn->lnc * psn->flnr * fnr * act * adp ;
	
	/* calculate Jmax = f(Vmax), reference:
	Wullschleger, S.D., 1993.  Biochemical limitations to carbon assimilation
		in C3 plants - A retrospective analysis of the A/Ci curves from
		109 species. Journal of Experimental Botany, 44:907-920.
	*/
	psn->Jmax = Jmax = 2.1*Vmax;
	
	/* calculate J = f(Jmax, ppfd), reference:
	de Pury and Farquhar 1997
	Plant Cell and Env.
	*/
	a = 0.7;
	b = -Jmax - (psn->ppfd*pabs/ppe);
	c = Jmax * psn->ppfd*pabs/ppe;
	psn->J = J = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
	
	/* solve for Av and Aj using the quadratic equation, substitution for Ci
	from A = g(Ca-Ci) into the equations from Farquhar and von Caemmerer:
	     
	       Vmax (Ci - gamma)
	Av =  -------------------   -   Rd
	      Ci + Kc (1 + O2/Ko)
	
	
	         J (Ci - gamma)
	Aj  =  -------------------  -   Rd
           4.5 Ci + 10.5 gamma  
    */

	/* quadratic solution for Av */    
    a = -1.0/g;
    b = Ca + (Vmax - Rd)/g + Kc*(1.0 + O2/Ko);
    c = Vmax*(gamma - Ca) + Rd*(Ca + Kc*(1.0 + O2/Ko));
    
    if ((det = b*b - 4.0*a*c) < 0.0)
    {
    	bgc_printf(BV_ERROR, "negative root error in psn routine\n");
    	ok=0;
    }
    
    psn->Av = Av = (-b + sqrt(det)) / (2.0*a);
    
    /* quadratic solution for Aj */
	a = -4.5/g;    
	b = 4.5*Ca + 10.5*gamma + J/g - 4.5*Rd/g;
	c = J*(gamma - Ca) + Rd*(4.5*Ca + 10.5*gamma);
		
	if ((det = b*b - 4.0*a*c) < 0.0)
	{
		bgc_printf(BV_ERROR, "negative root error in psn routine\n");
		ok=0;
	}
	
	psn->Aj = Aj = (-b + sqrt(det)) / (2.0*a);
	
	/* estimate A as the minimum of (Av,Aj) */
	if (Av < Aj) A = Av; 
	else         A = Aj;
	psn->A = A;
	bgc_printf(BV_DIAG, "psn->A: %f, A: %f\n", psn->A, A);
	psn->Ci = Ca - (A/g);
	/***************************************************/
	//if(mode==2&&sun==1)
		//fprintf(fp_ci,"%lf%12.4f%12.4f%\n",Ca,psn->Ci,Ca/psn->Ci);
	return (!ok);
}	


//  for high time resolution
int total_photosynthesisTimeRes(high_time_resolution* high_time_resolution, const wflux_struct* wf, std::vector<StationDataFlux*> &sfData, const cstate_struct* cs, const double albedo, const metvar_struct* metv, const epconst_struct* epc,
	epvar_struct* epv, cflux_struct* cf, psn_struct *psn_sun, psn_struct *psn_shade, const int yearS, const int daysS)//,FILE *fp_ci, int mode)
{

	/* This function is a wrapper and replacement for the photosynthesis code which used
		to be in the central bgc.c code.  At Mott Jolly's request, all of the science code
		is being moved into funtions. */
	int ok = 1;
	/* psn_struct psn_sun, psn_shade; */

	/* SUNLIT canopy fraction photosynthesis */
  /* set the input variables */
	psn_sun->c3 = epc->c3_flag;
	psn_sun->co2 = metv->co2;
	psn_sun->pa = metv->pa;
	psn_sun->t = metv->tday;
	psn_sun->lnc = 1.0 / (epv->sun_proj_sla * epc->leaf_cn);
	psn_sun->flnr = epc->flnr;
	psn_sun->ppfd = metv->ppfd_per_plaisun;
	/* convert conductance from m/s --> umol/m2/s/Pa, and correct
	for CO2 vs. water vapor */
	psn_sun->g = epv->gl_t_wv_sun * 1e6 / (1.6*R*(metv->tday + 273.15));
	psn_sun->dlmr = epv->dlmr_area_sun;
	if (ok && photosynthesisTimeRes(high_time_resolution, wf, sfData, epc, epv, cs, albedo, psn_sun, metv, yearS, daysS, 1))//,1,fp_ci,mode))
	{
		bgc_printf(BV_ERROR, "Error in photosynthesis() from bgc()\n");
		ok = 0;
	}

	//std::cout << cs->leafc << " " << yearS << " " << daysS << std::endl;
	bgc_printf(BV_DIAG, "\t\tdone sun psn\n");

	epv->assim_sun = psn_sun->A;

	/* for the final flux assignment, the assimilation output
		needs to have the maintenance respiration rate added, this
		sum multiplied by the projected leaf area in the relevant canopy
		fraction, and this total converted from umol/m2/s -> kgC/m2/d */
	cf->psnsun_to_cpool = (epv->assim_sun + epv->dlmr_area_sun) *
		epv->plaisun * 24*60*60 * 12.011e-9;

	/* SHADED canopy fraction photosynthesis */
	psn_shade->c3 = epc->c3_flag;
	psn_shade->co2 = metv->co2;
	psn_shade->pa = metv->pa;
	psn_shade->t = metv->tday;
	psn_shade->lnc = 1.0 / (epv->shade_proj_sla * epc->leaf_cn);
	psn_shade->flnr = epc->flnr;
	psn_shade->ppfd = metv->ppfd_per_plaishade;
	/* convert conductance from m/s --> umol/m2/s/Pa, and correct
	for CO2 vs. water vapor */
	psn_shade->g = epv->gl_t_wv_shade * 1e6 / (1.6*R*(metv->tday + 273.15));
	psn_shade->dlmr = epv->dlmr_area_shade;
	if (ok && photosynthesisTimeRes(high_time_resolution, wf, sfData, epc, epv, cs, albedo, psn_shade, metv, yearS, daysS, 0))//,0,fp_ci,mode))
	{
		bgc_printf(BV_ERROR, "Error in photosynthesis() from bgc()\n");
		ok = 0;
	}

	bgc_printf(BV_DIAG, "\t\tdone shade_psn\n");

	epv->assim_shade = psn_shade->A;

	/* for the final flux assignment, the assimilation output
		needs to have the maintenance respiration rate added, this
		sum multiplied by the projected leaf area in the relevant canopy
		fraction, and this total converted from umol/m2/s -> kgC/m2/d */
	cf->psnshade_to_cpool = (epv->assim_shade + epv->dlmr_area_shade) *
		epv->plaishade * 24 * 60 * 60 * 12.011e-9;
	return (!ok);
}

//  for high time resolution
int photosynthesisTimeRes(high_time_resolution* high_time_resolution, const wflux_struct* wf, std::vector<StationDataFlux*> &sfData, const epconst_struct* epc, epvar_struct* epv,
	const cstate_struct* cs, const double albedo, psn_struct *psn, const metvar_struct* metv, 
	const int yearS, const int daysS, const int sunorshade)//,int sun,FILE *fp_ci, int mode)
{
	int ok = 1;
	double meanT = 0;

	int dayCurr = yearS * 365 + daysS;

	/*if (dayCurr == 729)
		std::cout << yearS << " " << dayCurr << std::endl;*/	

	// time resolution
	int timeRes = sfData[1]->HRMIN - sfData[0]->HRMIN;

	// total final assimilation rate
	double totalA = 0;  // used for calculating temperature
	double totalPhotoNum = 0;  // used for calculating carbon or rad limitation

	int totalNum = 0;
	for (int i = 0; i < 24 * 60/ timeRes; ++i)
	{
		psn_struct psnT = *psn;
		metvar_struct metvT = *metv;		
		wflux_struct wfT = *wf;
		epvar_struct epvT = *epv;	

		// update temperature
		if (i == 0) 
		{
			psnT.t = metv->tday;
		}
		else
		{
			psnT.t = metv->tnight;
		}
		//psnT.t = sfData[dayCurr * 24 * 60 / timeRes + i]->TA;
		metvT.tday = psnT.t;
		metvT.tmin = psnT.t;
		meanT += metvT.tday;

		// update par and swavgfd
		metvT.swavgfd = sfData[dayCurr* 24 * 60 / timeRes + i]->Srad;
		metvT.par = metvT.swavgfd * RAD2PAR;	

		//std::cout << yearS << " " << dayCurr << " " << dayCurr * 24 * 60 / timeRes + i << " " << psnT.t << " " <<
		//	metvT.swavgfd << " " << sfData[dayCurr * 24 * 60 / timeRes + i]->VPD << std::endl;

		// update calculate ppfd 
		psnT.ppfd = calppfdT(cs, epc, &metvT, &epvT, albedo, sunorshade);
		/* for hightime stress output*/
		if (psnT.ppfd <= 0 )
		{
			if(sunorshade == 1 && high_time_resolution->output_stress) 
				high_time_resolution->tmpHighFile << ", ";
			continue;
		}
		// update calculate g
		metvT.vpd = sfData[dayCurr * 24 * 60 / timeRes + i]->VPD;
		double gl = calGl(&metvT, epc, &epvT, &wfT, 1, sunorshade);
		//psnT.g = gl * 1e6 / (1.6*R*(psnT.t + 273.15));	

		if (!photosynthesis(&psnT, &metvT))
		{
			totalA += psnT.A;
			//std::cout << psnT.t << ", " << psnT.g  << ", " << psnT.A << std::endl;
			++totalNum;
		}

		/* for hightime stress output*/
		if (psnT.Av > psnT.Aj )
		{
			if(high_time_resolution->output_stress)
				high_time_resolution->tmpHighFile << ", " << 0;	// rad limited
		}
		else 
		{
			if(high_time_resolution->output_stress)
				high_time_resolution->tmpHighFile << ", " << 1;	// carbon limited
			totalPhotoNum += 1;
		}
	}

	psn->A = totalA / totalNum;
	psn->t = meanT / (24 * 60 / timeRes);
	psn->O2 = totalPhotoNum / totalNum;  // carbon limitation percent
	return (!ok);
}


double simulationPar(const double inPar, const double timePer)
{
	if (timePer <= 0.5)
	{
		return inPar * timePer;
	}
	else
	{
		return inPar * fabs(1 - timePer);
	}
}

double simulationPar1(const double inPar)
{
	return inPar;
}

void replacePhotosynthesisResults(high_time_resolution* high_time_resolution, epvar_struct*epv, cflux_struct*cf, psn_struct*psn_sun, psn_struct*psn_shade,
	const epvar_struct*epvT, const cflux_struct*cfT, const psn_struct*psn_sunT, const psn_struct*psn_shadeT, const int yearS, const int daysS)
{
	// write to disk
	high_time_resolution->tmpHighFile << yearS << ", " << daysS << ", " << psn_sun->t << ", " << psn_sunT->t << ", " << 
		cf->psnsun_to_cpool << ", " << cfT->psnsun_to_cpool << ", " << cf->psnshade_to_cpool << ", " << cfT->psnshade_to_cpool;
	
	if (psn_sun->Av > psn_sun->Aj)
	{
		high_time_resolution->tmpHighFile << ", " << 0;	// rad limited
	}
	else
	{
		high_time_resolution->tmpHighFile << ", " << 1;	// carbon limited
	}

	// carbon limited percentage
	if (psn_sunT->O2 > 101)
	{
		high_time_resolution->tmpHighFile << ", " << -1;
	}
	else
	{ 
		high_time_resolution->tmpHighFile << ", " << psn_sunT->O2;
	}
	
	high_time_resolution->tmpHighFile << std::endl;	

	if (high_time_resolution->active && high_time_resolution->output_old_cpool == false)
	{
		cf->psnsun_to_cpool = cfT->psnsun_to_cpool;
		cf->psnshade_to_cpool = cfT->psnshade_to_cpool;
		psn_sun->A = psn_sunT->A;
		psn_shade->A = psn_shadeT->A;
	}

}


double calppfdT(const cstate_struct* cs, const epconst_struct* epc,
	metvar_struct* metv, epvar_struct* epv, double albedo, const int sunorshade)
{
	radtrans(cs, epc, metv, epv, albedo);

	if (sunorshade == 1)
	{
		return metv->ppfd_per_plaisun;
	}
	else
	{
		return metv->ppfd_per_plaishade;
	}
}

double calGl(const metvar_struct* metv, const epconst_struct* epc,
	epvar_struct* epv, wflux_struct* wf, const int mode, const int sunorshade)
{
	canopy_et(metv, epc, epv, wf, mode);

	if (sunorshade == 1)
	{
		// sun
		return epv->gl_t_wv_sun;
	}
	else
	{
		// shade
		return epv->gl_t_wv_shade;
	}
}

// read lai to vector 
int readLaiData(std::vector<float> &laiData, const char* laiDataFile, const int yearS)
{
	int ok = 1;
	//int ndays = 365 * yearS;
	std::ifstream infile(laiDataFile);
	//int i = 0;

	if (infile.is_open()) {
		std::string lineStr;
		while (getline(infile, lineStr)) {
			laiData.push_back(std::stof(lineStr));
		}
		infile.close();
	}
	else
	{
		bgc_printf(BV_ERROR, "Error in open lai data from bgc()\n");
		ok = 0;
	}

	return ok;
}

// read high time resolution data 
int readStationFluxData(std::vector<StationDataFlux*> &sfData, const char* fluxStationDataFile, const int yearS)
{
	int ok = 1;
	int ndays = 365 * yearS;
	std::ifstream infile(fluxStationDataFile);
	int i = 0;

	if (infile.is_open()) {
		std::string lineStr;
		while (getline(infile, lineStr)) {
			
			++i;
			//if (i == 1)continue;			

			// remove blank
			int index = 0;
			if (!lineStr.empty())
			{
				while ((index = lineStr.find(' ', index)) != std::string::npos) lineStr.erase(index, 1);
			}
			index = 0;
			if (!lineStr.empty())
			{
				while ((index = lineStr.find('"', index)) != std::string::npos) lineStr.erase(index, 1);
			}

			std::vector<std::string> resultStr;
			SplitString(lineStr, resultStr, ",");

			StationDataFlux *ff = new StationDataFlux();
			ff->YEAR = std::stoi(resultStr[0]);
			ff->DAY = std::stoi(resultStr[1]);
			ff->HRMIN = std::stoi(resultStr[2]);
			/*ff->Rd = std::stold(resultStr[3]);
			ff->GPP = std::stold(resultStr[4]);
			ff->NEE = std::stold(resultStr[5]);
			ff->TA = std::stold(resultStr[6]);
			ff->PREC = std::stold(resultStr[7]);
			ff->RH = std::stold(resultStr[8]);
			ff->VPD = std::stold(resultStr[9]);
			ff->Srad = std::stold(resultStr[10]);
			ff->SWC = std::stold(resultStr[11]);
			ff->PAR = std::stold(resultStr[12]);*/

			ff->TA = std::stold(resultStr[3]);
			ff->PREC = std::stold(resultStr[4]);
			ff->VPD = std::stold(resultStr[5]);
			ff->Srad = std::stold(resultStr[6]);

			ff->GPP = std::stold(resultStr[7]);
			ff->NEE = std::stold(resultStr[8]);
			

			sfData.push_back(ff);
		}
		infile.close();
	}
	else
	{
		bgc_printf(BV_ERROR, "Error in open high resolution data for photosynthesis() from bgc()\n");
		ok = 0;
	}

	return ok;
}

void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));

		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
		v.push_back(s.substr(pos1));
}

void analysisComm(const int argc, char **argv, high_time_resolution* highTM, lai_model* laiM)
{
	//char * tmp = nullptr;
	std::string inStationFile = "";
	laiM->active = false;

	highTM->active = false;
	highTM->output_stress = false;
	highTM->output_old_cpool = false;

	for (size_t i = 1; i < argc; i++)
	{
		if (std::string(argv[i]) == "-outOldC")
		{
			highTM->output_old_cpool = true;
		}

		if (std::string(argv[i]) == "-outStress")
		{
			highTM->output_stress = true;
		}
		// high time resolution
		if (std::string(argv[i]) == "-h")
		{
			if (i + 1 >= argc)
			{
				inStationFile = "wrong";
				//return tmp;
				std::cout << "please provide the high resoltion station file;" << std::endl;
			}
			else
			{
				inStationFile = std::string(argv[i + 1]);
				highTM->active = true;
				highTM->stationFile = (char *)malloc((inStationFile.length() + 1) * sizeof(char));
				strcpy(highTM->stationFile, inStationFile.c_str());

				//highTM->stationFile = inStationFile;
			}			
		}
		// ndvi
		if(std::string(argv[i]) == "-lai")
		{
			if (i + 1 >= argc)
			{
				inStationFile = "wrong";
				std::cout << "please provide the lai file;" << std::endl;
				//return tmp;
			}
			else
			{
				inStationFile = std::string(argv[i + 1]);
				laiM->active = true;
				laiM->laiFile = (char *)malloc((inStationFile.length() + 1) * sizeof(char));
				strcpy(laiM->laiFile, inStationFile.c_str());
			}
		}

	}
	//if(inStationFile == "")return tmp;

	//tmp = (char *)malloc((inStationFile.length() + 1) * sizeof(char));
	//strcpy(tmp, inStationFile.c_str());

	//return tmp;
}