#include "bgc.h"

/* Saturation vapor pressure Tsate (Pa), as a function of temperature (deg C),Temperature in Ceisiul degree???? 
  are calculated from the eighth-order polynomial fits of Flatau et al. (1992) */
double vaporpressure_sat(double temp)
{ 
	double a0,a1,a2,a3,a4,a5,a6,a7,a8;
    double t=temp;
    double pvap;
	if(temp >= 0)
	{
      a0=6.11213476;
      a1=0.444007856;
      a2=1.4306234e-2;
      a3=2.64461437e-4;
      a4=3.05903558e-6;
	  a5=1.96237241e-8;
	  a6=8.92344772e-11;
	  a7=-3.73208410e-13;
	  a8=2.09339997e-16;
	}
	else
	{
      a0=6.11123616;
      a1=0.503109514;
      a2=0.0188369801;
      a3=4.20547422e-4;
      a4=6.14396778e-6;
	  a5=6.02780717e-8;
	  a6=3.87940929e-10;
	  a7=1.49436277e-12;
	  a8=2.62655803e-15;
	}

  pvap = 100*(a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*pow(t,5)+a6*pow(t,6)+a7*pow(t,7)+a8*(pow(t,8)));

  /* A very good approximation can usually be made using the August-Roche-Magnus formula 
  (usually called the Magnus or Magnus-Tetens approximation, though this is historically inaccurate[5]):
  
  pvap=610.94*exp(17.625*t/(t+243.04*t);   //pvap(T) is the equilibrium or saturation vapor pressure in Pa, which is a function of temperature; T is in Celsius.
 经比较两种方法结果近似相等：t>0,yy方案1比方案2值稍大1~2pa;t<0,方案1比方案2小10~20pa.没有进行实测值验证。
 */
  return pvap;
}
/* the derivation of saturagtion vapour pressure,as a function of temperature T(centigrade),are calculated from the eighth-order polynomial fits of flatau etal.1992*/
double svp_dt(double temp)
{
	double a0,a1,a2,a3,a4,a5,a6,a7,a8;
    double t=temp;
    double pvap;
	if(temp>=0)
	{
      a0=0.444017302;
      a1=2.86064092e-2;
      a2=7.94683137e-4;
      a3=1.21211669e-5;
      a4=1.03354611e-7;
	  a5=4.04125005e-10;
	  a6=-7.88037859e-13;
	  a7=-1.14596802e-14;
	  a8=3.81294516e-17;
	}
	else
	{
      a0=5.03277922e-1;
      a1=3.77289173e-2;
      a2=1.26801703e-3;
      a3=2.49468427e-5;
      a4=3.13703411e-7;
	  a5=2.57180651e-9;
	  a6=1.33268878e-11;
	  a7=3.94116744e-14;
	  a8=4.98070196e-17;
	}
    pvap=100*(a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*pow(t,5)+a6*pow(t,6)+a7*pow(t,7)+a8*(pow(t,8)));
  return pvap;
}
double specifichumiditysat_dt(double tem,metvar_struct* metv)
{
	double Patm,esat,esat_dt,qsat_dt;
	Patm = metv->pa;
	esat=vaporpressure_sat(tem);
    esat_dt=svp_dt(tem);
	qsat_dt=0.622*Patm/((Patm-0.378*esat)*(Patm-0.378*esat))*esat_dt;
	return qsat_dt;
}
/* ground_heatflux.c for calcculating ground surface sensible heat  flux  W/m^2 */
//void sensible_heatflux(double* vaporheatflux,double* latenthf,double* rb)
void sensible_heatflux(double* vaporheatflux,double* latenthf,double* vhf_dt,double* lhf_dt,double* rb ,int snl,metvar_struct* metv,pmet_struct *pmetvar,soilpar_struct* vspar)
{
	double ra;
	double rs;                             /* (s/m) the aerodynamic resistance ;soil evapor resistance */
	double psi_m=0,psi_h=0;                /* psi_m is the stability correction for buoyancy effects on the momentum flux; psi_h is the stability correction for heat transportx */
	double k=0.41;                          /* karman constant */
	double z=2.0;                          /* (m)measured height,the reference height of temperature measurement  */
	double d=0.06;                         /* d is the zero-plane displacement height (m) d=0.066h  h is the height of canopy */
	double zm=0.02;                        /* the roughness length governing momentum (m), zm=0.13h h is the height of canopy; ln(zm/zh)=0.24 */
	double zh=0.02;                        /* (m)the surface roughness length for the heat flux (vapor transfer) (m) */
	double U;                            /* (m s-1) the friction velocity ,for simple,similar(approximate)to wind speed */
	double ts,ta;                            /* deg C  soil ,air temperature */
	double rou_a,tk;                          /* (kg/m3) the density of air */
	double ea_sat;                          /*  saturation vapor pressure of air */
	double es_sat;                          /*  saturation vapor pressure of ground */  
	double es,es_dt;                              /*  vapor pressure of ground */
	double ea;                              /*  vapor pressure of air */
	double Hs;
	double psi = svar[5].s_h;

	// zm=zh*exp(0.24);
	ta = pmetvar->ta;                         /* the temperature of daytime */
	U = metv->wspeed + 1;                    /*  wind speed */
	ts = svar[snl].s_t;//+273.15;                        /* the temperature of ground surface */
	tk = ts + 273.15;                          /* the unit transform to K from deg C */
	rou_a = 1.292-(0.00428*ta);              /* air density at temperature ta */
	if(snl<5)
	{
		d=0;
		zm=0.0002;
		Hs = exp(-Liv/Rw*(1/tk-1/Tk));
	}
	else
		Hs = exp(FSM_g * psi / (Rw * tk));                     /* relative vapor pressure */
	ra = 17.05;//((log((z-d)/zm)-psi_m)*(log((z-d)/zh)-psi_h))/(k*k*U);       /* By Campbell,1985  ra=17.05s/m */
	// rs=exp(8.2-4.225*svar[0].s_w/vspar->porosity);                   //Soil resistance is affected by soil type and soil water content and can be obtained from an empirical formula (Sun et al. 1998)
	rs = 3.5*pow(vspar->porosity/svar[5].s_w,2.3)+33.5;                /* By linjiading et al.,1983 */
	// rs=-805+4140*(vspar->porosity-svar[0].s_w);                      /* Camillo et al.,1986 */
	//饱和水汽压，单位：pa
	ea_sat = vaporpressure_sat(ta);
	es_sat = vaporpressure_sat(ts);
	*rb = ra;
	// 
	es_dt = svp_dt(ts);
	//空气水汽压,百帕
	ea = 0.01 * metv->rhumid * ea_sat;
	//ea = metv->vap;
	es = Hs * es_sat;
	/* evapor latent heat fluxes J/m2s = w/m2*/
	*latenthf = rou_a * CP * (es-ea)/(r*(ra+rs));

	/* counting air sensible heatflux (W/m2) */
	*vaporheatflux = rou_a * CP * (ts-ta)/(ra);       /* CP is specific heat of air. define constant */
	*vhf_dt = rou_a * CP / ra;
	*lhf_dt = rou_a * CP * Hs * es_dt/(r * (ra + rs));
}

/* sun net radiation Rn=LE+H+G (W/m^2) 
   Rn=(1-a)*Rs +R1-Ls ; 
   Li, S. B., and Zhao, C. Y. 2006. Estimating evapotranspiration based on energy balance in Guanchuan River Basin using remote sensing. Remote Sensing Technology and Application, 21(6), 521-526. (in Chinese)
   R1=Ea*SBC*pow(Ta,4)  Ta is the temperature of air (K)  
   SBC      5.67e-8         (W/(m2 K4)) Stefan-Boltzmann constant 
   Ea is long wave atmosphere albedo rato: Ea=1.24pow(Ed/Ta,1/7)
   Ed is vapour pressure (mb )
   
   Ls=Es*SBC*pow(Ts,4)  Ts is the temperature of ground (K).
   Es is long wave  ground albedo rato:Es=0.94+0.014*whi; whi is ground water content

   */
void sunnetradiation(double* Rn,double* Rn_dt, int snl,pmet_struct *pmetvar)
{ 
  //extern metvar_struct metvar;
  //extern soilvar_struct svar[10];
                                /* Net radiation of the accepted ground */
	//TODO 需要修改为bgc中的albedo
	double a = 0.15;                          /* ground Surface albedo rato */
	double Rs;                              /* total sun shortwave radiation */
	double Rl;                              /* atmosphere long wave radiation */
	double Ls;                              /* ground longwave radiation */
	double Ea;                              /* atmosphere longwave radiation albedo */
	double Es;                              /* ground longwave radiation albedo */
//	double Ed;                              /* (mb)atmosphere vapour pressure */
	double ea_sat;                           /* (pa)saturation vapor presssure */
	double ts,ta;                              /* deg C the temperature of the ground */
	double groundwater;                              /* the ground water content */

	//TODO dsr需要根据降雪的天数调整，在BGC中需要增加的量。
	int dsr = 5;   //连续未下雨天数，与BGC耦合时需考虑 /* number of days since rain*/

	if(snl==5)
		a = 0.20;
	else
		a = 0.60;
	//a = 0.85 * pow(0.92,pow(dsr,0.58));
	/* The snow surface albedo is assumed to decay with age following the functional form descrobed by laramie and schaake(1972)*/
	if(pmetvar->ta<0.0)
		a=0.85*pow(0.92,pow(dsr,0.58));  //accumulation season
	else
		a=0.85*pow(0.70,pow(dsr,0.58));  //melt season
	ts = svar[snl].s_t + 273.15;
	groundwater = svar[snl].s_w;
	ta = pmetvar->ta + 273.15;                   /* ta is the temperature of air */
	// Ed= 1013.25;                            /* 1013.25mb=1 atmospere pressure(760mmHg) */
	ea_sat = vaporpressure_sat(ta);    /* call vaporpressure_sat() */
	
	//Ed = ea_sat * metv->rhumid * 100;              /* unit pa->mb(毫巴=百巴=100pa）Ed=metv->vap/100;  */
	//Ea=1.24*pow(Ed/ta,1/7);
	Ea = 0.96;	
	//TODO BGC中，需要用向下的长波辐射代替。
	Rl = Ea * SBC * pow(ta,4);
	Es = 0.94 + 0.14 * groundwater;                        /* By zheng xiuqing,2001 */
	if(Es > 1.0)
		Es = 1.0;
	Ls = Es * SBC * pow(ts,4);
	Rs = pmetvar->irad;                       /* (W/m2)  daylight average shortwave flux */
	*Rn = (1-a) * Rs + Rl - Ls;
	/* count to the partial derivative of the longwave radiation absorbed by the ground(positive toward the atmosphere.*/
	*Rn_dt=-4 * Es * SBC * pow(ts,3);
}


/*气候数据************************************************************************************************/
/* 
daymet.c
transfer one day of meteorological data from metarr struct to metv struct

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions from version 4.1.1:
	changed the coefficient for Tday calculation from 0.212 to 0.45, to be
	consistent with the calculation of VPD in MTCLIM and Daymet.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int sdaymet(const metarr_struct* metarr, metvar_struct* metv, int metday)//,double* ytd_snd,double* ytd_snp)
{
	/* generates daily meteorological variables from the metarray struct */
	int ok=1;
	
	//*ytd_snd = metarr->snd[metday-1];
	//*ytd_snp = metarr->snp[metday-1];
	
	/* convert prcp from 0.1mm --> cm */
	metv->prcp = metarr->prcp[metday];
    //metv->snd = metarr->snd[metday];
   // metv->snp = metarr->snp[metday];
	
	//metv->den = metarr->den[metday];
	/* air temperature calculations (all temperatures from 0.1deg C to 1 deg C) */
	metv->tmax = metarr->tmax[metday];
	metv->tmin = metarr->tmin[metday];
	metv->tavg = metarr->tavg[metday];
	
	/*convert atmosphere(vapor) press from 0.1hpa to 1pa */
	//metv->atp = metarr->atp[metday];
	//metv->vap = metarr->vap[metday];
	
    /* convert windspeed from 0.1m/s to m/s */
	metv->wspeed =metarr->wspeed[metday];
	
	/* daylight average	shortwave flux density (W/m2) */
	metv->swavgfd =  metarr->swavgfd[metday];
	
	/* relative humidity */
	metv->rhumid = metarr->rhumid[metday];
	/* calculate VPD */
	//metv->vpd = metv->vap*(100/metv->rhumid - 1);
	metv->vpd = metarr->vpd[metday];
	/* daylength (s) from 0.1hour to s */
	metv->dayl = metarr->dayl[metday];

	return (!ok);
}

/*参数************************************************************************************************/

/* 
baresoil_evap.c
daily bare soil evaporation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
swht-BGC version 1.0 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int sparameters(soilpar_struct* vspar) 
{
	/*
	int typeNo;              soil texture types code numer: 1 sand,2 loamysand,3 sandyloam,4 loam,5 siltyloam,6 sandyclayloam,7 clayloam,8 siltyclayloam,9 sandyclay,10 siltyclay,11 clay,12 peat,13 moss,14 lichen 
	double psandy;           (%)percentage of sand 
	double psilt;           (%)percentage of slit 
	double pclay;            (%)percentage of clay 
	double k_solids;         (Wm-1k-1)the thermal conductivity of the soil solids 
	double c_solids;         (the heat capacity of the soil solids 
	double k_sat;            (mms-1)the saturated hydraulic conductivity 
	double psi_sat;          (mm)the saturated matric potential 
	double porosity;         the saturated volumetric water content(porosity) 
	double CH_b;             the Clap and Hornberger constant 
	*/

	int ok=1;
	
	/* calculating the specified parameters for each texture type by their original percentage of sand,silt, and clay. By Cosby et al.1984 */
	//vspar->k_solids=(8.80*vspar->psandy+2.92*vspar->pclay)/(vspar->psandy+vspar->pclay);            /* (Wm-1k-1)the thermal conductivity of the soil solids */
	vspar->c_solids=(2.128*vspar->psandy+2.385*vspar->pclay)/(vspar->psandy+vspar->pclay)*1000000;  /* (Jm-3k-1)the heat capacity of the soil solids */
	vspar->k_sat=0.0070556*pow(10,-0.884+0.0153*vspar->psandy)*0.001;                         /* (ms-1)the saturated hydraulic conductivity */
	vspar->psi_sat = -10 * pow(10,1.88-0.013*vspar->psandy)*0.001;                                /* (m)the saturated matric potential */
	vspar->porosity = 0.489 - 0.00126 * vspar->psandy;                                              /* the saturated volumetric water content(porosity) */
	vspar->CH_b = 2.91 + 0.159 * vspar->pclay;                                                     /* the Clap and Hornberger constant */    
	//vspar->c_solids = 2300000;
	vspar->k_solids = 1.0;
	
	//vspar->psi_sat = -(exp((1.54 - 0.0095*vspar->psandy + 0.0063* vspar->psilt)*log(10.0))*9.8e-5);
	/*Biome-BGC中的公式，vwc_sat即porosity*/
	//
	//vspar->CH_b = -(3.10 + 0.157*vspar->pclay - 0.003 * vspar->psandy);                                                      /* the Clap and Hornberger constant */ 
	//sitec->soil_b = -(3.10 + 0.157*clay - 0.003*sand);                                      /* (DIM) Clapp-Hornberger "b" parameter */
	//sitec->vwc_sat = (50.5 - 0.142*sand - 0.037*clay)/100.0;
	   

	return(!ok);





	/*
	switch(typeNo)
	{
	case 1 :
		vspar->psandy=92;
		vspar->psilt=5;
		vspar->pclay=3;
		vspar->k_solids=8.6143;
		vspar->c_solids=2136116;
		vspar->k_sat=0.000023558;
		vspar->psi_sat=-0.04729;
		vspar->porosity=0.3731;
		vspar->CH_b=3.39;
		break;
	case 2:
        vspar->psandy=82;
		vspar->psilt=12;
		vspar->pclay=6;
		vspar->k_solids=8.3991;
		vspar->c_solids=2145523;
		vspar->k_sat=0.000016563;
		vspar->psi_sat=-0.06394;
		vspar->porosity=0.3857;
		vspar->CH_b=3.86;
		break;
    case 3:
        vspar->psandy=58;
		vspar->psilt=32;
		vspar->pclay=10;
		vspar->k_solids=7.9353;
		vspar->c_solids=2165794;
		vspar->k_sat=0.000007111;
		vspar->psi_sat=-0.13188;
		vspar->porosity=0.4159;
		vspar->CH_b=4.50;
		break;
    case 4:
        vspar->psandy=43;
		vspar->psilt=39;
		vspar->pclay=18;
		vspar->k_solids=7.0649;
		vspar->c_solids=2203836;
		vspar->k_sat=0.000004912;
		vspar->psi_sat=-0.20734;
		vspar->porosity=0.4348;
		vspar->CH_b=5.77;
		break;
    case 5:
        vspar->psandy=17;
		vspar->psilt=70;
		vspar->pclay=13;
		vspar->k_solids=6.2520;
		vspar->c_solids=2239367;
		vspar->k_sat=0.000001677;
		vspar->psi_sat=-0.45425;
		vspar->porosity=0.4676;
		vspar->CH_b=4.98;
		break;
    case 6:
        vspar->psandy=58;
		vspar->psilt=15;
		vspar->pclay=27;
		vspar->k_solids=6.9323;
		vspar->c_solids=2209635;
		vspar->k_sat=0.000007111;
		vspar->psi_sat=-0.13188;
		vspar->porosity=0.4159;
		vspar->CH_b=7.20;
		break;
    case 7:
        vspar->psandy=32;
		vspar->psilt=34;
		vspar->pclay=34;
		vspar->k_solids=5.7709;
		vspar->c_solids=2260394;
		vspar->k_sat=0.000002845;
		vspar->psi_sat=-0.28893;
		vspar->porosity=0.4487;
		vspar->CH_b=8.32;
		break;
    case 8:
        vspar->psandy=10;
		vspar->psilt=56;
		vspar->pclay=34;
		vspar->k_solids=4.2564;
		vspar->c_solids=2326591;
		vspar->k_sat=0.000001311;
		vspar->psi_sat=-0.56104;
		vspar->porosity=0.4764;
		vspar->CH_b=8.32;
		break;
    case 9:
        vspar->psandy=52;
		vspar->psilt=6;
		vspar->pclay=42;
		vspar->k_solids=6.1728;
		vspar->c_solids=2242830;
		vspar->k_sat=0.000005756;
		vspar->psi_sat=-0.15805;
		vspar->porosity=0.4235;
		vspar->CH_b=9.59;
		break;
    case 10:
        vspar->psandy=6;
		vspar->psilt=47;
		vspar->pclay=47;
		vspar->k_solids=3.5856;
		vspar->c_solids=2355906;
		vspar->k_sat=0.000001139;
		vspar->psi_sat=-0.63299;
		vspar->porosity=0.4814;
		vspar->CH_b=10.38;
		break;
    case 11:
        vspar->psandy=22;
		vspar->psilt=20;
		vspar->pclay=58;
		vspar->k_solids=4.5370;
		vspar->c_solids=2314325;
		vspar->k_sat=0.000002;
		vspar->psi_sat=-0.39066;
		vspar->porosity=0.4613;
		vspar->CH_b=12.13;
		break;
	case 12:
        vspar->k_solids=0.2500;
		vspar->c_solids=2500000;
		vspar->k_sat=0.02;
		vspar->psi_sat=-0.1200;
		vspar->porosity=0.700;
		vspar->CH_b=4.00;
		break;
	case 13:
        vspar->k_solids=0.2500;
		vspar->c_solids=2500000;
		vspar->k_sat=0.15;
		vspar->psi_sat=-0.1200;
		vspar->porosity=0.900;
		vspar->CH_b=1.00;
		break;
	case 14:
        vspar->k_solids=0.2500;
		vspar->c_solids=2500000;
		vspar->k_sat=0.2;
		vspar->psi_sat=-0.0850;
		vspar->porosity=0.9500;
		vspar->CH_b=0.50;
		break;
	default:*/

}

/*初始化************************************************************************************************/
/* The variable of soil including temperature,water content,ice content,water head */
void svar_init(metarr_struct* metarr)
{
	FILE *fpp;
	int i;
	//int ok=1;
	//fpp = fopen("e:\\model\\@冻土BiomeBGC\\FSM-BGC\\ini\\初始化文件4.txt","r");

	//Uusing the data of run 2years result as the initization
	for(i = 0;i< N;i++)
	{
		svar[i].s_h = -0.05;
		svar[i].s_w = 0.15;
		svar[i].s_i = 0.05;
		svar[i].s_t = -0.31;
		//printf("soilt= %lf\n",svar[i].s_t);
	}
	/*fclose(fpp);
	if (ok && fscanf(fpp,"%*s%lf%lf%lf%lf",&svar[i].s_h,&svar[i].s_w,&svar[i].s_i,&svar[i].s_t)==EOF)
		{
			printf("Error reading s_init file, svar_init()\n");
			ok=0;
		}
	*/
}