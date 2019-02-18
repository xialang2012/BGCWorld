#include "bgc.h"

/* Saturation vapor pressure Tsate (Pa), as a function of temperature (deg C),Temperature in Ceisiul degree???? 
  are calculated from the eighth-order polynomial fits of Flatau et al. (1992) */
double vaporpressure_sat(double temp)
{ 
	double a0,a1,a2,a3,a4,a5,a6,a7,a8;
    double t=temp;
    double pvap;
	if(temp>=0)
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

  pvap=100*(a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*pow(t,5)+a6*pow(t,6)+a7*pow(t,7)+a8*(pow(t,8)));

  /* A very good approximation can usually be made using the August-Roche-Magnus formula 
  (usually called the Magnus or Magnus-Tetens approximation, though this is historically inaccurate[5]):
  
  pvap=610.94*exp(17.625*t/(t+243.04*t);   //pvap(T) is the equilibrium or saturation vapor pressure in Pa, which is a function of temperature; T is in Celsius.
 经比较两种方法结果近似相等：t>0,yy方案1比方案2值稍大1~2pa;t<0,方案1比方案2小10~20pa.没有进行实测值验证。
 */
  return pvap;
}

/* ground_heatflux.c for calcculating ground surface sensible heat  flux  W/m^2 */
int heatflux_suf(soilvar_struct ssvar[],double airvpd,double tair,double* vaporheatflux)
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
  double es;                              /*  vapor pressure of ground */
  double ea;                              /*  vapor pressure of air */
  double Hs;
  double psi=ssvar[0].s_h;
  int ok=1;
 // zm=zh*exp(0.24);
  ta=tair;                         /* the temperature of daytime */
  U = 0.5;                                 /*  wind speed */
  ts=ssvar[0].s_t;//+273.15;                        /* the temperature of ground surface */
  tk=ts+273.15;                          /* the unit transform to K from deg C */
  rou_a=1.292-(0.00428*ta);              /* air density at temperature ta */
  ra=((log((z-d)/zm)-psi_m)*(log((z-d)/zh)-psi_h))/(k*k*U);       /* By Campbell,1985 */
  // rs=exp(8.2-4.225*svar[0].s_w/spar->porosity);                   //Soil resistance is affected by soil type and soil water content and can be obtained from an empirical formula (Sun et al. 1998)
  rs=3.5*pow(vspar->porosity/ssvar[0].s_w,2.3)+33.5;                /* By linjiading et al.,1983 */
  // rs=-805+4140*(spar->porosity-svar[0].s_w);                      /* Camillo et al.,1986 */
  ea_sat=vaporpressure_sat(tk);
  es_sat=vaporpressure_sat(ts);
  Hs=exp(G_STD*psi/(Rw * tk));                     /* relative vapor pressure .R is universal gas constant(8.3143 J mole-1K-1 */
  if(Hs<1e-7) Hs=0;
  ea = ea_sat - airvpd;
  
  es = Hs*es_sat;
  /* evapor latent heat fluxes J/m2s = w/m2*/
 //  *latenthf = rou_a*CP*(es-ea)/(Hr*(ra+rs));
  
  /* counting air sensible heatflux (W/m2) */
  *vaporheatflux=rou_a*CP*(ts-ta)/(ra);       /* CP is specific heat of air. define constant */
 return (!ok);
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
int sunnetradiation(soilvar_struct ssvar[],double tair,double* Rn,double rad)
{ 
  //extern metvar_struct metvar;
  //extern soilvar_struct svar[10];
                                /* Net radiation of the accepted ground */
  double a=0.25;                          /* ground Surface albedo rato */
  double Rs;                              /* total sun shortwave radiation */
  double Rl;                              /* atmosphere long wave radiation */
  double Ls;                              /* ground longwave radiation */
  double Ea;                              /* atmosphere longwave radiation albedo */
  double Es;                              /* ground longwave radiation albedo */
 //double Ed;                              /* (mb)atmosphere vapour pressure */
 // double ea_sat;                           /* (pa)saturation vapor presssure */
  double ts,ta;                              /* deg C the temperature of the ground */
  double ws;                              /* the ground water content */
   int ok=1;
  ts=ssvar[0].s_t+273.15;
  ws=ssvar[0].s_w;
  ta=tair+273.15;                   /* ta is the temperature of air */
  // Ed= 1013.25;                            /* 1013.25mb=1 atmospere pressure(760mmHg) */
  //ea_sat=vaporpressure_sat(ta);    /* call vaporpressure_sat() */
  //Ed=(ea_sat-smetvar->vpd)/100;               /* unit pa->mb(毫巴=百巴=100pa） */
  //*Ea=1.24*pow(Ed/ta,1/7);
  Ea=0.96;
  Rl=Ea*SBC*pow(ta,4);

  Es=0.94+0.14*ws;                        /* By zheng xiuqing,2001 */
  Ls=Es*SBC*pow(ts,4);

  Rs=rad;                       /* (W/m2)  daylight average shortwave flux */

  *Rn=(1-a)*Rs+Rl-Ls;
  return (!ok);

}
