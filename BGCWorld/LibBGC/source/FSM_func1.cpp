#include "bgc.h"

//Thermal conductivity lamta_snow(Wm-1K-1)for snow is from Jordan(1991) at CLM
double lamta_snow(double rou_snow)
{
	double lamta;
	//lamta= lamta_a+(7.75e-5*rou_snow + 1.105e-6*rou_snow*rou_snow)*(lamta_i - lamta_a);
   // lamta=0.021+2.51*(rou_snow/rou_w)*(rou_snow/rou_w); //from Anderson E A,1976;郑秀清等，2002e
	lamta = 3.2217*pow(10.0,(-6))*(rou_snow*rou_snow);//Yen,1965
	return lamta;
}
/* The  volumetric heat capacities of snow(Jm-3k-1) (Verseghy,1991) */
double snow_c(double rou_snow)
{
	double sn_c;
	sn_c=1.9*pow(10.0,6)*rou_snow/rou_i;
	return sn_c;
}
/* the vapour density of snow derivative to Temperature */
double snow_vdensity_DT(double tem)
{
	double pv_sat,temp,tk;
	double snow_vden,snow_vden_dt;
	tk=tem+Tk;
	temp=52.57663-6790.4985/Tk-5.02808*log(Tk);
	pv_sat=(18/(R*Tk))*pow(FSM_e,temp);  //Campbell(1974)  saturated vapour density 

	temp=MW*Lf*tk/(R*Tk);
	snow_vden=pv_sat*pow(FSM_e,temp);
	snow_vden_dt=snow_vden*temp/tk;

	return snow_vden_dt;  
}
/* 水汽在雪层中的有效扩散系数(Effective diffusion coefficient),其值随温度而变化（Anderson E A,1976) */
double vapour_EDC(double tem,metvar_struct* metv)
{ 
	double De0, Nd,vedc; 
	De0=9*pow(15.0,-5);//De0 表示在0摄氏度，标准大汽压条件下，水汽在雪层中的有效扩散系数，取值9*15（-5）m2/s
	Nd=14.0 ; //经验常数
	vedc=De0*(101325/metv->pa)*pow(tem/Tk,Nd);

	return vedc;
}
/* soil water heat parameters function */
/* the volumetric heat capacity of the soil */
double Cv(double z_thw,double z_thi,soilpar_struct* vspar)
{ 
	double C = vspar->c_solids*(1-vspar->porosity)+C_w*z_thw+C_i*z_thi; return C;
}
/* the thermal conductivity By Peter-Lidard et al.,1998 */
double Kersten(double thw,double thi,soilpar_struct* vspar)
{
	double ke;
	double Sr=(thw+thi)/vspar->porosity;

	if(thi>0)
		ke=(thw+thi)/vspar->porosity;
	else 
		/*   if(Sr<0.05)
			 ke=0.7*log10(Sr)+1.0;
		 else
		     ke=log10(Sr)+1.0;*/

	   //2014年9月18日修正
	   if(Sr > 0.05)
       { 
		   ke = 0.7*log10(Sr)+1.0;
	       if( Sr > 0.1)
			  ke = log10(Sr)+1.0;
	   }
	   else
		   ke = 0;
	if(ke<0)
		ke=0;
	return ke;
}
double lamta_sat(double z_thw,double z_thi,soilpar_struct* vspar)
{   
	double xu=z_thw/(z_thw+z_thi)*vspar->porosity;
	double lamta=pow(lamta_i,vspar->porosity-xu)*pow(lamta_w,xu)*pow(vspar->k_solids,1-(vspar->porosity));return lamta;
}

double lamta(double z_thw,double z_thi,soilpar_struct* vspar)
{   
	double lamtaa;
	double Sr=(z_thw+z_thi)/vspar->porosity;
	double rou_s = (1 - vspar->porosity) * 2700;
	double lamta_dry=(0.135*rou_s+64.7)/(2700-0.947*rou_s);
	if(Sr>10e-7)
		lamtaa=Kersten(z_thw,z_thi,vspar)*(lamta_sat(z_thw,z_thi,vspar)-lamta_dry)+lamta_dry;
	else
		lamtaa=lamta_dry;
	return lamtaa;
}
/* hydraulic characteristic */
double K_water(double z_thw,double z_thi,soilpar_struct* vspar) 
{   
	double Ei=5/4*(vspar->k_sat*3600*100-3)*(vspar->k_sat*3600*100-3)+6;
	double K;
	if(z_thw<0.03||(vspar->porosity-z_thi)<0.05)
		K=0.0;
	else
		K=vspar->k_sat*pow(z_thw/(vspar->porosity-z_thi),3+2*vspar->CH_b)*pow(10,-Ei*z_thi);
	return K;
}
// Maximum liquid water when the soil temperature is below the freezing point
double Maxiwater(double T,soilpar_struct* vspar)
{   
	double maxiwater,temp;
	if(T>=0)
	{
		// printf("The temperature ought to be negative value");
		maxiwater=vspar->porosity ;
	}
	else
	{
		temp=Lf*T/(FSM_g*(T+273.15)*vspar->psi_sat);
		maxiwater=vspar->porosity*pow(temp,-1/vspar->CH_b);
	}
	return maxiwater;
}
//calculated Waterhead only by  watercontent
double Waterhead(double thw,soilpar_struct* vspar)
{
	double wh;
	wh = vspar->psi_sat * pow(thw/vspar->porosity,-vspar->CH_b);	
	return wh;
}
//calculated Waterhead by waterconten and icecontent
double Waterhead1(double thw,double thi,soilpar_struct* vspar)
{ 
	double wh;
	wh = vspar->psi_sat * pow(thw/vspar->porosity,-vspar->CH_b)*(1+Ck*thi)*(1+Ck*thi);
	return wh;
}
/**/double Waterhead2(double thw,double thi,soilpar_struct* vspar)
{
	double wh;
	wh = vspar->psi_sat*pow(thw/(vspar->porosity-thi),-vspar->CH_b);
	return wh;
}
//differential coefficient water content to water head
double derivationhead(double waterhead,double icecontent,soilpar_struct* vspar)
{ 
	double temp;
	double wwh;
	temp=(1+Ck*icecontent)*(1+Ck*icecontent);
	wwh=vspar->porosity*pow(temp,1/vspar->CH_b)*(-1/vspar->CH_b)*pow(waterhead/vspar->psi_sat,-1/vspar->CH_b-1)/vspar->psi_sat;
	// wwh=(vspar->porosity-icecontent)*(-1/vspar->CH_b)*pow(waterhead/vspar->psi_sat,-1/vspar->CH_b-1)/vspar->psi_sat;
	return wwh;
}
//differential coefficient water content to icecontent
double derivationice(double waterhead,double icecontent,soilpar_struct* vspar)
{ 
	double wwi;
	wwi=vspar->porosity*pow(vspar->psi_sat/waterhead,1/vspar->CH_b)*2*Ck*pow((1+Ck*icecontent),2/vspar->CH_b-1)/vspar->CH_b;//-1;
	//wwi=-pow(vspar->psi_sat/waterhead,1/vspar->CH_b);
	return wwi;
}
//differential coefficient waterhead to water content
double derivationw(double thw,double icecontent,soilpar_struct* vspar)
{
	double htow;
	htow=vspar->psi_sat*(1+Ck*icecontent)*(1+Ck*icecontent)*(-1/vspar->CH_b)*pow(thw/vspar->porosity,-vspar->CH_b-1)/vspar->porosity;
	return htow;
}
//differential maximum watercontent to negative temperature T 
double derivationt(double tem,double icecontent,soilpar_struct* vspar)
{   
	double temp, wtot;
	temp=FSM_g * Tk * (1 + Ck * icecontent) * (1 + Ck * icecontent) / Lf;
	if((tem>=0) || (icecontent<=0))
		wtot=0;
	else
		wtot=vspar->porosity*pow(temp,1/vspar->CH_b)*(-1/vspar->CH_b)*pow(tem/vspar->psi_sat,-1/vspar->CH_b-1)/vspar->psi_sat;
	return wtot;
}
//differential hydraulic conductivity K to unfrozen water content thw
double derivationk(double thw,double thi,soilpar_struct* vspar)
{
	double Ei=5/4*(vspar->k_sat*3600*100-3)*(vspar->k_sat*3600*100-3)+6;
	double KK;
	KK=vspar->k_sat*pow(10,-Ei*thi)*(3+2*vspar->CH_b)*pow(thw/(vspar->porosity-thi),2+2*vspar->CH_b)/(vspar->porosity-thi);
	return KK;
}


