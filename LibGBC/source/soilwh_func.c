#include "bgc.h"

/* soil water heat parameters function */
 
    //extern soilpar_struct vspar;      
	//extern soilvar_struct svar[N];  
   /* the volumetric heat capacity of the soil */
double Cv(double z_thw,double z_thi)
{ 
	double C=vspar->c_solids*(1-vspar->porosity)+C_w*z_thw+C_i*z_thi; return C;
}
  
   /* the thermal conductivity By Peter-Lidard et al.,1998 */

double Kersten(double thw,double thi)
{
	double ke;
	double Sr=(thw+thi)/vspar->porosity;

	if(thi>0)
		ke=(thw+thi)/vspar->porosity;
	else 
		if(Sr<0.05)
			ke=0.7*log10(Sr)+1.0;
		else
			ke=log10(Sr)+1.0;
	if(ke<0)
		ke=0;
	return ke;
}

double lamta_sat(double z_thw,double z_thi)
{   
	double xu=z_thw/(z_thw+z_thi)*vspar->porosity;
	double lamta=pow(lamta_i,vspar->porosity-xu)*pow(lamta_w,xu)*pow(vspar->k_solids,1-(vspar->porosity));return lamta;
} 

double lamta(double z_thw,double z_thi)
{   
	double lamtaa;
	double Sr=(z_thw+z_thi)/vspar->porosity;
	double rou_s=(1-vspar->porosity)*2700;
	double lamta_dry=(0.135*rou_s+64.7)/(2700-0.947*rou_s);
	if(Sr>10e-7)
		lamtaa=Kersten(z_thw,z_thi)*(lamta_sat(z_thw,z_thi)-lamta_dry)+lamta_dry;
	else
		lamtaa=lamta_dry;
	return lamtaa;
}

     /* hydraulic characteristic */

double K_water(double z_thw,double z_thi) 
{   
	double Ei=5/4*(vspar->k_sat*3600*100-3)*(vspar->k_sat*3600*100-3)+6;
	double K;
	if(z_thw<0.03||(vspar->porosity-z_thi)<0.05)
		K=0.0;
	else
		K=vspar->k_sat*pow(z_thw/(vspar->porosity-z_thi),3+2*vspar->CH_b)*pow(10,-Ei*z_thi);
	return K;
}

   /*  Unfrozen soil water content by frozen soil water head and ice content */
   double Watercontent1(double waterhead,double ice)
   {   
	   double wc,temp;
       temp=vspar->psi_sat*(1+Ck*ice)*(1+Ck*ice)/waterhead;
	   wc=vspar->porosity*pow(temp,1/vspar->CH_b);
	   return wc;
    }
   //Total water content by waterhead
   double Watercontent(double waterhead)
   {
	   double wc,temp;
	   temp=waterhead/vspar->psi_sat;
	   wc=vspar->porosity*pow(temp,-1/vspar->CH_b);
	   return wc;
   }

   // Maximum liquid water when the soil temperature is below the freezing point
   double Maxiwater(double T)
   {   
	   double maxiwater,temp;
	   if(T>=0)
	   {
           // printf("The temperature ought to be negative value");
	        maxiwater=0;
	   }
       else
       {
	        temp=Lf*T/(G_STD*(T+273.15)*vspar->psi_sat);
	        maxiwater=vspar->porosity*pow(temp,-1/vspar->CH_b);
	   }
	   return maxiwater;
   }
   //calculated Maximum supercooled soil liquid water by Koren et al.(1999)
   double Maxiwater1(double T,double icecontent)
   {
	   double maxiwater,temp;
	   temp=Lf*T/(G_STD*(T+273.15)*vspar->psi_sat*((1+Ck*icecontent)*(1+Ck*icecontent)));
	   maxiwater=vspar->porosity*pow(temp,-1/vspar->CH_b);
	   return maxiwater;
   }
   
   //calculated Waterhead only by  watercontent
   double Waterhead(double thw)
   {
	   double wh;
	   wh=vspar->psi_sat*pow(thw/vspar->porosity,-vspar->CH_b);
	   return wh;
   }
   //calculated Waterhead by waterconten and icecontent
   double Waterhead1(double thw,double thi)
   { 
	   double wh;
       wh=vspar->psi_sat*pow(thw/vspar->porosity,-vspar->CH_b)*(1+Ck*thi)*(1+Ck*thi);
	   return wh;
   }
   double Waterhead2(double thw,double thi)
   {
	   double wh;
	   wh = vspar->psi_sat * pow(thw/(vspar->porosity-thi),-vspar->CH_b);
	   return wh;
   }
   
   //differential coefficient water content to water head
   double derivationhead(double waterhead,double icecontent)
   { 
	   double temp;
       double wwh;
	   temp=(1+Ck*icecontent)*(1+Ck*icecontent);
	   wwh=vspar->porosity*pow(temp,1/vspar->CH_b)*(-1/vspar->CH_b)*pow(waterhead/vspar->psi_sat,-1/vspar->CH_b-1)/vspar->psi_sat;
	  // wwh=(vspar->porosity-icecontent)*(-1/vspar->CH_b)*pow(waterhead/vspar->psi_sat,-1/vspar->CH_b-1)/vspar->psi_sat;
	   return wwh;
   }
   //differential coefficient water content to icecontent
   double derivationice(double waterhead,double icecontent)
   { 
	   double wwi;
	   wwi=vspar->porosity*pow(vspar->psi_sat/waterhead,1/vspar->CH_b)*2*Ck*pow((1+Ck*icecontent),2/vspar->CH_b-1)/vspar->CH_b;//-1;
       //wwi=-pow(vspar->psi_sat/waterhead,1/vspar->CH_b);
	   return wwi;
   }
   //differential coefficient waterhead to water content
   double derivationw(double thw,double icecontent)
   {
	   double htow;
	   printf("derivationw %lf vspar->porosity= %lf -vspar->CH_b-1 = %lf, thw= %lf\n",vspar->CH_b,vspar->porosity,-vspar->CH_b-1,thw);
	   htow = vspar->psi_sat * (1 + Ck * icecontent) * (1 + Ck * icecontent) * (-1 / vspar->CH_b) * pow(thw/vspar->porosity,-vspar->CH_b-1)/vspar->porosity;
	   printf("htow= %lf, xxx=%lf£¬thw/vspar->porosity=%lf, yyy=%lf£¬zzz=%lf\n",htow,(-1 / vspar->CH_b),thw/vspar->porosity,pow(thw/vspar->porosity,-vspar->CH_b-1),1/vspar->porosity);
	   return htow;
   }

   //differential maximum watercontent to negative temperature T 
   double derivationt(double tem,double icecontent)
   {   
	   double temp, wtot;
	   temp=G_STD*Tk*(1+Ck*icecontent)*(1+Ck*icecontent)/Lf;
	   if((tem>=0) || (icecontent<=0))
		   wtot=0;
	   else
		   wtot=vspar->porosity*pow(temp,1/vspar->CH_b)*(-1/vspar->CH_b)*pow(tem/vspar->psi_sat,-1/vspar->CH_b-1)/vspar->psi_sat;
	   return wtot;
   }

   //differential hydraulic conductivity K to unfrozen water content thw
   double derivationk(double thw,double thi)
   {
       double Ei=5/4*(vspar->k_sat*3600*100-3)*(vspar->k_sat*3600*100-3)+6;
	   double KK;
	   KK=vspar->k_sat*pow(10,-Ei*thi)*(3+2*vspar->CH_b)*pow(thw/(vspar->porosity-thi),2+2*vspar->CH_b)/(vspar->porosity-thi);
	   return KK;

   }

   //caculate ice content by negative temperature and water content 
   double icecont(double tem,double thw)
   {
	   double temp,ice;
	   if(tem<0)
	   {
	   temp = tem*Lf*pow(thw/vspar->porosity,vspar->CH_b)/(vspar->psi_sat*G_STD*Tk);
       ice = (sqrt(temp)-1)/Ck;
	   }
	   else
		   ice =0;

	   if(ice>0)
		   return ice;
	   else
		   return 0;
   }
