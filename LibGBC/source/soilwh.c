#include "bgc.h"
/* here calculated water head,water content,ice content,and temperature in every layer 
   Using Shangsonghao method 2011.11.8 succeed 
*/

int soilwh(const siteconst_struct* sitec,soilvar_struct svar[],soilvar_struct soil_sum[],metvar_struct* metvar,wflux_struct* wf,FILE* fp)
{
  	/* function soilwh.c  variables  */
	
	soilvar_struct s0[N];                                             /* previous iterative value  */
	soilvar_struct s1[N];                                             /* value of beginning timestep */
	int ok = 1;
	int delta_t=60*30;                                                /* s time step length */
	int flag,i;
	int nn,t_index=0,temp; //Time 
	
	double evap,latenthf;
	double sensiblehf,groundhf,netradiation,precip,groundwf;
	
	double w,rr;
    double m,q1,qcm=0,T_amp=0;                                //qcm is the maximun infiltration capacity in precipitation day ;  tamp is the daily temperature amplitude
	double runoff=0;
	double lamtaii[N],Kii[N],Dii[N],derivationh[N],derivationi[N];            //note parameter:conduction of heat,diffuse of water,conductivity of water
	double K_s[N],K_n[N],D_s[N],D_n[N],lamta_s[N],lamta_n[N];                        //note boundary parameters
    double C1[N],Cvv[N],Ce[N],lamtas[N],lamtan[N],Ue[N],R3[N],R1n[N],R1s[N];   //solved temperature required varible
	
	double wm1,wm2,wm3,wm;
   
	double E[N],F[N],G[N],H[N];                                        //Water function coefficient
	double ET[N],FT[N],GT[N],HT[N];                                    //heat function coefficient 
	double dz[N]={0.02,0.08,0.10,0.20,0.20,0.20,0.20,0.30,0.30,0.40,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4}; //layers length(thickness)
	double dz_bgc[N]={2,10,20,40,60,80,100,130,160,200,240,280,320,360,400,440,480,520,560,600}; //layers depth(thickness)
	double Tn=15.0;
	double tair;
	double efrad;
	int N_bgc = 0;
	/*添加调整系数，计算每一层非冻土碳含量 2014-08-05*/
	double frost_soilc[N],nofrost_soilc[N],C_sum_ratio[N],Cdensity[N],Soil_ratio[N];
	//类型：森林k = -0.025
	//Cdensity = exp(-0.025);
	for(int i=0;i<N;i++)
	{
		C_sum_ratio[i] = Cdensity[i] * dz[i];//每一层的碳含量占总碳含量的比例，标准化为0~1之间
		if(sitec->soil_depth < dz_bgc[i])
			N_bgc = dz_bgc[i];
	}
	for(int i=0;i<N_bgc;i++)
	{
		//随着深度的变化，计算每一层中土壤碳密度占总碳库的比例
		C_sum_ratio[i] += exp(-0.025 * i);//需要求曲线积分之和
	}
	for(int i = 0;i < N_bgc; i++)
	{
		//随着深度的变化，计算每一层中土壤碳密度占总碳库的比例
		Cdensity[i] = exp(-0.025 * i) / C_sum_ratio[i];//需要求曲线积分之和
	}

	w=0.5;                                                              //relaxation factor
	rr=C_i/(Lf*rou_i);                                                  //the unit conversion factor(K)
	
	/* For counting every layer average temperaturre */
	for(i=0;i<N;i++)
	{
		soil_sum[i].s_t=0;
		soil_sum[i].s_tw=0;
		soil_sum[i].s_w=0;
		soil_sum[i].s_i=0;
		soil_sum[i].s_h=0;
	}
	precip = wf->snoww_to_soilw + wf->prcp_to_soilw + wf->canopyw_to_soilw - wf->soilw_evap - wf->soilw_trans;
	/* kg/m2 -> m */
	precip = precip/rou_w;
	latenthf=wf->soilw_evap*Lvap/(24*3600);
	if(metvar->prcp > 0.0)
	{
        q1=(vspar->k_sat*(vspar->psi_sat-svar[1].s_h)/(dz[0]+dz[1])/2)*24*3600;  //m/s -> m Top soil layer infiltration
        qcm=(dz[0]*(vspar->porosity-svar[0].s_w)+q1);  //+dz[1]*(vspar->porosity-svar[1].s_w)                //m 
        if(precip>qcm)
		{
			runoff=precip-qcm;
		    precip=(metvar->prcp-runoff)/(24*3600);              // m/s
		}
		else
			precip=precip/(24*3600);                      // m/s
	}
	else
		precip = precip/(24*3600);
	T_amp=(metvar->tmax-metvar->tmin)/2;  // the amplitude of the daliy temperature 		            
		   
	/* do loop timestep 30 min  within every day(24hour) */
	for(t_index = 0; t_index < 24*3600; t_index = t_index + delta_t)
	{
		nn=0; 
		/*
		        interpolated air temperature from the measured max and min values by a time step of 30 min
		        date:2011.09.10 modified the temperature of air every time step using cosine: minimun temperature 
			    at 03:00 hr and maximun at 15:00hr daily.before 03:00hr using the maximun temperature of the previous day ; 
			    after 15:00hr using the minimun pemperature of the next day in order to smoothly vary between days
		        */
				
        temp=t_index/3600; 
		  				
        tair=metvar->tavg+T_amp*cos(2*3.14159*(temp-15)/24);
		//given daytime swavgfd=metvar.swavgfd and night swavgfd=0
		temp=temp*3600;
		if(temp<7*3600||temp>7*3600+metvar->dayl)
			efrad = 0;
		else
			efrad = metvar->swtrans;
        /* call heatflux_suf() to calculate ground surface sensible and latent heat flux */       
        heatflux_suf(svar,metvar->vpd,tair,&sensiblehf);
	
		/* call sunnettradiation90 to calculate shortwave and longwave radation */
        sunnetradiation(svar,tair,&netradiation,efrad);

        /* entering ground surface heat flux (W/m2) */
		groundhf=netradiation-sensiblehf-latenthf;

	    /* entering ground surface water flux (m/s) */        
		groundwf=precip;
        for(i=0;i<N;i++)
		{  
			 s1[i]=svar[i];
		}
	    do   
	    {	
	    	flag=0;
	    	for(i=0;i<N;i++)
	        {  
		        s0[i]=svar[i];
		        Kii[i] = K_water(svar[i].s_w,svar[i].s_i);
				//vspar->CH_b = 6;
				
			    Dii[i] = Kii[i] * derivationw(svar[i].s_w,svar[i].s_i);
		        lamtaii[i] = lamta(svar[i].s_w,svar[i].s_i);
       
		        derivationh[i]=derivationhead(svar[i].s_h,svar[i].s_i);
		        derivationi[i]=derivationice(svar[i].s_h,svar[i].s_i);
            }
	        for(i=1;i<N-1;i++)
	        {
		        K_n[i]=(Kii[i]+Kii[i-1])/2;//sqrt(Kii[i]*Kii[i-1]);//2/(1/Kii[i]+1/Kii[i-1]);
                K_s[i]=(Kii[i]+Kii[i+1])/2;//sqrt(Kii[i+1]*Kii[i]);//2/(1/Kii[i]+1/Kii[i+1]);
				D_n[i]=(Dii[i]+Dii[i-1])/2;
                D_s[i]=(Dii[i]+Dii[i+1])/2;
		        lamta_n[i]=(lamtaii[i]+lamtaii[i-1])/2;//2/(1/lamtaii[i]+1/lamtaii[i-1]);
	            lamta_s[i]=(lamtaii[i]+lamtaii[i+1])/2;//2/(1/lamtaii[i]+1/lamtaii[i+1]);
		    }
            K_s[0]=K_n[1];//2/(1/Kii[0]+1/Kii[1]);K_n[0]=0;
			D_s[0]=D_n[1];
            lamta_s[0]=2/(1/lamtaii[0]+1/lamtaii[1]);lamta_n[0]=lamta_s[0];
            K_n[N-1]=K_s[N-2];//2/(1/Kii[N]+1/Kii[N-1]);K_s[N]=0;
            D_n[N-1]=D_s[N-2];
            lamta_n[N-1]=2/(1/lamtaii[N-1]+1/lamtaii[N-2]);lamta_s[N-1]=lamta_n[N-1];
					  
			//sloved temperature T[i] by equation (13)
            //coefficient of equation
		    
		    ET[0]=0;
			FT[0]=Cv(svar[0].s_w,svar[0].s_i)/delta_t+2*lamta_s[0]/(dz[0]*(dz[0]+dz[1]));
		    GT[0]=-2*lamta_s[0]/(dz[0]*(dz[0]+dz[1]));
		    HT[0]=Cv(svar[0].s_w,svar[0].s_i)*s1[0].s_t/delta_t+groundhf/dz[0];//+Lf*rou_i*(svar[0].s_i-s1[0].s_i)/delta_t;
		    /*
			FT[0]=lamta_s[0]/dz[0]+rou_a*CP/ra;
			GT[0]=-lamta_s[0]/dz[0];
			HT[0]=netradiation-latenthf2+rou_a*CP*metvar.tavg/ra;
			//HT[0]=rou_a*CP*metvar.tavg/ra;
			*/
		    for(i=1;i<N;i++)
		    {
			    //solved temperature required varible
			    C1[i]=Lf*rou_w*derivationt(svar[i].s_t,svar[i].s_i);
                Cvv[i]=Cv(svar[i].s_w,svar[i].s_i);
			    Ce[i]=Cvv[i]+C1[i];
			    lamtas[i]=lamta_s[i]+D_s[i]*C1[i];
			    lamtan[i]=lamta_n[i]+D_n[i]*C1[i];
			    Ue[i]=C1[i]*derivationk(svar[i].s_w,svar[i].s_i);

			    R3[i]=delta_t/dz[i];
			    R1n[i]=2*R3[i]/(dz[i]+dz[i-1]);
			    R1s[i]=2*R3[i]/(dz[i]+dz[i+1]);

			    ET[i]=-R1n[i]*lamtan[i]-R3[i]*Ue[i];
			    GT[i]=-R1s[i]*lamtas[i]+R3[i]*Ue[i];
			    FT[i]=Ce[i]-ET[i]-GT[i];
			    HT[i]=Ce[i]*s1[i].s_t;
			}
			
			 /*
			 HT[0]=HT[0]-ET[0]*T0;
			 ET[0]=0;
			 HT[N-1]=HT[N-1]-GT[N-1]*Tn;
			 GT[N-1]=0;
             */
			//Qn=0 heat fluxes equal to zero at lower boundary 
            /* ET[N-1]=-2*lamta_s[N-1]/(dz[N-1]*(dz[N-1]+dz[N-2]));
			FT[N-1]=Cv(svar[N-1].s_w,svar[N-1].s_i)/delta_t+2*lamta_s[N-1]/(dz[N-1]*(dz[N-1]+dz[N-2]));
			GT[N-1]=0;
			HT[N-1]=Cv(svar[N-1].s_w,svar[N-1].s_i)*s1[N-1].s_t/delta_t;//+Lf*rou_i*(svar[N-1].s_i-s1[N-1].s_i);
            */
		    //heat;counting upper triangle 
	        for(i=1;i<N;i++)
	        {
			    FT[i] = FT[i]-ET[i]*GT[i-1]/FT[i-1];
	            HT[i] = HT[i]-HT[i-1]*ET[i]/FT[i-1]; 
		    }
	        //heat;counting  temprature 
	        svar[N-1].s_t = 10;//HT[N-1]/FT[N-1];
	        for(i=N-2;i>=0;i--)
	        { 
	            svar[i].s_t = (HT[i]-GT[i]*svar[i+1].s_t)/FT[i];
				if(fabs(svar[i].s_t-s0[i].s_t)>0.01) 
				{
					 flag=1;
				}
	        }
            //solved  waterhead by equation (12)   
	    	//water head  implicit finite difference equation coefficient 
            E[0]=0;
		    F[0]=1+2*D_s[0]*delta_t/(dz[0]*(dz[0]+dz[1]));
		    G[0]=-2*D_s[0]*delta_t/(dz[0]*(dz[0]+dz[1]));
		    H[0]=s1[0].s_w+(groundwf-K_s[0])*delta_t/dz[0];//-rou_i/rou_w*(svar[0].s_i-s1[0].s_i);//-derivationh[0]*(WHEAD0[0]-WHEAD1[0])/delta_t;
     
		    for(i=1;i<N-1;i++)
	        {
	            E[i]=-2*D_n[i]*delta_t/(dz[i]*(dz[i]+dz[i-1]));
                G[i]=-2*D_s[i]*delta_t/(dz[i]*(dz[i]+dz[i+1]));
    	        F[i]=1-E[i]-G[i];
                H[i]=s1[i].s_w-2*(K_s[i]-K_n[i])*delta_t/(dz[i]+dz[i+1]);//-rou_i/rou_w*(svar[i].s_i-s1[i].s_i);//-derivationh[i]*(WHEAD0[i]-WHEAD1[i])/delta_t;
		    }
		 
		    E[N-1]=-2*D_n[N-1]*delta_t/(dz[N-1]*(dz[N-1]+dz[N-2]));
		    F[N-1]=1-E[N-1]; 
		    G[N-1]=0;
		    H[N-1]=s1[N-1].s_w-delta_t*(Kii[N-1]-K_n[N-1])/dz[N-1];//-rou_i/rou_w*(svar[N-1].s_i-s1[N-1].s_i);//-derivationh[N-1]*(WHEAD0[N-1]-WHEAD1[N-1])/delta_t;
	  
	        //counting upper triangle FOR Water content
	        for(i=1;i<N;i++)
	        { 
		        F[i]=F[i]-E[i]*G[i-1]/F[i-1];
	            H[i]=H[i]-H[i-1]*E[i]/F[i-1];
		    }
	        //counting  water content stu[i] 
	        svar[N-1].s_tw=H[N-1]/F[N-1];
			
	        for(i=N-2;i>=0;i--)
	        {
			    svar[i].s_tw=(H[i]-G[i]*svar[i+1].s_tw)/F[i];
		        if(fabs(svar[i].s_tw-s0[i].s_tw)>0.0001) 
				{
					flag=1;
				}
			}
			if(fabs(svar[0].s_t-s1[0].s_t) > 2)
			{
				svar[0].s_t=(svar[0].s_t+s1[0].s_t)*0.5;
			}
           	nn++;
			if(nn>10) break;
         // for(i=0;i<5;i++)
			//  fprintf(fp,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",nn,svar[i].s_t,svar[i].s_tw,svar[i].s_w,svar[i].s_i,svar[i].s_h);
					
        }while(flag && nn<6); //end loop one timestep
	    //update ice content,liquid water content
		//Update ice content(THI[i]) by equation (14)//||(svar[i].s_t>=-0.01)&&svar[i].s_i>=0.0001))
		       
	    for(i=0;i<N;i++)
	    {
	        if(svar[i].s_t<0)
		    {
			    svar[i].s_tt=(vspar->psi_sat*G_STD*Tk/Lf)*pow(svar[i].s_w/vspar->porosity,-vspar->CH_b)*(1+Ck*svar[i].s_i)*(1+Ck*svar[i].s_i);
                //  svar[i].s_t=273.16*(svar[i].s_h/(Lf/g-svar[i].s_h));
			    rr=Cv(svar[i].s_w,svar[i].s_i)/(Lf*rou_w);
                wm2=(svar[i].s_tt-svar[i].s_t)*rr;
											
			    wm1=Maxiwater1(svar[i].s_t,svar[i].s_i);
			    wm3=svar[i].s_tw-wm1;
			    wm=(fabs(wm2)>fabs(wm3))?wm3:wm2;
                if(wm3>0)
			    {
                    svar[i].s_w=svar[i].s_tw-wm;
				    svar[i].s_i=svar[i].s_i+wm*rou_w/rou_i;
			    }
		    }
		    else
		    	if(svar[i].s_i>0)
		    	{   
		    		rr=Cv(svar[i].s_w,svar[i].s_i)/(Lf*rou_w);
		    		wm2=svar[i].s_t*rr;
		    		if(svar[i].s_i>wm2)
		    	    {
                        svar[i].s_w=svar[i].s_w+wm2*rou_i/rou_w;
			    		svar[i].s_i=svar[i].s_i-wm2;
			    	}
		    		else
			    	{
				    	svar[i].s_w=svar[i].s_w+svar[i].s_i*rou_i/rou_w;
				        svar[i].s_i=0;
				    }
		        }
	    		else
		    		svar[i].s_w=svar[i].s_tw;
					
            if(svar[i].s_w<0.035) svar[i].s_w=0.035;
	    	if(svar[i].s_w>vspar->porosity) svar[i].s_w=vspar->porosity;
	    	svar[i].s_tw=svar[i].s_w + svar[i].s_i * rou_i/rou_w;
    	}
       //countint waterhead
	    for(i=0;i<N;i++)
	    {
	    	svar[i].s_h = Waterhead2(svar[i].s_w,svar[i].s_i);//+rou_i*svar[i].s_i/rou_w);
	    	soil_sum[i].s_t += svar[i].s_t;
		    soil_sum[i].s_tw += svar[i].s_tw;
		    soil_sum[i].s_w += svar[i].s_w;
		    soil_sum[i].s_i += svar[i].s_i;
		    soil_sum[i].s_h += svar[i].s_h;
	    }

	} //for end all timestep daily 
	/* every layer value summary translate to average */
	temp = 3600*24/delta_t;
	for(i=0;i<N;i++)
	{
        soil_sum[i].s_t /= temp;
		soil_sum[i].s_tw /= temp;
		soil_sum[i].s_w /= temp;
		soil_sum[i].s_i /= temp;
		soil_sum[i].s_h /= temp;
	}
	/*
	1 每一层计算一个含冰量，含水量，温度
	2 推算出每一层冻土的含量，和未冻的含量。冻土的含量，根据每一层的含冰量与所有水分的比值，估算出冻土占所有土的比例
	*/
	for(i=0;i<N_bgc;i++)
	{
		//计算每一层中，冻土中的碳所占的比例；非冻土中碳所占的比例；
		frost_soilc[i] = (soil_sum[i].s_tw - soil_sum[i].s_w) * Cdensity[i];
		nofrost_soilc[i] = soil_sum[i].s_w * Cdensity[i];		
	}
	




	return (!ok);
 
}//end function