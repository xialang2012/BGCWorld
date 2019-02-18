#include "bgc.h"
/* here calculated water head, water content, ice content, and temperature in every layer 
   Using Shangsonghao method 2011.11.8 succeed  . amended at 20140802 , mainly alter snow cover
*/
int FSM_core(metvar_struct* metv,int metday,pmet_struct *pmetvar,soilpar_struct* vspar,soilvar_struct soil_sum[],epvar_struct* epv,siteconst_struct* sitec,wstate_struct *ws)
{ 	   	
	int i = 0,delta_t = 1800,N_HalfHour = 48;                                                /* s time step length 30分钟1800秒*/
	int flag;
    int nn,t_index=0, soiltop = 5;//Time, 土壤第一层的索引
	int snl = 5;                     //snl:the number of snow layers
	int snow_flag = 0,newsnow_flag;
	//int temp_i = 0,temp_i1 = 0,temp_i2 = 0,temp_i3 = 0,temp_i4 = 0;
	int rain_start,rain_duration;
	double MinWater = 0.03;
	double temp;
	double evap,latenthf;
	double sensiblehf,groundhf,netradiation,precip,temp_precip,groundwf;
	double sumqn = 0.0,sumevap = 0.0;
	double ra;
    double m,qcm=0,T_amp=0;                                //qcm is the maximun infiltration capacity in precipitation day ;  tamp is the daily temperature amplitude
	double rou_a,runoff = 0.0;
	double lamtaii[N],Kii[N],Dii[N],derivationh[N],derivationi[N];            //note parameter:conduction of heat,diffuse of water,conductivity of water
	double K_s[N],K_n[N],D_s[N],D_n[N],lamta_s[N],lamta_n[N],sink[N];                        //note boundary parameters
    double C1[N],Cvv[N],Ce[N],lamtas[N],lamtan[N],Ue[N],R2[N],R3[N],R1n[N],R1s[N];   //solved temperature required varible
	
	double Tsum[N],Tlay[N],sumnr,sumshf,sumlhf,sumghf;
	double wm1;
	//snow_depth0:新雪厚度；snow_depth1：老雪厚度；
    double rou_snow_new = 0.0,rou_snow_old = 0.0,zsnow = 0,snow_depth0,snow_depth1 = 0,snow_density_delta1,snow_density_delta2,snow_density_delta3,temp0,temp1;                                        //the denisity of snow , the number of snow layers and  total snow thickness (depth)
	double E[N],F[N],G[N],H[N];                                        //Water function coefficient
	double ET[N],FT[N],GT[N],HT[N];                                    //heat function coefficient 
	double dz[N]={0.0,0.0,0.0,0.0,0.0,0.02,0.08,0.10,0.20,0.20,0.20,0.20,0.30,0.30,0.40,1.2,1.2,1.2,1.2,1.2}; // 前5层为雪盖次预留，soil layers length(thickness total 4m)
	double dz_sum[N]={0.0,0.0,0.0,0.0,0.0,0.02,0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.6,2.0,3.2,4.4,5.6,6.8,8.0};
	double HI[5];  //The excess or deficit of energy Hi(w/m2)
	double w_i[5],w_l[5],rou_sno[5],C_sno[5];
	double w_i1[5],w_l1[5];//w_i[i] and w_l[i] are the mass of ice and liquid water (kg/m2) in each snow/soil layer,respectively.
	double shf_dt,lhf_dt,Rn_dt,groundhf_dt;
	double melthf = 0,meltwf = 0,temp_w,melthigh = 0.0;//melthf:雪融化通量 kg/m2; meltwf:雪融化热通量W/m2
	double ps; //the weight of snow(expressed in terms of centimeters of water equivalent) above the layer of snow.
	soilvar_struct s0[N];                                             /* previous iterative value  */
	soilvar_struct s1[N];                                             /* value of beginning timestep */
	double FSM_psi = 0,FSM_water = 0,FSM_ice = 0,FSM_soilT = 0,FSM_snowT = 0,FSM_depth = 0,FSM_snowW=0;
	double FSM_lrad_down=0,FSM_lrad_up=0,bgc_lrad_input=0;
	//pa = 101300 * pow((1.0 - (0.0065 * 579.27)/288.15),5.25589);//plot1,elevation = 579.27m
	//降雨的日插值。算法：随机产生降雨开始的时间和降雨持续的时间。
	//rain_start = (int)(rand()%46); //0 - 47
	//rain_duration = (int)(rand()%24) + 1;  //1 - 48
	
	rain_start = 12; //0 - 47
	rain_duration = 10;  //1 - 48

	for(i = 0;i < N;i++)
	{
		Tsum[i] = 0;
		sink[i] = 0;
		//dz_sum[i] = 0;
		soil_sum[i].s_t = 0;
		soil_sum[i].s_tw = 0;
		soil_sum[i].s_w = 0;
		soil_sum[i].s_i = 0;
		soil_sum[i].s_h = 0;
	}
	/*长波辐射的计算 每半个小时计算一次，最后求全天的总平均，与输入的长波辐射进行对比
	FSM_lrad_down += 0.96 * SBC * pow(pmetvar->ta + 273.15,4);*/
	bgc_lrad_input = metv->lrad;


	/* Whether precipitation for snow */
	//判断是否有降雪
	if((metv->prcp> CRIT_PREC)&&((metv->tavg)<CRIT_PREC))
		newsnow_flag = 1;
	else
		newsnow_flag = 0;
	//rou_sno is the bulk density of newly fallen snow(kg.m3)(Anderson 1976)
	//TODO 雪的密度赋值有问题，没有考虑老雪。
	if(newsnow_flag)
	{
		if(metv->tavg > 2.0)
			rou_snow_new = 50 + 1.7 * pow(17,1.5);
		else if((metv->tavg <= 2.0) & (metv->tavg > -15.0))
			rou_snow_new = 50 + 1.7 * pow(metv->tavg +15,1.5);
		else 
			rou_snow_new = 50;
		/* when only given daily water equivalent of snow ,please notes: units conversion. cm to kg/m2 to m/m2*/
		snow_depth0 = 0.01 * metv->prcp * rou_w /rou_snow_new;//m
	}
	else
		snow_depth0 = 0;
	/*老雪的密度等于昨天雪的密度
	if(!snow_depth1)
	{
	for(int i=0;i<5;i++)
	rou_sno[i] = rou_snow_new;
	rou_snow_old = 0;
	}
	else
	{
	for(int i=0;i<5;i++)
	rou_sno[i] = rou_snow_old;
	}*/
	for(i = 0;i<5;i++)
		rou_sno[i] = rou_snow_new;
	//according to precipitation(cm) of everyday,transform to water fluex(m/s) of the ground within 2 hours
	//TODO 需要将降雨插值到逐半小时
	if(!snow_flag && metv->prcp>CRIT_PREC)
		precip = metv->prcp / (rain_duration * 180000);                      // m/s 厘米转化为米
	else
		precip = 0;
	
	/******************************For cycle****************************************/
	/* do loop timestep 30 min  within every day(24hour) */
	for(t_index = 0;t_index < 24*3600; t_index = t_index + delta_t)
	{
		nn = 0; 
		//降雨插值到半小时
		if(rain_start + rain_duration <= 48)
		{
			if(((metv->prcp > CRIT_PREC) && (metv->tavg > CRIT_PREC)) && (t_index >= rain_start * 1800) && (t_index < (rain_start + rain_duration) * 1800))//rain timeduring=3 hr
				temp_precip = precip;
			else
				temp_precip = 0;
		}
		else
		{
			if((t_index >= rain_start * 1800) && (t_index < (rain_start + rain_duration) * 1800) && ((metv->prcp > CRIT_PREC) && (metv->tavg > CRIT_PREC)))//rain timeduring=3 hr
				temp_precip = precip;
			else if(((metv->prcp > CRIT_PREC) && (metv->tavg > CRIT_PREC)) && (t_index >= 0 ) && (t_index < (rain_start + rain_duration - 48) * 1800))//rain timeduring=3 hr
				temp_precip = precip;
			else
				temp_precip = 0;
		}
		/*[BGC中删掉
		sumnr=0;
		sumshf=0;
		sumlhf=0;
		sumghf=0;
		temp_i1 = temp_i2 = temp_i3 = temp_i4 = 0;
		BGC中删掉]*/
		melthigh = 0.0;//每半小时重新初始化融雪量
		
		/*
		interpolated air temperature from the measured max and min values by a time step of 30 min
		date:2011.09.10 modified the temperature of air every time step using cosine: minimun temperature 
		at 03:00 hr and maximun at 15:00hr daily.before 03:00hr using the maximun temperature of the previous day ; 
		after 15:00hr using the minimun pemperature of the next day in order to smoothly vary between days
		*/
		//温度插值到逐半小时，T_amp表示温度的振幅
		T_amp = metv->tmax - metv->tmin;  // the amplitude of the daliy temperature 
		temp = t_index / 1800; 
		pmetvar->ta = metv->tavg + T_amp * 0.5 * cos(2 * 3.14159 * (temp - 30) / 48.0);	
		rou_a = 1.292 - (0.00428 * pmetvar->ta);              /* air density at temperature ta */
		//given daytime swavgfd=metv->swavgfd and night swavgfd=0
		//TODO 需要测试实际的日长
		if(t_index < 6*3600 || (t_index >= 6*3600 + metv->dayl))
			pmetvar->irad = 0;
		else if(t_index >= 6* 3600 && t_index < (6 * 3600 + ((int)(metv->dayl/1800)) * 1800))
			pmetvar->irad = metv->swavgfd;
		else
			pmetvar->irad = metv->swavgfd * (((int)metv->dayl % 1800)/1800.0);
		zsnow = snow_depth0 + snow_depth1;   //snow depth m
		// The number of snow layers and the thickness of each layer is a function of snow depth zno(m) as follows..
		// zsnow：单位m
		if(zsnow < 0.005)
		{
			snl = 5; 
			dz[4] = 0;
			if(snow_flag)
				meltwf = (0.01 * metv->prcp + snow_depth1 * rou_snow_old / rou_w) / delta_t;
			snow_flag = 0;
			snow_depth0 = 0;
			snow_depth1 = 0;
		}
		else if ((zsnow >= 0.005) && (zsnow <= 0.03))
		{
			snow_flag=1;
			snl=4;
			dz[4]=zsnow;
		}
		else if ((zsnow > 0.03) && (zsnow <= 0.04))
		{
			snow_flag=1;
			snl=3;
			dz[3]=zsnow/2;
			dz[4]=dz[3];
		}
		else if ((zsnow>0.04) && (zsnow<=0.07))
		{
			snow_flag=1;
			snl=3;
			dz[3]=0.02;
			dz[4]=zsnow-dz[3];
		}
		else if ((zsnow>0.07) && (zsnow<=0.12))
		{
			snow_flag=1;
			snl=2;
			dz[2]=0.02;
			dz[3]=(zsnow-0.02)/2;
			dz[4]=dz[3];
		}
		else if ((zsnow>0.12) && (zsnow<=0.18))
		{
			snow_flag=1;
			snl=2;
			dz[2]=0.02;
			dz[3]=0.05;
			dz[4]=zsnow-dz[2]-dz[3];
		}
		else if ((zsnow>0.18) && (zsnow<=0.29))
		{
			snow_flag=1;
			snl=1;
			dz[1]=0.02;
			dz[2]=0.05;
			dz[3]=(zsnow-dz[2]-dz[1])/2;
			dz[4]=dz[3];
		}
		else if ((zsnow>0.29) && (zsnow<=0.41))
		{
			snow_flag=1;
			snl=1;
			dz[1]=0.02;
			dz[2]=0.05;
			dz[3]=0.11;
			dz[4]=zsnow-dz[3]-dz[2]-dz[1];
		}
		else if ((zsnow>0.41) && (zsnow<=0.64))
		{
			snow_flag=1;
			snl=0;
			dz[0]=0.02;
			dz[1]=0.05;
			dz[2]=0.11;
			dz[3]=(zsnow-dz[0]-dz[2]-dz[1])/2;
			dz[4]=dz[3];
		}
		else if ((zsnow>0.64))
		{
			snow_flag=1;
			snl = 0;
			dz[0]=0.02;
			dz[1]=0.05;
			dz[2]=0.11;
			dz[3]=0.23;
			dz[4]=zsnow-dz[0]-dz[2]-dz[1]-dz[3];
		}
		//有雪的时候，每一层雪的质量应该按照新雪和老雪的密度
		for(i = snl; i < 5; i++)
			dz_sum[snl] += dz_sum[i] + dz[i]; 

		if(snow_flag)
		{
			for(i = snl; i<5; i++)
			{
				//雪的厚度小于第一层
				if((snow_depth0 - dz_sum[i] < 0) && (fabs(snow_depth0 - dz_sum[i]) < dz[i]))
				{
					//新雪和老雪混合的层
					//w_i[i] = (dz_sum[i] - snow_depth0) * rou_snow_old + (dz[i]-snow_depth0) * rou_snow_new;
					rou_sno[i] = rou_snow_new * (dz[i]- dz_sum[i] + snow_depth0) / dz_sum[i] + rou_snow_old * (dz_sum[i] - snow_depth0)/dz_sum[i];					
				}
				else if(snow_depth0 - dz_sum[i] < 0 && fabs(snow_depth0 - dz_sum[i]) >= dz[i])
				{
					//都是老雪的层
					//w_i[i] = rou_snow_old * dz[i];
					rou_sno[i] = rou_snow_old;
				}
				else if(snow_depth0 - dz_sum[i] >= 0)
				{
					//都是新雪的层
					//w_i[i] = rou_snow_new * dz[i];
					rou_sno[i] = rou_snow_new;
				}
				//雪的质量
				w_i[i] = rou_sno[i] * dz[i];
				w_l[i] = 0.0;
			}
		}
		else
		{
			for(i = 0;i < 5;i++)
				rou_sno[i] = 0;
		}
		//  value at the beginning of every time step
		if(snow_flag)
		{
			for(i = snl;i < 5;i++)
			{
				w_i1[i] = w_i[i];
				w_l1[i] = w_l[i];
			}	
		}
		for(i = 0;i < N;i++)
		{  
			s1[i] = svar[i];
		}
		sensible_heatflux(&sensiblehf,&latenthf,&shf_dt,&lhf_dt,&ra,snl,metv,pmetvar,vspar);
		sunnetradiation(&netradiation,&Rn_dt,snl,pmetvar);
		groundhf = (netradiation - sensiblehf - latenthf);
		groundhf_dt = Rn_dt - shf_dt - lhf_dt;
		/* calculate valour flux w/m2-->m/s*/ 
		evap = latenthf / (Lvap * rou_w);                 // Lvap is the latent heat of vaporization = 2.501e6 J/kg 
		if(evap<CRIT_PREC)
			evap=0.0;				 
		groundwf = temp_precip + meltwf;//-evap;
		meltwf = 0;
		/******************************Do...While****************************************/
		/* 求解一个时间步长内的温度、水含量do_1217edited  add, subtract and modify at 2014-08-05 */
		do
		{
			flag=0;
			for(i = snl;i < N;i++)
			{  
				s0[i]=svar[i];
				//土壤水势 add
				svar[i].s_h = Waterhead(svar[i].s_w,vspar);
				/**/
				if(svar[i].s_h >= 0)
					svar[i].s_h = -0.0001;
				if(svar[i].s_h< -5000)
					svar[i].s_h = -5000;
				if(i < soiltop)
				{
					lamtaii[i]=lamta_snow(rou_sno[i]);
					C_sno[i] = snow_c(rou_sno[i]);
					Dii[i]= vapour_EDC(svar[i].s_t,metv);
				}
				else
				{
					lamtaii[i] = lamta(svar[i].s_w,svar[i].s_i,vspar);
					Kii[i] = K_water(svar[i].s_w,svar[i].s_i,vspar);
					Dii[i] = Kii[i]*derivationw(svar[i].s_w,svar[i].s_i,vspar);
					derivationh[i] = derivationhead(svar[i].s_h,svar[i].s_i,vspar);
					derivationi[i] = derivationice(svar[i].s_h,svar[i].s_i,vspar);
				}
			}
			for(i = snl+1;i < N-1;i++)
			{
				K_n[i]=(Kii[i]+Kii[i-1])/2;//sqrt(Kii[i]*Kii[i-1]);//2/(1/Kii[i]+1/Kii[i-1]);
				K_s[i]=(Kii[i]+Kii[i+1])/2;//sqrt(Kii[i+1]*Kii[i]);//2/(1/Kii[i]+1/Kii[i+1]);
				D_n[i]=(Dii[i]+Dii[i-1])/2;
				D_s[i]=(Dii[i]+Dii[i+1])/2;
				lamta_n[i]=(lamtaii[i]+lamtaii[i-1])/2;//2/(1/lamtaii[i]+1/lamtaii[i-1]);
				lamta_s[i]=(lamtaii[i]+lamtaii[i+1])/2;//2/(1/lamtaii[i]+1/lamtaii[i+1]);
			}
			lamta_s[snl]=lamta_n[snl+1];
			D_s[snl]=D_n[snl+1];
			K_s[snl]=K_n[snl+1];
			K_n[N-1]=K_s[N-2];//2/(1/Kii[N]+1/Kii[N-1]);K_s[N]=0;
			D_n[N-1]=D_s[N-2];
			lamta_n[N-1]=2/(1/lamtaii[N-1]+1/lamtaii[N-2]);lamta_s[N-1]=lamta_n[N-1];
			//sloved temperature T[i] by equation 
			//coefficient of equation
			if(snow_flag)
			{
				//the top of snow or soil
				ET[snl]=0;
				FT[snl]=lamta_s[snl]/dz[snl]+rou_a*CP/ra;
				GT[snl]=-lamta_s[snl]/dz[snl];
				HT[snl]=netradiation-latenthf+rou_a*CP*pmetvar->ta/ra;
				/* CLM page 87  used at 2012-7-22 */
				/*  ET[snl]=0;
				FT[snl]=1+delta_t / (C_sno[snl]*dz[snl]) * (lamta_s[snl]/(dz[snl]+dz[snl+1])-groundhf_dt);
				GT[snl]=-delta_t/(C_sno[snl]*dz[snl])*(lamta_s[snl]/(dz[snl]+dz[snl+1]));
				HT[snl]=s1[snl].s_t+delta_t/(C_sno[snl]*dz[snl])*(groundhf-groundhf_dt*s1[snl].s_t-lamta_s[snl]*(s1[snl].s_t-s1[snl+1].s_t)/(dz[snl]+dz[snl+1]));
				*/  
				for(i = snl+1;i<soiltop;i++)
				{
					ET[i]=-2*lamta_n[i]/(dz[i]*(dz[i]+dz[i+1]));
					FT[i]=+2*(lamta_n[i]+lamta_s[i])/(dz[i]*(dz[i]+dz[i+1]))+C_sno[i]/delta_t;
					GT[i]=-2*lamta_s[i]/(dz[i]*(dz[i]+dz[i+1]));
					HT[i]=+sink[i]+C_sno[i]*s1[i].s_t/delta_t-rou_w*Lf*(svar[i].s_w-s1[i].s_w)/delta_t;
					/* 
					temp=snow_vdensity_DT(svar[i].s_t);
					ET[i]=2 /(dz[i]*(dz[i]+dz[i+1]))*(lamta_n[i]-Liv*temp*D_n[i]);
					FT[i]=2 * Liv * temp * (D_n[i]+lamta_s[i]) / (dz[i]*(dz[i]+dz[i+1]))-2*(lamta_n[i]+lamta_s[i])/(dz[i]*(dz[i]+dz[i+1]))-(rou_sno[i]*C_sno[i])/delta_t-Liv*temp/delta_t;
					GT[i]=(2 /(dz[i]*(dz[i]+dz[i+1])))*(lamta_s[i]-Liv*temp*D_s[i]);
					HT[i]=-((rou_sno[i]*C_sno[i])/delta_t+Liv*temp/delta_t)*s1[i].s_t;//-sink[i];
					*/ 
				}
				for(i = 5;i < N;i++)
				{
					//solved temperature required varible
					C1[i]=Lf*rou_w*derivationt(svar[i].s_t,svar[i].s_i,vspar);
					Cvv[i]=Cv(svar[i].s_w,svar[i].s_i,vspar);
					Ce[i]=Cvv[i]+C1[i];
					lamtas[i]=lamta_s[i]+D_s[i]*C1[i];
					lamtan[i]=lamta_n[i]+D_n[i]*C1[i];
					Ue[i]=C1[i]*derivationk(svar[i].s_w,svar[i].s_i,vspar);

					R3[i]=delta_t/dz[i];R2[i]=delta_t/(dz[i]+dz[i+1]);
					R1n[i]=2*R3[i]/(dz[i]+dz[i-1]);
					R1s[i]=2*R3[i]/(dz[i]+dz[i+1]);

					ET[i]=-R1n[i]*lamtan[i]-R2[i]*Ue[i];
					GT[i]=-R1s[i]*lamtas[i]+R2[i]*Ue[i];
					FT[i]=Ce[i]-ET[i]-GT[i];
					HT[i]=Ce[i]*s1[i].s_t;
				}
			}
			else
			{
				ET[5]=0;
				FT[5]=Cv(svar[5].s_w,svar[5].s_i,vspar)/delta_t+2*lamta_s[5]/(dz[5]*(dz[5]+dz[6]));
				GT[5]=-2*lamta_s[5]/(dz[5]*(dz[5]+dz[6]));
				HT[5]=Cv(svar[5].s_w,svar[5].s_i,vspar)*s1[5].s_t/delta_t+groundhf/dz[5];//+Lf*rou_i*(svar[0].s_i-s1[0].s_i)/delta_t;
				
				for(i=6;i<N;i++)
				{
					//solved temperature required varible
					C1[i]=Lf*rou_w*derivationt(svar[i].s_t,svar[i].s_i,vspar);
					Cvv[i]=Cv(svar[i].s_w,svar[i].s_i,vspar);
					Ce[i]=Cvv[i]+C1[i];
					lamtas[i]=lamta_s[i]+D_s[i]*C1[i];
					lamtan[i]=lamta_n[i]+D_n[i]*C1[i];
					Ue[i]=C1[i]*derivationk(svar[i].s_w,svar[i].s_i,vspar);

					R3[i]=delta_t/dz[i];R2[i]=delta_t/(dz[i]+dz[i+1]);
					R1n[i]=2*R3[i]/(dz[i]+dz[i-1]);
					R1s[i]=2*R3[i]/(dz[i]+dz[i+1]);

					ET[i]=-R1n[i]*lamtan[i]-R2[i]*Ue[i];
					GT[i]=-R1s[i]*lamtas[i]+R2[i]*Ue[i];
					FT[i]=Ce[i]-ET[i]-GT[i];
					HT[i]=Ce[i]*s1[i].s_t;
				}
			}
			//Qn=0 heat fluxes equal to zero at lower boundary 
			ET[N-1]=-2*lamta_s[N-1]/(dz[N-1]*(dz[N-1]+dz[N-2]));
			FT[N-1]=Cv(svar[N-1].s_w,svar[N-1].s_i,vspar)/delta_t+2*lamta_s[N-1]/(dz[N-1]*(dz[N-1]+dz[N-2]));
			GT[N-1]=0;
			HT[N-1]=Cv(svar[N-1].s_w,svar[N-1].s_i,vspar)*s1[N-1].s_t/delta_t+Lf*rou_i*(svar[N-1].s_i-s1[N-1].s_i);

			//heat;counting upper triangle 
			for(i = snl + 1;i < N;i++)
			{
				FT[i] = FT[i]-ET[i]*GT[i-1]/FT[i-1];
				HT[i] = HT[i]-HT[i-1]*ET[i]/FT[i-1]; 
			}				
			//heat;counting  temprature 
			svar[N-1].s_t= HT[N-1]/FT[N-1];
			for(i = N-2;i>=snl;i--)
			{ 
				svar[i].s_t = (HT[i] - GT[i] * svar[i+1].s_t) / FT[i];
				if(fabs(svar[i].s_t - s0[i].s_t) > 0.1) 
				{
					flag=1;
					if(fabs(svar[i].s_t - s0[i].s_t) > 2)
						svar[i].s_t = 0.5 * (svar[i].s_t - s0[i].s_t);
					}
				}	
			//solved  waterhead by equation (12)  
			//water head  implicit finite difference equation coefficient 
			E[5] = 0;
			F[5] = 1 + 2 * D_s[5] * delta_t / (dz[5] * (dz[5] + dz[6]));
			G[5] = -2*D_s[5]*delta_t/(dz[5]*(dz[5]+dz[6]));
			H[5] = s1[5].s_w+(groundwf-K_s[5])*delta_t/dz[5];//-rou_i/rou_w*(svar[0].s_i-s1[0].s_i);//-derivationh[0]*(WHEAD0[0]-WHEAD1[0])/delta_t;
				
			for(i=6;i<N-1;i++)
			{
				E[i]=-2*D_n[i]*delta_t/(dz[i]*(dz[i]+dz[i-1]));
				G[i]=-2*D_s[i]*delta_t/(dz[i]*(dz[i]+dz[i+1]));
				F[i]=1-E[i]-G[i];
				H[i]=s1[i].s_w-(K_s[i]-K_n[i])*delta_t/(dz[i]+dz[i+1]);//-(rou_i/rou_w)*(svar[i].s_i-s1[i].s_i);//-derivationh[i]*(WHEAD0[i]-WHEAD1[i])/delta_t;
			}
			E[N-1] = -2*D_n[N-1]*delta_t/(dz[N-1]*(dz[N-1]+dz[N-2]));
			F[N-1] = 1-E[N-1]; 
			G[N-1] = 0;
			H[N-1] = s1[N-1].s_w - delta_t * (Kii[N-1] - K_n[N-1]) / dz[N-1];//-(rou_i/rou_w)*(svar[N-1].s_i-s1[N-1].s_i);//-derivationh[N-1]*(WHEAD0[N-1]-WHEAD1[N-1])/delta_t;

			//counting upper triangle FOR Water content
			for(i = 6; i<N; i++)
			{ 
				F[i] = F[i] - E[i] * G[i-1] / F[i-1];
				H[i] = H[i] - H[i-1] * E[i] / F[i-1];
			}
			//counting  water content stu[i] 
			svar[N-1].s_w = H[N-1] / F[N-1];
			for(i = N-2; i >= 5; i--)
			{
				svar[i].s_w = (H[i] - G[i] * svar[i+1].s_w) / F[i];
				if((fabs(svar[i].s_w - s0[i].s_w) / svar[i].s_w) > 0.01) 
				{
					flag = 1;
				}
				if(svar[i].s_w < MinWater) 
					svar[i].s_w = MinWater;
				if(svar[i].s_w > vspar->porosity)
				{
					runoff = runoff + svar[i].s_w - vspar->porosity ;
					svar[i].s_w = vspar->porosity;
				}
			}

			nn++;
			for(i = 5;i < N;i++)
			{				
				if(abs(svar[i].s_t - s1[i].s_t) > 4.0) 
					svar[i].s_t = 0.5 * (svar[i].s_t + s1[i].s_t);//约束时间段间温度差值，防水数据出现波动影响收敛性
			}
		}while(flag && nn<4); //end loop one timestep
		/******************************Do...While****************************************/
		//update ice content,liquid water content
		//Update ice content(THI[i]) by equation (14)//||(svar[i].s_t>=-0.01)&&svar[i].s_i>=0.0001))
		// base on CLM P91 6.42
		if(snow_flag)
		{
			HI[snl] = groundhf+groundhf_dt*(0-s1[snl].s_t)-lamtaii[snl]*(s1[snl].s_t-s1[snl+1].s_t)/(dz[snl]+dz[snl+1])-lamtaii[snl]*(svar[snl].s_t-svar[snl+1].s_t)/(dz[snl]+dz[snl+1])-C_sno[snl]*dz[snl]*(0-s1[snl].s_t)/delta_t;
			for(i = snl+1;i<5;i++)
			{
				HI[i] = -lamtaii[i]*(s1[i].s_t-s1[i+1].s_t)/(dz[i]+dz[i+1])+lamtaii[i-1]*(s1[i-1].s_t-s1[i].s_t)/(dz[i]+dz[i-1])-lamtaii[i]*(svar[i].s_t-svar[i+1].s_t)/(dz[i]+dz[i+1])+lamtaii[i-1]*(svar[i-1].s_t-svar[i].s_t)/(dz[i]+dz[i-1])-C_sno[i]*dz[i]*(0-s1[i].s_t)/delta_t;
			}
			//fprintf(fp2,"\n");
			for(i = snl;i<5;i++)
			{   
				temp_w = HI[i] * delta_t / Lf; 
				if((svar[i].s_t > CRIT_PREC) && (w_i[i]>CRIT_PREC))
				{
					if(temp_w>CRIT_PREC)//melting
					{ 						
						w_i[i] = w_i1[i] - temp_w;
						if(w_i[i]<0) 
							w_i[i]=0;
					}
				}
				if((svar[i].s_t < CRIT_PREC) && (w_l[i]>CRIT_PREC) && (temp_w < 0.0))//freezing
				{
					if(w_l1[i] + temp_w > 0)
						w_i[i] = w_i1[i] - temp_w;
					else
						w_i[i] = w_l1[i] + w_i1[i];
				}
				w_l[i] = w_l1[i] + w_i1[i] - w_i[i];
				
				if(w_l[i]<0) 
					w_l[i]=0.0;
				//fprintf(fp2, "tmp_w=\t%.4f\tw_l[i]=%.2f\n",temp_w,w_l[i]); used to test
				if( svar[i].s_t > CRIT_PREC)
				{ 
					//融化的高度和融化的质量不对应，已修改
					//melthigh = melthigh + dz[i];
					if(w_i1[i] - w_i[i] > 0)
						melthigh = melthigh + (w_i1[i] - w_i[i])/rou_sno[i];
					meltwf +=  w_l[i] / (rou_w * delta_t);
				}
				// w_i[i] is ice mass(kg/m2) in NO.i layer turn to volumetric ice content;w_l[i] is liquid water mass 
				svar[i].s_w=w_l[i]/(dz[i]*rou_w);
				svar[i].s_i=w_i[i]/(dz[i]*rou_i);
				// rou_sno[i]=(w_l[i]+w_i[i])/dz[i];
			}
		}					

		//求冻结土壤的最大含水量及含冰量
		for(i = soiltop;i < N;i++)
		{
			if(abs(svar[i].s_t - s1[i].s_t) > 4.0) 
				svar[i].s_t = 0.5*(svar[i].s_t + s1[i].s_t);//约束时间段间温度差值，防水数据出现波动影响收敛性
			if(svar[i].s_t < CRIT_PREC)
			{
				wm1 = Maxiwater(svar[i].s_t,vspar);
				if(svar[i].s_w > wm1)
				{
					svar[i].s_i = s1[i].s_i+(svar[i].s_w-wm1)*rou_w/rou_i;
					svar[i].s_w = wm1;
				}
				else
				{
					if(svar[i].s_w < MinWater)
						svar[i].s_w = MinWater;
					svar[i].s_i = s1[i].s_i + (s1[i].s_w-svar[i].s_w) * rou_w/rou_i;
				}
				if((svar[i].s_i + svar[i].s_w) > vspar->porosity)
					svar[i].s_i = vspar->porosity - svar[i].s_w - 0.03;
			}
			else
			{
				if(s1[i].s_i>CRIT_PREC)
				{
					svar[i].s_i = s1[i].s_i-C_i*svar[i].s_t/(rou_i*Lf);
					if(svar[i].s_i>0.0)
						svar[i].s_w = svar[i].s_w+C_i*svar[i].s_t/(rou_i*Lf);
					else
					{
						svar[i].s_w = svar[i].s_w + s1[i].s_i * rou_i/rou_w;
						svar[i].s_i = 0;
					}
				}
				else
					svar[i].s_i=0;
			}
			if(svar[i].s_i<CRIT_PREC) 
				svar[i].s_i = 0.0;
		}

		if(!snow_flag)
		{				   
			for(i=0;i<5;i++)
			{   
				svar[i].s_h = 0;
				svar[i].s_w = 0;
				svar[i].s_i = 0;
				svar[i].s_t = 0;
				w_l[i] = 0;
				w_i[i] = 0;
			}
		}
		if(snow_flag)
		{
			for(i = 0;i < snl;i++)
			{  
				svar[i].s_h = 0;
				svar[i].s_w = 0;
				svar[i].s_i = 0;
				svar[i].s_t = 0;
				w_l[i] = 0;
				w_i[i] = 0;		
			}
		}
		//countint waterhead 
		for(i = 0;i < N;i++)
		{
			svar[i].s_h = Waterhead(svar[i].s_w,vspar);//add
			/**/
			if(svar[i].s_h >= 0)
					svar[i].s_h = -0.0001;
			if(svar[i].s_h< -5000)
					svar[i].s_h = -5000;
			//svar[i].s_h = vspar->psi_sat * pow(svar[i].s_w/vspar->porosity,-vspar->CH_b);
			//printf("psi_sat=%lf,water=%lf,porosity=%lf,CH_b=%lf,svar[i].s_h=%lf \n",vspar->psi_sat,svar[i].s_w,vspar->porosity,-vspar->CH_b,svar[i].s_h);
			Tsum[i] +=  svar[i].s_t;

			soil_sum[i].s_t += svar[i].s_t;
		    soil_sum[i].s_tw += svar[i].s_tw;
		    soil_sum[i].s_w += svar[i].s_w;
		    soil_sum[i].s_i += svar[i].s_i;
		    soil_sum[i].s_h += svar[i].s_h;
		}
		//printout the value of every layer at every timestep after of every year
		// if(metday==(simyr+1)*365-1)				
		m = 0.0;
	/*	
		fprintf(fp_temp,"NO.%4d day,%4d \t precip=%.4f\t snl=%d groundwf=%.10f,snow_depth1=%.4f,snow_depth0=%.4f\t,netradiation%.2f\tlatenthf%.2fmelthf%.8f ra=%.2f,snow_flag=%d\n",metday,t_index/3600,metv->prcp,snl,groundwf,snow_depth1,snow_depth0,netradiation,latenthf,melthf, ra,snow_flag);
		for(i = 0;i < N;i++)
		{	
			m=m+dz[i];
			fprintf(fp_temp,"%6.2fm",m); 
			fprintf(fp_temp,"\t%14.4f,%10.4f\t%10.4f\t%8.4f\t%.4f",svar[i].s_h,svar[i].s_w,svar[i].s_i,svar[i].s_t,rou_sno[i]);
			fprintf(fp_temp,"\n");
		}	   
		*/
		/*调整新雪和老雪的厚度*/
		if(newsnow_flag)
		{
			/*如果有新雪存在，则根据融化的雪量，先调整新雪，新雪化完，调整老雪*/
			if(snow_depth0 <= melthigh)
			{					  
				if(snow_depth1 - melthigh + snow_depth0 <= 0)
					snow_depth1 = 0;
				else
					snow_depth1 = snow_depth1 - melthigh + snow_depth0;
				snow_depth0 = 0;
			}
			else
				snow_depth0 = snow_depth0 - melthigh;
		}
		else
		{
			/*否则如果没有新雪，则直接从老雪的量中调整*/
			if(snow_depth1 <= melthigh)
				snow_depth1 = 0;
			else
				snow_depth1 = snow_depth1 - melthigh;
		}
		/* calculate sum of bottom infilation and top evapation */
		sumqn = sumqn + Kii[N-1] * delta_t * svar[N-1].s_w;
		sumevap = sumevap + evap * delta_t;
		sumnr = sumnr + netradiation;
		sumshf = sumshf + sensiblehf;
		sumlhf = sumlhf + latenthf;
		sumghf = sumghf + groundhf;
		//add
		FSM_lrad_down += 0.96 * SBC * pow(pmetvar->ta + 273.15,4);
		FSM_lrad_up += 0.96 * SBC * pow(svar[snl].s_t + 273.15,4);
	}//end loop all timestep(delta_t) within every day
	/******************************For cycle****************************************/
	/* adjust snow layer hight and snow density .Snow density incremental every step time (Anderson,1976)
	考虑雪的压力，融化，重力对雪的密度的影响，对雪的密度进行调整*/
	if(snow_flag)
	{
		snow_depth1 = 0.0;
		rou_snow_old = 0.0; //老雪的平均密度；
		ps = 0.0;
		for(i = snl;i<5;i++)
		{    /* notes: here ps is the weight of snow(unit in cm expressed water equivalent) above the layer of snow */
			ps = ps + 100 *(w_l[i]/rou_w + w_i[i]/rou_i);//the weight of snow(expressed in terms of centimeters of water equivalent) above the layer of snow. w_l[i] and w_i[i] unit is kg/m2
			snow_density_delta1=0.01*ps*pow(FSM_e,(0.08*svar[i].s_t - 21.0*rou_sno[i]/rou_w));  //compaction
			if(rou_sno[i]<150)
				snow_density_delta2=0.01*pow(FSM_e,0.04*svar[i].s_t);   //
			else
			{
				snow_density_delta2=0.01*pow(FSM_e,0.04*svar[i].s_t-46*(rou_sno[i]-150));
				if(svar[i].s_t>0)
					snow_density_delta2=2.0*snow_density_delta2;
			}
			/* the compaction rate due to melting  CLM p107 */
			temp0=w_i1[i]/(w_i1[i]+w_l1[i]);
			temp1=w_i[i]/(w_i[i]+w_l[i]);
			temp=(temp0-temp1)/temp0;
			if(temp>0.0)
				snow_density_delta3=temp;
			else
				snow_density_delta3=0;
			rou_sno[i] = rou_sno[i]*(1+snow_density_delta1+snow_density_delta2+snow_density_delta3);
			if(svar[i].s_t>0)//{ dz[i]=0;w_i[i]=0;w_l[i]=0;svar[i].s_t=0;snl=snl+1;}
				dz[i] = dz[i]*(1 - snow_density_delta1 - snow_density_delta2 - snow_density_delta3);
			snow_depth1 = snow_depth1 + dz[i];
			rou_snow_old += w_i[i] + w_l[i];
		}
		rou_snow_old /= snow_depth1;
	}
	//将半小时的结果求平均到逐日
	//temp = 3600*24/delta_t;
	for(i = 0;i < N;i++)
	{
        soil_sum[i].s_t /= N_HalfHour;
		soil_sum[i].s_tw /= N_HalfHour;
		soil_sum[i].s_w /= N_HalfHour;
		soil_sum[i].s_i /= N_HalfHour;
		soil_sum[i].s_h /= N_HalfHour;
	}
	// output daily average temprature of every layer
	//FSM_psi = (soil_sum[6].s_h + soil_sum[7].s_h + soil_sum[8].s_h + soil_sum[9].s_h + soil_sum[10].s_h + soil_sum[11].s_h)/600.0;
	//FSM_water = (soil_sum[6].s_w + soil_sum[7].s_w + soil_sum[8].s_w + soil_sum[9].s_w + soil_sum[10].s_w + soil_sum[11].s_w)/6.0;
	//FSM_ice = (soil_sum[6].s_i + soil_sum[7].s_i + soil_sum[8].s_i + soil_sum[9].s_i + soil_sum[10].s_i + soil_sum[11].s_i)/6.0;


	for(i = 5;i<= 11;i++)
	{
		//根据实际的土壤深度，如果实际的土壤深度在两层之间，则改层的温度也要加入求平均
		//土壤深度最小2cm，大于第一层。
		/*if(sitec->soil_depth >= dz_sum[i])
		{			
			FSM_soilT += soil_sum[i].s_t  * dz[i];
			FSM_water += soil_sum[i].s_w * dz[i];
			FSM_ice += soil_sum[i].s_i  * dz[i];
			FSM_psi += soil_sum[i].s_h * 0.01 * dz[i];//pa 转为Mpa
			FSM_depth = dz_sum[i];
			if((sitec->soil_depth < dz_sum[i+1]) && (sitec->soil_depth != dz_sum[i]))
			{
				FSM_soilT += soil_sum[i+1].s_t * dz[i+1];
				FSM_water += soil_sum[i+1].s_w * dz[i];
				FSM_ice += soil_sum[i+1].s_i * dz[i];
				FSM_psi += soil_sum[i+1].s_h * 0.01 * dz[i];
				FSM_depth = dz_sum[i+1];
			}
			//printf("i=%d,dz_sum[i]= %lf,FSM_depth = %lf,soil_depth = %lf\n",i,dz_sum[i],FSM_depth,sitec->soil_depth);
		}*/
		FSM_soilT += soil_sum[i].s_t  * dz[i];
		FSM_water += soil_sum[i].s_w * dz[i];
		FSM_ice += soil_sum[i].s_i  * dz[i];
		FSM_psi += soil_sum[i].s_h * 0.01 * dz[i];//pa 转为Mpa
		FSM_depth = 1;
	}
	
	//求全部土壤层的平均值
	FSM_lrad_down /= 48.0;
	FSM_lrad_up /= 48.0;
	FSM_soilT = FSM_soilT / FSM_depth;
	//转换为实际土壤深度的土壤水分重量，kgH20/m2
	FSM_water = FSM_water * sitec->soil_depth * 1000;
	FSM_ice = FSM_ice * sitec->soil_depth * 917;
	FSM_psi = FSM_psi / FSM_depth;
	FSM_snowT = (soil_sum[0].s_t + soil_sum[1].s_t + soil_sum[2].s_t + soil_sum[3].s_t + soil_sum[4].s_t) / 5.0;
	FSM_snowW = (soil_sum[0].s_w * dz[0] + soil_sum[1].s_w * dz[1] + soil_sum[2].s_w * dz[2] + soil_sum[3].s_w * dz[3] + soil_sum[4].s_w * dz[4]) * 1000;

	soil_sum[0].FSM_ice = FSM_ice;
	soil_sum[0].FSM_water = FSM_water;
	soil_sum[0].FSM_soilT = FSM_soilT;
	soil_sum[0].FSM_snowT = FSM_snowT;
	soil_sum[0].FSM_psi = FSM_psi;
	soil_sum[0].FSM_depth = FSM_depth;

	ws->FSM_icew = FSM_ice;
	ws->FSM_soilw = FSM_water;
	ws->FSM_snoww = FSM_snowW;
	//输出：DOY,空气温度,bgc土壤温度,FSM土壤温度,雪层温度,
	//bgc土壤水分,bgc土壤雪,FSM土壤水,FSM土壤冰,FSM冰+水,bgc水势,FSM水势
/*
	fprintf(fp1,"%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,"
		,metday,metv->tavg,metv->tsoil,FSM_soilT
		,FSM_snowT,ws->snoww,FSM_snowW,FSM_water,FSM_ice,(FSM_water + FSM_ice),ws->soilw,epv->psi,FSM_psi,FSM_lrad_up,FSM_lrad_down,bgc_lrad_input);
	for(i = 0;i < N; i++)
		fprintf(fp1,"%lf,",soil_sum[i].s_t);
	for(i = 0;i < N; i++)
		fprintf(fp1,"%lf,",soil_sum[i].s_w);
	fprintf(fp1,"\n");
	*/	
	return(0);
}//end function