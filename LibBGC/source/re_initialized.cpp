/* 
atm_pres.c
estimate atmospheric pressure as a function of elevation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
by qiushuai
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int re_init(siteconst_struct* sitec,const epconst_struct* epc, cstate_struct* cs,
cinit_struct* cinit, nstate_struct* ns)
{
	int ok=1;
	double r1[16],r2[16];//r1表示leafc不同类型的leafC2012/leafC1988，r2表示stemC的比率；
	r1[1]=1.370; r1[2]=1.872; r1[3]=2.429; r1[4]=2.429; r1[5]=3.535; r1[6]=1.316;r1[7]=1.774; r1[8]=1.546; 
	r1[9]=1.491; r1[10]=1.777; r1[11]=1.777; r1[12]=1.688; r1[13]=1.777; r1[14]=1.777; r1[15]=1.777;
	r2[1]=1.376; r2[2]=1.941; r2[3]=2.353; r2[4]=2.353; r2[5]=3.409; r2[6]=1.342;r2[7]=1.772; r2[8]=1.537; 
	r2[9]=1.433; r2[10]=1.780; r2[11]=1.780; r2[12]=1.813; r2[13]=1.780; r2[14]=1.780;  r2[15]=1.78;



	//printf("climate_id_site1=%d\n",&sitec->climate_id);
	//1 林地、灌木、竹林、建筑用地转化为草地，4种
	if(sitec->epctype<=17 && sitec->epctype_change==18)
	{
		cs->livestemc=0;
		cs->livestemc_storage=0;
		cs->livecrootc_transfer=0;
		cs->deadstemc=0;
		cs->deadstemc_storage=0;
		cs->deadstemc_transfer=0;
		//粗根转为粗木质部残体。
		cs->cwdc += cs->deadcrootc+cs->deadcrootc_storage+cs->deadcrootc_transfer;
		cs->livecrootc=0;
		cs->livecrootc_storage=0;
		cs->livecrootc_transfer=0;
		cs->deadcrootc=0;
		cs->deadcrootc_storage=0;
		cs->deadcrootc_transfer=0;
		//假如转化之前为建筑用地，则将土壤碳和凋落物初始化。
		if(sitec->epctype==0)
		{
			cs->litr1c=0.1 * 0.069;
			cs->litr2c=0.1 * 0.319;
			cs->litr3c=0.1 * 0.189;
			cs->litr4c=0.1 * 0.423;
			ns->litr2n = cs->litr2c / epc->leaflitr_cn;
			ns->litr3n = cs->litr3c / epc->leaflitr_cn;
			ns->litr4n = cs->litr4c / epc->leaflitr_cn;
			cs->soil1c=4.0 * 0.0037;
			cs->soil2c=4.0 * 0.0133;
			cs->soil3c=4.0 * 0.1334;
			cs->soil4c=4.0 * 0.8496;
			ns->soil1n = cs->soil1c/SOIL1_CN;
			ns->soil2n = cs->soil2c/SOIL2_CN;
			ns->soil3n = cs->soil3c/SOIL3_CN;
			ns->soil4n = cs->soil4c/SOIL4_CN;
		}
	}
	//2 草地转为林地，灌木，竹林2013-07-25，3种
	if(sitec->epctype ==18 && sitec->epctype_change<=17 && sitec->epctype_change>0)
	{
		/*方法一，草地转化为森林的，按照2010年林斑的生物量，反推到发生变化时的生物量。根据431个样地的树木的生长1988-2012
		的比值，按照该比值和发生*/
		if(sitec->epctype_change <= 15)
		{
			cinit->max_leafc = cinit->leafc2010 /r1[(int)sitec->epctype_change] + cinit->leafc2010 *(1- 1/r1[(int)sitec->epctype_change])*(cinit->change_date/8395.0);
			cinit->max_stemc = cinit->stemc2010 /r2[(int)sitec->epctype_change] + cinit->stemc2010 *(1- 1/r2[(int)sitec->epctype_change])*(cinit->change_date/8395.0);
			//printf("草地转林地后stemc=%lf\n",cinit->max_stemc);
		}
		else if(sitec->epctype_change == 16)
		{
			cinit->max_leafc = 0.092;
			cinit->max_stemc = 0.815;
		}
		else if(sitec->epctype_change == 17)
		{
			cinit->max_leafc = 0.036;
			cinit->max_stemc = 0.34;
		}

		/*利用样地中的平均值代替	方法二
		if(sitec->epctype_change == 2)
		{
			cinit->max_leafc = 0.237;
			cinit->max_stemc = 2.7;
		}
		else if(sitec->epctype_change == 3 ||sitec->epctype_change == 4)
		{
			cinit->max_leafc = 0.062;
			cinit->max_stemc = 0.475;
		}
		else if(sitec->epctype_change == 5)
		{
			cinit->max_leafc = 0.083;
			cinit->max_stemc = 1.04;
		}
		else if(sitec->epctype_change == 7)
		{
			cinit->max_leafc = 0.094;
			cinit->max_stemc = 1.19;
		}
		else if(sitec->epctype_change == 9)
		{
			cinit->max_leafc = 0.091;
			cinit->max_stemc = 1.18;
		}
		else if(sitec->epctype_change >=10 && sitec->epctype_change<= 15)
		{
			if(sitec->epctype_change == 12)
			{
				cinit->max_leafc = 0.054;
				cinit->max_stemc = 1.53;
			}
			else
			{
				cinit->max_leafc = 0.026;
				cinit->max_stemc = 0.58;
			}
		}
		
		else if(sitec->epctype_change == 16)
		{
			cinit->max_leafc = 0.092;
			cinit->max_stemc = 0.815;
		}
		else if(sitec->epctype_change == 17)
		{
			cinit->max_leafc = 0.036;
			cinit->max_stemc = 0.34;
		}
		else
		{		
			cinit->max_leafc = 0.12;
			cinit->max_stemc = 1.3;
		}
*/
	}
	//3 有林地，竹林、建筑用地转化为灌木 ，或有林地、建筑用地转化为竹林，5种
	if((sitec->epctype <=15 || sitec->epctype==17) && (sitec->epctype_change==16 ||sitec->epctype_change==17))
	{//前后类型不相同的判断在上个函数中
		if(sitec->epctype_change==16)//转为灌木
		{
			cinit->max_leafc = 0.092;
			cinit->max_stemc = 0.815;
		}
		if(sitec->epctype_change==17 && sitec->epctype <=15)//竹林
		{
			cinit->max_leafc = 0.036;
			cinit->max_stemc = 0.34;
		}
		cs->leafc=0;
		cs->leafc_storage=0;
		cs->leafc_transfer=0;
		cs->livestemc=0;
		cs->livestemc_storage=0;
		cs->livecrootc_transfer=0;
		cs->deadstemc=0;
		cs->deadstemc_storage=0;
		cs->deadstemc_transfer=0;
		//粗根转为粗木质部残体。
		cs->cwdc += cs->deadcrootc+cs->deadcrootc_storage+cs->deadcrootc_transfer;
		cs->livecrootc=0;
		cs->livecrootc_storage=0;
		cs->livecrootc_transfer=0;
		cs->deadcrootc=0;
		cs->deadcrootc_storage=0;
		cs->deadcrootc_transfer=0;
		//如果类型为0，则将土壤碳和凋落物碳初始化
		if(sitec->epctype==0)
		{
			cs->litr1c=0.1 * 0.069;
			cs->litr2c=0.1 * 0.319;
			cs->litr3c=0.1 * 0.189;
			cs->litr4c=0.1 * 0.423;
			ns->litr1n = cs->litr1c / epc->leaflitr_cn;
			ns->litr2n = cs->litr2c / epc->leaflitr_cn;
			ns->litr3n = cs->litr3c / epc->leaflitr_cn;
			ns->litr4n = cs->litr4c / epc->leaflitr_cn;		
			cs->soil1c=4.0 * 0.0037;
			cs->soil2c=4.0 * 0.0133;
			cs->soil3c=4.0 * 0.1334;
			cs->soil4c=4.0 * 0.8496;
			ns->soil1n = cs->soil1c/SOIL1_CN;
			ns->soil2n = cs->soil2c/SOIL2_CN;
			ns->soil3n = cs->soil3c/SOIL3_CN;
			ns->soil4n = cs->soil4c/SOIL4_CN;
		}
		
	}
	//4. 灌木、竹林、建筑用地转化为有林地； 灌木转化为或竹林； 4种
	if((sitec->epctype ==17 || sitec->epctype ==16|| (sitec->epctype ==0 && sitec->epctype_change!=17)) && ((sitec->epctype_change<=15 && sitec->epctype_change>0) || sitec->epctype_change==17))
	{
		if(sitec->epctype != sitec->epctype_change)
		{
			//利用读取的2010年的数据重新初始化。
			if(sitec->epctype_change <= 15)
			{
				cinit->max_leafc = cinit->leafc2010 /r1[(int)sitec->epctype_change] + cinit->leafc2010 *(1- 1/r1[(int)sitec->epctype_change])*(cinit->change_date/8395.0);
				cinit->max_stemc = cinit->stemc2010 /r2[(int)sitec->epctype_change] + cinit->stemc2010 *(1- 1/r2[(int)sitec->epctype_change])*(cinit->change_date/8395.0);
				//printf("灌木、竹林、其它转林地后stemc=%lf\n",cinit->max_stemc);
			}
			else if(sitec->epctype_change == 17)
			{
				cinit->max_leafc = 0.036;
				cinit->max_stemc = 0.34;
			}
			//如果类型为0，则将土壤碳和凋落物碳初始化
			if(sitec->epctype==0)
			{				
				//printf("stemc2010= %lf\n",cinit->stemc2010);
				cinit->max_leafc = cinit->leafc2010 /r1[(int)sitec->epctype_change] + cinit->leafc2010 *(1- 1/r1[(int)sitec->epctype_change])*(cinit->change_date/8395.0);
				cinit->max_stemc = cinit->stemc2010 /r2[(int)sitec->epctype_change] + cinit->stemc2010 *(1- 1/r2[(int)sitec->epctype_change])*(cinit->change_date/8395.0);
				cs->litr1c = cinit->litrc2010 * 0.069;
				cs->litr2c = cinit->litrc2010 * 0.319;
				cs->litr3c = cinit->litrc2010 * 0.189;
				cs->litr4c = cinit->litrc2010 * 0.423;
				ns->litr1n = cs->litr1c / epc->leaflitr_cn;
				ns->litr2n = cs->litr2c / epc->leaflitr_cn;
				ns->litr3n = cs->litr3c / epc->leaflitr_cn;
				ns->litr4n = cs->litr4c / epc->leaflitr_cn;		
				cs->soil1c = cinit->soilc2010 * 0.0037;
				cs->soil2c = cinit->soilc2010 * 0.0133;
				cs->soil3c = cinit->soilc2010 * 0.1334;
				cs->soil4c = cinit->soilc2010 * 0.8496;
				ns->soil1n = cs->soil1c/SOIL1_CN;
				ns->soil2n = cs->soil2c/SOIL2_CN;
				ns->soil3n = cs->soil3c/SOIL3_CN;
				ns->soil4n = cs->soil4c/SOIL4_CN;
			}
/*方法二
			if(sitec->epctype_change == 1 ||sitec->epctype_change == 6)
			{
				cinit->max_leafc = 0.071;
				cinit->max_stemc = 0.871;
			}
			if(sitec->epctype_change == 2)
			{
				cinit->max_leafc = 0.055;
				cinit->max_stemc = 0.933;
			}
			else if(sitec->epctype_change == 3)
			{
				cinit->max_leafc = 0.134;
				cinit->max_stemc = 1.57;
			}
			else if(sitec->epctype_change == 4)
			{
				cinit->max_leafc = 0.059;
				cinit->max_stemc = 0.34;
			}
			else if(sitec->epctype_change == 5)
			{
				cinit->max_leafc = 0.017;
				cinit->max_stemc = 0.2;
			}
			else if(sitec->epctype_change == 7)
			{
				cinit->max_leafc = 0.0499;
				cinit->max_stemc = 5.531;
			}
			else if(sitec->epctype_change == 8)
			{
				cinit->max_leafc = 0.026;
				cinit->max_stemc = 0.514;
			}
			else if(sitec->epctype_change == 9)
			{
				cinit->max_leafc = 0.071;
				cinit->max_stemc = 1.716;
			}
			else if(sitec->epctype_change >=10 && sitec->epctype_change<= 15)
			{
				cinit->max_leafc = 0.036;
				cinit->max_stemc = 0.684;
			}
			else if(sitec->epctype_change == 17)
			{
				cinit->max_leafc = 0.036;
				cinit->max_stemc = 0.34;
			}*/
		}
	}//结束4



	return (!ok);
}