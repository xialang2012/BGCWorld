#ifndef SOIL_FUNC_H
#define SOIL_FUNC_H
/*
bgc_func.h
header file for function prototypes

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/
#include "bgc.h"

#ifdef __cplusplus
extern "C"
{
#endif

//int snowmelt(const metvar_struct* metv, wflux_struct* wf, double snoww);
/*double Cv(double z_thw,double z_thi);
double Kersten(double thw,double thi);
double lamta_sat(double z_thw,double z_thi);
double lamta(double z_thw,double z_thi);
double K_water(double z_thw,double z_thi);
double Watercontent1(double waterhead,double ice);
double Watercontent(double waterhead);
double Waterhead(double thw);
double Waterhead1(double thw,double thi);
double Waterhead2(double thw,double thi);
double derivationhead(double waterhead,double icecontent);
double derivationice(double waterhead,double icecontent);
double derivationw(double thw,double icecontent);
double derivationt(double tem,double icecontent);
double derivationk(double thw,double thi);
double Maxiwater(double T);

double icecont(double tem,double thw);
double vaporpressure_sat(double temp);
//void sensible_heatflux(double* vaporheatflux,double* latenthf);
int heatflux_suf(soilvar_struct ssvar[],double airvpd,double tair,double* vaporheatflux);
int sunnetradiation(soilvar_struct ssvar[],double tair,double* Rn,double rad);

*/
double Maxiwater1(double T,double icecontent);
/*add 2014-07-30*/

double Cv(double z_thw,double z_thi,soilpar_struct* vspar);//the volumetric specific heat capacity of soil
double snow_c(double rou_snow);//the volumetric specific heat capacity of snow
double Kersten(double thw,double thi,soilpar_struct* vspar);

double lamta_sat(double z_thw,double z_thi,soilpar_struct* vspar);
double lamta(double z_thw,double z_thi,soilpar_struct* vspar);
double K_water(double z_thw,double z_thi,soilpar_struct* vspar);
double Maxiwater(double T,soilpar_struct* vspar);
double Waterhead(double thw,soilpar_struct* vspar);
double Waterhead1(double thw,double thi,soilpar_struct* vspar);
double Waterhead2(double thw,double thi,soilpar_struct* vspar);
double derivationhead(double waterhead,double icecontent,soilpar_struct* vspar);
double derivationice(double waterhead,double icecontent,soilpar_struct* vspar);
double derivationw(double thw,double icecontent,soilpar_struct* vspar);
double derivationt(double tem,double icecontent,soilpar_struct* vspar);
double derivationk(double thw,double thi,soilpar_struct* vspar);

double vaporpressure_sat(double temp);
double specifichumiditysat_dt(double tem,metvar_struct* metv);
double lamta_snow(double rou_snow);
double snow_vdensity_DT(double tem);
double vapour_EDC(double tem,metvar_struct* metv);

void sensible_heatflux(double* vaporheatflux,double* latenthf,double* vhf_dt,double* lhf_dt,double* rb ,int snl,metvar_struct* metv,pmet_struct *pmetvar,soilpar_struct* vspar);
void svar_init(metarr_struct* metarr);//TODO初始化函数，加入BGC模型中，需要修改。利用BGC模型的结果来初始化，不从文件中读取
void sunnetradiation(double* Rn,double* Rn_dt,int snl,pmet_struct *pmetvar);//TODO加入模型后需要修改，利用BGC计算的净辐射





#ifdef __cplusplus
}
#endif

#endif
