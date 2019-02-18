#ifndef POINTBGC_FUNC_H
#define POINTBGC_FUNC_H
/*
pointbgc_func.h
function prototypes for pointbgc
for use with pointbgc front-end to BIOME-BGC library

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#ifdef __cplusplus
extern "C"
{
#endif

int met_init(file init, point_struct* point);
int restart_init(file init, restart_ctrl_struct* restart, char* pointid);
int time_init(file init, control_struct *ctrl);
int scc_init(file init, climchange_struct* scc);
int co2_init(file init, co2control_struct* co2, int simyears);
int ndep_init(file ndepfile, ndepcontrol_struct* ndepctrl);
int sitec_init(file init, siteconst_struct* sitec);
//int sitec_init2(file sitefile,siteconst_struct* sitec);//new sitec initialization function. 10,2010
int sitec_init2(file sitefile, siteconst_struct* sitec,const epconst_struct* epc, 
				cstate_struct* cs,cinit_struct* cinit, nstate_struct* ns,co2control_struct* co2);//,int* climate_id);
//2013-03-26

int ramp_ndep_init(file init, ramp_ndep_struct* ramp_ndep);
//int epc_init(file init, epconst_struct* epc);
int epc_init(file epcfile, epconst_struct* epc, siteconst_struct* sitec);//new epc initialization function. 10,2010
int wstate_init(file init, const siteconst_struct* sitec, wstate_struct* ws);
int cnstate_init(file init, const epconst_struct* epc, cstate_struct* cs,
cinit_struct* cinit, nstate_struct* ns);
//int output_ctrl(file init, output_struct* output);
int output_ctrl(file init, output_struct* output, char* charcounter);//new output control function. 10,2010
int output_init(output_struct* output,fileopen_struct *fileopen);
int end_init(file init);
int metarr_init(file metf, metarr_struct* metarr, const climchange_struct* scc,
int nyears);
int presim_state_init(wstate_struct* ws, cstate_struct* cs, nstate_struct* ns,
cinit_struct* cinit);
int prephenology(const control_struct* ctrl, const epconst_struct* epc, 
const siteconst_struct* sitec, const metarr_struct* metarr,
phenarray_struct* phen);//add prephenology function 

#ifdef __cplusplus
}
#endif

#endif
