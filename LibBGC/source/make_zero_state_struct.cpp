/*2013-08-20 将状态量赋值为0，包括地上，地下*/


/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int make_zero_state_struct(wstate_struct *ws, cstate_struct* cs,
nstate_struct* ns)
{
	int ok=1;
	ws->canopyw = 0.0;
	ws->canopyevap_snk = 0.0;
	ws->trans_snk =0.0;
	ws->outflow_snk =0.0;
	ws->soilevap_snk =0.0;
	ws->snowsubl_snk =0.0;

	cs->cwdc=0.0;		
	cs->leafc=0.0; cs->leafc_storage=0.0; cs->leafc_transfer=0.0;
	cs->livestemc=0.0; cs->livestemc_storage=0.0; cs->livestemc_transfer=0.0;
	cs->deadstemc=0.0; cs->deadstemc_storage=0.0; cs->deadstemc_transfer=0.0;
	cs->livecrootc=0.0; cs->livecrootc_storage=0.0; cs->livecrootc_transfer=0.0;
	cs->deadcrootc=0.0; cs->deadcrootc_storage=0.0; cs->deadcrootc_transfer=0.0;
	cs->frootc=0.0; cs->frootc_storage=0.0; cs->frootc_transfer=0.0;	
	cs->litr1c=0.0; cs->litr2c=0.0; cs->litr3c=0.0; cs->litr4c=0.0;
	cs->soil1c=0.0; cs->soil2c=0.0; cs->soil3c=0.0; cs->soil4c=0.0;
	cs->cpool=0.0; 
	
	ns->cwdn=0.0;
	ns->leafn=0.0; ns->leafn_storage=0.0; ns->leafn_transfer=0.0;
	ns->livestemn=0.0; ns->livestemn_storage=0.0; ns->livestemn_transfer=0.0;
	ns->deadstemn=0.0; ns->deadstemn_storage=0.0; ns->deadstemn_transfer=0.0;
	ns->livecrootn=0.0; ns->livecrootn_storage=0.0; ns->livecrootn_transfer=0.0;
	ns->deadcrootn=0.0; ns->deadcrootn_storage=0.0; ns->deadcrootn_transfer=0.0;
	ns->frootn=0.0; ns->frootn_storage=0.0; ns->frootn_transfer=0.0;
	ns->litr1n=0.0; ns->litr2n=0.0; ns->litr3n=0.0; ns->litr4n=0.0;
	ns->soil1n=0.0; ns->soil2n=0.0; ns->soil3n=0.0; ns->soil4n=0.0;	
	ns->sminn=0.0;  ns->npool=0.0;  ns->retransn=0.0; 


	return(!ok);
}