/* 
output_init.c
Reads output control information from initialization file

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions from version 4.1.2:
Seperated init scanning and outputfile opening.
This function now opens output files. Output controls
are now read from output_ctrl.c

Revisions from version 4.1.1:
Fixed error in ascii output file that incorrectly gave the 
units for annual precipitation as cm/year - the real units are mm/yr.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "pointbgc.h"

int output_init(output_struct* output,fileopen_struct *fileopen)
{
  int ok = 1;
	
	/* open outfiles if specified */
	if (ok && output->dodaily)
	{
		strcpy(output->dayout.name,output->outprefix);
		strcat(output->dayout.name,".dayout");
		/*if (file_open(&(output->dayout),'w'))
		{
			bgc_printf(BV_ERROR, "Error opening daily outfile (%s) in output_ctrl()\n",output->dayout.name);
			ok=0;
		}
		else
		{
			bgc_printf(BV_WARN, "Opened binary daily output file in write mode\n");
		}
		*/
	}
	
	if (ok && output->domonavg)
	{
		strcpy(output->monavgout.name,output->outprefix);
		strcat(output->monavgout.name,".monavgout");
		/*if (file_open(&(output->monavgout),'w'))
		{
			bgc_printf(BV_ERROR, "Error opening monthly average outfile (%s) in output_ctrl()\n",output->monavgout.name);
			ok=0;
		}*/
	}
	if (ok && output->doannavg)
	{
		strcpy(output->annavgout.name,output->outprefix);
		strcat(output->annavgout.name,".annavgout");
		/*if (file_open(&(output->annavgout),'w'))
		{
			bgc_printf(BV_ERROR, "Error opening annual average outfile (%s) in output_ctrl()\n",output->annavgout.name);
			ok=0;
		}*/
	}
	if (ok && output->doannual)
	{
		strcpy(output->annout.name,output->outprefix);
		strcat(output->annout.name,".annout");
		/*if (file_open(&(output->annout),'w'))
		{
			bgc_printf(BV_ERROR, "Error opening annual outfile (%s) in output_ctrl()\n",output->annout.name);
			ok=0;
		}*/
	}
	/****************************************/
	/*					*/
	/* 		ASCII Outputs		*/
	/*					*/
	/****************************************/
		
	/* open daily ascii output files if specified */
	if (ok && output->bgc_ascii && output->dodaily)
	{
		strcpy(output->dayoutascii.name,output->outprefix);
		strcat(output->dayoutascii.name,".dayout.ascii");
		
		if(fileopen->day_pool_file)
		{
			/**/
			if (file_open(&(output->dayoutascii),'o'))
			{
				bgc_printf(BV_ERROR, "Error opening daily ascii outfile (%s) in output_ctrl()\n",output->dayoutascii.name);
				ok=0;
			}
			//fprintf(output->dayoutascii.ptr,"%s %s %s %s %s\n","NPP","NEE","HR","GPP","MR");
			fileopen->day_pool_file = 0;
		}
		
	}
	
	if (ok && output->bgc_ascii && output->domonavg)
	{
		strcpy(output->monoutascii.name,output->outprefix);
		strcat(output->monoutascii.name,".monavgout.ascii");
		/**/
		if(fileopen->month_pool_file)
		{
			if (file_open(&(output->monoutascii),'o'))
			{
				bgc_printf(BV_ERROR, "Error opening monthly ascii outfile (%s) in output_ctrl()\n",output->monoutascii.name);
				ok=0;
			}
			fileopen->month_pool_file = 0;
			//printf("monoutasciifile open\n");
		}
	}
	if (ok && output->bgc_ascii && output->doannual)
	{
		strcpy(output->annoutascii.name,output->outprefix);
		strcat(output->annoutascii.name,".annout.ascii");
		if(fileopen->annual_pool_file)
		{
			if (file_open(&(output->annoutascii),'o'))
			{
				bgc_printf(BV_ERROR, "Error opening annual ascii outfile (%s) in output_ctrl()\n",output->annoutascii.name);
				ok=0;
			}
			fprintf(output->annoutascii.ptr,"%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s%13s\n","plotID","soilWater","soilSnow","LAI","vegC",
					"litrC","soilC","totalC","leafC_sum","leafC_max","leafC","stemC","crootC","cwdc","FSM_soilw","FSM_icew","FSM_snoww");
			fileopen->annual_pool_file = 0;
			//printf("annoutasciifile open\n");
		}
	}
	
	/****************************************/
	/*					*/
	/* 	End ASCII Outputs		*/
	/*					*/
	/****************************************/
	
	/* Yeah, we need make this not happen on spinup */
	/* Spinngo is going to make this look like doodoo */
	if (ok && output->bgc_ascii && output->doannual)
	{
		/* simple text output */
		strcpy(output->anntext.name,output->outprefix);
		strcat(output->anntext.name,"_ann.txt");
		if(fileopen->annual_flux_file)
		{
			if (file_open(&(output->anntext),'o'))
			{
				bgc_printf(BV_ERROR, "Error opening annual text file (%s) in output_ctrl()\n",output->anntext.name);
				ok=0;
			}
			fprintf(output->anntext.ptr,"%10s%6s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%12s%12s%15s%15s\n",
				"plotID","year","PRCP","Tavg","max_LAI","ET","Trans","Evapor","CanopyW","OF","NPP","NBP","NEP","NEE","GPP","MR","GR","HR","root_MR","leaf_MR"
				,"bgc土壤水","bgc雪");
			/*fprintf(output->anntext.ptr,"%10s%6s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%12s%12s%15s%15s%15s%15s%15s%10s\n",
				"plotID","year","PRCP","Tavg","max_LAI","ET","Trans","Evapor","CanopyW","OF","NPP","NBP","NEP","NEE","GPP","MR","GR","HR","root_MR","leaf_MR"
				,"bgc土壤水","bgc雪","FSM土壤水","FSM冰","FSM雪","c价值");*/
			fileopen->annual_flux_file = 0;
			//printf("annoutasciifile open\n");
		}

		strcpy(output->fp_tsoil_error.name,output->outprefix);
		strcat(output->fp_tsoil_error.name,"_tsoil_Error.txt");
		if(fileopen->tsoil_error)
		{
			if (file_open(&(output->fp_tsoil_error),'o'))
			{
				bgc_printf(BV_ERROR, "Error opening annual text file (%s) in output_ctrl()\n",output->anntext.name);
				ok=0;
			}
			fprintf(output->fp_tsoil_error.ptr,"%10s%10s\n","plotID","TsoilErr");
			fileopen->tsoil_error = 0;
			//printf("soil T error recorded file open\n");
		}
		/* write the header info for simple text file 
		fprintf(output->anntext.ptr,"Annual summary output from Biome-BGC version %s\n",VERS);
		fprintf(output->anntext.ptr,"ann PRCP = annual total precipitation (mm/yr)\n");
		fprintf(output->anntext.ptr,"ann Tavg = annual average air temperature (deg C)\n");
		fprintf(output->anntext.ptr,"max LAI = annual maximum value of projected leaf area index (m2/m2)\n");
		fprintf(output->anntext.ptr,"ann ET = annual total evapotranspiration (mm/yr)\n");
		fprintf(output->anntext.ptr,"ann OF = annual total outflow (mm/yr)\n");
		fprintf(output->anntext.ptr,"ann NPP = annual total net primary production (gC/m2/yr)\n");
		fprintf(output->anntext.ptr,"ann NBP = annual total net biome production (gC/m2/yr)\n");
		fprintf(output->anntext.ptr,"ann NEP = annual total net ecosystem production (gC/m2/yr)\n");
		fprintf(output->anntext.ptr,"ann NEE = annual total net ecosystem  exchange (gC/m2/yr)\n");
		fprintf(output->anntext.ptr,"ann GPP = annual total gross primary production (gC/m2/yr)\n");
		fprintf(output->anntext.ptr,"ann MR = annual total maintenance respiration (gC/m2/yr)\n");
		fprintf(output->anntext.ptr,"ann GR = annual total growth respiration (gC/m2/yr)\n\n");
		fprintf(output->anntext.ptr,"ann HR = annual total heterotrophic respiration (gC/m2/yr)\n\n");*/
	}
	return (!ok);
}

