/*
pointbgc.c
front-end to BIOME-BGC for single-point, single-biome simulations
Uses BIOME-BGC function library

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "pointbgc.h"
#include "bgc.h"

#include "time.h"

char *argv_zero = NULL;
signed char cli_mode = MODE_INI;
//soilpar_struct* vspar;				/* add */
soilvar_struct svar[N];         /* add */
//pmet_struct* pmetvar;  
//metvar_struct metv;
//int soilT_error_mark = 0;
int main(int argc, char *argv[])
{
	int npoint = 10000;/*define number of point +1*/
	//static int climate_id;
	/* bgc input and output structures */
	bgcin_struct bgcin;
	bgcout_struct bgcout;

	/* local control information */
	point_struct point;
	restart_ctrl_struct restart;
	climchange_struct scc;
	output_struct output;
	
	/* initialization file */
	file init;
	file ndep_file;

	/*meteorology data initialization parameters. 10,2010*/
	file metfile;
	//char mpathname[100]="..\\..\\input";
	char mpathname[100]="./metdata/";
	//char met[100] = "e:\\model\\@冻土BiomeBGC\\FSM-BGC\\input\\mohesitemetdata1997_2000(snow).txt";
	char met[100];
	/*site condition initialization parameters. 10,2010*/
	file sitefile;
	file epctypefile;
	char spathname[100]=".\\site\\";
	char site[100];
	

	/*epc initialization parameters. 10,2010*/
	file epcfile;
	char epathname[100]=".\\epc\\";
	int epc_check=1;
	char epc[100];
	//int intepc=0;
	int intepc=0,epc_change=0;
	char charepc[100],junk[100],junk1[100];
	double x1,y1;/*used to read the junkhead of the file*/
	double junk_site=0;
	int mark_read = 1,N_begin=1,point_fid,point_id;
	/*output intermediate variables: soiltemprature, longterm average temprature, phenology arrays*/
	/*file lavgtem;
	file intervaule[20000];*/

	/* system time variables */
	struct tm *tm_ptr;
	time_t lt;

	extern signed char summary_sanity;

	int c; /* for getopt cli argument processing */
	extern signed char bgc_verbosity;
	extern int optind, opterr;
	unsigned char bgc_ascii = 1;
	//extern char *optarg;
	extern signed char cli_mode; /* What cli requested mode to run in.*/
	int readndepfile = 0;		/* Flag to tell the program to read an external NDEP file passed using getpopt -n */
	int pointcounter; /*control the point loop*/
	char charcounter[100],met_id1[30],met_id2[30];
	/*输出导度
	FILE *fp_write_gs;
	char gspathname0[100]=".\\daodu\\mco2andvmax_CI",gspathname[100];*/
	/******************************************************************************/
	/*调整平衡态——砍伐*/
	double litter_in=0,croot_in=0;
	//double ratio_cut=0.85,ratio_stem=0.24,ratio_croot=0.21;
	double ratio_cut=1.0,ratio_stem=1.0,ratio_croot=1.0;
	int typec_day,i=0;//tc is type change
	FILE *fp_site;
	fileopen_struct fileopen;
	fileopen.annual_flux_file = 1;
	fileopen.annual_pool_file = 1;
	fileopen.day_pool_file = 1;
	fileopen.month_pool_file = 1;
	fileopen.tsoil_error = 1;
	//FILE *fp_date;
	//动态分配内存
	//vspar = (soilpar_struct*) malloc(sizeof(soilpar_struct)); //add
	//pmetvar = (pmet_struct*) malloc(sizeof(pmet_struct)); //add
	//读取氮沉降数据
	bgcin.ndepctrl.varndep = 0;
	bgcin.cinit.change_date = 0;
	/* Store command name for use by bgc_print_usage() */
	argv_zero = (char *)malloc(strlen(argv[0])+1);
	strncpy(argv_zero, argv[0], strlen(argv[0])+1);

	/* Process command line arguments */
	opterr = 0;
	analysisComm(argc, argv, &bgcin.hModel, &bgcin.laiM, &bgcin.gsiM, &bgcin.pymcM);

	while((c = getopt(argc, argv, (char*)"pVsl:v:ugmn:a")) != -1)
	{
		switch(c)
		{
			case 'V':
				bgc_printf(BV_ERROR, "BiomeBGC version %s (built %s %s by %s on %s)\n", VERS, __DATE__, __TIME__, USER, HOST);
				exit(EXIT_SUCCESS);
				break;
			case 's':
				bgc_verbosity = BV_SILENT;      
				break;
			case 'v':
				bgc_verbosity = bgc_verbosity_decode(optarg);
				break;
			case 'l':
				bgc_logfile_setup(optarg);
				bgc_printf(BV_DIAG, "Using logfile for output.\n");
				break;
			case 'p':
				summary_sanity = SANE;
				break;
			case 'u':
				cli_mode = MODE_SPINUP;
				break;
			case 'm':
				cli_mode = MODE_MODEL;
				strcpy(ndep_file.name,"co2\\ndep.txt");
				//readndepfile = 0;
				//bgcin.ndepctrl.varndep = 0;
				break;
			case 'g':
				cli_mode = MODE_SPINNGO;
				break;
			case 'a':
				bgc_ascii = 1;
				break;
			case 'n':  /* Nitrogen deposition file */
				strcpy(ndep_file.name,optarg);
				bgc_printf(BV_DIAG,"Using annual NDEP file: %s\n",ndep_file.name);
				readndepfile = 1;
				bgcin.ndepctrl.varndep = 1;
				break;
			case '?':
				break;
			default:
				break;
			}
	}

	bgc_printf(BV_DIAG, "Verbosity Level Set To: %d\n", bgc_verbosity);
	
	if (summary_sanity == SANE)
		bgc_printf(BV_WARN, "Summary outputs will be calculated more sanely. See USAGE.TXT for details\n");

	if (cli_mode != MODE_INI)
	{
		//bgc_printf(BV_WARN, "Overridding ini mode. ");
		if (cli_mode == MODE_SPINUP)
			bgc_printf(BV_WARN, "Running in Spinup Mode.\n");
		if (cli_mode == MODE_MODEL)
			//bgc_printf(BV_WARN, "Running in Model mode.\n");
		if (cli_mode == MODE_SPINNGO)
			bgc_printf(BV_WARN, "Running in Spin-and-Go mode.\nThe spinup and model will both be run.\n");
	}
		
	bgc_printf(BV_DIAG, "Done processing CLI arguments.\n");

	/* get the system time at start of simulation */
	lt = time(NULL);
	tm_ptr = localtime(&lt);
	strcpy(point.systime,asctime(tm_ptr));
	/* Andrew tried this, you shouldn't. localtime returns a global extern. */
	/* free(tm_ptr); */
	output.anncodes = NULL;
	output.daycodes = NULL;
	output.bgc_ascii = bgc_ascii;
	
	/* initialize the bgcin state variable structures before filling with
	values from ini file */
	if (presim_state_init(&bgcin.ws, &bgcin.cs, &bgcin.ns, &bgcin.cinit))
	{
		bgc_printf(BV_ERROR, "Error in call to presim_state_init() from pointbgc()\n");
		exit(EXIT_FAILURE);
	}

	srand((unsigned int)time(NULL));
	/******************************
	**                           **
	**  BEGIN READING INIT FILE  **
	**                           **
	******************************/

	/* read the name of the main init file from the command line
	and store as init.name */

	/*打开记录转换日期的文件*/
	//fp_date=fopen(".\\地类转换日期1.txt","w");

	/*begin the point loop*/
	fp_site = fopen(".\\site\\site.csv", "r");
	while((fscanf(fp_site,"%[^,],%*[^\n]",junk1))!=EOF)
		i++;
	npoint = i-2;
	//printf("npoint = %d\n",npoint);
	fclose(fp_site);
	N_begin = 1;

for(pointcounter = N_begin; pointcounter<= npoint; pointcounter++)
{
	bgcin.cs.soilT_error_mark = 0;
	memset(charcounter,0,100);
	_itoa(pointcounter, charcounter, 10);
	if((pointcounter-1) % 100 == 0)
		//printf("第%d个点\n",pointcounter);

	if(pointcounter % 100 == 0)
	{	
		printf("*");		
	}
	for(i=0;i<N;i++)
	{
		svar[i].s_h = 0;
		svar[i].s_i = 0;
		svar[i].s_t = 0;
		svar[i].s_tw = 0;
		svar[i].s_w = 0;
		svar[i].FSM_depth = 0;
		svar[i].FSM_ice = 0;
		svar[i].FSM_psi = 0;
		svar[i].FSM_snowT = 0;
		svar[i].FSM_soilT = 0;
		svar[i].FSM_water =0;
		svar[i].FSM_water_total =0;
	}
	/*new parts to read meteorology file, site condition file, and epc file*/
	bgcin.sitec.pointcounter = pointcounter;
	//printf("转化的日期%d\n",typec_day);
	/*empty arrays*/
	//memset(met,0,100);
	memset(site,0,100);
	memset(epc,0,100);
	memset(charepc,0,100);
	memset(restart.out_restart.name,0,128);

	/*打开导度文件，开始写入
	strcpy(gspathname,gspathname0);
	strcat(gspathname,charcounter);
	strcat(gspathname,".txt");
	fp_write_gs=fopen(gspathname,"w");*/

	/*read init file*/
	if (optind >= argc)
	{
		bgc_print_usage();
		exit(EXIT_FAILURE);
	}
	strcpy(init.name, argv[optind]);
	/* open the main init file for ascii read and check for errors */
	if (file_open(&init,'i'))
	{
		bgc_printf(BV_ERROR, "Error opening init file, pointbgc.c\n");
		exit(EXIT_FAILURE);
	}

	/* read the header string from the init file */
	if (fgets(point.header, 100, init.ptr)==NULL)
	{
		bgc_printf(BV_ERROR, "Error reading header string: pointbgc.c\n");
		exit(EXIT_FAILURE);
	}
	//printf("读取init文件头\n");
	/* open met file, discard header lines *//*remove 10.3,2010*/
	/*if (met_init(init, &point))
	{
		bgc_printf(BV_ERROR, "Error in call to met_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}*/

	/* read restart control parameters */
	if (restart_init(init, &restart,charcounter))
	{
		bgc_printf(BV_ERROR, "Error in call to restart_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}

	/* read simulation timing control parameters */
	if (time_init(init, &(bgcin.ctrl)))
	{
		bgc_printf(BV_ERROR, "Error in call to epclist_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}

	/* read scalar climate change parameters */
	if (scc_init(init, &scc))
	{
		bgc_printf(BV_ERROR, "Error in call to scc_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}

	/* read CO2 control parameters */
	if (co2_init(init, &(bgcin.co2), bgcin.ctrl.simyears))
	{
		bgc_printf(BV_ERROR, "Error in call to co2_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}

	if(readndepfile)
	{
		if (ndep_init(ndep_file, &(bgcin.ndepctrl)))
		{
			bgc_printf(BV_ERROR, "Error in call to ndep_init() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
	}

	/* read site constants */
	if (sitec_init(init, &bgcin.sitec))
	{
		bgc_printf(BV_ERROR, "Error in call to sitec_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}

	/*get the name and open site condition file, only open once when i=0*/
	if(pointcounter == 1)
	{
		if(mark_read)
		{
			strcpy(site,spathname);
			strcat(site,"site.csv");
			//strcat(site,init.name);
			if ((sitefile.ptr = fopen(site,"r")) == NULL)
			{
				bgc_printf(BV_ERROR, "\n site文件打开错误 ... Exiting\n");
			}
			fopen_s(&(sitefile.ptr),site,"r");		
			fscanf(sitefile.ptr,"%[^,],%*[^\n]",junk1);//read the header of the file
			//printf("site读取成功\n");
			strcpy(site,spathname);
			strcat(site,"epc.txt");
			//fopen_s(&(epctypefile.ptr),site,"r");
			if((epctypefile.ptr = fopen(site,"r")) == NULL)
			{
				bgc_printf(BV_ERROR, "\n epctype文件打开错误 ... Exiting\n");
			}
			fscanf(epctypefile.ptr,"%s %*[^\n]",junk);//read the header of the file 2013-3-30
		}
		mark_read = 0;
	}
	else
	{
		if(mark_read)
		{
			strcpy(site,spathname);
			strcat(site,"site.csv");
			//strcat(site,init.name);
			if ((sitefile.ptr = fopen(site,"r")) == NULL)
			{
				bgc_printf(BV_ERROR, "\n site文件读取错误 ... Exiting\n");
			}
			fopen_s(&(sitefile.ptr),site,"r");		
			fscanf(sitefile.ptr,"%[^,],%*[^\n]",junk1);//read the header of the file
			//printf("site读取成功\n");
			strcpy(site,spathname);
			strcat(site,"epc.txt");
			//fopen_s(&(epctypefile.ptr),site,"r");
			if((epctypefile.ptr = fopen(site,"r")) == NULL)
			{
				bgc_printf(BV_ERROR, "\n epctype文件打开错误 ... Exiting\n");
			}
			fscanf(epctypefile.ptr,"%s %*[^\n]",junk);//read the header of the file 2013-3-30

			for(pointcounter = 1; pointcounter < N_begin; pointcounter++)
			{
				fscanf(sitefile.ptr,"%lf,%*[^\n]",&junk_site);//read the header of the file
				fscanf(epctypefile.ptr,"%lf%*[^\n]",&junk_site);//read the header of the file 2013-3-30
			}
			printf("跳过已算的记录\n");
		}
		mark_read = 0;
	}
	if(fscanf(epctypefile.ptr,"%d%d",&intepc,&epc_change)==NULL)
		intepc = -1;
	//fscanf(epctypefile.ptr,"%d%d",&intepc,&epc_change);
	//printf("epc=%d\n",intepc);
	bgcin.sitec.epctype = (double)intepc;
	bgcin.sitec.epctype_change = (double)epc_change;
	bgcin.epc.epctype = (double)intepc;
	bgcin.epc.epctype_change = (double)epc_change;
	//如果都是森林类型1-15，则按照当前的类型计算，而且不考虑类型的变化。2013-08-20.
	if(bgcin.sitec.epctype>0 && bgcin.sitec.epctype<16 && bgcin.sitec.epctype_change>0&& bgcin.sitec.epctype_change<16)
		bgcin.sitec.epctype = bgcin.sitec.epctype_change;
	if(bgcin.sitec.epctype == 19)
	{
		//如果转化之后的类型是19，则将其赋值为变化之前的数值，表示不变；
		bgcin.sitec.epctype = bgcin.sitec.epctype_change;
		intepc=(int)bgcin.sitec.epctype_change;
		//printf("类型为19,替换后为：%lf\n",bgcin.sitec.epctype);
	}
	bgcin.cinit.epctype = bgcin.sitec.epctype;
	bgcin.cinit.epctype_change = bgcin.sitec.epctype_change;

	//printf("epctype=%d\n",bgcin.sitec.epctype);
	/* read ramped nitrogen deposition block */
	if (ramp_ndep_init(init, &bgcin.ramp_ndep))
	{
		bgc_printf(BV_ERROR, "Error in call to ramp_ndep_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}
	//printf("读取soil n fixation\n");
	/* read ecophysiological constants */
	/*if (epc_init(init, &bgcin.epc))
	{
		bgc_printf(BV_ERROR, "Error in call to epc_init() from pointbgc.c... Exiting\n");
		exit(EXIT_FAILURE);
	}*/


	/*read epclist and epc parameter*/
	//2013-08-19修改，如果变化前后的类型都为0，则跳过该点，否则，仍然模拟。
	if(bgcin.sitec.epctype==0 && bgcin.sitec.epctype_change==0)
	{
		//intepc=epc_check;
		printf("第%d个点不读\n",pointcounter);
		fscanf(sitefile.ptr,"%[^,],%*[^\n]",junk);//跳过该行
		//都为0，类型不变化，则日期为-1
		//fprintf(fp_date,"%lf\n",bgcin.cinit.change_date);
	//	if(intepc==-1)//到达最后一个文件
	//		fclose(fp_date);
		if (bgcin.co2.varco2)
		{
			free(bgcin.co2.co2ppm_array);
			free(bgcin.co2.co2year_array);
		}
		if(restart.read_restart) fclose(restart.in_restart.ptr);
		if(restart.write_restart) fclose(restart.out_restart.ptr);
		fclose(init.ptr);

	}
	else
	{
		if(bgcin.sitec.epctype != bgcin.sitec.epctype_change)
		{
			/*产生随机数，发生地类变化的时间0~8395. 2013-08-20*/
			typec_day = (int)(rand()%8395);
			if(typec_day<2)
				typec_day=2;
			bgcin.cinit.change_date = (double)typec_day;
			//printf("地类变化的日期%d\n",typec_day);
		//	fprintf(fp_date,"%lf\n",bgcin.cinit.change_date);
		//	if(intepc==-1)//到达最后一个文件
		//		fclose(fp_date);
		}
		else
		{
		//	fprintf(fp_date,"%lf\n",bgcin.cinit.change_date);
		//	if(intepc==-1)//到达最后一个文件
		//		fclose(fp_date);
		}
		_itoa(intepc, charepc, 10);
		strcpy(epc,epathname);		
		strcpy(bgcin.sitec.epcfilepath,epathname);
		strcat(epc,charepc);
		strcat(epc,".epc");
		//printf("epc文件地址%s\n",epc);
		if((epcfile.ptr = fopen(epc,"r")) == NULL)
        {
			bgc_printf(BV_ERROR, "\n epc文件读取错误 ... Exiting\n");
			printf("point=%d\n",pointcounter);
			//printf("epc文件地址%s\n",epc);
		}
		//printf("epc文件地址%s\n",epc);
		if (epc_init(epcfile, &bgcin.epc, &bgcin.sitec))
		{
			bgc_printf(BV_ERROR, "Error in call to epc_init() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		fclose(epcfile.ptr);
	
		/*new site constants initialization function is added. 10,2010
		add initial pools,such as leafc,stemc,litrc,soilc,litrn,simn,2013.02.26*/	
	//	for(i=0;i<bgcin.co2.co2vals;i++)	
	//		printf("co2year = %d,co2ppm= %lf\n",bgcin.co2.co2year_array[i],bgcin.co2.co2ppm_array[i]);

		if (sitec_init2(sitefile, &bgcin.sitec, &bgcin.epc, &bgcin.cs, &bgcin.cinit,&bgcin.ns,&bgcin.co2))
		{
			bgc_printf(BV_ERROR, "Error in call to sitec_init2() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		bgcin.vspar.psandy = bgcin.sitec.sand;
		bgcin.vspar.pclay = bgcin.sitec.clay;
		bgcin.vspar.psilt = bgcin.sitec.silt;
		//printf("main_climate_i=%d\n",bgcin.sitec.climate_id);
	//	for(i=0;i<bgcin.co2.co2vals;i++)	
	//		printf("co2year = %d,co2ppm= %lf\n",bgcin.co2.co2year_array[i],bgcin.co2.co2ppm_array[i]);

		if(pointcounter == npoint)
		{
			fclose(sitefile.ptr);
			fclose(epctypefile.ptr);
		}

		/* initialize water state structure */
		if (wstate_init(init, &bgcin.sitec, &bgcin.ws))
		{
			bgc_printf(BV_ERROR, "Error in call to wstate_init() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		//printf("读取water structure\n");
		/* initialize carbon and nitrogen state structures 
		modified 2013-02-26
		if (cnstate_init(init, &bgcin.epc, &bgcin.cs, &bgcin.cinit,
			&bgcin.ns))
		{
			bgc_printf(BV_ERROR, "Error in call to cstate_init() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}*/
		//printf("读取cn structure\n");
		/* read the output control information */
		if (output_ctrl(init, &output,charcounter))
		{
			bgc_printf(BV_ERROR, "Error in call to output_ctrl() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		
		/* initialize output files. Does nothing in spinup mode*/
		/***************************************************************************/
		//如果变化前后
		//if(bgcin.sitec.epctype || bgcin.sitec.epctype_change)
		//{
		if (output_init(&output,&fileopen))
		{
			bgc_printf(BV_ERROR, "Error in call to output_init() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		//}
		/* read final line out of init file to test for proper file structure */
		if (end_init(init))
		{
			bgc_printf(BV_ERROR, "Error in call to end_init() from pointbgc.c... exiting\n");
			exit(EXIT_FAILURE);
		}
		fclose(init.ptr);

		/* get the name and open met file */
		/*"G:\\metdata\\input"*/
		/*point_fid = (int)floor((float)(bgcin.sitec.climate_id - 1)/30080.0) + 1;
		point_id = (bgcin.sitec.climate_id - 1)%30080 + 1;*/
		point_fid = 1;
		point_id = bgcin.sitec.climate_id;
		_itoa(point_fid,met_id1,10); //根据ID号查找气候数据的目录  input file 上级目录
		_itoa(point_id,met_id2,10);	//根据ID号查找气候数据文件
	    strcpy(met,mpathname);
		//strcat(met,met_id1);
		//strcat(met,"\\input");
		//添加
		strcat(met,met_id2);
		strcat(met,".txt");
		//tempt
		/*strcpy(met,".\\input1.txt");*/
		if ((metfile.ptr = fopen(met,"r")) == NULL)        
		{
			bgc_printf(BV_ERROR, "\n气候数据读取错误 %s ... Exiting\n",met);
			printf("soilc_海拔=%lf\n",bgcin.sitec.elev);
		}
		//printf("打开met文件\n");
		//printf(met);
		/* read meteorology file, build metarr arrays, compute running avgs */
		if (metarr_init(metfile, &bgcin.metarr, &scc, bgcin.ctrl.metyears))
		{
			bgc_printf(BV_ERROR, "Error in call to metarr_init() from pointbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		//printf("metarr_init complete \n");
		fclose(metfile.ptr);
		//printf("完成metarr_init\n");
		//}
		/* copy some of the info from input structure to bgc simulation control
		structure */
		bgcin.ctrl.onscreen = output.onscreen;
		bgcin.ctrl.dodaily = output.dodaily;
		bgcin.ctrl.domonavg = output.domonavg;
		bgcin.ctrl.doannavg = output.doannavg;
		bgcin.ctrl.doannual = output.doannual;
		bgcin.ctrl.ndayout = output.ndayout;
		bgcin.ctrl.nannout = output.nannout;
		bgcin.ctrl.daycodes = output.daycodes;
		bgcin.ctrl.anncodes = output.anncodes;
		bgcin.ctrl.read_restart = restart.read_restart;
		bgcin.ctrl.write_restart = restart.write_restart;
		bgcin.ctrl.keep_metyr = restart.keep_metyr;
		
		/* copy the output file structures into bgcout */
		if (output.dodaily) bgcout.dayout = output.dayout;
		if (output.domonavg) bgcout.monavgout = output.monavgout;
		if (output.doannavg) bgcout.annavgout = output.annavgout;
		if (output.doannual) bgcout.annout = output.annout;
		if (output.bgc_ascii && output.dodaily) bgcout.dayoutascii = output.dayoutascii;
		if (output.bgc_ascii && output.dodaily) bgcout.dayout_soil = output.dayout_soil;  //add
		if (output.bgc_ascii && output.domonavg) bgcout.monoutascii = output.monoutascii;
		if (output.bgc_ascii && output.doannual) bgcout.annoutascii = output.annoutascii;
		bgcout.anntext = output.anntext;
		bgcout.bgc_ascii = bgc_ascii;
		bgcout.fp_tsoil_error = output.fp_tsoil_error;
		/* if using ramped Ndep, copy preindustrial Ndep into ramp_ndep struct */
		if (bgcin.ramp_ndep.doramp)
		{
			bgcin.ramp_ndep.preind_ndep = bgcin.sitec.ndep;
		}
		//printf("开始读取restart文件\n");
		/* if using an input restart file, read a record */
		if (restart.read_restart)
		{
			/* 02/06/04
			 * The if statement gaurds against core dump on bad restart file.
			 * If spinup exits with error then the norm trys to use the restart,
			 * that has nothing in it, a seg fault occurs. Amac */
			if( fread(&(bgcin.restart_input),sizeof(restart_data_struct),1,restart.in_restart.ptr) == 0)
			{
				bgc_printf(BV_ERROR, "Error reading restart file! 0 bytes read. Aborting..\n");
				exit(EXIT_FAILURE);
			}
			/*此处修改文件的初始化，利用新的文件对restart中的变量进行赋值*/
		}
		//printf("读完restart文件\n");

	//	for(i=0;i<bgcin.co2.co2vals;i++)	
		//	printf("co2year = %d,co2ppm= %lf\n",bgcin.co2.co2year_array[i],bgcin.co2.co2ppm_array[i]);

		/*********************
		**                  **
		**  CALL BIOME-BGC  **
		**                  **
		*********************/

		/* all initialization complete, call model */
		/* either call the spinup code or the normal simulation code */

		if (bgcin.ctrl.spinup)
		{
			if (bgc(&bgcin, &bgcout, MODE_SPINUP))
			{
				bgc_printf(BV_ERROR, "Error in call to bgc()\n");
				exit(EXIT_FAILURE);
			}
			bgc_printf(BV_PROGRESS, "SPINUP: residual trend  = %.6lf\n",bgcout.spinup_resid_trend);
			bgc_printf(BV_PROGRESS, "SPINUP: number of years = %d\n",bgcout.spinup_years);
		}
		else
		{
			/**************************************************************************************************/
			//printf("CO2 value read: %i %lf\n",bgcin.co2.co2year_array[10],bgcin.co2.co2ppm_array[10]); 
			if (bgc(&bgcin, &bgcout, MODE_MODEL))
			{
				bgc_printf(BV_ERROR, "Error in call to bgc()\n");
				//break;
				exit(EXIT_FAILURE);
			}
			//printf("CO2 value read: %i %lf\n",bgcin.co2.co2year_array[10],bgcin.co2.co2ppm_array[10]); 
			//printf("bgc()主程序完成\n");
		}

		/* if using an output restart file, write a record */
		if (restart.write_restart)
		{
			/**************************************************************************************/
			/*添加砍伐，按照60%比例采伐*/
			/*调整平衡态的碳库，假设采伐之后立即达到平衡*/
			/*将计算细根和10%的叶子
			litter_in=0.25*ratio_cut*(bgcout.restart_output.frootc+bgcout.restart_output.frootc_storage+bgcout.restart_output.frootc_transfer+
					  0.1*(bgcout.restart_output.leafc+bgcout.restart_output.leafc_storage+bgcout.restart_output.leafc_transfer));*/
			/*将粗根移动到粗木质残体
			croot_in=ratio_cut*(bgcout.restart_output.livecrootc+bgcout.restart_output.livecrootc_storage+bgcout.restart_output.livecrootc_transfer+
					 bgcout.restart_output.deadcrootc+bgcout.restart_output.deadcrootc_storage+bgcout.restart_output.deadcrootc_transfer);
			bgcout.restart_output.cwdc+=croot_in;*/
			/*将细根和10%的叶子平均分配到掉落物中
			bgcout.restart_output.litr1c+=litter_in;
			bgcout.restart_output.litr2c+=litter_in;
			bgcout.restart_output.litr3c+=litter_in;
			bgcout.restart_output.litr4c+=litter_in;*/
			/*修改被砍伐的部分，叶子和细根以及粗根，茎=0*/
			if(bgcin.epc.woody)
			{
				bgcout.restart_output.leafc*=ratio_cut;
				bgcout.restart_output.leafc_storage*=ratio_cut;
				bgcout.restart_output.leafc_transfer*=ratio_cut;
				bgcout.restart_output.frootc*=ratio_cut;
				bgcout.restart_output.frootc_storage*=ratio_cut;
				bgcout.restart_output.frootc_transfer*=ratio_cut;
				bgcout.restart_output.livecrootc*=ratio_croot;
				bgcout.restart_output.livecrootc_storage*=ratio_croot;
				bgcout.restart_output.livecrootc_transfer*=ratio_croot;
				bgcout.restart_output.deadcrootc*=ratio_croot;
				bgcout.restart_output.deadcrootc_storage*=ratio_croot;
				bgcout.restart_output.deadcrootc_transfer*=ratio_croot;
				bgcout.restart_output.livestemc*=ratio_stem;
				bgcout.restart_output.livestemc_storage*=ratio_stem;
				bgcout.restart_output.livestemc_transfer*=ratio_stem;
				bgcout.restart_output.deadstemc*=ratio_stem;
				bgcout.restart_output.deadstemc_storage*=ratio_stem;
				bgcout.restart_output.deadstemc_transfer*=ratio_stem;
			}
			fwrite(&(bgcout.restart_output),sizeof(restart_data_struct),1,
				restart.out_restart.ptr);
		}
		/* Now do the Model part of Spin & Go. */
		if (cli_mode == MODE_SPINNGO)
		{
			bgc_printf(BV_PROGRESS, "Finished Spinup for Spin 'n Go. Now starting Model run ('Go' part of Spin'n Go)\n");
				
			bgc_printf(BV_PROGRESS, "Assigned bgcout struct to bgcin for spinngo model run\n");
			
			bgcin.ctrl.spinup = 0;
			output.doannavg = 1;
			output.doannual = 1;
			output.dodaily = 1;
			output.domonavg = 1;
			
			if (output_init(&output,&fileopen))
			{
				bgc_printf(BV_ERROR, "Error in call to output_init() from pointbgc.c... Exiting\n");
				exit(EXIT_FAILURE);
			}
			
			/* copy some of the info from input structure to bgc simulation control structure */
			bgcin.ctrl.dodaily = output.dodaily;
			bgcin.ctrl.domonavg = output.domonavg;
			bgcin.ctrl.doannavg = output.doannavg;
			bgcin.ctrl.doannual = output.doannual;

			/* copy the output file structures into bgcout */
			if (output.dodaily) bgcout.dayout = output.dayout;
			if (output.domonavg) bgcout.monavgout = output.monavgout;
			if (output.doannavg) bgcout.annavgout = output.annavgout;
			if (output.doannual) bgcout.annout = output.annout;
			if (output.bgc_ascii && output.dodaily) bgcout.dayoutascii = output.dayoutascii;
			if (output.bgc_ascii && output.domonavg) bgcout.monoutascii = output.monoutascii;
			if (output.bgc_ascii && output.doannual) bgcout.annoutascii = output.annoutascii;
			if (output.bgc_ascii && output.doannual) bgcout.anntext = output.anntext;
			
			/* initialize output files. Does nothing in spinup mode*/	
			bgcin.ctrl.read_restart = 1;
			bgcin.restart_input = bgcout.restart_output;
			
			if (bgc(&bgcin, &bgcout, MODE_MODEL))
			{
				bgc_printf(BV_ERROR, "Error in call to bgc()\n");
				exit(EXIT_FAILURE);
				//break;
			}
			restart.read_restart = 0;
			bgcin.ctrl.read_restart = 0;

			bgc_printf(BV_WARN, "Finished the bgc() Model call in spinngo\n");
			
		}

		/* post-processing output handling, if any, goes here */
		/* free memory *//*应该是每个点计算完了释放，还是所有点都计算完了再释放*/
		free(bgcin.metarr.tmax);
		free(bgcin.metarr.tmin);
		free(bgcin.metarr.prcp);
		free(bgcin.metarr.vpd);
		free(bgcin.metarr.tavg);
		free(bgcin.metarr.tavg_ra);
		free(bgcin.metarr.swavgfd);
		free(bgcin.metarr.par);
		free(bgcin.metarr.dayl);
		free(bgcin.metarr.wspeed);
		free(bgcin.metarr.lrad);
		free(bgcin.metarr.rhumid);
		//free(bgcin.metarr.tday);
		
		//free(vspar);    //add
		//free(pmetvar); //add
		if (bgcin.co2.varco2) free(bgcin.co2.co2ppm_array);
		if (bgcin.ndepctrl.varndep) free(bgcin.ndepctrl.ndepyear_array);
		if (bgcin.ndepctrl.varndep) free(bgcin.ndepctrl.ndep_array);
		if (output.anncodes != NULL) free(output.anncodes);
		if (output.daycodes != NULL) free(output.daycodes);
		/* close files */
		if (restart.read_restart) fclose(restart.in_restart.ptr);
		if (restart.write_restart) 
		{
			if (fclose(restart.out_restart.ptr) != 0)
			{
				bgc_printf(BV_WARN, "Warning, error closing restart file after write: %s\n", strerror(errno));
			}
		}

		//if (output.bgc_ascii && output.dodaily) fclose(output.dayout_soil.ptr);//add
		
		//if (output.dodaily) fclose(output.dayout.ptr);
		//if (output.domonavg) fclose(output.monavgout.ptr);
		//if (output.doannavg) fclose(output.annavgout.ptr);
		//if (output.doannual) fclose(output.annout.ptr);
		/* Close the ASCII output files */
		/*———————————————如果类型为0的话，不关闭文件—————————————————*/
		
		if(bgcin.sitec.epctype)
		{
		//	if (output.bgc_ascii && output.dodaily) fclose(output.dayoutascii.ptr);//日ascii文件
		//	if (output.bgc_ascii && output.domonavg) fclose(output.monoutascii.ptr);//月ascii文件
		}
			/**/
		/*write the soilT error file*/
		fprintf(output.fp_tsoil_error.ptr,"%10d%6d\n",bgcin.sitec.plotID,bgcin.cs.soilT_error_mark);
		//printf("soilT_error= %d\n",bgcin.cs.soilT_error_mark);
		if (output.bgc_ascii && output.doannual && (pointcounter == npoint))
		{
			fclose(output.annoutascii.ptr);//年ascii文件
			//printf("输出文件annual_pool关闭\n");
		}
		if (output.bgc_ascii && output.doannual && (pointcounter == npoint) && (fclose(output.anntext.ptr) != 0))
		{
			bgc_printf(BV_WARN, "Warning, error closing ascii annual output file: %s\n", strerror(errno));
			//printf("输出文件annual_flux关闭\n");
		}
		if (output.bgc_ascii && output.doannual && (pointcounter == npoint))
		{
			fclose(output.fp_tsoil_error.ptr);//年ascii文件
			//printf("输出文件annual_pool关闭\n");
		}
		if (output.bgc_ascii && output.dodaily && (pointcounter == npoint))
		{
			fclose(output.dayoutascii.ptr);//年ascii文件
			//printf("输出文件daily_pool关闭\n");
		}
		if (output.bgc_ascii && output.domonavg && (pointcounter == npoint))
		{
			fclose(output.monoutascii.ptr);//年ascii文件
			//printf("输出文件month_pool关闭\n");
		}
	//	}


	}/* 结束if(bgcin.sitec.epctype)，用于跳过epctype为0的点*/

	//printf("第%d个点完成\n",pointcounter);

}          /*end of the point loop,2010.10.03*/

	bgc_logfile_finish();
	free(argv_zero);
	return EXIT_SUCCESS;
} /* end of main */
	
