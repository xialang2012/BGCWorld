#ifndef BGC_STRUCT_H
#define BGC_STRUCT_H

/*
bgc_struct.h 
header file for structure definitions

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include <iostream>
#include <fstream>

#ifdef __cplusplus
extern "C"
{
#endif

/* simulation control variables */
typedef struct
{
	int metyears;          /* # years of met data */
    int simyears;          /* # years of simulation */
    int simstartyear;      /* first year of simulation */
	int spinup;            /* (flag) 1=spinup run, 0=normal run */
	int maxspinyears;      /* maximum number of years for spinup run */
	int dodaily;           /* flag for daily output */
	int domonavg;          /* flag for monthly average of daily outputs */
	int doannavg;          /* flag for annual average of daily outputs */
	int doannual;          /* flag for annual output */
	int ndayout;           /* number of daily outputs */
	int nannout;           /* number of annual outputs */
	int* daycodes;         /* array of indices for daily outputs */
	int* anncodes;         /* array of indices for annual outputs */
	int read_restart;      /* flag to read restart file */
	int write_restart;     /* flag to write restart file */
	int keep_metyr;        /* (flag) 1=retain restart metyr, 0=reset metyr */
	int onscreen;          /* (flag) 1=show progress on-screen 0=don't */
} control_struct;

/* a structure to hold information about ramped N-deposition scenario */
typedef struct
{
	int doramp;         /* (flag) 1=ramped Ndep, 0=constant Ndep */
	int ind_year;          /* (int)  reference year for indust. Ndep */
	double preind_ndep;    /* (double)  (kgN/m2/yr) preindustrial Ndep (at first metyear) */
	double ind_ndep;       /* (double)  (kgN/m2/yr) industrial Ndep at ref yr */
} ramp_ndep_struct;

/* a structure to hold information on the annual co2 concentration */
typedef struct
{
	int varco2;             /* (flag) 0=const 1=use file 2=const,file for Ndep */
    double co2ppm;          /* (ppm)  constant CO2 concentration */
	double* co2ppm_array;   /* (ppm)  annual CO2 concentration array */
	int* co2year_array;		/* (year) year corresponding to the concentration value in co2ppm_arry */
	int co2vals;			/* (num)  The number of CO2 concentration values in the co2ppm_array */
} co2control_struct;	

/* a structure to hold annual nitrogen deposition data */
typedef struct
{
	int varndep;            /* (flag) 0=const 1=use file  */
    double ndep;			/* (kgN m-2 yr-1)  constant ndep value */
	double* ndep_array;		/* (kgN m-2 yr-1)  annual ndep array */
	int* ndepyear_array;	/* (year) year corresponding to the ndep value in ndep_array */
	int ndepvals;			/* (num)  The number of ndep values in the ndep_array */
} ndepcontrol_struct;	

/* meteorological variable arrays */
/* inputs from mtclim, except for tavg and tavg_ra
which are used for an 11-day running average of daily average air T,
computed for the whole length of the met array prior to the 
daily model loop */

typedef struct
{
    double* tmax;          /* (0.1deg C) daily maximum air temperature */
    double* tmin;          /* (0.1deg C) daily minimum air temperature */
	double* vpd;           /* (Pa)    vapor pressure deficit */
	double* swavgfd;       /* (W/m2)  daylight avg shortwave flux density */
	double* prcp;          /* (0.1mm)    precipitation */
	double* par;           /* (W/m2)  photosynthetically active radiation */
	double* dayl;          /* (s)     daylength */
	double* tavg;          /* (0.1deg C) daily average temperature */
    double* tavg_ra;       /* (deg C) 11-day running avg of daily avg temp */
    //added
	//double* snd;           /* (cm)    snowdepth  */
	//double* snp;            /* 0.1g/cm2 snow press*/
	//double* den;           /*(kg/m3)  snow density */    
	//double* atp;           /* (0.1hPa)    atmospheric press      */
	//double* vap;           /* (0.1Pa)    vapour press           */
	double* rhumid;        /* (%)     relative humidity      */    
	double* wspeed;        /* (0.1m/s-1)  wind speed */    
	double* lrad;
	//double* tday;          /* (0.1deg C) daily average air temperature */
} metarr_struct;
/* daily values that are passed to daily model subroutines */
typedef struct
{
	double prcp;           /* (kg/m2) precipitation */
	double tmax;           /* (deg C) daily maximum air temperature */
	double tmin;           /* (deg C) daily minimum air temperature */
	double tavg;           /* (deg C) daily average air temperature */
	double tday;           /* (deg C) daylight average air temperature */
	double tnight;         /* (deg C) nightime average air temperature */
	double tsoil;          /* (deg C) daily soil temperature, avg, top 10 cm */ 
	double vpd;            /* (Pa)    vapor pressure deficit */
	double swavgfd;        /* (W/m2)  daylight average shortwave flux */
	double swabs;          /* (W/m2)  canopy absorbed shortwave flux */ 
	double swtrans;        /* (W/m2)  transmitted shortwave flux */
	double swabs_per_plaisun;   /* (W/m2) swabs per unit sunlit proj LAI */
	double swabs_per_plaishade; /* (W/m2) swabs per unit shaded proj LAI */
	double ppfd_per_plaisun;    /* (umol/m2/s) ppfd per unit sunlit proj LAI */
	double ppfd_per_plaishade;  /* (umol/m2/s) ppfd per unit shaded proj LAI */
	double par;            /* (W/m2)  photosynthetically active radiation */
	double parabs;         /* (W/m2)  PAR absorbed by canopy */
	double pa;             /* (Pa)    atmospheric pressure */
	double co2;            /* (ppm)   atmospheric concentration of CO2 */
	double dayl;           /* (s)     daylength */
	double wspeed;			/* (0.1m/s-1)  wind speed */ 
	double rhumid;			/*相对湿度*/
	double lrad;			/*长波辐射*/
} metvar_struct;
//extern metvar_struct metv;
/* water state variables (including sums for sources and sinks) */
typedef struct
{
    double soilw;          /* (kgH2O/m2) water stored in soil */
    double snoww;          /* (kgH2O/m2) water stored in snowpack */
    double canopyw;        /* (kgH2O/m2) water stored on canopy */
    double prcp_src;       /* (kgH2O/m2) SUM of precipitation */
    double outflow_snk;    /* (kgH2O/m2) SUM of outflow */
    double soilevap_snk;   /* (kgH2O/m2) SUM of soil water evaporation */
    double snowsubl_snk;   /* (kgH2O/m2) SUM of snow water sublimation */
    double canopyevap_snk; /* (kgH2O/m2) SUM of canopy water evaporation */
    double trans_snk;      /* (kgH2O/m2) SUM of transpiration */
	double FSM_soilw;      /* 冻土模块输出的土壤水分, KgH20/m2*/ 
	double FSM_icew;	   /* 冻土模块输出的土壤冰, KgH20/m2*/ 
	double FSM_snoww;	   /* 冻土模块输出的土壤雪, KgH20/m2*/ 
	double soilw_add;
} wstate_struct;
/*new added struct 2014-07-30*/
/* soiltype parameters */
typedef struct 
{
	int    typeNo;          /* soil texture types code numer: 1 sand,2 loamysand,3 sandyloam,4 loam,5 siltyloam,6 sandyclayloam,7 clayloam,8 siltyclayloam,9 sandyclay,10 siltyclay,11 clay,12 peat,13 moss,14 lichen */
	double psandy;          /* (%)percentage of sand */
	double psilt;           /* (%)percentage of slit */
	double pclay;           /* (%)percentage of clay */
	double k_solids;        /* (Wm-1K-1)the thermal conductivity of the soil solids */
	double c_solids;        /* (Jm-3K-1)the heat capacity of the soil solids */
	double k_sat;           /* (ms-1)the saturated hydraulic conductivity */
	double psi_sat;         /* (m)the saturated matric potential */
	double porosity;        /* the saturated volumetric water content(porosity) */
	double CH_b;            /* the Clap and Hornberger constant */
}soilpar_struct;
//extern soilpar_struct* vspar;
/*new added struct 2014-07-30*/
/* soil every layer variables */
typedef struct
{
    double s_t;           /* (*C)soil temperature */
	double s_tw;          /* soil total water conten */
	double s_w;           /* soil liquid water content */
	double s_i;           /* soil ice content */
	double s_h;           /* (m)soil water head  */
//	double s_tt;
//	double s_w_max;       /* the maxmum of the soil liquid water content */
//	double s_kt;          /* (Wm-3K-1)the thermal conductivity of the soil */
//	double s_ct;          /* (Jm-3K-1)the capacity of the soil */
//	double s_kw;          /* (ms-1)the hydralic conductivity of the soil */ 
	double FSM_psi;		  /* FSM模型计算的日平均水势，按照实际的土壤深度 */
	double FSM_water;	  /* FSM模型计算的日平均水分，按照实际的土壤深度 */
	double FSM_ice;
	double FSM_water_total;/* water + ice */
	double FSM_soilT;
	double FSM_snowT;
	double FSM_depth;
}soilvar_struct;
/*new added struct used for mark the status of a file which is open or not*/
typedef struct
{
	int annual_pool_file;
	int annual_flux_file;
	int month_pool_file;
	int month_flux_file;
	int day_pool_file;
	int day_flux_file;
	int tsoil_error;
}fileopen_struct;


typedef struct
{
	double ta;              /* (degC) air temperature */
	double pa;              /* (Pa)   air pressure */
	double vpd;             /* (Pa)   vapor pressure deficit */
	double irad;            /* (W/m2) incident shortwave flux density */
	double rv;              /* (s/m)  resistance to water vapor flux */
	double rh;              /* (s/m)  resistance to sensible heat flux */
	double windspeed;          /* (m/s)  wind speed */
	double lrad;
} pmet_struct;
//extern pmet_struct* pmetvar;

/* water flux variables */
typedef struct
{
    double prcp_to_canopyw;   /* (kgH2O/m2/d) interception on canopy */
    double prcp_to_soilw;     /* (kgH2O/m2/d) precip entering soilwater pool */
    double prcp_to_snoww;     /* (kgH2O/m2/d) snowpack accumulation */
    double canopyw_evap;      /* (kgH2O/m2/d) evaporation from canopy */
    double canopyw_to_soilw;  /* (kgH2O/m2/d) canopy drip and stemflow */
    double snoww_subl;        /* (kgH2O/m2/d) sublimation from snowpack */
    double snoww_to_soilw;    /* (kgH2O/m2/d) melt from snowpack */
    double soilw_evap;        /* (kgH2O/m2/d) evaporation from soil */
    double soilw_trans;       /* (kgH2O/m2/d) transpiration */     
    double soilw_outflow;     /* (kgH2O/m2/d) outflow */
	double et;                /* (kgH20/m2/d) evapotranspiration */
} wflux_struct;

/* carbon state initialization structure */
typedef struct
{
	double max_leafc;      /* (kgC/m2) first-year displayed + stored leafc */
	double max_stemc;      /* (kgC/m2) first-year total stem carbon */
	double litrc1988;	   /* (kgC/m2) litter carbon of 1988*/
	double litrc2010;	   /* (kgC/m2) litter carbon of 2010*/
	double soilc1988;	   /* (kgC/m2) soil carbon of 2010*/
	double soilc2010;	   /* (kgC/m2) soil carbon of 2010*/
	double leafc2010;      /* (kgC/m2) max leaf carbon of 2010*/
	double stemc2010;      /* (kgC/m2) max stem carbon of 2010*/
	double croot1988;
	double croot2010;
	double epctype;
	double epctype_change;
	double change_date;	   /* 发生类型变化的日期*/
	int frost;
} cinit_struct;

/* carbon state variables (including sums for sources and sinks) */
typedef struct 
{
    double leafc;            /* (kgC/m2) leaf C */
    double leafc_storage;    /* (kgC/m2) leaf C storage */
    double leafc_transfer;   /* (kgC/m2) leaf C transfer */
    double frootc;           /* (kgC/m2) fine root C */
    double frootc_storage;   /* (kgC/m2) fine root C storage */
    double frootc_transfer;  /* (kgC/m2) fine root C transfer */
    double livestemc;        /* (kgC/m2) live stem C */
    double livestemc_storage;  /* (kgC/m2) live stem C storage */
    double livestemc_transfer; /* (kgC/m2) live stem C transfer */
    double deadstemc;          /* (kgC/m2) dead stem C */
    double deadstemc_storage;  /* (kgC/m2) dead stem C storage */
    double deadstemc_transfer; /* (kgC/m2) dead stem C transfer */
    double livecrootc;         /* (kgC/m2) live coarse root C */
    double livecrootc_storage; /* (kgC/m2) live coarse root C storage */
    double livecrootc_transfer;/* (kgC/m2) live coarse root C transfer */
    double deadcrootc;         /* (kgC/m2) dead coarse root C */
    double deadcrootc_storage; /* (kgC/m2) dead coarse root C storage */
    double deadcrootc_transfer;/* (kgC/m2) dead coarse root C transfer */
	double gresp_storage;    /* (kgC/m2) growth respiration storage */
	double gresp_transfer;   /* (kgC/m2) growth respiration transfer */
    double cwdc;             /* (kgC/m2) coarse woody debris C */
    double litr1c;           /* (kgC/m2) litter labile C */
    double litr2c;           /* (kgC/m2) litter unshielded cellulose C */
    double litr3c;           /* (kgC/m2) litter shielded cellulose C */
    double litr4c;           /* (kgC/m2) litter lignin C */
    double soil1c;           /* (kgC/m2) microbial recycling pool C (fast) */
    double soil2c;           /* (kgC/m2) microbial recycling pool C (medium) */
    double soil3c;           /* (kgC/m2) microbial recycling pool C (slow) */
    double soil4c;           /* (kgC/m2) recalcitrant SOM C (humus, slowest) */
	double cpool;            /* (kgC/m2) temporary photosynthate C pool */
    double psnsun_src;       /* (kgC/m2) SUM of gross PSN from sulit canopy */
    double psnshade_src;     /* (kgC/m2) SUM of gross PSN from shaded canopy */
    double leaf_mr_snk;      /* (kgC/m2) SUM of leaf maint resp */
    double leaf_gr_snk;      /* (kgC/m2) SUM of leaf growth resp */
    double froot_mr_snk;     /* (kgC/m2) SUM of fine root maint resp */
    double froot_gr_snk;     /* (kgC/m2) SUM of fine root growth resp */
    double livestem_mr_snk;  /* (kgC/m2) SUM of live stem maint resp */
    double livestem_gr_snk;  /* (kgC/m2) SUM of live stem growth resp */
    double deadstem_gr_snk;  /* (kgC/m2) SUM of dead stem growth resp */
    double livecroot_mr_snk; /* (kgC/m2) SUM of live coarse root maint resp */
    double livecroot_gr_snk; /* (kgC/m2) SUM of live coarse root growth resp */
    double deadcroot_gr_snk; /* (kgC/m2) SUM of dead coarse root growth resp */
    double litr1_hr_snk;     /* (kgC/m2) SUM of labile litr microbial resp */
    double litr2_hr_snk;     /* (kgC/m2) SUM of cellulose litr microbial resp */
    double litr4_hr_snk;     /* (kgC/m2) SUM of lignin litr microbial resp */
    double soil1_hr_snk;     /* (kgC/m2) SUM of fast microbial respiration */
    double soil2_hr_snk;     /* (kgC/m2) SUM of medium microbial respiration */
    double soil3_hr_snk;     /* (kgC/m2) SUM of slow microbial respiration */
    double soil4_hr_snk;     /* (kgC/m2) SUM of recalcitrant SOM respiration */
	double fire_snk;         /* (kgC/m2) SUM of fire losses */
	int soilT_error_mark; /*标记土壤温度是否出错*/
} cstate_struct;

/* daily carbon flux variables */
typedef struct
{
	/* mortality fluxes */ 
	double m_leafc_to_litr1c;              /* (kgC/m2/d) */            
	double m_leafc_to_litr2c;              /* (kgC/m2/d) */
	double m_leafc_to_litr3c;              /* (kgC/m2/d) */
	double m_leafc_to_litr4c;              /* (kgC/m2/d) */
	double m_frootc_to_litr1c;             /* (kgC/m2/d) */
	double m_frootc_to_litr2c;             /* (kgC/m2/d) */
	double m_frootc_to_litr3c;             /* (kgC/m2/d) */
	double m_frootc_to_litr4c;             /* (kgC/m2/d) */
	double m_leafc_storage_to_litr1c;      /* (kgC/m2/d) */
	double m_frootc_storage_to_litr1c;     /* (kgC/m2/d) */
	double m_livestemc_storage_to_litr1c;  /* (kgC/m2/d) */
	double m_deadstemc_storage_to_litr1c;  /* (kgC/m2/d) */
	double m_livecrootc_storage_to_litr1c; /* (kgC/m2/d) */
	double m_deadcrootc_storage_to_litr1c; /* (kgC/m2/d) */
	double m_leafc_transfer_to_litr1c;     /* (kgC/m2/d) */
	double m_frootc_transfer_to_litr1c;    /* (kgC/m2/d) */
	double m_livestemc_transfer_to_litr1c; /* (kgC/m2/d) */
	double m_deadstemc_transfer_to_litr1c; /* (kgC/m2/d) */
	double m_livecrootc_transfer_to_litr1c;/* (kgC/m2/d) */
	double m_deadcrootc_transfer_to_litr1c;/* (kgC/m2/d) */
	double m_livestemc_to_cwdc;            /* (kgC/m2/d) */
	double m_deadstemc_to_cwdc;            /* (kgC/m2/d) */
	double m_livecrootc_to_cwdc;           /* (kgC/m2/d) */
	double m_deadcrootc_to_cwdc;           /* (kgC/m2/d) */
	double m_gresp_storage_to_litr1c;      /* (kgC/m2/d) */
	double m_gresp_transfer_to_litr1c;     /* (kgC/m2/d) */
	/* fire fluxes */
	double m_leafc_to_fire;                /* (kgC/m2/d) */
	double m_frootc_to_fire;               /* (kgC/m2/d) */
	double m_leafc_storage_to_fire;        /* (kgC/m2/d) */
	double m_frootc_storage_to_fire;       /* (kgC/m2/d) */
	double m_livestemc_storage_to_fire;    /* (kgC/m2/d) */
	double m_deadstemc_storage_to_fire;    /* (kgC/m2/d) */
	double m_livecrootc_storage_to_fire;   /* (kgC/m2/d) */
	double m_deadcrootc_storage_to_fire;   /* (kgC/m2/d) */
	double m_leafc_transfer_to_fire;       /* (kgC/m2/d) */
	double m_frootc_transfer_to_fire;      /* (kgC/m2/d) */
	double m_livestemc_transfer_to_fire;   /* (kgC/m2/d) */
	double m_deadstemc_transfer_to_fire;   /* (kgC/m2/d) */
	double m_livecrootc_transfer_to_fire;  /* (kgC/m2/d) */
	double m_deadcrootc_transfer_to_fire;  /* (kgC/m2/d) */
	double m_livestemc_to_fire;            /* (kgC/m2/d) */
	double m_deadstemc_to_fire;            /* (kgC/m2/d) */
	double m_livecrootc_to_fire;           /* (kgC/m2/d) */
	double m_deadcrootc_to_fire;           /* (kgC/m2/d) */
	double m_gresp_storage_to_fire;        /* (kgC/m2/d) */
	double m_gresp_transfer_to_fire;       /* (kgC/m2/d) */
	double m_litr1c_to_fire;               /* (kgC/m2/d) */
	double m_litr2c_to_fire;               /* (kgC/m2/d) */
	double m_litr3c_to_fire;               /* (kgC/m2/d) */
	double m_litr4c_to_fire;               /* (kgC/m2/d) */
	double m_cwdc_to_fire;                 /* (kgC/m2/d) */
	/* phenology fluxes from transfer pool */
	double leafc_transfer_to_leafc;          /* (kgC/m2/d) */
	double frootc_transfer_to_frootc;        /* (kgC/m2/d) */
	double livestemc_transfer_to_livestemc;  /* (kgC/m2/d) */
	double deadstemc_transfer_to_deadstemc;  /* (kgC/m2/d) */
	double livecrootc_transfer_to_livecrootc;/* (kgC/m2/d) */
	double deadcrootc_transfer_to_deadcrootc;/* (kgC/m2/d) */
	/* leaf and fine root litterfall */
	double leafc_to_litr1c;              /* (kgC/m2/d) */
	double leafc_to_litr2c;              /* (kgC/m2/d) */
	double leafc_to_litr3c;              /* (kgC/m2/d) */
	double leafc_to_litr4c;              /* (kgC/m2/d) */
	double frootc_to_litr1c;             /* (kgC/m2/d) */
	double frootc_to_litr2c;             /* (kgC/m2/d) */
	double frootc_to_litr3c;             /* (kgC/m2/d) */
	double frootc_to_litr4c;             /* (kgC/m2/d) */
	/* maintenance respiration fluxes */
	double leaf_day_mr;                  /* (kgC/m2/d) */
	double leaf_night_mr;                /* (kgC/m2/d) */
	double froot_mr;                     /* (kgC/m2/d) */
	double livestem_mr;                  /* (kgC/m2/d) */
	double livecroot_mr;                 /* (kgC/m2/d) */
	/* photosynthesis flux */
	double psnsun_to_cpool;              /* (kgC/m2/d) */
	double psnshade_to_cpool;            /* (kgC/m2/d) */
	/* litter decomposition fluxes */
	double cwdc_to_litr2c;               /* (kgC/m2/d) */
	double cwdc_to_litr3c;               /* (kgC/m2/d) */
	double cwdc_to_litr4c;               /* (kgC/m2/d) */
	double litr1_hr;                     /* (kgC/m2/d) */
	double litr1c_to_soil1c;             /* (kgC/m2/d) */
	double litr2_hr;                     /* (kgC/m2/d) */
	double litr2c_to_soil2c;             /* (kgC/m2/d) */
	double litr3c_to_litr2c;             /* (kgC/m2/d) */
	double litr4_hr;                     /* (kgC/m2/d) */
	double litr4c_to_soil3c;             /* (kgC/m2/d) */
	double soil1_hr;                     /* (kgC/m2/d) */
	double soil1c_to_soil2c;             /* (kgC/m2/d) */
	double soil2_hr;                     /* (kgC/m2/d) */
	double soil2c_to_soil3c;             /* (kgC/m2/d) */
	double soil3_hr;                     /* (kgC/m2/d) */
	double soil3c_to_soil4c;             /* (kgC/m2/d) */
	double soil4_hr;                     /* (kgC/m2/d) */
	/* daily allocation fluxes from current GPP */
	double cpool_to_leafc;               /* (kgC/m2/d) */
	double cpool_to_leafc_storage;       /* (kgC/m2/d) */
	double cpool_to_frootc;              /* (kgC/m2/d) */
	double cpool_to_frootc_storage;      /* (kgC/m2/d) */
	double cpool_to_livestemc;           /* (kgC/m2/d) */
	double cpool_to_livestemc_storage;   /* (kgC/m2/d) */
	double cpool_to_deadstemc;           /* (kgC/m2/d) */
	double cpool_to_deadstemc_storage;   /* (kgC/m2/d) */
	double cpool_to_livecrootc;          /* (kgC/m2/d) */
	double cpool_to_livecrootc_storage;  /* (kgC/m2/d) */
	double cpool_to_deadcrootc;          /* (kgC/m2/d) */
	double cpool_to_deadcrootc_storage;  /* (kgC/m2/d) */
	double cpool_to_gresp_storage;       /* (kgC/m2/d) */
	/* daily growth respiration fluxes */
	double cpool_leaf_gr;                /* (kgC/m2/d) */
	double cpool_leaf_storage_gr;        /* (kgC/m2/d) */
	double transfer_leaf_gr;             /* (kgC/m2/d) */
	double cpool_froot_gr;               /* (kgC/m2/d) */
	double cpool_froot_storage_gr;       /* (kgC/m2/d) */
	double transfer_froot_gr;            /* (kgC/m2/d) */
	double cpool_livestem_gr;            /* (kgC/m2/d) */
	double cpool_livestem_storage_gr;    /* (kgC/m2/d) */
	double transfer_livestem_gr;         /* (kgC/m2/d) */
	double cpool_deadstem_gr;            /* (kgC/m2/d) */
	double cpool_deadstem_storage_gr;    /* (kgC/m2/d) */
	double transfer_deadstem_gr;         /* (kgC/m2/d) */
	double cpool_livecroot_gr;           /* (kgC/m2/d) */
	double cpool_livecroot_storage_gr;   /* (kgC/m2/d) */
	double transfer_livecroot_gr;        /* (kgC/m2/d) */
	double cpool_deadcroot_gr;           /* (kgC/m2/d) */
	double cpool_deadcroot_storage_gr;   /* (kgC/m2/d) */
	double transfer_deadcroot_gr;        /* (kgC/m2/d) */
	/* annual turnover of storage to transfer pools */
	double leafc_storage_to_leafc_transfer;           /* (kgC/m2/d) */
	double frootc_storage_to_frootc_transfer;         /* (kgC/m2/d) */
	double livestemc_storage_to_livestemc_transfer;    /* (kgC/m2/d) */
	double deadstemc_storage_to_deadstemc_transfer;    /* (kgC/m2/d) */
	double livecrootc_storage_to_livecrootc_transfer; /* (kgC/m2/d) */
	double deadcrootc_storage_to_deadcrootc_transfer; /* (kgC/m2/d) */
	double gresp_storage_to_gresp_transfer;           /* (kgC/m2/d) */
	/* turnover of live wood to dead wood */
	double livestemc_to_deadstemc;        /* (kgC/m2/d) */
	double livecrootc_to_deadcrootc;      /* (kgC/m2/d) */
} cflux_struct;

/* nitrogen state variables (including sums for sources and sinks) */ 
typedef struct
{
    double leafn;              /* (kgN/m2) leaf N */
    double leafn_storage;      /* (kgN/m2) leaf N */
    double leafn_transfer;     /* (kgN/m2) leaf N */
    double frootn;             /* (kgN/m2) fine root N */ 
    double frootn_storage;     /* (kgN/m2) fine root N */ 
    double frootn_transfer;    /* (kgN/m2) fine root N */ 
    double livestemn;          /* (kgN/m2) live stem N */
    double livestemn_storage;  /* (kgN/m2) live stem N */
    double livestemn_transfer; /* (kgN/m2) live stem N */
    double deadstemn;          /* (kgN/m2) dead stem N */
    double deadstemn_storage;  /* (kgN/m2) dead stem N */
    double deadstemn_transfer; /* (kgN/m2) dead stem N */
    double livecrootn;         /* (kgN/m2) live coarse root N */
    double livecrootn_storage; /* (kgN/m2) live coarse root N */
    double livecrootn_transfer;/* (kgN/m2) live coarse root N */
    double deadcrootn;         /* (kgN/m2) dead coarse root N */
    double deadcrootn_storage; /* (kgN/m2) dead coarse root N */
    double deadcrootn_transfer;/* (kgN/m2) dead coarse root N */
    double cwdn;               /* (kgN/m2) coarse woody debris N */
    double litr1n;             /* (kgN/m2) litter labile N */
    double litr2n;             /* (kgN/m2) litter unshielded cellulose N */
    double litr3n;             /* (kgN/m2) litter shielded cellulose N */
    double litr4n;             /* (kgN/m2) litter lignin N */
    double soil1n;             /* (kgN/m2) microbial recycling pool N (fast) */
    double soil2n;             /* (kgN/m2) microbial recycling pool N (medium) */
    double soil3n;             /* (kgN/m2) microbial recycling pool N (slow) */
    double soil4n;             /* (kgN/m2) recalcitrant SOM N (humus, slowest) */
    double sminn;              /* (kgN/m2) soil mineral N */
	double retransn;           /* (kgN/m2) plant pool of retranslocated N */
	double npool;              /* (kgN/m2) temporary plant N pool */
    double nfix_src;           /* (kgN/m2) SUM of biological N fixation */
    double ndep_src;           /* (kgN/m2) SUM of N deposition inputs */
    double nleached_snk;       /* (kgN/m2) SUM of N leached */
    double nvol_snk;           /* (kgN/m2) SUM of N lost to volatilization */
	double fire_snk;           /* (kgN/m2) SUM of N lost to fire */
} nstate_struct;

/* daily nitrogen flux variables */
typedef struct
{
	/* mortality fluxes */
	double m_leafn_to_litr1n;              /* (kgN/m2/d) */
	double m_leafn_to_litr2n;              /* (kgN/m2/d) */
	double m_leafn_to_litr3n;              /* (kgN/m2/d) */
	double m_leafn_to_litr4n;              /* (kgN/m2/d) */
	double m_frootn_to_litr1n;             /* (kgN/m2/d) */
	double m_frootn_to_litr2n;             /* (kgN/m2/d) */
	double m_frootn_to_litr3n;             /* (kgN/m2/d) */
	double m_frootn_to_litr4n;             /* (kgN/m2/d) */
	double m_leafn_storage_to_litr1n;      /* (kgN/m2/d) */
	double m_frootn_storage_to_litr1n;     /* (kgN/m2/d) */
	double m_livestemn_storage_to_litr1n;  /* (kgN/m2/d) */
	double m_deadstemn_storage_to_litr1n;  /* (kgN/m2/d) */
	double m_livecrootn_storage_to_litr1n; /* (kgN/m2/d) */
	double m_deadcrootn_storage_to_litr1n; /* (kgN/m2/d) */
	double m_leafn_transfer_to_litr1n;     /* (kgN/m2/d) */
	double m_frootn_transfer_to_litr1n;    /* (kgN/m2/d) */
	double m_livestemn_transfer_to_litr1n; /* (kgN/m2/d) */
	double m_deadstemn_transfer_to_litr1n; /* (kgN/m2/d) */
	double m_livecrootn_transfer_to_litr1n;/* (kgN/m2/d) */
	double m_deadcrootn_transfer_to_litr1n;/* (kgN/m2/d) */
	double m_livestemn_to_litr1n;          /* (kgN/m2/d) */
	double m_livestemn_to_cwdn;            /* (kgN/m2/d) */
	double m_deadstemn_to_cwdn;            /* (kgN/m2/d) */
	double m_livecrootn_to_litr1n;         /* (kgN/m2/d) */
	double m_livecrootn_to_cwdn;           /* (kgN/m2/d) */
	double m_deadcrootn_to_cwdn;           /* (kgN/m2/d) */
	double m_retransn_to_litr1n;           /* (kgN/m2/d) */
	/* fire fluxes */
	double m_leafn_to_fire;                /* (kgN/m2/d) */
	double m_frootn_to_fire;               /* (kgN/m2/d) */
	double m_leafn_storage_to_fire;        /* (kgN/m2/d) */
	double m_frootn_storage_to_fire;       /* (kgN/m2/d) */
	double m_livestemn_storage_to_fire;    /* (kgN/m2/d) */
	double m_deadstemn_storage_to_fire;    /* (kgN/m2/d) */
	double m_livecrootn_storage_to_fire;   /* (kgN/m2/d) */
	double m_deadcrootn_storage_to_fire;   /* (kgN/m2/d) */
	double m_leafn_transfer_to_fire;       /* (kgN/m2/d) */
	double m_frootn_transfer_to_fire;      /* (kgN/m2/d) */
	double m_livestemn_transfer_to_fire;   /* (kgN/m2/d) */
	double m_deadstemn_transfer_to_fire;   /* (kgN/m2/d) */
	double m_livecrootn_transfer_to_fire;  /* (kgN/m2/d) */
	double m_deadcrootn_transfer_to_fire;  /* (kgN/m2/d) */
	double m_livestemn_to_fire;            /* (kgN/m2/d) */
	double m_deadstemn_to_fire;            /* (kgN/m2/d) */
	double m_livecrootn_to_fire;           /* (kgN/m2/d) */
	double m_deadcrootn_to_fire;           /* (kgN/m2/d) */
	double m_retransn_to_fire;             /* (kgN/m2/d) */
	double m_litr1n_to_fire;               /* (kgN/m2/d) */
	double m_litr2n_to_fire;               /* (kgN/m2/d) */
	double m_litr3n_to_fire;               /* (kgN/m2/d) */
	double m_litr4n_to_fire;               /* (kgN/m2/d) */
	double m_cwdn_to_fire;                 /* (kgN/m2/d) */
	/* phenology fluxes from transfer pool */
	double leafn_transfer_to_leafn;           /* (kgN/m2/d) */
	double frootn_transfer_to_frootn;         /* (kgN/m2/d) */
	double livestemn_transfer_to_livestemn;   /* (kgN/m2/d) */
	double deadstemn_transfer_to_deadstemn;   /* (kgN/m2/d) */
	double livecrootn_transfer_to_livecrootn; /* (kgN/m2/d) */
	double deadcrootn_transfer_to_deadcrootn; /* (kgN/m2/d) */
	/* litterfall fluxes */
	double leafn_to_litr1n;               /* (kgN/m2/d) */
	double leafn_to_litr2n;               /* (kgN/m2/d) */ 
	double leafn_to_litr3n;               /* (kgN/m2/d) */
	double leafn_to_litr4n;               /* (kgN/m2/d) */
	double leafn_to_retransn;             /* (kgN/m2/d) */
	double frootn_to_litr1n;              /* (kgN/m2/d) */
	double frootn_to_litr2n;              /* (kgN/m2/d) */
	double frootn_to_litr3n;              /* (kgN/m2/d) */
	double frootn_to_litr4n;              /* (kgN/m2/d) */
	/* deposition flux */
	double ndep_to_sminn;                 /* (kgN/m2/d) */
	double nfix_to_sminn;                 /* (kgN/m2/d) */
	/* litter and soil decomposition fluxes */
	double cwdn_to_litr2n;                /* (kgN/m2/d) */
	double cwdn_to_litr3n;                /* (kgN/m2/d) */
	double cwdn_to_litr4n;                /* (kgN/m2/d) */
	double litr1n_to_soil1n;              /* (kgN/m2/d) */
	double sminn_to_soil1n_l1;            /* (kgN/m2/d) */
	double litr2n_to_soil2n;              /* (kgN/m2/d) */
	double sminn_to_soil2n_l2;            /* (kgN/m2/d) */
	double litr3n_to_litr2n;              /* (kgN/m2/d) */
	double litr4n_to_soil3n;              /* (kgN/m2/d) */
	double sminn_to_soil3n_l4;            /* (kgN/m2/d) */
	double soil1n_to_soil2n;              /* (kgN/m2/d) */
	double sminn_to_soil2n_s1;            /* (kgN/m2/d) */
	double soil2n_to_soil3n;              /* (kgN/m2/d) */
	double sminn_to_soil3n_s2;            /* (kgN/m2/d) */
	double soil3n_to_soil4n;              /* (kgN/m2/d) */
	double sminn_to_soil4n_s3;            /* (kgN/m2/d) */
	double soil4n_to_sminn;               /* (kgN/m2/d) */
	/* denitrification (volatilization) fluxes */
	double sminn_to_nvol_l1s1;            /* (kgN/m2/d) */
	double sminn_to_nvol_l2s2;            /* (kgN/m2/d) */
	double sminn_to_nvol_l4s3;            /* (kgN/m2/d) */
	double sminn_to_nvol_s1s2;            /* (kgN/m2/d) */
	double sminn_to_nvol_s2s3;            /* (kgN/m2/d) */
	double sminn_to_nvol_s3s4;            /* (kgN/m2/d) */
	double sminn_to_nvol_s4;              /* (kgN/m2/d) */
	double sminn_to_denitrif;             /* (kgN/m2/d) */
	
	/* leaching flux */
	double sminn_leached;                 /* (kgN/m2/d) */
	/* daily allocation fluxes */
	double retransn_to_npool;             /* (kgN/m2/d) */
	double sminn_to_npool;                /* (kgN/m2/d) */
	double npool_to_leafn;                /* (kgN/m2/d) */
	double npool_to_leafn_storage;        /* (kgN/m2/d) */
	double npool_to_frootn;               /* (kgN/m2/d) */
	double npool_to_frootn_storage;       /* (kgN/m2/d) */
	double npool_to_livestemn;            /* (kgN/m2/d) */
	double npool_to_livestemn_storage;    /* (kgN/m2/d) */
	double npool_to_deadstemn;            /* (kgN/m2/d) */
	double npool_to_deadstemn_storage;    /* (kgN/m2/d) */
	double npool_to_livecrootn;           /* (kgN/m2/d) */
	double npool_to_livecrootn_storage;   /* (kgN/m2/d) */
	double npool_to_deadcrootn;           /* (kgN/m2/d) */
	double npool_to_deadcrootn_storage;   /* (kgN/m2/d) */
	/* annual turnover of storage to transfer */
	double leafn_storage_to_leafn_transfer;           /* (kgN/m2/d) */
	double frootn_storage_to_frootn_transfer;         /* (kgN/m2/d) */
	double livestemn_storage_to_livestemn_transfer;   /* (kgN/m2/d) */
	double deadstemn_storage_to_deadstemn_transfer;   /* (kgN/m2/d) */
	double livecrootn_storage_to_livecrootn_transfer; /* (kgN/m2/d) */
	double deadcrootn_storage_to_deadcrootn_transfer; /* (kgN/m2/d) */
	/* turnover of live wood to dead wood, with retranslocation */
	double livestemn_to_deadstemn;        /* (kgN/m2/d) */
	double livestemn_to_retransn;         /* (kgN/m2/d) */
	double livecrootn_to_deadcrootn;      /* (kgN/m2/d) */
	double livecrootn_to_retransn;        /* (kgN/m2/d) */
} nflux_struct;

/* temporary nitrogen variables for reconciliation of decomposition
immobilization fluxes and plant growth N demands */
typedef struct 
{
	double mineralized;
	double potential_immob;
	double plitr1c_loss;
	double pmnf_l1s1;
	double plitr2c_loss;
	double pmnf_l2s2;
	double plitr4c_loss;
	double pmnf_l4s3;
	double psoil1c_loss;
	double pmnf_s1s2;
	double psoil2c_loss;
	double pmnf_s2s3;
	double psoil3c_loss;
	double pmnf_s3s4;
	double psoil4c_loss;
	double kl4;
} ntemp_struct;
	
/* phenological control arrays */
typedef struct
{
	int* remdays_curgrowth; /* (nmetdays) days left in current growth season */
	int* remdays_transfer;  /* (nmetdays) number of transfer days remaining */
	int* remdays_litfall;   /* (nmetdays) number of litfall days remaining */
	int* predays_transfer;  /* (nmetdays) number of transfer days previous */
	int* predays_litfall;   /* (nmetdays) number of litfall days previous */

	/* for gsi */
	int** onday_arr;			/* (doy) ARRAY of first day of transfer period [year;yearday]*/
	int** offday_arr;		    /* (doy) ARRAY of last day of transfer period [year;yearday]*/
} phenarray_struct;

/* daily phenological data array */
typedef struct
{
	double remdays_curgrowth; /* days left in current growth season */
	double remdays_transfer;  /* number of transfer days remaining */
	double remdays_litfall;   /* number of litfall days remaining */
	double predays_transfer;  /* number of transfer days previous */
	double predays_litfall;   /* number of litfall days previous */
} phenology_struct;

/* ecophysiological variables */
typedef struct
{
	double day_leafc_litfall_increment;      /* (kgC/m2/d) rate leaf litfall */
	double day_frootc_litfall_increment;     /* (kgC/m2/d) rate froot litfall */
	double day_livestemc_turnover_increment; /* (kgC/m2/d) rate livestem turnover */
	double day_livecrootc_turnover_increment;/* (kgC/m2/d) rate livecroot turnover */
	double annmax_leafc;      /* (kgC/m2) annual maximum daily leaf C */
	double annmax_frootc;     /* (kgC/m2) annual maximum daily froot C */
	double annmax_livestemc;  /* (kgC/m2) annual maximum daily livestem C */
	double annmax_livecrootc; /* (kgC/m2) annual maximum daily livecroot C */
	double dsr;            /* (days) number of days since rain, for soil evap */
    double proj_lai;       /* (DIM) live projected leaf area index */
    double all_lai;        /* (DIM) live all-sided leaf area index */
    double plaisun;        /* (DIM) sunlit projected leaf area index */
    double plaishade;      /* (DIM) shaded projected leaf area index */
    double sun_proj_sla;   /* (m2/kgC) sunlit projected SLA */
    double shade_proj_sla; /* (m2/kgC) shaded projected SLA */
    double psi;            /* (MPa) water potential of soil and leaves */
    double vwc;            /* (DIM) volumetric water content */
    double dlmr_area_sun;  /* (umolC/m2projected leaf area/s) sunlit leaf MR */
    double dlmr_area_shade;/* (umolC/m2projected leaf area/s) shaded leaf MR */
    double gl_t_wv_sun;   /* (m/s) leaf-scale conductance to transpired water */
    double gl_t_wv_shade; /* (m/s) leaf-scale conductance to transpired water */
    double assim_sun;      /* (umol/m2/s) sunlit assimilation per unit pLAI */
    double assim_shade;    /* (umol/m2/s) shaded assimilation per unit pLAI */
    /* decomp variables */
    double t_scalar;       /* (DIM) decomp temperature scalar */
    double w_scalar;       /* (DIM) decomp water scalar */
    double rate_scalar;    /* (DIM) decomp combined scalar */
	double daily_gross_nmin; /* (kgN/m2/d) daily gross N mineralization */
	double daily_gross_nimmob; /* (kgN/m2/d) daily gross N immobilization */ 
	double daily_net_nmin; /* (kgN/m2/d) daily net N mineralization */
	double fpi;            /* (DIM) fraction of potential immobilization */

    /* the following are optional outputs, usually set if the appropriate
    functions are called with the flag verbose = 1 */
	double m_tmin;         /* (DIM) freezing night temperature multiplier */
	double m_psi;          /* (DIM) water potential multiplier */
	double m_co2;          /* (DIM) atmospheric [CO2] multiplier */
	double m_ppfd_sun;         /* (DIM) PAR flux density multiplier */
	double m_ppfd_shade;         /* (DIM) PAR flux density multiplier */
	double m_vpd;          /* (DIM) vapor pressure deficit multiplier */
	double m_final_sun;        /* (DIM) product of all other multipliers */
	double m_final_shade;        /* (DIM) product of all other multipliers */
	double gl_bl;          /* (m/s) leaf boundary layer conductance */
    double gl_c;           /* (m/s) leaf cuticular conductance */
    double gl_s_sun;           /* (m/s) leaf-scale stomatal conductance */
    double gl_s_shade;           /* (m/s) leaf-scale stomatal conductance */
	double gl_e_wv;        /* (m/s) leaf conductance to evaporated water */
	double gl_sh;          /* (m/s) leaf conductance to sensible heat */
    double gc_e_wv;        /* (m/s) canopy conductance to evaporated water */
    double gc_sh;          /* (m/s) canopy conductance to sensible heat */
	
	/* diagnostic variables for ouput purposes only */
	double ytd_maxplai;    /* (DIM) year-to-date maximum projected LAI */
} epvar_struct;

/* soil and site constants */
typedef struct
{
    double soil_depth;     /* (m)   effective depth of rooting zone */
    double soil_b;         /* (DIM) Clapp-Hornberger "b" parameter */
    double vwc_sat;        /* (DIM) volumetric water content at saturation */
    double vwc_fc;         /* (DIM) VWC at field capacity ( = -0.015 MPa) */
    double psi_sat;        /* (MPa) soil matric potential at saturation */
    double soilw_sat;      /* (kgH2O/m2) soilwater at saturation */
    double soilw_fc;       /* (kgH2O/m2) soilwater at field capacity */
    double elev;           /* (m) site elevation */
    double lat;            /* (degrees) site latitude (negative for south) */
    double sw_alb;         /* (DIM) surface shortwave albedo */
    double ndep;           /* (kgN/m2/yr) wet+dry atmospheric deposition of N */
    double nfix;           /* (kgN/m2/yr) symbiotic+asymbiotic fixation of N */
	double soiltype;       /* (DIM) soil type*/
	double epctype;        /* (DIM) vegetation type*/   
	double epctype_change; /* (DIM) vegetation type which were changed at a random day between a certain period 2013-06-26*/
	char epcfilepath[100];
	long int climate_id;
	long int pointcounter;
	double sand;
	double clay;
	double silt;
	long int plotID;
	
} siteconst_struct;

/* canopy ecophysiological constants */
typedef struct
{
	int woody;             /* (flag) 1=woody, 0=non-woody */
	int evergreen;         /* (flag) 1=evergreen, 0=deciduous */
	int c3_flag;           /* (flag) 1 = C3,  0 = C4 */
	int phenology_flag;    /* (flag) 1=phenology model, 0=user defined */
	int onday;             /* (yday) yearday leaves on */
	int offday;            /* (yday) yearday leaves off */
	double transfer_pdays;    /* (prop.) fraction of growth period for transfer */
	double litfall_pdays;     /* (prop.) fraction of growth period for litfall */
    double leaf_turnover;     /* (1/yr) annual leaf turnover fraction */
    double froot_turnover;    /* (1/yr) annual fine root turnover fraction */
    double livewood_turnover; /* (1/yr) annual live wood turnover fraction */
    double daily_mortality_turnover; /* (1/day) daily mortality turnover */
    double daily_fire_turnover;      /* (1/day) daily fire turnover */
    double alloc_frootc_leafc; /* (ratio) new fine root C to new leaf C */
	double alloc_newstemc_newleafc; /* (ratio) new stem C to new leaf C */
	double alloc_newlivewoodc_newwoodc; /* (ratio) new livewood C:new wood C */
    double alloc_crootc_stemc; /* (ratio) new live croot C to new live stem C */
	/*2014-11-20添加新的分配系数，树根和树叶比*/
	double alloc_crootc_leafc;
    double alloc_prop_curgrowth; /* (prop.) daily allocation to current growth */
    double avg_proj_sla;   /* (m2/kgC) canopy average proj. SLA */
    double sla_ratio;      /* (DIM) ratio of shaded to sunlit projected SLA */
    double lai_ratio;      /* (DIM) ratio of (all-sided LA / one-sided LA) */
    double int_coef;       /* (kg/kg/LAI/d) canopy precip interception coef */
    double ext_coef;       /* (DIM) canopy light extinction coefficient */
    double flnr;           /* (kg NRub/kg Nleaf) leaf N in Rubisco */
    double psi_open;       /* (MPa) psi at start of conductance reduction */
    double psi_close;      /* (MPa) psi at complete conductance reduction */
    double vpd_open;       /* (Pa)  vpd at start of conductance reduction */
	double vpd_close;      /* (Pa)  vpd at complete conductance reduction */
    double gl_smax;        /* (m/s) maximum leaf-scale stomatal conductance */
    double gl_c;           /* (m/s) leaf-scale cuticular conductance */
	double gl_bl;          /* (m/s) leaf-scale boundary layer conductance */
    double froot_cn;     /* (kgC/kgN) C:N for fine roots */
    double leaf_cn;      /* (kgC/kgN) C:N for leaves */
    double livewood_cn;  /* (kgC/kgN) C:N for live wood */
	double deadwood_cn;  /* (kgC/kgN) C:N for dead wood */
    double leaflitr_cn;      /* (kgC/kgN) constant C:N for leaf litter */
    double leaflitr_flab;    /* (DIM) leaf litter labile fraction */
    double leaflitr_fucel;   /* (DIM) leaf litter unshielded cellulose fract. */
    double leaflitr_fscel;   /* (DIM) leaf litter shielded cellulose fract. */
    double leaflitr_flig;    /* (DIM) leaf litter lignin fraction */
    double frootlitr_flab;   /* (DIM) froot litter labile fraction */
    double frootlitr_fucel;  /* (DIM) froot litter unshielded cellulose fract */
    double frootlitr_fscel;  /* (DIM) froot litter shielded cellulose fract */
    double frootlitr_flig;   /* (DIM) froot litter lignin fraction */
    double deadwood_fucel;   /* (DIM) dead wood unshileded cellulose fraction */
    double deadwood_fscel;   /* (DIM) dead wood shielded cellulose fraction */
    double deadwood_flig;    /* (DIM) dead wood lignin fraction */
	double decom_KS4;		 /* (DIM) 调整不同树种的分解速率，腐殖质*/
	double decom_KS3;		 /* (DIM) 调整不同树种的分解速率，主要对难分解碳库*/
	double decom_cwdc;		 /* (DIM) cwdc*/
	double decom_KL2;
	double decom_KL4;
	double epctype;       /* (DIM) epctype*/
	double epctype_change;        /* (DIM) vegetation type*/  

	/* for gsi */
	double snowcover_limit;		/* critical amount of snow (above: no vegetation period) */
	double heatsum_limit1;			/* lower limit of heatsum to calculate heatsum index */
	double heatsum_limit2;			/* upper limit of heatsum to calculate heatsum index */
	double tmin_limit1;				/* lower limit of tmin to calculate tmin index */
	double tmin_limit2;				/* upper limit of tmin to calculate tmin index  */
	double vpd_limit1;				/* lower limit of tmin to calculate vpd index  */
	double vpd_limit2;				/* upper limit of tmin to calculate vpd index  */
	double dayl_limit1;				/* lower limit of dayl to calculate vpd index  */
	double dayl_limit2;				/* upper limit of dayl to calculate vpd index  */
	int n_moving_avg;				/* moving average (calculated from indicatiors to avoid the effects of single extreme events)*/
	double GSI_limit_SGS;			/* when GSI first time greater that limit -> start of growing season */
	double GSI_limit_EGS;			/* when GSI first time less that limit -> end of growing season */
	/* end for gsi */

} epconst_struct;

/* structure for the photosynthesis routine */
typedef struct
{
	int c3;                 /* (flag) set to 1 for C3 model, 0 for C4 model */
	double pa;              /* (Pa) atmospheric pressure */
	double co2;             /* (ppm) atmospheric [CO2] */
	double t;               /* (deg C) temperature */
	double lnc;             /* (kg Nleaf/m2) leaf N per unit sunlit leaf area */
	double flnr;            /* (kg NRub/kg Nleaf) fract. of leaf N in Rubisco */
	double ppfd;            /* (umol/m2/s) PAR flux per unit sunlit leaf area */
	double g;               /* (umol/m2/s/Pa) conductance to CO2 */
	double dlmr;            /* (umol/m2/s) day leaf m. resp, proj. area basis */
	double Ci;              /* (Pa) intercellular [CO2] */
	double O2;              /* (Pa) atmospheric [O2] */
	double Ca;              /* (Pa) atmospheric [CO2] */
	double gamma;           /* (Pa) CO2 compensation point, no Rd */
	double Kc;              /* (Pa) MM constant carboxylation */
	double Ko;              /* (Pa) MM constant oxygenation */
	double Vmax;            /* (umol/m2/s) max rate carboxylation */
	double Jmax;            /* (umol/m2/s) max rate electron transport */
	double J;               /* (umol/m2/s) rate of RuBP regeneration */
	double Av;              /* (umol/m2/s) carboxylation limited assimilation */
	double Aj;              /* (umol/m2/s) RuBP regen limited assimilation */
	double A;               /* (umol/m2/s) final assimilation rate */
} psn_struct;


typedef struct
{
	double daily_npp;      /* kgC/m2/day = GPP - Rmaint - Rgrowth */
	double daily_nep;      /* kgC/m2/day = NPP - Rheterotroph */
	double daily_nee;      /* kgC/m2/day = NEP - fire losses */
	double daily_gpp;      /* kgC/m2/day  gross PSN source */
	double daily_mr;       /* kgC/m2/day  maintenance respiration */
	double daily_gr;       /* kgC/m2/day  growth respiration */
	double daily_hr;       /* kgC/m2/day  heterotrophic respiration */
	double daily_fire;     /* kgC/m2/day  fire losses */
	double daily_litfallc; /* kgC/m2/day  total litterfall */
	double daily_et;    /* kgW/m2/day daily evapotranspiration */
	double daily_evap;    /* kgW/m2/day daily evaporation */
	double daily_trans;   /* kgW/m2/day daily transpiration */
	double daily_outflow; /* kgW/m2/day daily outflow */
	double daily_soilw; /* kgW/m2/day daily soilw */
	double daily_snoww; /* kgW/m2/day daily snoww */
	double cum_npp;        /* kgC/m2  Summed over entire simulation */
	double cum_nep;        /* kgC/m2  Summed over entire simulation */
	double cum_nee;        /* kgC/m2  Summed over entire simulation */
	double cum_gpp;        /* kgC/m2  Summed over entire simulation */
	double cum_mr;         /* kgC/m2  Summed over entire simulation */
	double cum_gr;         /* kgC/m2  Summed over entire simulation */
	double cum_hr;         /* kgC/m2  Summed over entire simulation */
	double cum_fire;       /* kgC/m2  Summed over entire simulation */
	double vegc;           /* kgC/m2  total vegetation C */
	double litrc;          /* kgC/m2  total litter C */
	double soilc;          /* kgC/m2  total soil C */
	double totalc;         /* kgC/m2  total of vegc, litrc, and soilc */
	double leafc_sum;		   /*KgC/m2 leaf C 2013-03-27*/
	double stemc_sum;		   /*KgC/m2 stem C (include branch) 2013-03-27*/
	double daily_root_mr;	   /*KgC/m2 fine root and coarse root respiration 2013-07-21*/
	double daily_leaf_mr;	   /*KgC/m2 fine root and coarse root respiration 2013-08-12*/
	double crootc_sum;
	double leafc_max;
	double vegc_max;   //用于农田的生物量输出
	double crootc_max;//用于农田的生物量输出
	double frootc_max;//用于农田的生物量输出
	double litrc_max;
	double soilc_max;
	double totalc_max;
} summary_struct;

/* restart data structure */
typedef struct
{
	double soilw;
	double snoww;
	double canopyw;
	double leafc;
	double leafc_storage;
	double leafc_transfer;
	double frootc;
	double frootc_storage;
	double frootc_transfer;
	double livestemc;
	double livestemc_storage;
	double livestemc_transfer;
	double deadstemc;
	double deadstemc_storage;
	double deadstemc_transfer;
	double livecrootc;
	double livecrootc_storage;
	double livecrootc_transfer;
	double deadcrootc;
	double deadcrootc_storage;
	double deadcrootc_transfer;
	double gresp_storage;
	double gresp_transfer;
	double cwdc;
	double litr1c;
	double litr2c;
	double litr3c;
	double litr4c;
	double soil1c;
	double soil2c;
	double soil3c;
	double soil4c;
	double cpool;
	double leafn;
	double leafn_storage;
	double leafn_transfer;
	double frootn;
	double frootn_storage;
	double frootn_transfer;
	double livestemn;
	double livestemn_storage;
	double livestemn_transfer;
	double deadstemn;
	double deadstemn_storage;
	double deadstemn_transfer;
	double livecrootn;
	double livecrootn_storage;
	double livecrootn_transfer;
	double deadcrootn;
	double deadcrootn_storage;
	double deadcrootn_transfer;
	double cwdn;
	double litr1n;
	double litr2n;
	double litr3n;
	double litr4n;
	double soil1n;
	double soil2n;
	double soil3n;
	double soil4n;
	double sminn;
	double retransn;
	double npool;
	double day_leafc_litfall_increment;
	double day_frootc_litfall_increment;
	double day_livestemc_turnover_increment;
	double day_livecrootc_turnover_increment;
	double annmax_leafc;
	double annmax_frootc;
	double annmax_livestemc;
	double annmax_livecrootc;
	double dsr;
	int metyr;
} restart_data_struct;

/* add high time resolution data */
typedef struct
{
	int YEAR, DAY, HRMIN;  //time year, day of year, minutes
	double Rd;  //(umol / m2 / s)
	double GPP;	//(umol / m2 / s)
	double NEE;	//(umol / m2 / s)
	double TA;	//(deg C)
	double PREC;  //(mm)
	double RH;	//(%)
	double VPD;	//(kPa)
	double SWC;	//(%)
	double PAR;	//(umol / m2 / s)
	double ppfd;
} stationFlux_data_struct;

typedef struct
{
	bool output_stress;	// output stress
	bool output_old_cpool;
	bool active;		// main control of hightime
	bool output_carbon;	// output hightime carbon
	bool tempCorr; // high temp correction

	char * stationFile;
	std::ofstream tmpHighFile;
	std::ofstream highFile_carbon;
	std::ofstream highFile_stress;

	char * tempCorrFile;
	std::ofstream tempCorrStream;

}high_time_resolution;

typedef struct
{
	bool activate;	// 
	double a;
	double b;
}pymc;

typedef struct
{
	bool active;
	char * laiFile;
}
lai_model;

typedef struct
{
	bool active;
	char * gsiFile;
}
gsi_model;


#ifdef __cplusplus
}
#endif

#endif
