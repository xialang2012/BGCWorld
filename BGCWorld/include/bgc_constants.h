#ifndef BGC_CONSTANTS_H
#define BGC_CONSTANTS_H
/*
bgc_constants.h
Holds macro definitions for constants used in bgc()

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions from version 4.1.1:
Moved the heterotrophic respiration fractions and base decomposition rates 
out of the science subroutines and into the constants file so that it is
easier to make sure modifications are propagated correctly.
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#ifdef __cplusplus
extern "C"
{
#endif



/**************************************add***************************************/
/*[added */
#define FSM_e       2.71828
#define Liv      2.844e6        // latent heat of sublimation J/kg
#define FSM_g        9.8           //gravity  acceleration
#define r        66              // (pa/K)humiditor constant 
/*added 2014-09-26]*/
#define N        20             /* the number of layers of soil */ /*add 2014-07-30*/
#define Lf       3.336e5        /*  Latent heat of freezing J/kg */
#define Tk       273.15         // absolute zero temperature
#define Lvap     2.501e6        // latent heat of vaporization J/kg 
#define rou_w    1000           //Density of water kg/m3
#define rou_i    917            //Density of ice Kg/m3
//#define rou_s    1500           //Density of solid soil kg/m3
   
#define lamta_i  2.1597         //Conductivity of ice
#define lamta_w  0.5504         //Conductivity of water 0.58W/m.K
#define lamta_a  0.0209         //Conductivity of air
#define C_w      4.15e6         //Specific heat of water volume j/m3k
#define C_i      1.9228e6       //Specific heat of ice volume j/m3k
#define C_a      1.55818e3      //specific heat of air volume j/m3k

#define Tk       273.15         //absolution temprature 0C
#define Ck       8               //the specific surface of the prsence of ice
#define Hr        66              /* (pa)humiditor constant */
#define Rw       461.5           /* (J/kg.K) gas law constant of vapor */
/**************************************add***************************************/
/* atmospheric constants */
/* from the definition of the standard atmosphere, as established
by the International Civil Aviation Organization, and referenced in:

Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics. 2nd 
	Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
	(pp 10,167-168,245)
*/
#define G_STD    9.80665         /* (m/s2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define MA       28.9644e-3      /* (kg/mol) molecular weight of air */
#define MW       18.0148e-3      /* (kg/mol) molecular weight of water */
#define CP       1010.0          /* (J/kg K) specific heat of air */
#define LR_STD   0.0065          /* (-K/m) standard temperature lapse rate */
#define R        8.3143          /* (m3 Pa/ mol K) gas law constant */
#define SBC      5.67e-8         /* (W/(m2 K4)) Stefan-Boltzmann constant */
#define EPS      0.6219          /* (MW/MA) unitless ratio of molec weights */


/* ecosystem constants */
#define RAD2PAR     0.45     /* (DIM) ratio PAR / SWtotal  */
#define EPAR        4.55     /* (umol/J) PAR photon energy ratio */  
#define SOIL1_CN    12.0     /* C:N for fast microbial recycling pool */
#define SOIL2_CN    12.0     /* C:N for slow microbial recycling pool */
#define SOIL3_CN    10.0     /* C:N for recalcitrant SOM pool (humus) */
#define SOIL4_CN    10.0     /* C:N for recalcitrant SOM pool (humus) */
#define GRPERC      0.3      /* (DIM) growth resp per unit of C grown */
#define GRPNOW      1.0      /* (DIM) proportion of storage growth resp at fixation */
#define PPFD50      75.0     /* (umol/m2/s) PPFD for 1/2 stomatal closure */
#define DENITRIF_PROPORTION  0.01  /* fraction of mineralization to volatile */
#define MOBILEN_PROPORTION   0.1   /* fraction mineral N avail for leaching */

/* use this block of constants to include the dynamics for slowest soil pool (s4) */
/* respiration fractions for fluxes between compartments (unitless) */ 
//#define	RFL1S1		0.39	/* transfer from litter 1 to soil 1 */
//#define	RFL2S2		0.55	/* transfer from litter 2 to soil 2 */
//#define	RFL4S3		0.29	/* transfer from litter 4 to soil 3 */
//#define	RFS1S2		0.28	/* transfer from soil 1 to soil 2 */
//#define	RFS2S3		0.46    /* transfer from soil 2 to soil 3 */
//#define	RFS3S4		0.55	/* transfer from soil 3 to soil 4 */
//减少分解过程中的呼吸消耗 2013-07-19
#define	RFL1S1		0.39	/* transfer from litter 1 to soil 1 */
#define	RFL2S2		0.55	/* transfer from litter 2 to soil 2 */
#define	RFL4S3		0.29	/* transfer from litter 4 to soil 3 */
#define	RFS1S2		0.28	/* transfer from soil 1 to soil 2 */
#define	RFS2S3		0.46    /* transfer from soil 2 to soil 3 */
#define	RFS3S4		0.55	/* transfer from soil 3 to soil 4 */
/*默认 base decomposition rate constants (1/day) */ 
//#define KL1_BASE	0.7		/* labile litter pool */
//#define KL2_BASE	0.07	/* cellulose litter pool */
//#define KL4_BASE	0.014	/* lignin litter pool */
//#define KS1_BASE	0.07	/* fast microbial recycling pool */
//#define KS2_BASE	0.014	/* medium microbial recycling pool */
//#define KS3_BASE	0.0014	/* slow microbial recycling pool */
//#define KS4_BASE	0.0001	/* recalcitrant SOM (humus) pool */
//#define KFRAG_BASE	0.001	/* physical fragmentation of coarse woody debris */
//降低土壤分解速度，增加粗木质部残体分解速度
/*修改后——第一次修改*/
//#define KL1_BASE	0.7		/* labile litter pool */
//#define KL2_BASE	0.09	/* cellulose litter pool */
//#define KL4_BASE	0.05	/* lignin litter pool */
//#define KS1_BASE	0.07	/* fast microbial recycling pool 2013-07-18  */
//#define KS2_BASE	0.014	/* medium microbial recycling pool 2013-07-18  */
//#define KS3_BASE	0.0007	/* slow microbial recycling pool 2013-07-18 */
//#define KS4_BASE	0.00004	/* recalcitrant SOM (humus) pool 2013-07-13*/
//#define KFRAG_BASE	0.005	/* physical fragmentation of coarse woody debris 2013-07-13*/
/*修改后——第二次修改2013-08-23*/
//#define KL1_BASE	0.7		/* labile litter pool */
//#define KL2_BASE	0.09	/* cellulose litter pool */
//#define KL4_BASE	0.05	/* lignin litter pool */
//#define KS1_BASE	0.07	/* fast microbial recycling pool 2013-08-24  */
//#define KS2_BASE	0.014	/* medium microbial recycling pool 2013-08-24  */
//#define KS3_BASE	0.0007	/* slow microbial recycling pool 2013-08-24 */
//#define KS4_BASE	0.00004	/* recalcitrant SOM (humus) pool 2013-08-24*/
//#define KFRAG_BASE	0.005	/* physical fragmentation of coarse woody debris 2013-08-24*/
/*2014-11-05 青海土壤碳分解速率，根据各个土壤碳库的比例（spinup后的比例），推算出默认参数条件下土壤平均周转为5.3年*/
//#define KL1_BASE	0.7		/* labile litter pool */
//#define KL2_BASE	0.07	/* cellulose litter pool */
//#define KL4_BASE	0.014	/* lignin litter pool */

//#define KS1_BASE	0.0149	/* fast microbial recycling pool */
//#define KS2_BASE	0.00298	/* medium microbial recycling pool */
//#define KS3_BASE	0.000298	/* slow microbial recycling pool */
//#define KS4_BASE	0.0000213	/* recalcitrant SOM (humus) pool */
//#define KFRAG_BASE	0.002	/* physical fragmentation of coarse woody debris */

/*默认 2014-11-06 土壤分解回归默认 (1/day) */ 
#define KL1_BASE	0.7		/* labile litter pool */
#define KL2_BASE	0.07	/* cellulose litter pool */
#define KL4_BASE	0.014	/* lignin litter pool */
#define KS1_BASE	0.07	/* fast microbial recycling pool */
#define KS2_BASE	0.014	/* medium microbial recycling pool */
#define KS3_BASE	0.0014	/* slow microbial recycling pool */
#define KS4_BASE	0.0001	/* recalcitrant SOM (humus) pool */
#define KFRAG_BASE	0.001	/* physical fragmentation of coarse woody debris */



/* use this block of constants to exclude the dynamics for slowest soil pool (s4) */
/* respiration fractions for fluxes between compartments (unitless) */ 
/*#define	RFL1S1	0.39*/	/* transfer from litter 1 to soil 1 */
/*#define	RFL2S2	0.55*/	/* transfer from litter 2 to soil 2 */
/*#define	RFL4S3	0.29*/	/* transfer from litter 4 to soil 3 */
/*#define	RFS1S2	0.28*/	/* transfer from soil 1 to soil 2 */
/*#define	RFS2S3	0.46*/  /* transfer from soil 2 to soil 3 */
/*#define	RFS3S4	1.00*/	/* transfer from soil 3 to soil 4 */
/* base decomposition rate constants (1/day) */ 
/*#define KL1_BASE	0.7	*/	/* labile litter pool */
/*#define KL2_BASE	0.07*/	/* cellulose litter pool */
/*#define KL4_BASE	0.014*/	/* lignin litter pool */
/*#define KS1_BASE	0.07*/	/* fast microbial recycling pool */
/*#define KS2_BASE	0.014*/	/* medium microbial recycling pool */
/*#define KS3_BASE	0.0005*/	/* slow microbial recycling pool */
/*#define KS4_BASE	0.0000*/	/* recalcitrant SOM (humus) pool */

/* precision control */
/* This constant determines the lower limit of state variables before they
are set to 0.0 to control rounding and overflow errors */
#define CRIT_PREC 1e-20

/* spinup control */
/* maximum allowable trend in slow soil carbon at steady-state (kgC/m2/yr) */
#define SPINUP_TOLERANCE 0.0005
#define MODE_INI 0
#define MODE_SPINUP 1
#define MODE_MODEL 2
#define MODE_SPINNGO 3

/* output control constants */
#define NMAP 700

/* For modifying summary output as per pan-arctic bgc */
#define SANE 1
#define INSANE 0

#ifdef __cplusplus
}
#endif

#endif
