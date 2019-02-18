/* 
baresoil_evap.c
daily bare soil evaporation

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
swht-BGC version 1.0 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

int sparameters()//, int typeNo) 
{
	/*
	int typeNo;              soil texture types code numer: 1 sand,2 loamysand,3 sandyloam,4 loam,5 siltyloam,6 sandyclayloam,7 clayloam,8 siltyclayloam,9 sandyclay,10 siltyclay,11 clay,12 peat,13 moss,14 lichen 
	double psandy;           (%)percentage of sand 
	double psilt;           (%)percentage of slit 
	double pclay;            (%)percentage of clay 
	double k_solids;         (Wm-1k-1)the thermal conductivity of the soil solids 
	double c_solids;         (the heat capacity of the soil solids 
	double k_sat;            (mms-1)the saturated hydraulic conductivity 
	double psi_sat;          (mm)the saturated matric potential 
	double porosity;         the saturated volumetric water content(porosity) 
	double CH_b;             the Clap and Hornberger constant 
	*/

	int ok=1;
	
	/* calculating the specified parameters for each texture type by their original percentage of sand,silt, and clay. By Cosby et al.1984 */
	vspar->k_solids=(8.80*vspar->psandy+2.92*vspar->pclay)/(vspar->psandy+vspar->pclay);            /* (Wm-1k-1)the thermal conductivity of the soil solids */
	vspar->c_solids=(2.128*vspar->psandy+2.385*vspar->pclay)/(vspar->psandy+vspar->pclay)*1000000;  /* (Jm-3k-1)the heat capacity of the soil solids */
	vspar->k_sat=0.0070556*pow(10,-0.884+0.0153*vspar->psandy)*0.001;                         /* (ms-1)the saturated hydraulic conductivity */
	vspar->psi_sat=-10*pow(10,1.88-0.013*vspar->psandy)*0.001;                                /* (m)the saturated matric potential */
	vspar->porosity=0.489-0.00126*vspar->psandy;                                              /* the saturated volumetric water content(porosity) */
	vspar->CH_b=2.91+0.159*vspar->pclay;                                                      /* teh Clap and Hornberger constant */    

	return(!ok);





	/*
	switch(typeNo)
	{
	case 1 :
		vspar->psandy=92;
		vspar->psilt=5;
		vspar->pclay=3;
		vspar->k_solids=8.6143;
		vspar->c_solids=2136116;
		vspar->k_sat=0.000023558;
		vspar->psi_sat=-0.04729;
		vspar->porosity=0.3731;
		vspar->CH_b=3.39;
		break;
	case 2:
        vspar->psandy=82;
		vspar->psilt=12;
		vspar->pclay=6;
		vspar->k_solids=8.3991;
		vspar->c_solids=2145523;
		vspar->k_sat=0.000016563;
		vspar->psi_sat=-0.06394;
		vspar->porosity=0.3857;
		vspar->CH_b=3.86;
		break;
    case 3:
        vspar->psandy=58;
		vspar->psilt=32;
		vspar->pclay=10;
		vspar->k_solids=7.9353;
		vspar->c_solids=2165794;
		vspar->k_sat=0.000007111;
		vspar->psi_sat=-0.13188;
		vspar->porosity=0.4159;
		vspar->CH_b=4.50;
		break;
    case 4:
        vspar->psandy=43;
		vspar->psilt=39;
		vspar->pclay=18;
		vspar->k_solids=7.0649;
		vspar->c_solids=2203836;
		vspar->k_sat=0.000004912;
		vspar->psi_sat=-0.20734;
		vspar->porosity=0.4348;
		vspar->CH_b=5.77;
		break;
    case 5:
        vspar->psandy=17;
		vspar->psilt=70;
		vspar->pclay=13;
		vspar->k_solids=6.2520;
		vspar->c_solids=2239367;
		vspar->k_sat=0.000001677;
		vspar->psi_sat=-0.45425;
		vspar->porosity=0.4676;
		vspar->CH_b=4.98;
		break;
    case 6:
        vspar->psandy=58;
		vspar->psilt=15;
		vspar->pclay=27;
		vspar->k_solids=6.9323;
		vspar->c_solids=2209635;
		vspar->k_sat=0.000007111;
		vspar->psi_sat=-0.13188;
		vspar->porosity=0.4159;
		vspar->CH_b=7.20;
		break;
    case 7:
        vspar->psandy=32;
		vspar->psilt=34;
		vspar->pclay=34;
		vspar->k_solids=5.7709;
		vspar->c_solids=2260394;
		vspar->k_sat=0.000002845;
		vspar->psi_sat=-0.28893;
		vspar->porosity=0.4487;
		vspar->CH_b=8.32;
		break;
    case 8:
        vspar->psandy=10;
		vspar->psilt=56;
		vspar->pclay=34;
		vspar->k_solids=4.2564;
		vspar->c_solids=2326591;
		vspar->k_sat=0.000001311;
		vspar->psi_sat=-0.56104;
		vspar->porosity=0.4764;
		vspar->CH_b=8.32;
		break;
    case 9:
        vspar->psandy=52;
		vspar->psilt=6;
		vspar->pclay=42;
		vspar->k_solids=6.1728;
		vspar->c_solids=2242830;
		vspar->k_sat=0.000005756;
		vspar->psi_sat=-0.15805;
		vspar->porosity=0.4235;
		vspar->CH_b=9.59;
		break;
    case 10:
        vspar->psandy=6;
		vspar->psilt=47;
		vspar->pclay=47;
		vspar->k_solids=3.5856;
		vspar->c_solids=2355906;
		vspar->k_sat=0.000001139;
		vspar->psi_sat=-0.63299;
		vspar->porosity=0.4814;
		vspar->CH_b=10.38;
		break;
    case 11:
        vspar->psandy=22;
		vspar->psilt=20;
		vspar->pclay=58;
		vspar->k_solids=4.5370;
		vspar->c_solids=2314325;
		vspar->k_sat=0.000002;
		vspar->psi_sat=-0.39066;
		vspar->porosity=0.4613;
		vspar->CH_b=12.13;
		break;
	case 12:
        vspar->k_solids=0.2500;
		vspar->c_solids=2500000;
		vspar->k_sat=0.02;
		vspar->psi_sat=-0.1200;
		vspar->porosity=0.700;
		vspar->CH_b=4.00;
		break;
	case 13:
        vspar->k_solids=0.2500;
		vspar->c_solids=2500000;
		vspar->k_sat=0.15;
		vspar->psi_sat=-0.1200;
		vspar->porosity=0.900;
		vspar->CH_b=1.00;
		break;
	case 14:
        vspar->k_solids=0.2500;
		vspar->c_solids=2500000;
		vspar->k_sat=0.2;
		vspar->psi_sat=-0.0850;
		vspar->porosity=0.9500;
		vspar->CH_b=0.50;
		break;
	default:*/
}