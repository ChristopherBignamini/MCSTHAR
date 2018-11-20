#ifndef HERWIG6510WRAPPER_H
#define HERWIG6510WRAPPER_H

/**
* Herwig6510 maximum event record length (See herwig6510.inc file in Herwig6510 release sources)
*/  
const int nmxhep(4000);

/**
* Herwig6510 maximum number of defined particles (See herwig6510.inc file in Herwig6510 release sources)
*/
const int nmxres(500);

namespace MCSTHAR
{
/**
* @brief Herwig6510 HEPEVT common block wrapper structure
*
* @author Christopher Bignamini
*/
extern "C"
{
	extern struct
    {
        int NEVHEP;
        int NHEP;
		int ISTHEP[nmxhep];
        int IDHEP[nmxhep];
		int JMOHEP[nmxhep][2];
        int JDAHEP[nmxhep][2];
		double PHEP[nmxhep][5];
        double VHEP[nmxhep][4];
	} hepevt_;
}
#define hepevt hepevt_
}

/**
* @brief Herwig6510 HWEVNT common block wrapper structure
*
* @author Christopher Bignamini
*/
extern "C"
{
    extern struct
    {
		double AVWGT;
        double EVWGT;
        double GAMWT;
        double TLOUT;
        double WBIGST;
        double WGTMAX;
        double WGTSUM;
        double WSQSUM;
		int IDHW[nmxhep];
		int IERROR;
        int ISTAT;
        int LWEVT;
        int MAXER;
        int MAXPR;
		int NOWGT;
		int NRN[2];
		int NUMER;
        int NUMERU;
        int NWGTS;
		bool GENSOF; 
    } hwevnt_;
}
#define hwevnt hwevnt_

/**
* @brief Herwig6510 HWBMCH common block wrapper structure
*
* @author Christopher Bignamini
*/
extern "C"
{
    extern struct
    {
		char PART1[8];
        char PART2[8];
    } hwbmch_;
}
#define hwbmch hwbmch_

/**
* @brief Herwig6510 HWPROC common block wrapper structure
*
* @author Christopher Bignamini
*/
extern "C"
{
    extern struct
    {
		double EBEAM1;
        double EBEAM2;
        double PBEAM1;
        double PBEAM2;
        int IPROC;
        int MAXEV;
    } hwproc_;
}
#define hwproc hwproc_

/**
* @brief Herwig6510 HWPRAM common block wrapper structure
*
* @author Christopher Bignamini
*/
extern "C"
{
    extern struct
    {
		double AFCH[2][16];
        double ALPHEM;
        double B1LIM;
        double BETAF;
        double BTCLM;
        double CAFAC;
        double CFFAC;
	    double CLMAX;
        double CLPOW;
        double CLSMR[2];
        double CSPEED;
        double ENSOF;
        double ETAMIX;
        double F0MIX;
        double F1MIX;
        double F2MIX;
        double GAMH;
	    double GAMW;
        double GAMZ;
        double GAMZP;
        double GEV2NB;
        double H1MIX;
        double PDIQK;
        double PGSMX;
        double PGSPL[4];
        double PHIMIX;
        double PIFAC;
	    double PRSOF;
        double PSPLT[2];
        double PTRMS;
        double PXRMS;
        double QCDL3;
        double QCDL5;
        double QCDLAM;
        double QDIQK;
        double QFCH[16];
        double QG;
	    double QSPAC;
        double QV;
        double SCABI;
        double SWEIN;
        double TMTOP;
        double VFCH[2][16];
        double VCKM[3][3];
        double VGCUT;
        double VQCUT;
	    double VPCUT;
        double ZBINM;
        double EFFMIN;
        double OMHMIX;
        double ET2MIX;
        double PH3MIX;
        double GCUTME;
		int IOPREM;
        int IPRINT;
        int ISPAC;
        int LRSUD;
        int LWSUD;
        int MODPDF[2];
        int NBTRY;
        int NCOLO;
        int NCTRY;
	    int NDTRY;
        int NETRY;
        int NFLAV;
        int NGSPL;
        int NSTRU;
        int NSTRY;
        int NZBIN;
        int IOP4JT[2];
        int NPRFMT;
		int AZSOFT;
        int AZSPIN;
        int CLDIR[2];
        int HARDME;
        int NOSPAC;
        int PRNDEC;
        int PRVTX;
        int SOFTME;
	    int ZPRIME;
        int PRNDEF;
        int PRNTEX;
        int PRNWEB;
    } hwpram_;
}
#define hwpram hwpram_

/**
* @brief Herwig6510 HWPROP common block wrapper structure
*
* @author Christopher Bignamini
*/
extern "C"
{
	extern struct
    {
		double RLTIM[nmxres+1];
        double RMASS[nmxres+1];
        double RSPIN[nmxres+1];
		int ICHRG[nmxres+1];
        int IDPDG[nmxres+1];
        int IFLAV[nmxres+1];
		int NRES;
		bool VTOCDK[nmxres+1];
        bool VTORDK[nmxres+1];
        bool QORQQB[nmxres+1];
        bool QBORQQ[nmxres+1];
	} hwprop_;
}
#define hwprop hwprop_

/**
* @brief Herwig6510 subroutine and function wrapper declaration
*
* @author Christopher Bignamini
*/
#define hwigin hwigin_
#define hwuinc hwuinc_
#define hwusta hwusta_
#define hweini hweini_
#define hwuine hwuine_
#define hwepro hwepro_
#define hwbgen hwbgen_
#define hwdhob hwdhob_
#define hwdhvy hwdhvy_
#define hwdhad hwdhad_
#define hwmevt hwmevt_
#define hwufne hwufne_
#define hwefin hwefin_
extern "C"
{
	void hwigin(void);
	void hwuinc(void);
	void hwusta(const char[8]);
	void hweini(void);
	void hwuine(void);
	void hwepro(void);
	void hwbgen(void);
	void hwdhob(void);
    void hwdhad(void);
	void hwdhvy(void);
	void hwmevt(void);
	void hwufne(void);
	void hwefin(void);
}

#endif