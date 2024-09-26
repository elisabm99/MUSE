#ifndef constparams_h
#define constparams_h
#include <string>
#include <cmath>
using namespace std;

typedef enum {
  NFW,
  Ein,
  Iso,
  Bur,
  Gen,
  COUNTdm} DMprofile;

typedef enum {
  Plummer,
  Plummod,
  Sersic,
  King,
  Gaussian,
  Naive,
  COUNTst} starprof;

typedef enum {
  ConstB,
  ExpB,
  EquipB,
  COUNTB} Bprofile;

typedef enum {
  ConstD,
  ExpD,
  ConstConstD,
  BrokenPLD,
  COUNTD} Dprofile;

// constant parameters
const string data_path="./data/";
const string chi2listname = "chi2list_eri2.dat";
//const string mapname="LeoT-with-propagated-variance.fits", flagname="LeoT-with-propagated-variance.fits", rmsname="LeoT-with-propagated-variance.fits", maskname="mask-of-origin.fits",sextramaskname="Leo_sextractor_wpv_3s.fits";
const string mapname="DATACUBE_Eridanus2_field2.fits", flagname="", rmsname="DATACUBE_Eridanus2_field2.fits";
const string maskname="Eridanus2-Field2-SKY_MASK-m26.5.fits",sextramaskname="",expname="";
const int nummap=2,numrms=3,numflag=4,nummask=1; // position of the map in the fits file
const int userms=1; //0=with average rms; 1=rms from the map "rmsname"
const int useflag=0; //0=without flagging; 1=flagging from the map "flagname"
const int usemask=1, isEriMask=1; //0=without MUSE masking; 1=MUSE masking by using map maskname // at the moment it works only if also useflag is 1
const int usesextramask=0; //0=without sextractor masking; 1=sextractor masking by using map maskname // at the moment it works only if also useflag is 1
const int useexp=0; //0=without
const double explimit = 0.2; // minimum exposure
const int usephi=0; //0=perform integral of emissivity along los; 1=use some expression for phi=\int dl j
const double rmsave    = 2.e-4; // average rms in Jy (used only in the case userms=0)
const double ra_c = 56.08791810465522, dec_c = -43.53333651050871; // coordinates of center of dwarf. set to (0,0) if they are the same as CRVAL1,CRVAL2

const int numdSph=16;
const string namedSph[numdSph] = {"car","for","scu","boo","her","seg","ret","dra","uma","lmc","leo","eri","gru","hyd","ant","dum"};
const int field = 2; //set to 0 of there is only one field. For Eridanus2 we have 1,2,3,4,5
const int rundSph[numdSph]  = {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}; // 0=dSph not included; 1=dSph included;
const int nparams       = 3;
const double compprior[nparams][3] = {{-15.0,-10.0, 0.1}, {-1.0, 1.0, 0.1}, {-0.35, 0.505, 0.095}}; // {-2.0,2.0,0.05}
const double compprior_evi[nparams][3] = {{-23.0, -22.995, 0.1},{-1.0,1.0,0.1}, {-0.35, 0.505, 0.095}};



const double obs_beamx   = (120./60.0)/(2.*sqrt(log(4.))); // sigmax in arcmin: first number is fwhm beam in arcmin (40./60.0)
const double obs_beamy   = (120./60.0)/(2.*sqrt(log(4.))); // sigmay in arcmin: first number is fwhm beam in arcmin
// beams: car=(0.054,1.13,1.48)'; for=(0.068,1.28,1.72)'; scu=(0.068,0.82,1.71)'; boo=(0.128,1.12,1.53)'; her=(0.12,0.99,1.66); seg=(0.094,1.01,1.61);
const double obs_maxdist = 1.0 ; // consider only pixels at distances < obs_maxdist (in arcmin)
const double noncenter_ra=0;//0.1332;  // deg, shift to be applied to the coordinates in case the center of the map is not the center of the dSph (segue:0.1332 for the old maps)
const double noncenter_dec=0;//0.217; // deg, shift to be applied to the coordinates in case the center of the map is not the center of the dSph (segue:0.217 for the old maps)
const double obs_freq    = 0.888; // GHz //dummy at the moment in case cube=1
const double obs_bandwidth = 0.256; // GHz //dummy at the moment in case cube=1
const string obs_primbeamtype = "dummy"; // "GMRT610"; // primary beam model, dummy if you put a random name // "GMRT610" // "dummy"
const double obs_chunks = obs_bandwidth; // GHz
const int cube=1; // 0: continuum ; 1: spectral cube

const double distdSph[numdSph]={105.,147.,86.,42.,132.,35.,30.,76.,32.,51.,417.,366.,125.,151.,1350.,100.}; // kpc ! dSph distances from 1204.1562 and 2112.09374/2101.00253
// benchmark rho0 and r0: first six entries are for NFW, second six entries for Einasto and so on (they must respect the order of ..!
// from various papers, as 1111.1165 (BUR),  1001.4532, ...
//const double rho[6*5]={0.4,0.4,0.6,10.,1.,6.1,0.4,0.4,0.6,10.,1.,6.1,1.36,0.,0.913,0,0,0.25,1.36,0.,0.913,30.,0,0,0.4,0.4,0.6,10.,1.,6.1}; // 10^8 Msun kpc^-3
//const double r0[6*5]={0.63,1.,0.52,0.092,0.4,0.06,0.63,1.,0.52,0.092,0.4,0.06,0.22,0.,0.5,0,0,0,0.6,0.35,0.37,0.09,0,0,0.63,1.,0.52,0.092,0.4,0.06};// kpc
// from Martinez 1309.2641 and 2112.09374/2101.00253
const double rho[numdSph*5]      =  {3.0,  1.3,  1.7,  4.2,  3.6,  1.8,  0.0,  2.3,  4.2, 0.051, 0.227, 0.255, 0.304,  0.134, 0.255, 0.041,
0.37,0.22,0.28,0.53,0.49,0.54,0.7,0.36,0.49,1.0,0.0,0.0,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.57,32.9,9.1,13.7,39.5,34.4,43.3,0.0,13.7,31.4,0.301,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //10^8 Msun kpc^-3
const double r0[numdSph*5]       = {0.21, 0.47, 0.39, 0.17, 0.20, 0.16,  0.0, 0.35, 0.19, 9.8,   1.45,  1.09,  0.74,   5.85,  1.32, 5.1,
0.3,0.58,0.47,0.23,0.26,0.22,0.2,0.44,1.0,0.0,0.0,0.0,0.0,0.0,0.27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.063,0.19,0.13,0.052,0.06,0.048,0.0,0.13,0.064,4.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};// kpc//const double rho[numdSph*5]={3.0,1.3,1.7,4.2,3.6,1.8,0.0,2.3,4.2,0.051,0.231,0.259,0.293,0.102,0.245,0.041,0.37,0.22,0.28,0.53,0.49,0.54,0.7,0.36,0.49,1.0,0.0,0.0,0.0,0.0,0.0,0.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.57,32.9,9.1,13.7,39.5,34.4,43.3,0.0,13.7,31.4,0.301,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //10^8 Msun kpc^-3
//const double r0[numdSph*5]={0.21,0.47,0.39,0.17,0.20,0.16,0.0,0.35,0.19,9.8,1.41,1.07,0.76,6.40,1.34,5.1,0.3,0.58,0.47,0.23,0.26,0.22,0.2,0.44,1.0,0.0,0.0,0.0,0.0,0.0,0.27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1,0.063,0.19,0.13,0.052,0.06,0.048,0.0,0.13,0.064,4.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};// kpc
//const double rho[6*5]={11.5,3.2,5.0,4.2,13.2,27.7,1.87,0.81,0.98,0.53,2.14,3.24,0,0,0,0,0,0,57.1,17.2,27.3,39.5,65.6,99.2,0,0,0,0,0,0}; //10^8 Msun kpc^-3
//const double r0[6*5]={0.10,0.29,0.21,0.17,0.09,0.05,0.12,0.29,0.22,0.23,0.10,0.07,0,0,0,0,0,0,0.045,0.13,0.09,0.052,0.038,0.026,0,0,0,0,0,0};// kpc
// benchmark stellar profile size: first 6 entries are for Plummer, second 6 for Plummmod and so on (they must respect the order of ..!
const double rstelldSph[numdSph*6]={8.2,16.6,11.3,4.2,8.6,3.4,4.6,10.0,16.0,60.0,1.27,2.31,4.16,1.65,0.72,1.0,8.2,16.6,11.3,4.2,8.6,3.4,4.6,10.0,16.0,60.0,1.27,2.31,4.16,1.65,0.72,1.0,8.2,16.6,11.3,4.2,8.6,3.4,4.6,10.0,16.0,60.0,1.27,2.31,4.16,1.65,0.72,1.0,8.2,16.6,11.3,4.2,8.6,3.4,4.6,10.0,16.0,60.0,1.27,2.31,4.16,1.65,0.72,1.0,8.2,16.6,11.3,4.2,8.6,3.4,4.6,10.0,16.0,60.0,1.27,2.31,4.16,1.65,0.72,2.6,8.2,16.6,11.3,4.2,8.6,3.4,4.6,10.0,16.0,60.0,1.27,2.31,4.16,1.65,0.72,1.0};
// arcmin !  // rh from 1204.1562: implemented to be the same for all profiles - modify the rest of the code if they are changed!!!! for retII: sqrt(3.5*5.9) as in fig.3 of 1503.02584
const double alphaEin=0.4; // used only in case of Einasto profile // alphaEin=0.17 canonical; alphaEin=0.15 from Martinez 1309.2641 ; alphaEin=0.5 from MCMC of Bonnivard (1504.03309)
const double alphaGen=-1.0; // used only in case of Gen profile // alphaGen=-1.0 similar to NFW; alphaGen=0 similar to BUR
const string symmdSph="sph"; // sph=spherical; ell=elliptical (this case implemented only for optrunst=1,10)
const double AxialRatioProj[numdSph]={0.67,0.70,0.68,1.0,0.32,0.85,0.41,0.69,0.37,1.0,1.0,1.0,1.0,1.0,1.0,1.0}; // projected axial ratio of stellar comp.; from Table 1 in 1603.08046 (no data for boo)
const double QAxialRatioDM[numdSph]={0.6,1.1,0.8,1.0,1.0,1.0,1.1,1.4,1.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0};  //  axial ratio of dark matter halo; from Table 2 in 1603.08046 (no data for boo)
const double InclAngle[numdSph]={71.,72.,68.,90.,82.,54.,79.,59.,80.,90.,90.,90.0,90.0,90.0,90.0,90.0};  // deg // inclination angle; from Table 2 in 1603.08046 (no data for boo)
const double Uisrf[numdSph]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0};  // eV/cm^3 // ISRF density // not included in the case of dSphs
const double ngasneuHI[numdSph]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.8,0.0,0.0,0.0,0.0,0.0,0.0};  // cm^-3 // neutral atomic gas density // not included in the case of dSphs
const double ngasneuH2[numdSph]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  // cm^-3 // neutral molecular gas density // not included in the case of dSphs
const double ngasion[numdSph]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  // cm^-3// ionized gas density // not included in the case of dSphs

// parameters of choice for the various targets
/*0=Carina; 1=Fornax; 2=Sculptor; 3=BootesII; 4=Hercules; 5=Segue2; 6=ReticulumII; 7=Draco; 8=Ursa MajorII; 9=LMC; 10=generic;*/
const int optrun        = 2; // 0=single model;  1=MCMC;  2=grid-scan;
const int optrundm[numdSph] = {4,4,5,4,4,4,4,4,4,4,5,5,5,5,5,0}; // 0=DM not included, 1=DM in terms of emissivity j(j0,r0); 2=DM in terms of emissivity j(mchi,sv) from DS; 3=DM in terms of ne(mchi,sv) from DS after propag;  4=DM in terms of dN/dE from DS before propag; 5=DM decay-line (mchi,Gamma)
const int optdiff    = 2; // 0=no diffusion (radiate at inj place); 1= free-escape; 2=standard diffusion;
const int optselfconf    = 0; // 0=use the specified D model; 1= compute self-generated D;
const int optrunst[numdSph] = {0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,10}; // 0=CR not included, 1=CR in terms of emissivity j(j0,r0); 2=CR in terms of ne (power-law) after propag;  3=CR in terms of Q (power-law) before propag; 10=CR in terms of I(theta) from a phenomenological function
const double pcrdSph[numdSph] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5}; // spectral index for CR electron distribution
const DMprofile DMprofdSph[numdSph] = {Ein,  Bur,  NFW,  Ein,  Bur,  NFW,  Ein,  NFW,  NFW,  NFW, NFW, NFW,   NFW,    NFW,    NFW,   NFW};
const double rc = 0.08268, nbas = 0.8537, rt = 1.0862, deltabas = 3.3506; // rc, rt in kpc
const starprof stprofdSph[numdSph] = {Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian,Gaussian};
const Bprofile BprofdSph[numdSph] = {ExpB,ExpB,EquipB,ExpB,ExpB,ExpB,ExpB,ExpB,ExpB,ConstB,ConstB,ConstB,ConstB,ConstB,ConstB,ConstB};
const Dprofile DprofdSph[numdSph] = {ConstD,ConstD,ConstD,ConstD,ConstD,ConstD,ConstD,ConstD,ConstD,BrokenPLD,ConstD,ConstD,ConstD,ConstD,ConstD,ConstD};
const int ncomponents   = 3;
const string diffsource ="cosmic-rays";  // "dark-matter"; //  "both";
const string compname[ncomponents] = {"point-sources",diffsource,"extragalactic"};



const double B0dSph[numdSph]={0.7,1.2,1.2,0.4,0.4,0.4,1.0,1.0,1.0,4.3,1.0,1.0,1.0,1.0,1.0,1.0};//{0.9,2.0,1.6,1.2,0.4,1.4}; //{0.7,1.2,1.2,0.4,0.4,0.4}; // {1.,1.,1.,1.,1.,1.}; //muG
const double r0BdSph[numdSph]={rstelldSph[0],rstelldSph[1],rstelldSph[2],rstelldSph[3],rstelldSph[4],rstelldSph[5],rstelldSph[6],rstelldSph[7],rstelldSph[8],rstelldSph[9],rstelldSph[10],rstelldSph[11],rstelldSph[12],rstelldSph[13],rstelldSph[14],rstelldSph[15]}; // arcmin
const double D0dSph[numdSph]={3.e28,1.e30,1.e30,1.e30,3.e28,3.e28,3.e28,3.e28,3.e28,1.e28,3.e28,3.e28,3.e28,3.e28,3.e28,3.e28}; // cm^2/s
const double r0DdSph[numdSph]={rstelldSph[0],rstelldSph[1],rstelldSph[2],rstelldSph[3],rstelldSph[4],rstelldSph[5],rstelldSph[6],rstelldSph[7],rstelldSph[8],rstelldSph[9]*10.0,rstelldSph[10],rstelldSph[11],rstelldSph[12],rstelldSph[13],rstelldSph[14],rstelldSph[15]}; // arcmin
const double alphaDdSph[numdSph]={0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.51,0.3,0.3,0.3,0.3,0.3,0.3};
const string DMchannel="bb";

//other numerical constants
const double cc        = 2.99792458e10;  // cm/s
const double h_planck  = 6.6260755e-27;  // cm^2 g s^-1
const double year      = 365.25*24.*3600.;  // s
const double kpc       = 3.08568e21;  //cm
const double mp        = 0.938272029; // GeV
const double me        = 0.511e-3;  // GeV
const double Pi        = 3.14159265358979312;
const double eV_to_erg = 1.60217733e-12;
const double muG_GeV2  = 5.9157e-27;
const double GeV_to_K  = 1.1605e13;
const double GeV_to_g  = 1.7827e-24;
const double qe        = 0.30282;
const double invGeV_sec= 6.5822e-25;
const double GeV_erg   = 1.6022e-3;
const double aminrad   = Pi/180./60.;
const double MsunGeV   = 1.116e57; // GeV
const double Hmass     = 0.938; // GeV, hydrogen mass
const double H2mass    = 2.0*Hmass; // GeV, hydrogen mass
const double Hionizen  = 13.6e-9;  // GeV, ionization energy for H
const double H2ionizen = 4.52e-9;  // GeV, ionization energy for H2
const double radLH     = 62.8*5.60947e23; // GeV/cm^2, radiation length for H
const double radLH2    = 61.3*5.60947e23; // GeV/cm^2, radiation length for H2



#endif

