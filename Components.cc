#include "main.h"
#include "Components.h"
#include "Tools.h"

using namespace std;
/************************Class DM*********************************************/

DM::DM(int idSph, vector<double> parvar, double nuobs){ // The usual constructor
  if(optrundm[idSph]==1) {
    j0emDM=parvar[0]/(obs_beamx*obs_beamy*pow(Pi/180./60.,2)*Pi); // //divided by beam = Pi*sigmax*sigmay
    r0emDM=parvar[1];} // kpc
  else if(optrundm[idSph]==5) {
    double ma, g2Rho0;
    g2Rho0 = parvar[0]*parvar[0]; // g^2 \rho_0 in GeV^-2 10^8 M_\odot kpc^-3
    r0DM = parvar[2]; // kpc
    ma = 2 * nuobs * 2 * 3.14159 * 6.5822e-7;  // in eV
    GammaDec = ma*ma*ma / (64 * 3.14159) * g2Rho0 * 5.7664e-3; //decay rate * rho0 in s^-1 GeV cm^-3
//    cout<<"r0DM = "<<r0DM<<endl;
//    r0emDM=parvar[2];  // kpc // not used at the moment: it can be used after including a D-factor in constant and derivation of rho0
    }
  else {
    sigmav=parvar[0]; // cm^3/s 
    mchi=parvar[1];}  // GeV

  distDM=distdSph[idSph]; // kpc
  DMprof=DMprofdSph[idSph];
//  cout<<"DMprof = "<<DMprof<<" idSph = "<<idSph<<endl;
//  rho0DM=rho[idSph+DMprof*numdSph]*3.7988; // parvar[2] // GeV cm^-3
//  r0DM=r0[idSph+DMprof*numdSph]; // parvar[3] // kpc
  dsfileDM="outmodel_"; dsfileDM.append(namedSph[idSph]); dsfileDM.append(".dat");
  
}

double DM::getDMgprof(double xvar) {
/* dimensionless DM profile */
 double res, f;
 xvar=(xvar<1.e-5) ? 1.e-5 : xvar;
 switch (DMprof) {
  case 0:
    res=1./(xvar*pow(1+xvar,2)); // NFW
    break;
  case 1: 
    res=exp(-2./alphaEin*(pow(xvar,alphaEin)-1)); // Einasto
    break;
  case 2:
    res=1./(1+xvar*xvar); // isothermal
    break;
  case 3:
    res=1./((1+xvar)*(1+xvar*xvar)); // Burkert
    break;
  case 4:
    res=1./(pow(xvar,alphaGen)*pow(1+xvar*xvar,(alphaGen+3.0)/2.0)); // generalized profile
    break;
  case 5:
    if (xvar <= rt/r0DM){
        f = tanh(xvar * r0DM/rc);
        res = pow(f, nbas) * 1./(xvar*pow(1+xvar,2)) + nbas * pow(f, nbas -1) * (1 - f*f) / (xvar*xvar) *  r0DM / rc * (log(1+xvar) - xvar/(1 + xvar)); // NFWCoreTide}
    }else{
        f = tanh(rt/rc);
        res = (pow(f, nbas) * 1./(rt/r0DM * pow(1+rt/r0DM,2)) + nbas * pow(f, nbas -1) * (1 - f*f) / (rt/r0DM*rt/r0DM) *  r0DM / rc * (log(1+rt/r0DM) - rt/r0DM/(1 + rt/r0DM)) ) * pow(xvar * r0DM/rt, -deltabas);
    }
    break;
  default:
    cerr << "Problem in choosing halo model: "<< DMprof<<endl;
 }
 return res; 
};


double DM::darkmaint(double rvar)  {
/* DM-induced emissivity j in units of Jy/kpc (having defined I=int dl j)*/
//  cout<<"darkmaint"<<rvar<<" "<<r0emDM<<" "<<j0emDM<<" "<<pow(getDMgprof(rvar/r0emDM)/getDMgprof(1.),2)<<endl;
//  if(rvar>2.0) return 0.0;
  return j0emDM*pow(getDMgprof(rvar/r0emDM),2); ///getDMgprof(1.)
}; 

double DM::darkmajdec(double rvar,double deltafreqobs)  {
/* DM-induced emissivity j in units of Jy/kpc (having defined I=int dl j) in the case of spontaneous decay of DM (rvar in kpc, deltafreqobs in GHz)
   computed in the approximation of obs_bandwidth >> energy resolution >> DM velocity dispersion
   if you want to include stimulated emission, you need to add the proper factor */
//  cout<<"jdec= "<<GammaDec<<" "<<deltafreqobs<<" "<<rho0DM<<" "<<getDMgprof(rvar/r0DM)<<endl;
  return GammaDec/(deltafreqobs*1.e9) * getDMgprof(rvar/r0DM) * GeV_erg * kpc * 1.e23; // Jy/kpc
}; 



void DM::get_QDM_spectrum(double x,vector<double>& yQvec, vector<double>& Qvec,int imode)   {
/*total (e+ + e-) injection source: x in kpc or GeV, xQvec in log(GeV) or log(kpc), Qvec in log(GeV^-1 cm^-3 s^-1) */
  Qvec.clear();
  if(imode==0){
    yQvec.clear();
    for(int ie=0; ie<EQvecin.size(); ie++) {
      yQvec.push_back(EQvecin[ie]);
      Qvec.push_back(Qvecin[ie]+log(pow(rho0DM*getDMgprof(x/r0DM)/mchi,2)*sigmav/2.)); //annihilating
//      Qvec.push_back(Qvecin[ie]+log(rho0DM*getDMgprof(x/r0DM)/mchi*sigmav)); //decaying (sigmav is actually gamma) ; in DS mass should be taken to be mchi/2

// PBH case //
//      yQvec.push_back(log(0.1*pow(10,ie*4.0/(1.0*EQvecin.size()))));
//      double dNdEdot=3.82e-2*pow(mchi/1.e13*exp(yQvec[ie]),2)/(2.0*Pi*invGeV_sec*(exp(exp(yQvec[ie])/(1.06e13/mchi))+1)); Qvec.push_back(log(rho0DM*getDMgprof(x/r0DM)*GeV_to_g/mchi*sigmav*dNdEdot)); //PBH: sigmav is fPBH and mchi is mPBH in g
//    cout<<"QDM call "<<exp(yQvec[ie])<<" "<<exp(Qvec[ie])<<" "<<dNdEdot<<endl;
  } } 
  else if(imode==1){
    Tools* T;
    double res=T->linint(EQvecin,Qvecin,log(x));
    for(int ir=0; ir<yQvec.size(); ir++) Qvec.push_back(res+log(pow(rho0DM*getDMgprof(exp(yQvec[ir])/r0DM)/mchi,2)*sigmav/2.));}
  return;
};


/*****************************************************************************/

/************************Class star*******************************************/

star::star(int idSph, vector<double> parvar){ // The usual constructor
  if(optrunst[idSph]==1) {j0emstar=parvar[0]/(obs_beamx*obs_beamy*pow(Pi/180./60.,2)*Pi);} //divided by beam = Pi*sigmax*sigmay
  else{ n0star=parvar[0];}
  r0emstar=parvar[1];
  pcrstar=pcrdSph[idSph]; // spectral index for CR electron distribution
  diststar=distdSph[idSph]; // kpc
  stprof=stprofdSph[idSph];
  r0star=diststar*sin(rstelldSph[idSph]*Pi/180./60.); // parvar[3]
//  r0star=r0emstar;
}

double star::getstarprof(double xvar) {
 double res,msersic,xtking;
 switch (stprof) {
  case 0:
    res=1./pow((1.+xvar*xvar),2.5); // Plummer
    break;
  case 1:
    res=1./pow((1.+xvar*xvar),3); // modified Plummer
    break;
  case 2: 
    res=exp(-pow(xvar,1./msersic)); // Sersic ! implementare valori per msersic!!!
    break;
  case 3:
    res=pow(1./sqrt(1.+xvar*xvar)-1./sqrt(1.+xtking*xtking),2); // King ! implementare valori per xtking=rt/rc !!!!!
    break;
  case 4:
    res=exp(-xvar*xvar/2.); // Gaussian
    break;
  case 5:
    res=1./(1.+exp(xvar)); // naive case
    break;
  default:
    cout << "Problem in choosing stellar profile model: "<< stprof;
    abort();
 }
 return res; 
};


double star::crays(double rvar) {
/* CR-induced emissivity j in units of Jy/kpc (having defined I=int dl j)*/
double res=j0emstar*getstarprof(rvar/r0emstar);///getstarprof(1.); 
 return res; 
};

double star::getnestar(double rvar, double energy) {
/* CR electron density nel=n0*g(r/r0)*(E/GeV)^-p, where n0 is in cm^-3 GeV^-1 (if ne) or cm^-3 GeV^-1 s^-1 (if Q)*/
double res=n0star*getstarprof(rvar/r0star)*pow(energy,-pcrstar); 
 return res; 
};

void star::load_Qstar_spectrum(double r,vector<double>& EQvec, vector<double>& Qvec)   {
/*r in kpc, EQvec in log(GeV), Qvec in log(GeV^-1 cm^-3 s^-1) */
  int ne=50;
  Qvec.clear();
  EQvec.clear();
  for(int ie=0; ie<ne+1; ie++){
    EQvec.push_back(log(0.1*pow(10.,ie/10.))); // log(GeV), from 100 MeV to 10 TeV
    Qvec.push_back(log(max(getnestar(r,exp(EQvec[ie])),1.e-200)));  // log(GeV^-1 cm^-3 s^-1)
//    cout<<"load_Qstar_spectrum "<<exp(EQvec[ie])<<" "<<exp(Qvec[ie])<<endl;
}   
};

/*****************************************************************************/

