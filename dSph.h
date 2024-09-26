#ifndef dsph_h
#define dsph_h

#include"main.h"
#include"Components.h"
using namespace std;

class dSph {

public:
  dSph() { } // The default constructor
  dSph(int);// The constructor
  dSph(int,int,int,double,double,double,double);// The constructor for spectral line searches

  void setup(vector<double> parvar, double& rescale); // to assign params in the loop to dSph variables
  void compute(double,double&, int&, double&,double&); // Routine that performs computation and comparison with maps
  void GetData(int ); // routine to read (and print) main input
  void GetData(int , int, int); // routine to read (and print) main input from a 3D cube
  void GetOut(); // routine to print output maps
  double get_C(double , double , int );
  double get_p(double , double , int ); 
  double C_CRE (double , double );
  double p_CRE (double , double );
//  double intF (double x, void * params);
//  double intG (double x, void * params);
  double get_integral_Fpart(double fact,double rr);
  double get_integral_Gpart(double fact,double rr);
  void load_el_spectrum(string,int);
  void get_synchrotron_emissivity (int , double ,double& , double& );
  void emissivity (string , double, double& , double& );
  void emiss_tab (string , double, double& , double& );
  double diffint(string ,double);
  double diffintlos(string ,double);
  double diffintlosell(string ,double ,double);
  double phifact(double, double, double );
  double gaweint(double, double, double );
  double getnenodiff (vector<double> ,vector<double> ,double , double );
  double getnefree (vector<double> ,vector<double> , double );
  void getnediff (vector<vector<double> >, vector<vector<double> >& );
  double get_B_equip(double );
  double bloss(double , double);
  double ICbloss(double , double , double , double );
  double brembloss(double , double , double );
  void rundscode(int );
  void CRtabsph(string );
  void CRtabell(string );
  void GetROI(); // routine to select pixels which maximize S/N
  void Annuli();
  double getB(double rr) {return magBdSph->getBfield(rr);}
  double getD(double rvar, double Evar) {return DcoeffdSph->getDcoeff(rvar,Evar);}
  double getdDdr(double rvar,double Evar) {return DcoeffdSph->getdDdr(rvar,Evar);}
  inline double pbeamGMRT610(double xxvar) {return 1.0+(-3.486e-3)*xxvar*xxvar+(47.749e-7)*pow(xxvar,4)+(-35.203e-10)*pow(xxvar,6)+(10.399e-13)*pow(xxvar,8);}
  double getAverageEmission();

protected:
  vector<double> ra, dec, map, mapmodel, maprms, mapsource, mapexp; // vectors for map-coords, obs, th, and rms emissions in each pixel.
  vector<int> pixSN; // signal to noise in each pixel
  vector<unsigned short int> mapflag;
  double dist,ra0,dec0,egb,rintin,rintfin,beamx,beamy,nuobs,deltanuobs,maxdist;
  int eload,jload,idSphrec;
  string name;
  const char *dsfile;
  vector<vector<double> > specvec,specindex,crenorm,tabCRxyell; // matrix in (r,E) for spectrum, spectral index, and normalization
  vector<double> parrec, rcrevec, Ecrevec, tabCRx, tabCRy; // vectors for model: params, radius and energy of ne, theta and emission.
  vector<double> rjvec,jIvec,jPvec; // vectors for tabulated version of emissivity: radius, intensity, polarization.
  Bfield* magBdSph;
  Dcoeff* DcoeffdSph;
};



#endif
