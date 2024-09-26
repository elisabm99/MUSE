#ifndef Components_h
#define Components_h

#include"main.h"
using namespace std; 
  
class DM {
 public: 
   DM() { } // The default constructor
   DM(int, vector<double>, double);// The normal constructor
   double darkmaint(double);
   double darkmajdec(double, double);
   double getDMgprof(double );
   void load_QDM_spectrum( );
   void load_QDM_spectrumDS( );
   void get_QDM_spectrum(double,vector<double>& , vector<double>& ,int );
 protected:
  double distDM,r0DM,rho0DM,r0emDM,j0emDM,mchi,sigmav,GammaDec;
  DMprofile DMprof;
  string dsfileDM;
  vector<double> EQvecin,Qvecin;
};

class star
{
 public: 
  star() {} // The default constructor
  star(int,vector<double>); // The normal constructor
  double crays(double );
  double getstarprof(double );
  double getnestar(double ,double );
  void load_Qstar_spectrum(double,vector<double>& , vector<double>& );

 protected:
  double r0star,pcrstar,r0emstar,j0emstar,diststar,n0star;
  starprof stprof;
};

class Bfield
{
 public: 
  Bfield() {} // The default constructor
  Bfield(int,vector<double>); // The normal constructor
  double getBfield(double);

 protected:
  double B0,r0B;
  Bprofile Bprof;
};

class Dcoeff 
{
 public: 
   Dcoeff() {} // The default constructor
   Dcoeff(int,vector<double>); // The normal constructor
   double getDcoeff(double ,double);
   double getdDdr(double , double );
   double D0,alphaD,r0D;
   Dprofile Dprof;
 protected:
};

#endif

