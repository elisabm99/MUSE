#ifndef tools_h
#define tools_h

using namespace std;

class Tools {

public:
  Tools () {} // constructor
  double besselk23approx(double x);
  void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
  void bessik(float x, float xnu, float *ri, float *rk, float *rip, float *rkp);
  float chebev(float a, float b, float c[], int m, float x);
  void locate(vector<double> xx, double x, int &j);
  double linint(vector<double> xx, vector<double> yy, double x);
  double BiInterpolated (vector<vector<double> >, vector<double> , vector<double> , double , double );
  void tridag(const vector<double> , const vector<double> , const vector<double> , const vector<double> , vector<double>& );
};

#endif

