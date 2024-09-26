#include "main.h"
#include "Tools.h"
using namespace std;

double Tools::besselk23approx(double x)
{// This is a tabulated version of x * K_2/3(x), where K is the modified Bessel function
	static double cx[] ={0.,0.001,0.005,0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.,8.,9.,10.};
	static double cy[] ={0,0.107,0.184,0.231,0.312,0.388,0.438,0.475,0.527,0.56,0.582,0.596,0.607,0.603,0.590,0.57,0.547,0.521,0.494,0.439,0.386,0.336,0.290,0.250,0.168,0.111,0.0726,0.047,0.0298,0.0192,0.0077,0.0031,0.0012,0.00047,0.00018};
   int ic=0;
   if(x<=0.||x>=10.) return 0;
   while(x>cx[ic]) ic++;
   double x1=cx[ic];double x2=cx[ic-1];
   double y1=cy[ic];double y2=cy[ic-1];
   double coeffa=(y1-y2)/(x1-x2);
   double coeffb=y1-coeffa*x1;
   double res=coeffa*x+coeffb;
   return res;
}


#define EPS 1.0e-5
#define FPMIN 1.0e-30
#define MAXIT 1000
#define XMIN 2.0
#define PIg 3.141592653589793

void Tools::bessik(float x, float xnu, float *ri, float *rk, float *rip, float *rkp)
{
//	void beschb(double x, double *gam1, double *gam2, double *gampl,
//		double *gammi);
//	void nrerror(char error_text[]);
	int i,l,nl;
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
	if (x <= 0.0 || xnu < 0.0) {cout<<"bad arguments in bessik"<<endl;abort();}
	nl=(int)(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) {cout<<"x too large in bessik; try asymptotic expansion"<<endl;abort();}
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=PIg*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) {cout<<"bessk series failed to converge"<<endl;abort();}
		rkmu=sum;
		rk1=sum1*xi2;
	} else {
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=2;i<=MAXIT;i++) {
			a -= 2*(i-1);
			c = -a*c/i;
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) < EPS) break;
		}
		if (i > MAXIT) {cout<<"bessik: failure to converge in cf2"<<endl;abort();}
		h=a1*h;
		rkmu=sqrt(PIg/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	*ri=(rimu*ril1)/ril;
	*rip=(rimu*rip1)/ril;
	for (i=1;i<=nl;i++) {
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	*rk=rkmu;
	*rkp=xnu*xi*rkmu-rk1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN
#undef PIg
/* (C) Copr. 1986-92 Numerical Recipes Software 1>. */

#define NUSE1 5
#define NUSE2 5

void Tools::beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
//	float chebev(float a, float b, float c[], int m, float x);
	float xx;
	static float c1[] = {
		-1.142022680371172e0,6.516511267076e-3,
		3.08709017308e-4,-3.470626964e-6,6.943764e-9,
		3.6780e-11,-1.36e-13};
	static float c2[] = {
		1.843740587300906e0,-0.076852840844786e0,
		1.271927136655e-3,-4.971736704e-6,-3.3126120e-8,
		2.42310e-10,-1.70e-13,-1.0e-15};
	xx=8.0*x*x-1.0;
	*gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
}
#undef NUSE1
#undef NUSE2

float Tools::chebev(float a, float b, float c[], int m, float x)
{
//	void nrerror(char error_text[]);
	float d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) {cout<<"x not in range in routine chebev"<<endl;abort();}
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 1>. */


 void Tools::locate(vector <double> xx, double x, int &j) 
//modified for a zero-offset array
{
	int ju,jm,jl;
	int ascnd;
        int n=xx.size();
	jl=-1;
	ju=n;
	ascnd=(xx[n-1] >= xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 0;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
        if (x==xx[0]) j=0;
        else if (x==xx[n-1]) j=n-2;
	else j=jl;
}

 double Tools::linint(vector <double> xvec,vector <double> yvec, double x){
   if(x<xvec[0]||x>xvec[xvec.size()-1]) {
     cout<<" linint called out of range: "<<x<<" "<<xvec[0]<<" "<<xvec[xvec.size()-1]<<endl;
     return (x<xvec[0])? yvec[0]: yvec[xvec.size()-1];}
   int i;
   locate(xvec,x,i);
   double coeffa=(yvec[i]-yvec[i+1])/(xvec[i]-xvec[i+1]);
   double coeffb=yvec[i]-coeffa*xvec[i];
   return coeffa*x+coeffb;}     

double Tools::BiInterpolated (vector<vector<double> > vectmp, vector<double> vecx, vector<double> vecy, double xx, double yy)   {
 if(xx>vecx[vecx.size()-1]||xx<vecx[0]||yy>vecy[vecy.size()-1]||yy<vecy[0])
  {cerr<<"BiInterpolated out of range "<<xx<<" "<<vecx[0]<<" "<<vecx[vecx.size()-1]<<" "<<yy<<" "<<vecy[0]<<" "<<vecy[vecy.size()-1]<<endl;
   return 0;}
 else{
   int ix=0; int iy=0;
   locate(vecy, yy, iy); // locate gives i of your element is in between i and i+1
   locate(vecx, xx, ix);
   double t=(xx-vecx[ix])/(vecx[ix+1]-vecx[ix]);
   double u=(yy-vecy[iy])/(vecy[iy+1]-vecy[iy]);
   double res=(1.-t)*(1.-u)*vectmp[ix][iy]+t*(1.-u)*vectmp[ix+1][iy]+t*u*vectmp[ix+1][iy+1]+(1.-t)*u*vectmp[ix][iy+1];
   return res;}
 };

void Tools::tridag(const vector<double> a, const vector<double> b, const vector<double> c, const vector<double> r, vector<double>& u) {
 
   int j = 0;
   int n = a.size(); 
   double bet = 0.0; 
   static const int nmax = 1000;
   static double gam[nmax];
   //    One vector of workspace, gam, is needed. 
   if (b[0] == 0.0) cerr << "Error 1 in tridag" << endl; 
   //If this happens, then you should rewrite your equations as a set of order N-1, with u1 trivially eliminated. 
   u[0] = r[0] / (bet = b[0]); 
   for (j = 1; j < n; j++) {     //Decomposition and forward substitution.
     gam[j] = c[j-1]/bet; 
     bet = b[j] - a[j]*gam[j];
     if (bet == 0.0) cerr << "Error 2 in tridag" << endl; 
     u[j] = (r[j] - a[j]*u[j-1])/bet;
   } 
   for (j = (n-2); j >= 0; j--) u[j] -= gam[j+1]*u[j+1];    //Backsubstitution.
   
   return ;
 }


