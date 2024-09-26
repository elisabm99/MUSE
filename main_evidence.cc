
#include"main.h"
#include "dSph.h"
#include "Tools.h"

using namespace std;

#define STRING(num) #num


int getinpar(int ip) {return round((compprior_evi[ip][1]-compprior_evi[ip][0])/compprior_evi[ip][2])+1;}

int main() {

  time_t tin,tfin;
  time(&tin); // set initial time
  double time_compute1, time_compute2, time_compute1_tot, time_compute2_tot;
  time_compute1_tot = 0;   time_compute2_tot = 0;

  int idSph;
  for (int id=0; id<numdSph; id++) {
    if(rundSph[id]==1) {
      idSph = id;
      cout<<namedSph[idSph]<<endl;
    }
   }

  cout << "Starting the code" <<endl;
  vector<double> bestparams,paramsvar,chisquare; // best-fit params, auxiliary params, and series of chi2.
  for (int in=0; in<nparams; in++) bestparams.push_back(0.); // create the params-vector

  if(optrun==2){
    vector<dSph> targets; // vector of dSphs
    for (int in=0; in<numdSph; in++) {
      if(cube==1) {
        int ichlow,ichsup; double flat,deltanuvar,nuvar,angres;
        ifstream infilecube("inputfile_cube.dat"); infilecube>>ichlow>>ichsup>>flat>>nuvar>>deltanuvar>>angres;  infilecube.close();
//        cout<<"nuvar = "<<nuvar<<endl;
        dSph obj(in,ichlow,ichsup,flat,nuvar,deltanuvar,angres); targets.push_back(obj);}
    } // load obs data and initialize params
    cout << "Computing model predictions and comparing to map" << endl;
    vector<int> inpar,count;
    double chi2best=1.e30;
    int ntot=1;
    ofstream outfilebest("chi2best.dat");
    for(int ip=0; ip<nparams; ip++) {
      cout<<"getinpar(ip) "<<getinpar(ip)<<endl;
      inpar.push_back(getinpar(ip));
      cout<<"compprior_evi "<<ip<<" "<<inpar[ip]<<" "<<compprior_evi[ip][0]<<" "<<compprior_evi[ip][1]<<" "<<compprior_evi[ip][2]<<endl;
      ntot*=inpar[ip];
      paramsvar.push_back(0.);
      count.push_back(0);
    }
//    contour2D(0.0); abort(); // just for some checks!!
    for(int in=0; in<ntot; in++){
      int part=in,ichi2=0;
      double chi2=0.,chi2tmp=0.;
      for(int ip=nparams-1; ip>=0; ip--){
        count[ip]=(((part+1) % inpar[ip])==0) ? inpar[ip]-1 : ((part+1) % inpar[ip])-1;
        paramsvar[ip]=pow(10.,compprior_evi[ip][0]+compprior_evi[ip][2]*count[ip]);
        part=(part-count[ip])/inpar[ip];}
//      cout<<" Count "<<in<<" of "<<ntot-1<<" "<<count[0]<<" "<<count[1]<<" "<<endl;
      for (int id=0; id<numdSph; id++) {
        if(rundSph[id]==1) {
          double rescale=0;
          targets[id].setup(paramsvar,rescale); // redefine dSph-params according to paramsvar
          if (rescale != 1) cout<<"After setup. rescale = "<<rescale<<endl;
//          cout<<"Computing"<<endl;
          targets[id].compute(rescale,chi2tmp,ichi2, time_compute1, time_compute2); // compute map and chi-2 (if rescale.ne.1 just rescale the map)
          time_compute1_tot += time_compute1;
          time_compute2_tot += time_compute2;
          chi2+=chi2tmp;
          }
      }
      if(namedSph[idSph] == "eri"){
      	      outfilebest <<paramsvar[0]<<" "<<paramsvar[1]<<" "<<paramsvar[2]<<" "<<chi2<<endl;
      }else{
          if(chi2<chi2best) {
//            cout<<" Best-Params "<<chi2<<" "<<ichi2<<" "<<chi2/(ichi2*1.0)<<" "<<paramsvar[0]<<" "<<paramsvar[1]<<" "<<paramsvar[2]<<endl;
            bestparams=paramsvar;
            chi2best=chi2;
          }
          paramsvar=bestparams; // to use best-fit in GetOut
      }
    }
    if(namedSph[idSph] != "eri"){
          cout<<"writing outfilebest. chi2best = "<<chi2best<<endl;
          outfilebest<<bestparams[0]<<" "<<bestparams[1]<<" "<<bestparams[2]<<" "<<chi2best<<endl; // writes to file logGamma_best, Sflat_best, r0_best chi2null
    }
    outfilebest.close();


  }


  time(&tfin); // set final time
  cout << "Time length for this run " << (double)(tfin-tin) << " s." << endl;
  cout << "Time for integration " << time_compute1_tot << " s." << endl;
  cout << "Time for calculating chi2 " << time_compute2_tot << " s." << endl;

  cout << "Done.\n" << endl;

  return 0;
}



