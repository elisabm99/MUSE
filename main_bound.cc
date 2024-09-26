
#include"main.h"
#include "dSph.h"
#include "Tools.h"

using namespace std;

#define STRING(num) #num


int getinpar(int ip) {return round((compprior[ip][1]-compprior[ip][0])/compprior[ip][2])+1;}

void contour2D(double chi2null) {
  cout<<"Computing 2D contour for variables: chi2null = "<<chi2null<<endl;
  vector<double> parvar,vecpar0,vecpar1,vecchi2;
  /* list of chi2 for the two main parameters, profiling-out the other params (ie taking the best-fit) */
  ofstream chi2file2D("chi2list2D.dat");
  ifstream chi2file (chi2listname);
  for(int ip=0; ip<nparams; ip++) parvar.push_back(0.); // parvar has nparams=2 elements for us
  if(nparams<=1) {cout<<"Nothing to do for 2D contour"<<endl;}
  else{ // for us nparams=2
    int nmarg=nparams-1;
    for(int ip=2; ip<nparams; ip++) nmarg*=getinpar(ip);
    for(int in1=0; in1<getinpar(0); in1++) { // loops over values of Gamma
      for(int in2=0; in2<getinpar(1); in2++){ // loops ove values of Sflat
        double chi2var=1.e30,chi2tmp; // assign a high initial value to chi2var
        for(int in=0; in<nmarg; in++){ // for us nmarg=1, so in=0
          for(int ip=0; ip<nparams; ip++) chi2file>>parvar[ip]; // parvar[0] = Gamma, parvar[1] = ~Sflat
          chi2file>>chi2tmp;
          if(chi2tmp<chi2var) chi2var=chi2tmp;} // update chi2var if a lower value is found
        chi2file2D<<log10(parvar[1])<<" "<<log10(parvar[0])<<" "<<chi2var<<endl;
        vecpar1.push_back(log10(parvar[1])); vecpar0.push_back(log10(parvar[0])); vecchi2.push_back(chi2var); //vecpar1=log10(~Sflat),vecpar0=log10(Gamma)
     // cout<<"chi2list "<<parvar[1]<<" "<<parvar[0]<<" "<<chi2var<<endl;
    }}
    chi2file2D.close();  chi2file.close();}
    //{vecpar0,vecpar1,vecchi2} are {{Gamma0, Sflat0, chi2_00},{Gamma0, Sflat1, chi2_01}... {Gamma1, Sflat0, chi2_10},{Gamma1, Sflat1, chi2_11}....}

//  /* regions of discovery/bounds for two main parameters */
//  ofstream chi2filelim("chi2bound.dat"); // this will have 2 columns: Sflat, Gamma. If a value of Sflat is repeated, it means there are both lower and upper bounds
//  ofstream chi2filedisc("chi2disc.dat");
//  double sigmadisc=3.0, Pcl=0.95, disccl=chi2null-sigmadisc*sigmadisc, discpar=-1.e30, p0, res, p0best=0, p1best=0;
//  Tools* T;
//  parvar.clear(); for(int in2=0; in2<getinpar(0); in2++) parvar.push_back(0.); //initialize array parvar of length the number of values of Gamma
//  vector<double> chi2tmp; for(int in2=0; in2<getinpar(0); in2++) chi2tmp.push_back(0.);//initialize array chi2tmp of length the number of values of Gamma
//  cout<<"getinpar(1) = "<<getinpar(1)<<endl;
//  for(int in1=0; in1<getinpar(1); in1++) { // loops over values of Sflat
//    double chi2best=1.e30,chi2worst=0.0;
//    // For given Sflat, finds best and worst chi2 among all Gammas. For best chi2, saves Sflat and Gamma (bzw. p0best, p2best)
//    for(int in0=0; in0<getinpar(0); in0++){ // loops over values of Gamma
//      p0=vecpar1[in1]; // vecpar1 contains log10(~Sflat). Notice the meaning of 0,1 is changed!! now p0->Sflat
//      parvar[in0]=vecpar0[in0*getinpar(1)]; // picks up all the different values of Gamma. vecpar0 contains log10(Gamma)
//      chi2tmp[in0]=vecchi2[in1+in0*getinpar(1)]; // finds the chi2 that corresponds to Gamma, Sflat above. See line 33
//      if(chi2tmp[in0]>chi2worst) chi2worst=chi2tmp[in0]; // finds largest chi2 for given Sflat
//      if(chi2tmp[in0]<chi2best) {chi2best=chi2tmp[in0]; p0best=p0; p1best=parvar[in0];} // finds smallest chi2 for given Sflat. p0best->Sflat, p1best->Gamma
////      cout<<"chi2list2 "<<p0<<" "<<parvar[in0]<<" "<<chi2tmp[in0]<<" "<<chi2best<<" "<<chi2worst-chi2best<<endl;
//    }
//
//    // Find discovery
//    if(chi2null-chi2best>discpar) discpar=chi2null-chi2best; //sqrt(discpar) = # of sigmas of discrepancy b/w best fit and null hypothesis
//    double oneminusP0=1.0;//(chi2best>chi2null) ? 1.0 : erfc(sqrt((chi2null-chi2best)/2.0));
////    cout<<" chi2best-chi2null "<<chi2best-chi2null<<" "<<chi2best<<" "<<oneminusP0<<endl;
//    if(chi2best>disccl) chi2filedisc<<pow(10.0,p0)<<" "<<1.e100<<endl; // no discovery at >= 3sigma. disccl=chi2null-sigmadisc*sigmadisc
//    else if(chi2worst<disccl) chi2filedisc<<pow(10.0,p0)<<" "<<1.e-100<<endl; //
//    else {
//      for(int in0=1; in0<getinpar(0); in0++) {
//        if((chi2tmp[in0-1]<disccl && chi2tmp[in0]>disccl) || (chi2tmp[in0-1]>disccl && chi2tmp[in0]<disccl)) {
//          res=(parvar[in0]-parvar[in0-1])/(chi2tmp[in0]-chi2tmp[in0-1])*(disccl-chi2tmp[in0])+parvar[in0];// parvar contains all values of log10(Gamma).
//          chi2filedisc<<pow(10.0,p0)<<" "<<pow(10.0,res)<<endl; }}
//     }
//
//    // Find upper and lower limits
//    if(chi2null<chi2best && (1.0-0.5*erfc(sqrt(chi2best-chi2null)/2.0))/oneminusP0>Pcl) chi2filelim<<pow(10.0,p0)<<" "<<1.e-100<<endl;
//    else{
//      chi2best=(chi2best<chi2null) ? chi2best : chi2null;
//      if((1.0-0.5*erfc(sqrt(chi2worst-chi2best)/2.0))/oneminusP0<Pcl) chi2filelim<<pow(10.0,p0)<<" "<<1.e100<<endl;
//      else {
//        for(int in0=1; in0<getinpar(0); in0++){ // loop over Gammas. chi2tmp contains all the chi2 for different Gammas, for given Sflat
//          double P1=1.0-0.5*erfc(sqrt((chi2tmp[in0]-chi2best)/2.0))/oneminusP0, P0=1.0-0.5*erfc(sqrt((chi2tmp[in0-1]-chi2best)/2.0))/oneminusP0;
////            cout<<"bounds "<<pow(10.0,p0)<<" "<<chi2tmp[in0]-chi2best<<" "<<P0<<" "<<P1<<endl; // p0 is log10(~Sflat)
//          if((P1<Pcl&&P0>Pcl)||(P0<Pcl&&P1>Pcl)){ // if only one of P0, P1 is greater than the confidence level. This works because the Gammas are ordered in ascending order, so the chi2s should also be in order, acending or descending
//            res=(parvar[in0]-parvar[in0-1])/(P1-P0)*(Pcl-P1)+parvar[in0]; // parvar contains all values of log10(Gamma). Here we interpolate to find which value of Gamma corresponds to the confidence level
////            cout<<"Bound. i_Sflat"<<in0<<"  Sflat = "<<pow(10, p0)<<"   P0 = "<<P0<<"   P1 ="<<P1<<"  chi2best = "<<chi2best<<"  chi2high = "<<chi2tmp[in0]<<"  chi2low = "<<chi2tmp[in0-1]<<" Gammabound ="<<pow(10.0,res)<<endl;
////            cout<<"log_Gamma_high ="<<parvar[in0]<<"  log_Gamma_low ="<<parvar[in0-1]<<endl;
////            cout<<pow(10.0,p0)<<" "<<pow(10.0,res)<<endl;
//            chi2filelim<<pow(10.0,p0)<<" "<<pow(10.0,res)<<endl;}}}}
//  }
//  cout<<" Discovery at sigma = "<<sqrt(discpar)<<" "<<pow(10.,p0best)<<" "<<pow(10.,p1best)<<endl;
//  chi2filedisc.close(); chi2filelim.close();
  return;}



int main() {

  time_t tin,tfin;
  time(&tin); // set initial time
  double time_compute1, time_compute2, time_compute1_tot, time_compute2_tot;
  time_compute1_tot = 0;   time_compute2_tot = 0;

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
//    cout<<"\nABORTING"<<endl; abort();
    cout << "\nComputing model predictions and comparing to map" << endl;
    vector<int> inpar,count;
    double chi2best=1.e30;
    int ntot=1;
    ofstream outfile(chi2listname); // this is the file of the chi2_lists
    ofstream outfilebest("chi2best.dat"); // this is the file of the chi2_nulls
    for(int ip=0; ip<nparams; ip++) {
      cout<<"ip = "<<ip<<" getinpar(ip) = "<<getinpar(ip)<<endl;
      inpar.push_back(getinpar(ip));
      cout<<"compprior. ip = "<<ip<<" inpar(ip) = "<<inpar[ip]<<" compprior = "<<compprior[ip][0]<<" "<<compprior[ip][1]<<" "<<compprior[ip][2]<<endl;
      ntot*=inpar[ip]; // total number of elements in the parameter grid
      paramsvar.push_back(0.);
      count.push_back(0);}
    for(int in=0; in<ntot; in++){
//      cout<<"in = "<<in<<endl;
      int part=in,ichi2=0;
      double chi2=0.,chi2tmp=0.;
      for(int ip=nparams-1; ip>=0; ip--){
        count[ip]=(((part+1) % inpar[ip])==0) ? inpar[ip]-1 : ((part+1) % inpar[ip])-1;
        paramsvar[ip]=pow(10.,compprior[ip][0]+compprior[ip][2]*count[ip]); // contains current values of the parameters
//        cout<<"count = "<<count[ip]<<endl;
//        cout<<"ip = "<<ip<<", count(ip) = "<<count[ip]<<", paramsvar(ip) = "<<paramsvar[ip]<<endl;
        part=(part-count[ip])/inpar[ip];}
//      cout<<" Count "<<in<<" of "<<ntot-1<<" "<<count[0]<<" "<<count[1]<<" "<<endl;
//        cout<<"paramsvar = "<<paramsvar[0]<<" "<<paramsvar[1]<<" "<<paramsvar[2]<<" "<<endl;
      for (int id=0; id<numdSph; id++) {
        if(rundSph[id]==1) {
          double rescale=0;
          targets[id].setup(paramsvar,rescale); // redefine dSph-params according to paramsvar
          if (rescale != 1) cout<<"After setup. rescale = "<<rescale<<endl;
          targets[id].compute(rescale,chi2tmp, ichi2, time_compute1, time_compute2); // compute map and chi-2 (if rescale.ne.1 just rescale the map)
          time_compute1_tot += time_compute1;
          time_compute2_tot += time_compute2;
          chi2+=chi2tmp;}}
//      std::cout.rdbuf(old_buffer);
//        cout<<"number of good pixels = "<<ichi2<<endl;
//      cout<<" Params "<<in<<" "<<paramsvar[0]<<" "<<paramsvar[1]<<" "<<paramsvar[2]<<" "<<chi2<<endl;
      outfile <<paramsvar[0]<<" "<<paramsvar[1]<<" "<<paramsvar[2]<<" "<<chi2<<endl;
      chisquare.push_back(chi2);

        }
            outfile.close();
  }


  time(&tfin); // set final time
  cout << "Time length for this run " << (double)(tfin-tin) << " s." << endl;
  cout << "Time for integration " << time_compute1_tot << " s." << endl;
  cout << "Time for calculating chi2 " << time_compute2_tot << " s." << endl;

  cout << "Done.\n" << endl;

  return 0;
}



