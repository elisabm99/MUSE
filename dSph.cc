#include "main.h"
#include "dSph.h"
#include "Tools.h"
#include "Components.h"

using namespace std;


/***********************   CLASS dSph ***********************************/


dSph::dSph(int idSph, int ichlow, int ichsup, double flat, double nuvar, double deltanuvar, double angres) {
/* routine to initialize a dSph */ //flat in the same units as the map, nuobs in GHz, deltanuvar in A, angres in arcmin
    if(rundSph[idSph]!=1) return;
    cout <<namedSph[idSph] <<" included in the analysis with spectral channels from "<<ichlow<<" to "<<ichsup<< endl;
// Do not include here params which may be looped on (because of rescale in setup)
    name=namedSph[idSph];
    string f = to_string(field);
    if (field !=0) name.append(f);
    dist=distdSph[idSph]; // kpc
    egb=flat;
    nuobs=nuvar;
    deltanuobs=deltanuvar;
//    double mchi=2.0*nuobs*2.0*Pi*6.5822e-16; // GeV, if decay into two photons
    beamx=angres; 
    beamy=angres; // elliptic case to be implemented!!!!!
    maxdist= obs_maxdist;
    GetData(idSph,ichlow,ichsup); //abort();
    idSphrec=idSph;
    rintin=1.e-3;
    rintfin=dist*sin(maxdist*aminrad)*1.2;//0.2; // setting 200 pc for LeoT // max(5.*r0[idSphrec+numdSph*DMprofdSph[idSphrec]],dist*sin(30.0*60.*aminrad));  // setting rfin=3*r0DM
    for(int in=0; in<nparams; in++) parrec.push_back(1.e30);
}


void dSph::GetData(int idSph, int ichlow, int ichsup) {
/* routine to load one dSph image from a 3D cube*/
   int nmapsinput=1+userms+useflag+usemask; // number of maps/cubes from minimum 1=only data to maximum 4=data+err+flags+maskcube
   cout<<"Number of maps = "<<nmapsinput<<endl;
   for(int ichannel=ichlow; ichannel<=ichsup; ichannel++){ //this numbering is zero based

     int binsize = ichsup - ichlow + 1;
     cout << "\nObservational Data. ichannel = "<<ichannel <<endl;
     for (int imap=0; imap<nmapsinput;imap++){
       cout <<endl<< "Loading map number  imap = "<<imap <<endl;
       int status=0,CRPIX1,CRPIX2,NAXIS1,NAXIS2,anynul,hdutype=0,hdunum;
       long long felement=1LL;
       char comment[100];
       float CRVAL1,CRVAL2,CDELT1,CDELT2,nulval=0;
       fitsfile *fptr; /* FITS file pointer */
       string filevar = data_path; filevar.append("images/");

       int maptype=0; // 0=flux; 1=rms; 2=flags; 3=mask;
       if(imap==1){ if(userms==1) maptype=1; else if(useflag==1) maptype=2; else maptype=3;}
       if(imap==2){ if(useflag==1) maptype=2; else maptype=3;}
       if(imap==3) maptype=3;

       if(maptype==0) filevar.append(mapname);
       if(maptype==1)  filevar.append(rmsname);
       if(maptype==2)  filevar.append(flagname);
       if(maptype==3) filevar.append(maskname);
       cout<<"maptype = "<<maptype<<" filevar = "<<filevar<<endl;

       const char *filename; filename=filevar.c_str();
       if( fits_open_image(&fptr,filename,READONLY,&status) ) cout<<"read observational map: open status= "<<status<<" "<<filevar<<endl;
//       fits_get_num_hdus(fptr,&hdunum,&status); // get hdu of last table

       if(maptype==1) {fits_movabs_hdu(fptr,numrms,&hdutype,&status); cout<<"rms. Moving to header number "<<numrms<<" "<<filename<<endl;}
       if(maptype==2) {fits_movabs_hdu(fptr,numflag,&hdutype,&status); cout<<"flag. Moving to header number "<<numflag<<" "<<filename<<endl;}// move to last table
       if(maptype==3) {fits_movabs_hdu(fptr,nummask,&hdutype,&status); cout<<"mask. Moving to header number "<<nummask<<" "<<filename<<endl;}// move to last table

//       cout<<imap<<"hdulast "<<hdunum<<" "<<hdutype<<endl; //fits_get_hdu_num(fptr,&hdunum); cout<<imap<<"hducurr "<<hdunum<<" "<<hdutype<<endl;
       if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) cout<<"NAXIS1 = "<<NAXIS1<<" status= "<<status<<endl;
       if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) cout<<"NAXIS2 = "<<NAXIS2<<" status= "<<status<<endl;
       if(maptype<3) if( fits_read_key(fptr,TINT,"CRPIX1",&CRPIX1,comment,&status) ) cout<<"CRPIX1 = "<<CRPIX1<<" status= "<<status<<endl;
       if(maptype<3) if( fits_read_key(fptr,TINT,"CRPIX2",&CRPIX2,comment,&status) ) cout<<"CRPIX2 = "<<CRPIX2<<" status= "<<status<<endl;
       if(maptype<3) if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) cout<<"CRVAL1 = "<<CRVAL1<<" status= "<<status<<endl;
       if(maptype<3) if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) cout<<"CRVAL2 = "<<CRVAL2<<" status= "<<status<<endl;


/*       if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) cout<<"CDELT1 = "<<CDELT1<<" status= "<<status<<endl;*/ CDELT1=0.2/3600.0;
/*       if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) cout<<"CDELT2 = "<<CDELT2<<" status= "<<status<<endl; */ CDELT2=0.2/3600.0;
       cout<<"CDELT1 = "<<CDELT1<<" CDELT2 = "<<CDELT2<<" CRPIX1 = "<<CRPIX1<<" CRPIX2 = "<<CRPIX2<<" NAXIS1 = "<<NAXIS1<<" NAXIS2 = "<<NAXIS2<<" status = "<<status<<endl;
       ra0 = ra_c - CRVAL1; // central position RA in degree
       dec0 = -dec_c + CRVAL2; // central position DEC in degree
       status=0;
       int nelements=NAXIS1*NAXIS2;
       if(imap>0&&nelements!=map.size()) {cout<<"Dimensions not matching in input tables "<<imap<<" "<<nelements<<" "<<map.size()<<endl; abort();}
       float *map_in=new float[nelements];
       unsigned short int *map_in2=new unsigned short int[nelements];
       long int *map_insou=new long int[nelements];
       felement = static_cast<long long> (nelements) * static_cast<long long> (ichannel) + 1;
//       cout<<"felement = "<<felement<<"  "<<(nelements*ichannel+1)<<" "<<nelements<<" "<<ichannel<<endl;
       if(maptype==2) if( fits_read_img(fptr,TUSHORT,felement,nelements,&nulval,map_in2,&anynul,&status) ) cout<<"read table status = "<<status<<endl;
       if(maptype<2||(maptype==3 && isEriMask==0)) if( fits_read_img(fptr,TFLOAT,felement,nelements,&nulval,map_in,&anynul,&status) ) cout<<"read table status= "<<status<<endl;
       if(maptype==3 && isEriMask==1) if( fits_read_img(fptr,TLONGLONG,1,nelements,&nulval,map_insou,&anynul,&status) ) cout<<"read table status= "<<status<<endl;
       fits_close_file(fptr, &status);  // close the file */

       // Print the values of some pixels to check we are reading correctly
       if(maptype==2) cout<<map_in2[NAXIS1*100+106]<<endl;
       if(maptype<2) cout<<map_in[NAXIS1*100+106]<<endl;
       if(maptype==3 && isEriMask==1) cout<<map_insou[NAXIS1*100+106]<<endl;
       if(maptype==3 && isEriMask==0) cout<<map_in[NAXIS1*100+106]<<endl;

// we use (l,m) because the source is spherical in (l,m) coordinates, and not in (RA,DEC) if the latter are treated in an euclidian space (eg computing distances with (RA*cosDEC)^2+DEC^2) 
// instead of in a spherical space (the sky)
       for(int idec=0; idec<NAXIS2; idec++){
         if(maptype==0&&ichannel==ichlow) dec.push_back(CDELT2*((CRPIX2-1)-idec));  // degree 
// above is actually m not dec // m=(sin(dec)*cos(dec0)-cos(dec)*sin(dec0)*cos(ra-ra0)))] however m-|dec-dec0|< arcsec in our maps (and we are not interested in DEC)
         for(int ira=0; ira<NAXIS1; ira++){
           if(idec==0&&maptype==0&&ichannel==ichlow) ra.push_back(CDELT1*(ira-(CRPIX1-1))); // degree
//above is actually l not ra // to get RA one should take arcsin(l/(15*cos(dec))) // note that ra and l may differ by few tens of arcsec at the boundary of the maps (but we are not interested in RA) 
           double flux= (map_in[idec*NAXIS1+ira]==map_in[idec*NAXIS1+ira]) ? map_in[idec*NAXIS1+ira]/(binsize*1.0): 0.; 
//           if(abs(flux)>1.e-30&&ichannel==ichlow&&imap<2) cout<<ichannel<<" "<<imap<<" "<<flux<<endl;
//         flux*=1.5/(nuobs*1.e9)*1.e3; // to go from 10^-20 erg/cm2/s/A to Jy (assuming a 1.5 A binning)
           if(maptype==0) {
             if(ichannel==ichlow) {map.push_back(flux); mapmodel.push_back(0.); maprms.push_back(rmsave); mapflag.push_back(0); mapsource.push_back(0.0); mapexp.push_back(0.0);}
             else {map[idec*NAXIS1+ira]+=flux;}  }
           if(maptype==1) {
             if(abs(flux)<1.e-100||isinf(flux)>0) mapflag[idec*ra.size()+ira]+=1;
//             if(isinf(flux)>0) mapflag[idec*ra.size()+ira]+=1;


             maprms[idec*NAXIS1+ira]+=flux;}
//         if(map[map.size()-1]>1.e-5) cout<< "Map "<<idec<<" "<<ira<<" "<<map_in[map.size()-1]<<endl;
           if(maptype==2) {mapflag[idec*NAXIS1+ira]+=map_in2[idec*NAXIS1+ira]; 
//             if(map_in2[idec*NAXIS1+ira]>0) cout<<"flag "<<map_in2[idec*NAXIS1+ira]<<" "<<mapflag[idec*NAXIS1+ira]<<endl;
           }
           if(maptype==3 && isEriMask==1) {
                if(map_insou[idec*NAXIS1+ira]==0){
                    mapsource[idec*NAXIS1+ira]+=1;
                }
//                cout<<map_insou[idec*NAXIS1+ira]<<" "<<mapsource[idec*NAXIS1+ira]<<endl;
//             if(map_insou[idec*NAXIS1+ira]>0) cout<<"mask "<<map_insou[idec*NAXIS1+ira]<<" "<<mapsource[idec*NAXIS1+ira]<<endl;
           }
           if(maptype==3 && isEriMask==0) {mapsource[idec*NAXIS1+ira]+=map_in[idec*NAXIS1+ira];
//             if(map_insou[idec*NAXIS1+ira]>0) cout<<"mask "<<map_insou[idec*NAXIS1+ira]<<" "<<mapsource[idec*NAXIS1+ira]<<endl;
           }

       }}

//       cout<<"About to delete these maps: "<<map_in<<" "<<&map_in<<" "<<map_in2<<" "<<&map_in2<<" "<<map_insou<<" "<<&map_insou<<endl;
       delete [] map_in; delete [] map_in2; delete [] map_insou;

       // Read our mask made with sextractor
       if(ichannel==ichlow && maptype==0 && usesextramask==1) {
         cout<<"Reading sextractor map"<<endl;
         float *map_ins = new float[nelements];
//       string filemask=data_path+"images/LeoT-approximate-segmentation.fits"; const char *filenamemask; filename=filemask.c_str();
//         string filemask=data_path+"images/Leo_sextractor_fin.fits"; const char *filenamemask; filenamemask=filemask.c_str();
         string filemask=data_path+"images/"+sextramaskname; const char *filenamemask; filenamemask=filemask.c_str();
         fitsfile *fptrmask;
         if( fits_open_image(&fptrmask,filenamemask,READONLY,&status) ) cout<<"read observational map: open status= "<<status<<" filenamemask = "<<filenamemask<<endl;
         felement=1;
         if( fits_read_img(fptrmask,TFLOAT,felement,nelements,&nulval,map_ins,&anynul,&status) ) cout<<"read table status= "<<status<<endl;
         for(int in=0; in<nelements; in++) mapsource[in]+=map_ins[in];
         delete [] map_ins;
         fits_close_file(fptrmask, &status);
         cout<<mapsource[NAXIS1*100+106]<<endl;
       }

       // Read exposure
       if(ichannel==ichlow && maptype==0 && useexp==1) {
         cout<<"Reading exposure map"<<endl;
         float *map_ins=new float[nelements];
         string fileexp = data_path+"images/"+expname; const char *filenameexp; filenameexp=fileexp.c_str();
         fitsfile *fptrexp;
         if( fits_open_image(&fptrexp,filenameexp,READONLY,&status) ) cout<<"read observational map: open status= "<<status<<" filenameexp = "<<filenameexp<<endl;
         felement=1;
         if( fits_read_img(fptrexp,TFLOAT,felement,nelements,&nulval,map_ins,&anynul,&status) ) cout<<"read table status= "<<status<<endl;

         for(int in=0; in<nelements; in++) {mapexp[in]+=map_ins[in];}
         delete [] map_ins;
         fits_close_file(fptrexp, &status);
         cout<<mapexp[NAXIS1*100+106]<<endl;
       }
  }
   cout<<"Maps loaded "<<(ra[1]-ra[0])*3600<<" "<<(dec[1]-dec[0])*3600<<" "<<ra[0]<<" "<<ra[ra.size()-1]<<" "<<dec[0]<<" "<<dec[dec.size()-1]<<endl;
//          abort();
//  Annuli(); abort();
 }
//  abort();

  // Exclude 20 pixels at the edges of the image and apply mask, flag, exposure
  int njump=20;
  int counter=0, counter_mapsource=0, counter_mapflag=0;
  for(int idec=0; idec<dec.size(); idec++)  for(int ira=0; ira<ra.size(); ira++) { 
    int in=idec*ra.size()+ira;
    double distvar=sqrt(ra[ira]*ra[ira] + dec[idec]*dec[idec])*60.;
//    if(distvar>maxdist) cout<<"distvar>maxdist"<<endl;
    if(idec<njump||ira<njump||idec>dec.size()-njump||ira>ra.size()-njump||abs(mapsource[in])>1.e-50 ||distvar>maxdist||mapflag[in]>=1 ||(useexp==1 && mapexp[in]<=explimit)) {
        map[in]=0.0;
        if (abs(mapsource[in])>1.e-50) counter_mapsource++;
        if (mapflag[in]>=1) counter_mapflag++;
        counter++;}
    else if(namedSph[idSph] == "eri" && field==2 && dec[idec] <=-0.006833333333333334 && ra[ira]>=0.00044444444444444447){
        map[in]=0.0; counter++;}
    else if(namedSph[idSph] == "ant" && idec>200){
//    cout<<idec<<" ";
        map[in]=0.0; counter++;}
    }
    cout<<"Number of pixels set to zero = "<<counter<<endl;
    cout<<"Number of abs(mapsource[in])>1.e-50 = "<<counter_mapsource<<endl;
    cout<<"Number ofmapflag[in]>=1 = "<<counter_mapflag<<endl;


  // Alternative derivation of the noise - Standard deviation computed in 10x10 boxes (or modify for different size)
/*  int ndec=dec.size(),nra=ra.size(), nbox=100;
  cout<<"Alternative derivation of the noise, nbox = "<<nbox<<endl;
  for(int idec=0; idec<ndec; idec++){
      int ndecmin=max(0,idec-nbox),ndecmax=min(idec+nbox,ndec);
      for(int ira=0; ira<nra; ira++){
        int nramin=max(0,ira-nbox),nramax=min(ira+nbox,nra);
        int ngrid=0;
        for(int in3=ndecmin; in3<ndecmax; in3++) for(int in4=nramin; in4<nramax; in4++) if(abs(map[in3*nra+in4])>1.e-50) ngrid++;
        double fluxave=0.0, rmsvar=0.0;
        for(int in3=ndecmin; in3<ndecmax; in3++) for(int in4=nramin; in4<nramax; in4++) if(abs(map[in3*nra+in4])>1.e-50) fluxave+=map[in3*nra+in4]/(ngrid*1.0);
        for(int in3=ndecmin; in3<ndecmax; in3++) for(int in4=nramin; in4<nramax; in4++) if(abs(map[in3*nra+in4])>1.e-50) rmsvar+=pow(fluxave-map[in3*nra+in4],2)/(ngrid*1.0);
        maprms[idec*nra+ira]=rmsvar;
      }
  } */

  // Initial estimates of Sflat. Media spaziale della mappa: aritmetica e pesata su 1/sigma^2
//  ofstream outfile(name+"_spectrum_ave_tmp.dat");
  double avespect=0.0, avevar=0.0, avestd=0.0, avevarinv=0.0, avespectw=0.0;
  int ncountpix=0;
  cout<<"map.size() = "<<map.size()<<endl;

  for(int in=0; in<map.size(); in++) if (abs(map[in])>1.e-50) {
    ncountpix++;
    avespect+=map[in];
    avevar+=maprms[in];
    avespectw+=map[in]/maprms[in];
    avevarinv+=1.0/maprms[in];
  }
  for(int in=0; in<map.size(); in++) if(abs(map[in])>1.e-50) avestd+=pow(map[in]-avespect/(ncountpix*1.0),2);
  cout<< "Map "<<(ichlow+ichsup)/2.0<<" "<<ncountpix<<endl;
  cout<<"avespect = "<<avespect<<" avevar = "<<avevar<<" sqrt(avestd/avevar) = "<<sqrt(avestd/avevar)<<endl;
  double aveave=avespectw/avevarinv; //avespect/(ncountpix*1.0); //
  cout<<"avespectw = "<<avespectw<<" aveave = "<<aveave<<endl;

  double aveerr=1.0/sqrt(avevarinv); //sqrt(avevar)/(ncountpix*1.0); //    sqrt(avestd)/(ncountpix*1.0);  //
//  double ratiob=(20.7*1.116e57*1.6022e-3/pow(3.0856e21,2)*1.e20*0.762*1.e-25/deltanuobs/(4.0*Pi)/pow(maxdist/(beamx),2))/(aveave+2.0*aveerr); // taking the D-factor to be 20.7 Msun/kpc^2 at 30 arcsec // the factor 0.762 is because we take the integral over the energy fwhm and not between -inf and +inf
//  outfile<<4699.57+(ichlow+ichsup)/2.0*1.25<<" "<<avespect/(ncountpix*1.0)<<" "<<avespectw/avevarinv<<" "<<sqrt(avevar)/(ncountpix*1.0)<<" "<<sqrt(avestd)/(ncountpix*1.0)<<" "<<1.0/sqrt(avevarinv)<<" "<<2.0*nuobs*2.0*Pi*6.5822e-7<<" "<<ratiob<<endl;
//  outfile.close();
  egb=avespectw/avevarinv; //avespect/(ncountpix*1.0); //  // Comment this line to use weighted average as initial guess for Sflat (when best fit is read from inputfile_cube.dat. Otherwise we calculate it on the fly setting egb=avespectw/avevarinv;). used as a first guess for egb
   
  for(int in=0; in<maprms.size(); in++) maprms[in]=sqrt(maprms[in]);
  return;
}

double dSph::getAverageEmission() {

  double avg = 0;
  int Npix = 0;

  for(int in=0; in<map.size(); in++) if (abs(map[in])>1.e-50) {
    Npix++;
    avg+=map[in];
  }

  avg = avg/Npix;
  cout<<"Npix = "<<Npix<<", avg = "<<avg<<endl;
  return avg;
}

void dSph::setup(vector<double> parvar, double& rescale) {
  /* routine to assign params of the scan and decide whether to do the computation or just rescale results of previous step */
  // parvar contains the value of the parameters, rescale is initialized to 0
  // parrec has dimension = number of parameters (defined in ConstParams), and all elements are initialized to 1e30 in the constructor
  rescale=parvar[0]/parrec[0]; // Ratio between current and previous value of the first parameter in compprior. Currently, Gamma
  if(abs(1-parrec[1]/parvar[1])>1.e-5) rescale=1; // 1-Ratio between previous and current values on the second parameter in compprior. Currently Sflat. If the ratio is not close enough to 1, we run the computation
  // rescale=1 implies to run the computation (par[1] is mchi or r0)
//  cout<<"In dSph::setup. rescale = "<<rescale<<endl;
  parrec=parvar;
  return;
}


void dSph::compute(double rescale, double& chi2, int& ichi2, double& time1, double& time2) {
/* Chi-sqare computation on either all pixels or ROI */
//  cout << "Computing model. egb = "<<egb << endl;

  ::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
  double egbmodel=egb,pbeam=1.0;
//  cout<<"egbmodel = "<<egbmodel<<endl;
  if(cube==1) egbmodel=egb*parrec[1]; // multiplied by Sflat
  if(rescale<1.e-10) egbmodel=egb; // when the function is called to evaluate the null case 
  if(abs(1-rescale)<1.e-5) {CRtabsph(diffsource);}
  else {for(int in=0; in<tabCRy.size(); in++) tabCRy[in]+=log(rescale);}
  Tools* T;
  int npixx=int(beamx*2.355/2.0*sqrt(Pi)/60.0/abs(ra[1]-ra[0])); // we need an approx to go from elliptical to square: take ellipse area=l1*l2 with l=sqrt(pi)*r, r=FWHM/2
  int npixy=int(beamy*2.355/2.0*sqrt(Pi)/60.0/abs(dec[1]-dec[0]));
//  cout<<" Number of independent points "<<map.size()/(1.0*npixx*npixy)<<" out of "<<map.size()<<endl;
  int waydof=2; // 0->all pixels; 1->"average over beam", ie all pixels but divide chi2 by npix; 2->only 1 pixel per beam 
  int ROIon=0; // 0->full map; 1->selected ROI;
  int npixgood=0,npixcons=0;
  double pixsize=abs(ra[1]-ra[0])*60.;
  double totem=0.0,totrms=0.0;
  double chi2var=0.0;
  ::chrono::steady_clock::duration elapsedTime = ::chrono::steady_clock::now() - startTime;

  ::chrono::steady_clock::time_point startTime1 = std::chrono::steady_clock::now();
  int ioffset, idx;
  for(int idec=0; idec<dec.size(); idec++) {
      ioffset = ra.size()*idec;
      if(idec%npixy==0||waydof!=2){
      for(int ira=0; ira<ra.size(); ira++){
        idx = ioffset + ira;
        if(ira%npixx==0||waydof!=2){
          double thetavar;
          if (ra_c != 0.0 && dec_c != 0.0) {
              thetavar=sqrt( (ra[ira]-ra0)*(ra[ira]-ra0) + (dec[idec]-dec0)*(dec[idec]-dec0))*60.; // arcmin
          }
          else{thetavar=sqrt(ra[ira]*ra[ira] + dec[idec]*dec[idec])*60.;} // arcmin       if(thetavar<=maxdist&&inuisance==0) npixcons+=1;
          if(abs(map[idx])>1.e-20 && maprms[idx]>1.e-20 && thetavar<=maxdist && abs(mapsource[idx])<1.e-30 && map[idx]>-3.* abs(maprms[idx])){
            ichi2+=1; npixgood+=1;
            double thlog=(thetavar>pixsize) ? log(thetavar) : log(pixsize);
            mapmodel[idx] = 0.0 + exp(T->linint(tabCRx,tabCRy,thlog)) + egbmodel; // put the PS-model instead of 0 if you add it
            chi2var+=pow((map[idx]-mapmodel[idx])/maprms[idx],2);
            if(isnan(chi2var)) {cout<<"Problem in chi2 "<<map[idx]<<" "<<maprms[idx]<<" "<<mapmodel[idx]<<endl; abort();}
            if(chi2var>1.e30){cout<<"chi2>1.e30 "<<mapmodel[idx]<<" "<<map[idx]<<" "<<maprms[idx]<<" "<<chi2var<<" "<<idec<<" "<<ira<<" "<<egbmodel<<endl;abort();}
            }
          }
    }}}
  chi2=chi2var;
  ::chrono::steady_clock::duration elapsedTime1 = ::chrono::steady_clock::now() - startTime1;

  if(waydof==1) chi2/=(1.0*npixx*npixy);  

  time1 = ::chrono::duration_cast< ::chrono::duration< double > >(elapsedTime).count();
  time2 = ::chrono::duration_cast< ::chrono::duration< double > >(elapsedTime1).count();

  return;
}



void dSph::CRtabsph(string namec) {
/* tabulation of diffuse intensity at different angles assuming spherical symmetry*/
  tabCRx.clear();
  tabCRy.clear(); 
  int ntheta=50;
  double thin=log(beamx/10.), thfin=log(asin(rintfin*0.9/dist)/aminrad); //arcmin
  for(int ith=0; ith<ntheta+1;ith++) tabCRx.push_back(thin+(thfin-thin)/(ntheta*1.)*ith); //log (arcmin)
  int nchunks=1;//floor(obs_bandwidth/obs_chunks+0.5);
  double res=0.;
  for(int inu=0; inu<nchunks; inu++){
    if(cube!=1) nuobs=obs_freq-obs_bandwidth/2.+obs_chunks*(inu+0.5);
    jload=0; // to reset emissivity-table
    for(int ith=0; ith<ntheta+1;ith++) {
      if(inu==0) {tabCRy.push_back(max(diffint(namec,exp(tabCRx[ith])),1.e-200)/(nchunks*1.));}
//      if(inu==0) {tabCRy.push_back(1.e-2*parrec[0]*exp(-pow(exp(tabCRx[ith])/parrec[1],2)));} // just for checking !!!
//      if(inu==0) {tabCRy.push_back(parrec[0]/1.e-25/(pow(180.,2)/4.0)*exp(-pow(exp(tabCRx[ith])/(180./2.355),2)/2));} // just for checking !!!
      else {tabCRy[ith]+=diffint(namec,exp(tabCRx[ith]))/(nchunks*1.);}
//      cout<<"Flux "<<ith<<" "<<exp(tabCRx[ith])<<" "<<diffint(namec,exp(tabCRx[ith]))/(nchunks*1.)<<" "<<tabCRy[ith]<<endl;
//      cout<<exp(tabCRx[ith])<<" "<<tabCRy[ith]*1.e3<<endl;//arcmin vs mJy
  }} 
  for(int ith=0; ith<ntheta+1;ith++) tabCRy[ith]=log(max(tabCRy[ith],1.e-200)); //abort();
  return; 
};


double dSph::diffint(string namec,double th0){
/* diffuse intensity I in units of Jy (having defined I= int dphi dtheta sintheta dl j/(4*Pi) );  th0 in arcmin */
  return Pi*pow(beamx*aminrad*2.355/2.0,2)*diffintlos(namec,th0*aminrad)/(4.*Pi);

  double resth=0, thmax=5.*beamx, thvec[1000]; 
  int ntheta; 
  ntheta=min(int(beamx/th0*101),999); // to increase the precision towards the center
  if(ntheta<11) ntheta=11; // must use ntheta > (thfin-thin)/beam (ie >10 if thmax=5*beam) !!!!
  if(beamx<2.e-3||beamy<2.e-3) cout<<"Warning: l.o.s. integral may not be accurate for such small angle "<<beamx<<endl;
  double sinthetacut=rintin/dist; // cut angle
  double thetasplit=asin(3.*sinthetacut); // inner angle at which the integral is split
  double onemcosthetasplit=(thetasplit<1.e-4) ? pow(thetasplit,2)/2.-pow(thetasplit,4)/24. : 1-cos(thetasplit);
  double thetasup=min(Pi,(th0+thmax)*aminrad); // upper boundary of integration in theta
  double thetainf=max(0.,(th0-thmax)*aminrad); // lower boundary of integration in theta
  double onemcosinf= (thetainf<1.e-4) ? pow(thetainf,2)/2.-pow(thetainf,4)/24.: 1-cos(thetainf);
  double onemcossup= (thetasup<1.e-4) ? pow(thetasup,2)/2.-pow(thetasup,4)/24.: 1-cos(thetasup);

  if(onemcosthetasplit>onemcosinf) {
    thvec[0]=onemcosinf;
    thvec[ntheta]=min(onemcosthetasplit,onemcossup);
    for(int ith=1; ith<ntheta+1;ith++){
      thvec[ith]=thvec[0]+(thvec[ntheta]-thvec[0])/(ntheta*1.)*ith;
      double thvar=(thvec[ith]+thvec[ith-1])/2.;
      resth+=(thvec[ith]-thvec[ith-1])*phifact(thvar,thmax*aminrad,th0*aminrad)*diffintlos(namec,acos(1-thvar));}} //Jy
  if(onemcossup>onemcosthetasplit) {
    thvec[0]=max(onemcosthetasplit,onemcosinf);
    thvec[ntheta]=onemcossup;
    for(int ith=1; ith<ntheta+1;ith++){
      thvec[ith]=thvec[0]+(thvec[ntheta]-thvec[0])/(ntheta*1.)*ith;
      double thvar=(thvec[ith]+thvec[ith-1])/2.;
//      cout<<"Integration "<<th0<<" "<<acos(1.-thvar)/aminrad<<" "<<thvar<<endl;
      resth+=(thvec[ith]-thvec[ith-1])*phifact(thvar,thmax*aminrad,th0*aminrad)*diffintlos(namec,acos(1-thvar));}} //Jy 
//      resth+=(thvec[ith]-thvec[ith-1])*parrec[0]*exp(-pow(thvar/(parrec[1]*aminrad),2)/2.0)*1.e3;}} //Jy // just for checking !!!
   return resth/(4.*Pi);};  


double dSph::phifact(double onemcostheta,double thmax, double th0) {
/* integral in dphi of the angular part of the integral for diffuse intensity I; output in rad; thmax,th0 in rad */
  double cosphi,phisup;
  if(onemcostheta==0||th0<1.e-15) {phisup=Pi;}        
  else{
    double stheta=sqrt(2*onemcostheta-onemcostheta*onemcostheta);
    double sintheta=(stheta>1.e-16) ? stheta : 1.e-16;
    double theta=acos(1-onemcostheta);
    double sinth0=(sin(th0)>1.e-16) ? sin(th0) : 1.e-16;
    if(abs(th0)<1.e-4&&abs(thmax)<1.e-4) {
      cosphi=((th0*th0-thmax*thmax)/2.*(1-(th0*th0+thmax*thmax)/12.)+onemcostheta*cos(th0))/sintheta/sinth0;}
    else if(abs(theta-th0)<1.e-4&&abs(thmax)<1.e-4) {
      cosphi=1-((thmax-theta+th0)/2.*(1.-(thmax*thmax+pow(theta-th0,2))/12.)*(thmax+theta-th0))/sintheta/sinth0;}
    else{ 
      cosphi=(cos(thmax)-(1-onemcostheta)*cos(th0))/sintheta/sinth0;}
    if(cosphi>1||cosphi<-1) {phisup=Pi;}
    else {phisup=abs(acos(cosphi));}}

  double phivec[100],resint=0.;
  int nphi=10;
  phivec[0]=-phisup;
  for(int ith=1; ith<nphi+1;ith++){
    phivec[ith]=-phisup+(2*phisup)/(nphi*1.)*ith;
    resint+=(phivec[ith]-phivec[ith-1])*gaweint((phivec[ith]+phivec[ith-1])/2.,onemcostheta,th0);}
  return resint;}


double dSph::gaweint(double phi, double onemcth, double th0) {
/* Gaussian (or elliptical) angular response for the telescope */
  double sintheta=sqrt(2.*onemcth-onemcth*onemcth);
  double onemcthnew= (th0<1.e-4) ? onemcth+(th0*th0/2.-pow(th0,4)/24.)*(1-onemcth)-sin(th0)*sintheta*cos(phi) : (1-cos(th0))+cos(th0)*onemcth-sin(th0)*sintheta*cos(phi);
  if(onemcthnew<0&&onemcthnew>-1.e-15) onemcthnew=0;
  if(onemcthnew>2.||onemcthnew<0) cerr<< "wrong onemcosthetanew in gaweint "<<onemcthnew<<endl;
  double tanphinew=sin(th0)*sin(phi)/(cos(phi)*sintheta*cos(th0)+sin(th0)*(1-onemcth));  // (or is it a minus sign???)
  double sinphinew2 = (abs(tanphinew)<1.) ? tanphinew*tanphinew/(1.+tanphinew*tanphinew) : 1./(1.+1./(tanphinew*tanphinew));
  double tanthetanew2=onemcthnew*(2.-onemcthnew)/(1.-onemcthnew*(2.-onemcthnew));
  double cosphinew2=1-sinphinew2;
/* choose here gaussian or elliptical response function */
  double weight=exp(-0.5*tanthetanew2/pow(tan(beamx*aminrad),2));
//  double weight=exp(-0.5*(tanthetanew2/pow(tan(beamx*aminrad),2)*cosphinew2+tanthetanew2/pow(tan(beamy*aminrad),2)*sinphinew2));
  return weight;}


double dSph::diffintlos(string namec,double thin){
/* diffuse intensity I along a given los specified by thin in units of Jy or similar (having defined I=int dl j)
   thin is in rad and is th0+thvar */
  if(usephi==0) {
//    cout<<"\nparrec.size() =" <<parrec.size()<<"  "<<"parrec[0] = "<<parrec[0]<<"  "<<"parrec[1] = "<<parrec[1]<<endl;
    DM targ(idSphrec, parrec, nuobs);
    double res=0.0,rrvec[1000];
    int nr=50;
    double rdistvar=dist*sin(thin); //kpc, distance of the l.o.s. from the center of the dSph 
    rrvec[0]=rintin;
    rrvec[nr]=sqrt(rintfin*rintfin-rdistvar*rdistvar);
    for(int ir=1; ir<nr+1; ir++){
      rrvec[ir]=exp(log(rrvec[0])+(log(rrvec[nr])-log(rrvec[0]))/(nr*1.)*ir); //kpc, distance on the line of sight from the max point
      double rrvar=sqrt(rrvec[ir]*rrvec[ir-1]+rdistvar*rdistvar); // kpc, distance from the center of the dSph
      double jj=targ.darkmajdec(rrvar, deltanuobs/1.e9) / 1.e3; // in units of 10^-20 erg/cm2/s/A/kpc
//      cout<<"diffintlos in loop "<<thin<<" "<<rrvar<<" "<<jj<<" "<<deltanuobs/1.e9<<endl;
      res+=jj*(rrvec[ir]-rrvec[ir-1]);}
//    cout<<"\ndiffintlos "<<rdistvar<<" "<<thin<<" "<<endl;
//    cout<<res<<" ";
//    double ress = parrec[0]/(4.0*Pi*deltanuobs)*1.116e57*1.6022e-3/pow(3.0856e21,2)*1.e20*1.048e5*1.782*pow(thin/aminrad/60.0,1.782)/(2.0*Pi*thin*thin)*4*Pi;
//    cout/*<<ress<<" "*/<<(ress/res)<<endl;

    return res;}

 else return parrec[0]/(4.0*Pi*deltanuobs)*1.116e57*1.6022e-3/pow(3.0856e21,2)*1.e20*1.048e5*1.782*pow(thin/aminrad/60.0,1.782)/(2.0*Pi*thin*thin)*4*Pi; //the factor 4*Pi cancels with the one in diffint

 // LeoT: signal = Gamma/(4Pi*deltanu)*Phi*4Pi (here in units of 10^-20 erg/cm2/s/A) with Phi~D0*a*theta^(a-2)/(2*Pi) 
// D-factor from extrapolation of Bonnivard et al. + analytic derivation // bf: D0=1.048e5, a=1.782; low: D0=4.346e4, a=1.696; up: D0=2.619e5, a=1.823; [D0]=Msun/kpc^2
// alternative D-factor derived from best-fit profiles in 1309.0815: D0=7.943e4, a=1.95 (Burkert) or a=1.73 (NFW)

};  

void dSph::GetOut() {

   cout<<" GetOut outputs "<<endl;
   if(optrun==1) {
     CRtabsph(diffsource);
     string filevar="outmodel_"; filevar.append(name); filevar.append(".dat");
     const char *outfilename; outfilename=filevar.c_str();
     ofstream outfile (outfilename);
     int npoints=100;
     double decmin=beamx/3./60.,decmax=asin(rintfin*0.9/dist)/aminrad/60.;
     outfile<<npoints<<endl;
     Tools* T;
     for(int idec=0; idec<100; idec++) {
       double dectmp=decmin*pow(10,idec/(npoints*1./log10(decmax/decmin)));
       double res=exp(T->linint(tabCRx,tabCRy,log(dectmp*60.)));
//       cout<<dectmp<<" "<<res<<endl;
       outfile<<dectmp<<" "<<res<<endl;}
     abort();
     outfile.close();}
   else{
   cout<<" Starting Remake "<<endl;
   long nelements=ra.size()*dec.size();
/* questa parte e' una ripetizione: se si puo' eliminarla!!*/ 
   double *array;
   array=new double[nelements];
   int icount=0;
   CRtabsph(diffsource);
   Tools* T;
   for(int idec=0; idec<dec.size(); idec++) {
     for(int ira=0; ira<ra.size(); ira++){
       double thetavar=sqrt(ra[ira]*ra[ira]+dec[idec]*dec[idec])*60.; // arcmin
       double thlog=(thetavar>1.e-1) ? log(thetavar) : log(1.e-1);
       array[ra.size()*idec+ira]=exp(T->linint(tabCRx,tabCRy,thlog));//+egb;
/*              double discdec=abs(idec-1.96*ira+3520),discra=abs(ira-(idec+3520)/1.96),ddisc=2./3600.*discdec*discra/sqrt(discdec*discdec+discra*discra);
              double A0disc=5.037e-3, d0disc=0.184,discmodel=A0disc*exp(-pow(ddisc/d0disc,2)/2.0);
              array[ra.size()*idec+ira]+=discmodel;*/
//       array[ra.size()*idec+ira]= (pixSN[ra.size()*idec+ira]==1) ? 1 : 0;
   }}
   cout<<" Remake done "<<endl;
  
   for(int in=0; in<2; in++) {
     int status=0;
     long fpixel=1,naxis = 2, naxes[2];
     double CRVAL1=ra0,CRVAL2=dec0;
     fitsfile *fptr; /* FITS file pointer */
     double CDELT1=ra[1]-ra[0];
     double CDELT2=-dec[1]+dec[0]; //opposite sign because table reflects image but dec grows in the opposite direction
     long CRPIX1=abs(ra[0]/CDELT1)+1;
     long CRPIX2=abs(dec[0]/CDELT2)+1;
     naxes[0]=ra.size();
     naxes[1]=dec.size();
     string filevar="!map_";
     filevar.append(namedSph[idSphrec]);
     if(in==0) filevar.append("_model.fits");
     if(in==1) filevar.append("_resid.fits");
     const char *outfile;
     outfile=filevar.c_str();
     if(in==1) for(long ip=0; ip<nelements; ip++) array[ip]=map[ip]-array[ip];
     fits_create_file(&fptr, outfile, &status); // create new file or overwrite existing one
     fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status); // Create the primary array image (32-bit float pixels)
     fits_write_img(fptr, TDOUBLE, fpixel, nelements, array, &status);
     fits_update_key(fptr, TLONG, "CRPIX1", &CRPIX1,"central pixel for RA", &status);
     fits_update_key(fptr, TLONG, "CRPIX2", &CRPIX2,"central pixel for DEC", &status);
     fits_update_key(fptr, TDOUBLE, "CRVAL1", &CRVAL1,"central value in RA", &status);
     fits_update_key(fptr, TDOUBLE, "CRVAL2", &CRVAL2,"central value in DEC", &status);
     fits_update_key(fptr, TDOUBLE, "CDELT1", &CDELT1,"step in RA", &status);
     fits_update_key(fptr, TDOUBLE, "CDELT2", &CDELT2,"step in DEC", &status);
     fits_close_file(fptr, &status);     // close the file
     fits_report_error(stderr, status);  // print out any error messages
   }
// aggiungere output distribuzione elettroni!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   cout<<" Output fits files produced "<<endl;}
  return;
}


void dSph::GetROI() {
    for(int in=0; in<map.size(); in++) pixSN.push_back(0);
    double Rston=0.;
// select the pixel with maximum S/N:
    int imax;
    for(int in=0; in<map.size(); in++) if(abs(map[in])>1.e-10&&maprms[in]>1.e-6) {
      double Rstonvar=mapmodel[in]/maprms[in];
      if(Rstonvar>Rston) {Rston=Rstonvar; imax=in;}}
    double sumsign=mapmodel[imax];
    double sumrms=pow(maprms[imax],2);
    pixSN[imax]=1;   
    int count1=1,count2=1,counttot=0;

    while(count1>0||count2>0){ // loop until the regions of pixels involved do not change
      count1=0;count2=0;
// add pixels which improve overall S/N
      for(int in=0; in<map.size(); in++) if(abs(map[in])>1.e-10&&maprms[in]>1.e-6) {
        if(pixSN[in]==0){
          double Rstonvar=(mapmodel[in]+sumsign)/sqrt(sumrms+pow(maprms[in],2));
          if(Rstonvar>Rston) {
            pixSN[in]=1;
            count1+=1;}}}
      sumsign=0.;
      sumrms=0.;
       for(int in=0; in<map.size(); in++) {
         if(pixSN[in]==1) {
           sumsign+=mapmodel[in];
           sumrms+=pow(maprms[in],2);}}
      Rston=sumsign/sqrt(sumrms);

// subract added pixels which degrade overall S/N
      for(int in=0; in<map.size(); in++){
        if(pixSN[in]==1){
          double Rstonvar=(sumsign-mapmodel[in])/sqrt(sumrms-pow(maprms[in],2));
          if(Rstonvar>Rston) {
            pixSN[in]=0;
            count2+=1;}}}
      sumsign=0.;
      sumrms=0.;
      counttot=0;
      for(int in=0; in<map.size(); in++) {
        if(pixSN[in]==1) {
          counttot+=1;
          sumsign+=mapmodel[in];
          sumrms+=pow(maprms[in],2);}}
     Rston=sumsign/sqrt(sumrms);
}

  cout<<"Total ratio signal to noise from ROI method "<<Rston<<" with pixels involved "<<counttot<<endl;
//  for(int in=0; in<map.size(); in++) if(pixSN[in]==1) cout<<in<<" pixel within ROI "<<mapmodel[in]<<" "<<maprms[in]<<endl;

  return;
}


void dSph::Annuli() {
// routine to compute flux density as a function of radial distance, making the average in spherical annuli
  cout<<"Computing flux density as a function of radial distance"<<endl;
  double radin=0.0, radout=100.0, radstep=5.0; // inner and outer radii and step size in arcmin
  int nann=int((radout-radin)/radstep+0.5);
  double beamvar=Pi*beamx*beamy*pow(2.355/2.0,2),normbeam=beamvar/abs((ra[1]-ra[0])*(dec[1]-dec[0]))/3600.0; // number of pixels in a synthesized beam
  double Iprof[100],rmsprof[100];
  ofstream annfile("annfile.dat");
  cout<<"Starting Annuli "<<normbeam<<endl;
  for (int in=0; in<nann; in++) {
    Iprof[in]=0.0; rmsprof[in]=0.0;
    double area=Pi*((1.0+2.0*in)*radstep*radstep+2.0*radin*radstep) ; // arcmin^2 // area annulus= Pi*(r2^2-r1^2)
    double npixarea=area/abs((ra[1]-ra[0])*(dec[1]-dec[0]))/3600.0; // number of pixels in the annulus area
    int intmp=0;
    for(int idec=0; idec<dec.size(); idec++) for(int ira=0; ira<ra.size(); ira++){
      double dist=sqrt(ra[ira]*ra[ira]+dec[idec]*dec[idec])*60.0; // dist from center // arcmin
      if(dist>=(radin+in*radstep)&&dist<(radin+(in+1)*radstep)&&abs(mapsource[ra.size()*idec+ira])<1.e-30&&map[ra.size()*idec+ira]>-3.*abs(maprms[ra.size()*idec+ira])&&abs(map[ra.size()*idec+ira])>1.e-20&&maprms[ra.size()*idec+ira]>1.e-20) {
        intmp++;
//        Iprof[in]+=map[ra.size()*idec+ira]/normbeam/area;  // Jy/arcmin^2
        if(isnan(mapmodel[ra.size()*idec+ira])) {cout<<"Problem in Annuli "<<mapmodel[ra.size()*idec+ira]<<endl; mapmodel[ra.size()*idec+ira]=0.0;}
        Iprof[in]+=mapmodel[ra.size()*idec+ira]/normbeam/area;  // Jy/arcmin^2
        rmsprof[in]+=maprms[ra.size()*idec+ira]/normbeam/area; }   // Jy/arcmin^2
    }
// we add rms linearly over the beam (and then take the average), but what about the sum over the area?
// is it quadratic (and if the rms ~ constant we can sum everything up linearly and then correct by a factor of sqrt(number of beams in the area), see below) or linear?
    double maskfact=intmp/(npixarea*1.0);
    Iprof[in]/=maskfact;
    rmsprof[in]/=max(1.0,sqrt(area*maskfact/beamvar)); 
    double stddev=0.;
// the variance should be computed comparing the global mean with the mean in each beam but summing over all pixels should be ok (assuming I ~ constant within the beam)
    for(int idec=0; idec<dec.size(); idec++) for(int ira=0; ira<ra.size(); ira++){
      double dist=sqrt(ra[ira]*ra[ira]+dec[idec]*dec[idec])*60.0; // dist from center // arcmin
      if(dist>=(radin+in*radstep)&&dist<(radin+(in+1)*radstep)&&abs(mapsource[ra.size()*idec+ira])<1.e-30&&map[ra.size()*idec+ira]>-3.*abs(maprms[ra.size()*idec+ira])&&abs(map[ra.size()*idec+ira])>1.e-20&&maprms[ra.size()*idec+ira]>1.e-20)  stddev+=pow(normbeam*(map[ra.size()*idec+ira]/normbeam/(area*maskfact)-Iprof[in]/(intmp*1.0)),2)/(intmp-1.0);} // variance   

    stddev=sqrt(stddev*intmp/normbeam+rmsprof[in]*rmsprof[in]); // adding the rms in quadrature to Nbeam * variance
    cout<<"Annuli "<<radin+in*radstep<<" "<<Iprof[in]<<" "<<rmsprof[in]<<" "<<stddev<<" "<<intmp<<endl;
    annfile<<radin+(in*1.0+0.5)*radstep<<" "<<Iprof[in]*1.e3<<" "<<stddev*1.e3<<endl; // arcmin, mJy/arcmin^2, mJy/arcmin^2
  }
  
  return;
}


/*****************************************************************************/






