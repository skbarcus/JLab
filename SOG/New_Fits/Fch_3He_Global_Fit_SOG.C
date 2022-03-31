#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <math.h>
#include "TMinuit.h"
#include "TCanvas.h" //Needed to compile canvas objects.
#include "TStopwatch.h" //Needed to compile stopwatch objects.
#include "TRandom.h" //Needed to compile random objects.
#include "TMarker.h" //Needed to compile marker objects.
#include "TLegend.h" //Needed to compile legend objects.
#include "TStyle.h" //Needed to compile gstyle objects.
#include "TSystem.h" //Needed to compile gsystem objects.

#include <TMath.h>
#include "Math/IFunction.h"
#include <cmath>
#include "Math/SpecFunc.h" //Access special functions like spherical Bessel.

#include "Math/Functor.h"
#include "Math/RichardsonDerivator.h"

using namespace std;

Double_t xyz = 5.5;
//Double_t truncate = 5.5;

Double_t pi = 3.141592654;
//Double_t deg2rad = pi/180.0; //For some odd reason this equals zero always when I'm working on the new fits on 3/29/22.
Double_t deg2rad = 0.017453293;
Double_t GeV2fm = 1./0.0389;            //Convert Q^2 units from GeV^2 to fm^-2.
Double_t hbar = 6.582*pow(10.0,-16.0);   //hbar in [eV*s].
Double_t C = 299792458.0;                //Speed of light [m/s]. 
Double_t e = 1.60217662E-19;             //Electron charge [C].
Double_t e2_nuclear = 1.4399643929E-3;             //Electron charge squared in nuclear units [GeV * fm].
Double_t alpha = 0.0072973525664;//1.0/137.0;              //Fine structure constant.
Double_t muHe3 = -2.1275*(3.0/2.0); //Diens has this 3/2 factor for some reason, but it fits the data much better.  //2*2.793-1.913 is too naive.
//Double_t muHe3 = 2.9788*(3.0/1.0); //3H
//Double_t mu3H = 2.9788*(3.0/1.0); //Magnetic moment of trinucleon (H3 or He3). NIST: http://physics.nist.gov/cgi-bin/cuu/Results?search_for=magnet+moment   //MCEEP Code for H3 and He3 eleastic FFs has magnetic moments multiplied by 3.0/Z. I don't know why but it works. Maybe it's a factor of A/Z?

Int_t loops = 1;
Int_t userand = 2;                       //0 = use predetermined Ri from Amroun. 1 = use random Ri in generated in a range around Amroun's. 2 = use random Ri, ngaus=12, generated in increments of 0.1 with larger possible spacing at greater radii. 3 = use predetermined Ri for the purposes of trying to tune the fit by hand. 4 = ngaus=8. 5 = ngaus=9. 6 = ngaus=10. 7 = ngaus=11.
Int_t usedifmin = 1;                     //0 = Remove some of the points in the diffractive minimum. 
Int_t showgaus = 1;
Int_t fitvars = 0;                       //0 = fit only Qi, 1 = fit R[i] and Qi, 2 = Fit R[i], Qi, and gamma.
Int_t fft = 0;                           //0 = don't use FFT to try to get a charge radii. 1 = do use FFT to extract a charge radii.
Int_t Amroun_Qi = 0;                     //1 = Override fitted Qi and use Amroun's values.
Int_t showplots = 1;
Int_t useFB = 0;                         //Turn on Fourier Bessel fit.
Int_t useFB_FM = 1;                      //0 = Turn on Fourier Bessel fit just for FC. 1 = Turn on Fourier Bessel fit attempting FC and FM.
Int_t improve = 0;                       //1 = use mnimpr() to check for other minima around the one MIGRAD finds.
Int_t MINOS = 0;                         //1 = use MINOS to calculate parameter errors. With ERRordef=30, npar=24, 10000 calls took about 1.5 hours and gave results only slightly different from intial parameter errors given. Several pars were hitting limits. 
Int_t optimize_Ri = 1;                   //1 = Have code loop over each Ri value shifting it 0.1 higher and 0.1 lower until chi2 stops improving.
Int_t bootstrap = 0;                     //0 = No bootstrapping. 1 = Using a fixed Ri set randomly select points in the dataset a number of times equal to the number of points in the dataset and then use those points for a fit.
Int_t npar = 48;                         //Number of parameters in fit.
Int_t ngaus = 12;                        //Number of Gaussians used to fit data.
Int_t ngaus_Amroun = 12;
Int_t nFB = 12;                          //Number of Fourrier-Bessel sums to use.
Double_t Z = 2.;                         //Atomic number He3.
//Double_t Z = 1.;                         //Atomic number H3.
Double_t A = 3.;                        //Mass number He3.
Double_t MtHe3 = 3.0160293*0.9315;         //Mass of He3 in GeV.
//Double_t MtHe3 = 3.0160492*0.9315;         //Mass of H3 in GeV.
//Double_t Mt3H = 3.0160492*0.9315;         //Mass of H3 in GeV.
//Double_t Gamma = 0.8*pow(2.0/3.0,0.5);   //Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.//For some odd reason this equals zero always when I'm working on the new fits on 3/29/22.
Double_t Gamma = 0.653197;
//Double_t E0 = 0.5084;                    //Initial e- energy GeV.
Double_t Ef = 0.;                        //Final e- energy GeV.
Double_t ymin = 30.;//30
Double_t ymax = 100.;//100
Double_t yminFF = 0.0001;//30
Double_t ymaxFF = 6.;
Double_t range = fabs(ymaxFF - yminFF);
Int_t n = 10000;
Int_t ndim = n+1;
Int_t npdraw = 10001;                     //Number of points to be used when drawing a function.
Double_t Truncate = 100.;                 //Truncate the histogram before inverse FFT. [fm^-2]
Int_t skip = 2.;                          //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                        //Counts number of lines in the data file. 
Int_t ncols;                             //Set how many columns of data we have in the data file.
char str[1000];                          //Variable to read lines of the data file..
Float_t thetatemp,qefftemp,sigexptemp,uncertaintytemp,E0temp;
Float_t theta[1000];                     //Angle in degrees.
Float_t qeff[1000];                      //q effective in fm^-1.
Float_t sigexp[1000];                    //Sigma experimental (cross section). Not sure on units yet.
Float_t uncertainty[1000];
Float_t E0[1000];

Int_t Amroun_pts = 57;                 //Dropped two points with no energy value give.
Int_t Collard_pts = 118;
Int_t Szlata_pts = 22;
Int_t Dunn_pts = 27;
Int_t Camsonne_pts = 18;               //Dropped two points with crazy Chi^2 values. Should reevaluate eventually.
Int_t Nakagawa_pts = 5;
Int_t my_pts = 1;
Int_t Arnold_pts = 11;                 //These XSs had to be calculated from A^1/2 function.
const Int_t datapts = 115;//115

Double_t m = 2.;

Double_t R[12] = {0.3,0.7,0.9,1.1,1.5,1.9,2.2,2.7,3.3,4.2,4.3,4.8};//Final 3He representative fit.

//Double_t R[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_Amroun[12] = {0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2}; //Amroun Fit
Double_t R_init[12] = {};
Double_t R_best[12] = {};
Double_t R_best_chi2 = 0;

//3He
Double_t Qich[12] = {0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.98962e-12};//3He final #30.
Double_t Qim[12] = {0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.77602e-12,0.0240927,8.94934e-12};//3He final #30.

Double_t Qich_Amroun[12] = {0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338};//3He
Double_t Qim_Amroun[12] = {0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.};//3He

//3H
//Double_t Qich_Amroun[12] = {0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121, 0.0, 0.0};//3H
//Double_t Qim_Amroun[12] = {0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077, 0.0, 0.0};//3H
Double_t Qich_best[12] = {};
Double_t Qim_best[12] = {};
Double_t av[24] = {9.9442E-3, 2.0829E-2, 1.8008E-2, 8.9117E-3, 2.3151E-3, 2.3263E-3, 2.5850E-3, 1.9014E-3, 1.2746E-3, 7.0446E-4, 3.0493E-4, 1.1389E-4};
Double_t averr[24] = {};
Double_t Qicherr[12]={}; 
Double_t Qimerr[12]={};
Float_t Q2[datapts];
Double_t Chi2[datapts]={};
Double_t residual[datapts]={};
Double_t xsfit[datapts]={};
Double_t fchfit[datapts]={};
Double_t Chi2_FB[datapts]={};
Double_t residual_FB[datapts]={};
Double_t FBfit[datapts]={};
Double_t E0_bs[datapts] = {};
Double_t theta_bs[datapts] = {};
Double_t sigexp_bs[datapts] = {};
Double_t uncertainty_bs[datapts] = {};
Double_t Q2_bs[datapts] = {};

Double_t  Qichtot = 0.;
Double_t  Qimtot = 0.;
Double_t amin = 0.;
Double_t maxQ2 = 0.;
TMarker *m1,*m2,*m3,*m4,*m5,*m6,*m7,*m8;

Float_t fch[datapts];
Float_t dfch[datapts];
Float_t Q2temp, fchtemp, dfchtemp;

//Define Fch with free pars for chi2 minimization fitting.

Double_t Fch(float Q2, Double_t *par)
{ 
  Double_t sumchtemp = 0.;
  Double_t fitch = 0.;
  Double_t Q2eff = Q2;//Need to add back E0 to make this correction again.
  
   //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 
      //cout<<"par["<<i<<"] = "<<par[i]<<endl;   //Starts with Amroun's pars I believe.
      //cout<<"R["<<i<<"] = "<<R[i]<<endl;
      //Fit just the Qi values using predetermined R[i] values.
      sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
      //cout<<"sumchtemp = "<<sumchtemp<<endl; //Fails here 3/29/2022.
      fitch =  fitch + sumchtemp;
      //cout<<"fitch["<<i<<"] = "<<fitch<<endl;
    }

  fitch =  fitch * exp(-0.25*Q2eff*pow(Gamma,2.0));
  return fitch;
}

//Define SOG FFs and XS.
Double_t XS(float E0, float theta, Double_t *par)
{
  //Double_t value=( (par[0]*par[0])/(x*x)-1)/ ( par[1]+par[2]*y-par[3]*y*y);
  //Double_t value = par[0] * x*x + par[1];
  Double_t val = 0.;
  Double_t mottxs = 0.;
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
  Double_t fitm = 0.;
  Double_t summtemp = 0.;
  
  Ef = E0/(1.0+2.0*E0*pow(sin(theta*deg2rad/2.0),2.0)/MtHe3);
  Double_t Q2 = 4.0*E0*Ef*pow(sin(theta*deg2rad/2.0),2.0) * GeV2fm;
  Double_t Q2eff = pow( pow(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*pow(GeV2fm,0.5)*1.12*pow(A,1.0/3.0))) ,2.0);   //Z=2 A=3
  
  Double_t W = E0 - Ef;
  //wHe3 = (Q2*1.0/GeV2fm)/(2.0*MtHe3);
  Double_t q2_3 = fabs(  pow(W,2.0)*GeV2fm - Q2eff  );        //Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
  Double_t eta = 1.0 + Q2eff/(4.0*pow(MtHe3,2.0)*GeV2fm);       //Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2.
  //cout<<"Q2 = "<<Q2<<"   Q2eff = "<<Q2eff<<"   W = "<<W<<"   q2_3 = "<<q2_3<<"   eta = "<<eta<<endl;
  Double_t Qtot = 1.0;
  Double_t Qtemp = 0.;

  /*
    for(Int_t i=0;i++,ngaus)
    {
    Qtemp = Qtemp + par[i];
    }*/
  
  //Calculate Mott XS.
  mottxs = (  (pow(Z,2.)*(Ef/E0)) * (pow(alpha,2.0)/(4.0*pow(E0,2.0)*pow(sin(theta*deg2rad/2.0),4.0)))*pow(cos(theta*deg2rad/2.0),2.0)  ) * 1.0/25.7;    //Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
  //cout<<"mottxs = "<<mottxs<<endl;

  //if(par[0]+par[1]+par[2]+par[3]+par[4]+par[5]+par[6]+par[7]+par[8]+par[9]+par[10]+par[11] == 1.)
  //{
  //cout<<"***************************************************************************"<<endl;
  /*
    for(Int_t i=0;i<2*ngaus;i++)
    {
    cout<<"par["<<i<<"] = "<<par[i]<<endl;
    }*/

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 
      //cout<<"par["<<i<<"] = "<<par[i]<<endl;   //Starts with Amroun's pars I believe.
      //cout<<"R["<<i<<"] = "<<R[i]<<endl;
      //Fit just the Qi values using predetermined R[i] values.
      sumchtemp = (par[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );
      //cout<<"sumchtemp = "<<sumchtemp<<endl; //Fails here 3/29/2022.
      fitch =  fitch + sumchtemp;
      //cout<<"fitch["<<i<<"] = "<<fitch<<endl;
    }

  //}
  //fitch =  fitch * exp(-0.25*Q2eff*pow(gamma,2.0));
  fitch =  fitch * exp(-0.25*Q2eff*pow(Gamma,2.0));
  
  //if(par[ngaus+0]+par[ngaus+1]+par[ngaus+2]+par[ngaus+3]+par[ngaus+4]+par[ngaus+5]+par[ngaus+6]+par[ngaus+7]+par[ngaus+8]+par[ngaus+9]+par[ngaus+10]+par[ngaus+11] == 1.)
  //{
  //Define SOG for magnetic FF.
  for(Int_t i=0; i<ngaus; i++)
    {
      //Fit just the Qi values using predetermined R[i] values.
      summtemp = (par[ngaus+i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) );	
      
      fitm = fitm + summtemp;
      //cout<<"fitm["<<i<<"] = "<<fitm<<endl;
    }
  //}
  fitm = fitm * exp(-0.25*Q2eff*pow(Gamma,2.0));   //For some reason had fabs(fitm).

  /*
    cout<<"E0 = "<<E0<<"   theta = "<<theta<<endl;
    cout<<"Ef = "<<Ef<<"   Q2 = "<<Q2<<"   Q2eff = "<<Q2eff<<"   W = "<<W<<"   q2_3 = "<<q2_3<<"   eta = "<<eta<<endl;
    cout<<"Mott XS = "<<mottxs<<endl;
    cout<<"fitch = "<<fitch<<"   fitm = "<<fitm<<endl;
  */
  /*  
      for(Int_t i=0;i<2*ngaus;i++)
      {
      cout<<"par["<<i<<"] = "<<par[i]<<endl;
      }
  */
  /*
    for(Int_t i=0;i<ngaus;i++)
    {
    cout<<"R["<<i<<"] = "<<R[i]<<endl;
    }
  */

  //Weird attempt to anchor Fch,m(0)=1. Relace first line in Amroun_3He_Data.txt. Previously the correct line was 0.3144	30	0.001923	0.000063459
  /*
  if(E0 == 0.3144 && theta == 30)
    {
      fitch = 1.0;
      fitm = 1.0;
    }
  */
  val = mottxs * (1./eta) * ( (Q2eff/q2_3)*pow(fitch,2.) + (pow(muHe3,2.0)*Q2eff/(2*pow(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + pow(tan(theta*deg2rad/2),2))*pow(fitm,2.) ); 
  //cout<<"XS = "<<val<<endl;
  return val;
}

//Create a Chi^2 function to minimize. 
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //const Int_t nbins = datapts;//177
  //Int_t i;
  //calculate chisquare
  Double_t chisq = 0;
  Double_t delta;
  Double_t res;
  if(bootstrap == 0)
    {
      for(Int_t i=0;i<datapts+1;i++) 
	{
	  if(i!=datapts)
	    {
	      //cout<<"XS["<<i<<"] = "<<XS(E0[i],theta[i],par)<<endl;
	      delta  = (fch[i]-Fch(Q2[i],par))/dfch[i];
	      chisq += delta*delta;
	      Chi2[i] = delta*delta;
	      //residual[i] = (sigexp[i] - XS(E0[i],theta[i],par))/sigexp[i]; 
	      //residual[i] = fabs(sigexp[i] - XS(E0[i],theta[i],par))/XS(E0[i],theta[i],par);
	      residual[i] = (fch[i] - Fch(Q2[i],par))/Fch(Q2[i],par); 
	      fchfit[i] = Fch(Q2[i],par);
	      //cout<<"xsfit["<<i<<"] = "<<xsfit[i]<<endl;
	    }
	  if(i==datapts)
	    {
	      delta  = (1-Fch(0.000001,par))/0.0001;
	      chisq += delta*delta;
	      Chi2[i] = delta*delta;
	      //residual[i] = (1 - Fch(0,par))/Fch(0,par); 
	      //fchfit[i] = Fch(0,par);
	    }
	  
	}
    }
  f = chisq;
}

//Define FFs for plotting purposes.
//Plot Charge FF Fch(Q) fm^-1.
Double_t ChFF(Double_t *Q, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      //Use SOG fit for C12 Qi coefficients and R[i] values. 
      //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

      //Convert to fm. Not sure I need to do this.
      //sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(Q[0]*pow(GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(Q[0]*pow(GeV2fm,0.5)*R[i])/(Q[0]*pow(GeV2fm,0.5)*R[i])) );
      sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(Q[0]*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(Q[0]*R[i])/(Q[0]*R[i])) );
	
      fitch = fitch + sumchtemp;
    }
  //Convert to fm. Not sure I need to do this.
  //fitch = fitch * exp(-0.25*pow(Q[0]*pow(GeV2fm,0.5),2.)*pow(gamma,2.0));
  fitch = fitch * exp(-0.25*pow(Q[0],2.)*pow(Gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}

//Plot Charge FF Fch(Q^2) fm^-2.
Double_t ChFF_Q2(Double_t *Q2, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    { 	
      //Use SOG fit for C12 Qi coefficients and R[i] values. 
      //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

      //Convert to fm. Not sure I need to do this.
      //sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0]*GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0]*GeV2fm,0.5)*R[i])/(pow(Q2[0]*GeV2fm,0.5)*R[i])) );
      sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );
	
      fitch = fitch + sumchtemp;
    }
  //Convert to fm. Not sure I need to do this.
  //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
  fitch = fitch * exp(-0.25*Q2[0]*pow(Gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}


//Plot Amroun's charge FF. No idea whay I can't just redefine Qi from ChFF_Q2.
Double_t ChFF_Q2_Amroun(Double_t *Q2, Double_t *par)
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;

  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus_Amroun; i++)
    { 	
      //Use SOG fit for C12 Qi coefficients and R[i] values. 
      //sumchtemp = (Qi[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0],0.5)*R[i])/(pow(Q2[0],0.5)*R[i])) );

      //Convert to fm. Not sure I need to do this.
      //sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( cos(pow(Q2[0]*GeV2fm,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (sin(pow(Q2[0]*GeV2fm,0.5)*R[i])/(pow(Q2[0]*GeV2fm,0.5)*R[i])) );
      sumchtemp = (Qich_Amroun[i]/(1.0+2.0*pow(R_Amroun[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2[0],0.5)*R_Amroun[i]) + (2.0*pow(R_Amroun[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2[0],0.5)*R_Amroun[i])/(pow(Q2[0],0.5)*R_Amroun[i])) );
	
      fitch = fitch + sumchtemp;
    }
  //Convert to fm. Not sure I need to do this.
  //fitch = fitch * exp(-0.25*Q2[0]*GeV2fm*pow(gamma,2.0));
  fitch = fitch * exp(-0.25*Q2[0]*pow(Gamma,2.0));
  fitch = fabs(fitch);
  return fitch;
}


Double_t fitg(Double_t *Q, Double_t *par)
{
  Double_t val = 0.;
  
  //Show Gaussian Part of FFs.
  val = (par[0]/(1.0+2.0*pow(par[1],2.0)/pow(Gamma,2.0))) * ( cos(Q[0]*par[1]) + (2.0*pow(par[1],2.0)/pow(Gamma,2.0)) * (sin(Q[0]*par[1])/(Q[0]*par[1])) );
  
  val = val * exp(-0.25*pow(Q[0],2.)*pow(Gamma,2.0));

  return val;
}

Double_t fitg_rho(Double_t *r, Double_t *par)
{
  Double_t val = 0.;

  //Show Gaussian part of rho.
  val = par[0]/( 1+2*pow(par[1],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-par[1]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+par[1]),2.)/pow(Gamma,2.) )  );

  val = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * val;

  return val;
}

//Define the charge density from I. Sick. 
Double_t rho_ch(Double_t *r, Double_t *par)
{
  Double_t rho = 0;
  Double_t rho_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(Gamma,2.) )  );
      rho = rho + rho_temp;
    }

  rho = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * rho; //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho;
}

//Create a function that can be integrated to check that the normilaization to Ze is correct.
Double_t rho_ch_int(Double_t *r, Double_t *par)
{
  Double_t rho_int = 0;
  Double_t rho_int_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_int_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(Gamma,2.) )  );
      rho_int = rho_int + rho_int_temp;
    }

  rho_int = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * rho_int * 4*pi*pow(r[0],2.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho_int;
}

//Create a function to calculate rms radius.
Double_t rho_rms(Double_t *r, Double_t *par)
{
  Double_t rho_rms = 0;
  Double_t rho_rms_temp = 0;
   
  for(Int_t i=0;i<ngaus;i++)
    {
      rho_rms_temp = Qich[i]/( 1+2*pow(R[i],2.)/pow(Gamma,2.) ) * (  exp( -pow((r[0]-R[i]),2.)/pow(Gamma,2.) ) + exp( -pow((r[0]+R[i]),2.)/pow(Gamma,2.) )  );
      rho_rms = rho_rms + rho_rms_temp;
    }

  rho_rms = Z/(2*pow(pi,1.5)*pow(Gamma,3.)) * rho_rms * 4*pi*pow(r[0],4.); //Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.

  return rho_rms;
}

Double_t ChFF_Deriv(Double_t Q2) 
{
  Double_t fitch = 0.;
  Double_t sumchtemp = 0.;
   
  //Define SOG for charge FF.
  for(Int_t i=0; i<ngaus; i++)
    {
      sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(Gamma,2.0))) * ( cos(pow(Q2,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(Gamma,2.0)) * (sin(pow(Q2,0.5)*R[i])/(pow(Q2,0.5)*R[i])) );
      fitch = fitch + sumchtemp;
    }
  fitch = fitch * exp(-0.25*Q2*pow(Gamma,2.0));
  //fitch = fabs(fitch);
  return fitch;
}


void Fch_3He_Global_Fit_SOG()
{

  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  FILE *fp;
  fp = fopen("/home/skbarcus/JLab/SOG/New_Fits/3He_Fch_Data.txt","r");

   //Read in data.
  while (1) {
    //Skips the first 5 lines of the file. 
    if (nlines < skip)
      {
	fgets(str,1000,fp);
	nlines++;
      }
    //Reads the two columns of data into x and y.
    else
      {
	//Read in the number of columns of data in your data file. 
	ncols = fscanf(fp,"%f %f %f",&Q2temp, &fchtemp, &dfchtemp);
	if (ncols < 0) break;   
	//cout<<"ncols = "<<ncols<<endl;
	//cout<<thetatemp<<"   "<<qefftemp<<"   "<<sigexptemp<<"   "<<uncertaintytemp<<endl;
	Q2[nlines-skip] = Q2temp;
	fch[nlines-skip] = fchtemp;
	dfch[nlines-skip] = dfchtemp; 

	//Q2[nlines-skip] = 4 * E0[nlines-skip] * (E0[nlines-skip]/(1.0+2.0*E0[nlines-skip]*pow(sin(theta[nlines-skip]*deg2rad/2.0),2.0)/MtHe3)) * pow(sin(theta[nlines-skip]*deg2rad/2.0),2.) * GeV2fm;

	cout<<"Q2["<<nlines-skip<<"] = "<<Q2[nlines-skip]<<"   fch["<<nlines-skip<<"] = "<<fch[nlines-skip]<<"   dfch["<<nlines-skip<<"] = "<<dfch[nlines-skip]<<endl;

	nlines++;
      }
  }

  fclose(fp);

 //Create an output file to store fit results.
  std::ofstream output ("Fch_Fit_Pars.txt", std::ofstream::out);
  output<<"Chi2   rChi2   BIC   AIC    Qichtot    R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch"<<endl;

  //Begin loop over fit with different Ri values each time.
  for(Int_t q=0;q<loops;q++)
    {

      if(userand == 2)
	{
	  //Generate random R[i] values. 
	  Double_t d = 0.49;
	  Double_t step = 0.5;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand = new TF1("rand","x",0.,.01);
	  R[0] = 0.1;
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand1 = new TF1("rand1","x",3.,4.);
	  R[1] = TMath::Nint(rand1->GetRandom())/10.+R[0];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand2 = new TF1("rand2","x",3.,4.);
	  R[2] = TMath::Nint(rand2->GetRandom())/10.+R[1];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand3 = new TF1("rand3","x",3.,4.);
	  R[3] = TMath::Nint(rand3->GetRandom())/10.+R[2];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand4 = new TF1("rand4","x",3.,4.);
	  R[4] = TMath::Nint(rand4->GetRandom())/10.+R[3];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand5 = new TF1("rand5","x",3.,4.);
	  R[5] = TMath::Nint(rand5->GetRandom())/10.+R[4];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand6 = new TF1("rand6","x",3.,4.);
	  R[6] = TMath::Nint(rand6->GetRandom())/10.+R[5];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand7 = new TF1("rand7","x",5.,6.);
	  R[7] = TMath::Nint(rand7->GetRandom())/10.+R[6];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand8 = new TF1("rand8","x",5.,6.);
	  R[8] = TMath::Nint(rand8->GetRandom())/10.+R[7];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand9 = new TF1("rand9","x",5.,6.);
	  R[9] = TMath::Nint(rand9->GetRandom())/10.+R[8];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand10 = new TF1("rand10","x",5.,6.);
	  R[10] = TMath::Nint(rand10->GetRandom())/10.+R[9];
	  gRandom->SetSeed(0);                    //Sets new random seed.
	  TF1 *rand11 = new TF1("rand11","x",5.,6.);
	  R[11] = TMath::Nint(rand11->GetRandom())/10.+R[10];
	}

      //Initiate Minuit for minimization.
      TMinuit *gMinuit = new TMinuit(24);  //initialize TMinuit with a maximum of 24 params
      gMinuit->SetFCN(fcn);

      Double_t arglist[10];
      Int_t ierflg = 0;

      arglist[0] = 30.; //1 is for simple chi^2. For multiparameter errors this needs to be increased. 30 adds ~ 1 min.
      gMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //Set the ERRordef or UP. 
  
      //Set step sizes.
      static Double_t stepsize[4] = {0.001 , 0.1 , 0.01 , 0.00001};
  
      //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
      for(Int_t i=0;i<ngaus;i++)
	{
	  gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], 0.,1.,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], 0.,0.,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.001,Qich[i]+0.001,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.000000000001,Qich[i]+0.000000000001,ierflg);
	  //gMinuit->mnparm(i, Form("Qich%d",i+1), Qich[i], stepsize[0], Qich[i]-0.05,Qich[i]+0.05,ierflg);
	}
  
      // Now ready for minimization step
      arglist[0] = 10000.;//Max calls. 50000.
      arglist[1] = 0.1;//Tolerance for convergance. 1 seems to give the same results as 0.1.
      //cout<<"Sup1"<<endl;
      gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
      if(improve == 1)
	{
	  gMinuit->mnimpr(); //Check for other minima to see if we're trapped in a local minima.
	}

      //Implement MINOS to calculate multiparameter errors. arglist[0] is giving # of calls again.
      if(MINOS == 1)
	{
	  gMinuit->mnexcm("MINOS", arglist, 1, ierflg);
	}

      //cout<<"Sup2"<<endl;
      // Print results
      Double_t edm,errdef;//Moved amin to global.
      Int_t nvpar,nparx,icstat;
      gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
      //gMinuit->mnprin(3,amin);

//Attempt to optimize the Ri by varying them systematically.
      if(optimize_Ri == 1)
	{
	  Int_t chi2_init = amin;
	  R_best_chi2 = amin;
	  cout<<"Initial Chi^2 = "<<amin<<endl;
	  R_init[0] = R[0];
	  R_best[0] = R[0];

	  //Set initial best Ri and Qi values to the initial Ri values.
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      gMinuit->GetParameter(i,Qich[i],Qicherr[i]);
	    }
	  for(Int_t i=0;i<ngaus;i++)
	    {
	      R_best[i] = R[i];
	      Qich_best[i] = Qich[i];
	    }

	  for(Int_t i=0;i<ngaus;i++)
	    {
	      Int_t chi2_better = 1;  //Test if the change improved chi2.
	      R_init[i] = R[i];
	      cout<<"*********************************************************"<<endl;
	      cout<<"Optimizing Initial R["<<i<<"] = "<<R[i]<<endl;
	  
	      //Check for better Ri values below the initial Ri value.
	      while(chi2_better == 1)
		{
		  //amin = 0.;
		  //Protect against R[i]=0. while still basically testing R[i]=0.
		  if(R[i]==0.1)
		    {
		      R[i] = 0.0001; 
		    }
		  else
		    {
		      R[i] = R[i] - 0.1;
		    }
		  cout<<"R["<<i<<"] set to "<<R[i]<<endl;
	      
		  //arglist[0] = 30.; //1 is for simple chi^2. For multiparameter errors this needs to be increased. 30 adds ~ 1 min.
		  //gMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //Set the ERRordef or UP.
	      
		  //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //gMinuit->mnparm(j, Form("Qich%d",j+1), Qich[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are Amroun's Qi.
		      gMinuit->mnparm(j, Form("Qich%d",j+1), Qich_best[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are best fit Qi values so far.
		    }
	      
		  arglist[0] = 10000.;//Max calls. 50000.
		  arglist[1] = 0.1;//Tolerance for convergance. 1 seems to give the same results as 0.1.
		  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
		  if(improve == 1)
		    {
		      gMinuit->mnimpr();
		    }
		  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		  cout<<"- Updated R["<<i<<"] = "<<R[i]<<"   New Chi^2 = "<<amin<<endl;
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //Set Qi values equal to the fitted parameters.
		      gMinuit->GetParameter(j,Qich[j],Qicherr[j]);
		      cout<<"R["<<j<<"] = "<<R[j]<<"   Qich["<<j<<"] = "<<Qich[j]<<endl;
		    }
		  if(amin<R_best_chi2)
		    {
		      for(Int_t j=0;j<ngaus;j++)
			{
			  Qich_best[j] = Qich[j];
			}
		      R_best[i] = R[i]; //Store best Ri value thus far.
		      R_best_chi2 = amin; //Store the best updated chi2 so far.
		      chi2_better = 1;
		    }
		  else
		    {
		      R[i] = R_best[i]; //No improvement so reset Ri to the best previous value.
		      cout<<"No improvement found. Setting R["<<i<<"] to best value = "<<R_best[i]<<".   With best Chi^2 = "<<R_best_chi2<<"."<<endl;
		      chi2_better = 0;
		    }
		} 
	      //Check for better Ri values above the initial Ri value.
	      chi2_better = 1;
	      R[i] = R_init[i];  //Start from initial Ri and move up now.
	      while(chi2_better == 1)
		{
		  R[i] = R[i] + 0.1;
		  cout<<"R["<<i<<"] set to "<<R[i]<<endl;

		  //Set starting guesses for parameters. (Use Amroun's SOG parameters.)
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //gMinuit->mnparm(j, Form("Qich%d",j+1), Qich[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are Amroun's Qi.
		      gMinuit->mnparm(j, Form("Qich%d",j+1), Qich_best[j], stepsize[0], 0.,1.,ierflg); //Starting guesses are best fit Qi values so far.
		    }

		  arglist[0] = 10000.;//Max calls. 50000.
		  arglist[1] = 0.1;//Tolerance for convergance. 1 seems to give the same results as 0.1.
		  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
		  if(improve == 1)
		    {
		      gMinuit->mnimpr();
		    }
		  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		  cout<<"+ Updated R["<<i<<"] = "<<R[i]<<"   New Chi^2 = "<<amin<<endl;
		  for(Int_t j=0;j<ngaus;j++)
		    {
		      //Set Qi values equal to the fitted parameters.
		      gMinuit->GetParameter(j,Qich[j],Qicherr[j]);
		      cout<<"R["<<j<<"] = "<<R[j]<<"   Qich["<<j<<"] = "<<Qich[j]<<endl;
		    }
		  if(amin<R_best_chi2)
		    {
		      for(Int_t j=0;j<ngaus;j++)
			{
			  Qich_best[j] = Qich[j];
			}
		      R_best[i] = R[i]; //Store best Ri value thus far.
		      R_best_chi2 = amin; //Store the best updated chi2 so far.
		      chi2_better = 1;
		    }
		  else
		    {
		      R[i] = R_best[i]; //No improvement so reset Ri to the best previous value.
		      cout<<"No improvement found. Setting R["<<i<<"] to final best value = "<<R_best[i]<<".   With best Chi^2 = "<<R_best_chi2<<"."<<endl;
		      chi2_better = 0;
		    }
		} 
	    }
	  //Final fit using all of the Ri_best values.
	  cout<<"Preforming final fit with optimized R[i] values."<<endl;
	  for(Int_t j=0;j<ngaus;j++)
	    {
	      cout<<"R["<<j<<"] = "<<R[j]<<endl;
	    }
	  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
	  if(improve == 1)
	    {
	      gMinuit->mnimpr();
	    }
	  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	}//End optimize Ri.

      Qichtot = 0; //Reset to zero in case looping over several fits.
      //Calculate the sum of the Qich parameters.
      for(Int_t i=0;i<ngaus;i++)
	{
	  //Qich[i] = fxs0->GetParameter(i);
	  //Qich[i] = gMinuit->GetParameter(i,1.,1.);
	  gMinuit->GetParameter(i,Qich[i],Qicherr[i]);
	  Qichtot = Qichtot + Qich[i];
	}

      //Fill text file with fcn (chi2), Qichtot, Qimtot, Ri, Qich, Qim.
      output<<amin<<" "<<amin/(datapts-2*ngaus-1)<<" "<<datapts*TMath::Log(amin/datapts)+TMath::Log(datapts)*ngaus<<" "<<datapts*TMath::Log(amin/datapts)+2*ngaus<<" "<<Qichtot<<" "<<R[0]<<" "<<R[1]<<" "<<R[2]<<" "<<R[3]<<" "<<R[4]<<" "<<R[5]<<" "<<R[6]<<" "<<R[7]<<" "<<R[8]<<" "<<R[9]<<" "<<R[10]<<" "<<R[11]<<" "<<Qich[0]<<" "<<Qich[1]<<" "<<Qich[2]<<" "<<Qich[3]<<" "<<Qich[4]<<" "<<Qich[5]<<" "<<Qich[6]<<" "<<Qich[7]<<" "<<Qich[8]<<" "<<Qich[9]<<" "<<Qich[10]<<" "<<Qich[11]; 
    }

  output.close();


  st->Stop();
  cout<<"*********************************************"<<endl;
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}

int main()
{
  //Double_t xyz = 5.5;
  Fch_3He_Global_Fit_SOG();
}
