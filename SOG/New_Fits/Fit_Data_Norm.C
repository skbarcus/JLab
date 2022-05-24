
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

/*
const int nvar=6; //Number of free parameters in the fit.
const int n7 = 33; //This definition lets you use the variable as the number of array elements so you don't need to hard code it in each time.
Double_t Q[n7]  ={1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 0.290, 0.386, 0.450, 0.481, 0.600, 0.670, 0.775, 0.850, 0.900, 0.951, 1.000, 0.051, 0.098, 0.258, 0.571, 1.05, 1.51, 2.21, 2.98, 0.3, 0.5, 0.7, 0.9};
Double_t F[n7]  ={0.622, 0.503,0.387, 0.312, 0.267, 0.215, 0.175, 0.137, 0.118, 0.0758, 0.885, 0.863, 0.769, 0.784, 0.79, 0.789, 0.769, 0.743, 0.603, 0.682, 0.579, 0.967, 1.039, 0.95, 0.826, 0.697, 0.548, 0.387, 0.282, 0.872, 0.782, 0.785, 0.659};
Int_t Dataset[n7] = {1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1};//Array marking which dataset each data point belongs too. These dataset associations are made up just for an example/test.
//Int_t Dataset[n7] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1};
// create the error arrays
Double_t eQ[n7] = {0};
Double_t eF[n7] ={0.014, 0.012, 0.010, 0.008, 0.007, 0.006, 0.005, 0.004, 0.005, 0.004, 0.062, 0.053, 0.082, 0.071, 0.064, 0.0234, 0.088, 0.044, 0.071, 0.054, 0.075, 0.047, 0.046, 0.029, 0.028, 0.027, 0.025, 0.028, 0.032, 0.047, 0.046, 0.038, 0.035};
*/


const int nvar=6; //Number of free parameters in the fit.
const int n7 = 12; //This definition lets you use the variable as the number of array elements so you don't need to hard code it in each time.
Double_t Q[n7]  ={1,2,3,4,5,6,7,8,9,10,11,12};
//Double_t F[n7]  ={1,4,9,16,25,36};
Double_t dA = 0.2;
Double_t dB = -0.1;
Double_t NormA = 1+dA;
Double_t NormB = 1+dB;
//Double_t F[n7]  ={1*NormA,4*NormB,9*NormA,16*NormB,25*NormA,36*NormB};
Double_t F[n7] = {pow(Q[0],2)*NormA,pow(Q[1],2)*NormB,pow(Q[2],2)*NormA,pow(Q[3],2)*NormB,pow(Q[4],2)*NormA,pow(Q[5],2)*NormB,pow(Q[6],2)*NormA,pow(Q[7],2)*NormB,pow(Q[8],2)*NormA,pow(Q[9],2)*NormB,pow(Q[10],2)*NormA,pow(Q[11],2)*NormB};
Int_t Dataset[n7] = {1,2,1,2,1,2,1,2,1,2,1,2};
Double_t eQ[n7] = {0};
//Double_t eF[n7] = {0.2,0.4,1.8,1.6,5,3.6};
Double_t eF[n7] = {pow(Q[0],2)*abs(dA*0.3),pow(Q[1],2)*abs(dB*0.3),pow(Q[2],2)*abs(dA*0.3),pow(Q[3],2)*abs(dB*0.3),pow(Q[4],2)*abs(dA*0.3),pow(Q[5],2)*abs(dB*0.3),pow(Q[6],2)*abs(dA*0.3),pow(Q[7],2)*abs(dB*0.3),pow(Q[8],2)*abs(dA*0.3),pow(Q[9],2)*abs(dB*0.3),pow(Q[10],2)*abs(dA*0.3),pow(Q[11],2)*abs(dB*0.3)};


// fixed chi2 parameters with prof. Long
Double_t poly3(Double_t *X, Double_t *par)
{
  Double_t val = 0.;

  val = par[0] + par[1]*X[0] + par[2]*pow(X[0],2) + par[3]*pow(X[0],3);

  return val;
}

Double_t poly3_norm(Double_t *X, Int_t Dataset, Double_t *par)
{
  //When adding floating normalizations for each individual dataset you'll need an if statement to perform the polynomial fit uniquely for each dataset. Basically each data point in a dataset is multiplied by a unique constant shared by that whole dataset. You may want to pass this function a number representing the dataset. Something like:
  //if dataset == 1
  //val = par[4] * ( par[0] + par[1]*X[0] + par[2]*pow(X[0],2) + par[3]*pow(X[0],3) );
  //if dataset == 2
  //val = par[5] * ( par[0] + par[1]*X[0] + par[2]*pow(X[0],2) + par[3]*pow(X[0],3) );
  //This way the polynomial coefficients are shared between data sets, but each form factor value is multiplied by a floating normalization parameter.
  //You'll also need to change the number of free parameters in TF1 and Minuit to represent the new free parameters.
  Double_t val = 0.;
  //Double_t c0;
  //Double_t c1;
  
  if(Dataset==1)
    {
      //c0 = 1;
      //c1 = 0;
      val = par[4] * poly3(X,par);
    }
  else if(Dataset==2)
    {
      //c0 = 0;
      //c1 = 1;
      val = par[5] * poly3(X,par);
    }
  
  //val =  c0 * par[4] * ( par[0] + par[1]*X[0] + par[2]*pow(X[0],2) + par[3]*pow(X[0],3) ) + c1 * par[5] * ( par[0] + par[1]*X[0] + par[2]*pow(X[0],2) + par[3]*pow(X[0],3) ) ;

  //val =  c0 * par[4] * poly3(X,par) + c1 * par[5] * poly3(X,par) ;
  
  return val;
}

Double_t poly3_norm1(Double_t *X, Int_t Dataset, Double_t *par)
{
  Double_t val = 0.;

  val = par[4] * poly3(X,par);
  
  //cout<<"par[4] = "<<par[4]<<"    par[2] = "<<par[2]<<endl;

  return val;
}

Double_t poly3_norm2(Double_t *X, Int_t Dataset, Double_t *par)
{
  Double_t val = 0.;

  val = par[5] * poly3(X,par);

  //cout<<"par[5] = "<<par[5]<<"    par[2] = "<<par[2]<<endl;

  return val;
}

void FCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t iflag)
//Define chi2 function to minimize.
{
    Double_t  chi2[n7] = {0};
    Double_t  A[n7]= {0};
    Double_t  A1[n7]= {0};
    Double_t   B[n7]= {0};
    Double_t   P[n7]= {0};
    Double_t  totcount= 0;
    Double_t count;
    for (Int_t i=0; i<n7; i++)
      {
	P[i]= poly3_norm(&Q[i],Dataset[i],par);//Attempts to find dataset normalizations.
	/*
	if(Dataset[i]==1)
	  {
	    P[i]= poly3_norm1(&Q[i],Dataset[i],par);
	  }
	else if(Dataset[i]==2)
	  {
	    P[i]= poly3_norm2(&Q[i],Dataset[i],par);
	  }
	*/
	A[i]= P[i]- F[i];
	A1[i]= pow(A[i],2);
	B[i]= pow(eF[i],2);
	chi2[i] = A1[i] /B[i];
	count= chi2[i];
	totcount= totcount+count;  
      }
    
    f = totcount;
    //cout<<"fcn = "<<f<<endl;
}
   
void Fit_Data_Norm()
{
  //  TMinuit //
  TMinuit *gMinuit = new TMinuit(nvar);  //initialize TMinuit with a maximum number of parameters. 
  gMinuit->SetFCN(FCN);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  //static Double_t vstart[nvar] = {1, -4.21077e-01 , 0.0713488  , -0.00447744, 1, 1};
  //static Double_t vstart[nvar] = {1, -0.4, 0.07, -0.005, 1.0, 1.0};
  static Double_t vstart[nvar] = {0.0, 0.0, 1., 0.0, 1.3, 0.8};
  //static Double_t vstart[nvar] = {0.0, 0.0, 1.05056, 0.0, 1.20472, 0.856687};
  //static Double_t vstart[nvar] = {-0.89894, -0.150002, 0.974971, -0.00416685, 1, 1};
  //static Double_t vstart[nvar] = {0.966523, -0.400214, 0.0697626, -0.00526376, 1.0, 1.0};
  //static Double_t vstart[nvar] = {1, -0.400214, 0.0697626, -0.00526376, 1.0, 1.0};
  //static Double_t vstart[nvar] = {1, 1, 1, 1};
  static Double_t step[4] = {0.001 , 0.01 , 0.1 , 1.0};
  gMinuit->mnparm(0, "a1", vstart[0], step[1], -10,10,ierflg);
  gMinuit->mnparm(1, "a2", vstart[1], step[1], -10,10,ierflg);
  gMinuit->mnparm(2, "a3", vstart[2], step[1], 0.7,1.3,ierflg);
  gMinuit->mnparm(3, "a4", vstart[3], step[1], -10,10,ierflg);
  gMinuit->mnparm(4, "a5", vstart[4], step[1], 1.,1.4,ierflg);
  gMinuit->mnparm(5, "a6", vstart[5], step[1], 0.7,1.,ierflg);
  
  //gMinuit->FixParameter(0);
  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);
  //gMinuit->FixParameter(4);
  //gMinuit->FixParameter(5);

  // Now ready for minimization step
  arglist[0] = 10000;
  arglist[1] = 0.1;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  //    Double_t *p1= gMinuit->Getmnparm();
  Double_t Fit_Pars[nvar] = {};
  Double_t Fit_Pars_Err[nvar] = {};
  for(Int_t i=0; i<nvar; i++)
    {
      gMinuit->GetParameter(i,Fit_Pars[i],Fit_Pars_Err[i]);
      //poly3_fit->GetParameter(i,Fit_Pars[i],Fit_Pars_Err[i]);
      cout<<"Fit_Pars["<<i<<"] = "<<Fit_Pars[i]<<endl;
    }
   
  // find chi2_min //
  Double_t  chi2_min[n7] = {0};
  Double_t  A_min[n7]= {0};
  Double_t  A1_min[n7]= {0};
  Double_t   B_min[n7]= {0};
  Double_t   P_min[n7]= {0};
  Double_t  totcount_min= 0;
  Double_t count_min;
  cout<<"********************************************"<<endl;
  for (Int_t i=0; i<n7; i++) {
    //P_min[i]= poly3(&Q[i], p1);
    P_min[i]=poly3(&Q[i],Fit_Pars);
    cout<<"P_min = "<<P_min[i]<<endl;
    A_min[i]= P_min[i]- F[i];
    A1_min[i]= pow(A_min[i],2);
    B_min[i]= pow(eF[i],2);
    chi2_min[i] = A1_min[i] /B_min[i];
    count_min= chi2_min[i];
    totcount_min= totcount_min+count_min;
    cout <<"( Q[i],F[i] ,P_min[i])" << endl;
    cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P_min[i] << ")"<< endl;
    cout<< "chi2_min=" << totcount_min << endl;
    cout<<"********************************************"<<endl;
    
  }

  for(Int_t i=0; i<nvar; i++)
    {
      cout<<"Fit_Pars["<<i<<"] = "<<Fit_Pars[i]<<endl;
    }
  
  //Normalize the data and errors in the original arrays using their floating normalization fit parameters.
  Double_t F_norm[n7] = {};
  Double_t eF_norm[n7] = {};

  //Loop over all data.
  for(Int_t i=0;i<n7;i++)
    {
      //Check which dataset the point belongs to.
      if(Dataset[i]==1)
	{
	  //Multiply the form factor and its uncertainty by the fitted floating normalization constant.
	  F_norm[i] = F[i]/Fit_Pars[4];
	  eF_norm[i] = eF[i]/Fit_Pars[4];
	}
      else if(Dataset[i]==2)
	{
	  F_norm[i] = F[i]/Fit_Pars[5];
	  eF_norm[i] = eF[i]/Fit_Pars[5];
	}
      //cout<<"Data point "<<i<<": F_norm = "<<F_norm[i]<<" eF_norm = "<<eF_norm[i]<<endl;
    }

  //Plot the overall fit and the data points with the floating normalizations applied.
  TCanvas* c_fit=new TCanvas("c_fit");

  TGraph *gr7= new TGraphErrors(n7,Q,F_norm,eQ,eF_norm); //with error
  gr7->SetTitle("3H Charge Form Factor ");
  gr7->GetXaxis()->SetTitle("Q^{2}(fm^{-2})");
  gr7->GetYaxis()->SetTitle("Charge Form Factor");
  gr7->SetMarkerStyle(8);
  gr7->SetMarkerSize(1);
  gr7->SetLineColor(1);
  gr7->SetLineWidth(1);
  gr7->GetYaxis()->SetRangeUser(0.01,200);
  gr7->GetXaxis()->SetRangeUser(0,12.5);
  gr7->GetXaxis()->SetTitleSize(0.04);
  gr7->GetYaxis()->SetTitleSize(0.04);
  gr7->GetXaxis()->CenterTitle(true);
  gr7->GetYaxis()->CenterTitle(true);
  gr7->Draw("a*");

  //TGraph *gr8= new TGraphErrors(n7,Q,F,eQ,eF);
  //gr8->SetMarkerColor(4);
  //gr8->Draw("* same");
  
  //Plot Polynomial fit.
  gStyle->SetOptFit(11111);
  TF1 * fpoly3 = new TF1("fpoly3", poly3,0,12,4);
  fpoly3->SetParameter(0,Fit_Pars[0]);
  fpoly3->SetParameter(1,Fit_Pars[1]);
  fpoly3->SetParameter(2,Fit_Pars[2]);
  fpoly3->SetParameter(3,Fit_Pars[3]);
  fpoly3->Draw("same");
  
}
  
