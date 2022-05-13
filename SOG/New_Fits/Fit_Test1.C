
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

//Int_t n7 = 6;
const int n7 = 5; //This definition lets you use the variable as the number of array elements so you don't need to hard code it in each time.
//Double_t Q[n7]  ={0,1,1.5,2,2.5,3};
//Double_t F[n7]  ={1,0.622, 0.503,0.387, 0.312, 0.267};

Double_t Q[n7]  ={1,1.5,2,2.5,3};
Double_t F[n7]  ={0.622, 0.503,0.387, 0.312, 0.267};
// create the error arrays
Double_t eQ[n7] = {0};
//Double_t eF[n7] ={0.00001, 0.014, 0.012, 0.010, 0.008, 0.007};
Double_t eF[n7] ={0.014, 0.012, 0.010, 0.008, 0.007};

const int nvar = 4; //Number of free parameters in the fit.

// fixed chi2 parameters with prof. Long
Double_t poly3(Double_t *X, Double_t *par)
{
  //When adding floating normalizations for each individual dataset you'll need an if statement to perform the polynomial fit uniquely for each dataset. Basically each data point in a dataset is multiplied by a unique constant shared by that whole dataset. This constant will also need to be raised to a power when multiplied into the poly fit I believe. You may want to pass this function a number representing the dataset. Something like:
  //if dataset == 1
  //val = par[0] + par[1]*X[0]*par[4] + par[2]*pow(X[0]*par[4],2) + par[3]*pow(X[0]*par[4],3);
  //if dataset == 2
  //val = par[0] + par[1]*X[0]*par[5] + par[2]*pow(X[0]*par[5],2) + par[3]*pow(X[0]*par[5],3);
  //This way the polynomial coefficients are shared between data sets, but each form factor value is multiplied by a floating normalization parameter.
  //You'll also need to change the number of free parameters in TF1 and Minuit to represent the new free parameters.
    Double_t val = 0.;

  val = par[0] + par[1]*X[0] + par[2]*pow(X[0],2) + par[3]*pow(X[0],3);

  return val;
}

void FCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t iflag)
// ===========new canvas===========//
// plot the old data 3He form factor by using data points
{
  /*
  TCanvas* c1;
  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  
          
// create the coordinate arrays 3He for  all data

Int_t n7 = 6;
    Double_t Q[6]  ={0,1,1.5,2,2.5,3};
    Double_t F[6]  ={1,0.622, 0.503,0.387, 0.312, 0.267};
            // create the error arrays
    Double_t eQ[5] = {0};
    Double_t eF[6] ={0.00001,0.014,
        0.012,
        0.010,
        0.008,
        0.007};

    TGraph *gr7= new TGraphErrors(n7,Q,F,eQ,eF); //with error
     gr7->SetTitle("3H Charge Form Factor ");
     gr7->GetXaxis()->SetTitle("Q^{2}(fm^{-2})");
     gr7->GetYaxis()->SetTitle("Charge Form Factor");
     gr7->SetMarkerStyle(8);
     gr7->SetMarkerSize(1);
     gr7->SetLineColor(1);
     gr7->SetLineWidth(1);
     gr7->GetYaxis()->SetRangeUser(0.01,1.01);
     gr7->GetXaxis()->SetRangeUser(0,6);
     gr7->GetXaxis()->SetTitleSize(0.04);
     gr7->GetYaxis()->SetTitleSize(0.04);
     gr7->GetXaxis()->CenterTitle(true);
     gr7->GetYaxis()->CenterTitle(true);
     gr7->Draw("a*");
    
         
          // fiting //
          
    gStyle->SetOptFit(11111);
    TF1 *poly3_fit = new TF1("poly3_fit",poly3,0,10,nvar);
    //poly3_fit->FixParameter (0,1);
    gr7->Fit(poly3_fit);
    //Access the fit resuts
    TF1 * f5 = gr7->GetFunction("f5");
   // f5->SetLineWidth(3);
   // f5->SetLineColor(2);
   // f5->Draw("same");
   // Double_t par2[4];
    Double_t *p0= poly3_fit->GetParameters();
    // find chi2 //

*/
    Double_t  chi2[n7] = {0};
    Double_t  A[n7]= {0};
    Double_t  A1[n7]= {0};
    Double_t   B[n7]= {0};
    Double_t   P[n7]= {0};
    Double_t  totcount= 0;
    Double_t count;
    for (Int_t i=0; i<n7; i++) {
             P[i]= poly3(&Q[i], par);
             A[i]= P[i]- F[i];
             A1[i]= pow(A[i],2);
             B[i]= pow(eF[i],2);
             chi2[i] = A1[i] /B[i];
             count= chi2[i];
             totcount= totcount+count;
             //cout <<"( Q[i],F[i] ,P[i])" << endl;
             //cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P[i] << ")"<< endl;
             //cout<< "chi2=" << totcount << endl;
         
    }
    
    /*
    TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
    leg->AddEntry(gr7,"All Data with Current Data");
    leg->Draw("same");
    */
    
    f = count;
}
   
void Fit_Test1()
{
//  TMinuit //
    TMinuit *gMinuit = new TMinuit(nvar);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(FCN);
    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    // Set starting values and step sizes for parameters
    //static Double_t vstart[nvar] = {1, -4.21077e-01 , 0.0713488  , -0.00447744};
    static Double_t vstart[nvar] = {1, -0.4, 0.07, -0.005};
    //static Double_t vstart[nvar] = {1, 1, 1, 1};
    static Double_t step[nvar] = {0.001 , 1.0 , 1.0 , 1.0};
    gMinuit->mnparm(0, "a1", vstart[0], step[0], -1,1,ierflg);
    gMinuit->mnparm(1, "a2", vstart[1], step[0], -1,1,ierflg);
    gMinuit->mnparm(2, "a3", vstart[2], step[0], -1,1,ierflg);
    gMinuit->mnparm(3, "a4", vstart[3], step[0], -1,1,ierflg);
    
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
}
  
