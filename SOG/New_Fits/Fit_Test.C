
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
  TCanvas* c1;
  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  
          
// create the coordinate arrays 3He for  all data

 Int_t n7 = 6;
 Float_t Q[6]  ={0,1,1.5,2,2.5,3};
 Float_t F[6]  ={1,0.622, 0.503,0.387, 0.312, 0.267};
            // create the error arrays
    Float_t eQ[5] = {0};
    Float_t eF[6] ={0.00001,0.014,
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
    //TF1 * f4 = new TF1("pol3", "[0]+ [1]*x+ [2]*x^2 + [3]*x^3 ");
    TF1 *poly3_fit = new TF1("poly3",poly3,0,10,4);
    poly3_fit->FixParameter (0,1);
    gr7->Fit("pol3");
    //Access the fit resuts
    TF1 *f5 = gr7->GetFunction("pol3");
    f5->SetLineWidth(3);
    f5->SetLineColor(2);
    f5->Draw("same");
 

    // find chi2 //
    Float_t  chi2[5] = {0};
    Float_t  A[5]= {0};
    Float_t  A1[5]= {0};
    Float_t   B[5]= {0};
    Float_t   P[5]= {0};
    for (Int_t i=0; i<n7; i++) {
        Double_t  totcount, count, P[i] ;
        P[i]= 1 - (4.21077e-01 * Q[i])+ (0.0713488 * Q[i]* Q[i])-(0.00447744* Q[i]* Q[i]* Q[i]);
        
             A[i]= P[i]- F[i];
             A1[i]= pow(A[i],2);
             B[i]= pow(eF[i],2);
             chi2[i] = A1[i] /B[i];
             count= chi2[i];
             totcount= totcount+count;
             cout <<"( Q[i],F[i] ,P[i])" << endl;
             cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P[i] << ")"<< endl;
             cout<< "chi2=" << totcount << endl;
         
    }
    
    
    TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
    leg->AddEntry(gr7,"All Data with Current Data");
    leg->Draw("same");
   
}
   
void Fit_Test()
{
//  TMinuit //
    TMinuit *gMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(FCN);
    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
// Set starting values and step sizes for parameters
   static Double_t vstart[4] = {1, -4.21077e-01 , 0.0713488  , -0.00447744};
   static Double_t step[4] = {1.0 , 1.0 , 1.0 , 1.0};
   gMinuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
   gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);
   gMinuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);

// Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

// Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);

}
  
