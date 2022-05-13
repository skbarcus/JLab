
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


const int n7 = 33; //This definition lets you use the variable as the number of array elements so you don't need to hard code it in each time.
Double_t Q[n7]  ={1,1.5,2,2.5,3, 3.5, 4,4.5, 5, 6,0.290,0.386,0.450,0.481,0.600,0.670,0.775,0.850
    ,0.900,0.951,1.000,0.051,0.098,0.258,0.571,1.05,1.51,2.21,2.98, 0.3,0.5,0.7,0.9};
Double_t F[n7]  ={0.622, 0.503,0.387, 0.312, 0.267, 0.215, 0.175, 0.137, 0.118, 0.0758,0.885,0.863,0.769,0.784,0.79,0.789,0.769,0.743,0.603,0.682,0.579,0.967,1.039
    ,0.95,0.826,0.697,0.548,0.387,0.282, 0.872,0.782,0.785,0.659};
// create the error arrays
Double_t eQ[n7] = {0};
Double_t eF[n7] ={0.014,0.012,0.010,0.008,0.007,0.006,0.005,0.004,0.005,0.004,0.062
,0.053,0.082,0.071,0.064,0.0234,0.088,0.044,0.071,0.054,0.075, 0.047,0.046,0.029,0.028,0.027
,0.025,0.028,0.032,0.047,0.046,0.038,0.035};



const int nvar = 5; //Number of free parameters in the fit.

// fixed chi2 parameters with prof. Long
Double_t poly3(Double_t *X, Double_t *par)
{
  //When adding floating normalizations for each individual dataset you'll need an if statement to perform the polynomial fit uniquely for each dataset. so what we did we add par[4] to the y-axis which is form factor to move the data point up and down. Then we want the fitting eqaution pass throug 1 so we need to divided all the equation by par[0]*par[4]
    Double_t val = 0.;

 // val = 1 + par[0]*par[4]*X[0] + par[1]*par[4]*pow(X[0],2) + par[2]*par[4]*pow(X[0],3);
    val = (par[0]+ par[1]*X[0] + par[2]*pow(X[0],2)+ par[3]* pow(X[0],3))*par[4];

  return val;
}

void FCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t iflag)
// ===========new canvas===========//
// plot the old data 3He form factor by using data points
{
  
    Double_t  chi2[n7] = {0};
    Double_t  A[n7]= {0};
    Double_t  A1[n7]= {0};
    Double_t   B[n7]= {0};
    Double_t   P[n7]= {0};
    Double_t  totcount= 0;
    Double_t count= 0;
    for (Int_t i=0; i<n7; i++) {
             P[i]= poly3(&Q[i], par);
             A[i]= P[i]- F[i];
             A1[i]= pow(A[i],2);
             B[i]= pow(eF[i],2);
             chi2[i] = A1[i] /B[i];
             count= chi2[i];
             totcount= totcount+count;
           //  cout <<"( Q[i],F[i] ,P[i])" << endl;
           //  cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P[i] << ")"<< endl;
           //  cout<< "chi2=" << totcount << endl;
         
    }
    
  //  cout <<"( Q[i],F[i] ,P[i])" << endl;
  //  cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P[i] << ")"<< endl;
    cout<< "chi2=" << totcount << endl;
    f = totcount;
}
   
void Fit_Test4()
{
//  TMinuit //
    TMinuit *gMinuit = new TMinuit(nvar);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(FCN);
    Double_t arglist[10];
    Int_t ierflg = 0;
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    // Set starting values and step sizes for parameters
   // static Double_t vstart[nvar] = { -0.421077 , 0.0713488  , -0.00447744};
    static Double_t vstart[nvar] = {0.9, 0.9, 0.9,0.9,0.9};
    static Double_t step[nvar] = {0.001 , 1.0 , 1.0,1.0,1.0};
    gMinuit->mnparm(0, "a0", vstart[0], step[0], -10,10,ierflg);
    gMinuit->mnparm(1, "a1", vstart[1], step[0], -1,1,ierflg);
    gMinuit->mnparm(2, "a2", vstart[2], step[0], -1,1,ierflg);
    gMinuit->mnparm(3, "a3", vstart[3], step[0], -1,1,ierflg);
    gMinuit->mnparm(4, "a4", vstart[4], step[0], -10,10,ierflg);
    
    // Now ready for minimization step
    arglist[0] = 10000;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    Double_t Fit_Pars[nvar] = {};
    Double_t Fit_Pars_Err[nvar] = {};
    for(Int_t i=0; i<nvar; i++)
      {
       gMinuit->GetParameter(i,Fit_Pars[i],Fit_Pars_Err[i]);
       //poly3_fit->GetParameter(i,Fit_Pars[i],Fit_Pars_Err[i]);
       cout<<"Fit_Pars["<<i<<"] = "<<Fit_Pars[i]<<endl;
       cout<<"vstart["<<i<<"] = "<<vstart[i]<<endl;
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
    for (Int_t j=0; j<n7; j++) {
//	P_min[i]= poly3(&Q[i], p1);
	P_min[j]=poly3(&Q[j],Fit_Pars);
	cout<<"P_min = "<<P_min[j]<<endl;
	A_min[j]= P_min[j]- F[j];
	A1_min[j]= pow(A_min[j],2);
	B_min[j]= pow(eF[j],2);
	chi2_min[j] = A1_min[j] /B_min[j];
	count_min= chi2_min[j];
	totcount_min= totcount_min+count_min;
	//cout <<"( Q[j],F[j] ,P_min[j])" << endl;
	//cout <<"( "<< Q[j]<< ","<< F[j]<< ","<< P_min[j] << ")"<< endl;
	//cout<< "chi2_min=" << totcount_min << endl;
	//cout<<"********************************************"<<endl;
        
    }
    cout<< "chi2_min=" << totcount_min << endl;
    cout<<"********************************************"<<endl;


TCanvas* c1;
c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

        
// create the coordinate arrays 3He for  all data

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
    TF1 * poll3 = new TF1("poll3", "[0]+ [1]*x + [2]* x^2 + [3]*x^3 ");
    gr7->Fit("poll3");
   // poll3->FixParameter (4,1);
    //Access the fit resuts
    TF1 *f5 = gr7->GetFunction("poll3");
    f5->SetLineWidth(3);
    f5->SetLineColor(2);
    f5->Draw("same");;
    
   /* Int_t n = 6;
    
    
   
    for (Int_t i=0; i<n; i++) {
        
        Double_t  q[i]={0} ;
        Double_t  f[i];
        q[i] = i*0.001;
        
        f[i] = 0.995  -0.421077 * q[i] + 0.0713488*pow(q[i],2)+ -0.00447744*pow(q[i],3);
    }
  
        TGraph *gr2 = new TGraph (n, q[i], f[i]);
        TMultiGraph* mg  = new TMultiGraph();
        mg->Add(gr2, "AC*");
        mg->Add(gr7, "AL*");
        mg->Draw("AC*");
TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
leg->AddEntry(gr7,"All Data with Current Data");
leg->Draw("same");*/
}
