
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

  Int_t  Z =8; // number of points in each set
  Int_t  Z1 =8;
  const Int_t npts1 = 8;
  const Int_t npts2 = 8;
  Int_t n7[2] ={npts1,npts2}; //2 number of data sets.
  Double_t Q[npts1]  = {0,1,2,3,4,5,6,7};
  Double_t F[npts1] ={3,12,45,120,255,468, 777,1200};

  Double_t eQ[npts1] = {0};
  Double_t eF[npts1] = {0.046,0.048,0.028,0.023,0.019,0.014,0.011,0.009};
     
  Double_t Q1[npts2]  = {0,1,2,3,4,5,6,7};
  Double_t F1[npts2] ={0.5,2,7.5 ,20, 42.5 ,78 ,129.5,200};
  // create the ,error arrays
  Double_t eQ1[npts2] = {0};
  Double_t eF1[npts2] = {0.046,0.048,0.028,0.023,0.019,0.014,0.011,0.009};

  const Int_t n9 = npts1+npts2;
  Double_t Qtot[n9]= {0};
  Double_t Ftot[n9]={0};
  Double_t eQtot[n9]= {0};
  Double_t eFtot[n9]= {0};
  Int_t dataset[n9] = {0};
  const int nvar = 6; //Number of free parameters in the fit.

// fixed chi2 parameters with prof. Long
Double_t poly3(Double_t *X, Double_t *par)
{
  //When adding floating normalizations for each individual dataset you'll need an if statement to perform the polynomial fit uniquely for each dataset. so what we did we add par[4] to the y-axis which is form factor to move the data point up and down. Then we want the fitting eqaution pass throug 1 so we need to divided all the equation by par[0]*par[4]
  Double_t val = 0.;

  // val = 1 + par[0]*par[4]*X[0] + par[1]*par[4]*pow(X[0],2) + par[2]*par[4]*pow(X[0],3);
  val = (par[0]+ par[1]*X[0] + par[2]*pow(X[0],2) + par[3]* pow(X[0],3));

  return val;
}

void FCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t iflag)
// ===========new canvas===========//
// plot the old data 3He form factor by using data points
{
  Double_t  totcount= 0;
  Double_t  chi22;
  Double_t  chi2[n9] = {0};
  Double_t  A[n9]= {0};
  Double_t  A1[n9]= {0};
  Double_t   B[n9]= {0};
  Double_t   P[n9]= {0};
  Double_t  count= 0;
  Double_t  N1= par[4];
  Double_t  N2= par[5];

  for (Int_t i=0; i<n9; i++) {

    if(dataset[i]==1)
      {
	P[i]= (poly3(&Qtot[i], par)/N1);
      }
    else if(dataset[i]==2)
      {
	P[i]= (poly3(&Qtot[i], par)/N2);
      }
      A[i]= P[i]- Ftot[i];
      A1[i]= pow(A[i],2);
      B[i]= pow(eFtot[i],2);
      chi2[i] = A1[i] /B[i];
      count= chi2[i];
      totcount= totcount+count;
             
      //  cout <<"( Q[i],F[i] ,P[i])" << endl;
      //  cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P[i] << ")"<< endl;
      //  cout<< "chi2=" << totcount << endl;
      //chi22= pow(((N-1)/0.1),2);
    if(dataset[i]==1)
      {
        chi22= pow(((N1-1)/0.1),2);
      }
    else if(dataset[i]==2)
      {
        chi22= pow(((N2-1)/0.1),2);
      }
    totcount= totcount+chi22;
  }
    
  //  cout <<"( Q[i],F[i] ,P[i])" << endl;
  //  cout <<"( "<< Q[i]<< ","<< F[i]<< ","<< P[i] << ")"<< endl;
  //cout<< "chi2=" << totcount << endl;
  f = totcount;
}

void Fit_lk1()
{

  for (Int_t i=0; i<n7[0]; i++) {
    //cout<<i<<": Q[i] = "<<Q[i]<<"   F[i] = "<<F[i]<<"   eQ[i] = "<<eQ[i]<<"   eF[i] = "<<eF[i]<<endl;
    //cout<<i<<": Q1[i] = "<<Q1[i]<<"   F1[i] = "<<F1[i]<<"   eQ1[i] = "<<eQ1[i]<<"   eF1[i] = "<<eF1[i]<<endl;
    Qtot[i]= Q[i];
    Ftot[i]= F[i];
    eQtot[i]= eQ[i];
    eFtot[i]= eF[i];
    dataset[i] = 1;
    //cout<<i<<endl;
  }
  for (Int_t i=0; i<n7[1]; i++) {
    //cout<<i<<": Q[i] = "<<Q[i]<<"   F[i] = "<<F[i]<<"   eQ[i] = "<<eQ[i]<<"   eF[i] = "<<eF[i]<<endl;
    //cout<<i<<": Q1[i] = "<<Q1[i]<<"   F1[i] = "<<F1[i]<<"   eQ1[i] = "<<eQ1[i]<<"   eF1[i] = "<<eF1[i]<<endl;
    //cout<<F1[i]<<endl;
    Qtot[i+n7[0]]= Q1[i];
    Ftot[i+n7[0]]= F1[i];
    eQtot[i+n7[0]]= eQ1[i];
    eFtot[i+n7[0]]= eF1[i];
    dataset[i+n7[0]] = 2;
    //cout<<i<<endl;
  }
    
  for (Int_t i=0; i<n9; i++) 
    {
      cout<<"Qtot["<<i<<"] = "<<Qtot[i]<<endl;
    }

  //definsion();
  //  TMinuit //
  TMinuit *gMinuit = new TMinuit(nvar);  //initialize TMinuit with a maximum of 4 params
  gMinuit->SetFCN(FCN);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  // static Double_t vstart[nvar] = { -0.421077 , 0.0713488  , -0.00447744};
  static Double_t vstart[nvar] = { 1, 1 ,1,1,3,0.5};
  static Double_t step[nvar] = {0.01,1,1,1,1,1};
  gMinuit->mnparm(0, "a0", vstart[0], step[0], 0.8,1.2,ierflg);
  gMinuit->mnparm(1, "a1", vstart[1], step[0], 0.8,1.2,ierflg);
  gMinuit->mnparm(2, "a2", vstart[2], step[0], 0.8,1.2,ierflg);
  gMinuit->mnparm(3, "a3", vstart[3], step[0], 0.8,1.2,ierflg);
  gMinuit->mnparm(4, "a4", vstart[4], step[0], 2.5,3.5,ierflg);
  gMinuit->mnparm(5, "a5", vstart[5], step[0], 0.1,1,ierflg);
  
  //gMinuit->FixParameter(0);
  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);
  //gMinuit->FixParameter(4);
  //gMinuit->FixParameter(5);
    
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 0.01;
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
      //cout<<"vstart["<<i<<"] = "<<vstart[i]<<endl;
    }
  /*
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
  */

  TCanvas* c1;
  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

        
  // create the coordinate arrays 3He for  all data

  TGraph *gr7= new TGraphErrors(n9,Qtot,Ftot,eQtot,eFtot); //with error
  gr7->SetTitle("3H Charge Form Factor ");
  gr7->GetXaxis()->SetTitle("Q^{2}(fm^{-2})");
  gr7->GetYaxis()->SetTitle("Charge Form Factor");
  gr7->SetMarkerStyle(8);
  gr7->SetMarkerSize(1);
  gr7->SetLineColor(1);
  gr7->SetLineWidth(1);
  gr7->GetYaxis()->SetRangeUser(0,2000);
  gr7->GetXaxis()->SetRangeUser(0,10);
  gr7->GetXaxis()->SetTitleSize(0.04);
  gr7->GetYaxis()->SetTitleSize(0.04);
  gr7->GetXaxis()->CenterTitle(true);
  gr7->GetYaxis()->CenterTitle(true);
  gr7->Draw("a*");
    
  // fiting //
  gStyle->SetOptFit(11111);
  TF1 * poll3 = new TF1("poll3", "[0]+ [1]*x + [2]* x^2 + [3]*x^3 ");
  //  gr7->Fit("poll3");
  // poll3->FixParameter (4,1);
  //Access the fit resuts
  // TF1 *f5 = gr7->GetFunction("poll3");
  // f5->SetLineWidth(3);
  // f5->SetLineColor(2);
  //  f5->Draw("same");
    
  
  TF1 * fpoly3 = new TF1("fpoly3", poly3,0,6,4);
  fpoly3->SetParameter(0,Fit_Pars[0]);
  fpoly3->SetParameter(1,Fit_Pars[1]);
  fpoly3->SetParameter(2,Fit_Pars[2]);
  fpoly3->SetParameter(3,Fit_Pars[3]);
  //fpoly3->SetParameter(4,Fit_Pars[4]);
  //fpoly3->SetLineColor(3);
  
  fpoly3->Draw("same");
    
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
