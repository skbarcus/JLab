
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

  
  const Int_t npts1 = 6; // number of points in each set
  const Int_t npts2 = 6;
  Int_t n7[2] ={npts1,npts2}; //2 number of data sets.

/*
  Double_t Q[npts1]  = {1,2,3,4,5,6,7,8};
  Double_t F[npts1] ={1.15,4.6,17.25,46,97.75,179.4,297.85,460};
  
  Double_t UN = 0.01;
  Double_t eQ[npts1] = {0};
  Double_t eF[npts1] = {Q[0]*UN,Q[1]*UN,Q[2]*UN,Q[3]*UN,Q[4]*UN,Q[5]*UN,Q[6]*UN,Q[7]*UN};

  Double_t Q1[npts2]  = {1,2,3,4,5,6,7,8};
  Double_t F1[npts2]= {0.95,3.8,14.25,38,80.75,148.2,246.05,380};
  // create the ,error arrays
  Double_t UN1 = 0.01;
  Double_t eQ1[npts2] = {0};
  Double_t eF1[npts2] = {Q[0]*UN1,Q[1]*UN1,Q[2]*UN1,Q[3]*UN1,Q[4]*UN1,Q[5]*UN1,Q[6]*UN1,Q[7]*UN1};
*/

Double_t dA = 0.2;
Double_t dB = -0.1;
Double_t NormA = 1+dA;
Double_t NormB = 1+dB;

Double_t Q[npts1]  = {1,3,5,7,9,11};
Double_t F[npts1] ={pow(Q[0],2)*NormA,pow(Q[1],2)*NormA,pow(Q[2],2)*NormA,pow(Q[3],2)*NormA,pow(Q[4],2)*NormA,pow(Q[5],2)*NormA};

Double_t eQ[npts1] = {0};
Double_t eF[npts1] = {pow(Q[0],2)*abs(dA*0.3),pow(Q[1],2)*abs(dA*0.3),pow(Q[2],2)*abs(dA*0.3),pow(Q[3],2)*abs(dA*0.3),pow(Q[4],2)*abs(dA*0.3),pow(Q[5],2)*abs(dA*0.3)};

Double_t Q1[npts2]  = {2,4,6,8,10,12};
Double_t F1[npts2] ={pow(Q1[0],2)*NormB,pow(Q1[1],2)*NormB,pow(Q1[2],2)*NormB,pow(Q1[3],2)*NormB,pow(Q1[4],2)*NormB,pow(Q1[5],2)*NormB};
// create the ,error arrays
Double_t eQ1[npts2] = {0};
Double_t eF1[npts2] = {pow(Q1[0],2)*abs(dB*0.3),pow(Q1[1],2)*abs(dB*0.3),pow(Q1[2],2)*abs(dB*0.3),pow(Q1[3],2)*abs(dB*0.3),pow(Q1[4],2)*abs(dB*0.3),pow(Q1[5],2)*abs(dB*0.3)};



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

  for (Int_t i=0; i<n9; i++) {

    if(dataset[i]==1)
      {
	P[i]= (poly3(&Qtot[i], par)*par[4]);
      }
   else if(dataset[i]==2)
    //  {
	P[i]= (poly3(&Qtot[i], par)*par[5]);
   //   }
      A[i]= P[i]- Ftot[i];
      A1[i]= pow(A[i],2);
      B[i]= pow(eFtot[i],2);
      chi2[i] = A1[i] /B[i];
      count= chi2[i];
      totcount= totcount+count;
             
    if(dataset[i]==1)
    {
      chi22= pow(((par[4]-1)/0.1),2);
      //cout <<"(cc"<< Qtot[i]<< ","<< Ftot[i]<< ","<< P[i] << ")"<< endl;
    }
      
    else if(dataset[i]==2)
      {
       chi22= pow(((par[5]-1)/0.1),2);
       //cout <<"(cc"<< Qtot[i]<< ","<< Ftot[i]<< ","<< P[i] << ")"<< endl;
      }
      totcount= totcount+chi22;
      
  }
  //cout<< "chi2=" << totcount << endl;
  //cout<< "chi2/DF=" << totcount/ n9 <<endl;
  //cout<< "Par[4]=" << par[4] << "  Par[5]= " << par[5] << endl;
  f = totcount;
}

void Fit_lk3()
{
    cout<<"n7[0] "<<n7[0]<< endl;
    cout<<"n7[1] "<<n7[1]<< endl;
  for (Int_t i=0; i<n7[0]; i++) {

    Qtot[i]= Q[i];
    Ftot[i]= F[i];
    eQtot[i]= eQ[i];
    eFtot[i]= eF[i];
    dataset[i] = 1;
    //cout<<i<<endl;
  }
    
  for (Int_t i=0; i<n7[1]; i++) {
    Qtot[i+n7[0]]= Q1[i];
    Ftot[i+n7[0]]= F1[i];
    eQtot[i+n7[0]]= eQ1[i];
    eFtot[i+n7[0]]= eF1[i];
    dataset[i+n7[0]] = 2;
  }
    
  for (Int_t i=0; i<n9; i++) 
    {
      cout<<"Qtot["<<i<<"] = "<<Qtot[i]<<endl;
    }

  //  TMinuit //
  TMinuit *gMinuit = new TMinuit(nvar);  //initialize TMinuit with a maximum of 4 params
  gMinuit->SetFCN(FCN);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  //static Double_t vstart[nvar] = { 0.9, 1.1 ,0.95,0.99,1,1};
  static Double_t vstart[nvar] = {0.0, 0.0, 1., 0.0, 1.3, 0.8};
  static Double_t step[nvar] = {0.0001 , 0.0001 , 0.0001 , 0.0001, 0.0001,0.0001};
  gMinuit->mnparm(0, "a0", vstart[0], step[0], -10,10,ierflg);
  gMinuit->mnparm(1, "a1", vstart[1], step[1], -10,10,ierflg);
  gMinuit->mnparm(2, "a2", vstart[2], step[2], -10,10,ierflg);
  gMinuit->mnparm(3, "a3", vstart[3], step[3], -10,10,ierflg);
  gMinuit->mnparm(4, "a4", vstart[4], step[4], -10,10,ierflg);
  gMinuit->mnparm(5, "a5", vstart[5], step[5], -10,10,ierflg);
  
  //gMinuit->FixParameter(0);
  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  //gMinuit->FixParameter(3);
  //gMinuit->FixParameter(4);
  //gMinuit->FixParameter(5);
    
  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 0.0001;
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
      cout<<"Fit_Pars["<<i<<"] = "<<Fit_Pars[i]<<endl;
      
    }
  

  TCanvas* c1;
  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
    Double_t Ff[npts1] ={F[0]/Fit_Pars[4],F[1]/Fit_Pars[4],F[2]/Fit_Pars[4],F[3]/Fit_Pars[4],F[4]/Fit_Pars[4],F[5]/Fit_Pars[4]};
    Double_t eFf[npts1] ={eF[0]/Fit_Pars[4],eF[1]/Fit_Pars[4],eF[2]/Fit_Pars[4],eF[3]/Fit_Pars[4],eF[4]/Fit_Pars[5],eF[5]/Fit_Pars[4]};
  // create the coordinate arrays 3H for  all data
  TGraph *gr7= new TGraphErrors(npts1,Q,F,eQ,eF); //with error
  gr7->SetMarkerStyle(33);
  gr7->SetMarkerColor(2);
  gr7->SetMarkerSize(2);
  gr7->SetLineColor(1);
  gr7->SetLineWidth(1);
 
    Double_t Ff1[npts2] ={F1[0]/Fit_Pars[5],F1[1]/Fit_Pars[5],F1[2]/Fit_Pars[5],F1[3]/Fit_Pars[5],F1[4]/Fit_Pars[5],F1[5]/Fit_Pars[5]};
    Double_t eFf1[npts2] ={eF1[0]/Fit_Pars[5],eF1[1]/Fit_Pars[5],eF1[2]/Fit_Pars[5],eF1[3]/Fit_Pars[5],eF1[4]/Fit_Pars[5],eF1[5]/Fit_Pars[5]};
    
    TGraph *gr8= new TGraphErrors(npts2,Q1,F1,eQ1,eF1); //with error
    gr8->SetMarkerStyle(33);
    gr8->SetMarkerSize(1.5);
    gr8->SetLineColor(1);
    gr8->SetLineWidth(1);
    auto g = new TMultiGraph();
         g->Add(gr7);
         g->Add(gr8);
         g->Draw("AP");
    gPad->Modified();
    g->GetYaxis()->SetRangeUser(0,200);
    g->GetXaxis()->SetRangeUser(0,12.5);
    g->SetTitle("3H Charge Form Factor ");
    g->GetXaxis()->SetTitle("Q^{2}(fm^{-2})");
    g->GetYaxis()->SetTitle("Charge Form Factor");
    g->GetXaxis()->SetTitleSize(0.04);
    g->GetYaxis()->SetTitleSize(0.04);
    g->GetXaxis()->CenterTitle(true);
    g->GetYaxis()->CenterTitle(true);

  
  TF1 * fpoly3 = new TF1("fpoly3", poly3,0,12,4);
  fpoly3->SetParameter(0,Fit_Pars[0]);
  fpoly3->SetParameter(1,Fit_Pars[1]);
  fpoly3->SetParameter(2,Fit_Pars[2]);
  fpoly3->SetParameter(3,Fit_Pars[3]);
//fpoly3->SetLineColor(3);
  fpoly3->Draw("same");
    
      TLegend *leg=new TLegend(0.2,0.2,0.5,0.5);
      leg->AddEntry(gr7,"Data set1");
      leg->AddEntry(gr8,"Data set2");
      leg->AddEntry(fpoly3,"Fitting");
      leg->Draw("same");

}
