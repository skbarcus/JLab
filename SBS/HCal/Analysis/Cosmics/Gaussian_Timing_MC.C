#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"
#include <vector>
#include <TStopwatch.h>
#include <TMath.h>
#include <TF1.h>
#include<TStyle.h>
#include<TList.h>
#include<TChain.h>
#include<TRandom.h>
using namespace std;

const Int_t ngaus = 3;
Int_t nevt = 10000000;
Int_t bins = 1000;
Int_t time_min = -10;
Int_t time_max = 10;
Int_t sig1 = 2, sig2 = 4, sig3 = 3;

//Create Gaussian to fit the timing resolution.
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

void Gaussian_Timing_MC()
{

  TH1F **hGaus = new TH1F*[ngaus];
  for(Int_t i = 0; i<ngaus; i++)
    {
      hGaus[i] = new TH1F(Form("hGaus%d",i),Form("Gaussian %d",i),bins,-10.,10.);
    }

  TH1F *hTime = new TH1F("hTime","hTime",bins,time_min,time_max);
  TH1F *hRef_Time = new TH1F("hRef_Time","hRef_Time",bins,time_min,time_max);

  for(Int_t i=0;i<nevt;i++)
    {
      //hGaus[0]->Fill(gRandom->Gaus(0.,sig1));
      //hGaus[1]->Fill(gRandom->Gaus(0.,sig2));
      //hGaus[2]->Fill(gRandom->Gaus(0.,sig3));
      hTime->Fill(gRandom->Gaus(0.,sig2) - (gRandom->Gaus(0.,sig1)+gRandom->Gaus(0.,sig3))/2.);
      hRef_Time->Fill((gRandom->Gaus(0.,sig1)+gRandom->Gaus(0.,sig3))/2.);
    }

  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  hTime->Draw();
  hTime->SetLineColor(1);
  hRef_Time->Draw("same");
  //hRef_Time->SetLineColor(2);
  /*
  hGaus[0]->Draw("same");
  hGaus[0]->SetLineColor(3);
  hGaus[1]->Draw("same");
  hGaus[1]->SetLineColor(4);
  hGaus[2]->Draw("same");
  hGaus[2]->SetLineColor(6);
  */

  TF1 *func_gaus = new TF1("func_gaus",fit_gaus, time_min, time_max, 3);
  func_gaus->SetParameter(0,hTime->GetMaximum());
  func_gaus->SetParameter(1,hTime->GetMean());
  func_gaus->SetParameter(2,hTime->GetStdDev());
  hTime->Fit("func_gaus","R M");

  TF1 *func_gaus_ref = new TF1("func_gaus_ref",fit_gaus, time_min, time_max, 3);
  func_gaus_ref->SetLineColor(4);
  func_gaus_ref->SetParameter(0,hRef_Time->GetMaximum());
  func_gaus_ref->SetParameter(1,hRef_Time->GetMean());
  func_gaus_ref->SetParameter(2,hRef_Time->GetStdDev());
  hRef_Time->Fit("func_gaus_ref","R M");
}
