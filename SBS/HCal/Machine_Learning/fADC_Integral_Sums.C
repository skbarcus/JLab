#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include <vector>
#include "hcal.h"
#include <TStopwatch.h>
#include <TMath.h>
#include <TF1.h>
#include<TStyle.h>
#include<TList.h>
#include<TChain.h>
using namespace std;

Double_t entry;
Int_t gCurrentEntry = 0;
Int_t DISP_MIN_SAMPLE = 0;
Int_t DISP_MAX_SAMPLE = 20;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);

Int_t limit_evts = 0;
Int_t loop_max = 0;
Int_t max_evts = 10000;

const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;

Double_t f1_res = 0.112;//0.112
Double_t fadc_res = 4.;

Int_t nevt;
Int_t r,c,n,idx;
TChain *T = 0;
Double_t fadc_int_sum = 0;

Double_t adc[kNrows][kNcols];
Double_t tdc[kNrows][kNcols];

TH1F *histos[kNrows][kNcols];

Int_t max_fadc_int_sum = 0;
Int_t max_fadc_int_sum_temp = 0;
Int_t max_fadc_int_sum_evt;

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

void fADC_Integral_Sums(Int_t run = 820)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);
  
  if(!T) 
    { 
      T = new TChain("T");
      //T->Add(TString::Format("rootfiles/fadc_f1tdc_%d.root",run));//325
      //T->Add(TString::Format("rootfiles/fadc_f1tdc_%d_1.root",run));
      
      
      //TChain chain("T");
      //chain.Add(TString::Format("rootfiles/fadc_f1tdc_%d.root",run));
      //chain.Add(Form("rootfiles/fadc_f1tdc_%d.root",run));
      
      
      //============  Reading the Rootfile =======================//
      
      const TString rootfilePath = "/home/skbarcus/JLab/SBS/HCal/Analysis/Cosmics/rootfiles/";
      std::ostringstream string;
      string << rootfilePath<<"fadc_f1tdc_"<<run;
      TString basename = string.str().c_str();
      TString rootfile = basename + ".root";
      cout<<basename<<endl;
      //TChain* T;
      //T = new TChain("T");
      
      Long_t split=0;
      char* file = 0;
      
      //====adding splits rootfiles =======================//
      
      Long_t u=0;
      while ( !gSystem->AccessPathName(rootfile.Data()) ) 
	{
	  T->Add(rootfile.Data());
	  cout << "ROOT file " << rootfile << " added to TChain." << endl;
	  u++;
	  rootfile = basename + "__" + u + ".root";
	}
      
      if(!T->GetEntries())
	{
	  cerr<< "No root file was found" << endl;
	  return;
	}
      //==finish adding splits rootfiles=====================//
      
      //T->Add(TString::Format("rootfiles/fadc_f1tdc_%d.root",run));
      T->SetBranchStatus("*",0);
      T->SetBranchStatus("sbs.hcal.*",1);
      T->SetBranchAddress("sbs.hcal.nsamps",hcalt::nsamps);       //Number of samples for given row-col
      T->SetBranchAddress("sbs.hcal.a",hcalt::a);                 //Raw ADC amplitudes
      T->SetBranchAddress("sbs.hcal.tdc",hcalt::tdc);             //Raw TDC value
      //T->SetBranchAddress("sbs.hcal.ledbit",&hcalt::ledbit);
      //T->SetBranchAddress("sbs.hcal.ledcount",&hcalt::ledcount);
      T->SetBranchAddress("sbs.hcal.samps",hcalt::samps);         //RAW ADC samples
      T->SetBranchAddress("sbs.hcal.samps_idx",hcalt::samps_idx); //Index in samples vector for given row-col module 
      T->SetBranchAddress("sbs.hcal.row",hcalt::row);             //Row for block in data vectors
      T->SetBranchAddress("sbs.hcal.col",hcalt::col);             //Col for block in data vectors
      T->SetBranchStatus("Ndata.sbs.hcal.row",1);
      T->SetBranchAddress("Ndata.sbs.hcal.row",&hcalt::ndata);    //???
      std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    }
  
  nevt = T->GetEntries();
  cout<<"Number of events = "<<nevt<<endl;

  //Make histos for fADC samples.
  for(Int_t r = 0; r < kNrows; r++) 
    {
      for(Int_t c = 0; c < kNcols; c++) 
	{
	  histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
	}
    }

  TH1F *hfadc_int_sums = new TH1F("hfadc_int_sums","fADC Integral Sums per Event",1000,7e5,8.5e6);

  //Loop over all events.
  if(limit_evts==1)
    {
      loop_max = max_evts;
    }
  else
    {
      loop_max = nevt;
    }

  for(Int_t i=0; i<loop_max ;i++)
    {
      T->GetEntry(gCurrentEntry);

      //Set arrays to zero;
      for(r  = 0; r < kNrows; r++) 
	{
	  for(c  = 0; c < kNcols; c++) 
	    {
	      adc[r][c] = 0.0;
	      tdc[r][c] = 0.0;
	    }
	}

      Double_t fadc_int[channels] = {};
      //Loop over all channels and fill the arrays. Do this seperately so reference time is set for all channels.
      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];
	  adc[r][c] = hcalt::a[j];
	  tdc[r][c] = hcalt::tdc[j];
	  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) //Takes about 1 minute for 300k evts.
	    {
	      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	      fadc_int[j] = fadc_int[j] + histos[r][c]->GetBinContent(s+1-DISP_MIN_SAMPLE);
	    }
	  //No pedestal subtraction.
	  fadc_int[j] = fadc_int[j];
	  //cout<<"fADC integral (no pedestal subtraction) = "<<fadc_int[j]<<endl;
	  
	  //Sum fADC integrals for histo.
	  fadc_int_sum = fadc_int_sum + adc[r][c];
	}
      hfadc_int_sums->Fill(fadc_int_sum);
      //cout<<"Event "<<i<<": fadc_int_sum = "<<fadc_int_sum<<endl;
      max_fadc_int_sum_temp = fadc_int_sum;
      if(max_fadc_int_sum_temp>max_fadc_int_sum)
	{
	  max_fadc_int_sum = max_fadc_int_sum_temp;
	  max_fadc_int_sum_evt = i;
	}
      fadc_int_sum = 0;
      gCurrentEntry++;
    }

  TCanvas *cfadc_int_sums = new TCanvas("cfadc_int_sums");
  cfadc_int_sums->cd(1);
  hfadc_int_sums->GetXaxis()->SetTitle("fADC Integral Sum per Event");
  hfadc_int_sums->GetXaxis()->CenterTitle(true);
  hfadc_int_sums->GetXaxis()->SetLabelSize(0.05);
  hfadc_int_sums->GetXaxis()->SetTitleSize(0.06);
  hfadc_int_sums->GetXaxis()->SetTitleOffset(0.75);
  hfadc_int_sums->GetYaxis()->SetTitle("Counts");
  hfadc_int_sums->GetYaxis()->CenterTitle(true);
  hfadc_int_sums->GetYaxis()->SetLabelSize(0.05);
  hfadc_int_sums->GetYaxis()->SetTitleSize(0.06);
  hfadc_int_sums->GetYaxis()->SetTitleOffset(0.85);
  gPad->SetLogy();
  hfadc_int_sums->Draw("");

  cout<<"Max fADC Integral Sum Event = "<<max_fadc_int_sum_evt<<" with a value of "<<max_fadc_int_sum<<"."<<endl;

  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
