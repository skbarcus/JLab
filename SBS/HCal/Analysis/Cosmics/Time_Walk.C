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
using namespace std;

Double_t entry;
Int_t gCurrentEntry = 0;
Int_t limit_evts = 1;                  //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t loop_max = 0;                    //Dummy variable set equal to max_evts or nevt for the loop.
Int_t max_evts = 1;                    //Maximum number of events to analyze if limit_evts = 1.
Int_t plot_time_walks =0;              //0 = don't plot individual time walk plots. 1 = do plot individual time walk plots for each PMT.
Int_t calc_corr_time_res = 0;          //0 = Just plots the time walk fadc integral vs threshold crossing time plots. 1 = also use the time walk plot fits to get a corrected timing resolution that `removes' the time walk effect. 
Int_t ref_trig_plots = 0;              //0 = plot nothing for the ref ch(s). 1 = Plot the individual histos and exp fits to the trigger for the final event analyzed. 2 = Plot the relative timing resolution between two reference channels.
Int_t ref_fit_type = 1;                //1 = exponential fir of ref ch. 2 = linear fit of ref ch.
Int_t use_cfd = 2;                     //Use a constant fraction threshold to measure the crossing time defined as 1/4 the current peak height. 2 = Use CFD technique of scaling the signal by a constant and inverting and delaying a copied signal and summing it with the scale signal to get a histo with a zero crossing for time. (See notes from CAEN). 
Double_t cfd_thr = 0.25;              //Threshold setting for constant threshold fraction operation. 
Int_t cfd_delay = 2;                //Number of fADC units (4ns each) to delay the inverted fADC signal for CFD timing.
Int_t plot_ind_landau = 0;             //0 = Don't plot the individual landau fit for mod_row and mod_col. Note this landau function fit only updates when there is a TDC signal so plotting a histo that didn't have a TDC hit will show the last event of that module that did. 1 = plot individual PMT fADC Landau fit for mod_row mod_col PMT of the final event analyzed. 
Int_t mod_row = 1;                     //Row and column for plotting a single fADC with Landau and threshold displayed. Row and columns can be entered counting from 1 as on the actual HCal.
Int_t mod_col = 1;
Int_t module = 0;
Int_t DISP_MIN_SAMPLE = 0;
Int_t DISP_MAX_SAMPLE = 20;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
Double_t fadc_int_max_cut = 500;
Double_t fadc_int_min_cut = 500;
Int_t min_time = -300;
Int_t max_time = 50;
Int_t bins = 100;
Int_t ref_ch = 143;
Int_t ref_row = 11;
Int_t ref_col = 11;
Int_t ref_ch_2 = 142;                  //Second reference channel if want to see relative timing resolutions between same triggers.
Int_t ref_row_2 = 11;
Int_t ref_col_2 = 10;
Int_t fit_width = 2;                   //Number of bins between lower bound and upper bound of exp ref ch fits.
Int_t ref_th = 1250;                   //The reference ch fADC threshold chosen to pick the reference time.
Int_t lower_fit_th = 400;              //Threshold above which first fADC bin is used as the lower bound for the exponential ref ch fit.
Int_t r,c,n,idx;
Int_t nevt;
Int_t fit_status;
const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
Double_t ind_res[144];
Double_t f1_res = 0.112;//0.112
Double_t fadc_res = 4.;
TChain *T = 0;
std::string user_input;

Int_t skip = 3;                        //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                      //Counts number of lines in the data file. 
Int_t ncols;                           //Set how many columns of data we have in the data file.
char str[1000];                        //Variable to read lines of the data file.
Float_t mod[channels], pedestal[channels], threshold[channels];
Float_t mod_temp, pedestal_temp, threshold_temp;

TH1F *histos[kNrows][kNcols];
TH1F *histos_vert[kNrows][kNcols];
TH1F *histos_cfd_k[kNrows][kNcols];
TH1F *histos_cfd_d[kNrows][kNcols];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

TH1F* MakeHisto_vert(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d_vert",row,col),
		     TString::Format("%d-%d_vert",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

TH1F* MakeHisto_cfd_k(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d_cfd_k",row,col),
		     TString::Format("%d-%d_cfd_k",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

TH1F* MakeHisto_cfd_d(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d_cfd_d",row,col),
		     TString::Format("%d-%d_cfd_d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Create Landau to fit fADC pulses.
Double_t fit_landau(Double_t *X, Double_t *par) 
{
  Double_t landau = par[0]*TMath::Landau(X[0],par[1],par[2])+par[3];
  return landau;
}

//Create function to fit time walk correction.
Double_t fit_log(Double_t *X, Double_t *par) 
{
  Double_t log = par[0] * TMath::Log(par[1] * (X[0] - par[2])) + par[3];
  return log;
}

//Create function to fit time walk correction.
Double_t fit_poly3(Double_t *X, Double_t *par) 
{
  Double_t poly3 = par[0] + par[1]*X[0] + par[2]*pow(X[0],2.) + par[3]*pow(X[0],3.);
  return poly3;
}

Double_t fit_poly4(Double_t *X, Double_t *par) 
{
  Double_t poly4 = par[0] + par[1]*X[0] + par[2]*pow(X[0],2.) + par[3]*pow(X[0],3.) + par[4]*pow(X[0],4.); //+ par[5]*pow(X[0],5.) + par[6]*pow(X[0],6.) + par[7]*pow(X[0],7.) + par[8]*pow(X[0],8.) + par[9]*pow(X[0],9.) + par[10]*pow(X[0],10.);
  return poly4;
}

//Create a linear fit.
Double_t fit_lin(Double_t *X,Double_t *par) 
{
  Double_t lin = par[0]*X[0]+par[1];
  return lin;
}

Double_t fit_exp(Double_t *X, Double_t *par) 
{
  //Double_t exp = par[0] * TMath::Exp(par[1]*X[0]) + par[4] + par[2] * TMath::Exp(par[3]*X[0]);
  //Double_t exp = par[0] * TMath::Exp(par[1]*(X[0]-par[2])) + par[3];
  Double_t exp = par[0]*TMath::Exp(par[1]+par[2]*X[0]) + par[3];
  return exp;
}

//Exponential.
/*
TF1 **func_exp = new TF1*[channels];
for(Int_t i = 0; i<channels; i++)
  {
    func_exp[i] = new TF1("func_exp",fit_exp, 0., 25000, 4);
  }
*/
TF1 *func_exp = new TF1("func_exp",fit_exp, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);
TF1 *func_exp_2 = new TF1("func_exp_2",fit_exp, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);

TF1 *func_lin = new TF1("func_lin",fit_lin, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 2);
TF1 *func_lin_2 = new TF1("func_lin_2",fit_lin, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 2);

Double_t fit_pwr(Double_t *X, Double_t *par) 
{
  Double_t pwr = par[0] * pow(X[0],par[1]) + par[2];
  return pwr;
}

Double_t fit_exp_ln(Double_t *X, Double_t *par) 
{
  Double_t exp_ln = par[0] * TMath::Exp(par[1]*(X[0]-par[2])) + par[3]*TMath::Log(par[4]*(X[0]-par[2])) + par[5];
  return exp_ln;
}

//Create Gaussian to fit the timing resolution.
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

void Time_Walk(Int_t run = 501)
//Int_t Test_Display(Int_t run = 290)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Create a date object.
  TDatime time;

  FILE *fp;
  fp = fopen("/home/skbarcus/JLab/SBS/HCal/Analysis/Cosmics/Pedestals.txt","r");
 //Read in data.
  while (1) 
    {
      //Skips the first skip lines of the file. 
      if (nlines < skip)
	{
	  fgets(str,1000,fp);
	  nlines++;
	}
      //Reads the two columns of data into x and y.
      else
	{
	  //Read in the number of columns of data in your data file. 
	  ncols = fscanf(fp,"%f %f %f", &mod_temp, &pedestal_temp, &threshold_temp);
	  
	  if (ncols < 0) break;    
	  
	  mod[nlines-skip] = mod_temp;
	  pedestal[nlines-skip] = pedestal_temp;
	  //Threshold is defined as 1/4 of the average peak height minus the pedestal during a real event. Then add back in the pedestal value since the Landau fit I use doesn't subtract the pedestal.
	  threshold[nlines-skip] = (threshold_temp-pedestal[nlines-skip])/4.0+pedestal[nlines-skip];

	  nlines++;
	}
    }

  for(Int_t i=0;i<channels;i++)
    {
      cout<<"Module "<<i<<": pedestal = "<<pedestal[i]<<" threshold = "<<threshold[i]<<endl;
    }
  
  if(!T) 
    { 
      T = new TChain("T");
      //T->Add(TString::Format("rootfiles/fadc_f1tdc_%d.root",run));//325

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
	  rootfile = basename + "_" + u + ".root";
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

  for(Int_t r = 0; r < kNrows; r++) 
    {
      for(Int_t c = 0; c < kNcols; c++) 
	{
	  histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);        //fADC histos for time walk fit calculation.
	  histos_vert[r][c] = MakeHisto_vert(r,c,DISP_FADC_SAMPLES);   //fADC histos for vertical cosmics time walk corrections.
	  histos_cfd_k[r][c] = MakeHisto_cfd_k(r,c,DISP_FADC_SAMPLES);   //fADC histos for CFD fADC scaled by cfd_thr.
	  histos_cfd_d[r][c] = MakeHisto_cfd_d(r,c,DISP_FADC_SAMPLES);   //fADC histos for CFD fADC inverted and delayed by cfd_delay.
	}
    }

  nevt = T->GetEntries();
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];
  Double_t pedestals[channels] = {};
  Int_t no_tdc[channels] = {};
  Int_t yes_tdc[channels] = {};
  Double_t fadc_int[channels] = {};
  Double_t peak_height[channels] = {};

  //Define Landau functions to fit the fADC peaks and find the time over threshold.
  TF1 **func_landau = new TF1*[channels];
  for(Int_t i=0;i<channels;i++)
    {
      func_landau[i] = new TF1("func_landau",fit_landau, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);
    }

  //Define a line to draw where the threshold is on an fADC plot.
  Double_t threshold_line = 0;

  TLine **line = new TLine*[channels];
  for(Int_t i=0;i<channels;i++)
    {
      line[i] = new TLine(threshold_line,0,threshold_line,10000);
    }

  // Create a 2D histogram to hold the time walk correction results.
  //TH2F *htime_walk = new TH2F("htime_walk","",1000,0.,25000.,1000,0.,80.);
  TH2F **htime_walk = new TH2F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      htime_walk[i] = new TH2F(Form("htime_walk%d",i),Form("Time Walk Correction for Module %d",i),1000,-5000.,25000.,1000,-80.,80.);//,1000,0.,25000.,1000,0.,80.);
    }
  
  //Time walk corrected timing resolution plots.
  TH1F **htdc_uncorrected = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      htdc_uncorrected[i] = new TH1F(Form("htdc_uncorrected%d",i),Form("TDC Times with Ref Time (No Time Walk Correction) for Module %d",i),1000,-50,30);
    }

  //Time walk corrected timing resolution plots.
  TH1F **htdc_corrected = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      htdc_corrected[i] = new TH1F(Form("htdc_corrected%d",i),Form("TDC Times with Ref Time and Time Walk Corrected for Module %d",i),1000,-50,30);
    }

  //Timing resolution of the ref ch.
  TH1F *hrefch_timing = new TH1F("hrefch_timing","Reference Channel Timing Resolution",1000,0.,40.);

  //Create histo to contain the timing resolution between the two ref channels.
  TH1F *Ref_Timing_Resolution = new TH1F("Ref_Timing_Resolution","Double Reference Channel Timing Resolution",1000,-10.,10.);

  TFitResultPtr fit;

  //Loop over all events.
  if(limit_evts==1)
    {
      loop_max = max_evts;
    }
  else
    {
      loop_max = nevt;
    }
  //cout<<"*************"<<endl;
  for(Int_t i=0; i<loop_max ;i++)
    {
      // cout<<"*************"<<endl;
      T->GetEntry(gCurrentEntry);
      //      cout<<"!!!!!!!!!!!!!"<<endl;
      //Set arrays to zero;
      for(r  = 0; r < kNrows; r++) 
	{
	  for(c  = 0; c < kNcols; c++) 
	    {
	      adc[r][c] = 0.0;
	      tdc[r][c] = 0.0;
	    }
	}

      //Find the reference time by fitting the fADC pulse in the reference channel.
      idx = hcalt::samps_idx[ref_ch];
      n = hcalt::nsamps[ref_ch];
      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
	{
	  histos[ref_row][ref_col]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	  //cout<<"s"<<s<<" = "<<hcalt::samps[idx+s]<<endl;
	}

      //Draw the reference channel fADC.
      Double_t lower_bound = 0;
      Double_t upper_bound = DISP_MAX_SAMPLE;

      if(ref_fit_type == 1)
	{
	  //Find the lower bound of the fit range. Defined as when ref pulse amplitude passes a given threshold.
	  for(Int_t m=1;m<DISP_MAX_SAMPLE;m++)
	    {
	      if(histos[ref_row][ref_col]->GetBinContent(m) > lower_fit_th)
		//if(histos[ref_row][ref_col]->GetBinContent(m) > histos[ref_row][ref_col]->GetBinContent(m-1)*1.1 && m > 1)
		{
		  lower_bound = histos[ref_row][ref_col]->GetBinCenter(m);
		  break;
		}
	    }
	  
	  //Find the upper bound of the fit range.
	  //upper_bound = histos[ref_row][ref_col]->GetBinCenter(histos[ref_row][ref_col]->GetMaximumBin());
	  upper_bound = lower_bound + fit_width;
	  
	  //cout<<"Lower range bound = "<<lower_bound<<". Upper range bound = "<<upper_bound<<"."<<endl;
	  func_exp->SetRange(lower_bound,upper_bound);
	  func_exp->SetParameter(0,-2224);//-2224
	  func_exp->SetParameter(1,2.15);//2.15
	  func_exp->SetParameter(2,-0.4258);//-0.4258
	  func_exp->SetParameter(3,2682);//2682
	  func_exp->SetLineColor(1);
	  histos[ref_row][ref_col]->Fit(func_exp,"q 0 r");
	}

      //Linear fit to ref channel.
      Double_t lower_bound_lin = 0;
      Double_t upper_bound_lin = DISP_MAX_SAMPLE;

      if(ref_fit_type == 2)
	{
	  //Find the lower bound of the fit range. Find when height > 1200 as upper range and bin below that as lower (just fit two bins).
	  for(Int_t m=1;m<DISP_MAX_SAMPLE;m++)
	    {
	      if(histos[ref_row][ref_col]->GetBinContent(m) > ref_th)
		//if(histos[ref_row][ref_col]->GetBinContent(m) > histos[ref_row][ref_col]->GetBinContent(m-1)*1.1 && m > 1)
		{
		  upper_bound_lin = histos[ref_row][ref_col]->GetBinCenter(m);
		  lower_bound_lin = upper_bound_lin - 1;
		  //cout<<"lower_lin = "<<lower_bound_lin<<". upper_lin = "<<upper_bound_lin<<endl;
		  break;
		}
	    }
	  func_lin->SetRange(lower_bound_lin,upper_bound_lin);
	  func_lin->SetParameter(0,750);//750
	  func_lin->SetParameter(1,-2000);//-2000
	  histos[ref_row][ref_col]->Fit(func_lin,"q 0 r");
	}
      
      if(ref_trig_plots == 1)
	{
	  histos[ref_row][ref_col]->Draw();
	  if(ref_fit_type == 1)
	    {
	      func_exp->Draw("same");
	    }
	  if(ref_fit_type == 2)
	    {
	      func_lin->Draw("same");
	    }
	}
      //cout<<"func_exp: Chi^2 = "<<func_exp->GetChisquare()<<", NDF = "<<func_exp->GetNDF()<<", rChi^2 = "<<func_exp->GetChisquare()/func_exp->GetNDF()<<"."<<endl;
      //Find a point on the fit to use as a reference. The shape and fit should be consistent enough from event to event that a reference y-value should be sufficient. 

      Double_t ref_time = 0;
      if(ref_fit_type == 1)
	{
	  ref_time = func_exp->GetX(ref_th,lower_bound,upper_bound) * fadc_res;
	  //cout<<"Ref time is "<<ref_time<<"."<<endl;
	}
      if(ref_fit_type == 2)
	{
	  //linear ref time.
	  ref_time = func_lin->GetX(ref_th,lower_bound_lin,upper_bound_lin) * fadc_res;
	}


      //Repeat reference time calculation if using two reference channels. (Rare.)
      //Find the reference time by fitting the fADC pulse in the reference channel.
      idx = hcalt::samps_idx[ref_ch_2];
      n = hcalt::nsamps[ref_ch_2];
      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
	{
	  histos[ref_row_2][ref_col_2]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	  //cout<<"s"<<s<<" = "<<hcalt::samps[idx+s]<<endl;
	}

      if(ref_fit_type == 1)
	{
	  //Draw the reference channel fADC.
	  lower_bound = 0;
	  upper_bound = DISP_MAX_SAMPLE;
	  
	  //Find the lower bound of the fit range. Defined as when ref pulse amplitude passes 1000.
	  for(Int_t m=1;m<DISP_MAX_SAMPLE;m++)
	    {
	      if(histos[ref_row_2][ref_col_2]->GetBinContent(m) > lower_fit_th)
		//if(histos[ref_row_2][ref_col_2]->GetBinContent(m) > histos[ref_row_2][ref_col_2]->GetBinContent(m-1)*1.1 && m > 1)
		{
		  lower_bound = histos[ref_row_2][ref_col_2]->GetBinCenter(m);
		  break;
		}
	    }
	  
	  
	  //Find the upper bound of the fit range.
	  //upper_bound = histos[ref_row_2][ref_col_2]->GetBinCenter(histos[ref_row_2][ref_col_2]->GetMaximumBin());
	  upper_bound = lower_bound + fit_width;
	  
	  //cout<<"Lower range bound = "<<lower_bound<<". Upper range bound = "<<upper_bound<<"."<<endl;
	  func_exp_2->SetRange(lower_bound,upper_bound);
	  func_exp_2->SetParameter(0,-2224);//-2224
	  func_exp_2->SetParameter(1,2.15);//2.15
	  func_exp_2->SetParameter(2,-0.4258);//-0.4258
	  func_exp_2->SetParameter(3,2682);//2682
	  histos[ref_row_2][ref_col_2]->Fit(func_exp_2,"q 0 r");
	  histos[ref_row_2][ref_col_2]->SetLineColor(2);
	  func_exp_2->SetLineColor(1);
	}

      if(ref_fit_type == 2)
	{
	  //Linear fit to 2nd ref channel.
	  lower_bound_lin = 0;
	  upper_bound_lin = DISP_MAX_SAMPLE;
	  
	  //Find the lower bound of the fit range. Find when height > 1200 as upper range and bin below that as lower (just fit two bins).
	  for(Int_t m=1;m<DISP_MAX_SAMPLE;m++)
	    {
	      if(histos[ref_row_2][ref_col_2]->GetBinContent(m) > ref_th)
		//if(histos[ref_row_2][ref_col_2]->GetBinContent(m) > histos[ref_row_2][ref_col_2]->GetBinContent(m-1)*1.1 && m > 1)
		{
		  upper_bound_lin = histos[ref_row_2][ref_col_2]->GetBinCenter(m);
		  lower_bound_lin = upper_bound_lin - 1;
		  //cout<<"lower_lin = "<<lower_bound_lin<<". upper_lin = "<<upper_bound_lin<<endl;
		  break;
		}
	    }
	  func_lin_2->SetRange(lower_bound_lin,upper_bound_lin);
	  func_lin_2->SetParameter(0,750);//750
	  func_lin_2->SetParameter(1,-2000);//-2000
	  func_lin_2->SetLineColor(1);
	  histos[ref_row_2][ref_col_2]->Fit(func_lin_2,"q 0 r");
	  histos[ref_row_2][ref_col_2]->SetLineColor(2);
	}

      if(ref_trig_plots == 1)
	{
	  histos[ref_row_2][ref_col_2]->Draw("same");
	  if(ref_fit_type == 1)
	    {
	      func_exp_2->Draw("same");
	    }
	  if(ref_fit_type == 2)
	    {
	      func_lin_2->Draw("same");
	    }
	}
      //cout<<"func_exp_2: Chi^2 = "<<func_exp_2->GetChisquare()<<", NDF = "<<func_exp_2->GetNDF()<<", rChi^2 = "<<func_exp_2->GetChisquare()/func_exp_2->GetNDF()<<"."<<endl;

      //Find a point on the fit to use as a reference. The shape and fit should be consistent enough from event to event that a reference y-value should be sufficient. 
      Double_t ref_time_2 = 0;
      if(ref_fit_type == 1)
	{
	  //exp ref time.
	  ref_time_2 = func_exp_2->GetX(ref_th,lower_bound,upper_bound) * fadc_res;
	}
      if(ref_fit_type == 2)
	{
	  //linear ref time.
	  ref_time_2 = func_lin_2->GetX(ref_th,lower_bound_lin,upper_bound_lin) * fadc_res;
	}
      //cout<<"Ref time 2 is "<<ref_time_2<<"."<<endl;

      //Plot histo to contain the timing resolution between the two ref channels.
      Ref_Timing_Resolution->Fill(ref_time-ref_time_2);


      //Loop over all channels and fill the arrays. Do this seperately so reference time is set for all channels.
      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];
	  adc[r][c] = hcalt::a[j];
	  tdc[r][c] = hcalt::tdc[j];
	  fadc_int[j] = 0;
	  Double_t ped_val = 0;
	  
	  //Make histograms of fADC pulses when the TDC also fired.
	  if(tdc[r][c] != 0)
	    {
	      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
		{
		  if(use_cfd == 2)
		    {
		      //Fill histo with fADC scaled by cfd_thr constant value.
		      histos_cfd_k[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]*cfd_thr-pedestal[j]*cfd_thr);
		      //Fill histo with inverted fADC signal delayed by cfd_delay fADC time steps.
		      if(s > (cfd_delay - 1))
			{
			  //Fill delayed and inverted bins with fADC signal.
			  histos_cfd_d[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,-1*hcalt::samps[idx+s-cfd_delay]+pedestal[j]);
			}
		      else
			{
			  //Fill first few bins with zeros so the histogram has all samples filled.
			  histos_cfd_d[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,0);
			}
		    }
		  else
		    {
		      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
		      fadc_int[j] = fadc_int[j] + histos[r][c]->GetBinContent(s+1-DISP_MIN_SAMPLE);
		    }
		}
	      //sum the scaled and delayed histos for CFD zero crossing.
	      histos_cfd_k[r][c]->Add(histos_cfd_k[r][c],histos_cfd_d[r][c],1.,1.);

	      //Subtract off avg pedestals.
	      fadc_int[j] = fadc_int[j] - pedestal[j]*DISP_MAX_SAMPLE;

	      //Fit fADC pulse with a Landau function.
	      //func_landau[0]->GetParameter(0);
	      func_landau[j]->SetParameter(0,histos[r][c]->GetMaximum()*4.0);//Parameter setting the height scale (not 1:1). ~4800
	      //func_landau[j]->SetParLimits(0,0,1000000);                     //Force parameter 0 to be positive.
	      func_landau[j]->SetParameter(1,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()));//Parameter finding the approximate location of the peak. ~15
	      //func_landau[j]->SetParLimits(1,0,25);
	      func_landau[j]->SetParameter(2,1);//Parameter setting the width. ~1
	      func_landau[j]->SetParLimits(2,0,5);                             //Force parameter 2 to be positive.
	      func_landau[j]->SetParameter(3,pedestal[j]);//Parameter measuring the pedestal. ~250
	      fit = histos[r][c]->Fit(func_landau[j],"0 Q S opt");//Don't use M or get some issues that make code take forever.//opt in options was here for some reason I don't recall?

	      fit_status = histos[r][c]->Fit(func_landau[j],"Q");
	      //cout<<"Event "<<i<<" Module "<<j<<": fit status = "<<fit_status<<endl;

	      //func_landau[j]->Result();
	      //Fitter::Result() = fit_result;
	      //cout<<"Status = "<<TFitResultPtr<<endl;

	      //Int_t fitstatus = fit;                  //Gives zero when Landau fit works or fails.
	      //cout<<"Convergence = "<<fitstatus<<endl;

	      //Fill histogram with the fadc integral and the threshold times.
	      //htime_walk[j]->Fill(fadc_int[j],func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res);
	      //Combination of fADC and TDC. Must move into separate loop below so all tdc[r][c] are filled first.
	      //htime_walk[j]->Fill(fadc_int[j],(func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res) + (tdc[r][c] - (tdc[r+1][c]+tdc[r-1][c])/2.)*f1_res)/2.;
	      if(j==0)
		{
		  //cout<<"Event "<<i<<": Approximate peak location = "<<histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin())<<endl;
		  //cout<<"module "<<j<<" fadc_int = "<<fadc_int[j]<<" threshold x value (ns) = "<<func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res<<endl;
		}
	    }
	}

      //Draw CFD method 2 single fADC histogram.
      histos_cfd_k[mod_row-1][mod_col-1]->Draw("");
      //histos_cfd_d[mod_row-1][mod_col-1]->Draw("");

      //Draw a single fADC peak with Landau fit. run 820
      if(plot_ind_landau == 1)
	{
	  if(i == max_evts - 1)
	    {
	      module = (mod_row-1)*kNcols + mod_col-1;
	      cout<<"Event "<<i<<" mod_row-1 = "<<mod_row-1<<" mod_col-1 = "<<mod_col-1<<" module = "<<module<<endl;
	      cout<<"threshold[j] crossing = "<<func_landau[module]->GetX(threshold[module],0.,histos[mod_row-1][mod_col-1]->GetBinCenter(histos[mod_row-1][mod_col-1]->GetMaximumBin()))<<endl;
	      threshold_line = func_landau[module]->GetX(threshold[module],0.,histos[mod_row-1][mod_col-1]->GetBinCenter(histos[mod_row-1][mod_col-1]->GetMaximumBin()));
	      line[module]->SetX1(threshold_line);
	      line[module]->SetY1(0.);
	      line[module]->SetX2(threshold_line);
	      line[module]->SetY2(10000.);
	      histos[mod_row-1][mod_col-1]->Draw();
	      func_landau[module]->Draw("same");
	      line[module]->Draw("same");
	    }
	}

      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  //if(tdc[r][c] != 0 && func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res > 1.) //&& r != 0 && r != 11)
	  if(tdc[r][c] != 0)
	    {
	      //Fill time walk histograms.
	      //fADC and TDC.
	      //htime_walk[j]->Fill(fadc_int[j],(func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res) + ((tdc[r][c] - (tdc[r+1][c]+tdc[r-1][c])/2.)*f1_res)/2.);
	      //fADC int vs TDC with ref.
	      //	      htime_walk[j]->Fill(fadc_int[j],(tdc[r][c] - (tdc[r+1][c] + tdc[r-1][c])/2.)*f1_res);
	      //fADC only.
	      //No ref times.
	      //cout<<"PMT: "<<j<<" r = "<<r<<" r-1 = "<<r-1<<endl;
	      //	      htime_walk[j]->Fill(fadc_int[j],func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res);
	      //Ref times.
	      //htime_walk[j]->Fill(fadc_int[j],(   func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin())) - (func_landau[j-12]->GetX(threshold[j-12],0.,histos[r-1][c]->GetBinCenter(histos[r-1][c]->GetMaximumBin())) + func_landau[j+12]->GetX(threshold[j+12],0.,histos[r+1][c]->GetBinCenter(histos[r+1][c]->GetMaximumBin())))/2.   )*fadc_res);

	      //fADC time over threshold vs TDC time with ref time.
	      //htime_walk[j]->Fill(func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res,(tdc[r][c] - (tdc[r+1][c] + tdc[r-1][c])/2.)*f1_res);

	      //Using proper fADC trigger copy reference time.
	      //Condition to make sure the fADC peak amplitude is above the 1/4 average peak height threshold.
	      //cout<<"Event "<<i<<" Module "<<j<<": Histo_vert["<<r<<"]["<<c<<"] Apmlitude = "<<histos_vert[r][c]->GetMaximum()<<" Threshold = "<<threshold[j]<<endl;
	      if(use_cfd == 1 && histos[r][c]->GetMaximumBin() > 2 && histos[r][c]->GetMaximum() > (pedestal[j] + 10) && fit_status == 0 && TMath::IsNaN(func_landau[j]->GetX((histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res) == 0)
		{
		  //cout<<"Event "<<i<<" Module "<<j<<" (histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j] "<<(histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j]<<endl;
		  //cout<<"Event "<<i<<" Module "<<j<<": histo max = "<<histos[r][c]->GetMaximum()<<" Landau fit max = "<<func_landau[j]->GetMaximum()<<" histo max - landau fit max = "<<histos[r][c]->GetMaximum() - func_landau[j]->GetMaximum()<<endl;
		  //cout<<"Event "<<i<<" Module "<<j<<": IsNaN = "<<TMath::IsNaN(func_landau[j]->GetX((histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res)<<endl;		  
		  
		  htime_walk[j]->Fill(fadc_int[j],func_landau[j]->GetX((histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res - ref_time);
		}
	      else
		{
		  if(histos[r][c]->GetMaximum() > threshold[j]*1.5 && histos[r][c]->GetMaximumBin() > 3 && histos[r][c]->GetBinContent(1) < pedestal[j]*1.5)
		    {//cout<<"Hello "<<endl;
		      //cout<<"Event "<<i<<" Module "<<j<<": Histo["<<r<<"]["<<c<<"] Apmlitude = "<<histos[r][c]->GetMaximum()<<" Threshold = "<<threshold[j]<<endl;
		      htime_walk[j]->Fill(fadc_int[j],func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res - ref_time);
		      //cout<<"there!"<<endl;
		    }
		}
	    }
	}
      if(calc_corr_time_res == 1)
	{
	  if(gCurrentEntry%5000==0)
	    {
	      cout<<gCurrentEntry<<" events processed. "<<((double)gCurrentEntry/(double)loop_max)*50<<" % complete."<<endl;
	    }
	}
      else
	{
	  if(gCurrentEntry%5000==0)
	    {
	      cout<<gCurrentEntry<<" events processed. "<<((double)gCurrentEntry/(double)loop_max)*100<<" % complete."<<endl;
	    }
	}
      //cout<<"loop # = "<<gCurrentEntry<<endl;
      gCurrentEntry++;
    }

  //Draw the relative timing resolution between the two reference channels.
  if(ref_trig_plots == 2)
    {
      Ref_Timing_Resolution->Draw();
    }

  //Define functions to fit the time walk corrections.
  //log
  TF1 **func_log = new TF1*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      func_log[i] = new TF1("func_log",fit_log, 0., 25000, 4);
    }

  //Third order polynomial.
  TF1 **func_poly3 = new TF1*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      func_poly3[i] = new TF1("func_poly3",fit_poly3, 0., 7000, 4);
    }
  //Fourth order polynomial.
  TF1 **func_poly4 = new TF1*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      func_poly4[i] = new TF1("func_poly4",fit_poly4, 0., 25000, 5);
    }
  
  //Exponential.
  TF1 **func_exp_time_walk = new TF1*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      func_exp_time_walk[i] = new TF1("func_exp_time_walk",fit_exp, 0., 25000, 4);
    }
  
  //Exponential + ln.
  TF1 **func_exp_ln = new TF1*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      func_exp_ln[i] = new TF1("func_exp_ln",fit_exp_ln, 0., 25000, 6);
    }

  //Power.
  TF1 **func_pwr = new TF1*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      func_pwr[i] = new TF1("func_pwr",fit_pwr, 0., 25000, 3);
    }
  
  if(plot_time_walks == 1)
    {
      //Plot the time walk correction plots for each pmt.
      TCanvas *ctime_walk[kNrows];
      for(Int_t i=0; i<kNrows; i++)
	{
	  ctime_walk[i] = new TCanvas(Form("ctime_walk_Row_%d",i));
	  ctime_walk[i]->SetGrid();
	  ctime_walk[i]->Divide(4,3);
	  
	  for(Int_t j=0;j<12;j++)
	    {
	      ctime_walk[i]->cd(j+1);
	      //gPad->SetLogx();
	      htime_walk[i*12+j]->Draw("");
	      //htime_walk[i*12+j]->Fit(func_log[i*12+j],"Q");
	      //htime_walk[i*12+j]->Fit(func_poly3[i*12+j],"QR");
	      //htime_walk[i*12+j]->Fit(func_poly4[i*12+j],"Q");
	      //	  func_exp[i*12+j]->SetParameter(0,10);
	      //	  func_exp[i*12+j]->SetParLimits(0,1,20);
	      //	  func_exp[i*12+j]->SetParameter(1,-0.0001);
	      //	  func_exp[i*12+j]->SetParLimits(1,-0.001,0);
	      //	  func_exp[i*12+j]->SetParameter(2,1);
	      //	  func_exp[i*12+j]->SetParLimits(2,-10000,10000);
	      //	  func_exp[i*12+j]->SetParameter(3,10);
	      //	  func_exp[i*12+j]->SetParLimits(3,0,30);
	      //htime_walk[i*12+j]->Fit(func_exp[i*12+j],"Q");
	      //func_exp[i*12+j]->SetParameter(4,10);
	      //func_exp[i*12+j]->Draw("same");\\Pars (10,-0.0001,10) close start.
	      //func_pwr[i*12+j]->SetParameter(0,10000);
	      //func_pwr[i*12+j]->SetParameter(1,-1);
	      //func_pwr[i*12+j]->SetParameter(2,0);
	      //htime_walk[i*12+j]->Fit(func_pwr[i*12+j],"Q");
	      
	      //Fit with an exponential and a ln combined.
	      /*
		func_exp_ln[i*12+j]->SetParameter(0,10);
		func_exp_ln[i*12+j]->SetParLimits(0,1,30);
		func_exp_ln[i*12+j]->SetParameter(1,-0.0001);
		func_exp_ln[i*12+j]->SetParLimits(1,-0.001,0);
		func_exp_ln[i*12+j]->SetParameter(2,0);
		func_exp_ln[i*12+j]->SetParLimits(2,-500,500);
		func_exp_ln[i*12+j]->SetParameter(3,-1);
		func_exp_ln[i*12+j]->SetParLimits(3,0,-10);
		func_exp_ln[i*12+j]->SetParameter(4,0.05);
		func_exp_ln[i*12+j]->SetParLimits(4,0,0.1);
		func_exp_ln[i*12+j]->SetParameter(5,0);
		func_exp_ln[i*12+j]->SetParLimits(5,-30,30);
		//func_exp_ln[i*12+j]->Draw("same");
		cout<<"Fit of PMT "<<i*12+j<<"."<<endl;
		htime_walk[i*12+j]->Fit(func_exp_ln[i*12+j],"0");
		func_exp_ln[i*12+j]->Draw("same");
		func_exp_ln[i*12+j]->SetLineColor(2);
	      */
	      
	      //Fit with just an exponential.
	      cout<<"Fit of PMT "<<i*12+j<<"."<<endl;
	      func_exp_time_walk[i*12+j]->SetParameter(0,5.2);//5.2
	      func_exp_time_walk[i*12+j]->SetParLimits(0,0,20);
	      func_exp_time_walk[i*12+j]->SetParameter(1,1);//1
	      func_exp_time_walk[i*12+j]->SetParameter(2,-0.00049);//-0.00049
	      func_exp_time_walk[i*12+j]->SetParLimits(2,-0.01,0);
	      func_exp_time_walk[i*12+j]->SetParameter(3,-27);//-27
	      func_exp_time_walk[i*12+j]->SetParLimits(3,-40,40);
	      htime_walk[i*12+j]->Fit(func_exp_time_walk[i*12+j],"0");
	      func_exp_time_walk[i*12+j]->Draw("same");
	      func_exp_time_walk[i*12+j]->SetLineColor(2);
	    }	  
	}
    }
  gPad->SetLogx(0);
  
  if(calc_corr_time_res == 1)
    {
      //Loop to find vertical cosmics and get TDC times using a reference time. Then corrects for time walk calculated above.
      //Reset gCurrentEntry for next loop.
      gCurrentEntry = 0;
      //Loop over events.
      for(Int_t i=0; i<loop_max ;i++)
	{
	  T->GetEntry(gCurrentEntry);
	  
	  //Set arrays to zero;
	  for(r = 0; r < kNrows; r++) 
	    {
	      for(c = 0; c < kNcols; c++) 
		{
		  //peak[r][c] = 0.0;
		  adc[r][c] = 0.0;
		  tdc[r][c] = 0.0;
		}
	    }
	  
	  //Find the reference time by fitting the fADC pulse in the reference channel.
	  idx = hcalt::samps_idx[ref_ch];
	  n = hcalt::nsamps[ref_ch];
	  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
	    {
	      histos_vert[ref_row][ref_col]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	      //cout<<"s"<<s<<" = "<<hcalt::samps[idx+s]<<endl;
	    }
	  
	  //Draw the reference channel fADC.
	  Double_t lower_bound_vert = 0;
	  Double_t upper_bound_vert = DISP_MAX_SAMPLE;
	  
	  //Find the lower bound of the fit range. Defined as when ref pulse amplitude passes 1000.
	  for(Int_t m=1;m<DISP_MAX_SAMPLE;m++)
	    {
	      if(histos_vert[ref_row][ref_col]->GetBinContent(m) > lower_fit_th)
		{
		  lower_bound_vert = histos_vert[ref_row][ref_col]->GetBinCenter(m);
		  break;
		}
	    }
	  
	  //Find the upper bound of the fit range.
	  //upper_bound_vert = histos_vert[ref_row][ref_col]->GetBinCenter(histos_vert[ref_row][ref_col]->GetMaximumBin());
	  upper_bound_vert = lower_bound_vert + fit_width;
	  
	  //cout<<"Lower range bound = "<<lower_bound<<". Upper range bound = "<<upper_bound<<"."<<endl;
	  func_exp->SetRange(lower_bound_vert,upper_bound_vert);
	  func_exp->SetParameter(0,-2224);//-2224
	  func_exp->SetParameter(1,2.15);//2.15
	  func_exp->SetParameter(2,-0.4258);//-0.4258
	  func_exp->SetParameter(3,2682);//2682
	  histos_vert[ref_row][ref_col]->Fit(func_exp,"q 0 r");
	  //histos_vert[ref_row][ref_col]->Draw();
	  //Find a point on the fit to use as a reference. The shape and fit should be consistent enough from event to event that a reference y-value should be sufficient. 
	  Double_t ref_time_vert = func_exp->GetX(ref_th,lower_bound_vert,upper_bound_vert) * fadc_res;
	  //cout<<"Ref time is "<<ref_time_vert<<"."<<endl;
	  //Fill ref ch timing resolution histo.
	  hrefch_timing->Fill(ref_time_vert);
	  
	  //Loop over tubes. Do this first so all PMTs have TDC times etc. 
	  for(Int_t j=0; j<hcalt::ndata; j++)
	    {
	      //Reset the fADC integral.
	      fadc_int[j] = 0;
	      r = hcalt::row[j]-1;
	      c = hcalt::col[j]-1;
	      idx = hcalt::samps_idx[j];
	      n = hcalt::nsamps[j];
	      adc[r][c] = hcalt::a[j];
	      tdc[r][c] = hcalt::tdc[j];
	      //Fill fADC histos for each PMT in event.
	      if(tdc[r][c] != 0)
		{
		  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
		    {
		      histos_vert[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
		      fadc_int[j] = fadc_int[j] + histos_vert[r][c]->GetBinContent(s+1-DISP_MIN_SAMPLE);
		    }
		  //Subtract off avg pedestals.
		  fadc_int[j] = fadc_int[j] - pedestal[j]*DISP_MAX_SAMPLE;
		  
		  //Fit fADC pulse with a Landau function.
		  //func_landau[0]->GetParameter(0);
		  func_landau[j]->SetParameter(0,histos_vert[r][c]->GetMaximum()*4.0);//Parameter setting the height scale (not 1:1). ~4800
		  func_landau[j]->SetParameter(1,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin()));//Parameter finding the location of the peak. ~15
		  func_landau[j]->SetParameter(2,1);//Parameter setting the width. ~1
		  func_landau[j]->SetParameter(3,250);//Parameter measuring the pedestal. ~250
		  histos_vert[r][c]->Fit(func_landau[j],"Q 0");//Don't use M or get some issues that make code take forever.
		}
	      //cout<<"Event "<<i<<" PMT "<<j<<": fadc_int = "<<fadc_int[j]<<endl;
	    }
	  //Loop over PMTs again now that adc/tdc values are filled and we have fADC integrals. We will do time walk corrections here.
	  for(Int_t j=0; j<hcalt::ndata; j++)
	    {
	      r = hcalt::row[j]-1;
	      c = hcalt::col[j]-1;
	      idx = hcalt::samps_idx[j];
	      n = hcalt::nsamps[j];
	      //Find three PMTs firing vertically and apply time walk correction to TDC time using reference time of average times of PMTs directly above and below the PMT being studied.
	      //cout<<"Event "<<i<<" PMT "<<j<<": TDC[r][c] = "<<tdc[r][c]<<" TDC[r+1][c] = "<<tdc[r+1][c]<<" TDC[r-1][c] = "<<tdc[r-1][c]<<" fit(0) = "<<func_poly3[j]->Eval(0)<<" fit(fadc_int) = "<<func_poly3[j]->Eval(fadc_int[j])<<endl;
	      //if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && fadc_int[j]>fadc_int_max_cut && fadc_int[j-12]>fadc_int_max_cut && fadc_int[j+12]>fadc_int_max_cut && fadc_int[j-1]<fadc_int_min_cut && fadc_int[j-13]<fadc_int_min_cut && fadc_int[j+11]<fadc_int_min_cut && fadc_int[j+1]<fadc_int_min_cut && fadc_int[j-11]<fadc_int_min_cut && fadc_int[j+13]<fadc_int_min_cut)
	      //Three vertical tubes.
	      if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0)
		//Six vertical tubes.
		//if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0&& tdc[r-2][c] != 0 && tdc[r-3][c] != 0 && tdc[r-4][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r-2][c-1] == 0 && tdc[r-3][c-1] == 0 && tdc[r-4][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0 && tdc[r-2][c+1] == 0 && tdc[r-3][c+1] == 0 && tdc[r-4][c+1] == 0)
		{
		  //cout<<"Event "<<i<<" PMT "<<j<<": TDC[r][c] = "<<tdc[r][c]<<" TDC[r+1][c] = "<<tdc[r+1][c]<<" TDC[r-1][c] = "<<tdc[r-1][c]<<" fit(0) = "<<func_poly3[j]->Eval(0)<<" fit(fadc_int) = "<<func_poly3[j]->Eval(fadc_int[j])<<" tdc corrected = "<<( (tdc[r][c] - (tdc[r+1][c]+tdc[r-1][c])/2.) + func_poly3[j]->Eval(0) - func_poly3[j]->Eval(fadc_int[j]) ) *f1_res<<endl;
		  //htdc_corrected[j]->Fill( (tdc[r][c] - (tdc[r+1][c]+tdc[r-1][c])/2.)*f1_res );

		  //htdc_uncorrected[j]->Fill((   func_landau[j]->GetX(threshold[j],0.,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin())) - (func_landau[j-12]->GetX(threshold[j-12],0.,histos_vert[r-1][c]->GetBinCenter(histos_vert[r-1][c]->GetMaximumBin())) + func_landau[j+12]->GetX(threshold[j+12],0.,histos_vert[r+1][c]->GetBinCenter(histos_vert[r+1][c]->GetMaximumBin())))/2.   )*fadc_res);
		  //htdc_corrected[j]->Fill((   func_landau[j]->GetX(threshold[j],0.,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin())) - (func_landau[j-12]->GetX(threshold[j-12],0.,histos_vert[r-1][c]->GetBinCenter(histos_vert[r-1][c]->GetMaximumBin())) + func_landau[j+12]->GetX(threshold[j+12],0.,histos_vert[r+1][c]->GetBinCenter(histos_vert[r+1][c]->GetMaximumBin())))/2.   )*fadc_res - func_exp_ln[j]->Eval(fadc_int[j]));
		  
		  //Perform time walk corrections using the correct trigger copy as a ref channel.
		  //Condition to make sure the fADC peak amplitude is above the 1/4 average peak height threshold.
		  //cout<<"Event "<<i<<" Module "<<j<<": Histo_vert["<<r<<"]["<<c<<"] Apmlitude = "<<histos_vert[r][c]->GetMaximum()<<" Threshold = "<<threshold[j]<<endl;
		  if(use_cfd == 1 && histos[r][c]->GetMaximumBin() > 2 && histos[r][c]->GetMaximum() > (pedestal[j] + 10) && fit_status == 0 && TMath::IsNaN(func_landau[j]->GetX((histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res) == 0)
		    {
		      //cout<<"Event "<<i<<" Module "<<j<<" (histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j] "<<(histos[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j]<<endl;
		      //cout<<"Event "<<i<<" Module "<<j<<": histo max = "<<histos[r][c]->GetMaximum()<<" Landau fit max = "<<func_landau[j]->GetMaximum()<<" histo max - landau fit max = "<<histos[r][c]->GetMaximum() - func_landau[j]->GetMaximum()<<endl;		  
		      
		      htdc_uncorrected[j]->Fill(func_landau[j]->GetX((histos_vert[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j],0.,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin()))*fadc_res - ref_time_vert);
		      htdc_corrected[j]->Fill(func_landau[j]->GetX((histos_vert[r][c]->GetMaximum()-pedestal[j])*cfd_thr+pedestal[j],0.,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin()))*fadc_res - ref_time_vert - func_exp_time_walk[j]->Eval(fadc_int[j]));
		    }
		  else
		    {
		      if(histos_vert[r][c]->GetMaximum() > threshold[j]*1.5)
			{
			  htdc_uncorrected[j]->Fill((func_landau[j]->GetX(threshold[j],0.,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin()))*fadc_res) - ref_time_vert);
			  htdc_corrected[j]->Fill((func_landau[j]->GetX(threshold[j],0.,histos_vert[r][c]->GetBinCenter(histos_vert[r][c]->GetMaximumBin()))*fadc_res) - ref_time_vert - func_exp_time_walk[j]->Eval(fadc_int[j]));
			}
		    }
		}
	    }
	  
	  if(gCurrentEntry%5000==0)
	    {
	      cout<<gCurrentEntry<<" events processed twice. "<<50 + ((double)gCurrentEntry/(double)loop_max)*50<<" % complete."<<endl;
	    }
	  gCurrentEntry++;
	}
      
      TF1 **func_gaus_fit_uncorrected = new TF1*[channels];
      for(Int_t i = 0; i<hcalt::ndata; i++)
	{
	  func_gaus_fit_uncorrected[i] = new TF1("func_gaus_fit_uncorrected",fit_gaus, 0, 50, 3);
	}
      
      TF1 **func_gaus_fit = new TF1*[channels];
      for(Int_t i = 0; i<hcalt::ndata; i++)
	{
	  func_gaus_fit[i] = new TF1("func_gaus_fit",fit_gaus, -20, 20, 3);
	}

      TF1 *func_gaus_fit_refch = new TF1("func_gaus_fit_refch",fit_gaus, 0, 40, 3);
      
      gStyle->SetOptFit(1111);

      //Make plot of ref ch timing resolution.
      TCanvas* crefch_timing=new TCanvas("crefch_timing");
      crefch_timing->SetGrid();
      hrefch_timing->Draw();
      func_gaus_fit_refch->SetLineColor(4);
      func_gaus_fit_refch->SetNpx(1000);
      func_gaus_fit_refch->SetParameter(0,hrefch_timing->GetMaximum());
      func_gaus_fit_refch->SetParameter(1,hrefch_timing->GetMean());
      func_gaus_fit_refch->SetParameter(2,hrefch_timing->GetStdDev());
      hrefch_timing->Fit(func_gaus_fit_refch,"Q 0");
      func_gaus_fit_refch->Draw("same");
      cout<<"Reference channel = r"<<ref_row+1<<"c"<<ref_col+1<<": Initial Max "<<hrefch_timing->GetMaximum()<<" Initial Mean "<<hrefch_timing->GetMean()<<" Initial StdDev "<<hrefch_timing->GetStdDev()<<" Final Max "<<func_gaus_fit_refch->GetParameter(0)<<" Final Mean "<<func_gaus_fit_refch->GetParameter(1)<<" Final StdDev "<<func_gaus_fit_refch->GetParameter(2)<<endl;

  Int_t x1 = -20;
  Int_t y1 = 20;
  Int_t x2 = -10;
  Int_t y2 = 0;

  //Double_t sig_refch = func_gaus_fit_refch->GetParameter(2);  //Timing resolution of the ref ch.
  Double_t sig_refch = hrefch_timing->GetStdDev();
  Double_t sig_pmt = 0;                                       //Timing resolution of a single PMT.

  TPaveText **text = new TPaveText*[channels];
      
      //Plot the time walk correction plots for each pmt.
      TCanvas *ctime_corrected[kNrows];
      for(Int_t i=0; i<kNrows; i++)
	{
	  ctime_corrected[i] = new TCanvas(Form("ctime__corrected_Row_%d",i));
	  ctime_corrected[i]->SetGrid();
	  ctime_corrected[i]->Divide(4,3);
	  
	  for(Int_t j=0;j<12;j++)
	    {
	      ctime_corrected[i]->cd(j+1);
	      htdc_corrected[i*12+j]->Draw();
	      func_gaus_fit[i*12+j]->SetLineColor(4);
	      func_gaus_fit[i*12+j]->SetNpx(1000);
	      func_gaus_fit[i*12+j]->SetParameter(0,htdc_corrected[i*12+j]->GetMaximum());
	      func_gaus_fit[i*12+j]->SetParameter(1,htdc_corrected[i*12+j]->GetMean());
	      func_gaus_fit[i*12+j]->SetParameter(2,htdc_corrected[i*12+j]->GetStdDev());
	      htdc_corrected[i*12+j]->Fit(func_gaus_fit[i*12+j],"Q");
	      func_gaus_fit[i*12+j]->Draw("same");
	      cout<<"Module "<<i*12+j<<" Initial Max "<<htdc_corrected[i*12+j]->GetMaximum()<<" Initial Mean "<<htdc_corrected[i*12+j]->GetMean()<<" Initial StdDev "<<htdc_corrected[i*12+j]->GetStdDev()<<" Final Max "<<func_gaus_fit[i*12+j]->GetParameter(0)<<" Final Mean "<<func_gaus_fit[i*12+j]->GetParameter(1)<<" Final StdDev "<<func_gaus_fit[i*12+j]->GetParameter(2)<<endl;
	      
	      htdc_uncorrected[i*12+j]->Draw("sames");
	      htdc_uncorrected[i*12+j]->SetLineColor(2);
	      func_gaus_fit_uncorrected[i*12+j]->SetLineColor(2);
	      func_gaus_fit_uncorrected[i*12+j]->SetNpx(1000);
	      func_gaus_fit_uncorrected[i*12+j]->SetParameter(0,htdc_uncorrected[i*12+j]->GetMaximum());
	      func_gaus_fit_uncorrected[i*12+j]->SetParameter(1,htdc_uncorrected[i*12+j]->GetMean());
	      func_gaus_fit_uncorrected[i*12+j]->SetParameter(2,htdc_uncorrected[i*12+j]->GetStdDev());
	      htdc_uncorrected[i*12+j]->Fit(func_gaus_fit_uncorrected[i*12+j],"Q sames");
	      func_gaus_fit_uncorrected[i*12+j]->Draw("same");
	      cout<<"Module "<<i*12+j<<" Initial Max "<<htdc_uncorrected[i*12+j]->GetMaximum()<<" Initial Mean "<<htdc_uncorrected[i*12+j]->GetMean()<<" Initial StdDev "<<htdc_uncorrected[i*12+j]->GetStdDev()<<" Final Max "<<func_gaus_fit_uncorrected[i*12+j]->GetParameter(0)<<" Final Mean "<<func_gaus_fit_uncorrected[i*12+j]->GetParameter(1)<<" Final StdDev "<<func_gaus_fit_uncorrected[i*12+j]->GetParameter(2)<<endl;

	      ctime_corrected[i]->Update();
	      x1 = -50;
	      y1 =  gPad->GetUymax();
	      //y1 = func_gaus_fit[i*12+j]->GetMaximum(-DISP_MAX_SAMPLE,DISP_MAX_SAMPLE)*0.7;
	      x2 = -15;
	      y2 = gPad->GetUymax()*0.333;
	      text[i*12+j] = new TPaveText(x1,y1,x2,y2);
	      text[i*12+j]->AddText("Uncorrected StdDev: ");  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText(Form("%0.3f",func_gaus_fit_uncorrected[i*12+j]->GetParameter(2)));  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText("Corrected StdDev: ");  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText(Form("%0.3f",func_gaus_fit[i*12+j]->GetParameter(2)));  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText("Ref Ch StdDev: ");  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText(Form("%0.3f",sig_refch));  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText("PMT StdDev: ");  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText(Form("%0.3f",pow(fabs(pow(func_gaus_fit[i*12+j]->GetParameter(2),2.) - pow(sig_refch,2.)),0.5)));  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->Draw("same");
	    }
	}
    }
  
  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
