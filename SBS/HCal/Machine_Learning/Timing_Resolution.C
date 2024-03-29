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

Int_t timing_method = 2; //0 = TDC timing using trigger copy as ref ch. 1 = fADC timing using a Landau to find the peak location for each hit (require a TDC hit too). 2 = TDC timing with fADC cuts above and below target PMT requiring an fADC peak and the 6 surrounding tubes requiring no fADC peak. Ref ch is avg of above and below TDC times. 3 = fADC integrals. Takes integral of each fADC hit (requiring TDC hit too) and makes a histo for each PMT. 4 = Plot just one module fADC integral, TDC time spectrum with channel below as ref time, and same fADC integral with optional cut on TDC time. 5 = TDC time with reference time subtracted when ref is avg of above and below tubes vs. fadc integral. 6 = Top PMT of three vertically - bottom PMT's TDC times vs. fADC integrals. 7 = Calculate time over threshold using only fADCs and a Landau fit of the signal.
Double_t entry;
Int_t gCurrentEntry = 0;
Int_t DISP_MIN_SAMPLE = 0;
Int_t DISP_MAX_SAMPLE = 20;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
Int_t min_time = -300;
Int_t max_time = 50;
Int_t bins = 100;
Double_t fadc_int_max_cut = 1000;
Double_t fadc_int_min_cut = 500;
Int_t r,c,n,idx;
Int_t nevt;
Int_t refrow = 1000;     //Counting from 1. 
Int_t refcol = 1000;     //Counting from 1.
Int_t timerow = 1;    //Counting from 1.
Int_t timecol = 9;    //Counting from 1.
const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
Int_t limit_evts = 0;                     //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t loop_max = 0;                       //Dummy variable set equal to max_evts or nevt for the loop.
Int_t max_evts = 1000;                   //Maximum number of events to analyze if limit_evts = 1.
//const Int_t maxevts = 10000;
Double_t ind_res[144];
//vector<double> temp(0,0);
//vector<vector<double>> timing(144,temp);
//Double_t times[1000000][144] = {};
Double_t f1_res = 0.112;//0.112
Double_t fadc_res = 4.;
TChain *T = 0;
std::string user_input;
Int_t skip = 3;                            //Gives number of lines to skip at top of data file. 
Int_t nlines = 0;                          //Counts number of lines in the data file. 
Int_t ncols;                               //Set how many columns of data we have in the data file.
char str[1000];                           //Variable to read lines of the data file.
Float_t mod[channels], pedestal[channels], threshold[channels];
Float_t mod_temp, pedestal_temp, threshold_temp;

//Create Gaussian to fit the timing resolution.
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

TH1F *histos[kNrows][kNcols];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
		     TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Create Gaussian to fit the timing resolution.
Double_t fit_landau(Double_t *X, Double_t *par) 
{
  //Double_t landau =0;
  //Double_t landau = TMath::Gaus(par[0],par[1],par[2]);
  Double_t landau = par[0]*TMath::Landau(X[0],par[1],par[2])+par[3];
  return landau;
}

//TF1 **func_landau = new TF1*[channels];
//TF1 *func_landau[channels] = new TF1("func_landau",fit_landau, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 3);

void Timing_Resolution(Int_t run = 501)
//Int_t Test_Display(Int_t run = 290)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

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

  for(Int_t r = 0; r < kNrows; r++) 
    {
      for(Int_t c = 0; c < kNcols; c++) 
	{
	  histos[r][c] = MakeHisto(r,c,DISP_FADC_SAMPLES);
	}
    }

  //idx = hcalt::samps_idx[m];
  //n = hcalt::nsamps[m];
  //gCurrentEntry = entry;
  cout<<"Hello"<<endl;
  //T->GetEntry(0);
  //T->GetEntry(345345);
  T->GetEntry(gCurrentEntry);
  cout<<"There!"<<endl;
  cout << "Displaying event " << gCurrentEntry << endl;
  cout<<"ndata = "<<hcalt::ndata<<endl;                                  //Number of channels.
  cout<<"row[] = "<<hcalt::row[16]<<"   col[] = "<<hcalt::col[16]<<endl;
  cout<<"a[] = "<<hcalt::a[1]<<endl;
  cout<<"tdc[] = "<<hcalt::tdc[1]<<endl;

  n = hcalt::nsamps[0];
  cout<<"n = "<<n<<endl;

  idx = hcalt::samps_idx[0];
  cout<<"idx = "<<idx<<endl;

  nevt = T->GetEntries();
  cout<<"Number of events = "<<nevt<<endl;

  Float_t peak[kNrows][kNcols];
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];

  TF1 **func_landau = new TF1*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      func_landau[i] = new TF1("func_landau",fit_landau, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);
    }

  TH1F **htimes = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htimes[i] = new TH1F(Form("htimes%d",i),Form("TDC Timing for Module %d",i),bins,min_time,max_time);
    }
  
  TH1F **htimes_fadc = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htimes_fadc[i] = new TH1F(Form("htimes_fadc%d",i),Form("fADC Timing for Module %d",i),bins,DISP_MIN_SAMPLE+10,DISP_MAX_SAMPLE*fadc_res);
    }

  TH1F **hfadc_integrals = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      //Allow separate x range for ref chs. 
      if(i == 142 || i ==143)
	{
	  hfadc_integrals[i] = new TH1F(Form("hfadc_integrals%d",i),Form("Integral of fADC for Module %d",i),8000,-4000,4000);
	}
      else
	{
	  hfadc_integrals[i] = new TH1F(Form("hfadc_integrals%d",i),Form("Integral of fADC for Module %d",i),250,-1000,20000);//1000,-3000,20000
	}
    }

  TH1F **htdc_times = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htdc_times[i] = new TH1F(Form("htdc_times%d",i),Form("TDC Times with Ref Channel %d",i),61,-30,30);
    }

  TH1F **hfadc_integrals_tdc_cuts = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      hfadc_integrals_tdc_cuts[i] = new TH1F(Form("hfadc_integrals_tdc_cuts%d",i),Form("Integral of fADC for Module %d with Cut on TDC Time",i),1000,-3000,20000);
    }

  TH1F **htdc_fadc_cuts = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htdc_fadc_cuts[i] = new TH1F(Form("htdc_fadc_cuts%d",i),Form("TDC times with avg of above and below TDCs for ref time with cuts on integral of fADC for Module %d",i),bins,-DISP_MAX_SAMPLE,DISP_MAX_SAMPLE);
    }

  TH1F **htdc_ref = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htdc_ref[i] = new TH1F(Form("htdc_ref%d",i),Form("TDC Reference Time Jitter for Module %d",i),bins,-20,20);
    }

TH2F **htdc_time_vs_fadc_int = new TH2F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htdc_time_vs_fadc_int[i] = new TH2F(Form("htdc_time_vs_fadc_int%d",i),Form("TDC times with avg of above and below TDCs for ref time vs. fADC Integrals %d",i),1000,-1000,20000,1000,-30,30);
    }

TH2F **htdc_verticallity = new TH2F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      htdc_verticallity[i] = new TH2F(Form("htdc_verticallity%d",i),Form("TDC times of Top Minus Bottom PMTs vs. fADC Integrals %d",i),1000,-1000,20000,1000,-100,100);
    }

  TH1F **hfadc_thr_times = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      hfadc_thr_times[i] = new TH1F(Form("hfadc_thr_times%d",i),Form("fADC Time Over Threshold for Module %d",i),bins,-20,20);
    }

  /*
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  
  TH1F *h1 = new TH1F("h1","h1",100,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h1->SetStats(0);
  h1->SetLineWidth(2);
  //return h1;

  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++)  
    {
      //cout<<"hello"<<endl;
      h1->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
      //cout<<hcalt::samps[idx+s]<<endl;
    }
  h1->Draw();
  */

  /*
  TH1F *htiming = new TH1F("htiming",Form("Timing Resolution for Ch r%d-c%d (ref ch r%d-c%d time - raw time)",timerow,timecol,refrow,refcol),bins,min_time,max_time);
  //htiming->SetStats(0);
  htiming->SetLineWidth(2);
  */
  Int_t hits[144] = {};
  Double_t pedestals[144] = {};
  Int_t no_tdc[144] = {};
  
  //Loop over all events.
  if(limit_evts==1)
    {
      loop_max = max_evts;
    }
  else
    {
      loop_max = nevt;
    }
  //for(Int_t i=0; i<1 ;i++)
  for(Int_t i=0; i<loop_max ;i++)
    {
      T->GetEntry(gCurrentEntry);
      //T->GetEntry(11);
      //T->GetEntry(665);

      //Set arrays to zero;
      for(r  = 0; r < kNrows; r++) 
	{
	  for(c  = 0; c < kNcols; c++) 
	    {
	      peak[r][c] = 0.0;
	      adc[r][c] = 0.0;
	      tdc[r][c] = 0.0;
	    }
	}

      Double_t fadc_int[channels] = {};
      //Find reference time.
      //Loop over all channels and fill the arrays. Do this seperately so reference time is set for all channels.
      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];
	  adc[r][c] = hcalt::a[j];
	  tdc[r][c] = hcalt::tdc[j];
	  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
	    {
	      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	      fadc_int[j] = fadc_int[j] + histos[r][c]->GetBinContent(s+1-DISP_MIN_SAMPLE);
	    }
	  //Subtract off avg pedestals.
	  //cout<<"Event "<<gCurrentEntry<<" PMT from 0 "<<j<<": fADC integral with pedestal "<<fadc_int[j]<<" pedestal value "<<pedestal[j]<<" DISP_MAX_SAMPLE "<<DISP_MAX_SAMPLE<<endl;
	  fadc_int[j] = fadc_int[j] - pedestal[j]*DISP_MAX_SAMPLE;
	  //cout<<"fADC integral pedestal subtracted "<<fadc_int[j]<<endl;
	  //func_landau[j] = new TF1("func_landau",fit_landau, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);
	}

      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];
	  //adc[r][c] = hcalt::a[j];
	  //tdc[r][c] = hcalt::tdc[j];
	  
	  //	  Double_t fadc_int[channels] = {};
	  //Fill each PMT channel with its fADC data.
	  /*
	  for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
	    {
	      histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
	      fadc_int[j] = fadc_int[j] + histos[r][c]->GetBinContent(s+1-DISP_MIN_SAMPLE);
	    }
	  */
	  Double_t ped_val = 0;
	  //Get the pedestals.
	  if(tdc[r][c] == 0)
	    {
	      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++)
		{
		  ped_val = ped_val + histos[r][c]->GetBinContent(s+1);
		}
	      ped_val = ped_val/DISP_FADC_SAMPLES;
	      pedestals[j] = pedestals[j] + ped_val;
	      no_tdc[j] = no_tdc[j]+1;
	    }
	  //	  func_landau[j] = new TF1("func_landau",fit_landau, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 4);//This kills code if many evts.
	  if(timing_method == 1)
	    {
	      func_landau[j]->SetParameter(0,histos[r][c]->GetMaximum()*4.0);//Parameter setting the height scale (not 1:1). ~4800
	      func_landau[j]->SetParameter(1,15);//Parameter finding the location of the peak. ~15
	      func_landau[j]->SetParameter(2,1);//Parameter setting the width. ~1
	      func_landau[j]->SetParameter(3,250);//Parameter measuring the pedestal. ~250
	    }
	  //cout<<"histos->GetMaximum() = "<<histos[r][c]->GetMaximum()<<endl;

	  //func_landau[j]->SetParameter(0,histos[r][c]->GetMaximum());
	  //func_landau[j]->SetParameter(1,histos[r][c]->GetMean());
	  //func_landau[j]->SetParameter(2,histos[r][c]->GetStdDev());
	  //	  histos[r][c]->Fit("func_landau","QM");

	  //Fill the timing vector with the corrected times for each channel.
	  if(tdc[r][c] !=  0)   //Remove this TDC firing requirement for mode 4.
	  //Require 3 vertical PMTs.
	  //if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0)
	    {
	      //	      timing[j].push_back((tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res);
	      //	      times[i][j] = (tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res;
	      if((r != (refrow - 1) || c != (refcol - 1)))
		{
		  //cout<<"Event = "<<i<<"Module = "<<j<<endl;
		  //		  histos[r][c]->Fit("func_landau","M");
		  if(timing_method == 1)
		    {
		      histos[r][c]->Fit(func_landau[j],"Q");//Don't use M or get some issues that make code take forever.
		    }
		  if(timing_method == 3)
		    {	      
		      hfadc_integrals[j]->Fill(fadc_int[j]);
		      /*
		      if(j==142 || j==143)
			{
			  cout<<"Event # "<<i<<": PMT "<<j<<" has fadc_int[j] = "<<fadc_int[j]<<"."<<endl;
			}
		      */
		      /*
		      if(j==76 && fadc_int[j]>100)
			{
			  cout<<"Event # "<<i<<" fadc_int[j] = "<<fadc_int[j]<<endl;
			}
		      */
		    }
		  if(timing_method == 4)
		    {
		      if(tdc[r+1][c] != 0)
			{
			  hfadc_integrals[j]->Fill(fadc_int[j]);
			  htdc_times[j]->Fill((tdc[r][c]-tdc[r+1][c])*f1_res);
			  //Apply a cut on the TDC time to the fADC integral.
			  if((tdc[r][c]-tdc[r+1][c])*f1_res > -10 && (tdc[r][c]-tdc[r+1][c])*f1_res < 10 )
			    {
			      hfadc_integrals_tdc_cuts[j]->Fill(fadc_int[j]);
			    }
			}
		    }
		  //cout<<"Event "<<i<<" module "<<j<<". Peak position = "<<func_landau[j]->GetMaximumX()<<endl;
		  //cout<<"par[0] = "<<func_landau[j]->GetParameter(0)<<", par[1] = "<<func_landau[j]->GetParameter(1)<<", par[2] = "<<func_landau[j]->GetParameter(2)<<", par[3] = "<<func_landau[j]->GetParameter(3)<<endl;
		}
 	    } 
	  //Fill the histogram with the TDC times for the chosen channel (ref time - raw time).
	  //if(tdc[r][c] !=  0 && (r != (refrow - 1) || c != (refcol - 1)) && (r == (timerow - 1) && c == (timecol - 1)) && tdc[r+1][c] != 0 && tdc[r-1][c] != 0)
	  if(tdc[r][c] !=  0 && (r != (refrow - 1) || c != (refcol - 1)) && tdc[r+1][c] != 0 && tdc[r-1][c] != 0)
	    {
	      //	      htiming->Fill((tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res);
	      if(timing_method == 0)
		{
		  htimes[j]->Fill((tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res);
		}
	      if(timing_method == 1)
		{
		  htimes_fadc[j]->Fill((func_landau[j]->GetMaximumX())*fadc_res);
		}
	      hits[j] = hits[j] + 1;
	      //timing[j].push_back((tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res);
	      //cout<<"Event "<<gCurrentEntry<<": Corrected time for hit in r"<<r+1<<"-c"<<c+1<<" = "<<"(tdc["<<refrow<<"]["<<refcol<<"] - tdc["<<r+1<<"]["<<c+1<<"])*"<<f1_res<<" ns = "<<(tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res<<"."<<endl;
	      //cout<<"**** "<<(tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res<<" ****"<<endl;
	      if(timing_method == 5)
		{
		  htdc_time_vs_fadc_int[j]->Fill(fadc_int[j],(tdc[r][c] - (tdc[r+1][c]+tdc[r-1][c])/2.)*f1_res);
		}
	      if(timing_method == 6)
		{
		  htdc_verticallity[j]->Fill(fadc_int[j],(tdc[r+1][c]-tdc[r-1][c])*f1_res);
		}
	    }
	  //Cut on fADC integral, three vertical PMTs, and none of the six surounding PMTs. 
	  
	  //if(fadc_int[j]>5500. && fadc_int[j-12]>5500. && fadc_int[j+12]>5500. && tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r-1][c]!=0 && fadc_int[j-1]<5100 && fadc_int[j-13]<5100 && fadc_int[j+11]<5100 && fadc_int[j+1]<5100 && fadc_int[j-11]<5100 && fadc_int[j+13]<5100)
	  if(timing_method == 2)
	    {
	      //if(fadc_int[j]>fadc_int_max_cut && fadc_int[j-12]>fadc_int_max_cut && fadc_int[j+12]>fadc_int_max_cut && tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r-1][c]!=0 && fadc_int[j-1]<fadc_int_min_cut && fadc_int[j-13]<fadc_int_min_cut && fadc_int[j+11]<fadc_int_min_cut && fadc_int[j+1]<fadc_int_min_cut && fadc_int[j-11]<fadc_int_min_cut && fadc_int[j+13]<fadc_int_min_cut)
	      //3 vertical hits.
	      //if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0)
		//4 vertial hits.
	      //if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && tdc[r-2][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r-2][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0 && tdc[r-2][c+1] == 0)
		//5 vertial hits.
	      //if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && tdc[r-2][c] != 0 && tdc[r-3][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r-2][c-1] == 0 && tdc[r-3][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0 && tdc[r-2][c+1] == 0 && tdc[r-3][c+1] == 0)
	      //6 vertial hits.
	      if(tdc[r][c] !=  0 && tdc[r+1][c] != 0 && tdc[r-1][c] != 0 && tdc[r-2][c] != 0 && tdc[r-3][c] != 0 && tdc[r-4][c] != 0 && tdc[r][c-1] == 0 && tdc[r+1][c-1] == 0 && tdc[r-1][c-1] == 0 && tdc[r-2][c-1] == 0 && tdc[r-3][c-1] == 0 && tdc[r-4][c-1] == 0 && tdc[r][c+1] == 0 && tdc[r+1][c+1] == 0 && tdc[r-1][c+1] == 0 && tdc[r-2][c+1] == 0 && tdc[r-3][c+1] == 0 && tdc[r-4][c+1] == 0)
		{
		  //Fill histo with tdc times minus the reference time.
		  htdc_fadc_cuts[j]->Fill((tdc[r][c] - (tdc[r+1][c] + tdc[r-1][c])/2.)*f1_res);
		  //Fill histo with timing jitter of ref times to subtract from the final timing resolution plots in quadrature.
		  htdc_ref[j]->Fill( ((tdc[r+1][c] - tdc[r-1][c])/2.)*f1_res );
		  //cout<<"&&&&&&&&&&&&&&&&&&&& "<<(tdc[r][c] - (tdc[r+1][c]+tdc[r-1][c])/2.)*f1_res<<endl;
		}
	    }

	  //Fit the fADC signals with a Landau and find the time over threshold.
	  if(timing_method == 7)
	    {
	      if(fadc_int[j]>fadc_int_max_cut && fadc_int[j-12]>fadc_int_max_cut && fadc_int[j+12]>fadc_int_max_cut && tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r-1][c]!=0 && fadc_int[j-1]<fadc_int_min_cut && fadc_int[j-13]<fadc_int_min_cut && fadc_int[j+11]<fadc_int_min_cut && fadc_int[j+1]<fadc_int_min_cut && fadc_int[j-11]<fadc_int_min_cut && fadc_int[j+13]<fadc_int_min_cut)
		{
		  func_landau[j]->SetParameter(0,histos[r][c]->GetMaximum()*4.0);//Parameter setting the height scale (not 1:1). ~4800
		  func_landau[j]->SetParameter(1,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()));//Parameter finding the location of the peak. ~15
		  func_landau[j]->SetParameter(2,1);//Parameter setting the width. ~1
		  func_landau[j]->SetParameter(3,250);//Parameter measuring the pedestal. ~250
		  histos[r][c]->Fit(func_landau[j],"Q");//Don't use M or get some issues that make code take forever.
		  
		  //hfadc_thr_times[j]->Fill( func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))*fadc_res);
		}
	    }
	  
	}
      
      //Fill the fADC time over threshold histogram. Now that all tubes have had their fADCs fit we can use the above and below tubes for reference times.
      if(timing_method == 7)
	{
	  //Loop over PMTs
	  for(Int_t j=0; j<hcalt::ndata; j++)
	    {
	      r = hcalt::row[j]-1;
	      c = hcalt::col[j]-1;
	      if(fadc_int[j]>fadc_int_max_cut && fadc_int[j-12]>fadc_int_max_cut && fadc_int[j+12]>fadc_int_max_cut && tdc[r][c]!=0 && tdc[r+1][c]!=0 && tdc[r-1][c]!=0 && fadc_int[j-1]<fadc_int_min_cut && fadc_int[j-13]<fadc_int_min_cut && fadc_int[j+11]<fadc_int_min_cut && fadc_int[j+1]<fadc_int_min_cut && fadc_int[j-11]<fadc_int_min_cut && fadc_int[j+13]<fadc_int_min_cut)
		{
		  hfadc_thr_times[j]->Fill(    (  func_landau[j]->GetX(threshold[j],0.,histos[r][c]->GetBinCenter(histos[r][c]->GetMaximumBin()))  -  ( func_landau[j-12]->GetX(threshold[j-12],0.,histos[r-1][c]->GetBinCenter(histos[r-1][c]->GetMaximumBin()))  +  func_landau[j+12]->GetX(threshold[j+12],0.,histos[r+1][c]->GetBinCenter(histos[r+1][c]->GetMaximumBin())) )/2.0  )*fadc_res    );
		}
	    }
	}
      if(i%5000==0)
	{
	  cout<<i<<" events processed. "<<((double)i/(double)loop_max)*100.<<" % complete."<<endl;
	}
      gCurrentEntry++;
    }
  
  /*
  htiming->GetXaxis()->SetTitle("Time (ns)");
  htiming->GetXaxis()->CenterTitle(true);
  htiming->GetXaxis()->SetLabelSize(0.05);
  htiming->GetXaxis()->SetTitleSize(0.06);
  htiming->GetXaxis()->SetTitleOffset(0.75);
  htiming->GetYaxis()->SetTitle("Number of Measurements");
  htiming->GetYaxis()->CenterTitle(true);
  htiming->GetYaxis()->SetLabelSize(0.05);
  htiming->GetYaxis()->SetTitleSize(0.06);
  htiming->GetYaxis()->SetTitleOffset(0.85);
  */

  if(timing_method == 1)
    {
      TCanvas* c1=new TCanvas("c1");
      c1->SetGrid();
      //func_landau[0]->Draw();
      //histos[0][0]->Draw("same");
      histos[timerow-1][timecol-1]->Draw();
      
      //func_landau[15]->Draw("same");
      cout<<"func_landau par[0] = "<<func_landau[0]->GetParameter(0)<<", func_landau par[1] = "<<func_landau[0]->GetParameter(1)<<", func_landau par[2] = "<<func_landau[0]->GetParameter(2)<<endl;
      cout<<"Peak Max = "<<func_landau[0]->GetMaximum()<<", Peak Height = "<<func_landau[0]->GetMaximum()-func_landau[0]->GetParameter(3)<<", Peak Location = "<<func_landau[0]->GetMaximumX()<<endl;
    }

  /*
  TCanvas* ctiming=new TCanvas("ctiming");
  ctiming->SetGrid();

  htiming->Draw();
*/

  //Calculate and print timing resolutions for each channel. 
  //TF1 **func_gaus = new TF1*[nfunc];
  TF1 *func_gaus = new TF1("func_gaus",fit_gaus, min_time, max_time, 3);
  /*
  for(Int_t i=0; i<hcalt::ndata; i++)
    {
      TH1F *h = new TH1F(TString::Format("h%02d",i),TString::Format("h%02d",i),bins,min_time,max_time);
      cout<<"Gaussian fit for timing in module "<<i<<" (from 0)."<<endl;
      for(Int_t j=0; j<timing[i].size(); j++)
	{
	  h->Fill(times[i][j]);
	}
      func_gaus->SetParameter(0,h->GetMaximum());
      func_gaus->SetParameter(1,h->GetMean());
      func_gaus->SetParameter(2,h->GetStdDev());
      h->Fit("func_gaus","R 0 M");
      ind_res[i] = fabs(func_gaus->GetParameter(2));
      cout<<"The timing resolution is "<<fabs(func_gaus->GetParameter(2))<<" ns."<<endl;
      }*/

  /*
  //TDC times for each tube using ref ch that is a copy of the trigger.
  for(Int_t i=0; i<hcalt::ndata; i++)
    {
      TH1F *h = new TH1F(TString::Format("h%02d",i),TString::Format("h%02d",i),bins,min_time,max_time);
      //      cout<<"Gaussian fit for timing in module "<<i<<" (from 0)."<<endl;
      //func_gaus[i] = new TF1("func_gaus",fit_gaus, min_time, max_time, 3);
      for(Int_t j=0; j<timing[i].size(); j++)
	{
	  h->Fill(timing[i][j]);
	}
      func_gaus->SetParameter(0,h->GetMaximum());
      func_gaus->SetParameter(1,h->GetMean());
      func_gaus->SetParameter(2,h->GetStdDev());
      h->Fit("func_gaus","R 0 M Q");
      ind_res[i] = fabs(func_gaus->GetParameter(2));
      //      cout<<"The timing resolution is "<<fabs(func_gaus->GetParameter(2))<<" ns."<<endl;
    }
  */
  
  /*
  for(Int_t i=0; i<hcalt::ndata; i++)
    {
            cout<<"Timing resolution for module "<<i<<" (from 0) = "<<ind_res[i]<<" with "<<hits[i]<<" hits used for timing calculation."<<endl;
    }
  */

  /*
  //Fit the histo containing the timing data with a Gaussian fit for the channel chosen to be displayed..
  //TF1 *func_gaus = new TF1("func_gaus",fit_gaus, min_time, max_time, 3);
  cout<<"**** Single module to be displayed: Ch r"<<timerow<<"-c"<<timecol<<". ****"<<endl;
  func_gaus->SetLineColor(2);
  func_gaus->SetNpx(1000);
  func_gaus->SetParameter(0,htiming->GetMaximum());
  func_gaus->SetParameter(1,htiming->GetMean());
  func_gaus->SetParameter(2,htiming->GetStdDev());
  htiming->Fit("func_gaus","R 0 M");
  func_gaus->Draw("same");
  cout<<"The timing resolution is "<<fabs(func_gaus->GetParameter(2))<<" ns."<<endl;
  */

  /*
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  //func_landau[48]->Draw();
  c2->Divide(3,2);
  c2->cd(1);
  htimes[24]->Draw();
  c2->cd(2);
  htimes[25]->Draw();
  c2->cd(3);
  htimes[26]->Draw();
  c2->cd(4);
  htimes[27]->Draw();
  c2->cd(5);
  htimes[28]->Draw();
  c2->cd(6);
  htimes[29]->Draw();
  */

  /*
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  //func_landau[48]->Draw();
  c3->Divide(3,2);
  c3->cd(1);
  htimes_fadc[24]->Draw();
  c3->cd(2);
  htimes_fadc[25]->Draw();
  c3->cd(3);
  htimes_fadc[26]->Draw();
  c3->cd(4);
  htimes_fadc[27]->Draw();
  c3->cd(5);
  htimes_fadc[28]->Draw();
  c3->cd(6);
  htimes_fadc[29]->Draw();
*/

  /*
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();
  //func_landau[48]->Draw();
  c4->Divide(3,2);
  c4->cd(1);
  hfadc_integrals[24]->Draw();
  c4->cd(2);
  hfadc_integrals[25]->Draw();
  c4->cd(3);
  hfadc_integrals[26]->Draw();
  c4->cd(4);
  hfadc_integrals[27]->Draw();
  c4->cd(5);
  hfadc_integrals[28]->Draw();
  c4->cd(6);
  hfadc_integrals[29]->Draw();
*/

  /*
  TCanvas* c5=new TCanvas("c5");
  c5->SetGrid();
  //func_landau[48]->Draw();
  c5->Divide(3,2);
  c5->cd(1);
  htdc_fadc_cuts[24]->Draw();
  c5->cd(2);
  htdc_fadc_cuts[25]->Draw();
  c5->cd(3);
  htdc_fadc_cuts[26]->Draw();
  c5->cd(4);
  htdc_fadc_cuts[27]->Draw();
  c5->cd(5);
  htdc_fadc_cuts[28]->Draw();
  c5->cd(6);
  htdc_fadc_cuts[29]->Draw();
  */

  TCanvas *cTimes[12];

  gStyle->SetOptFit(1111);
  /*
  // Retrieve the stat box
  TPaveStats *ps = (TPaveStats*)cTimes[0]->GetPrimitive("stats");
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  ps->SetAllWith("....","size",1.04);
  */

  TF1 **func_gaus_fit = new TF1*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      func_gaus_fit[i] = new TF1("func_gaus_fit",fit_gaus, -20, 20, 3);
    }

  TF1 **func_gaus_fit_ref = new TF1*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      func_gaus_fit_ref[i] = new TF1("func_gaus_fit_ref",fit_gaus, -20, 20, 3);
    }

  Int_t x1 = -20;
  Int_t y1 = 20;
  Int_t x2 = -10;
  Int_t y2 = 0;

  Double_t sig_tdc = 0;

  TPaveText **text = new TPaveText*[channels];
  /*
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      text[i] = new TPaveText(x1,y1,x2,y2);
    }
  */

  //TPaveText *t = new TPaveText(.05,.3,.95,.6);
  //t->AddText("This line is blue"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
  //t->AddText("This line is red");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kRed);
  //t->Draw();

  if(timing_method != 4)
    {
  for(Int_t i=0; i<12; i++)
    {
      cTimes[i] = new TCanvas(Form("cTimes_Row_%d",i));
      cTimes[i]->SetGrid();
      cTimes[i]->Divide(4,3);
      /*
      cTimes[i]->cd(1);
      htdc_fadc_cuts[0]->Draw();
      cTimes[i]->cd(2);
      htdc_fadc_cuts[0]->Draw();
      cTimes[i]->cd(7);
      htdc_fadc_cuts[0]->Draw();
      */

       for(Int_t j=0;j<12;j++)
	{
	  cTimes[i]->cd(j+1);

	  if(timing_method == 0)
	    {
	      htimes[i*12+j]->Draw();
	    }
	  if(timing_method == 1)
	    {
	      htimes_fadc[i*12+j]->Draw();
	    }
	  if(timing_method == 2)
	    {
	      htdc_fadc_cuts[i*12+j]->Draw();
	      htdc_fadc_cuts[i*12+j]->SetLineColor(2);
	      func_gaus_fit[i*12+j]->SetLineColor(2);
	      func_gaus_fit[i*12+j]->SetNpx(1000);
	      func_gaus_fit[i*12+j]->SetParameter(0,htdc_fadc_cuts[i*12+j]->GetMaximum());
	      func_gaus_fit[i*12+j]->SetParameter(1,htdc_fadc_cuts[i*12+j]->GetMean());
	      func_gaus_fit[i*12+j]->SetParameter(2,htdc_fadc_cuts[i*12+j]->GetStdDev());
	      //func_gaus_fit[i*12+j]->SetParameter(0,htdc_fadc_cuts[i*12+j]->GetMaximum());
	      //func_gaus_fit[i*12+j]->SetParameter(1,htdc_fadc_cuts[i*12+j]->GetMean());
	      //func_gaus_fit[i*12+j]->SetParameter(2,htdc_fadc_cuts[i*12+j]->GetStdDev());
	      htdc_fadc_cuts[i*12+j]->Fit(func_gaus_fit[i*12+j],"Q");
	      func_gaus_fit[i*12+j]->Draw("same");
	      cout<<"Module "<<i*12+j<<" Initial Max "<<htdc_fadc_cuts[i*12+j]->GetMaximum()<<" Initial Mean "<<htdc_fadc_cuts[i*12+j]->GetMean()<<" Initial StdDev "<<htdc_fadc_cuts[i*12+j]->GetStdDev()<<" Final Max "<<func_gaus_fit[i*12+j]->GetParameter(0)<<" Final Mean "<<func_gaus_fit[i*12+j]->GetParameter(1)<<" Final StdDev "<<func_gaus_fit[i*12+j]->GetParameter(2)<<endl;

	      //Plot the reference time jitter on the same histograms.
	      htdc_ref[i*12+j]->Draw("sames");
	      func_gaus_fit_ref[i*12+j]->SetLineColor(4);
	      func_gaus_fit_ref[i*12+j]->SetNpx(1000);
	      func_gaus_fit_ref[i*12+j]->SetParameter(0,htdc_ref[i*12+j]->GetMaximum());
	      func_gaus_fit_ref[i*12+j]->SetParameter(1,htdc_ref[i*12+j]->GetMean());
	      func_gaus_fit_ref[i*12+j]->SetParameter(2,htdc_ref[i*12+j]->GetStdDev());
	      htdc_ref[i*12+j]->Fit(func_gaus_fit_ref[i*12+j],"Q + sames");
	      func_gaus_fit_ref[i*12+j]->Draw("sames");
	      cout<<"Module "<<i*12+j<<" Initial Max "<<htdc_ref[i*12+j]->GetMaximum()<<" Initial Mean "<<htdc_ref[i*12+j]->GetMean()<<" Initial StdDev "<<htdc_ref[i*12+j]->GetStdDev()<<" Final Max "<<func_gaus_fit_ref[i*12+j]->GetParameter(0)<<" Final Mean "<<func_gaus_fit_ref[i*12+j]->GetParameter(1)<<" Final StdDev "<<func_gaus_fit_ref[i*12+j]->GetParameter(2)<<endl;
	      cout<<"Standard deviation for module "<<i*12+j<<" = sqrt(|(TDC time -ref time StdDev)^2 - (TDC ref top - TDC ref bottom StdDev)^2|) = sqrt(|"<<func_gaus_fit[i*12+j]->GetParameter(2)<<"^2 - "<<func_gaus_fit_ref[i*12+j]->GetParameter(2)<<"^2|) = "<<" sqrt(|"<<pow(func_gaus_fit[i*12+j]->GetParameter(2),2.)<<" - "<<pow(func_gaus_fit_ref[i*12+j]->GetParameter(2),2.)<<"|) = "<<pow(fabs(pow(func_gaus_fit[i*12+j]->GetParameter(2),2.) - pow(func_gaus_fit_ref[i*12+j]->GetParameter(2),2.)),0.5)<<" ns."<<endl;

	      sig_tdc = pow(fabs(pow(func_gaus_fit[i*12+j]->GetParameter(2),2.) - pow(func_gaus_fit_ref[i*12+j]->GetParameter(2),2.)),0.5);
	      x1 = -DISP_MAX_SAMPLE;
	      y1 = func_gaus_fit_ref[i*12+j]->GetMaximum(-DISP_MAX_SAMPLE,DISP_MAX_SAMPLE)*0.7;
	      x2 = -DISP_MAX_SAMPLE*.25;
	      y2 = func_gaus_fit_ref[i*12+j]->GetMaximum(-DISP_MAX_SAMPLE,DISP_MAX_SAMPLE)*0.5;
	      text[i*12+j] = new TPaveText(x1,y1,x2,y2);
	      text[i*12+j]->AddText("TDC StdDev: ");  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->AddText(Form("%0.3f ns",sig_tdc));  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	      text[i*12+j]->Draw("same");
	    }
	  if(timing_method == 3)
	    {
	      gPad->SetLogy();
	      hfadc_integrals[i*12+j]->Draw();
	    }
	  if(timing_method == 5)
	    {
	      htdc_time_vs_fadc_int[i*12+j]->Draw();
	    }
	  if(timing_method == 6)
	    {
	      htdc_verticallity[i*12+j]->Draw();
	    }
	  if(timing_method == 7)
	    {
	      hfadc_thr_times[i*12+j]->Draw();
	      func_gaus_fit[i*12+j]->SetLineColor(2);
	      func_gaus_fit[i*12+j]->SetNpx(1000);
	      func_gaus_fit[i*12+j]->SetParameter(0,hfadc_thr_times[i*12+j]->GetMaximum());
	      func_gaus_fit[i*12+j]->SetParameter(1,hfadc_thr_times[i*12+j]->GetMean());
	      func_gaus_fit[i*12+j]->SetParameter(2,hfadc_thr_times[i*12+j]->GetStdDev());
	      hfadc_thr_times[i*12+j]->Fit(func_gaus_fit[i*12+j],"Q");
	      func_gaus_fit[i*12+j]->Draw("same");
	      cout<<"Module "<<i*12+j<<" Initial Max "<<hfadc_thr_times[i*12+j]->GetMaximum()<<" Initial Mean "<<hfadc_thr_times[i*12+j]->GetMean()<<" Initial StdDev "<<hfadc_thr_times[i*12+j]->GetStdDev()<<" Final Max "<<func_gaus_fit[i*12+j]->GetParameter(0)<<" Final Mean "<<func_gaus_fit[i*12+j]->GetParameter(1)<<" Final StdDev "<<func_gaus_fit[i*12+j]->GetParameter(2)<<endl;
	    }
	  //htdc_fadc_cuts[0]->Draw();
	  //	  htdc_fadc_cuts[i*12+j]->Draw();    //
	  //       	  htimes_fadc[i*12+j]->Draw();       //
       	  //hfadc_integrals[i*12+j]->Draw();   //
	  //	  htimes[i*12+j]->Draw();            //
	}

       if(timing_method == 0)
	 {
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/TDC_Time_Trig/TDC_Time_Run%d_Row%d.png",run,i));
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/TDC_Time_Trig/TDC_Time_Run%d_Row%d.C",run,i));
	 }
       if(timing_method == 1)
	 {
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/fADC_Time/fADC_Time_Run%d_Row%d.png",run,i));
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/fADC_Time/fADC_Time_Run%d_Row%d.C",run,i));
	 }
       if(timing_method == 2)
	 {
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/TDC_Time_fADC_Cuts/TDC_Time_fADC_Cuts_Run%d_Row%d.png",run,i));
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/TDC_Time_fADC_Cuts/TDC_Time_fADC_Cuts_Run%d_Row%d.C",run,i));
	 }
       if(timing_method == 3)
	 {
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/fADC_Integrals/fADC_Integrals_Run%d_Row%d.png",run,i));
	   cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/fADC_Integrals/fADC_Integrals_Run%d_Row%d.C",run,i));
	 }
    }
    }

  /*
  TCanvas* cint=new TCanvas("cint");
  cint->SetGrid();

  TCanvas* ctdc=new TCanvas("ctdc");
  ctdc->SetGrid();

  TCanvas* cint_tdc_cut=new TCanvas("cint_tdc_cut");
  cint_tdc_cut->SetGrid();

  if(timing_method == 4)
    {
      hfadc_integrals_tdc_cuts[0]->SetLineColor(2);
      cint->cd();
      cout<<"fADC integrals with no cut = "<<hfadc_integrals[0]->GetEntries()<<endl;
      hfadc_integrals[0]->Draw();
      hfadc_integrals_tdc_cuts[0]->Draw("same");
      ctdc->cd();
      htdc_times[0]->Draw();
      cint_tdc_cut->cd();
      cout<<"fADC integrals with TDC cut = "<<hfadc_integrals_tdc_cuts[0]->GetEntries()<<endl;
      hfadc_integrals_tdc_cuts[0]->Draw();
    }
  */

  //Print the pedestal values.
  /*
  for(Int_t i=0; i<hcalt::ndata; i++)
    {
      pedestals[i] = pedestals[i]/no_tdc[i];
      cout<<"Pedestals["<<i<<"] = "<<pedestals[i]<<endl;
      pedestal[i] = pedestal[i];
      cout<<"Pedestal["<<i<<"] = "<<pedestal[i]<<endl;
    }
  */

  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
