#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "hcal.h"
#include <vector>
#include <TStopwatch.h>
#include <TDatime.h>
using namespace std;

Double_t entry;
Int_t gCurrentEntry = 0;
Int_t DISP_MIN_SAMPLE = 0;
Int_t DISP_MAX_SAMPLE = 20;
Int_t update_pedestals_file = 1;          //DOesn't work! 0-> Doesn't update the Pedestals.txt file with new values. 1-> Updates the Pedestals.txt file with the newly calculated values.
Int_t limit_evts = 0;                     //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t loop_max = 0;                       //Dummy variable set equal to max_evts or nevt for the loop.
Int_t max_evts = 5000;                   //Maximum number of events to analyze if limit_evts = 1.
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
Int_t min_time = -300;
Int_t max_time = 50;
Int_t bins = 250;
Int_t r,c,n,idx;
Int_t nevt;
const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
Int_t ref_ch = 143;
Int_t ref_row = 11;
Int_t ref_col = 11;
Int_t ref_ch_2 = 142;               //Second reference channel if want to see relative timing resolutions between same triggers.
Int_t ref_row_2 = 11;
Int_t ref_col_2 = 10;
Int_t ped_correlation_1 = 0;        //Modules (from zero) to check pedestal correlations.
Int_t ped_correlation_2 = 6;
Double_t ped_correlation_1_fadc_int = 0;          //fADC integrals for the modules being checked for pedestal correlations.
Double_t ped_correlation_2_fadc_int = 0;
Double_t ind_res[144];
TChain *T = 0;
std::string user_input;

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
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

void Pedestals(Int_t run = 501)
//Int_t Test_Display(Int_t run = 290)
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Create a date object.
  TDatime time;

  //Open file to write output to.
  std::ofstream output1 (Form("./Pedestals.txt"), std::ofstream::out);
  if(update_pedestals_file == 1)
    {
      output1<<"These pedestals were generated from run "<<run<<"."<<endl;
      output1<<"They were generated at "<<time.GetHour()<<":"<<time.GetMinute()<<" on "<<time.GetMonth()<<"/"<<time.GetDay()<<"/"<<time.GetYear()<<"."<<endl;
      output1<<"Module,"<<" Avg Pedestal Value"<<" Avg Peak Height"<<endl;
    }

  if(!T) 
    { 
      T = new TChain("T");
      //T->Add(TString::Format("rootfiles/fadc_f1tdc_%d.root",run));//325
 
      //============  Reading the Rootfile =======================//
      
      const TString rootfilePath = "/home/skbarcus/JLab/SBS/HCal/Analysis/Cosmics/rootfiles/";
      std::ostringstream str;
      str << rootfilePath<<"fadc_f1tdc_"<<run;
      TString basename = str.str().c_str();
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

  nevt = T->GetEntries();
  Double_t adc[kNrows][kNcols];
  Double_t tdc[kNrows][kNcols];
  Double_t pedestals[channels] = {};
  Int_t no_tdc[channels] = {};
  Int_t yes_tdc[channels] = {};
  Double_t fadc_int[channels] = {};
  Double_t peak_height[channels] = {};

  //Histograms to store maximum amplitudes for fADC hits.
  TH1F **hamplitudes = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hamplitudes[i] = new TH1F(Form("hamplitudes%d",i),Form("Maximum Amplitudes for Module %d",i),bins,0,2000);
    }

  //Histograms for showing the fADC integrals for each pedestal event (i.e. not TDC hit).
  TH1F **hped_int = new TH1F*[channels];
  for(Int_t i = 0; i<channels; i++)
    {
      hped_int[i] = new TH1F(Form("hped_int%d",i),Form("fADC Integrals for Pedestal Events for Module %d",i),bins,0,7000);
    }

  TH2F *hped_correlation = new TH2F("hped_correlation",Form("Corelation Between Module %d and %d (from 0) Pedestals",ped_correlation_1,ped_correlation_2),1000,0,7000,1000,0,7000);


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
      //T->GetEntry(11);
      //T->GetEntry(665);
      
      //Set arrays to zero;
      for(r  = 0; r < kNrows; r++) 
	{
	  for(c  = 0; c < kNcols; c++) 
	    {
	      adc[r][c] = 0.0;
	      tdc[r][c] = 0.0;
	    }
	}

      //Loop over all channels and fill the arrays.
      ped_correlation_1_fadc_int = 0;
      ped_correlation_2_fadc_int = 0;
      for(Int_t j=0; j<hcalt::ndata; j++)
	{
	  r = hcalt::row[j]-1;
	  c = hcalt::col[j]-1;
	  idx = hcalt::samps_idx[j];
	  n = hcalt::nsamps[j];
	  adc[r][c] = hcalt::a[j];
	  tdc[r][c] = hcalt::tdc[j];
	  Double_t ped_val = 0;
	  fadc_int[j] = 0;             //Reset fadc integrals.
	  
	  //Find the pedestal values when the TDC didn't fire (should be no fADC signal either). 
	  if(tdc[r][c] == 0 || j==ref_ch || j==ref_ch_2)// || j==36 || j==96)
	    {
	      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
		{
		  histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
		  fadc_int[j] = fadc_int[j] + histos[r][c]->GetBinContent(s+1-DISP_MIN_SAMPLE);
		  //ped_val = ped_val + histos[r][c]->GetBinContent(s+1);
		}
	      //cout<<"Event "<<i<<" PMT "<<j<<": fadc_int = "<<fadc_int[j]<<endl;

	      ped_val = fadc_int[j]/DISP_FADC_SAMPLES;
	      pedestals[j] = pedestals[j] + ped_val;
	      no_tdc[j] = no_tdc[j]+1;
	      //Fill the histogram of fADC pedestal event integrals.
	      hped_int[j]->Fill(fadc_int[j]);
	      //Fill the histogram for the two modules to check their pedestal correlation.
	      if(j == ped_correlation_1)
		{
		  ped_correlation_1_fadc_int = fadc_int[j];
		}
	      if(j == ped_correlation_2)
		{
		  ped_correlation_2_fadc_int = fadc_int[j];
		}
	    }

	  //Find average fadc integral peak height for time over threshold calculations.
	  if(tdc[r][c] != 0)
	    {
	      //Make fadc histo for each event and each pmt.
	      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
		{
		  histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE,hcalt::samps[idx+s]);
		  //peak_height[j] = histos[r][c]->GetMaximum() + peak_height[j];
		  //yes_tdc[j] = yes_tdc[j]+1;
		}
	      //Find and store the maximum peak height in the histo.
	      hamplitudes[j]->Fill(histos[r][c]->GetMaximum());
	      peak_height[j] = histos[r][c]->GetMaximum() + peak_height[j];
	      yes_tdc[j] = yes_tdc[j]+1;
	    }
	}
      if(i%5000==0)
	{
	  cout<<i<<" events processed. "<<((double)i/(double)loop_max)*100.<<" % complete."<<endl;
	}

      //Fill the pedestal correlation histogram.
      hped_correlation->Fill(ped_correlation_1_fadc_int,ped_correlation_2_fadc_int);

      gCurrentEntry++;
    }

  TCanvas *cAmplitudes[12];

  gStyle->SetOptFit(1111);

  Int_t fit_max = 1000;
  TF1 **func_gaus_fit = new TF1*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
    {
      func_gaus_fit[i] = new TF1("func_gaus_fit",fit_gaus, hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin())+50, fit_max, 3);
      //cout<<"Location of peak "<<hamplitudes[i]->GetMaximumBin()<<" or "<<hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin())<<" or "<<hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximum())<<endl;
      //cout<<hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin(hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin())+50,1500))<<endl;
      //cout<<hamplitudes[i]->GetMaximumBin(hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin())+50,1500,"")<<endl;
      //Int_t x = 300;
      //Int_t y = 1500;
      //Int_t z = 10000;
      //cout<<hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin(x,y,z))<<endl;
      /*
      cout<<"Module "<<i<<": Bin Range ("<<hamplitudes[i]->GetMaximumBin()+5<<","<<bins<<") = ("<<hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin()+2)<<","<<hamplitudes[i]->GetBinCenter(bins)<<")"<<endl;
      hamplitudes[i]->GetXaxis()->SetRange(hamplitudes[i]->GetMaximumBin()+5,bins);
      cout<<hamplitudes[i]->GetBinCenter(hamplitudes[i]->GetMaximumBin())<<endl;
      hamplitudes[i]->GetXaxis()->SetRange(0,bins);
      */
    }

  Double_t ped_peak = 0;
  Double_t spe_peak = 0;
  Double_t spe_bin = 0;
  Double_t chi1 = 0;
  Double_t chi2 = 0;

  //Make plots of fADC amplitudes.
  for(Int_t i=0; i<12; i++)
    {
      cAmplitudes[i] = new TCanvas(Form("cAmplitudes_%d",i));
      cAmplitudes[i]->SetGrid();
      cAmplitudes[i]->Divide(4,3);
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
	  spe_bin = 0;
	  spe_peak = 0;
	  ped_peak = 0;
	  chi1 = 0;
	  chi2 = 0;
	  fit_max = 1000;

	  //Find pedestal location.
	  ped_peak = hamplitudes[i*12+j]->GetMaximumBin();
	  //Restrict range to above pedestal.
	  hamplitudes[i*12+j]->GetXaxis()->SetRange(ped_peak+5,bins);
	  //Find SPE peak bin above the pedestal.
	  spe_bin = hamplitudes[i*12+j]->GetMaximumBin();
	  //Find SPE peak value above the pedestal.
	  spe_peak = hamplitudes[i*12+j]->GetBinCenter(spe_bin);
	  //Return range to full.
	  hamplitudes[i*12+j]->GetXaxis()->SetRange(0,bins);

	  cAmplitudes[i]->cd(j+1);
	  hamplitudes[i*12+j]->Draw();
	  func_gaus_fit[i*12+j]->SetLineColor(2);
	  func_gaus_fit[i*12+j]->SetNpx(1000);
	  //func_gaus_fit[i*12+j]->SetParameter(0,hamplitudes[i*12+j]->GetMaximum());
	  //func_gaus_fit[i*12+j]->SetParameter(0,50);
	  func_gaus_fit[i*12+j]->SetParameter(0,hamplitudes[i*12+j]->GetBinContent(spe_bin));
	  //cout<<hamplitudes[i*12+j]->GetBinContent(spe_bin)<<endl;
	  func_gaus_fit[i*12+j]->SetParLimits(0,hamplitudes[i*12+j]->GetBinContent(spe_bin)*0.9,hamplitudes[i*12+j]->GetBinContent(spe_bin)*1.1);
	  //func_gaus_fit[i*12+j]->SetParameter(1,hamplitudes[i*12+j]->GetMean());
	  //func_gaus_fit[i*12+j]->SetParameter(1,500);
	  func_gaus_fit[i*12+j]->SetParameter(1,spe_peak);
	  func_gaus_fit[i*12+j]->SetParLimits(1,spe_peak-100,spe_peak+100);
	  //func_gaus_fit[i*12+j]->SetParLimits(1,300,900);
	  //func_gaus_fit[i*12+j]->SetParameter(2,hamplitudes[i*12+j]->GetStdDev());
	  func_gaus_fit[i*12+j]->SetParameter(2,400);
	  //func_gaus_fit[i*12+j]->SetParameter(0,htdc_fadc_cuts[i*12+j]->GetMaximum());
	  //func_gaus_fit[i*12+j]->SetParameter(1,htdc_fadc_cuts[i*12+j]->GetMean());
	  //func_gaus_fit[i*12+j]->SetParameter(2,htdc_fadc_cuts[i*12+j]->GetStdDev());
	  hamplitudes[i*12+j]->Fit(func_gaus_fit[i*12+j],"QR");

	  //Try to improve Gaussian fitting range by decreaing the upper fit limit until chi^2 is minimized.
	  //Calculate initial chi^2.
	  chi1 = func_gaus_fit[i*12+j]->GetChisquare()/func_gaus_fit[i*12+j]->GetNDF();
	  //Reset fit maximum by 50 units.
	  fit_max = fit_max - 50;
	  func_gaus_fit[i*12+j]->SetRange(hamplitudes[i*12+j]->GetBinCenter(ped_peak)+50,fit_max);
	  //Redo the fit with the new range.
	  hamplitudes[i*12+j]->Fit(func_gaus_fit[i*12+j],"QR");
	  //Calculate the updated chi^2.
	  chi2 = func_gaus_fit[i*12+j]->GetChisquare()/func_gaus_fit[i*12+j]->GetNDF();
	  //While the new fit has a "better" (lower reduced chi^2) than the original keep decreasing the fit range maximum.
	  for(;chi2<chi1;)
	    {
	      chi1 = chi2;
	      fit_max = fit_max - 50;
	      //If we get to close to passing over the SPE peak stop decreasing the max fit range.
	      if(fit_max < spe_peak + 50)
		{
		  break;
		}
	      func_gaus_fit[i*12+j]->SetRange(hamplitudes[i*12+j]->GetBinCenter(ped_peak)+50,fit_max);
	      hamplitudes[i*12+j]->Fit(func_gaus_fit[i*12+j],"QR");
	      chi2 = func_gaus_fit[i*12+j]->GetChisquare()/func_gaus_fit[i*12+j]->GetNDF();
	    }
	  //Go up in fit range back to where the ch1^2 was best and recalculate the fits.
	  fit_max = fit_max + 50;
	  func_gaus_fit[i*12+j]->SetRange(hamplitudes[i*12+j]->GetBinCenter(ped_peak)+50,fit_max);
	  hamplitudes[i*12+j]->Fit(func_gaus_fit[i*12+j],"QR");

	  func_gaus_fit[i*12+j]->Draw("same");
	  //cout<<"Module "<<i*12+j<<" Initial Max "<<hamplitudes[i*12+j]->GetMaximum()<<" Initial Mean "<<hamplitudes[i*12+j]->GetMean()<<" Initial StdDev "<<hamplitudes[i*12+j]->GetStdDev()<<" Final Max "<<func_gaus_fit[i*12+j]->GetParameter(0)<<" Final Mean "<<func_gaus_fit[i*12+j]->GetParameter(1)<<" Final StdDev "<<func_gaus_fit[i*12+j]->GetParameter(2)<<endl;
	  cout<<"Module "<<i*12+j<<": SPE Peak Counts = "<<func_gaus_fit[i*12+j]->GetParameter(0)<<", Mean fADC Amplitude = "<<func_gaus_fit[i*12+j]->GetParameter(1)<<", StdDev fADC Amplitude = "<<func_gaus_fit[i*12+j]->GetParameter(2)<<", Avg NPE = "<<pow(func_gaus_fit[i*12+j]->GetParameter(1)/func_gaus_fit[i*12+j]->GetParameter(2),2.)<<" with a rChi^2 = "<<func_gaus_fit[i*12+j]->GetChisquare()/func_gaus_fit[i*12+j]->GetNDF()<<"."<<endl;
	}
    }

  //Plot the pedestal events fADC integrals for each module. 
  TCanvas *cped_int[12];

  Double_t ped_int_std_dev[channels] = {};
  Double_t ped_int_mean[channels] = {};
  Double_t ratio[channels] = {};
  Double_t JLab_pmt_std_dev_avg = 0;
  Double_t JLab_pmt_mean_avg = 0;
  Double_t CMU_pmt_std_dev_avg = 0;
  Double_t CMU_pmt_mean_avg = 0;

  for(Int_t i=0; i<12; i++)
    {
      cped_int[i] = new TCanvas(Form("cped_int_%d",i));
      cped_int[i]->SetGrid();
      cped_int[i]->Divide(4,3);

       for(Int_t j=0;j<12;j++)
	{
	  cped_int[i]->cd(j+1);
	  hped_int[i*12+j]->Draw();
	  //Find the avg ratio of std dev to mean for each tube.
	  ped_int_std_dev[i*12+j] = hped_int[i*12+j]->GetStdDev();
	  ped_int_mean[i*12+j] = hped_int[i*12+j]->GetMean();
	  ratio[i*12+j] = ped_int_std_dev[i*12+j]/ped_int_mean[i*12+j];
	  if(j>3 && j<8)
	    {
	      JLab_pmt_std_dev_avg = JLab_pmt_std_dev_avg + hped_int[i*12+j]->GetStdDev();
	      JLab_pmt_mean_avg = JLab_pmt_mean_avg + hped_int[i*12+j]->GetMean();
	    }
	  else
	    {
	      CMU_pmt_std_dev_avg = CMU_pmt_std_dev_avg + hped_int[i*12+j]->GetStdDev();
	      CMU_pmt_mean_avg = CMU_pmt_mean_avg + hped_int[i*12+j]->GetMean();
	    }
	  cout<<"Standard deviation of the mean pedestal event fADC integral as a percent of the mean pedestal event fADC integal value for PMT "<<i*12+j<<" = "<<ratio[i*12+j]*100.<<"%."<<endl;
	}
    }
  cout<<"JLab PMTs average mean pedestal event fADC integral = "<<JLab_pmt_mean_avg/(channels/3)<<"."<<endl;
  cout<<"JLab PMTs average standard deviation of the mean pedestal event fADC integral = "<<JLab_pmt_std_dev_avg/(channels/3)<<"."<<endl;
  cout<<"Standard deviation of the mean pedestal event fADC integral as a percent of the mean pedestal event fADC integal value for JLab PMTs = "<<(JLab_pmt_std_dev_avg/(channels/3))/(JLab_pmt_mean_avg/(channels/3))*100.<<"%."<<endl;
  cout<<"CMU PMTs average mean pedestal event fADC integral = "<<CMU_pmt_mean_avg/(2*channels/3)<<"."<<endl;
  cout<<"CMU PMTs average standard deviation of the mean pedestal event fADC integral = "<<CMU_pmt_std_dev_avg/(2*channels/3)<<"."<<endl;
  cout<<"Standard deviation of the mean pedestal event fADC integral as a percent of the mean pedestal event fADC integal value for CMU PMTs = "<<(CMU_pmt_std_dev_avg/(2*channels/3))/(CMU_pmt_mean_avg/(2*channels/3))*100.<<"%."<<endl;

  TCanvas* cped_correlations = new TCanvas("cped_correlations");
  cped_correlations->SetGrid();
  hped_correlation->Draw();
  
  //Print the pedestal values and average peak height values.
  for(Int_t i=0; i<hcalt::ndata; i++)
    {
      cout<<"Module "<<i<<" pedestal[i] = "<<pedestals[i]<<" no_tdc[i] = "<<no_tdc[i]<<endl;
      pedestals[i] = pedestals[i]/no_tdc[i];
      peak_height[i] = peak_height[i]/yes_tdc[i];
      cout<<"Module "<<i<<" has an average pedestal of "<<pedestals[i]<<" and an average peak height of "<<peak_height[i]<<endl;
      if(update_pedestals_file == 1)
	{
	  output1<<i<<"   "<<pedestals[i]<<"   "<<peak_height[i]<<endl;
	}
    }

  if(update_pedestals_file == 1)
    {
      output1.close();
    }

  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
