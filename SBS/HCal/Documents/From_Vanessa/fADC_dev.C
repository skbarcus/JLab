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


// compile me by setting COMPILE to 1, then
// g++ `root-config --libs --cflags` fADC.C -o fADC.exe -Wall
// you can execute me via
// ./fADC.exe runnumber
// where runnumber is the run number to analyze (see variable "run" for current default)

//#define COMPILE 1
// otherwise (interactive mode) just comment the line above

// you might want to change here...
Int_t limit_evts = 0;                     //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t max_evts = 10000;                   //Maximum number of events to analyze if limit_evts = 1.

// ...and here
// path Vanessa
const TString pedestal_file = "/Users/vanessabrio/Pedestals_run820.txt";
const TString rootfilePath = "/Users/vanessabrio/rootfiles/";

Int_t loop_max = 0;                       //Dummy variable set equal to max_evts or nevt for the loop.
Double_t entry;
Int_t gCurrentEntry = 0;
Int_t DISP_MIN_SAMPLE = 0;
Int_t DISP_MAX_SAMPLE = 20;
const Int_t DISP_FADC_SAMPLES = (DISP_MAX_SAMPLE-DISP_MIN_SAMPLE);
Int_t min_time = -300;
Int_t max_time = 50;
Int_t bins = 100;
Double_t fadc_int_max_cut = 0;
Double_t fadc_int_min_cut = 100;
Int_t r,c,n,idx;
Int_t nevt;
Int_t refrow = 1000;     //Counting from 1. 
Int_t refcol = 1000;     //Counting from 1.
Int_t timerow = 1;    //Counting from 1.
Int_t timecol = 9;    //Counting from 1.
const Int_t kNrows = 12;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;
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
Double_t rms[kNrows][kNcols];

TH1F* MakeHisto(Int_t row, Int_t col, Int_t bins)
{
  TH1F *h = new TH1F(TString::Format("h%02d%02d",row,col),
      TString::Format("%d-%d",row+1,col+1),bins,DISP_MIN_SAMPLE,DISP_MAX_SAMPLE);
  h->SetStats(0);
  h->SetLineWidth(2);
  return h;
}

//Create landau to fit the timing resolution.
Double_t fit_landau(Double_t *X, Double_t *par) 
{
  Double_t landau = par[0]*TMath::Landau(X[0],par[1],par[2])+par[3];
  return landau;
}

Double_t fit_exp(Double_t *x, Double_t *par)
{
  Double_t expo = par[0] + par[1] * TMath::Exp(- x[0] / par[2]);
  return expo;
}

#ifndef COMPILE
void Timing_Resolution(Int_t run = 820)
{
#else
int main(int argc, char **argv)
{
  Int_t run = 820;
  if (argc == 2){
    run = atoi(argv[1]);
    cout << "run "<< run << endl;
  }
#endif
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  FILE *fp;
  fp = fopen(pedestal_file,"r");
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

    //============  Reading the Rootfile =======================//

    std::ostringstream str;
    str << rootfilePath<<"fadc_f1tdc_"<<run;
    TString basename = str.str().c_str();
    TString rootfile = basename + ".root";
    cout<<basename<<endl;

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
#ifndef COMPILE
      return;
#else
      return -1;
#endif
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

  cout<<"Hello"<<endl;
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

  TF1 **func_exp = new TF1*[2];
  for(Int_t i = 142; i<hcalt::ndata; i++)
  {
    func_exp[i-142] = new TF1("func_expo",fit_exp, DISP_MIN_SAMPLE, DISP_MAX_SAMPLE, 3);
  }

  TH1F **htimes_fadc = new TH1F*[channels];
  TH1F **histos_cut  = new TH1F*[channels];
  TH1F **htimes_at_threshold = new TH1F*[channels];
  for(Int_t i = 0; i<hcalt::ndata; i++)
  {
    htimes_fadc[i] = new TH1F(Form("htimes_fadc%d",i),Form("fADC Timing for Module %d",i),bins,DISP_MIN_SAMPLE+10,DISP_MAX_SAMPLE*fadc_res);
    histos_cut[i] = new TH1F(Form("htimes_fadc_cut%d",i),Form("fADC Timing for Module %d cut",i),bins,DISP_MIN_SAMPLE+10,DISP_MAX_SAMPLE*fadc_res);
    histos_cut[i]->SetLineColor(kRed);
    htimes_at_threshold[i] = new TH1F(Form("htimes_threshold%d",i),Form("fADC Time at threshold for Module %d",i),bins, 0, DISP_MAX_SAMPLE*fadc_res/2);
  }

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

  /////////////////////////////////////////////////////////
  // test samps
  TCanvas *test = new TCanvas("test","test samples",1400,1000);
  test->Divide(12,12);
  test->Update();
  for (int entrynum = 0; entrynum<loop_max ; entrynum++)
  {
    T->GetEntry(entrynum);
    memset(peak,     0, channels * sizeof(Int_t));
    memset(adc,      0, channels * sizeof(Int_t));
    memset(tdc,      0, channels * sizeof(Int_t));
    Double_t ped_val = 0;
    Double_t time_at_threshold = 0;

    for(Int_t j=0; j<hcalt::ndata; j++)
    {
      r = hcalt::row[j]-1;
      c = hcalt::col[j]-1;
      idx = hcalt::samps_idx[j];
      n = hcalt::nsamps[j];
      adc[r][c] = hcalt::a[j];
      Bool_t at_threshold = false;
      for(Int_t s = DISP_MIN_SAMPLE; s < DISP_MAX_SAMPLE && s < n; s++) 
      {
        Double_t value = hcalt::samps[idx+s];
        histos[r][c]->SetBinContent(s+1-DISP_MIN_SAMPLE, value);
        ped_val = ped_val + histos[r][c]->GetBinContent(s+1);
        if ( at_threshold == false && value >= threshold[j] ) {
          at_threshold = true;
          time_at_threshold = (s-DISP_MIN_SAMPLE)*f1_res ;
          htimes_at_threshold[j]->Fill( time_at_threshold );
        }
      }

      test->cd(12*(r)+c+1);
      histos[r][c]->Draw();

      func_landau[j]->SetParameter(0,histos[r][c]->GetMaximum()*4.0);//Parameter setting the height scale (not 1:1). ~4800
      func_landau[j]->SetParameter(1,15);//Parameter finding the location of the peak. ~15
      func_landau[j]->SetParameter(2,1);//Parameter setting the width. ~1
      func_landau[j]->SetParameter(3,250);//Parameter measuring the pedestal. ~250
      histos[r][c]->Fit(func_landau[j],"Q");//Don't use M or get some issues that make code take forever.
      rms[r][c] = func_landau[j]->GetParameter(2);

      //tdc[r][c] = hcalt::tdc[j];

      //cout<<"** event "<< entrynum<<" row "<< r+1 <<" col "<< c+1 <<" RMS: "<< rms[r][c] << " bin at threshold: "<< time_at_threshold / f1_res << " time: "<< time_at_threshold <<endl;

    }
    test->Update();
 //   test->SaveAs(Form("./Timing_Resolution_Plots/fADC_Time/samples_event_%d.png",entrynum));
  }

  Double_t time_at_threshold_rms_ref = (htimes_at_threshold[hcalt::ndata -2]->GetRMS() + htimes_at_threshold[hcalt::ndata -1]->GetRMS())/2.;


  TCanvas *ctimes_at_threshold = new TCanvas("test1","times at threshold",1400,1000);
  ctimes_at_threshold->Divide(12,12);
  cout <<endl;
  for (int chan=0; chan<hcalt::ndata; chan++){
    Double_t time_at_threshold_rms = htimes_at_threshold[chan]->GetRMS();
    cout << "Threshold[ "<<chan<< " ]: "<< threshold[chan]
      << " Entries: "<< htimes_at_threshold[chan]->GetEntries()
      << " Mean: "<< htimes_at_threshold[chan]->GetMean()
      <<" RMS: "<< htimes_at_threshold[chan]->GetRMS()
      <<" time res (ns): " << 4*sqrt( abs(pow(time_at_threshold_rms,2) - pow(time_at_threshold_rms_ref,2)) )
      << endl;
    ctimes_at_threshold->cd( chan+1 );
    ctimes_at_threshold->SetLogy();
    htimes_at_threshold[chan]->Draw();
  }
  ctimes_at_threshold->SaveAs(Form("./Timing_Resolution_Plots/fADC_Time/time_res_run_%d.png",run));
  /////////////////////////////////////////////////////////
#ifndef COMPILE
  return; // just quit now, no further analysis
#else
  return 0;
#endif



  //for(Int_t i=0; i<1 ;i++)
  for(Int_t i=0; i<loop_max ;i++)
  {
    T->GetEntry(gCurrentEntry);

    //Set arrays to zero;
    memset(peak,     0, channels * sizeof(Int_t));
    memset(adc,      0, channels * sizeof(Int_t));
    memset(tdc,      0, channels * sizeof(Int_t));

    Double_t fadc_int[channels];
    memset(fadc_int, 0, channels * sizeof(Int_t));
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

      func_landau[j]->SetParameter(0,histos[r][c]->GetMaximum()*4.0);//Parameter setting the height scale (not 1:1). ~4800
      func_landau[j]->SetParameter(1,15);//Parameter finding the location of the peak. ~15
      func_landau[j]->SetParameter(2,1);//Parameter setting the width. ~1
      func_landau[j]->SetParameter(3,250);//Parameter measuring the pedestal. ~250

      if(tdc[r][c] !=  0)   //Remove this TDC firing requirement for mode 4. 
      {
        //	      timing[j].push_back((tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res);
        //	      times[i][j] = (tdc[refrow-1][refcol-1] - tdc[r][c])*f1_res;
        if((r != (refrow - 1) || c != (refcol - 1)))
        {
          histos[r][c]->Fit(func_landau[j],"Q");//Don't use M or get some issues that make code take forever.
          rms[r][c] = func_landau[j]->GetParameter(2);

          if(tdc[r+1][c] != 0 && tdc[r-1][c] != 0)
          {
            htimes_fadc[j]->Fill((func_landau[j]->GetMaximumX())*fadc_res);
            hits[j] = hits[j] + 1;
          }
        }
      }
    }

    // apply cuts
    // normal case: module not in first/last row/column
    for(Int_t r=1; r<10; r++)
    {
      for (Int_t c=1; c<10; c++)
      {
        if (htimes_fadc[r*12+c]->GetEntries() > fadc_int_max_cut && 
            // values in same column
            htimes_fadc[(r-1)*12+c]->GetEntries() > fadc_int_max_cut && 
            htimes_fadc[(r+1)*12+c]->GetEntries() > fadc_int_max_cut &&
            // values in adjacent columns
            htimes_fadc[(r-1)*12+c-1]->GetEntries() < fadc_int_min_cut &&
            htimes_fadc[(r-1)*12+c+1]->GetEntries() < fadc_int_min_cut &&
            htimes_fadc[r  *12+c-1]->GetEntries() < fadc_int_min_cut &&
            htimes_fadc[r  *12+c+1]->GetEntries() < fadc_int_min_cut &&
            htimes_fadc[(r+1)*12+c-1]->GetEntries() < fadc_int_min_cut &&
            htimes_fadc[(r+1)*12+c+1]->GetEntries() < fadc_int_min_cut)
        {
          for (Int_t i=0; i<channels; i++)
          {
            Double_t oldcounts = histos_cut[r*12+c]->GetBinContent(i);
            histos_cut[r*12+c]->SetBinContent(i, oldcounts + htimes_fadc[r*12+c]->GetBinContent(i));
            //cout << "XXX module "<< r*12+c <<" bin "<< i <<" oldcounts "<< oldcounts <<" htimes_fadc "<< htimes_fadc[r*12+c]->GetBinContent(i) <<" newcounts "<< histos_cut[r*12+c]->GetBinContent(i) <<endl;
          }
        }
      }
    }
    // module in first/last row, except corners
    for (Int_t c=1; c<10; c++)
    {
      if (htimes_fadc[0*12+c]->GetEntries() > fadc_int_max_cut && 
          htimes_fadc[1*12+c]->GetEntries() > fadc_int_max_cut &&
          htimes_fadc[0*12+c-1]->GetEntries() < fadc_int_min_cut &&
          htimes_fadc[0*12+c+1]->GetEntries() < fadc_int_min_cut &&
          htimes_fadc[1*12+c-1]->GetEntries() < fadc_int_min_cut &&
          htimes_fadc[1*12+c+1]->GetEntries() < fadc_int_min_cut)
      {
        for (Int_t i=0; i<channels; i++)
        {
          Double_t oldcounts = histos_cut[0*12+c]->GetBinContent(i);
          histos_cut[0*12+c]->SetBinContent(i, oldcounts + htimes_fadc[0*12+c]->GetBinContent(i));
        }
      }
      if (htimes_fadc[11*12+c]->GetEntries() > fadc_int_max_cut && 
          htimes_fadc[10*12+c]->GetEntries() > fadc_int_max_cut &&
          htimes_fadc[11*12+c-1]->GetEntries() < fadc_int_min_cut &&
          htimes_fadc[11*12+c+1]->GetEntries() < fadc_int_min_cut &&
          htimes_fadc[10*12+c-1]->GetEntries() < fadc_int_min_cut &&
          htimes_fadc[10*12+c+1]->GetEntries() < fadc_int_min_cut)
      {
        for (Int_t i=0; i<channels; i++)
        {
          Double_t oldcounts = histos_cut[11*12+c]->GetBinContent(i);
          histos_cut[11*12+c]->SetBinContent(i, oldcounts + htimes_fadc[11*12+c]->GetBinContent(i));
        }
      }
    }
    // module in a corner
    if (htimes_fadc[0*12+0]->GetEntries() > fadc_int_max_cut && 
        htimes_fadc[1*12+0]->GetEntries() > fadc_int_max_cut &&
        htimes_fadc[0*12+1]->GetEntries() < fadc_int_min_cut &&
        htimes_fadc[1*12+1]->GetEntries() < fadc_int_min_cut)
    {
      for (Int_t i=0; i<channels; i++)
      {
        Double_t oldcounts = histos_cut[0*12+0]->GetBinContent(i);
        histos_cut[0*12+0]->SetBinContent(i, oldcounts + htimes_fadc[0*12+0]->GetBinContent(i));
      }
    }
    if (htimes_fadc[11*12+0]->GetEntries() > fadc_int_max_cut && 
        htimes_fadc[10*12+0]->GetEntries() > fadc_int_max_cut &&
        htimes_fadc[11*12+1]->GetEntries() < fadc_int_min_cut &&
        htimes_fadc[10*12+1]->GetEntries() < fadc_int_min_cut)
    {
      for (Int_t i=0; i<channels; i++)
      {
        Double_t oldcounts = histos_cut[11*12+0]->GetBinContent(i);
        histos_cut[11*12+0]->SetBinContent(i, oldcounts + htimes_fadc[11*12+0]->GetBinContent(i));
      }
    }
    if (htimes_fadc[11*12+11]->GetEntries() > fadc_int_max_cut && 
        htimes_fadc[10*12+11]->GetEntries() > fadc_int_max_cut &&
        htimes_fadc[11*12+10]->GetEntries() < fadc_int_min_cut &&
        htimes_fadc[10*12+10]->GetEntries() < fadc_int_min_cut)
    {
      for (Int_t i=0; i<channels; i++)
      {
        Double_t oldcounts = histos_cut[11*12+11]->GetBinContent(i);
        histos_cut[11*12+11]->SetBinContent(i, oldcounts + htimes_fadc[11*12+11]->GetBinContent(i));
      }
    }
    if (htimes_fadc[0*12+11]->GetEntries() > fadc_int_max_cut && 
        htimes_fadc[1*12+11]->GetEntries() > fadc_int_max_cut &&
        htimes_fadc[0*12+10]->GetEntries() < fadc_int_min_cut &&
        htimes_fadc[1*12+10]->GetEntries() < fadc_int_min_cut)
    {
      for (Int_t i=0; i<channels; i++)
      {
        Double_t oldcounts = histos_cut[0*12+11]->GetBinContent(i);
        histos_cut[0*12+11]->SetBinContent(i, oldcounts + htimes_fadc[0*12+11]->GetBinContent(i));
      }
    }

    if(i%5000==0)
    {
      cout<<i<<" events processed. "<< (i+0.)/loop_max*100.<<" % complete."<<endl;
      //cout<<i<<" events processed. "<<((double)i/(double)loop_max)*100.<<" % complete."<<endl;
    }
    gCurrentEntry++;
  }

  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  //func_landau[0]->Draw();
  //histos[0][0]->Draw("same");
  histos[timerow-1][timecol-1]->Draw();
  //func_landau[15]->Draw("same");
  cout<<"func_landau par[0] = "<<func_landau[0]->GetParameter(0)<<", func_landau par[1] = "<<func_landau[0]->GetParameter(1)<<", func_landau par[2] = "<<func_landau[0]->GetParameter(2)<<endl;
  cout<<"Peak Max = "<<func_landau[0]->GetMaximum()<<", Peak Height = "<<func_landau[0]->GetMaximum()-func_landau[0]->GetParameter(3)<<", Peak Location = "<<func_landau[0]->GetMaximumX()<<endl;


  //Calculate and print timing resolutions for each channel. 

  // last two channels are used as trigger, get their RMS and use their average as reference
  Double_t rms_ref_1 = htimes_fadc[142]->GetRMS();
  Double_t rms_ref_2 = htimes_fadc[143]->GetRMS();
  Double_t rms_ref = (rms_ref_1 + rms_ref_2)/2;
  Double_t rms_ref2 = rms_ref * rms_ref; // square here, just once, to save time in loop
  Double_t rms_current = -1;
  TString title;

  cout <<"rms 1: "<< rms_ref_1 <<" rms 2: "<< rms_ref_2 <<" rms ref: "<< rms_ref <<endl;

  TCanvas *cTimes[12];

  gStyle->SetOptFit(1111);

  TPaveText **text = new TPaveText*[channels];

  for(Int_t i=0; i<12; i++)
  {
    cTimes[i] = new TCanvas(Form("cTimes_Row_%d",i), Form("cTimes_Row_%d",i), 1400, 1000);
    cTimes[i]->SetGrid();
    cTimes[i]->Divide(4,3);

    for(Int_t j=0;j<12;j++)
    {
      cTimes[i]->cd(j+1);
      htimes_fadc[i*12+j]->Draw();
      histos_cut[i*12+j]->Draw("same");
      rms_current = rms[i][j];
      if ( i*12+j < 142 && htimes_fadc[i*12+j]->GetEntries() > 0) // add text pad with time resolution, only to modules with large enough rms, and excluding trigger channels
        //if (rms_current > rms_ref && i*12+j < 142 && htimes_fadc[i*12+j]->GetEntries() > 0) // add text pad with time resolution, only to modules with large enough rms, and excluding trigger channels
      {
        Double_t rms_red = 4*sqrt(abs(rms_current*rms_current - rms_ref2)); // ns
        text[i*12+j] = new TPaveText(.05,.8,.4,.9,"NDC");
        text[i*12+j]->AddText("TDC StdDev: ");  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
        text[i*12+j]->AddText(Form("%0.3f",rms_red));  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
        text[i*12+j]->Draw("same");

        cout << "module: "<< i*12+j << " rms (current) "<< rms_current*4 <<" (reduced) "<< rms_red <<" ns"<<endl;
      }
    }

    cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/fADC_Time/fADC_Time_Run%d_Row%d.png",run,i));
    cTimes[i]->SaveAs(Form("./Timing_Resolution_Plots/fADC_Time/fADC_Time_Run%d_Row%d.C",run,i));
  }

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

#ifndef COMPILE
  return;
#else
  return 0;
#endif
}
