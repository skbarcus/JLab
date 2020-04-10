#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
//#include "hcal.h"
#include <vector>
#include <TStopwatch.h>
#include <TDatime.h>

//#define gmn_tree_cxx
//#include "gmn_tree.h"
#include <TH2.h>
#include <TStyle.h>
using namespace std;

Double_t entry;
Int_t gCurrentEntry = 0;
Int_t limit_evts = 0;                     //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t loop_max = 0;                       //Dummy variable set equal to max_evts or nevt for the loop.
Int_t max_evts = 2;                   //Maximum number of events to analyze if limit_evts = 1.
Int_t bins = 250;
Int_t nevt;
const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;

Double_t thr_eng = 0.010;              //Minimal energy threshold GeV. Default 10 MeV theshold in G4sbs unless otherwise specified.
Double_t max_edep = 0.;                //Maximum energy deposited in a module GeV.
Int_t max_edep_evt = 0;                //Event containing maximum edep module hit. 
Int_t max_edep_row = 0;                //Row of maximum energy deposited in a module.
Int_t max_edep_col = 0;                //Col of maximum energy deposited in a module.
Double_t max_edep_tot = 0.;            //Maximum edep for a whole event GeV.
Int_t max_edep_tot_evt = 0;              //Event with total maximum edep.
Double_t edep_tot = 0.;                //Total energy deposited for a single event GeV.
Double_t x[1],y[1];                    //Arrays to hold points for tgraph used to draw a mark on the max edep module.
Int_t n = 1;                           //Number of points on the graph

TChain *T = 0;
std::string user_input;

void Test()
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Create a date object.
  TDatime time;

  if(!T) 
    { 

      T = new TChain("T");
      T->Add("/lustre19/expphy/volatile/halla/sbs/skbarcus/rootfiles/gmn_3.5GeV2_th22_5to42_5ph-30to30.root");

      if(!T->GetEntries())
	{
	  cerr<< "No root file was found" << endl;
	  return;
	}
 
      std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    }
 
  vector<int> *hcal_row = 0;
  vector<int> *hcal_col = 0;
  vector<double> *hcal_sumedep = 0;
  vector<double> *hcal_xhit = 0;
  vector<double> *hcal_yhit = 0;

  T->SetBranchStatus("Harm.HCalScint.hit.row",1);
  T->SetBranchStatus("Harm.HCalScint.hit.col",1);
  T->SetBranchStatus("Harm.HCalScint.hit.sumedep",1);
  T->SetBranchStatus("Harm.HCalScint.hit.xhit",1);
  T->SetBranchStatus("Harm.HCalScint.hit.yhit",1);

  T->SetBranchAddress("Harm.HCalScint.hit.row",&hcal_row);
  T->SetBranchAddress("Harm.HCalScint.hit.col",&hcal_col);
  T->SetBranchAddress("Harm.HCalScint.hit.sumedep",&hcal_sumedep);  //GeV
  T->SetBranchAddress("Harm.HCalScint.hit.xhit",&hcal_xhit);  //m
  T->SetBranchAddress("Harm.HCalScint.hit.yhit",&hcal_yhit);  //m
  
  //Create 2D histogram for cell hits.
  TH2F *hcell_hits = new TH2F("hcell_hits","Number of Hits per Cell",12,0.,12.,24,0.,24.);

  //Create 2D histogram for (X,Y) hits.
  TH2F *hxy_hits = new TH2F("hxy_hits","(X,Y) Location of Hits",100,-0.1,0.1,100,-0.1,0.1);

  //Create 1D histogram for the total energy deposited in the scintillators.
  TH1F *hsumedep = new TH1F("hsumedep","Total Energy Deposited in Scintillators",100,0.,1.);

  nevt = T->GetEntries();

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

      //cout<<"*** "<<(*hcal_row)[i]<<endl;     //de-referenced pointer
      //cout<<"*** "<<hcal_row->at(i)<<endl;      //The at function

      /*
      for(Int_t j=0; j<(*hcal_row).size(); j++)
	{
	  cout<<"Event "<<i<<": nhits = "<<(*hcal_row).size()<<" row = "<<(*hcal_row)[j]<<" col = "<<(*hcal_col)[j]<<endl;
	}
      */

      //Reset the total energy deposited in a single event.
      edep_tot = 0.;
      
      //Loop over all individual hits in all events.
      for(Int_t j=0; j<(*hcal_row).size(); j++)
	{
	  //Add threshold cut for energy deposited.
	  if((*hcal_sumedep)[j] > thr_eng)
	    {
	      //Fill (row,col) hits, (X,Y) hits, and edep histos.
	      hcell_hits->Fill((*hcal_col)[j],(*hcal_row)[j]);
	      hxy_hits->Fill((*hcal_xhit)[j],(*hcal_yhit)[j]);
	      hsumedep->Fill((*hcal_sumedep)[j]);
	      
	      //Sum total energy deposited for all hits in event.
	      edep_tot = edep_tot + (*hcal_sumedep)[j];

	      //Find maximal energy deposited in a module.
	      if((*hcal_sumedep)[j] > max_edep)
		{
		  max_edep = (*hcal_sumedep)[j];
		  max_edep_row = (*hcal_row)[j];
		  max_edep_col = (*hcal_col)[j];
		  max_edep_evt = i;
		}
	    }
	  //Find event with maximum total edep for all hits combined.
	  if(edep_tot > max_edep_tot)
	    {
	      max_edep_tot = edep_tot;
	      max_edep_tot_evt = i;
	    }
	}

      if(i%5000==0)
	{
	  cout<<i<<" events processed. "<<((double)i/(double)loop_max)*100.<<" % complete."<<endl;
	}

      gCurrentEntry++;
    }
      
  TCanvas* cmod_hits=new TCanvas("cmod_hits");
  cmod_hits->SetGrid();
  hcell_hits->Draw("colz");

  x[0] = max_edep_col - 0.5;
  y[0] = max_edep_row - 0.5;
  TGraph *gmax_edep = new TGraph (n, x, y);
  gmax_edep->Draw("same *");

  TCanvas* cxy_hits=new TCanvas("cxy_hits");
  cxy_hits->SetGrid();
  hxy_hits->Draw("colz");

  TCanvas* csumedep=new TCanvas("csumedep");
  csumedep->SetGrid();
  hsumedep->Draw("");

  //Print the module with the maximal edep and its location.
  cout<<"The maximum energy deposition of "<<max_edep*1000.<<" MeV was deposited in row "<<max_edep_row<<" col "<<max_edep_col<<" during event "<<max_edep_evt<<"."<<endl;

  //Print max edep for whole event and the event.
  cout<<"Maximum energy deposited for all hits in a single event = "<<max_edep_tot*1000<<" MeV for event "<<max_edep_tot_evt<<"."<<endl;

  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
