#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
//#include "hcal.h"
#include <vector>
#include <algorithm> 
#include <TStopwatch.h>
#include <TDatime.h>

//#define gmn_tree_cxx
//#include "gmn_tree.h"
#include <TH2.h>
#include <TStyle.h>
using namespace std;

Double_t entry;
Int_t gCurrentEntry = 0;
Int_t limit_evts = 0;                 //0-> analyze all events. 1-> only analyze events up until max_evts.
Int_t loop_max = 0;                   //Dummy variable set equal to max_evts or nevt for the loop.
Int_t max_evts = 2;                   //Maximum number of events to analyze if limit_evts = 1.
Int_t hard_thr = 1;                   //0 = use thresholds calculated as a % of max edep in a single event for BB and HCal. 1 = use the hard coded threshold values below.
Double_t hard_hcal_thr = 0.15;        //Hard coded energy threshold for HCal GeV.
Double_t hard_bb_thr = 2.5;            //Hard coded energy threshold for BB GeV.
Int_t bins = 250;
Int_t nevt;
const Int_t kNrows = 24;
const Int_t kNcols = 12;
const Int_t channels = kNrows * kNcols;

Double_t theta = 14.8;                 //SBS theta in degrees.
Double_t max_edep = 0.;                //Maximum energy deposited in a module GeV.
Int_t max_edep_evt = 0;                //Event containing maximum edep module hit. 
Int_t max_edep_row = 0;                //Row of maximum energy deposited in a module.
Int_t max_edep_col = 0;                //Col of maximum energy deposited in a module.
Double_t max_edep_tot = 0.;            //Maximum edep for a whole event GeV.
Int_t max_edep_tot_evt = 0;              //Event with total maximum edep.
Double_t edep_tot = 0.;                //Total energy deposited for a single event GeV.
Double_t x[1],y[1];                    //Arrays to hold points for tgraph used to draw a mark on the max edep module.
Int_t n = 1;                           //Number of points on the graph
Double_t npe_mev = 5.5;                //Average number of photoelectrons per MeV.

Double_t tot_hit_eng = 0.;             //Total energy in all hits accepted. Used to get average energy per hit.
Double_t tot_evt_eng = 0.;             //Total energy in all events accepted. Used to get average energy per event.
Int_t tot_hits = 0;                    //Counts total number of hits accepted.
Int_t tot_evts = 0;                    //Counts total number of events accepted. 
Double_t edep_tot_bb_ps = 0.;          //Total energy deposited in a single BB preshower event GeV.
Double_t edep_tot_bb_sh = 0.;          //Total energy deposited in a single BB shower event GeV.
Double_t edep_tot_bb = 0.;             //Total energy deposited in a single BB event (preshower + shower) GeV.
Double_t edep_tot_bb_max = 0.;         //Maximum energy deposited in BB preshower and shower for a single event.
Double_t bb_thr = 0.;                  //Energy threshold for accepting HCal events based on a coincident BB hit.
Double_t bb_thr_lvl = 0.5;             //What fraction of the maximum BB energy for a single event defines the BB threshold.
Double_t hcal_thr = 0.;                //Energy threshold for accepting HCal events based on removing low eng HCal events.
Double_t hcal_thr_lvl = 0.1;           //What fraction of the maximum HCal energy for a single event defines the HCal threshold.
vector<Double_t> hcal_hit_eng;         //Vector to store the energies of each individual hcal hit.
vector<Double_t> hcal_evt_eng;         //Vector to store the energies of each hcal event.

TChain *T = 0;
std::string user_input;

//Create Gaussian to fit the timing resolution.
Double_t fit_gaus(Double_t *X,Double_t *par) 
{
  Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((X[0]-par[1])/par[2]),2));
  return fitval;
}

void Energy_Deposition()
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Create a date object.
  TDatime time;

  if(!T) 
    { 

      T = new TChain("T");
      T->Add("/home/skbarcus/JLab/SBS/HCal/Analysis/Simulation/G4SBS/rootfiles/gmn_13.5GeV2_th23_0to43_0ph-30to30.root");

      if(!T->GetEntries())
	{
	  cerr<< "No root file was found" << endl;
	  return;
	}
 
      std::cerr << "Opened up tree with nentries=" << T->GetEntries() << std::endl;
    }
 
  vector<int> *hcal_row = 0;
  vector<int> *hcal_col = 0;
  Int_t hcal_nhits = 0;
  vector<double> *hcal_sumedep = 0;
  vector<double> *hcal_xhit = 0;
  vector<double> *hcal_yhit = 0;
  vector<double> *hcal_xhitg = 0;
  vector<double> *hcal_yhitg = 0;
  vector<double> *hcal_zhitg = 0;
  Int_t bb_gem_nhits = 0;
  Int_t bb_hodo_nhits = 0;
  Int_t bb_ps_nhits = 0;
  vector<int> *bb_ps_row = 0;
  vector<double> *bb_ps_sumedep = 0;
  Int_t bb_sh_nhits = 0;
  vector<int> *bb_sh_row = 0;
  vector<double> *bb_sh_sumedep = 0;
  //vector<int> *bb_gem_nhits = 0;
  //vector<int> *bb_hodo_nhits = 0;
  //vector<int> *bb_ps_nhits = 0;
  //vector<int> *bb_sh_nhits = 0;


  T->SetBranchStatus("Harm.HCalScint.hit.row",1);
  T->SetBranchStatus("Harm.HCalScint.hit.col",1);
  T->SetBranchStatus("Harm.HCalScint.hit.sumedep",1);
  T->SetBranchStatus("Harm.HCalScint.hit.nhits",1);
  T->SetBranchStatus("Harm.HCalScint.hit.xhit",1);
  T->SetBranchStatus("Harm.HCalScint.hit.yhit",1);
  T->SetBranchStatus("Harm.HCalScint.hit.xhitg",1);
  T->SetBranchStatus("Harm.HCalScint.hit.yhitg",1);
  T->SetBranchStatus("Harm.HCalScint.hit.zhitg",1);
  T->SetBranchStatus("Earm.BBGEM.hit.nhits",1);
  T->SetBranchStatus("Earm.BBHodoScint.hit.nhits",1);
  T->SetBranchStatus("Earm.BBPSTF1.hit.nhits",1);
  T->SetBranchStatus("Earm.BBPSTF1.hitrow",1);
  T->SetBranchStatus("Earm.BBPSTF1.hitsumedep",1);
  T->SetBranchStatus("Earm.BBSHTF1.hit.nhits",1);
  T->SetBranchStatus("Earm.BBSHTF1.hitrow",1);
  T->SetBranchStatus("Earm.BBSHTF1.hitsumedep",1);

  T->SetBranchAddress("Harm.HCalScint.hit.row",&hcal_row);
  T->SetBranchAddress("Harm.HCalScint.hit.col",&hcal_col);
  T->SetBranchAddress("Harm.HCalScint.hit.sumedep",&hcal_sumedep);  //GeV
  T->SetBranchAddress("Harm.HCalScint.hit.nhits",&hcal_nhits);  //GeV
  T->SetBranchAddress("Harm.HCalScint.hit.xhit",&hcal_xhit);        //m
  T->SetBranchAddress("Harm.HCalScint.hit.yhit",&hcal_yhit);        //m
  T->SetBranchAddress("Harm.HCalScint.hit.xhitg",&hcal_xhitg);      //m
  T->SetBranchAddress("Harm.HCalScint.hit.yhitg",&hcal_yhitg);      //m
  T->SetBranchAddress("Harm.HCalScint.hit.zhitg",&hcal_zhitg);      //m
  T->SetBranchAddress("Earm.BBGEM.hit.nhits",&bb_gem_nhits);
  T->SetBranchAddress("Earm.BBHodoScint.hit.nhits",&bb_hodo_nhits);
  T->SetBranchAddress("Earm.BBPSTF1.hit.nhits",&bb_ps_nhits);
  T->SetBranchAddress("Earm.BBPSTF1.hit.row",&bb_ps_row);
  T->SetBranchAddress("Earm.BBPSTF1.hit.sumedep",&bb_ps_sumedep);
  T->SetBranchAddress("Earm.BBSHTF1.hit.nhits",&bb_sh_nhits);
  T->SetBranchAddress("Earm.BBSHTF1.hit.row",&bb_sh_row);
  T->SetBranchAddress("Earm.BBSHTF1.hit.sumedep",&bb_sh_sumedep);
  
  //Create 2D histogram for cell hits.
  TH2F *hcell_hits = new TH2F("hcell_hits","Number of Hits per HCal Module",12,0.,12.,24,0.,24.);
  hcell_hits->GetXaxis()->SetTitle("Columns");
  hcell_hits->GetXaxis()->CenterTitle(1);
  hcell_hits->GetYaxis()->SetTitle("Rows");
  hcell_hits->GetYaxis()->CenterTitle(1);

  //Create 2D histogram for (X,Y) hits.
  TH2F *hxy_hits = new TH2F("hxy_hits","Energy Weighted (X,Y) Location of Hits on HCal",30,-1.5,1.5,40,-2,2);//50,-3,3,50,-2,3 looks very nice and rectangular.
  hxy_hits->GetXaxis()->SetTitle("X-Position (m)");
  hxy_hits->GetXaxis()->CenterTitle(1);
  hxy_hits->GetYaxis()->SetTitle("Y-Position (m)");
  hxy_hits->GetYaxis()->CenterTitle(1);

  //Create 1D histogram for the total energy deposited in the scintillators of individual PMTs.
  TH1F *hsumedep = new TH1F("hsumedep","Total Energy Deposited in Scintillators for Each Individual PMT",100,0.,1.);
  hsumedep->GetXaxis()->SetTitle("Energy (GeV)");
  hsumedep->GetXaxis()->CenterTitle(1);
  hsumedep->GetYaxis()->SetTitle("Number of Occurences");
  hsumedep->GetYaxis()->CenterTitle(1);
  hsumedep->GetYaxis()->SetTitleOffset(1.3);
  hsumedep->SetLineWidth(2);

  //Create 1D histogram for the total energy deposited in all HCal scintillators during a single event for the entire run.
  TH1F *hsumedep_evt = new TH1F("hsumedep_evt","Total Energy Deposited in Scintillators for Entire Events",100,0.,1.);
  hsumedep_evt->GetXaxis()->SetTitle("Energy (GeV)");
  hsumedep_evt->GetXaxis()->CenterTitle(1);
  hsumedep_evt->GetYaxis()->SetTitle("Number of Occurences");
  hsumedep_evt->GetYaxis()->CenterTitle(1);
  hsumedep_evt->GetYaxis()->SetTitleOffset(1.3);
  hsumedep_evt->SetLineWidth(2);

  //Create 1D histogram for the total energy deposited in all HCal scintillators during a single event passing energy threshold cuts for HCal and BB.
  TH1F *hsumedep_evt_cuts = new TH1F("hsumedep_evt_cuts","Total Energy Deposited in Scintillators for Entire Events Passing HCal and BB Energy Cuts",100,0.,1.);
  hsumedep_evt_cuts->GetXaxis()->SetTitle("Energy (GeV)");
  hsumedep_evt_cuts->GetXaxis()->CenterTitle(1);
  hsumedep_evt_cuts->GetYaxis()->SetTitle("Number of Occurences");
  hsumedep_evt_cuts->GetYaxis()->CenterTitle(1);
  hsumedep_evt_cuts->GetYaxis()->SetTitleOffset(1.3);
  hsumedep_evt_cuts->SetLineWidth(2);

  //Create 1D histogram for the total energy deposited in the BB preshower.
  TH1F *hsumedep_bb_ps = new TH1F("hsumedep_bb_ps","Total Energy Deposited in BB Preshower",100,0.,3.);
  hsumedep_bb_ps->GetXaxis()->SetTitle("Energy (GeV)");
  hsumedep_bb_ps->GetXaxis()->CenterTitle(1);
  hsumedep_bb_ps->GetYaxis()->SetTitle("Number of Occurences");
  hsumedep_bb_ps->GetYaxis()->CenterTitle(1);
  hsumedep_bb_ps->GetYaxis()->SetTitleOffset(1.3);
  hsumedep_bb_ps->SetLineWidth(2);

  //Create 1D histogram for the total energy deposited in the BB shower.
  TH1F *hsumedep_bb_sh = new TH1F("hsumedep_bb_sh","Total Energy Deposited in BB Shower",100,0.,5.);
  hsumedep_bb_sh->GetXaxis()->SetTitle("Energy (GeV)");
  hsumedep_bb_sh->GetXaxis()->CenterTitle(1);
  hsumedep_bb_sh->GetYaxis()->SetTitle("Number of Occurences");
  hsumedep_bb_sh->GetYaxis()->CenterTitle(1);
  hsumedep_bb_sh->GetYaxis()->SetTitleOffset(1.3);
  hsumedep_bb_sh->SetLineWidth(2);

  //Create 1D histogram for the total energy deposited in the BB preshower and shower combined.
  TH1F *hsumedep_bb_edep_tot = new TH1F("hsumedep_bb_edep_tot","Total Energy Deposited in BB Preshower and Shower",100,0.,5.);
  hsumedep_bb_edep_tot->GetXaxis()->SetTitle("Energy (GeV)");
  hsumedep_bb_edep_tot->GetXaxis()->CenterTitle(1);
  hsumedep_bb_edep_tot->GetYaxis()->SetTitle("Number of Occurences");
  hsumedep_bb_edep_tot->GetYaxis()->CenterTitle(1);
  hsumedep_bb_edep_tot->GetYaxis()->SetTitleOffset(1.3);
  hsumedep_bb_edep_tot->SetLineWidth(2);

  nevt = T->GetEntries();

  //Switch SBS angle from degrees to radians.
  theta = theta * TMath::Pi()/180.;

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
      edep_tot_bb_ps = 0.;
      edep_tot_bb_sh = 0.;

      //cout<<"hcal_row size = "<<(*hcal_row).size()<<". hcal_col size = "<<(*hcal_col).size()<<"."<<endl;

      //Loop over all individual hits in all events purely to get edep for each total event. This lets us cut on this in the histogram filling loop.
      for(Int_t j=0; j<(*hcal_row).size(); j++)
	{
	  //Sum total energy deposited for all hits in event.
	  edep_tot = edep_tot + (*hcal_sumedep)[j];
	}

      //Find HCal event with maximum total edep for all hits combined.
      if(edep_tot > max_edep_tot)
	{
	  max_edep_tot = edep_tot;
	  max_edep_tot_evt = i;
	  hcal_thr = hcal_thr_lvl * max_edep_tot;
	}
      
      //Loop over all BB preshower hits to sum energy.
      for(Int_t j=0; j<(*bb_ps_row).size(); j++)
	{
	  //Sum total energy deposited for all hits in BB preshower event.
	  edep_tot_bb_ps = edep_tot_bb_ps + (*bb_ps_sumedep)[j];
	}

      //Fill BB preshower events energies histo.
      if(bb_ps_nhits > 0)
      	{
	  hsumedep_bb_ps->Fill(edep_tot_bb_ps);
	}

      //Loop over all BB shower hits to sum energy.
      for(Int_t j=0; j<(*bb_sh_row).size(); j++)
	{
	  //Sum total energy deposited for all hits in BB preshower event.
	  edep_tot_bb_sh = edep_tot_bb_sh + (*bb_sh_sumedep)[j];
	}

      //Fill BB shower events energies histo.
      if(bb_sh_nhits > 0)
	{
	  hsumedep_bb_sh->Fill(edep_tot_bb_sh);
	}

      //Sum total BB preshower and shower energies deposited.
      edep_tot_bb = edep_tot_bb_ps + edep_tot_bb_sh;

      //Find maximum energy deposited in BB shower and preshower. Have a given % of this be the threshold.
      if(edep_tot_bb > edep_tot_bb_max)
	{
	  edep_tot_bb_max = edep_tot_bb;
	  bb_thr = bb_thr_lvl * edep_tot_bb_max;
	}

      //Overrides setting the HCal and BB energy thresholds by percentage of maximum edep for an event and instead uses hardcoded values.
      if(hard_thr == 1)
	{
	  hcal_thr = hard_hcal_thr;
	  bb_thr = hard_bb_thr;
	}

      //Fill BB preshower + shower energies histo.
      if(bb_ps_nhits > 0 || bb_sh_nhits > 0)    //Require a hit in either preshower or shower.
      	{
	  hsumedep_bb_edep_tot->Fill(edep_tot_bb);
	}

      //Fill HCal energies histo. No cuts other than ignoring events with no hits.
      if(hcal_nhits > 0)
      	{
	  hsumedep_evt->Fill(edep_tot);
	}

      //Apply an energy threshold to HCal hits and events based on BB PS+SH energy and HCal energy. .
      if(edep_tot > hcal_thr && edep_tot_bb > bb_thr)
	{
	  //Fill histo with total event energies.
	  hsumedep_evt_cuts->Fill(edep_tot);

	  //Fill HCal event energies vector.
	  hcal_evt_eng.push_back(edep_tot);

	  //Count total HCal events passing threshold and sum total energy for average energy calculation.
	  tot_evt_eng = tot_evt_eng + edep_tot;
	  tot_evts++;
	  
	  //Loop over all individual hits in all HCal events.
	  for(Int_t j=0; j<(*hcal_row).size(); j++)
	    {
	      //Fill (row,col) hits, (X,Y) hits, and edep histos.
	      hcell_hits->Fill((*hcal_col)[j],(*hcal_row)[j]);
	      //hxy_hits->Fill((*hcal_xhitg)[j],(*hcal_yhitg)[j]);
	      hxy_hits->Fill((*hcal_xhitg)[j]*cos(theta)+(*hcal_zhitg)[j]*sin(theta),(*hcal_yhitg)[j]-0.45);
	      hsumedep->Fill((*hcal_sumedep)[j]);
	  
	      //Fill HCal hit energies vector.
	      hcal_hit_eng.push_back((*hcal_sumedep)[j]);

	      //Sum total hit energy (should match total event energy) and count number of hits passing threshold.
	      tot_hit_eng = tot_hit_eng + (*hcal_sumedep)[j];
	      tot_hits++;
	      
	      //Find maximal energy deposited in a module.
	      if((*hcal_sumedep)[j] > max_edep)
		{
		  max_edep = (*hcal_sumedep)[j];
		  max_edep_row = (*hcal_row)[j];
		  max_edep_col = (*hcal_col)[j];
		  max_edep_evt = i;
		}
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

  x[0] = max_edep_col + 0.5;
  y[0] = max_edep_row + 0.5;

  TGraph *gmax_edep = new TGraph (n, x, y);
  gmax_edep->Draw("same *");

  TCanvas* cxy_hits=new TCanvas("cxy_hits");
  cxy_hits->SetGrid();
  hxy_hits->Draw("colz");

  TCanvas* csumedep=new TCanvas("csumedep");
  csumedep->SetGrid();
  csumedep->SetLogy();
  hsumedep->Draw("");

  TCanvas* csumedep_evt=new TCanvas("csumedep_evt");
  csumedep_evt->SetGrid();
  hsumedep_evt->Draw("");
  csumedep_evt->Update();//Need update to get gPad->GetUymax() to work.
  TLine *line_hcal_thr = new TLine(hcal_thr,0,hcal_thr,gPad->GetUymax());
  line_hcal_thr->SetLineColor(2);
  line_hcal_thr->SetLineWidth(2);
  line_hcal_thr->Draw("same");

  TCanvas* csumedep_evt_cuts=new TCanvas("csumedep_evt_cuts");
  csumedep_evt_cuts->SetGrid();
  hsumedep_evt_cuts->Draw("");

  TF1 *func_gaus_fit_edep_evt = new TF1("func_gaus_fit_edep_evt",fit_gaus, 0., 1., 3);
  func_gaus_fit_edep_evt->SetLineColor(2);
  func_gaus_fit_edep_evt->SetNpx(1000);
  func_gaus_fit_edep_evt->SetParameter(0,hsumedep_evt_cuts->GetMaximum());
  func_gaus_fit_edep_evt->SetParameter(1,hsumedep_evt_cuts->GetMean());
  func_gaus_fit_edep_evt->SetParameter(2,hsumedep_evt_cuts->GetStdDev());
  hsumedep_evt_cuts->Fit(func_gaus_fit_edep_evt,"");
  func_gaus_fit_edep_evt->Draw("same");

  TCanvas* csumedep_bb_ps=new TCanvas("csumedep_bb_ps");
  csumedep_bb_ps->Divide(1,3);
  csumedep_bb_ps->cd(1);
  csumedep_bb_ps->SetGrid();
  hsumedep_bb_ps->Draw("");

  //TCanvas* csumedep_bb_sh=new TCanvas("csumedep_bb_sh");
  //csumedep_bb_sh->SetGrid();
  csumedep_bb_ps->cd(2);
  hsumedep_bb_sh->Draw("");

  //TCanvas* csumedep_bb_edep_tot=new TCanvas("csumedep_bb_edep_tot");
  //csumedep_bb_edep_tot->SetGrid();
  csumedep_bb_ps->cd(3);
  hsumedep_bb_edep_tot->Draw("");
  //csumedep_bb_edep_tot->Update();//Need update to get gPad->GetUymax() to work.
  csumedep_bb_ps->Update();
  TLine *line_bb_thr = new TLine(bb_thr,0,bb_thr,gPad->GetUymax());
  line_bb_thr->SetLineColor(2);
  line_bb_thr->SetLineWidth(2);
  line_bb_thr->Draw("same");

  TCanvas* csumedep_bb_edep_tot=new TCanvas("csumedep_bb_edep_tot");
  csumedep_bb_edep_tot->SetGrid();
  hsumedep_bb_edep_tot->Draw("");
  line_bb_thr->Draw("same");
   

  //Print the module with the maximal edep and its location.
  cout<<"The maximum energy deposition of "<<max_edep*1000.<<" MeV ("<<max_edep*1000.*npe_mev<<" PE) was deposited in row "<<max_edep_row<<" col "<<max_edep_col<<" during event "<<max_edep_evt<<"."<<endl;

  //Print max edep for whole event and the event.
  cout<<"Maximum energy deposited for all hits in a single event = "<<max_edep_tot*1000.<<" MeV ("<<max_edep_tot*1000.*npe_mev<<" PE) for event "<<max_edep_tot_evt<<"."<<endl;

  //Print number of events that passed the threshold.
  cout<<"Number of events passing threshold = "<<tot_evts<<endl;

  //Print the average energy deposited in a PMT's scintillator and the average energy deposited in an event.
  cout<<"Total hits accepted = "<<tot_hits<<". Total energy deposited in scintillators of accepted hits = "<<tot_hit_eng<<" GeV. Average energy deposited in a PMT's scintillator per hit = "<<(tot_hit_eng/tot_hits)*1000.<<" MeV."<<endl;
  cout<<"Total events accepted = "<<tot_evts<<". Total energy deposited in scintillators of accepted events = "<<tot_evt_eng<<" GeV. Average energy deposited in scintillator's per event = "<<(tot_evt_eng/tot_evts)*1000.<<" MeV."<<endl;

  cout<<"The standard deviation of the energy deposited in events is "<<func_gaus_fit_edep_evt->GetParameter(2)*1000.<<" MeV."<<endl;

  cout<<"Max energy deposited in BB preshower and shower for a single event = "<<edep_tot_bb_max<<". BB threshold = "<<bb_thr*1000<<" MeV."<<endl;

  //Sort the HCal accepted event energies in ascending order.
  sort(hcal_evt_eng.begin(), hcal_evt_eng.end());
  /*
  for(Int_t i=0;i<hcal_evt_eng.size();i++)
    {
      cout<<hcal_evt_eng[i]<<" ";
    }
  */

  //Sort the HCal accepted hit energies in ascending order.
  sort(hcal_hit_eng.begin(), hcal_hit_eng.end());
  /*
  for(Int_t i=0;i<hcal_hit_eng.size();i++)
    {
      cout<<hcal_hit_eng[i]<<" ";
    }
  */

  cout<<"Max energy for a single hit passing thresholds = "<<hcal_hit_eng.back()*1000.<<" MeV."<<endl;
  Int_t hcal_hit_995;
  hcal_hit_995 = tot_hits * 0.995;
  cout<<"hcal_hit_995 = "<<hcal_hit_995<<" has an energy of "<<hcal_hit_eng[hcal_hit_995]*1000.<<" MeV."<<endl;

  cout<<"Max energy for a single event passing thresholds = "<<hcal_evt_eng.back()*1000.<<" MeV."<<endl;
  Int_t hcal_evt_995;
  hcal_evt_995 = tot_evts * 0.995;
  cout<<"hcal_evt_995 = "<<hcal_evt_995<<" has an energy of "<<hcal_evt_eng[hcal_evt_995]*1000.<<" MeV."<<endl;
 
  //Stop and print out stopwatch information.
  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
