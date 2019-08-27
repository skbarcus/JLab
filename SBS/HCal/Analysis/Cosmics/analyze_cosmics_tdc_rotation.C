const Int_t kNrows=1;
const Int_t kNcols=4;
const Int_t kN=kNcols*kNrows;

void analyze_cosmics_tdc_rotation(Int_t run = 302)
{
  TFile *file = new TFile(TString::Format("%s/cosmics/hcal_tdc_%d.root",
    getenv("HCAL_ROOTFILES"),run),"READ");
  TTree *T = (TTree*)file->Get("TInt");

  TCanvas *canvas = new TCanvas("canvas","canvas",500*kNcols,400*kNrows);
  canvas->Divide(kNcols,kNrows);
  for(Int_t m = 1; m <= kN; m++) {
    canvas->cd(m);
    T->Draw(Form("c%d>>h%d(100,-100,4000)",m,m),Form("smax%d>sped%d+20&&c%d<3000000",m,m,m));
    TH1F *h = (TH1F*)gDirectory->Get(Form("h%d",m));
    h->Fit("landau");
    gPad->SetLogy(1);
  }
}
