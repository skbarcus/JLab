
//Script to output several histograms with fits for final sampling fraction plot HCAL - 10.20.20 sseeds
void HCAL_sampFrac_plotgen(){
  TFile *HarmEsumFile = new TFile("HarmEsumFile.root","RECREATE");
  TTree *D;
  TH1D *HarmEsum;
  vector<double> FitMin = {0.0,0.01,0.06,0.08,0.12,0.2,0.3,0.35,0.4,0.5};  //Fit window min params
  vector<double> FitMax = {0.08,0.2,0.35,0.45,0.55,0.65,0.7,0.85,0.95,1.1};  //Fit window max params
  vector<double> emin = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};  //Plot window min params
  vector<double> emax = {0.14,0.25,0.35,0.5,0.6,0.7,0.8,0.9,1.0,1.1};  //Plot window max params
  TVectorD meanValues(10);
  TVectorD meanErrors(10);
  TVectorD IncHadE(10);
  TVectorD IncHadEErr(10);
  TF1 *fFit;
  for(int i=1;i<=10;i++){
    TFile *fOpenedFile = new TFile(Form("gen_HCAL_sampFrac_%dGeV.root",i));
    D = (TTree*)fOpenedFile->Get("T");
    D->Draw(Form("Harm.HCalScint.det.esum>>HarmEsum(100,%f,%f)",emin[i-1],emax[i-1]));
    HarmEsum = (TH1D*)fOpenedFile->Get("HarmEsum");
    HarmEsum->SetTitle(Form("Harm Detector ESum, Incident Hadron %d GeV",i));
    cout << Form("Fit for %dGeV incoming hadrons:",i) <<endl;
    //fFit = new TF1("f1","gaus",FitMin[i-1],FitMax[i-1]);
    //HarmEsum->Fit(fFit);
    HarmEsum->Fit("gaus","","",FitMin[i-1],FitMax[i-1]);  //Fit histogram with gaussian
    fFit=HarmEsum->GetFunction("gaus");
    //cout << "hi" << endl;
    //cout << fFit->GetParameter(1)<<endl;
    meanValues[i-1] = fFit->GetParameter(1);
    meanErrors[i-1] = fFit->GetParError(1);
    IncHadE[i-1] = i;
    IncHadEErr[i-1] = 0;
    //cout << "hi again" << endl;
    HarmEsumFile->cd();
    HarmEsum->Write(Form("HarmEsum%dGeV",i));
    cout << Form("The mean value is: %f", meanValues[i-1]) << endl;
    cout << Form("The error on the mean is: %f", meanErrors[i-1]) << endl;
  }
  TGraphErrors *g1 = new TGraphErrors(IncHadE,meanValues,IncHadEErr,meanErrors);
  TCanvas *C1 = new TCanvas();
  C1->cd();
  g1->SetTitle("HCAL Sampling Fraction");
  g1->Draw("AP");
  g1->Fit("pol1","N");
  TF1 *fUp = new TF1("fUp","[1]*x+[0]",0,10);
  fUp->SetParameters(-0.0351703-0.00125294, 0.0714149+0.000511056);
  TF1 *fDown = new TF1("fDown","[1]*x+[0]",0,10);
  fDown->SetParameters(-0.0351703+0.00125294, 0.0714149-0.000511056);
  TGraph *gfUp = new TGraph(fUp);
  TGraph *gfDown = new TGraph(fDown);
  gfUp->Draw("same");
  gfDown->Draw("same");  //Centroid via wolfram alpha x = 2.45167, y = 0.13992
  TGraph *Centroid = new TGraph();
  Centroid->SetPoint(0,2.45167,0.13991546854652);
  Centroid->SetMarkerSize(3);
  Centroid->SetMarkerStyle(82);
  g1->GetXaxis()->SetTitle("E Incident Hadron [GeV]");
  g1->GetXaxis()->CenterTitle();
  g1->GetYaxis()->SetTitle("E Deposition in HCAL Scintillator [GeV]");
  g1->GetYaxis()->CenterTitle();
  Centroid->Draw("P same");
  cout << "The Centroid is at x = 2.45167, y = 0.13992" << endl;
  return;
  }
