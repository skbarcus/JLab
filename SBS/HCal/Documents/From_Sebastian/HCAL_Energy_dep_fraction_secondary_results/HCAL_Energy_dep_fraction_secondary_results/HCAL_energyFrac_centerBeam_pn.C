//Script to output several histograms with fits for proton and neutron energy deposition fraction in HCAL plot - 10.20.23 sseeds
void HCAL_energyFrac_centerBeam_pn(){
  TFile *HarmEsumFile = new TFile("HarmEsumFile.root","RECREATE");
  TTree *D;
  TH1D *HarmEsum_p;
  TH1D *HarmEsum_n;
  //All window and fit parameters determined QUALITATIVELY with effort to capture full distribution and cut low amp events. Some low inc E neutron histos should probably be fit with a poisson.
  //Low amplitude events correspond to losses of secondary events out of the HCAL scaling with inc had E? Seems to make sense

  //Protons
  vector<double> FitMin_p = {0.0,0.01,0.12,0.18,0.22,0.3,0.38,0.42,0.52,0.58};  //Fit window min params
  vector<double> FitMax_p = {0.08,0.2,0.35,0.45,0.55,0.65,0.7,0.85,0.95,1.1};  //Fit window max params
  vector<double> emin_p = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};  //Plot window min params
  vector<double> emax_p = {0.14,0.25,0.35,0.5,0.6,0.7,0.8,0.9,1.0,1.1};  //Plot window max params
  TVectorD meanValues_p(10); //Vector to contain means from guassian fits for the proton
  TVectorD meanErrors_p(10); //Vector to conatin mean errors from fits for the proton
  TVectorD IncHadE_p(10); //Same for neutron
  TVectorD IncHadEErr_p(10); //Same for neutron
  TF1 *fFit_p;
  //TCanvas *C1 = new TCanvas();  //Canvas to save fits to E dep histos
  //C1->Divide(5,2);  //Not working 10.23.20
  for(int i=1; i<=10; i++){
    TFile *fOpenedFile = new TFile(Form("gen_HCAL_sampFrac_%dGeVbeam_p.root",i));  //Tailored to protons
    D = (TTree*)fOpenedFile->Get("T");
    D->Draw(Form("Harm.HCalScint.det.esum>>HarmEsum_p(100,%f,%f)",emin_p[i-1],emax_p[i-1]));
    HarmEsum_p = (TH1D*)fOpenedFile->Get("HarmEsum_p");
    HarmEsum_p->SetTitle(Form("Harm Detector ESum, Incident Proton %d GeV",i));
    HarmEsum_p->SetMinimum(0.0);
    cout << Form("Fit for %dGeV incoming protons:",i) <<endl;
    HarmEsum_p->Fit("gaus","","",FitMin_p[i-1],FitMax_p[i-1]);  //Fit histogram with gaussian
    fFit_p=HarmEsum_p->GetFunction("gaus");
    meanValues_p[i-1] = fFit_p->GetParameter(1);
    meanErrors_p[i-1] = fFit_p->GetParError(1);
    IncHadE_p[i-1] = i;
    IncHadEErr_p[i-1] = 0;
    HarmEsumFile->cd();
    HarmEsum_p->Write(Form("HarmEsum_p%dGeV",i));
    HarmEsum_p->Draw("AP");
    //C1->cd(i);
    //HarmEsum_p->Draw("AP");
    cout << Form("The mean value is: %f", meanValues_p[i-1]) << endl;
    cout << Form("The error on the mean is: %f", meanErrors_p[i-1]) << endl;
  }
  //cout << "Saving fit to E dep by proton E plots.." << endl;
  //C1->SaveAs("HCAL_EDepPlots_p.pdf");
  
  //Neutrons
  vector<double> FitMin_n = {0.0,0.01,0.12,0.18,0.22,0.3,0.38,0.42,0.52,0.58};  //Fit window min params
  vector<double> FitMax_n = {0.05,0.13,0.35,0.45,0.55,0.65,0.7,0.85,0.95,1.1};  //Fit window max params
  vector<double> emin_n = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};  //Plot window min params
  vector<double> emax_n = {0.14,0.25,0.35,0.5,0.6,0.7,0.8,0.9,1.0,1.1};  //Plot window max params
  TVectorD meanValues_n(10);
  TVectorD meanErrors_n(10);
  TVectorD IncHadE_n(10);
  TVectorD IncHadEErr_n(10);
  TF1 *fFit_n;
  //TCanvas *C2 = new TCanvas();  //Canvas to save fits to E dep histos
  //C2->Divide(5,2);
  for(int i=1; i<=10; i++){
    TFile *fOpenedFile = new TFile(Form("gen_HCAL_sampFrac_%dGeVbeam_n.root",i));  //Tailored to protons
    D = (TTree*)fOpenedFile->Get("T");
    D->Draw(Form("Harm.HCalScint.det.esum>>HarmEsum_n(100,%f,%f)",emin_n[i-1],emax_n[i-1]));
    HarmEsum_n = (TH1D*)fOpenedFile->Get("HarmEsum_n");
    HarmEsum_n->SetTitle(Form("Harm Detector ESum, Incident Neutron %d GeV",i));
    HarmEsum_n->SetMinimum(0.0);
    cout << Form("Fit for %dGeV incoming neutrons:",i) <<endl;
    HarmEsum_n->Fit("gaus","","",FitMin_n[i-1],FitMax_n[i-1]);  //Fit histogram with gaussian
    fFit_n=HarmEsum_n->GetFunction("gaus");
    meanValues_n[i-1] = fFit_n->GetParameter(1);
    meanErrors_n[i-1] = fFit_n->GetParError(1);
    IncHadE_n[i-1] = i;
    IncHadEErr_n[i-1] = 0;
    HarmEsumFile->cd();
    HarmEsum_n->Write(Form("HarmEsum_n%dGeV",i));
    HarmEsum_n->Draw("AP");
    //C2->cd(i);
    //HarmEsum_n->Draw("AP");
    cout << Form("The mean value is: %f", meanValues_n[i-1]) << endl;
    cout << Form("The error on the mean is: %f", meanErrors_n[i-1]) << endl;
  }
  //cout << "Saving fit to E dep by neutron E plots.." << endl;
  //C2->SaveAs("HCAL_EDepPlots_n.pdf");


  //Protons
  TCanvas *C3 = new TCanvas(); //Canvas to contain proton centroid fit and fit lines
  C3->cd();
  TGraphErrors *g1_p = new TGraphErrors(IncHadE_p,meanValues_p,IncHadEErr_p,meanErrors_p); //Construct E dep frac plot
  TF1 *gFit_p;
  g1_p->SetTitle("HCAL Proton Energy Deposition");
  g1_p->SetMinimum(0.0);
  g1_p->SetMarkerColor(2);
  g1_p->Draw("AP");
  g1_p->Fit("pol1");
  gFit_p = g1_p->GetFunction("pol1");
  gFit_p->SetLineColor(0);
  double p0_p = gFit_p->GetParameter("p0");
  double p0err_p = gFit_p->GetParError(0);
  double p1_p = gFit_p->GetParameter(1);
  double p1err_p = gFit_p->GetParError(1);
  cout << p0_p << endl;
  cout << p1_p << endl;
  TF1 *fUp_p = new TF1("fUp","[1]*x+[0]",0,10);
  TF1 *fDown_p = new TF1("fDown","[1]*x+[0]",0,10);
  fUp_p->SetParameters(p0_p-p0err_p, p1_p+p1err_p); 
  fDown_p->SetParameters(p0_p+p0err_p, p1_p-p1err_p);
  fUp_p->SetLineColor(2);
  fDown_p->SetLineColor(2);
  TGraph *gfUp_p = new TGraph(fUp_p);
  TGraph *gfDown_p = new TGraph(fDown_p);
  gfUp_p->Draw("same");
  gfDown_p->Draw("same");  
  TGraph *Centroid_p = new TGraph();
  double centX_p = ((p0_p+p0err_p)-(p0_p-p0err_p))/((p1_p+p1err_p)-(p1_p-p1err_p));  //Determine crossing point x
  double centY_p = (p1_p+p1err_p)*centX_p+(p0_p-p0err_p);  //Determine crossing point y

  Centroid_p->SetPoint(0,centX_p,centY_p);
  Centroid_p->SetMarkerSize(3);
  Centroid_p->SetMarkerColor(2);
  Centroid_p->SetMarkerStyle(82);
  g1_p->GetXaxis()->SetTitle("E Incident Proton [GeV]");
  g1_p->GetXaxis()->CenterTitle();
  g1_p->GetYaxis()->SetTitle("E Deposition in HCAL Scintillator [GeV]");
  g1_p->GetYaxis()->CenterTitle();
  Centroid_p->Draw("P same");
  HarmEsumFile->cd();
  
  cout << Form("The proton centroid is at x = %f, y = %f", centX_p, centY_p) << endl;

  //Neutrons
  TCanvas *C4 = new TCanvas(); //Canvas to contain neutron centroid fit and fit lines
  C4->cd();
  TGraphErrors *g1_n = new TGraphErrors(IncHadE_n,meanValues_n,IncHadEErr_n,meanErrors_n); //Construct E dep frac plot
  TF1 *gFit_n;
  g1_n->SetTitle("HCAL Neutron Energy Deposition");
  g1_n->SetMinimum(0.0);
  g1_n->SetMarkerColor(4);
  g1_n->Draw("AP");
  g1_n->Fit("pol1");
  gFit_n = g1_n->GetFunction("pol1");
  gFit_n->SetLineColor(0);
  double p0_n = gFit_n->GetParameter("p0");
  double p0err_n = gFit_n->GetParError(0);
  double p1_n = gFit_n->GetParameter(1);
  double p1err_n = gFit_n->GetParError(1);
  cout << p0_n << endl;
  cout << p1_n << endl;
  TF1 *fUp_n = new TF1("fUp","[1]*x+[0]",0,10);
  TF1 *fDown_n = new TF1("fDown","[1]*x+[0]",0,10);
  fUp_n->SetParameters(p0_n-p0err_n, p1_n+p1err_n); 
  fDown_n->SetParameters(p0_n+p0err_n, p1_n-p1err_n);
  fUp_n->SetLineColor(4);
  fDown_n->SetLineColor(4);
  TGraph *gfUp_n = new TGraph(fUp_n);
  TGraph *gfDown_n = new TGraph(fDown_n);
  gfUp_n->Draw("same");
  gfDown_n->Draw("same");  
  TGraph *Centroid_n = new TGraph();
  double centX_n = ((p0_n+p0err_n)-(p0_n-p0err_n))/((p1_n+p1err_n)-(p1_n-p1err_n));  //Determine crossing point x
  double centY_n = (p1_n+p1err_n)*centX_n+(p0_n-p0err_n);  //Determine crossing point y

  Centroid_n->SetPoint(0,centX_n,centY_n);
  Centroid_n->SetMarkerSize(3);
  Centroid_n->SetMarkerColor(4);
  Centroid_n->SetMarkerStyle(84);
  g1_n->GetXaxis()->SetTitle("E Incident Neutron [GeV]");
  g1_n->GetXaxis()->CenterTitle();
  g1_n->GetYaxis()->SetTitle("E Deposition in HCAL Scintillator [GeV]");
  g1_n->GetYaxis()->CenterTitle();
  Centroid_n->Draw("P same");
  cout << Form("The neutron centroid is at x = %f, y = %f", centX_n, centY_n) << endl;

  //Save canvas containing final plot
  //C4->SaveAs("HCAL_EdepVSEhad_pn.pdf");

  //Make final plot with proton energy signatures (red) and neutron energy signatures (blue)
  cout << "Drawing E Dep vs E Inc plot" << endl;
  TVectorD meanValues(20);
  TVectorD meanErrors(20);
  TVectorD IncHadE(20);
  TVectorD IncHadEErr(20);

  //TCanvas *C5 = new TCanvas();
  auto *C5 = new TCanvas();
  C5->cd();
  
  for(int i=0; i<10; i++){  //Add data together for final plot
    meanValues[i] = meanValues_p[i];
    meanErrors[i] = meanErrors_p[i];
    IncHadE[i] = IncHadE_p[i];
    IncHadEErr[i] = IncHadEErr_p[i];
    meanValues[i+10] = meanValues_n[i];
    meanErrors[i+10] = meanErrors_p[i];
    IncHadE[i+10] = IncHadE_p[i];
    IncHadEErr[i+10] = IncHadEErr_p[i];
  }
  //TGraphErrors *EvE_pn = new TGraphErrors(IncHadE,meanValues,IncHadEErr,meanErrors);
  auto EvE_pn = new TGraphErrors(IncHadE,meanValues,IncHadEErr,meanErrors);
  
  EvE_pn->SetTitle("Proton and Neutron Energy Deposition in HCAL Scintillator vs Incident Hadron Energy");
  EvE_pn->GetXaxis()->SetTitle("E Incident Hadron [GeV]");
  EvE_pn->GetXaxis()->CenterTitle();
  EvE_pn->GetYaxis()->SetTitle("E Deposition in HCAL Scintillator [GeV]");
  EvE_pn->GetYaxis()->CenterTitle();
  EvE_pn->SetMinimum(0.0);
  EvE_pn->Draw("AP");

  //TGraph *ProtonGraph = new TGraph();
  auto ProtonGraph = new TGraph();
  for(int i=0; i<10; i++){
    ProtonGraph->SetPoint(i,IncHadE_p[i],meanValues_p[i]);
  }
  ProtonGraph->SetMarkerColor(2);
  ProtonGraph->SetMarkerStyle(4);
  ProtonGraph->Draw("P same");

  //TGraph *NeutronGraph = new TGraph();
  auto NeutronGraph = new TGraph();
  for(int i=0; i<10; i++){
    NeutronGraph->SetPoint(i,IncHadE_n[i],meanValues_n[i]);
  }
  NeutronGraph->SetMarkerColor(4);
  NeutronGraph->SetMarkerStyle(25);
  NeutronGraph->Draw("P same");

  //Include centroids
  Centroid_p->Draw("P same");
  Centroid_n->Draw("P same");

  //Add legend
  //TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetHeader("Incident Hadrons","C"); 		        
  legend->AddEntry(ProtonGraph,"Protons","p");
  legend->AddEntry(NeutronGraph,"Neutrons","p");
  legend->Draw();
  
  return;
}
