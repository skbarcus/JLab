void FF_Comparison()
{
//=========Macro generated from canvas: c5/
//=========  (Wed Apr  6 15:00:14 2022) by ROOT version 6.18/04
   TCanvas *c5 = new TCanvas("c5", "",65,52,700,500);
   gStyle->SetOptStat(0);
   c5->Range(-7.499888,-6.16875,67.49999,55.51875);
   c5->SetFillColor(0);
   c5->SetBorderMode(0);
   c5->SetBorderSize(2);
   c5->SetGridx();
   c5->SetGridy();
   c5->SetFrameBorderMode(0);
   c5->SetFrameBorderMode(0);
   
   TH1F *h1__6 = new TH1F("h1__6","^{3}H Cross Section World Data Distribution",50,0.0001,60);
   h1__6->SetBinContent(1,30);
   h1__6->SetBinContent(2,42);
   h1__6->SetBinContent(3,47);
   h1__6->SetBinContent(4,35);
   h1__6->SetBinContent(5,25);
   h1__6->SetBinContent(6,14);
   h1__6->SetBinContent(7,8);
   h1__6->SetBinContent(8,3);
   h1__6->SetBinContent(9,4);
   h1__6->SetBinContent(10,5);
   h1__6->SetBinContent(11,3);
   h1__6->SetBinContent(12,3);
   h1__6->SetBinContent(13,1);
   h1__6->SetBinContent(14,2);
   h1__6->SetBinContent(16,3);
   h1__6->SetBinContent(17,1);
   h1__6->SetBinContent(18,2);
   h1__6->SetBinContent(20,2);
   h1__6->SetBinContent(22,1);
   h1__6->SetBinContent(23,1);
   h1__6->SetBinContent(25,1);
   h1__6->SetBinContent(27,1);
   h1__6->SetEntries(234);
   h1__6->SetStats(0);
   h1__6->SetFillColor(4);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h1__6->SetLineColor(ci);
   h1__6->GetXaxis()->SetTitle("Q^{2} (fm^{-2})");
   h1__6->GetXaxis()->CenterTitle(true);
   h1__6->GetXaxis()->SetLabelFont(42);
   h1__6->GetXaxis()->SetLabelSize(0.05);
   h1__6->GetXaxis()->SetTitleSize(0.06);
   h1__6->GetXaxis()->SetTitleOffset(0.75);
   h1__6->GetXaxis()->SetTitleFont(42);
   h1__6->GetYaxis()->SetTitle("Number of Measurements");
   h1__6->GetYaxis()->CenterTitle(true);
   h1__6->GetYaxis()->SetLabelFont(42);
   h1__6->GetYaxis()->SetLabelSize(0.05);
   h1__6->GetYaxis()->SetTitleSize(0.06);
   h1__6->GetYaxis()->SetTitleOffset(0.75);
   h1__6->GetYaxis()->SetTitleFont(42);
   h1__6->GetZaxis()->SetLabelFont(42);
   h1__6->GetZaxis()->SetLabelSize(0.035);
   h1__6->GetZaxis()->SetTitleSize(0.035);
   h1__6->GetZaxis()->SetTitleOffset(1);
   h1__6->GetZaxis()->SetTitleFont(42);
   h1__6->Draw("");
   
   TPaveText *pt = new TPaveText(0.1547564,0.94,0.8452436,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("^{3}H Cross Section World Data Distribution");
   pt->Draw();
   c5->Modified();
   c5->cd();
   c5->SetSelected(c5);
}
