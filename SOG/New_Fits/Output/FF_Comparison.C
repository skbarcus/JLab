void FF_Comparison()
{
//=========Macro generated from canvas: c5/
//=========  (Thu Mar 31 18:45:21 2022) by ROOT version 6.18/04
   TCanvas *c5 = new TCanvas("c5", "",65,52,700,500);
   gStyle->SetOptStat(0);
   c5->Range(-8.749888,-4.9875,78.74999,44.8875);
   c5->SetFillColor(0);
   c5->SetBorderMode(0);
   c5->SetBorderSize(2);
   c5->SetGridx();
   c5->SetGridy();
   c5->SetFrameBorderMode(0);
   c5->SetFrameBorderMode(0);
   
   TH1F *h1__6 = new TH1F("h1__6","^{3}He Cross Section World Data Distribution",70,0.0001,70);
   h1__6->SetBinContent(1,31);
   h1__6->SetBinContent(2,30);
   h1__6->SetBinContent(3,38);
   h1__6->SetBinContent(4,33);
   h1__6->SetBinContent(5,25);
   h1__6->SetBinContent(6,16);
   h1__6->SetBinContent(7,7);
   h1__6->SetBinContent(8,5);
   h1__6->SetBinContent(9,9);
   h1__6->SetBinContent(10,3);
   h1__6->SetBinContent(11,4);
   h1__6->SetBinContent(12,3);
   h1__6->SetBinContent(13,3);
   h1__6->SetBinContent(14,3);
   h1__6->SetBinContent(15,1);
   h1__6->SetBinContent(16,3);
   h1__6->SetBinContent(18,1);
   h1__6->SetBinContent(19,1);
   h1__6->SetBinContent(20,3);
   h1__6->SetBinContent(21,1);
   h1__6->SetBinContent(23,3);
   h1__6->SetBinContent(25,4);
   h1__6->SetBinContent(26,2);
   h1__6->SetBinContent(28,1);
   h1__6->SetBinContent(30,3);
   h1__6->SetBinContent(31,2);
   h1__6->SetBinContent(32,1);
   h1__6->SetBinContent(34,2);
   h1__6->SetBinContent(35,2);
   h1__6->SetBinContent(36,1);
   h1__6->SetBinContent(37,1);
   h1__6->SetBinContent(41,2);
   h1__6->SetBinContent(42,1);
   h1__6->SetBinContent(43,1);
   h1__6->SetBinContent(45,3);
   h1__6->SetBinContent(47,1);
   h1__6->SetBinContent(49,1);
   h1__6->SetBinContent(52,1);
   h1__6->SetBinContent(55,1);
   h1__6->SetBinContent(56,1);
   h1__6->SetBinContent(58,1);
   h1__6->SetBinContent(61,2);
   h1__6->SetBinContent(65,1);
   h1__6->SetBinContent(71,1);
   h1__6->SetEntries(259);
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
   
   TPaveText *pt = new TPaveText(0.15,0.94,0.85,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("^{3}He Cross Section World Data Distribution");
   pt->Draw();
   c5->Modified();
   c5->cd();
   c5->SetSelected(c5);
}
