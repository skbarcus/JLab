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

void TPaveText_Test()
{
  TH1F **htimes = new TH1F*[24];
  for(Int_t i = 0; i<24; i++)
    {
      htimes[i] = new TH1F(Form("htimes%d",i),Form("Module %d",i),100,0,100);
    }

  for(Int_t i=0; i<24; i++)
    {
      for(Int_t j=0;j<10;j++)
	{
	  htimes[i]->Fill(pow(j,2.));
	}
    }

  TCanvas *cTimes[2];

  TPaveText **text = new TPaveText*[24];
  for(Int_t i = 0; i<24; i++)
    {
      text[i] = new TPaveText(10.0,.70,40.,0.5);//.05,.3,.95,.6 //top left (x1,y1) then lower right (x2,y2).
    }

  for(Int_t i=0; i<2; i++)
    {
      cTimes[i] = new TCanvas(Form("cTimes_Row_%d",i));
      cTimes[i]->SetGrid();
      cTimes[i]->Divide(4,3);
      
      for(Int_t j=0;j<12;j++)
	{
	  cTimes[i]->cd(j+1);
	  htimes[i*12+j]->Draw();
	  text[i*12+j]->AddText("This line is blue"); 
	  ((TText*)text[i*12+j]->GetListOfLines()->Last())->SetTextColor(kBlue);
	  text[i*12+j]->Draw("same");
	}
    }
}
