#include <iostream>
#include <fstream>
#include <string>

#include "TGraph.h"

using namespace std;

int plotHV2()
{
  const int header = 1; // header lines to skip
  const int npmt = 12 * 12; // number of pmt's
  const float peakAmp = 600; // peak amplitude we want
  const float nevents_threshold = 600; // min number of events to consider the current point


  // files to read
  const char* files[]={
    //"run978.txt",
    //"Salvuccio.txt",
    //"Bogdan.txt",
    "plot_2020_11_10/run978.txt",
    "plot_2020_11_10/run980.txt",
    "plot_2020_11_10/run987.txt",
    "plot_2020_11_10/run988.txt",
    "plot_2020_11_10/run989.txt",
  };
  const int nfiles = sizeof( files ) / sizeof( *files );

  // variables to read lines of the data file.
  string str;      
  float HVArray  [ nfiles ][ npmt ];
  float PeakArray[ nfiles ][ npmt ];
  int channel;
  int ndata[ npmt ];
  int nevents [ nfiles ][ npmt ];

  // variables for plots
  float HV       [ nfiles ];
  float Peak     [ nfiles ];
  TGraph *plots  [ npmt ];
  TString title;
  TCanvas *canvas [ npmt ];
  int npoints;

  // variables for fit
  float minx = 10000;
  float maxx = 0;
  TF1 *powerlaw;
  float y[ npmt ]; // HV to set the pmt's in order to have amplitude = peakAmp
  gStyle->SetOptFit(1111);

  cout << nfiles << endl;

  for (int ifile = 0; ifile < nfiles; ifile++)
  {
    ifstream fp ( files[ ifile ] );
    cout << "Opened file "<< files[ ifile ]<<endl;

    // skip header
    for (int h = 0; h < header; h++)
      getline(fp, str);

    for (int pmt = 0; pmt < npmt; pmt++)
    {
      fp >> channel 
        >> HVArray[ ifile ][ pmt ] 
        >> PeakArray[ ifile ][ pmt ]
        >> nevents [ ifile ][ pmt ];
      if ( PeakArray [ ifile ][ pmt ] < minx ) minx = PeakArray [ ifile ][ pmt ];
      if ( PeakArray [ ifile ][ pmt ] < maxx ) maxx = PeakArray [ ifile ][ pmt ];
    }

    fp.close();
  }
  powerlaw = new TF1("powerlaw","[0] * x**[1]", minx - 100, maxx + 100);
  powerlaw->SetParameter(0,100);
  powerlaw->SetParameter(1, 1);
  powerlaw->SetParName(0,"Normalization");
  powerlaw->SetParName(1,"Exponent");
  for ( int pmt = 0; pmt < npmt; pmt++)
  {
    title.Form("PMT_%d", pmt);
    ndata [ pmt ] = 0;
    for ( int index = 0; index < nfiles; index++)
    {
      // check if we have enough good events
      if ( nevents[ index ][ pmt ] > nevents_threshold )
      {
        HV [ ndata [pmt] ] = HVArray[ index ][ pmt ];
        Peak [ ndata [pmt] ] = PeakArray[ index ][ pmt ];
        ndata[pmt]++;
        //cout << pmt <<" "<< index <<" "<< HV[index] <<" "<< Peak[index] <<endl;
      }
    }
    if ( ndata[pmt] > 2 ) // fit has 2 params, so we need at least 3 points
    {
      canvas [ pmt ] = new TCanvas(title);
      plots [ pmt ] = new TGraph( ndata[pmt], Peak, HV);
      plots [ pmt ]->SetTitle( title );
      plots [ pmt ]->GetXaxis()->SetTitle( "Peak amplitude" );
      plots [ pmt ]->GetYaxis()->SetTitle( "HV (V)");
      plots [ pmt ]->GetYaxis()->SetTitleOffset(1.3);
      plots [ pmt ]->SetMarkerStyle(3);
      plots [ pmt ]->SetMarkerSize(3);
      plots [ pmt ]->Draw("AP");
      plots [ pmt ]->Fit( "powerlaw","Q" );
      y [ pmt ] = powerlaw->Eval( peakAmp );
      canvas[ pmt ]->SaveAs(Form("fit/Fit_PMT_%d.png", pmt));
    } else {
      y [ pmt ] = 0;
    }
  }
  for ( int pmt = 0; pmt < npmt; pmt++)
  {
    cout << "To have a peak amplitude of "<< peakAmp << " set HV of PMT "<< pmt <<" at "<< y [pmt] << " V, #points (fit better if > 3):  "<< ndata[pmt] << endl; 
  }

  return 0;
}
