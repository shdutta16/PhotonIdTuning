#include <cstdlib>
#include <iostream>


#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TString.h"

void ROCcomp(){


  gStyle->SetOptStat(0);

  TString f1_name = "../output_1M200K/Loose_w_scaled.root";   // BAR
  TFile *f1 = new TFile( f1_name );
  if ( !f1 || !f1->IsOpen() ) {
    cout << "\nERROR! Could not open root file " << f1_name
         << endl;
    exit(0);
  }


  f1->cd("Method_Cuts/Cut_Loose_r");
  TH1F *SC = (TH1F*)MVA_Cut_Loose_r_trainingRejBvsS->Clone();


  
  TString f2_name = "../output_1M200K/Loose_w_scaled.root"      // BAR_REL
  TFile *f2 = new TFile( f2_name );
  f2->cd("Method_Cuts/Cut_Loose_r");
  TH1F *REL = (TH1F*)MVA_Cut_Loose_r_trainingRejBvsS->Clone();

  SC->SetLineColor(kGray +3);
  REL->SetLineColor(kRed);
  


  TCanvas *c1  = new TCanvas("c1","",600,600);
  c1->cd()->SetGrid();

  SC->Draw();
  SC->GetXaxis()->SetTitle("Signal Efficiency"); 
  SC->GetYaxis()->SetTitle("Background Rejection");
  REL->Draw("same");
  

  c1->SaveAs("ROCcomP.png");
  c1->SaveAs("ROCcomP.C");
  c1->SaveAs("ROCcomP.root");


}

