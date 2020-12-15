#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"


#include <sstream>
#include <string.h>
#include <fstream>

void Etas(){


  TH1F *Letas = new TH1F("Letas","Loose Cut Efficiency Eta",100,-3,3);
  TH1F *Metas = new TH1F("Metas","Medium Cut Efficiency Eta",100,-3,3);
  TH1F *Tetas = new TH1F("Tetas","Tight Cut Efficiency Eta",100,-3,3);

  TH1F *Letab = new TH1F("Letab","Loose Cut b Efficiency Eta",100,-3,3);
  TH1F *Metab = new TH1F("Metab","Medium Cut b Efficiency Eta",100,-3,3);
  TH1F *Tetab = new TH1F("Tetab","Tight Cut b Efficiency Eta",100,-3,3);


  TString fname = "Eff.root";
  TFile *input = TFile::Open( fname );
  if (!input || !input->IsOpen()) {
    cout << "\nERROR! Could not open root file " << fname
         << endl;
    exit(0);
  }

  TH1F *EffETA0 = (TH1F*)input->Get("EffETA0");
  TH1F *EffETAL = (TH1F*)input->Get("EffETAL");
  TH1F *EffETAM = (TH1F*)input->Get("EffETAM");
  TH1F *EffETAT = (TH1F*)input->Get("EffETAT");

  TH1F *EffETA0b = (TH1F*)input->Get("EffETA0b");
  TH1F *EffETALb = (TH1F*)input->Get("EffETALb");
  TH1F *EffETAMb = (TH1F*)input->Get("EffETAMb");
  TH1F *EffETATb = (TH1F*)input->Get("EffETATb");



  Letas->Divide(EffETAL,EffETA0,1.,1.,"B");
  Metas->Divide(EffETAM,EffETA0,1.,1.,"B");
  Tetas->Divide(EffETAT,EffETA0,1.,1.,"B");

  Letab->Divide(EffETALb,EffETA0b,1.,1.,"B");
  Metab->Divide(EffETAMb,EffETA0b,1.,1.,"B");
  Tetab->Divide(EffETATb,EffETA0b,1.,1.,"B");



  TCanvas *ceta = new TCanvas("ceta","Eta Eff",500,500);
  ceta->cd();
  Letas->GetXaxis()->SetTitle("#eta");
  Letas->GetXaxis()->SetRangeUser(-3,3);                                                                
  //  Letas->GetYaxis()->SetRangeUser(-5,5);
  Letas->SetLineColor(kGray + 3);
  Letas->SetMarkerColor(kGray +3);
  Letas->SetMarkerStyle(20);
  Letas->Draw();

  Letab->SetLineColor(kAzure + 3);
  Letab->SetMarkerColor(kAzure +3);
  Letab->SetMarkerStyle(20);
  Letab->Draw("esame");
  
  Metab->SetLineColor(kAzure + 5);
  Metab->SetMarkerColor(kAzure +5);
  Metab->SetMarkerStyle(20);
  Metab->Draw("esame");

  Tetab->SetLineColor(kAzure + 10);
  Tetab->SetMarkerColor(kAzure +10);
  Tetab->SetMarkerStyle(20);
  Tetab->Draw("esame");

  Metas->Draw("esame");
  Metas->SetMarkerColor(kOrange -3);
  Metas->SetLineColor(kOrange -3);
  Metas->SetMarkerStyle(20);

  Tetas->Draw("esame");
  Tetas->SetMarkerColor(kYellow );
  Tetas->SetLineColor(kYellow );
  Tetas->SetMarkerStyle(20);
  ceta->Update();

  ceta->SaveAs("EfETA_all.png");
  ceta->SaveAs("EfETA_all.C");
  ceta->SaveAs("EfETA_all.root");




}
