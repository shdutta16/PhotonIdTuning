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


#include "tdrstyle.C"


void Plotter(){


  gROOT->ProcessLine(" .L tdrstyle.C");
  setTDRStyle();
  gStyle->SetOptStat(0);

  TH1F *Letas = new TH1F("Letas","Loose Cut Efficiency Eta",100,-3,3);
  TH1F *Metas = new TH1F("Metas","Medium Cut Efficiency Eta",100,-3,3);
  TH1F *Tetas = new TH1F("Tetas","Tight Cut Efficiency Eta",100,-3,3); 

  TH1F *Lpts = new TH1F("Lpts","Loose Cut Efficiency pt",100,0,200);
  TH1F *Mpts = new TH1F("Mpts","Medium Cut Efficiency pt",100,0,200);
  TH1F *Tpts = new TH1F("Tpts","Tight Cut Efficiency pt",100,0,200); 

  TH1F *Lnvtxs = new TH1F("Lnvtxs","Loose Cut Efficiency vertices",100,0,100);
  TH1F *Mnvtxs = new TH1F("Mnvtxs","Medium Cut Efficiency vertices",100,0,100);
  TH1F *Tnvtxs = new TH1F("Tnvtxs","Tight Cut Efficiency vertices",100,0,100); 


  TH1F *Letab = new TH1F("Letab","Loose Cut b Efficiency Eta",100,-3,3);
  TH1F *Metab = new TH1F("Metab","Medium Cut b Efficiency Eta",100,-3,3);
  TH1F *Tetab = new TH1F("Tetab","Tight Cut b Efficiency Eta",100,-3,3); 


  TH1F *Lptb = new TH1F("Lptb","Loose Cut b Efficiency pt",100,0,200);
  TH1F *Mptb = new TH1F("Mptb","Medium Cut b Efficiency pt",100,0,200);
  TH1F *Tptb = new TH1F("Tptb","Tight Cut b Efficiency pt",100,0,200); 


  
  TH1F *Lnvtxb = new TH1F("Lnvtxb","Loose Cut b  Efficiency vertices",100,0,100);
  TH1F *Mnvtxb = new TH1F("Mnvtxb","Medium Cut b Efficiency vertices",100,0,100);
  TH1F *Tnvtxb = new TH1F("Tnvtxb","Tight Cut b Efficiency vertices",100,0,100); 


  TH1F *Sieaft  = new TH1F("Sieaft","Sieie cut only",100,0,200); 
  TH1F *ToEaft  = new TH1F("ToEaft","HoE cut only",100,0,200); 
  TH1F *IsoPaft = new TH1F("IsoPaft","IsoP cut only",100,0,200); 
  TH1F *IsoCaft = new TH1F("IsoCaft","IsoC cut only",100,0,200); 
  TH1F *IsoNaft = new TH1F("IsoNaft","IsoN cut only",100,0,200); 



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
  


  TH1F *EffPT0 = (TH1F*)input->Get("EffPT0");
  TH1F *EffPTL = (TH1F*)input->Get("EffPTL");
  TH1F *EffPTM = (TH1F*)input->Get("EffPTM");
  TH1F *EffPTT = (TH1F*)input->Get("EffPTT");

  TH1F *EffPT0b = (TH1F*)input->Get("EffPT0b");
  TH1F *EffPTLb = (TH1F*)input->Get("EffPTLb");
  TH1F *EffPTMb = (TH1F*)input->Get("EffPTMb");
  TH1F *EffPTTb = (TH1F*)input->Get("EffPTTb");



  TH1F *EffNVTX0 = (TH1F*)input->Get("EffNVTX0");
  TH1F *EffNVTXL = (TH1F*)input->Get("EffNVTXL");
  TH1F *EffNVTXM = (TH1F*)input->Get("EffNVTXM");
  TH1F *EffNVTXT = (TH1F*)input->Get("EffNVTXT");

  TH1F *EffNVTX0b = (TH1F*)input->Get("EffNVTX0b");
  TH1F *EffNVTXLb = (TH1F*)input->Get("EffNVTXLb");
  TH1F *EffNVTXMb = (TH1F*)input->Get("EffNVTXMb");
  TH1F *EffNVTXTb = (TH1F*)input->Get("EffNVTXTb");




  Letas->Divide(EffETAL,EffETA0,1.,1.,"B");
  Metas->Divide(EffETAM,EffETA0,1.,1.,"B");
  Tetas->Divide(EffETAT,EffETA0,1.,1.,"B");
  
  Letab->Divide(EffETALb,EffETA0b,1.,1.,"B");
  Metab->Divide(EffETAMb,EffETA0b,1.,1.,"B");
  Tetab->Divide(EffETATb,EffETA0b,1.,1.,"B");
  


  Lpts->Divide(EffPTL,EffPT0,1.,1.,"B");
  Mpts->Divide(EffPTM,EffPT0,1.,1.,"B");
  Tpts->Divide(EffPTT,EffPT0,1.,1.,"B");
  
  Lptb->Divide(EffPTLb,EffPT0b,1.,1.,"B");
  Mptb->Divide(EffPTMb,EffPT0b,1.,1.,"B");
  Tptb->Divide(EffPTTb,EffPT0b,1.,1.,"B");



  Lnvtxs->Divide(EffNVTXL,EffNVTX0,1.,1.,"B");
  Mnvtxs->Divide(EffNVTXM,EffNVTX0,1.,1.,"B");
  Tnvtxs->Divide(EffNVTXT,EffNVTX0,1.,1.,"B");
  
  Lnvtxb->Divide(EffNVTXL,EffNVTX0b,1.,1.,"B");
  Mnvtxb->Divide(EffNVTXM,EffNVTX0b,1.,1.,"B");
  Tnvtxb->Divide(EffNVTXT,EffNVTX0b,1.,1.,"B");



  TH1F *EffPTs = (TH1F*)input->Get("EffPTs");
  TH1F *EffPTt = (TH1F*)input->Get("EffPTt");
  TH1F *EffPTp = (TH1F*)input->Get("EffPTp");
  TH1F *EffPTc = (TH1F*)input->Get("EffPTc");
  TH1F *EffPTn = (TH1F*)input->Get("EffPTn");



  // the branch  out cuts 


  Sieaft->Divide(EffPTs,EffPT0,1.,1.,"B"); 
  ToEaft->Divide(EffPTt,EffPT0,1.,1.,"B");
  IsoPaft->Divide(EffPTp,EffPT0,1.,1.,"B");
  IsoCaft->Divide(EffPTc,EffPT0,1.,1.,"B");
  IsoNaft->Divide(EffPTn,EffPT0,1.,1.,"B");




  
  TCanvas *c1  = new TCanvas("c1","Medium",600,600);
  c1->Divide(2,2);

  c1->cd(1);
  Mnvtxs->Draw();
  Mnvtxs->GetYaxis()->SetRangeUser(0,1.0);
  Mnvtxs->GetXaxis()->SetRangeUser(0,50);
  Mnvtxs->GetXaxis()->SetTitle("# Nvtx");
  Mnvtxs->SetLineColor(kRed);
  Mnvtxs->SetMarkerColor(kRed);
  Mnvtxs->SetMarkerSize(0.5);
  Mnvtxb->SetMarkerSize(0.5);
  Mnvtxb->Draw("same");
  
  c1->cd(2);
  Mpts->Draw();
  Mpts->GetYaxis()->SetRangeUser(0,1.0);
  Mpts->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  Mpts->SetLineColor(kRed);
  Mpts->SetMarkerColor(kRed);
  Mpts->SetMarkerSize(0.5);
  Mptb->SetMarkerSize(0.5);
  Mptb->Draw("same");

  c1->cd(3);
  Metas->Draw();
  Metas->GetYaxis()->SetRangeUser(0,1.0);
  Metas->GetXaxis()->SetRangeUser(1.5,3.0);
  Metas->GetXaxis()->SetTitle("#eta");
  Metas->SetLineColor(kRed);
  Metas->SetMarkerColor(kRed);
  Metas->SetMarkerSize(0.5);
  Metab->SetMarkerSize(0.5);
  Metab->Draw("same");
  

  c1->SaveAs("MediumEffBck.png");
  c1->SaveAs("MediumEffBck.C");
  c1->SaveAs("MediumEffBck.root");



  TCanvas *c2  = new TCanvas("c2","Loose",600,600);
  c2->Divide(2,2);

  c2->cd(1);
  Lnvtxs->Draw();
  Lnvtxs->GetXaxis()->SetTitle("# Nvtx");
  Lnvtxs->GetYaxis()->SetRangeUser(0,1.0);
  Lnvtxs->GetXaxis()->SetRangeUser(0,50);
  Lnvtxs->SetLineColor(kRed);
  Lnvtxs->SetMarkerColor(kRed);
  Lnvtxs->SetMarkerSize(0.5);
  Lnvtxb->SetMarkerSize(0.5);
  Lnvtxb->Draw("esame");
  
  c2->cd(2);
  Lpts->Draw();
  Lpts->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  Lpts->GetYaxis()->SetRangeUser(0,1.0);
  Lpts->SetLineColor(kRed);
  Lpts->SetMarkerColor(kRed);
  Lpts->SetMarkerSize(0.5);
  Lptb->SetMarkerSize(0.5);
  Lptb->Draw("esame");

  c2->cd(3);
  Letas->Draw();
  Letas->GetXaxis()->SetTitle("#eta");
  Letas->GetYaxis()->SetRangeUser(0,1.0);
  Letas->GetXaxis()->SetRangeUser(1.5,3.0);
  Letas->SetLineColor(kRed);
  Letas->SetMarkerColor(kRed);
  Letas->SetMarkerSize(0.5);
  Letab->SetMarkerSize(0.5);
  Letab->Draw("esame");

  
  c2->SaveAs("LooseEffBck.png");
  c2->SaveAs("LooseEffBck.C");
  c2->SaveAs("LooseEffBck.root");




  TCanvas *c3  = new TCanvas("c3","Tight",600,600);
  c3->Divide(2,2);

  c3->cd(1);
  Tnvtxs->Draw();
  Tnvtxs->GetYaxis()->SetRangeUser(0,1.0);
  Tnvtxs->GetXaxis()->SetRangeUser(0,50);
  Tnvtxs->GetXaxis()->SetTitle("# Nvtx");
  Tnvtxs->SetLineColor(kRed);
  Tnvtxs->SetMarkerColor(kRed);
  Tnvtxs->SetMarkerSize(0.5);
  Tnvtxb->SetMarkerSize(0.5);
  Tnvtxb->Draw("esame");
  
  c3->cd(2);
  Tpts->Draw();
  Tpts->GetYaxis()->SetRangeUser(0,1.0);
  Tpts->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  Tpts->SetLineColor(kRed);
  Tpts->SetMarkerColor(kRed);
  Tpts->SetMarkerSize(0.5);
  Tptb->SetMarkerSize(0.5);
  Tptb->Draw("esame");

  c3->cd(3);
  Tetas->Draw();
  Tetas->GetYaxis()->SetRangeUser(0,1.0);
  Tetas->GetXaxis()->SetTitle("#eta");
  Tetas->GetXaxis()->SetRangeUser(1.5,3.0);
  Tetas->SetLineColor(kRed);
  Tetas->SetMarkerColor(kRed);
  Tetas->SetMarkerSize(0.5);
  Tetab->SetMarkerSize(0.5);
  Tetab->Draw("same");


  c3->SaveAs("TightEffBck.png");
  c3->SaveAs("TightEffBck.C");
  c3->SaveAs("TightEffBck.root");

 


  TCanvas *c10 = new TCanvas("c10","Branch Out Cuts",900,900);
  c10->Divide(3,2);
  
  c10->cd(1);  
  Sieaft->Draw();
  Sieaft->SetMarkerSize(0.5);
  Sieaft->GetYaxis()->SetTitle("Only Sieie Cut Efficiency");
  Sieaft->GetXaxis()->SetTitle("Pt GeVc^{-1}");
 
  c10->cd(2);
  ToEaft->Draw();
  Sieaft->SetMarkerSize(0.5);
  ToEaft->GetYaxis()->SetTitle("Only HOE Cut Efficiency");
  ToEaft->GetXaxis()->SetTitle("Pt GeVc^{-1}");

  c10->cd(3);
  IsoPaft->Draw();
  IsoPaft->SetMarkerSize(0.5);
  IsoPaft->GetYaxis()->SetTitle("Only iso p Cut Efficiency");
  IsoPaft->GetXaxis()->SetTitle("Pt GeVc^{-1}");

  c10->cd(4);
  IsoCaft->Draw();
  IsoCaft->SetMarkerSize(0.5);
  IsoCaft->GetYaxis()->SetTitle("Only iso c Cut Efficiency");
  IsoCaft->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  
  c10->cd(5);
  IsoNaft->Draw();
  IsoNaft->SetMarkerSize(0.5);
  IsoNaft->GetYaxis()->SetTitle("Only iso n Cut Efficiency");
  IsoNaft->GetXaxis()->SetTitle("Pt GeVc^{-1}");
  
  c10->SaveAs("BranchOutCuts.png");
  c10->SaveAs("BranchOutCuts.C");
  c10->SaveAs("BranchOutCuts.root");


}
