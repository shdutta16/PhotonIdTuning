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

void PlotBuilder(){


  //Read the cuts: 


  ifstream tight;
  ifstream medium;
  ifstream loose;

  TString tight_name = "TightR.txt";
  TString medium_name = "MediumR.txt";
  TString loose_name = "LooseR.txt";

  tight.open(tight_name);
  if( !tight || !tight.is_open() ){
    cout << "\nERROR! Could not open txt file " 
	 << tight_name << endl;
    exit(0);
  }

  medium.open(medium_name);
  if( !medium || !medium.is_open() ){
    cout << "\nERROR! Could not open txt file "
         << medium_name << endl;
    exit(0);
  }

  loose.open(loose_name);
  if( !loose || !loose.is_open() ){
    cout << "\nERROR! Could not open txt file "
         << loose_name << endl;
    exit(0);
  }


  double isoCL,isoCM,isoCT;
  double isoPL,isoPM,isoPT;
  double isoNL,isoNM,isoNT;
  double sieL,sieM,sieT; 
  double toeL,toeM,toeT;

  if(tight.is_open()){
    while(!tight.eof()){
      tight>>sieT>>isoCT>>isoNT>>isoPT>>toeT;
      break;
    }
  }

  
  if(medium.is_open()){
    while(!medium.eof()){
      medium>>sieM>>isoCM>>isoNM>>isoPM>>toeM;
      break;
    }
  }


  if(loose.is_open()){
    while(!loose.eof()){
      loose>>sieL>>isoCL>>isoNL>>isoPL>>toeL;
      break;
    }
  }



  // PT

  // Signal
  TH1F *EffPT0   = new TH1F("EffPT0","Signal Eff vs PT 0",100,0,250);
  TH1F *EffPTL   = new TH1F("EffPTL","Signal Eff vs PT L",100,0,250);
  TH1F *EffPTM   = new TH1F("EffPTM","Signal Eff vs PT M",100,0,250);
  TH1F *EffPTT   = new TH1F("EffPTT","Signal Eff vs PT T",100,0,250);

  // Background
  TH1F *EffPT0b   = new TH1F("EffPT0b","Background Eff vs PT 0",100,0,250);
  TH1F *EffPTLb   = new TH1F("EffPTLb","Background Eff vs PT L",100,0,250);
  TH1F *EffPTMb   = new TH1F("EffPTMb","Background Eff vs PT M",100,0,250);
  TH1F *EffPTTb   = new TH1F("EffPTTb","Background Eff vs PT T",100,0,250);



  // NVTX

  // Signal
  TH1F *EffNVTX0 = new TH1F("EffNVTX0","Signal Eff vs NVTX 0",180,0,180);
  TH1F *EffNVTXL = new TH1F("EffNVTXL","Signal Eff vs NVTX L",180,0,180);
  TH1F *EffNVTXM = new TH1F("EffNVTXM","Signal Eff vs NVTX M",180,0,180);
  TH1F *EffNVTXT = new TH1F("EffNVTXT","Signal Eff vs NVTX T",180,0,180);

  // Background
  TH1F *EffNVTX0b = new TH1F("EffNVTX0b","Background Eff vs NVTX 0",180,0,180);
  TH1F *EffNVTXLb = new TH1F("EffNVTXLb","Background Eff vs NVTX L",180,0,180);
  TH1F *EffNVTXMb = new TH1F("EffNVTXMb","Background Eff vs NVTX M",180,0,180);
  TH1F *EffNVTXTb = new TH1F("EffNVTXTb","Background Eff vs NVTX T",180,0,180);



  // ETA

  // Signal
  TH1F *EffETA0  = new TH1F("EffETA0","Signal Eff vs Eta 0",100,-3,3);
  TH1F *EffETAL  = new TH1F("EffETAL","Signal Eff vs Eta L",100,-3,3);
  TH1F *EffETAM  = new TH1F("EffETAM","Signal Eff vs Eta M",100,-3,3);
  TH1F *EffETAT  = new TH1F("EffETAT","Signal Eff vs Eta T",100,-3,3);

  // Background
  TH1F *EffETA0b  = new TH1F("EffETA0b","Background Eff vs Eta 0",100,-3,3);
  TH1F *EffETALb  = new TH1F("EffETALb","Background Eff vs Eta L",100,-3,3);
  TH1F *EffETAMb  = new TH1F("EffETAMb","Background Eff vs Eta M",100,-3,3);
  TH1F *EffETATb  = new TH1F("EffETATb","Background Eff vs Eta T",100,-3,3);



  // PHI

  // Signal
  TH1F *EffPHI0  = new TH1F("EffPHI0","Signal Eff vs PHI 0",100,-4,4);
  TH1F *EffPHIL  = new TH1F("EffPHIL","Signal Eff vs PHI L",100,-4,4);
  TH1F *EffPHIM  = new TH1F("EffPHIM","Signal Eff vs PHI M",100,-4,4);
  TH1F *EffPHIT  = new TH1F("EffPHIT","Signal Eff vs PHI T",100,-4,4);

  // Background
  TH1F *EffPHI0b  = new TH1F("EffPHI0b","Background Eff vs PHI 0",100,-4,4);
  TH1F *EffPHILb  = new TH1F("EffPHILb","Background Eff vs PHI L",100,-4,4);
  TH1F *EffPHIMb  = new TH1F("EffPHIMb","Background Eff vs PHI M",100,-4,4);
  TH1F *EffPHITb  = new TH1F("EffPHITb","Background Eff vs PHI T",100,-4,4);



  // Branch out Cuts


  // PT

  // Signal
  TH1F *EffPTs  = new TH1F("EffPTs","Signal Eff vs PT L s",100,0,250);
  TH1F *EffPTt  = new TH1F("EffPTt","Signal Eff vs PT L t",100,0,250);
  TH1F *EffPTc  = new TH1F("EffPTc","Signal Eff vs PT L c",100,0,250);
  TH1F *EffPTp  = new TH1F("EffPTp","Signal Eff vs PT L p",100,0,250);
  TH1F *EffPTn  = new TH1F("EffPTn","Signal Eff vs PT L n",100,0,250);

  TH1F *EffPTMs  = new TH1F("EffPTMs","Signal Eff vs PT M s",100,0,250);
  TH1F *EffPTMt  = new TH1F("EffPTMt","Signal Eff vs PT M t",100,0,250);
  TH1F *EffPTMc  = new TH1F("EffPTMc","Signal Eff vs PT M c",100,0,250);
  TH1F *EffPTMp  = new TH1F("EffPTMp","Signal Eff vs PT M n",100,0,250);
  TH1F *EffPTMn  = new TH1F("EffPTMn","Signal Eff vs PT M p",100,0,250);

  TH1F *EffPTTs  = new TH1F("EffPTTs","Signal Eff vs PT T s",100,0,250);
  TH1F *EffPTTt  = new TH1F("EffPTTt","Signal Eff vs PT T t",100,0,250);
  TH1F *EffPTTc  = new TH1F("EffPTTc","Signal Eff vs PT T c",100,0,250);
  TH1F *EffPTTp  = new TH1F("EffPTTp","Signal Eff vs PT T n",100,0,250);
  TH1F *EffPTTn  = new TH1F("EffPTTn","Signal Eff vs PT T p",100,0,250);

  // Background
  TH1F *EffPTbs  = new TH1F("EffPTbs","Background Eff vs PT L s",100,0,250);
  TH1F *EffPTbt  = new TH1F("EffPTbt","Background Eff vs PT L t",100,0,250);
  TH1F *EffPTbc  = new TH1F("EffPTbc","Background Eff vs PT L c",100,0,250);
  TH1F *EffPTbp  = new TH1F("EffPTbp","Background Eff vs PT L p",100,0,250);
  TH1F *EffPTbn  = new TH1F("EffPTbn","Background Eff vs PT L n",100,0,250);

  TH1F *EffPTMbs  = new TH1F("EffPTMbs","Background Eff vs PT M s",100,0,250);
  TH1F *EffPTMbt  = new TH1F("EffPTMbt","Background Eff vs PT M t",100,0,250);
  TH1F *EffPTMbc  = new TH1F("EffPTMbc","Background Eff vs PT M c",100,0,250);
  TH1F *EffPTMbp  = new TH1F("EffPTMbp","Background Eff vs PT M p",100,0,250);
  TH1F *EffPTMbn  = new TH1F("EffPTMbn","Background Eff vs PT M n",100,0,250);

  TH1F *EffPTTbs  = new TH1F("EffPTTbs","Background Eff vs PT T s",100,0,250);
  TH1F *EffPTTbt  = new TH1F("EffPTTbt","Background Eff vs PT T t",100,0,250);
  TH1F *EffPTTbc  = new TH1F("EffPTTbc","Background Eff vs PT T c",100,0,250);
  TH1F *EffPTTbp  = new TH1F("EffPTTbp","Background Eff vs PT T p",100,0,250);
  TH1F *EffPTTbn  = new TH1F("EffPTTbn","Background Eff vs PT T n",100,0,250);



  // NVTX

  // Signal
  TH1F *EffNVTXs = new TH1F("EffNVTXs","Signal Eff vs NVTX L s",180,0,180);
  TH1F *EffNVTXt = new TH1F("EffNVTXt","Signal Eff vs NVTX L t",180,0,180);
  TH1F *EffNVTXc = new TH1F("EffNVTXc","Signal Eff vs NVTX L c",180,0,180);
  TH1F *EffNVTXp = new TH1F("EffNVTXp","Signal Eff vs NVTX L p",180,0,180);
  TH1F *EffNVTXn = new TH1F("EffNVTXn","Signal Eff vs NVTX L n",180,0,180);

  TH1F *EffNVTXMs  = new TH1F("EffNVTXMs","Signal Eff vs NVTX M s",180,0,180);
  TH1F *EffNVTXMt  = new TH1F("EffNVTXMt","Signal Eff vs NVTX M t",180,0,180);
  TH1F *EffNVTXMc  = new TH1F("EffNVTXMc","Signal Eff vs NVTX M c",180,0,180);
  TH1F *EffNVTXMp  = new TH1F("EffNVTXMp","Signal Eff vs NVTX M n",180,0,180);
  TH1F *EffNVTXMn  = new TH1F("EffNVTXMn","Signal Eff vs NVTX M p",180,0,180);

  TH1F *EffNVTXTs  = new TH1F("EffNVTXTs","Signal Eff vs NVTX T s",180,0,180);
  TH1F *EffNVTXTt  = new TH1F("EffNVTXTt","Signal Eff vs NVTX T t",180,0,180);
  TH1F *EffNVTXTc  = new TH1F("EffNVTXTc","Signal Eff vs NVTX T c",180,0,180);
  TH1F *EffNVTXTp  = new TH1F("EffNVTXTp","Signal Eff vs NVTX T n",180,0,180);
  TH1F *EffNVTXTn  = new TH1F("EffNVTXTn","Signal Eff vs NVTX T p",180,0,180);


  // Background
  TH1F *EffNVTXbs = new TH1F("EffNVTXbs","Background Eff vs NVTX L s",180,0,180);
  TH1F *EffNVTXbt = new TH1F("EffNVTXbt","Background Eff vs NVTX L t",180,0,180);
  TH1F *EffNVTXbc = new TH1F("EffNVTXbc","Background Eff vs NVTX L c",180,0,180);
  TH1F *EffNVTXbn = new TH1F("EffNVTXbn","Background Eff vs NVTX L n",180,0,180);
  TH1F *EffNVTXbp = new TH1F("EffNVTXbp","Background Eff vs NVTX L p",180,0,180);

  TH1F *EffNVTXMbs = new TH1F("EffNVTXMbs","Background Eff vs NVTX M s",180,0,180);
  TH1F *EffNVTXMbt = new TH1F("EffNVTXMbt","Background Eff vs NVTX M t",180,0,180);
  TH1F *EffNVTXMbc = new TH1F("EffNVTXMbc","Background Eff vs NVTX M c",180,0,180);
  TH1F *EffNVTXMbn = new TH1F("EffNVTXMbn","Background Eff vs NVTX M n",180,0,180);
  TH1F *EffNVTXMbp = new TH1F("EffNVTXMbp","Background Eff vs NVTX M p",180,0,180);

  TH1F *EffNVTXTbs = new TH1F("EffNVTXTbs","Background Eff vs NVTX T s",180,0,180);
  TH1F *EffNVTXTbt = new TH1F("EffNVTXTbt","Background Eff vs NVTX T t",180,0,180);
  TH1F *EffNVTXTbc = new TH1F("EffNVTXTbc","Background Eff vs NVTX T c",180,0,180);
  TH1F *EffNVTXTbn = new TH1F("EffNVTXTbn","Background Eff vs NVTX T n",180,0,180);
  TH1F *EffNVTXTbp = new TH1F("EffNVTXTbp","Background Eff vs NVTX T p",180,0,180);



  // ETA

  // Signal
  TH1F *EffETAs  = new TH1F("EffETAs","Signal Eff vs ETA L s",100,-3,3);
  TH1F *EffETAt  = new TH1F("EffETAt","Signal Eff vs ETA L t",100,-3,3);
  TH1F *EffETAc  = new TH1F("EffETAc","Signal Eff vs ETA L c",100,-3,3);
  TH1F *EffETAn  = new TH1F("EffETAn","Signal Eff vs ETA L n",100,-3,3);
  TH1F *EffETAp  = new TH1F("EffETAp","Signal Eff vs ETA L p",100,-3,3);

  TH1F *EffETAMs  = new TH1F("EffETAMs","Signal Eff vs ETA M s",100,-3,3);
  TH1F *EffETAMt  = new TH1F("EffETAMt","Signal Eff vs ETA M t",100,-3,3);
  TH1F *EffETAMc  = new TH1F("EffETAMc","Signal Eff vs ETA M c",100,-3,3);
  TH1F *EffETAMp  = new TH1F("EffETAMp","Signal Eff vs ETA M n",100,-3,3);
  TH1F *EffETAMn  = new TH1F("EffETAMn","Signal Eff vs ETA M p",100,-3,3);

  TH1F *EffETATs  = new TH1F("EffETATs","Signal Eff vs ETA T s",100,-3,3);
  TH1F *EffETATt  = new TH1F("EffETATt","Signal Eff vs ETA T t",100,-3,3);
  TH1F *EffETATc  = new TH1F("EffETATc","Signal Eff vs ETA T c",100,-3,3);
  TH1F *EffETATp  = new TH1F("EffETATp","Signal Eff vs ETA T n",100,-3,3);
  TH1F *EffETATn  = new TH1F("EffETATn","Signal Eff vs ETA T p",100,-3,3);


  // Background
  TH1F *EffETAbs  = new TH1F("EffETAbs","Background Eff vs ETA L s",100,-3,3);
  TH1F *EffETAbt  = new TH1F("EffETAbt","Background Eff vs ETA L t",100,-3,3);
  TH1F *EffETAbc  = new TH1F("EffETAbc","Background Eff vs ETA L c",100,-3,3);
  TH1F *EffETAbn  = new TH1F("EffETAbn","Background Eff vs ETA L n",100,-3,3);
  TH1F *EffETAbp  = new TH1F("EffETAbp","Background Eff vs ETA L p",100,-3,3);

  TH1F *EffETAMbs  = new TH1F("EffETAMbs","Background Eff vs ETA M s",100,-3,3);
  TH1F *EffETAMbt  = new TH1F("EffETAMbt","Background Eff vs ETA M t",100,-3,3);
  TH1F *EffETAMbc  = new TH1F("EffETAMbc","Background Eff vs ETA M c",100,-3,3);
  TH1F *EffETAMbn  = new TH1F("EffETAMbn","Background Eff vs ETA M n",100,-3,3);
  TH1F *EffETAMbp  = new TH1F("EffETAMbp","Background Eff vs ETA M p",100,-3,3);

  TH1F *EffETATbs  = new TH1F("EffETATbs","Background Eff vs ETA T s",100,-3,3);
  TH1F *EffETATbt  = new TH1F("EffETATbt","Background Eff vs ETA T t",100,-3,3);
  TH1F *EffETATbc  = new TH1F("EffETATbc","Background Eff vs ETA T c",100,-3,3);
  TH1F *EffETATbn  = new TH1F("EffETATbn","Background Eff vs ETA T n",100,-3,3);
  TH1F *EffETATbp  = new TH1F("EffETATbp","Background Eff vs EAT T p",100,-3,3);


  
  //Setting the Tree Branches

  TString finput_name = "../../CutTMVABarrel90_test.root"; 
  TFile *finput = new TFile(finput_name);
  if (!finput || !finput->IsOpen() ) {
    cout << "\nERROR! Could not open root file " << finput_name
         << endl;
    exit(0);
  }


  float Ppt,Peta,Pphi,isoP,isoC,isoN,sieie,toe,weight;
  int   nvtx; 
  
  //  weight = 1.0; 

  finput->cd();

  //Signal Tree
  TTree *t_S = (TTree*)finput->Get("t_S");
  
  t_S->SetBranchAddress("Sieie",&sieie);
  t_S->SetBranchAddress("isoP",&isoP);
  t_S->SetBranchAddress("isoC",&isoC);
  t_S->SetBranchAddress("isoN",&isoN);
  t_S->SetBranchAddress("ToE",&toe);
  t_S->SetBranchAddress("weighT",&weight);
  t_S->SetBranchAddress("Nvtx",&nvtx);
  t_S->SetBranchAddress("Peta",&Peta);
  t_S->SetBranchAddress("Ppt",&Ppt);

  //Background Tree
  TTree *t_B = (TTree*)finput->Get("t_B");
                                          
  t_B->SetBranchAddress("Sieie",&sieie);
  t_B->SetBranchAddress("isoP",&isoP);
  t_B->SetBranchAddress("isoC",&isoC);
  t_B->SetBranchAddress("isoN",&isoN);
  t_B->SetBranchAddress("ToE",&toe);
  t_B->SetBranchAddress("weighT",&weight);
  t_B->SetBranchAddress("Nvtx",&nvtx);
  t_B->SetBranchAddress("Peta",&Peta);
  t_B->SetBranchAddress("Ppt",&Ppt);



  EffPT0->Sumw2();
  EffPTL->Sumw2();
  EffPTM->Sumw2();
  EffPTT->Sumw2();

  EffPT0b->Sumw2();
  EffPTLb->Sumw2();
  EffPTMb->Sumw2();
  EffPTTb->Sumw2();


  EffNVTX0->Sumw2();
  EffNVTXL->Sumw2();
  EffNVTXM->Sumw2();
  EffNVTXT->Sumw2();

  EffNVTX0b->Sumw2();
  EffNVTXLb->Sumw2();
  EffNVTXMb->Sumw2();
  EffNVTXTb->Sumw2();


  EffETA0->Sumw2();
  EffETAL->Sumw2();
  EffETAM->Sumw2();
  EffETAT->Sumw2();

  EffETA0b->Sumw2();
  EffETALb->Sumw2();
  EffETAMb->Sumw2();
  EffETATb->Sumw2();


  EffPHI0->Sumw2();
  EffPHIL->Sumw2();
  EffPHIM->Sumw2();
  EffPHIT->Sumw2();
  
  EffPHI0b->Sumw2();
  EffPHILb->Sumw2();
  EffPHIMb->Sumw2();
  EffPHITb->Sumw2();

  

  EffPTs->Sumw2();
  EffPTt->Sumw2();
  EffPTc->Sumw2();
  EffPTn->Sumw2();
  EffPTp->Sumw2();

  EffPTMs->Sumw2();
  EffPTMt->Sumw2();
  EffPTMc->Sumw2();
  EffPTMp->Sumw2();
  EffPTMn->Sumw2();

  EffPTTs->Sumw2();
  EffPTTt->Sumw2();
  EffPTTc->Sumw2();
  EffPTTp->Sumw2();
  EffPTTn->Sumw2();


  EffPTbs->Sumw2();
  EffPTbt->Sumw2();
  EffPTbc->Sumw2();
  EffPTbn->Sumw2();
  EffPTbp->Sumw2();

  EffPTMbs->Sumw2();
  EffPTMbt->Sumw2();
  EffPTMbc->Sumw2();
  EffPTMbn->Sumw2();
  EffPTMbp->Sumw2();

  EffPTTbs->Sumw2();
  EffPTTbt->Sumw2();
  EffPTTbc->Sumw2();
  EffPTTbn->Sumw2();
  EffPTTbp->Sumw2();



  EffNVTXs->Sumw2();
  EffNVTXt->Sumw2();
  EffNVTXc->Sumw2();
  EffNVTXn->Sumw2();
  EffNVTXp->Sumw2();

  EffNVTXMs->Sumw2();
  EffNVTXMt->Sumw2();
  EffNVTXMc->Sumw2();
  EffNVTXMp->Sumw2();
  EffNVTXMn->Sumw2();

  EffNVTXTs->Sumw2();
  EffNVTXTt->Sumw2();
  EffNVTXTc->Sumw2();
  EffNVTXTp->Sumw2();
  EffNVTXTn->Sumw2();


  EffNVTXbs->Sumw2();
  EffNVTXbt->Sumw2();
  EffNVTXbc->Sumw2();
  EffNVTXbn->Sumw2();
  EffNVTXbp->Sumw2();

  EffNVTXMbs->Sumw2();
  EffNVTXMbt->Sumw2();
  EffNVTXMbc->Sumw2();
  EffNVTXMbn->Sumw2();
  EffNVTXMbp->Sumw2();

  EffNVTXTbs->Sumw2();
  EffNVTXTbt->Sumw2();
  EffNVTXTbc->Sumw2();
  EffNVTXTbn->Sumw2();
  EffNVTXTbp->Sumw2();



  EffETAs->Sumw2();
  EffETAt->Sumw2();
  EffETAc->Sumw2();
  EffETAn->Sumw2();
  EffETAp->Sumw2();

  EffETAMs->Sumw2();
  EffETAMt->Sumw2();
  EffETAMc->Sumw2();
  EffETAMp->Sumw2();
  EffETAMn->Sumw2();

  EffETATs->Sumw2();
  EffETATt->Sumw2();
  EffETATc->Sumw2();
  EffETATp->Sumw2();
  EffETATn->Sumw2();


  EffETAbs->Sumw2();
  EffETAbt->Sumw2();
  EffETAbc->Sumw2();
  EffETAbn->Sumw2();
  EffETAbp->Sumw2();

  EffETAMbs->Sumw2();
  EffETAMbt->Sumw2();
  EffETAMbc->Sumw2();
  EffETAMbn->Sumw2();
  EffETAMbp->Sumw2();

  EffETATbs->Sumw2();
  EffETATbt->Sumw2();
  EffETATbc->Sumw2();
  EffETATbn->Sumw2();
  EffETATbp->Sumw2();




  for(int i = 0; i < t_S->GetEntries(); i++){
    t_S->GetEntry(i);
    
    EffNVTX0->Fill(nvtx,weight);
    EffPT0->Fill(Ppt,weight);
    EffETA0->Fill(Peta,weight);

  
    //Loose Cut: 
    if((sieie  < sieL)&&
       (toe    < toeL)&&
       (isoP-0.002544*Ppt < isoPL)&&
       (isoC   < isoCL)&&
       (isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNL)){
      
      EffNVTXL->Fill(nvtx,weight);
      EffPTL->Fill(Ppt,weight);
      EffETAL->Fill(Peta,weight);
    }
    
    //Medium Cut:
    if((sieie  < sieM)&&
       (toe    < toeM)&&
       (isoP -0.002544*Ppt < isoPM)&&
       (isoC   < isoCM)&&
       (isoN- (0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNM)){
      EffNVTXM->Fill(nvtx,weight);
      EffPTM->Fill(Ppt,weight);
      EffETAM->Fill(Peta,weight);
    }

    //Tight Cut:
    if((sieie  < sieT)&&
       (toe    < toeT)&&
       (isoP -0.002544*Ppt < isoPT)&&
       (isoC   < isoCT)&&
       (isoN - (0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNT)){
      EffNVTXT->Fill(nvtx,weight);
      EffPTT->Fill(Ppt,weight);
      EffETAT->Fill(Peta,weight);
    }

    //Branch out Cuts Loose
    if(sieie  < sieL){
       EffNVTXs->Fill(nvtx,weight);
       EffPTs->Fill(Ppt,weight);
       EffETAs->Fill(Peta,weight);
       
    }
    if(toe    < toeL){
      EffNVTXt->Fill(nvtx,weight);
      EffPTt->Fill(Ppt,weight);
      EffETAt->Fill(Peta,weight);
     
    }    
    if(isoP -0.002544*Ppt < isoPL){    
      EffNVTXp->Fill(nvtx,weight);
      EffPTp->Fill(Ppt,weight);
      EffETAp->Fill(Peta,weight);
     
    }    
    if(isoC   < isoCL){
      EffNVTXc->Fill(nvtx,weight);
      EffPTc->Fill(Ppt,weight);
      EffETAc->Fill(Peta,weight);
     
    }
    if(isoN-(0.01556*Ppt - 0.000001129*Ppt*Ppt) < isoNL){
      EffNVTXn->Fill(nvtx,weight);
      EffPTn->Fill(Ppt,weight);
      EffETAn->Fill(Peta,weight);
     
    }


    //Branch out Cuts Medium
    if(sieie  < sieM){
       EffPTMs->Fill(Ppt,weight);
       EffNVTXMs->Fill(nvtx,weight);
       EffETAMs->Fill(Peta,weight);
    }
    if(toe    < toeM){
      EffPTMt->Fill(Ppt,weight);
      EffNVTXMt->Fill(nvtx,weight);
      EffETAMt->Fill(Peta,weight);
    }    
    if(isoP -0.002544*Ppt < isoPM){    
      EffPTMp->Fill(Ppt,weight);
      EffNVTXMp->Fill(nvtx,weight);
      EffETAMp->Fill(Peta,weight);
    }    
    if(isoC   < isoCM){
      EffPTMc->Fill(Ppt,weight);
      EffNVTXMc->Fill(nvtx,weight);
      EffETAMc->Fill(Peta,weight);
    }
    if(isoN-(0.01556*Ppt - 0.000001129*Ppt*Ppt) < isoNM){
      EffPTMn->Fill(Ppt,weight);
      EffNVTXMn->Fill(nvtx,weight);
      EffETAMn->Fill(Peta,weight);
    }


    //Branch out Cuts Tight
    if(sieie  < sieT){
       EffPTTs->Fill(Ppt,weight);
       EffNVTXTs->Fill(nvtx,weight);
       EffETATs->Fill(Peta,weight);
    }
    if(toe    < toeT){
      EffPTTt->Fill(Ppt,weight);
      EffNVTXTt->Fill(nvtx,weight);
      EffETATt->Fill(Peta,weight);
    }    
    if(isoP -0.002544*Ppt < isoPT){    
      EffPTTp->Fill(Ppt,weight);
      EffNVTXTp->Fill(nvtx,weight);
      EffETATp->Fill(Peta,weight);
    }    
    if(isoC   < isoCT){
      EffPTTc->Fill(Ppt,weight);
      EffNVTXTc->Fill(nvtx,weight);
      EffETATc->Fill(Peta,weight);
    }
    if(isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNT){
      EffPTTn->Fill(Ppt,weight);
      EffNVTXTn->Fill(nvtx,weight);
      EffETATn->Fill(Peta,weight);
    }


  }  



  //plots for Background 

  for(int i = 0; i < t_B->GetEntries(); i++){
    t_B->GetEntry(i);
    EffNVTX0b->Fill(nvtx,weight);
    EffPT0b->Fill(Ppt,weight);
    EffETA0b->Fill(Peta,weight);
  
    //Loose Cut: 
    if((sieie  < sieL)&&
       (toe    < toeL)&&
       (isoP -0.002544*Ppt < isoPL)&&
       (isoC   < isoCL)&&
       (isoN-((0.01556*Ppt-0.000001129*Ppt*Ppt)) < isoNL)){
      
      EffNVTXLb->Fill(nvtx,weight);
      EffPTLb->Fill(Ppt,weight);
      EffETALb->Fill(Peta,weight);
    }
    
    //Medium Cut:
    if((sieie  < sieM)&&
       (toe    < toeM)&&
       (isoP -0.002544*Ppt < isoPM)&&
       (isoC   < isoCM)&&
       (isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNM)){
      EffNVTXMb->Fill(nvtx,weight);
      EffPTMb->Fill(Ppt,weight);
      EffETAMb->Fill(Peta,weight);
    }

    //Tight Cut:
    if((sieie  < sieT)&&
       (toe    < toeT)&&
       (isoP -0.002544*Ppt < isoPT)&&
       (isoC   < isoCT)&&
       (isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNT)){
      EffNVTXTb->Fill(nvtx,weight);
      EffPTTb->Fill(Ppt,weight);
      EffETATb->Fill(Peta,weight);
    }


    //Branch out Cuts

    //Loose
    if(sieie  < sieL){
       EffNVTXbs->Fill(nvtx,weight);
       EffPTbs->Fill(Ppt,weight);
       EffETAbs->Fill(Peta,weight);
    }
    if(toe    < toeL){
      EffNVTXbt->Fill(nvtx,weight);
      EffPTbt->Fill(Ppt,weight);
      EffETAbt->Fill(Peta,weight);
    }    
    if(isoP -0.002544*Ppt < isoPL){    
      EffNVTXbp->Fill(nvtx,weight);
      EffPTbp->Fill(Ppt,weight);
      EffETAbp->Fill(Peta,weight);
     }    
    if(isoC   < isoCL){
      EffNVTXbc->Fill(nvtx,weight);
      EffPTbc->Fill(Ppt,weight);
      EffETAbc->Fill(Peta,weight);
 
    }
    if(isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNL){
      EffNVTXbn->Fill(nvtx,weight);
      EffPTbn->Fill(Ppt,weight);
      EffETAbn->Fill(Peta,weight);
      
    }

    //Medium
    if(sieie  < sieM){
       EffNVTXMbs->Fill(nvtx,weight);
       EffPTMbs->Fill(Ppt,weight);
       EffETAMbs->Fill(Peta,weight);
    }
    if(toe    < toeM){
      EffNVTXMbt->Fill(nvtx,weight);
      EffPTMbt->Fill(Ppt,weight);
      EffETAMbt->Fill(Peta,weight);
    }    
    if(isoP -0.002544*Ppt < isoPM){    
      EffNVTXMbp->Fill(nvtx,weight);
      EffPTMbp->Fill(Ppt,weight);
      EffETAMbp->Fill(Peta,weight);
     }    
    if(isoC   < isoCM){
      EffNVTXMbc->Fill(nvtx,weight);
      EffPTMbc->Fill(Ppt,weight);
      EffETAMbc->Fill(Peta,weight);
 
    }
    if(isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNM){
      EffNVTXMbn->Fill(nvtx,weight);
      EffPTMbn->Fill(Ppt,weight);
      EffETAMbn->Fill(Peta,weight);
      
    }


    //Tight
    if(sieie  < sieT){
       EffNVTXTbs->Fill(nvtx,weight);
       EffPTTbs->Fill(Ppt,weight);
       EffETATbs->Fill(Peta,weight);
    }
    if(toe    < toeT){
      EffNVTXTbt->Fill(nvtx,weight);
      EffPTTbt->Fill(Ppt,weight);
      EffETATbt->Fill(Peta,weight);
    }    
    if(isoP -0.002544*Ppt < isoPT){    
      EffNVTXTbp->Fill(nvtx,weight);
      EffPTTbp->Fill(Ppt,weight);
      EffETATbp->Fill(Peta,weight);
     }    
    if(isoC   < isoCT){
      EffNVTXTbc->Fill(nvtx,weight);
      EffPTTbc->Fill(Ppt,weight);
      EffETATbc->Fill(Peta,weight);
 
    }
    if(isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) < isoNT){
      EffNVTXTbn->Fill(nvtx,weight);
      EffPTTbn->Fill(Ppt,weight);
      EffETATbn->Fill(Peta,weight);
      
    }

  }



  TFile *f1 = new TFile("Eff1.root","recreate");
  f1->cd();

  EffPT0->Write();
  EffNVTX0->Write();
  EffETA0->Write();
  EffPHI0->Write();
    
  EffPTL->Write();
  EffNVTXL->Write();
  EffETAL->Write();
  EffPHIL->Write();
  
  EffPTM->Write();
  EffNVTXM->Write();
  EffETAM->Write();
  EffPHIM->Write();
  
  EffPTT->Write();
  EffNVTXT->Write();
  EffETAT->Write();
  EffPHIT->Write();


  EffPT0b->Write();
  EffNVTX0b->Write();
  EffETA0b->Write();
  EffPHI0b->Write();

  EffPTLb->Write();
  EffNVTXLb->Write();
  EffETALb->Write();
  EffPHILb->Write();

  EffPTMb->Write();
  EffNVTXMb->Write();
  EffETAMb->Write();
  EffPHIMb->Write();

  EffPTTb->Write();
  EffNVTXTb->Write();
  EffETATb->Write();
  EffPHITb->Write(); 



  // Loose

  EffPTs->Write();
  EffNVTXs->Write();
  EffETAs->Write();

  EffPTt->Write();
  EffNVTXt->Write();
  EffETAt->Write();

  EffPTc->Write();
  EffNVTXc->Write();
  EffETAc->Write();

  EffPTp->Write();
  EffNVTXp->Write();
  EffETAp->Write();

  EffPTn->Write();
  EffNVTXn->Write();
  EffETAn->Write();

    
  EffPTbs->Write();
  EffNVTXbs->Write();
  EffETAbs->Write();
  
  EffPTbt->Write();
  EffNVTXbt->Write();
  EffETAbt->Write();
  
  EffPTbc->Write();
  EffNVTXbc->Write();
  EffETAbc->Write();

  EffPTbp->Write();
  EffNVTXbp->Write();
  EffETAbp->Write();

  EffPTbn->Write();
  EffNVTXbn->Write();
  EffETAbn->Write();


  // Medium

  EffPTMs->Write();
  EffNVTXMs->Write();
  EffETAMs->Write();

  EffPTMt->Write();
  EffNVTXMt->Write();
  EffETAMt->Write();

  EffPTMc->Write();
  EffNVTXMc->Write();
  EffETAMc->Write();

  EffPTMp->Write();
  EffNVTXMp->Write();
  EffETAMp->Write();

  EffPTMn->Write();
  EffNVTXMn->Write();
  EffETAMn->Write();


  EffPTMbs->Write();
  EffNVTXMbs->Write();
  EffETAMbs->Write();

  EffPTMbt->Write();
  EffNVTXMbt->Write();
  EffETAMbt->Write();

  EffPTMbc->Write();
  EffNVTXMbc->Write();
  EffETAMbc->Write();

  EffPTMbp->Write();
  EffNVTXMbp->Write();
  EffETAMbp->Write();

  EffPTMbn->Write();
  EffNVTXMbn->Write();
  EffETAMbn->Write();


  // Tight

  EffPTTs->Write();
  EffNVTXTs->Write();
  EffETATs->Write();

  EffPTTt->Write();
  EffNVTXTt->Write();
  EffETATt->Write();

  EffPTTc->Write();
  EffNVTXTc->Write();
  EffETATc->Write();

  EffPTTp->Write();
  EffNVTXTp->Write();
  EffETATp->Write();

  EffPTTn->Write();
  EffNVTXTn->Write();
  EffETATn->Write();


  EffPTTbs->Write();
  EffNVTXTbs->Write();
  EffETATbs->Write();

  EffPTTbt->Write();
  EffNVTXTbt->Write();
  EffETATbt->Write();

  EffPTTbc->Write();
  EffNVTXTbc->Write();
  EffETATbc->Write();

  EffPTTbp->Write();
  EffNVTXTbp->Write();
  EffETATbp->Write();

  EffPTTbn->Write();
  EffNVTXTbn->Write();
  EffETATbn->Write();



  TFile *feta = new TFile("Eff1etaB.root","recreate");
  feta->cd();

  EffETA0->Write();
  EffETA0b->Write();

  EffETAs->Write();
  EffETAt->Write();
  EffETAc->Write();
  EffETAp->Write();
  EffETAn->Write();

  EffETAbs->Write();
  EffETAbt->Write();
  EffETAbc->Write();
  EffETAbp->Write();
  EffETAbn->Write();


}

