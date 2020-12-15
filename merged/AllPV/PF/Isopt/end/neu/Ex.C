#include <TMath.h>
#include <cstdlib>
#include <TRandom.h>
#include <TGraph.h>
#include <vector>
#include <new>
#include "sstream"
#include <string>
#include <fstream>
#include <iostream>



void ErrCalc(TH1F*, int, double, double&, double&, double&);

void Ex(){


  TString f1_name = "../../CutTMVAEndCap90_HPT.root";
  TFile *f1 = new TFile(f1_name);
  if (!f1 || !f1->IsOpen()) {
    cout << "\nERROR! Could not open root file " 
	 << f1_name << endl;
    exit(0);
  }

  float genPt,ppt,peta,Sie_ie,iso_P,iso_C,iso_N,to_e,weighT;
  int nvtx;
  gStyle->SetOptStat(0);

  //Signal Tree
  TTree *t_S = (TTree*)f1->Get("t_S");  
      
  t_S->SetBranchAddress("Sieie",&Sie_ie);
  t_S->SetBranchAddress("isoP",&iso_P);
  t_S->SetBranchAddress("isoC",&iso_C);
  t_S->SetBranchAddress("isoN",&iso_N);
  t_S->SetBranchAddress("ToE",&to_e);
  t_S->SetBranchAddress("weighT",&weighT);

  t_S->SetBranchAddress("Nvtx",&nvtx);
  t_S->SetBranchAddress("Peta",&peta);
  t_S->SetBranchAddress("Ppt",&ppt);
  t_S->SetBranchAddress("genPt",&genPt);

  //Background Tree
  TTree *t_B = (TTree*)f1->Get("t_B");
     
  t_B->SetBranchAddress("Sieie",&Sie_ie);
  t_B->SetBranchAddress("isoP",&iso_P);
  t_B->SetBranchAddress("isoC",&iso_C);
  t_B->SetBranchAddress("isoN",&iso_N);
  t_B->SetBranchAddress("ToE",&to_e);
  t_B->SetBranchAddress("weighT",&weighT);

  t_B->SetBranchAddress("Nvtx",&nvtx);
  t_B->SetBranchAddress("Peta",&peta);
  t_B->SetBranchAddress("Ppt",&ppt);

  TH2F *isoPptS = new TH2F("isoPptS","Iso Neutral vs Pt",250,0,1000,1000,0,100);


  for(int i  = 0; i < t_S->GetEntries();i++){
    t_S->GetEntry(i);

    if(ppt < 15 || ppt > 1000) continue;
    if(iso_P == 0) continue;

    isoPptS->Fill(ppt,iso_N);
  }

  cout<<"Builded the 2d HISTOGRAM"<<endl;



  TH2F *his2 = (TH2F*)isoPptS->Clone();

  int dim = his2->GetXaxis()->GetNbins(); 
  
  double * cutV; 
  double * errVH; 
  double * errVL; 
  double * binc; 
  double * bincerH; 
  double * bincerL; 
  cutV   = new double[dim];
  errVH   = new double[dim];
  errVL   = new double[dim];
  binc   = new double[dim]; 
  bincerH = new double[dim];
  bincerL = new double[dim];
  
  for(int i  = 1; i < dim ; i++){
    double xval = 0; 
    double errXH = 0; 
    double errXL = 0; 
    
    TH1F *r22 = (TH1F*)his2->ProjectionY(" ",i,i+1," ");
    
    TH1F *h1 = (TH1F*)r22->Clone();
    ErrCalc(h1,i,0.950,xval,errXL,errXH);
    
    cout << "bin: " << i << " "
	 << "x-val: " << xval << " " 
	 << "-" << errXL << " " 
	 << "+" << errXH << endl;
    
    cutV[i-1]   = xval; 
    
    errVL[i-1]   = errXL; 
    errVH[i-1]   = errXH;
    
    binc[i-1]   = his2->GetXaxis()->GetBinCenter(i);
    bincerL[i-1] = 0;
    bincerH[i-1] = 0;
  }
  TGraphAsymmErrors * IsoptScaling = new TGraphAsymmErrors(dim,binc,cutV,
							   bincerL,bincerH,
							   errVL,errVH);

  TGraphAsymmErrors * IsoptScaling2 = new TGraphAsymmErrors(dim,binc,cutV,
							    bincerL,bincerH,
							    errVL,errVH);

  TGraphAsymmErrors * IsoptScalingLin = new TGraphAsymmErrors(dim,binc,cutV,
							      bincerL,bincerH,
							      errVL,errVH);

  // double DownL = his2->GetXaxis()->GetBinCenter(1);
  // double UpperL = his2->GetXaxis()->GetBinCenter(dim);
  
  TF1 *fn1 = new TF1("fn1","exp([0]*x + [1])",20,400);
  TF1 *fn2 = new TF1("fn2","[1]*x + [2]*x*x + [0]",20,400);
  TF1 *fnlin = new TF1("fnlin","[1]*x + [0]",20,400);

  IsoptScaling->Fit("fn1","R");
  IsoptScaling2->Fit("fn2","R");
  IsoptScalingLin->Fit("fnlin","R");
  

  gStyle->SetOptFit(1);
  
  TCanvas *c3 = new TCanvas("c3","Iso Pt",1200,600);
  c3->Divide(2,1);
  
  c3->cd(1);
  IsoptScaling->SetMarkerStyle(24); 
  IsoptScaling->SetMarkerSize(0.4);
  IsoptScaling->Draw("AP");
  IsoptScaling->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScaling->GetYaxis()->SetTitle("Isolation Contour PF::h0 rho corr.");
  c3->cd(2); 
  c3->cd(2);
  c3->SetRightMargin(0.3);
  c3->SetLogz();
  his2->Draw("colz");
  his2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  his2->GetYaxis()->SetTitle("Isolation PF::h0  rho corr.");

  c3->SaveAs("./IsoPt_end_neu/final_400.png");
  c3->SaveAs("./IsoPt_end_neu/final_400.C");
  

  
  TCanvas *c4 = new TCanvas("c4","Iso Pt",1200,600);
  c4->Divide(2,1);
  
  c4->cd(1);
  IsoptScaling2->SetMarkerStyle(24); 
  IsoptScaling2->SetMarkerSize(0.4);
  IsoptScaling2->Draw("AP");
  IsoptScaling2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScaling2->GetYaxis()->SetTitle("Isolation Contour PF::h0 rho corr.");

  c4->cd(2);
  IsoptScalingLin->SetMarkerStyle(24); 
  IsoptScalingLin->SetMarkerSize(0.4);
  IsoptScalingLin->Draw("AP");
  IsoptScalingLin->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScalingLin->GetYaxis()->SetTitle("Isolation Contour PF::h0 rho corr.");

  c4->SaveAs("./IsoPt_end_neu/POl.png");
  c4->SaveAs("./IsoPt_end_neu/POl.C");



  TCanvas *c6 = new TCanvas("c6","Iso Pt",1200,600);
  c6->Divide(2,1);  
  c6->cd(1);
  IsoptScaling2->SetMarkerStyle(24); 
  IsoptScaling2->SetMarkerSize(0.4);
  IsoptScaling2->Draw("AP");
  IsoptScaling2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  IsoptScaling2->GetYaxis()->SetTitle("Isolation Contour PF::h0 rho corr.");
  c6->cd(2);
  c6->SetRightMargin(0.3);
  c6->SetLogz();
  his2->Draw("colz"); 
  his2->GetXaxis()->SetTitle("Photon Pt GeVc^{-1}");
  his2->GetYaxis()->SetTitle("Isolation PF::h0  rho corr.");

  c6->SaveAs("./IsoPt_end_neu/POL2.png");
  c6->SaveAs("./IsoPt_end_neu/POL2.C");


}



void ErrCalc(TH1F *HIST,int binxn,double perc,double & X_val, double & errXL,double & errXH){
  
  // double perc = 0.95;  // efficiency contour value
  //int binxn = 0;       // name of the bin to be scanned

  // Naming the output 
  ostringstream gn; 
  gn << binxn; 
  string gnm_png = "./IsoPt_end_neu_eff/Efficiency_CutVal"+gn.str()+".png";
  string gnm_C   = "./IsoPt_end_neu_eff/Efficiency_CutVal"+gn.str()+".C";
  const char *graphname_png = gnm_png.c_str();
  const char *graphname_C   = gnm_C.c_str();

  /*   TEST CODE
  TH1F *Gaus1 = new TH1F("Gaus1","Gaussian for test",100,0,100);
  
  TRandom r1; 
  
  for(int i  = 0 ; i < 1000; i++){
    double x = r1.Gaus(30,4);
    Gaus1->Fill(x);
  }
  */ 

  TH1F *h1 = (TH1F*)HIST->Clone();
  int arsize = h1->GetXaxis()->GetNbins(); 
  
  if(h1->GetEntries() != 0 ){// goto endd; 

  double *eff; 
  double *eff_err; 
  double *cutV; 
  double *cutV_err; 
  
  eff      = new double[arsize];
  cutV     = new double[arsize];
  eff_err  = new double[arsize];
  cutV_err = new double[arsize];
  h1->Draw();
  int tot = h1->GetEntries();
  float tot_nz = tot - h1->GetBinContent(1);
  // cout<<tot<<endl;
    int integ = 0; 
    for(int i  = 1; i < (h1->GetXaxis()->GetNbins() + 1); i++){ 
      double xCut = h1->GetXaxis()->GetBinLowEdge(i); 
      integ += h1->GetBinContent(i); 
     

      if(integ != 0 && tot != 0 ){ 
	eff[i -1] = (integ*1.0/tot);
	//eff_err[i -1] = ((integ*1.0/tot)*sqrt(pow(sqrt(tot)/tot,2) + pow(sqrt(integ)/integ,2) )); 
	eff_err[i -1] = sqrt ( eff[i-1] * ( 1 - eff[i-1] ) / tot_nz ) ;
      }else{
	eff_err[i -1] = 0; 
	eff[i -1] = 0;
      }
      
      cutV[i -1] = (xCut);
      cutV_err[i - 1] = 0;
      //  cout<<"bin :"<<i<<"  % :"<<integ*1.0/tot<<" cut val :"<<xCut <<endl;
    }

    
    gStyle->SetOptStat(1);

    //Now draw the resulting curve 
    TCanvas *c1 = new TCanvas("c1","The Eff- cut value plot",1200,600);
    
    c1->Divide(2,1);

    c1->cd(1);
    double x[2] = {h1->GetXaxis()->GetBinCenter(1),h1->GetXaxis()->GetBinCenter(arsize)};
    double y[2] = {perc,perc};

    TGraph *lineP =  new TGraph(2,x,y); 
    lineP->SetLineColor(kRed);

    TGraphErrors *efC = new TGraphErrors(arsize ,cutV,eff,cutV_err,eff_err);
    efC->SetMarkerStyle(20);
    efC->SetMarkerSize(0.5);
    
    //   efC->Draw("AP");
    efC->GetXaxis()->SetTitle("Cut Value");
    efC->GetYaxis()->SetTitle("Efficiency");
    



    TMultiGraph *GRPhs = new TMultiGraph(); 
    GRPhs->Add(lineP,"l");
    GRPhs->Add(efC,"p");

    GRPhs->Draw("APL");

    c1->cd(2);

    h1->Draw();

    c1->SaveAs(graphname_png);
    c1->SaveAs(graphname_C);

    double * err_up; 
    double * err_down; 

    err_up = new double [arsize];
    err_down = new double [arsize];
    
    for(int i = 0;  i < arsize ; i++){
      err_up[i] = eff[i] + eff_err[i];
      err_down[i] = eff[i] - eff_err[i];
      // cout<< "i: "<<i <<" "<<eff[i] <<" "<<err_up[i]<<"  "<<err_down[i]<<endl;
    }
   


    //Extrapolation method to find the CutValue errors
    int up = -99; 
    int down = -99;  
    for( int i  = 0; i < arsize ; i++){
      if( err_up [i] > perc){ up = i; 
      	break;
      }
    }
    for( int i  = 0; i < arsize ; i++){
      if( err_down [i] > perc){
	down = i; 
	break;
      }
    }
   
    // here it returns the extrapolated result

    double Usl   = (err_up[up] - err_up[up-1])/(cutV[up] - cutV[up-1]); 
    double Ustom = err_up[up] - Usl*cutV[up]; 
    double err_1 = (perc - Ustom)/Usl; 


    double Dsl   = (err_down[down] - err_down[down-1])/(cutV[down] - cutV[down-1]); 
    double Dstom = err_down[down] - Dsl*cutV[down]; 
    double err_2 = (perc - Dstom)/Dsl; 
    
    
    // cout<<"down err:"<<err_1<<endl;
    // cout<<"UP err:"<<err_2<<endl;
    //cout<<"cut value"<<(err_1+err_2)/2.; 
    //cout<<"sym err:"<< (err_2-err_1)/2.;

    int xc = 0; 
   
    for(int i = 0; i < arsize; i++){
      if( eff[i] > perc ){ 
	xc = i; 
	break;
      }
    }
    
   

    if(fabs(eff[xc] - perc) > fabs(eff[xc - 1] - perc) ){
      X_val = cutV[xc -1 ]; 
   
    }else{
      X_val = cutV[xc ];
   
    }
    
    // X_val = (err_1+err_2)/2.; 

    //checking case that the lines do not cross ! 

    if(up == -99){
      errXL = X_val - h1->GetXaxis()->GetBinCenter(1);
    }else{
      errXL = X_val - err_1;
    }

    if(down == -99){
      
      // errXL = ( h1->GetXaxis()->GetBinCenter(arsize) - X_val );
      int dif = 0; 
      for(int i = 0; i < arsize; i++){
	if(eff[i] >  0.98){
	  dif = i; 
	  break; 
	}
      }
      errXH = (cutV[dif] - X_val );


    }else{
       errXH = err_2 - X_val;
    }

    
    


    //Deleting all the dynamical arrays

    
    delete[] eff;
    delete[] cutV;

    delete[] eff_err;
    delete[] cutV_err;

    delete[] err_up; 
    delete[] err_down; 
  }  
 
  //endd:
  if(h1->GetEntries() == 0)  cout<<"Empty histo"<<endl;

}
