#include "sstream"
#include <string>
#include <fstream>
#include <iostream>
#include <ErrScanning.C>



void FitterAL(double e_min,double e_max,int bin,double &isoP,double &isoC,double &isoN){
  //------------Describtion of file functions---------------
  // This Macro calculates and returns the effective areas 
  // and the the profile plots  
  //
  // Output of this file : 1)  A canvas saved as Bin_X.png 
  // with the Isolation profiles and the Rho fitted for sanity checks 
  // where X is the specific eta bin 
  // 
  // 2) A file called Bin_X_ProfPlots.root containing the fitted profile 
  // plots of the Rho and three isolations.
  //
  // 3)It returns (as a function) the three effective areas under the names :
  // isoP, isoC, isoN, by changing the variables to these values that are 
  // fed to it.  
  //
  //
  // Required Input to this Macro: 3 variables that must be checked to 
  // correspond to the correct isolation, to be changed and given the 
  // effective areas values. 
  // The Bin number of eta, because with that we open the file 
  // containing the 2d plots to create the profles. And thus these files 
  // must exist. These files are generated by the AreaCalc.C macro.
  //
  //
  //-------------------------------------------------------------------------
  
  



  ostringstream biin; 
  biin << bin; 
  string fop = "AreaCalc_"+biin.str()+".root";
  string fou = "BinL_"+biin.str()+"ProfPlots.root";
  string Pngg = "BinL_"+biin.str()+".png";
  string ftxt  = "EfAL"+biin.str()+".txt";


  const char *fnmm = fop.c_str();
  const char *outP = fou.c_str();
  const char *PNG = Pngg.c_str();
  const char *fttx  = ftxt.c_str(); 
  
  ofstream offf; 
  offf.open(fttx);


  TFile *f1 = new TFile(fnmm);
  
  f1->cd();
  cout<<"openinng file"<<endl;
                                                                                                      
  TProfile *IsoPvsRhop = IsoPvsRho->ProfileX();
  TProfile *IsoCvsRhop = IsoCvsRho->ProfileX();
  TProfile *IsoNvsRhop = IsoNvsRho->ProfileX();

  //Calling the Contour and Err calculator

  ErrScan(IsoPvsRho);
  ErrScan(IsoNvsRho);
  ErrScan(IsoNvsRho);




  

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);



  TF1 *fn2 = new TF1("fn2","[0]*x + [1]",5,30);
  TF1 *fn3 = new TF1("fn3","[0]*x + [1]",5,30);
  TF1 *fn4 = new TF1("fn4","[0]*x + [1]",5,30);

  //Fitting the Profile plots


  IsoPvsRhop->Fit("fn2","R");
  IsoNvsRhop->Fit("fn3","R");
  IsoCvsRhop->Fit("fn4","R");


  
  TCanvas *c1  = new TCanvas("c1","Iso and Rhos",1200,400);
  c1->Divide(3,1);
  c1->cd(1);
  IsoPvsRhop->Draw();                                                           
  IsoPvsRhop->GetXaxis()->SetTitle(" Rho ");
  IsoPvsRhop->GetYaxis()->SetTitle(" Iso Profile ");


  c1->cd(2);
  IsoNvsRhop->Draw();    
  IsoNvsRhop->GetXaxis()->SetTitle(" Rho ");
  IsoNvsRhop->GetYaxis()->SetTitle(" Iso Profile ");
  
  c1->cd(3);
  IsoCvsRhop->Draw();    
  IsoCvsRhop->GetXaxis()->SetTitle(" Rho ");
  IsoCvsRhop->GetYaxis()->SetTitle(" Iso Profile ");
  c1->SaveAs(PNG);

  //Retrieving the fit parameters to derive effective areas  

  TF1 *fgm  = IsoPvsRhop->GetFunction("fn2");
  TF1 *fne  = IsoNvsRhop->GetFunction("fn3");
  TF1 *fch1 = IsoCvsRhop->GetFunction("fn4");

  //  double Rho_a  = fR->GetParameter(0);
  double Isog_a = fgm->GetParameter(0);
  double Isog_aer = fgm->GetParError(0);

  double Ison_a = fne->GetParameter(0);
  double Ison_aer = fne->GetParError(0);

  double Isoc_a = fch1->GetParameter(0);
  double Isoc_aer = fch1->GetParError(0);
  
  cout<<"-------------effective areas -----"<<endl;
  cout<<" gammas  iso EA :"<< Isog_a<<" "<<Isog_aer <<endl;
  cout<<" neutral iso EA :"<< Ison_a<<" "<<Ison_aer <<endl;
  cout<<" charge  iso EA :"<< Isoc_a<<" "<<Isoc_aer <<endl;
  

  offf<<"-------------effective areas -----"<<endl;
  offf<<" gammas  iso EA :"<< Isog_a<<" "<<Isog_aer <<endl;
  offf<<" neutral iso EA :"<< Ison_a<<" "<<Ison_aer <<endl;
  offf<<" charge  iso EA :"<< Isoc_a<<" "<<Isoc_aer <<endl;

  offf.close();


  isoP = (Isog_a >= 0 ) ? Isog_a : 0;
  isoN = (Ison_a >= 0 ) ? Ison_a : 0;
  isoC = (Isoc_a >= 0 ) ? Isoc_a : 0;

  
  TFile *f2 = new TFile(outP,"RECREATE");
  f2->cd();


  IsoPvsRhop->Write();
  IsoCvsRhop->Write();
  IsoNvsRhop->Write();

  f2->Close();
}
