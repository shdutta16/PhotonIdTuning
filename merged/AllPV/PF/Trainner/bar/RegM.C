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
#include "TMVA/DataLoader.h"
//#include "TMVARegGui.C"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"

using namespace TMVA;


void RegM(){
  
  TMVA::Tools::Instance();
  std::cout << "==> Start TMVARegression" << std::endl;
    
  ifstream myfile; TString myfile_name;
  myfile_name = "LooseR.txt";
  myfile.open(myfile_name);

  if (!myfile || !myfile.is_open()) {
    cout << "\nERROR! Could not open txt file " << myfile_name
         << endl;
    exit(0);
  }


  ostringstream xcS,xcH,xcP,xcC,xcN;  
  double xS,xC,xN,xP,xH;

  if(myfile.is_open()){
    while(!myfile.eof()){
      myfile>>xS>>xC>>xN>>xP>>xH;
    }
  }

  xcS<<xS;
  xcC<<xC;
  xcN<<xN;
  xcP<<xP;
  xcH<<xH;

  //Output file 
  TString outfileName( "Medium_w_BAR_scaled.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  //Declaring the factory
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", 
					      outputFile, 
					      "!V:!Silent:Color:DrawProgressBar" );

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");


  //Declaring Input Variables 
                                                                                  
  dataloader->AddVariable( "Sieie",'F');
  dataloader->AddVariable( "isoC",'F' );
  dataloader->AddVariable( "(isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) > 0 ) ? isoN-(0.01556*Ppt-0.000001129*Ppt*Ppt) : 0.0 ",'F' );
  dataloader->AddVariable( "(isoP-0.002544*Ppt > 0 ) ? isoP-0.002544*Ppt : 0.0 ",'F' );
  dataloader->AddVariable( "(ToE > 0.0) ? ToE : 0.0 ",'F' );


  //factory->AddVariable("(isoN-(0.014*Ppt+0.000019*Ppt*Ppt) > 0 ) ? isoN-(0.014*Ppt+0.000019*Ppt*Ppt) : 0.0 ",'F' );
  //factory->AddVariable( "(isoP-0.0053*Ppt > 0 ) ? isoP-0.0053*Ppt : 0.0 ",'F' );

  dataloader->AddSpectator( "Ppt",'F' );
  //dataloader->AddSpectator( "ToE",'F' );

  
  TString fname = "../CutTMVABarrel90_test.root";
  TFile *input = TFile::Open( fname );
  if (!input || !input->IsOpen()) {
    cout << "\nERROR! Could not open root file " << fname
         << endl;
    exit(0);
  }

  
  // --- Register the regression tree
  TTree *signal = (TTree*)input->Get("t_S");
  TTree *background = (TTree*)input->Get("t_B");
  
  //Just Some more settings
   Double_t signalWeight      = 1.0; 
   Double_t backgroundWeight  = 1.0; 

   // You can add an arbitrary number of regression trees
   dataloader->AddSignalTree( signal, signalWeight );
   dataloader->AddBackgroundTree( background , backgroundWeight );
 
   //TCut mycuts ="Ppt>15 && Ppt < 200 && ToE < 0.05";
   //TCut mycutb ="Ppt>15 && Ppt < 200 && ToE < 0.05";
   TCut mycuts ="Ppt>15  && Ppt < 200";
   TCut mycutb ="Ppt>15  && Ppt < 200";
   //factory->PrepareTrainingAndTestTree(mycuts,mycutb,"");
   dataloader->PrepareTrainingAndTestTree(mycuts,mycutb,"nTrain_Signal=1000000:nTrain_Background=1000000:nTest_Signal=200000:nTest_Background=200000");
   //factory->PrepareTrainingAndTestTree(mycuts,mycutb,"nTrain_Signal=300000:nTrain_Background=300000:nTest_Signal=300000:nTest_Background=300000");
   //factory->PrepareTrainingAndTestTree(mycuts,mycutb,"nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=1000:nTest_Background=1000");
   //"nTrain_Signal=9000:nTrain_Background=9000:nTest_Signal=10000:nTest_Background=10000");x

   dataloader->SetBackgroundWeightExpression("weighT");
   dataloader->SetSignalWeightExpression("weighT");

   TString methodName = "Cut_Medium_r";
   TString methodOptions ="!H:!V:FitMethod=GA:EffMethod=EffSEl"; 
   methodOptions +=":VarProp[0]=FMin:VarProp[1]=FMin:VarProp[2]=FMin:VarProp[3]=FMin:VarProp[4]=FMin";
  
   methodOptions +=":CutRangeMax[0]="+xcS.str(); 
   methodOptions +=":CutRangeMax[1]="+xcC.str();
   methodOptions +=":CutRangeMax[2]="+xcN.str();
   methodOptions +=":CutRangeMax[3]="+xcP.str();
   methodOptions +=":CutRangeMax[4]="+xcH.str();
   //   methodOptions +=":CutRangeMax[4]="+xcP.str();

   //************
   factory->BookMethod(dataloader,
		       TMVA::Types::kCuts,
		       methodName,
		       methodOptions);

   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();    
   
   // --------------------------------------------------------------
   // Save the output
   outputFile->Close();



   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      
   delete factory;

}
