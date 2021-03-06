void Iso_vsRho5()
{
//=========Macro generated from canvas: c3/Iso vs Pt
//=========  (Thu Nov 19 20:55:52 2020) by ROOT version 6.22/02
   TCanvas *c3 = new TCanvas("c3", "Iso vs Pt",0,0,1200,400);
   gStyle->SetOptFit(1);
   c3->Range(0,0,1,1);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c3_1
   TPad *c3_1 = new TPad("c3_1", "c3_1",0.01,0.01,0.3233333,0.99);
   c3_1->Draw();
   c3_1->cd();
   c3_1->Range(-13.675,-0.4675,123.075,4.2075);
   c3_1->SetFillColor(0);
   c3_1->SetBorderMode(0);
   c3_1->SetBorderSize(2);
   c3_1->SetFrameBorderMode(0);
   c3_1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx3037[100] = {
   0.5,
   1.5,
   2.5,
   3.5,
   4.5,
   5.5,
   6.5,
   7.5,
   8.5,
   9.5,
   10.5,
   11.5,
   12.5,
   13.5,
   14.5,
   15.5,
   16.5,
   17.5,
   18.5,
   19.5,
   20.5,
   21.5,
   22.5,
   23.5,
   24.5,
   25.5,
   26.5,
   27.5,
   28.5,
   29.5,
   30.5,
   31.5,
   32.5,
   33.5,
   34.5,
   35.5,
   36.5,
   37.5,
   38.5,
   39.5,
   40.5,
   41.5,
   42.5,
   43.5,
   44.5,
   45.5,
   46.5,
   47.5,
   48.5,
   49.5,
   50.5,
   51.5,
   52.5,
   53.5,
   54.5,
   55.5,
   56.5,
   57.5,
   58.5,
   59.5,
   60.5,
   61.5,
   62.5,
   63.5,
   64.5,
   65.5,
   66.5,
   67.5,
   68.5,
   69.5,
   70.5,
   71.5,
   72.5,
   73.5,
   74.5,
   75.5,
   76.5,
   77.5,
   78.5,
   79.5,
   80.5,
   81.5,
   82.5,
   83.5,
   84.5,
   85.5,
   86.5,
   87.5,
   88.5,
   89.5,
   90.5,
   91.5,
   92.5,
   93.5,
   94.5,
   95.5,
   96.5,
   97.5,
   98.5,
   99.5};
   Double_t Graph0_fy3037[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   3.3,
   3.4,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.3,
   3.2,
   3.2,
   3.2,
   3.1,
   3.2,
   3.2,
   3,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_felx3037[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fely3037[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0.09533871,
   0.1591513,
   0.05885087,
   0.1061852,
   0.09443806,
   0.07969357,
   0.06790836,
   0.06129924,
   0.07222418,
   0.08103475,
   0.05351062,
   0.05447748,
   0.04176058,
   0.04821061,
   0.07937684,
   0.08586103,
   0.07884152,
   0.076702,
   0.07949686,
   0.06678794,
   0.08215473,
   0.09235877,
   0.09877897,
   0.08146918,
   0.08869003,
   0.1285111,
   0.134034,
   0.1345003,
   0.1183767,
   0.1287844,
   0.1544047,
   0.1510313,
   0.09553282,
   0.09689239,
   0.1219722,
   0.1392407,
   0.2072679,
   0.1994103,
   0.1714731,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fehx3037[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fehy3037[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0.048604,
   -0.03836531,
   0.05085856,
   -0.008550587,
   -0.009146206,
   0.004321207,
   0.01239114,
   0.00779862,
   -0.01128678,
   -0.02035159,
   0.009464091,
   0.002577673,
   0.01370095,
   0.01162782,
   -0.01403825,
   -0.02435012,
   -0.0217394,
   -0.01606167,
   -0.01656452,
   -0.004300485,
   -0.01223082,
   -0.02365216,
   -0.0305853,
   -0.009358903,
   -0.01202622,
   -0.04225911,
   -0.04501596,
   -0.04047776,
   -0.005693968,
   -0.01032374,
   -0.02519331,
   -0.01476594,
   0.05056153,
   0.06433347,
   0.05234076,
   0.06254018,
   0.03132605,
   0.04485698,
   0.1175542,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(100,Graph0_fx3037,Graph0_fy3037,Graph0_felx3037,Graph0_fehx3037,Graph0_fely3037,Graph0_fehy3037);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(24);
   grae->SetMarkerSize(0.4);
   
   TH1F *Graph_Graph03037 = new TH1F("Graph_Graph03037","Graph",100,0,109.4);
   Graph_Graph03037->SetMinimum(0);
   Graph_Graph03037->SetMaximum(3.74);
   Graph_Graph03037->SetDirectory(0);
   Graph_Graph03037->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph03037->SetLineColor(ci);
   Graph_Graph03037->GetXaxis()->SetTitle("#rho");
   Graph_Graph03037->GetXaxis()->SetLabelFont(42);
   Graph_Graph03037->GetXaxis()->SetTitleOffset(1);
   Graph_Graph03037->GetXaxis()->SetTitleFont(42);
   Graph_Graph03037->GetYaxis()->SetTitle("Neutral Isolation 95% Contour ");
   Graph_Graph03037->GetYaxis()->SetLabelFont(42);
   Graph_Graph03037->GetYaxis()->SetTitleFont(42);
   Graph_Graph03037->GetZaxis()->SetLabelFont(42);
   Graph_Graph03037->GetZaxis()->SetTitleOffset(1);
   Graph_Graph03037->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03037);
   
   
   TF1 *fnn3038 = new TF1("fnn","[0]*x + [1]",20,70, TF1::EAddToList::kNo);
   fnn3038->SetFillColor(19);
   fnn3038->SetFillStyle(0);
   fnn3038->SetLineColor(2);
   fnn3038->SetLineWidth(2);
   fnn3038->SetChisquare(12.79577);
   fnn3038->SetNDF(37);
   fnn3038->GetXaxis()->SetLabelFont(42);
   fnn3038->GetXaxis()->SetTitleOffset(1);
   fnn3038->GetXaxis()->SetTitleFont(42);
   fnn3038->GetYaxis()->SetLabelFont(42);
   fnn3038->GetYaxis()->SetTitleFont(42);
   fnn3038->SetParameter(0,-0.002945898);
   fnn3038->SetParError(0,0.0005726374);
   fnn3038->SetParLimits(0,0,0);
   fnn3038->SetParameter(1,3.385068);
   fnn3038->SetParError(1,0.01914807);
   fnn3038->SetParLimits(1,0,0);
   fnn3038->SetParent(grae);
   grae->GetListOfFunctions()->Add(fnn3038);
   
   TPaveStats *ptstats = new TPaveStats(0.59,0.04000001,0.95,0.2,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf =  12.8 / 37");
   ptstats_LaTex = ptstats->AddText("p0       = -0.002946 #pm 0.0005726 ");
   ptstats_LaTex = ptstats->AddText("p1       = 3.385 #pm 0.01915 ");
   ptstats->SetOptStat(0);
   ptstats->SetOptFit(111);
   ptstats->Draw();
   grae->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(grae->GetListOfFunctions());
   grae->Draw("ap");
   
   TPaveText *pt = new TPaveText(0.4166171,0.9291672,0.5833829,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Graph");
   pt->Draw();
   c3_1->Modified();
   c3->cd();
  
// ------------>Primitives in pad: c3_2
   TPad *c3_2 = new TPad("c3_2", "c3_2",0.3433333,0.01,0.6566667,0.99);
   c3_2->Draw();
   c3_2->cd();
   c3_2->Range(-13.675,-0.9644555,123.075,8.680099);
   c3_2->SetFillColor(0);
   c3_2->SetBorderMode(0);
   c3_2->SetBorderSize(2);
   c3_2->SetFrameBorderMode(0);
   c3_2->SetFrameBorderMode(0);
   
   Double_t Graph0_fx3040[100] = {
   0.5,
   1.5,
   2.5,
   3.5,
   4.5,
   5.5,
   6.5,
   7.5,
   8.5,
   9.5,
   10.5,
   11.5,
   12.5,
   13.5,
   14.5,
   15.5,
   16.5,
   17.5,
   18.5,
   19.5,
   20.5,
   21.5,
   22.5,
   23.5,
   24.5,
   25.5,
   26.5,
   27.5,
   28.5,
   29.5,
   30.5,
   31.5,
   32.5,
   33.5,
   34.5,
   35.5,
   36.5,
   37.5,
   38.5,
   39.5,
   40.5,
   41.5,
   42.5,
   43.5,
   44.5,
   45.5,
   46.5,
   47.5,
   48.5,
   49.5,
   50.5,
   51.5,
   52.5,
   53.5,
   54.5,
   55.5,
   56.5,
   57.5,
   58.5,
   59.5,
   60.5,
   61.5,
   62.5,
   63.5,
   64.5,
   65.5,
   66.5,
   67.5,
   68.5,
   69.5,
   70.5,
   71.5,
   72.5,
   73.5,
   74.5,
   75.5,
   76.5,
   77.5,
   78.5,
   79.5,
   80.5,
   81.5,
   82.5,
   83.5,
   84.5,
   85.5,
   86.5,
   87.5,
   88.5,
   89.5,
   90.5,
   91.5,
   92.5,
   93.5,
   94.5,
   95.5,
   96.5,
   97.5,
   98.5,
   99.5};
   Double_t Graph0_fy3040[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   1.7,
   1.7,
   0.9,
   0.4,
   0.8,
   0.8,
   0.9,
   1,
   1,
   1.1,
   1.1,
   1.1,
   1.2,
   1.3,
   1.3,
   1.3,
   1.3,
   1.4,
   1.4,
   1.4,
   1.5,
   1.5,
   1.6,
   1.6,
   1.6,
   1.7,
   1.7,
   1.8,
   1.8,
   1.8,
   1.9,
   1.9,
   2,
   2,
   2,
   2.1,
   2.1,
   2.2,
   2.2,
   2.2,
   2.3,
   2.3,
   2.4,
   2.4,
   2.4,
   2.5,
   2.5,
   2.5,
   2.6,
   2.6,
   2.7,
   2.8,
   2.7,
   2.7,
   2.8,
   2.9,
   2.9,
   3,
   3,
   3.1,
   3.2,
   3.2,
   3.4,
   3.4,
   3.4,
   3.3,
   3.3,
   3.9,
   4,
   3.9,
   3.6,
   3.6,
   2.9,
   2.3,
   2.7,
   3,
   3.6,
   3.6,
   3.6,
   5.8,
   5.9,
   1.4,
   3.4,
   3.4,
   1,
   0,
   0,
   3.1,
   3.1,
   0};
   Double_t Graph0_felx3040[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fely3040[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0.06828427,
   0.830718,
   0.6770784,
   0.4,
   0.3477429,
   0.1833126,
   0.1090676,
   0.1554968,
   0.1065177,
   0.1483456,
   0.06771351,
   0.07097948,
   0.1145069,
   0.1273987,
   0.1227761,
   0.09319167,
   0.02185072,
   0.08651344,
   0.05013262,
   0.01511183,
   0.0842587,
   0.02960907,
   0.09060291,
   0.06293989,
   0.03433083,
   0.0927358,
   0.0330989,
   0.08066456,
   0.05787239,
   0.02077712,
   0.07588828,
   0.05063133,
   0.1053052,
   0.05522235,
   0.02524958,
   0.09575548,
   0.04141906,
   0.09410166,
   0.07239876,
   0.03521228,
   0.08759624,
   0.04072447,
   0.07087071,
   0.02996331,
   0.02606817,
   0.1184066,
   0.09602403,
   0.03498432,
   0.09022928,
   0.06155564,
   0.09542736,
   0.1402005,
   0.05857803,
   0.06729859,
   0.1074392,
   0.08540703,
   0.1374775,
   0.1939607,
   0.1351204,
   0.1897434,
   0.2214504,
   0.1924253,
   0.2201242,
   0.370887,
   0.570785,
   0.4866977,
   0.4271746,
   0.7398359,
   0.7396918,
   1.072251,
   1.00427,
   1.043327,
   0.6567172,
   1.035701,
   0.8029571,
   0.3361901,
   0.8687423,
   1.516004,
   2.700855,
   5.397818,
   5.625861,
   1.225861,
   2.064495,
   2.43,
   1,
   0,
   0,
   3.1,
   3.1,
   0};
   Double_t Graph0_fehx3040[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fehy3040[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   -0.01171573,
   -0.02143594,
   -0.01797959,
   0.4390431,
   0.02355202,
   0.1277971,
   0.05134184,
   -0.02362255,
   0.06139071,
   -0.01113528,
   0.0344181,
   0.005162786,
   -0.05099167,
   -0.06864505,
   -0.07106693,
   -0.04567449,
   0.01895858,
   -0.04890459,
   -0.01781062,
   0.01394854,
   -0.05617447,
   -0.005635937,
   -0.0660773,
   -0.04018153,
   -0.01352076,
   -0.07037827,
   -0.01218111,
   -0.05824139,
   -0.03677469,
   -6.769145e-05,
   -0.05297366,
   -0.02732106,
   -0.08154341,
   -0.03094301,
   -0.0004776532,
   -0.06911144,
   -0.01332549,
   -0.0655577,
   -0.04187054,
   -0.002965272,
   -0.05155549,
   -0.003421953,
   -0.03166681,
   0.01185662,
   0.02045473,
   -0.06734446,
   -0.04437977,
   0.02709672,
   -0.02508823,
   0.00623657,
   -0.007131983,
   -0.04077226,
   0.03373957,
   0.02982016,
   0.05012835,
   0.06077577,
   0.02799689,
   -0.04291355,
   0.02799827,
   0.05007816,
   0.05175556,
   0.1745764,
   0.1014626,
   0.06774783,
   0.1252445,
   0.2561007,
   0.464999,
   0.3118766,
   0.2143788,
   0.1806226,
   0.4552478,
   0.4909969,
   0.5504978,
   0.6793738,
   4.314222,
   0.825849,
   0.2939176,
   1.932922,
   2.137868,
   0.04900348,
   -0.0261895,
   -0.0261895,
   -0.025,
   -0.015,
   -1,
   0,
   0,
   -3.1,
   -3.1,
   0};
   grae = new TGraphAsymmErrors(100,Graph0_fx3040,Graph0_fy3040,Graph0_felx3040,Graph0_fehx3040,Graph0_fely3040,Graph0_fehy3040);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(24);
   grae->SetMarkerSize(0.4);
   
   TH1F *Graph_Graph03040 = new TH1F("Graph_Graph03040","Graph",100,0,109.4);
   Graph_Graph03040->SetMinimum(0);
   Graph_Graph03040->SetMaximum(7.715644);
   Graph_Graph03040->SetDirectory(0);
   Graph_Graph03040->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph03040->SetLineColor(ci);
   Graph_Graph03040->GetXaxis()->SetTitle("#rho");
   Graph_Graph03040->GetXaxis()->SetLabelFont(42);
   Graph_Graph03040->GetXaxis()->SetTitleOffset(1);
   Graph_Graph03040->GetXaxis()->SetTitleFont(42);
   Graph_Graph03040->GetYaxis()->SetTitle("Photon Isolation 95% Contour ");
   Graph_Graph03040->GetYaxis()->SetLabelFont(42);
   Graph_Graph03040->GetYaxis()->SetTitleFont(42);
   Graph_Graph03040->GetZaxis()->SetLabelFont(42);
   Graph_Graph03040->GetZaxis()->SetTitleOffset(1);
   Graph_Graph03040->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03040);
   
   
   TF1 *fnp3041 = new TF1("fnp","[0]*x + [1]",20,80, TF1::EAddToList::kNo);
   fnp3041->SetFillColor(19);
   fnp3041->SetFillStyle(0);
   fnp3041->SetLineColor(2);
   fnp3041->SetLineWidth(2);
   fnp3041->SetChisquare(29.31228);
   fnp3041->SetNDF(58);
   fnp3041->GetXaxis()->SetLabelFont(42);
   fnp3041->GetXaxis()->SetTitleOffset(1);
   fnp3041->GetXaxis()->SetTitleFont(42);
   fnp3041->GetYaxis()->SetLabelFont(42);
   fnp3041->GetYaxis()->SetTitleFont(42);
   fnp3041->SetParameter(0,0.03962235);
   fnp3041->SetParError(0,0.0004727209);
   fnp3041->SetParLimits(0,0,0);
   fnp3041->SetParameter(1,0.2349185);
   fnp3041->SetParError(1,0.01867267);
   fnp3041->SetParLimits(1,0,0);
   fnp3041->SetParent(grae);
   grae->GetListOfFunctions()->Add(fnp3041);
   
   ptstats = new TPaveStats(0.59,0.04000001,0.95,0.2,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf = 29.31 / 58");
   ptstats_LaTex = ptstats->AddText("p0       = 0.03962 #pm 0.0004727 ");
   ptstats_LaTex = ptstats->AddText("p1       = 0.2349 #pm 0.01867 ");
   ptstats->SetOptStat(0);
   ptstats->SetOptFit(111);
   ptstats->Draw();
   grae->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(grae->GetListOfFunctions());
   grae->Draw("ap");
   
   pt = new TPaveText(0.4166171,0.9291672,0.5833829,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt_LaTex = pt->AddText("Graph");
   pt->Draw();
   c3_2->Modified();
   c3->cd();
  
// ------------>Primitives in pad: c3_3
   TPad *c3_3 = new TPad("c3_3", "c3_3",0.6766667,0.01,0.99,0.99);
   c3_3->Draw();
   c3_3->cd();
   c3_3->Range(-13.675,-0.496792,123.075,4.471128);
   c3_3->SetFillColor(0);
   c3_3->SetBorderMode(0);
   c3_3->SetBorderSize(2);
   c3_3->SetFrameBorderMode(0);
   c3_3->SetFrameBorderMode(0);
   
   Double_t Graph0_fx3043[100] = {
   0.5,
   1.5,
   2.5,
   3.5,
   4.5,
   5.5,
   6.5,
   7.5,
   8.5,
   9.5,
   10.5,
   11.5,
   12.5,
   13.5,
   14.5,
   15.5,
   16.5,
   17.5,
   18.5,
   19.5,
   20.5,
   21.5,
   22.5,
   23.5,
   24.5,
   25.5,
   26.5,
   27.5,
   28.5,
   29.5,
   30.5,
   31.5,
   32.5,
   33.5,
   34.5,
   35.5,
   36.5,
   37.5,
   38.5,
   39.5,
   40.5,
   41.5,
   42.5,
   43.5,
   44.5,
   45.5,
   46.5,
   47.5,
   48.5,
   49.5,
   50.5,
   51.5,
   52.5,
   53.5,
   54.5,
   55.5,
   56.5,
   57.5,
   58.5,
   59.5,
   60.5,
   61.5,
   62.5,
   63.5,
   64.5,
   65.5,
   66.5,
   67.5,
   68.5,
   69.5,
   70.5,
   71.5,
   72.5,
   73.5,
   74.5,
   75.5,
   76.5,
   77.5,
   78.5,
   79.5,
   80.5,
   81.5,
   82.5,
   83.5,
   84.5,
   85.5,
   86.5,
   87.5,
   88.5,
   89.5,
   90.5,
   91.5,
   92.5,
   93.5,
   94.5,
   95.5,
   96.5,
   97.5,
   98.5,
   99.5};
   Double_t Graph0_fy3043[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   1.4,
   1.4,
   1.2,
   0.6,
   0.9,
   1,
   1,
   1.3,
   1.5,
   1.5,
   1.4,
   1.4,
   1.4,
   1.5,
   1.5,
   1.6,
   1.6,
   1.7,
   1.7,
   1.7,
   1.8,
   1.8,
   1.8,
   1.8,
   1.9,
   1.9,
   1.9,
   1.9,
   2,
   2,
   2,
   2,
   2,
   2.1,
   2.1,
   2.1,
   2.1,
   2.1,
   2.2,
   2.2,
   2.2,
   2.2,
   2.2,
   2.2,
   2.2,
   2.3,
   2.3,
   2.3,
   2.3,
   2.4,
   2.4,
   2.3,
   2.3,
   2.4,
   2.4,
   2.4,
   2.4,
   2.4,
   2.5,
   2.5,
   2.5,
   2.5,
   2.6,
   2.6,
   2.7,
   2.6,
   2.5,
   2.7,
   2.6,
   2.6,
   2.9,
   2.8,
   2.5,
   1.2,
   1.4,
   2,
   2.1,
   1.9,
   1.9,
   1.9,
   1.8,
   1.6,
   1.6,
   0.6,
   0.5,
   0,
   0,
   0.3,
   0.3,
   0};
   Double_t Graph0_felx3043[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fely3043[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   1.4,
   0.2770784,
   1.2,
   0.6,
   0.3728915,
   0.1408852,
   0.1660202,
   0.3177835,
   0.1201945,
   0.2688016,
   0.1529269,
   0.05966357,
   0.05565208,
   0.06157035,
   0.03592119,
   0.1290271,
   0.05747002,
   0.1226633,
   0.09558477,
   0.04880421,
   0.1063006,
   0.06997755,
   0.04878371,
   0.01974023,
   0.09420384,
   0.09259557,
   0.07282666,
   0.02663628,
   0.1045265,
   0.08267808,
   0.03711166,
   0.03636049,
   0.03479102,
   0.1045984,
   0.07981584,
   0.07843525,
   0.0695614,
   0.03027511,
   0.103378,
   0.08111601,
   0.08722395,
   0.06453219,
   0.04905742,
   0.03962727,
   0.02996445,
   0.1174569,
   0.07230915,
   0.07305564,
   0.05115168,
   0.0979135,
   0.1097788,
   0.09248778,
   0.06541695,
   0.1155903,
   0.1314493,
   0.103414,
   0.07290649,
   0.1156536,
   0.1909057,
   0.1460681,
   0.1648852,
   0.136684,
   0.2271854,
   0.3061979,
   0.3359125,
   0.2251728,
   0.142218,
   0.3465071,
   0.3293812,
   0.3191381,
   0.530057,
   0.6113362,
   0.6341186,
   0.165644,
   0.3827746,
   0.6487299,
   0.6509387,
   0.8677458,
   0.523919,
   0.4038399,
   0.8007107,
   0.3816228,
   0.3493662,
   0.13,
   0.5,
   0,
   0,
   0.3,
   0.3,
   0};
   Double_t Graph0_fehx3043[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fehy3043[100] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   -1.4,
   -0.01797959,
   -0.01464102,
   0.04333333,
   0.05440163,
   0.1846899,
   0.2183673,
   0.1149012,
   0.1110457,
   0.0312548,
   0.01345907,
   0.06705825,
   0.05570184,
   0.03039247,
   0.03915578,
   -0.06142848,
   0.007390324,
   -0.06817786,
   -0.04631704,
   -0.007660915,
   -0.06170267,
   -0.0281555,
   -0.0135285,
   0.01283891,
   -0.06144882,
   -0.06057986,
   -0.04231033,
   0.001986174,
   -0.07541584,
   -0.05472547,
   -0.01074873,
   -0.009633664,
   -0.007349935,
   -0.07477951,
   -0.04899875,
   -0.04756839,
   -0.03861519,
   0.000634486,
   -0.0700327,
   -0.04582368,
   -0.04944537,
   -0.0249613,
   -0.00909529,
   0.001637864,
   0.0132031,
   -0.0656795,
   -0.02213166,
   -0.01651678,
   0.01327238,
   -0.01986884,
   -0.01040208,
   -0.0001473817,
   0.02553415,
   -0.02396584,
   -0.02254083,
   0.003570538,
   0.04758852,
   0.05254296,
   0.01366249,
   0.07660736,
   0.06448244,
   0.1136294,
   0.06846339,
   0.09355522,
   0.1336927,
   0.08863387,
   0.1705679,
   0.2158553,
   0.4668446,
   0.4148011,
   0.08085798,
   0.2614575,
   1.113033,
   0.5775435,
   0.5029834,
   0.7137828,
   0.1886763,
   0.3364171,
   0.831294,
   0.8379607,
   -0.02832816,
   -0.03,
   -0.0263932,
   -0.015,
   -0.5,
   0,
   0,
   -0.3,
   -0.3,
   0};
   grae = new TGraphAsymmErrors(100,Graph0_fx3043,Graph0_fy3043,Graph0_felx3043,Graph0_fehx3043,Graph0_fely3043,Graph0_fehy3043);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(24);
   grae->SetMarkerSize(0.4);
   
   TH1F *Graph_Graph03043 = new TH1F("Graph_Graph03043","Graph",100,0,109.4);
   Graph_Graph03043->SetMinimum(0);
   Graph_Graph03043->SetMaximum(3.974336);
   Graph_Graph03043->SetDirectory(0);
   Graph_Graph03043->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph03043->SetLineColor(ci);
   Graph_Graph03043->GetXaxis()->SetTitle("#rho");
   Graph_Graph03043->GetXaxis()->SetLabelFont(42);
   Graph_Graph03043->GetXaxis()->SetTitleOffset(1);
   Graph_Graph03043->GetXaxis()->SetTitleFont(42);
   Graph_Graph03043->GetYaxis()->SetTitle("Charge Isolation 95% Contour ");
   Graph_Graph03043->GetYaxis()->SetLabelFont(42);
   Graph_Graph03043->GetYaxis()->SetTitleFont(42);
   Graph_Graph03043->GetZaxis()->SetLabelFont(42);
   Graph_Graph03043->GetZaxis()->SetTitleOffset(1);
   Graph_Graph03043->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03043);
   
   
   TF1 *fnc3044 = new TF1("fnc","[0]*x + [1]",20,80, TF1::EAddToList::kNo);
   fnc3044->SetFillColor(19);
   fnc3044->SetFillStyle(0);
   fnc3044->SetLineColor(2);
   fnc3044->SetLineWidth(2);
   fnc3044->SetChisquare(150.1212);
   fnc3044->SetNDF(58);
   fnc3044->GetXaxis()->SetLabelFont(42);
   fnc3044->GetXaxis()->SetTitleOffset(1);
   fnc3044->GetXaxis()->SetTitleFont(42);
   fnc3044->GetYaxis()->SetLabelFont(42);
   fnc3044->GetYaxis()->SetTitleFont(42);
   fnc3044->SetParameter(0,0.02012828);
   fnc3044->SetParError(0,0.0002289483);
   fnc3044->SetParLimits(0,0,0);
   fnc3044->SetParameter(1,1.062124);
   fnc3044->SetParError(1,0.01407823);
   fnc3044->SetParLimits(1,0,0);
   fnc3044->SetParent(grae);
   grae->GetListOfFunctions()->Add(fnc3044);
   
   ptstats = new TPaveStats(0.59,0.04000001,0.95,0.2,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf = 150.1 / 58");
   ptstats_LaTex = ptstats->AddText("p0       = 0.02013 #pm 0.0002289 ");
   ptstats_LaTex = ptstats->AddText("p1       = 1.062 #pm 0.01408 ");
   ptstats->SetOptStat(0);
   ptstats->SetOptFit(111);
   ptstats->Draw();
   grae->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(grae->GetListOfFunctions());
   grae->Draw("ap");
   
   pt = new TPaveText(0.4166171,0.9291672,0.5833829,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt_LaTex = pt->AddText("Graph");
   pt->Draw();
   c3_3->Modified();
   c3->cd();
   c3->Modified();
   c3->cd();
   c3->SetSelected(c3);
}
