#include "AlgoZero.h"

AlgoZero::AlgoZero(){

  ptype  = ".gif";

  fitInput = new MassAnaInput();

  string decayType;
  fitInput->GetParameters( "DecayType", &decayType );
  fitInput->GetParameters( "MassLBound", &mL);
  fitInput->GetParameters( "MassHBound", &mH);

  cname = decayType;

  fitTools = new MassAna();
  fitInput->Initialize( &hfolder );

}

AlgoZero::~AlgoZero(){

 delete fitInput;
 delete fitTools;

}

// separate background : tt-wrong permutation,  wjets + qcd
void AlgoZero::MoreCombinedFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag ){

  FILE* logfile = fopen(hfolder+"/Outputf.log","a"); 

  gStyle->SetOptStat("i");
  gStyle->SetOptFit(111);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatTextColor(1);
  c7 = new TCanvas("c7","", 900, 800);
  c7->SetFillColor(10);
  c7->SetFillColor(10);
  c7->Divide(2,2);

  int nbin = ( mH - mL )/rbin ;

  // Get fake data information
  THStack* ttstk = new THStack("ttstk", "Combined Fitting"); 
  TH1D* fakedata = new TH1D("fakedata","", nbin, mL, mH );
  /*
  TH1D* dth0 = new TH1D("dth0","", nbin, mL, mH );      // tt-signal
  TH1D* dth1 = new TH1D("dth1","", nbin, mL, mH );      // tt-wrong combination
  TH1D* dth2 = new TH1D("dth2","", nbin, mL, mH );      // w+jets
  TH1D* dth3 = new TH1D("dth3","", nbin, mL, mH );      // single top t-channel
  TH1D* dth4 = new TH1D("dth4","", nbin, mL, mH );      // single top tW-channel
  */
  THStack* dthb1 = new THStack("dthb1", "background group1 ");
  THStack* dthb2 = new THStack("dthb2", "background group2 ");
  TH1D* dtadd1 = new TH1D("dtadd1","", nbin, mL, mH );      
  TH1D* dtadd2 = new TH1D("dtadd2","", nbin, mL, mH );     

  //getFakeData( mName, fakedata, ttstk, dth0, dth1, dth2, dth3, dth4 );
  vector<TH1D*> hlist;
  fitInput->getFakeData( rbin, fakedata, ttstk, hlist );

  // pre-fit background
  double statErr = 0;
  double bestMass = fitTools->Chi2Test(mName, fakedata, lowBound, upBound, 12, NBTag, &statErr, nbin );
  // Fit the data
  TF1 *func1 = new TF1("func1",MassFitFunction::fitData1, lowBound, upBound,12);
  if ( cname == "lep" ) func1 = new TF1("func1",MassFitFunction::fitData2, lowBound, upBound,12);

  // get the extrapolated parameters
  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 12, NBTag, rbin );

  func1->SetParLimits(0,  5., 100.);
  func1->FixParameter(1, fpar[1] );
  func1->FixParameter(2, fpar[2] );
  func1->SetParLimits(3, fpar[3]- 0.1*fpar[3], fpar[3]+0.1*fpar[3] );
  func1->FixParameter(4, fpar[4] );
  func1->FixParameter(5, fpar[5] );
  func1->FixParameter(6, fpar[6] );
  func1->FixParameter(7, fpar[7] );
  func1->FixParameter(8, fpar[8] );
  func1->FixParameter(9, fpar[9] );
  func1->FixParameter(10, fpar[10] );
  func1->SetParLimits(11, fpar[11]- 0.1*fpar[11], fpar[11]+0.1*fpar[11] );

  c7->cd(1);
  ttstk->Draw();
  func1->SetLineColor(1);
  func1->SetLineWidth(2);
  fakedata->Fit( func1, "R","sames",110,330);

  TLegend *leg = new TLegend(.65, .4, .95, .75);
  vector<string> channelNames ;
  fitInput->GetParameters("channel", &channelNames);
  leg->AddEntry(hlist[0], "tt-correct", "f");
  for (size_t i=1; i< hlist.size(); i++) {
      TString legName = channelNames[i] ;
      leg->AddEntry( hlist[i], legName , "f");
  }
  /*
  leg->AddEntry(dth0, "tt-correct", "f");
  leg->AddEntry(dth1, "tt-wrong", "f");
  leg->AddEntry(dth2, "wjets", "f");
  leg->AddEntry(dth3, "SingleTop_t", "f");
  leg->AddEntry(dth4, "SingleTop_tW", "f");
  */
  leg->Draw("same");

  c7->Update();

  Double_t apars[12];
  double m1 = fitTools->MassDigi(mName);
  for (int i=0; i< 12; i++) { 
      apars[i]   = func1->GetParameter(i);
  }
  fprintf(logfile," %.2f %.2f %.2f \n", m1, apars[1], statErr );

  // Draw the expected signal
  c7->cd(2);
  TF1 *func2 = new TF1("func2",MassFitFunction::fitSG, 100, 360, 6 );
  if ( cname == "lep" ) func2 = new TF1("func2",MassFitFunction::fitSG1, 100, 360, 6 );
  for (int j=0; j<6; j++ ) {
      func2->FixParameter( j, apars[j] );
  }
  func2->SetLineColor(1);
  func2->SetLineWidth(3);
  //dth0->Draw();
  hlist[0]->Draw();
  func2->Draw("sames");
  c7->Update();

  // Draw the expected background , background group 1
  c7->cd(3);
  TF1 *func3 = new TF1("func3",MassFitFunction::fitLD, 100, 360, 3);
  func3->FixParameter(0, apars[0]*apars[10] );
  func3->FixParameter(1, apars[6] );
  func3->FixParameter(2, apars[7] );
  func3->SetLineColor(1);
  func3->SetLineWidth(3);

  dthb1->Add( hlist[0] );
  dtadd1->Add( hlist[0]);
  
  /*
  if ( comp[4] ) dthb1->Add( dth4 );
  if ( comp[3] ) dthb1->Add( dth3 );
  if ( comp[2] ) dthb1->Add( dth2 );
  if ( comp[1] ) dthb1->Add( dth1 );
  if ( comp[4] ) dtadd1->Add( dth4 );
  if ( comp[3] ) dtadd1->Add( dth3 );
  if ( comp[2] ) dtadd1->Add( dth2 );
  if ( comp[1] ) dtadd1->Add( dth1 );
  */
  dthb1->Draw();
  dtadd1->Draw("sames");

  func3->Draw("sames");

  c7->Update();

  // Draw the expected background , background group 2
  c7->cd(4);
  TF1 *func4 = new TF1("func4",MassFitFunction::fitLD, 100, 360, 3);
  func4->FixParameter(0, apars[0]*apars[11] );
  func4->FixParameter(1, apars[8] );
  func4->FixParameter(2, apars[9] );
  func4->SetLineColor(1);
  func4->SetLineWidth(3);
  
  for (size_t i= hlist.size(); i>1; i--) {
      dthb2->Add( hlist[i-1] );
      dtadd2->Add( hlist[i-1]);
  }
  /*
  if ( !comp[4] ) dthb2->Add( dth4 );
  if ( !comp[3] ) dthb2->Add( dth3 );
  if ( !comp[2] ) dthb2->Add( dth2 );
  if ( !comp[1] ) dthb2->Add( dth1 );
  if ( !comp[4] ) dtadd2->Add( dth4 );
  if ( !comp[3] ) dtadd2->Add( dth3 );
  if ( !comp[2] ) dtadd2->Add( dth2 );
  if ( !comp[1] ) dtadd2->Add( dth1 );
  */
  dthb2->Draw();
  dtadd2->Draw("sames");
  func4->Draw("sames");
  c7->Update();

  c7->Print(hfolder+mName+"Zero"+ptype);

  fclose(logfile);

  delete func1;
  delete func2;
  delete func3;
  delete func4;
  delete leg;
  delete ttstk;
  /*
  delete dth0;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth4;
  */
  delete dthb1;
  delete dthb2;
  delete dtadd1;
  delete dtadd2;
  delete fakedata;

}

// combined all background : tt-wrong permutation, wjets and qcd
void AlgoZero::CombinedFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag ){

  FILE* logfile = fopen(hfolder+"/Outputf.log","a"); 

  gStyle->SetOptStat("i");
  gStyle->SetOptFit(111);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatTextColor(1);
  c3 = new TCanvas("c3","", 900, 800);
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->Divide(2,2);

  // Get fake data information
  int nbin = ( mH - mL )/ rbin ;
  THStack* ttstk = new THStack("ttstk", "Combined Fitting"); 
  TH1D* fakedata = new TH1D("fakedata","", nbin, mL, mH );
  /*
  TH1D* dth0 = new TH1D("dth0","", nbin, mL, mH );
  TH1D* dth1 = new TH1D("dth1","", nbin, mL, mH );
  TH1D* dth2 = new TH1D("dth2","", nbin, mL, mH );
  TH1D* dth3 = new TH1D("dth3","", nbin, mL, mH );
  TH1D* dth4 = new TH1D("dth4","", nbin, mL, mH );
  getFakeData( mName, fakedata, ttstk, dth0, dth1, dth2, dth3, dth4 );
  */
  TH1D* dth123 = new TH1D("dt123","", nbin, mL, mH );
  vector<TH1D*> hlist;
  fitInput->getFakeData( rbin, fakedata, ttstk, hlist );

  double statErr = 0;
  double bestMass = fitTools->Chi2Test(mName, fakedata, lowBound, upBound ,9 , NBTag, &statErr, rbin );
  // Fit the data
  TF1 *func7 = new TF1("func7",MassFitFunction::fitData1, lowBound, upBound,12);
  if ( cname == "lep" ) func7 = new TF1("func7",MassFitFunction::fitData2, lowBound, upBound, 12);

  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 9, NBTag, rbin );

  func7->SetParLimits(0,  5., 200.);
  func7->FixParameter(1, fpar[1] );
  func7->FixParameter(2, fpar[2] );
  func7->SetParLimits(3, fpar[3]- 0.1*fpar[3], fpar[3]+0.1*fpar[3] );
  func7->FixParameter(4, fpar[4] );
  func7->FixParameter(5, fpar[5] );
  func7->FixParameter(6, fpar[6] );
  func7->FixParameter(7, fpar[7] );
  func7->FixParameter(8, fpar[8] );
  func7->FixParameter(9, fpar[9] );
  func7->SetParLimits(10, fpar[10]- 0.1*fpar[10], fpar[10]+0.1*fpar[10] );
  func7->FixParameter(11, fpar[11] );

  c3->cd(1);
  ttstk->Draw();

  func7->SetLineColor(kBlue+2);
  func7->SetLineStyle(1);
  func7->SetLineWidth(2);
  fakedata->Fit( func7, "R","sames",lowBound,upBound);
  func7->Draw("same");

  TLegend *leg = new TLegend(.65, .4, .95, .75);
  vector<string> channelNames ;
  fitInput->GetParameters("channel", &channelNames);
  leg->AddEntry(hlist[0], "tt-correct", "f");
  for (size_t i=1; i< hlist.size(); i++) {
      TString legName = channelNames[i] ;
      leg->AddEntry( hlist[i], legName , "f");
  }
  /*
  leg->AddEntry(dth0, "tt-correct", "f");
  leg->AddEntry(dth1, "tt-wrong", "f");
  leg->AddEntry(dth2, "wjets", "f");
  leg->AddEntry(dth3, "SingleTop_t", "f");
  leg->AddEntry(dth4, "SingleTop_tW", "f");
  */
  leg->Draw("same");

  c3->Update();

  // Draw the expected signal
  Double_t apars[12];
  double m1 = fitTools->MassDigi(mName);
  for (int i=0; i< 12; i++) { 
      apars[i]   = func7->GetParameter(i);
  }
  fprintf(logfile," %.2f %.2f %.2f \n", m1, apars[1], statErr );

  TF1 *func2 = new TF1("func2",MassFitFunction::fitSG, 100, 360, 6 );
  if ( cname == "lep" ) func2 = new TF1("func2",MassFitFunction::fitSG1, 100, 360, 6 );
  for (int j=0; j<6; j++ ) {
      func2->FixParameter( j, apars[j] );
  }
  func2->SetLineColor(1);
  func2->SetLineWidth(3);

  c3->cd(2);
  //dth0->Draw();
  hlist[0]->Draw();
  func2->Draw("sames");
  c3->Update();

  // Draw the expected background
  TF1 *func3 = new TF1("func3",MassFitFunction::fitLD, 100, 360, 3);
  func3->FixParameter(0, apars[10]*apars[0] );
  func3->FixParameter(1, apars[6] );
  func3->FixParameter(2, apars[7] );
  func3->SetLineColor(1);
  //func3->SetLineStyle(2);
  func3->SetLineWidth(3);

  c3->cd(3);
  for (size_t i= hlist.size(); i>1; i--) {
      dth123->Add( hlist[i-1] );
  }
  /*
  dth123->Add(dth1, 1);
  dth123->Add(dth2, 1);
  dth123->Add(dth3, 1);
  dth123->Add(dth4, 1);
  */
  dth123->SetFillColor(kOrange+7);
  dth123->Draw();
  func3->Draw("sames");
  c3->Update();

  c3->Print(plot2+mName+ptype);

  delete func7;
  delete func2;
  delete func3;
  delete leg;
  delete ttstk;
  /*
  delete dth0;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth4;
  */
  delete dth123;
  delete fakedata;
 
  fclose(logfile);
}

// combined the fake data
/*
void AlgoZero::getFakeData( TString mName,  TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, TH1D* dth4, TH1D* dth5 ){

  // get the file names of fake data
  TString fNameList[5];
  fitInput->GetFileName( mName, 0, fNameList );
  // get all fake data 
  if ( cname == "had" ) {
     if (dth0 != NULL && dth1 != NULL ) { 
        fitInput->get_h1Obj( fNameList[0],  "mcmTt", "hadTM", dth0 );                  // tt-signal
        fitInput->getHadPermutation( fNameList[0], "hadTM", dth1 ) ;
        dth1->Add( dth0, -1 ) ;
     }
     if (dth2 != NULL ) fitInput->getHadPermutation( fNameList[1], "hadTM", dth2 );  // w+jets 
     if (dth3 != NULL ) fitInput->getHadPermutation( fNameList[2], "hadTM", dth3, 0.195 );  // single top t channel
     if (dth4 != NULL ) fitInput->getHadPermutation( fNameList[3], "hadTM", dth4, 0.183 );  // single top tW channel 
     if (dth5 != NULL ) fitInput->getHadPermutation( fNameList[4], "hadTM", dth5, 2.01  );  // QCD
  }
  if ( cname == "lep" ) {
     if (dth0 != NULL && dth1 != NULL ) { 
        fitInput->get_h1Obj( fNameList[0],  "mcmTt", "lepTM", dth0 );                  // tt-signal
        fitInput->getLepPermutation( fNameList[0], "lepTM", dth1 ) ;
        dth1->Add( dth0, -1 ) ;
     }
     if (dth2 != NULL ) fitInput->getLepPermutation( fNameList[1], "hadTM", dth2 );  // w+jets 
     if (dth3 != NULL ) fitInput->getLepPermutation( fNameList[2], "hadTM", dth3, 0.195 );  // single top t channel
     if (dth4 != NULL ) fitInput->getLepPermutation( fNameList[3], "hadTM", dth4, 0.183 );  // single top tW channel 
     if (dth5 != NULL ) fitInput->getLepPermutation( fNameList[4], "hadTM", dth5, 2.01  );  // QCD
  }

  // mix all fake samples
  if (dth0 != NULL ) ttadd->Add(dth0, 1);
  if (dth1 != NULL ) ttadd->Add(dth1, 1);
  if (dth2 != NULL ) ttadd->Add(dth2, 1);
  if (dth3 != NULL ) ttadd->Add(dth3, 1);
  if (dth4 != NULL ) ttadd->Add(dth4, 1);
  if (dth5 != NULL ) ttadd->Add(dth5, 1);

  // give different color for samples
  if (dth5 != NULL ) dth5->SetFillColor(5);  // QCD
  if (dth4 != NULL ) dth4->SetFillColor(6);  // single top tw channel
  if (dth3 != NULL ) dth3->SetFillColor(4);  // single top t channel
  if (dth2 != NULL ) dth2->SetFillColor(2);  // w+jets
  if (dth1 != NULL ) dth1->SetFillColor(7);
  if (dth0 != NULL ) dth0->SetFillColor(3);

  // stack them
  if (dth5 != NULL ) ttstk->Add( dth5 );
  if (dth4 != NULL ) ttstk->Add( dth4 );
  if (dth3 != NULL ) ttstk->Add( dth3 );
  if (dth2 != NULL ) ttstk->Add( dth2 );
  if (dth1 != NULL ) ttstk->Add( dth1 );
  if (dth0 != NULL ) ttstk->Add( dth0 );

}
*/

void AlgoZero::SetFitParameters( double mass, Double_t* para, int nPara, int NBTag, int rbin ) {

     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau Wjets , 10: tt-ttwrong ratio , 11: tt-Wjets ratio 
     //                      p0      p1      p2     p3      p4     p5       p6      p7       p8     p9      p10      p11
     Double_t a0[12] =  { 12.389, 15.65, -10.689, 53.002,  4.379,  5.061, 92.660, 23.041, 196.402, 60.891, 29.205,  7.971 };
     Double_t a1[12] =  {  0.288,  0.91,    0.19, -0.099,  0.006, -0.001,  0.556,  0.131,  0.000,   0.000, -0.061, -0.021 };
     // Lep: value for No BTagging
     Double_t f0[12] = {  41.246,  -1., -6.505,  2.336,   0.0,  31.796, 107.098, -0.886,  151.034,  38.422,  38.950,  13.727 };
     Double_t f1[12] = {   0.092,   1.,  0.196, -0.007,   1.0,  -0.012,   0.316,  0.204,    0.023,  -0.003,  -0.002,  -0.020 };

     // Had: value for All BTagging
     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau/all Bg , 10: tt-ttwrong ratio , 11: tt-allBg ratio 
     //                      p0     p1      p2     p3     p4     p5      p6      p7      p8      p9     p10      p11
     Double_t b0[12] = { -3.495, 7.507, -13.302, 66.885, 4.379,  5.177, 97.443, 17.985,   0.,    0.,  29.546,  0. };
     Double_t b1[12] = {  0.284, 0.960,   0.202, -0.184, 0.006, -0.001,  0.522,  0.156,   0.,    0.,  -0.076,  0. };
     // Lep: value for All BTagging
     Double_t g0[12] = {  -4.506,  -1., 12.513,  1.185,  0.0,  31.356, 110.244, -5.392,   0.,    0.,  62.543,  0. };
     Double_t g1[12] = {   0.268,   1.,  0.084,  0.000,  1.0,  -0.009,   0.308,  0.225,   0.,    0.,  -0.170,  0. };

     // Had: value for 1 BTagging
     Double_t c0[12] = {  8.225, 10.816, -13.0, 43.413, 4.379,  5.285, 101.584, 16.580,   0.,    0.,  29.127,  0. };
     Double_t c7[12] = {  0.128,  0.939,   0.2, -0.049, 0.006, -0.002,   0.500,  0.162,   0.,    0.,  -0.054,  0. };
     // Lep: value for 1 BTagging
     Double_t h0[12] = {  12.558,  -1.498, 11.048,  1.185,  0.0,  32.059, 121.560, -1.993,   0.,    0.,  32.521,  0. };
     Double_t h1[12] = {   0.096,   1.003,  0.214,  0.000,  1.0,  -0.013,   0.232,  0.209,   0.,    0.,   0.042,  0. };

     // value for 2 BTagging
     // *** Had:
     // 10 GeV Bin
     Double_t d0[12] = { -14.213, -0.377, -10.394, 76.862, 4.379,  5.327, 42.436, 18.694,  0.,  0., 23.793,  0. };
     Double_t d1[12] = {  -0.170,  1.005,   0.175, -0.259, 0.006, -0.002,  0.812,  0.142,  0.,  0., -0.083,  0. };
     // 20 GeV Bin
     Double_t j0[12] = { -20.639,  6.181, -11.187, 15.988, 4.374,  4.452, 33.826, -12.159,  0.,  0.,  9.555,   0. };
     Double_t j1[12] = {   0.291,  0.967,   0.185,  0.093, 0.006,  0.003,  0.870,   0.333,  0.,  0.,  0.000,   0. };
     // 20 GeV Bin, 12 parameters
     Double_t e0[12] = { -20.639,  6.181, -11.187, 15.988, 4.374,  4.452, 16.027, -21.647, 99.958, -1.391, 22.120, 2.501 };
     Double_t e1[12] = {   0.291,  0.967,   0.185,  0.093, 0.006,  0.003,  0.957,   0.376,  0.626,  0.299, -0.077, -0.01 };
     // *** Lep: 
     // 10 GeV Bin
     Double_t i0[12] = {  -0.791, -1.388, 41.876,  1.185,   0.0,  30.738, 59.454, -18.817,   0.,    0.,  17.506,  0. };
     Double_t i1[12] = {   0.085,  1.007, -0.114,  0.000,   1.0,  -0.004,  0.637,   0.282,   0.,    0.,   0.000,  0. };
     //20 GeV Bin
     Double_t k0[12] = {   6.194,  -5.040,   4.708,  1.185,  0.0,  25.131, 58.172, -20.989,   0.,    0.,  17.190,  0. };
     Double_t k1[12] = {   0.134,   1.026,   0.108,  0.000,  1.0,   0.027,  0.644,   0.302,   0.,    0.,   0.000,  0. };

     // *** for kinematic had: 
     // 10 GeV Bin - hadW Constrain
     //Double_t kk0[12] = {  27.869,  13.838, 22.286, 0.155, 148.285, 53.366, 184.916, 52.355,   0.,    0.,  2.422,  0. };
     //Double_t kk1[12] = {   0.071,   0.922,  0.025, 0.017,   0.303,  0.018,   0.000,  0.000,   0.,    0., -0.003,  0. };
     for ( int i =0; i < 12; i++ ) {
         if ( cname == "had" ) {
            if ( nPara == 12 && NBTag == -1 ) para[i] = a0[i] + a1[i]*mass ;
            if ( nPara ==  9 && NBTag ==  0 ) para[i] = b0[i] + b1[i]*mass ;
            if ( nPara ==  9 && NBTag ==  1 ) para[i] = c0[i] + c7[i]*mass ;
            if ( nPara ==  9 && NBTag ==  2 ) para[i] = d0[i] + d1[i]*mass ;
            if ( nPara ==  9 && NBTag ==  2 && rbin == 20 ) para[i] = j0[i] + j1[i]*mass ;
            if ( nPara == 12 && NBTag ==  2 ) para[i] = e0[i] + e1[i]*mass ;
         }
         if ( cname == "lep" ) {
            if ( nPara ==12 && NBTag == -1 ) para[i] = f0[i] + f1[i]*mass ;
            if ( nPara == 9 && NBTag ==  0 ) para[i] = g0[i] + g1[i]*mass ;
            if ( nPara == 9 && NBTag ==  1 ) para[i] = h0[i] + h1[i]*mass ;
            if ( nPara == 9 && NBTag ==  2 ) para[i] = i0[i] + i1[i]*mass ;
            if ( nPara == 9 && NBTag ==  2 && rbin == 20 ) para[i] = k0[i] + k1[i]*mass ;
         }
     }

}

