#include "AlgoKcon.h"

AlgoKcon::AlgoKcon( TString channel, int NBTag, double massL, double massH  ){

  ptype  = ".gif";

  cname = channel;
  hname = channel+"TM";

  mL = massL ;
  mH = massH ;

  //fitFunc  = new MassFitFunction();
  fitInput = new MassAnaInput( channel, NBTag, massL, massH );
  fitTools = new MassAna( channel, NBTag, massL, massH );
  fitInput->Initialize( &hfolder );
}

AlgoKcon::~AlgoKcon(){

 //delete fitFunc;
 delete fitInput;
 delete fitTools;

}


// Fit the data with kinematic constrain
void AlgoKcon::ConstrainFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag ){

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
  TH1D* dth1     = new TH1D("dth1","", nbin, mL, mH );
  TH1D* dth2     = new TH1D("dth2","", nbin, mL, mH );
  TH1D* dth3     = new TH1D("dth3","", nbin, mL, mH );
  TH1D* dth4     = new TH1D("dth4","", nbin, mL, mH );
  TH1D* dth123   = new TH1D("dt123","", nbin, mL, mH );

  getFakeData( mName, fakedata, ttstk, NULL, dth1, dth2, dth3, dth4 );

  double statErr = 0;
  double bestMass = fitTools->Chi2Test(mName, fakedata, lowBound, upBound ,9 , NBTag, &statErr, rbin, true );
  // Fit the data
  TF1 *func7 = new TF1("func7",MassFitFunction::fitData2, lowBound, upBound,12);

  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 9, NBTag );

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
  leg->AddEntry(dth1, "ttbar-all", "f");
  leg->AddEntry(dth2, "wjets", "f");
  leg->AddEntry(dth3, "singleTop_t", "f");
  leg->AddEntry(dth4, "singleTop_tW", "f");
  leg->Draw("same");

  c3->Update();
 
  // Draw the expected signal
  Double_t apars[12];
  double m1 = fitTools->MassDigi(mName);
  for (int i=0; i< 12; i++) { 
      apars[i]   = func7->GetParameter(i);
  }
  fprintf(logfile," %.2f %.2f %.2f \n", m1, apars[1], statErr );

  TF1 *func2 = new TF1("func2",MassFitFunction::fitSG1, 100, 360, 6 );
  for (int j=0; j<6; j++ )  func2->FixParameter( j, apars[j] );
  
  func2->SetLineColor(1);
  func2->SetLineWidth(3);

  c3->cd(2);
  dth1->Draw();
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
  dth123->Add(dth2, 1);
  dth123->Add(dth3, 1);
  dth123->Add(dth4, 1);
  dth123->SetFillColor(kOrange+7);
  dth123->Draw();
  func3->Draw("sames");
  c3->Update();

  //plot2 = "BTag_";
  c3->Print(hfolder+plot2+mName+ptype);

  delete func7;
  delete func2;
  delete func3;
  delete leg;
  delete ttstk;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth4;
  delete dth123;
  delete fakedata;
 
  fclose(logfile);
}

// combined the fake data
void AlgoKcon::getFakeData( TString mName, TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, TH1D* dth4, TH1D* dth5 ){

  // get the file names of fake data
  TString fNameList[5]; 
  fitInput->GetFileName( mName, 0, fNameList );

  // get all fake data 
  if (dth1 != NULL && dth0 == NULL ) fitInput->get_h1Obj( fNameList[0], "solTt", hname, dth1 ); // for kinematic constrain  
  if (dth2 != NULL ) fitInput->get_h1Obj( fNameList[1], "solTt", hname, dth2 );            // w+jets
  if (dth3 != NULL ) fitInput->get_h1Obj( fNameList[2], "solTt", hname, dth3,  0.195 );    // single top t-channel
  if (dth4 != NULL ) fitInput->get_h1Obj( fNameList[3], "solTt", hname, dth4,  0.183 );    // single top tW-channel
  if (dth5 != NULL ) fitInput->get_h1Obj( fNameList[4], "solTt", hname, dth5,  2.01  );    // QCD

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

void AlgoKcon::SetFitParameters( double mass, Double_t* para, int nPara, int NBTag ) {

     // Lep: value for No BTagging
     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau Wjets , 10: tt-ttwrong ratio , 11: tt-Wjets ratio 
     //                      p0      p1      p2     p3      p4     p5       p6      p7       p8     p9      p10      p11
     Double_t f0[12] = {  41.246,  -1., -6.505,  2.336,   0.0,  31.796, 107.098, -0.886,  151.034,  38.422,  38.950,  13.727 };
     Double_t f1[12] = {   0.092,   1.,  0.196, -0.007,   1.0,  -0.012,   0.316,  0.204,    0.023,  -0.003,  -0.002,  -0.020 };

     // *** for kinematic had: 
     // 10 GeV Bin - hadW Constrain
     //Double_t kk0[12] = {  27.869,  13.838, 22.286, 0.155, 148.285, 53.366, 184.916, 52.355,   0.,    0.,  2.422,  0. };
     //Double_t kk1[12] = {   0.071,   0.922,  0.025, 0.017,   0.303,  0.018,   0.000,  0.000,   0.,    0., -0.003,  0. };
     // 10 GeV Bin - dM Constrain
     //Double_t kk0[12] = {  1.904, 12.051, 5.438,  4.202, 117.285, 21.724, 195.073, 39.88,   0.,    0.,  3.820,  0. };
     //Double_t kk1[12] = {  0.258,  0.933, 0.139, -0.010,   0.653,  0.115,   0.000,  0.00,   0.,    0., -0.011,  0. };
     // 10 GeV Bin - dM*hadW Constrain
     Double_t kk0[12] = {  62.167, 4.508, 2.097, -3.413, 161.463, 10.687, 183.757, 37.475,   0.,    0.,  0.919,  0. };
     Double_t kk1[12] = {  -0.099, 0.969, 0.142,  0.038,   0.205,  0.134,   0.000,  0.000,   0.,    0.,  0.007,  0. };
     for ( int i =0; i < 12; i++ ) {
         if ( cname == "had" ) {
            if ( nPara == 9 && NBTag == -1 ) para[i] = kk0[i] + kk1[i]*mass ;
         }
         if ( cname == "lep" ) {
            if ( nPara == 9 && NBTag == -1 ) para[i] = f0[i] + f1[i]*mass ;
         }
     }

}

