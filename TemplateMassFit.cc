#include "TemplateMassFit.h"

TemplateMassFit::TemplateMassFit( TString aname, TString channel ){


  // For different mass &/ different algorithm
  fname = aname ;
  // "Lep" or "Had"
  cname = channel ;
  hname = channel+"_selTmass_pt0";
  sname = "mcTtMass0";
  hfolder = "AllTags0/";


  plot1 = hfolder+"/"+channel+"_CombinedFit2_"+aname+".gif";
  plot2 = hfolder+"/"+channel+"_CombinedFit_"+aname+".gif";
  plot4 = hfolder+"/"+channel+"_TemplateFit_"+aname+".gif";
  plot5 = hfolder+"/"+channel+"_Norm2_"+aname+".gif";
  plot6 = hfolder+"/"+channel+"_TemplateFit2_"+aname+".gif";
  plot9 = hfolder+"/"+channel+"_X2Test"+aname+".gif";

  fitFunc = new MassFitFunction();

  gSystem->mkdir(hfolder);
  
}

TemplateMassFit::~TemplateMassFit(){
 
  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
  delete c6;
  delete c7;
  delete c8;
  delete c9; 
  delete fitFunc;

}

// separate background : tt-wrong permutation,  wjets + qcd
void TemplateMassFit::MoreCombinedFitting( int rbin, int lowBound, int upBound ){

  logfile = fopen(hfolder+"/Outputf.log","a"); 

  gStyle->SetOptStat("i");
  gStyle->SetOptFit(111);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatTextColor(1);
  c1 = new TCanvas("c1","", 900, 800);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);

  int nbin = 480./rbin ;

  // Get fake data information
  THStack* ttstk = new THStack("ttstk", "Combined Fitting"); 
  TH1D* fakedata = new TH1D("fakedata","", nbin, 0., 480.);
  TH1D* dth0 = new TH1D("dth0","", nbin, 0., 480.);
  TH1D* dth1 = new TH1D("dth1","", nbin, 0., 480.);
  TH1D* dth2 = new TH1D("dth2","", nbin, 0., 480.);
  TH1D* dth3 = new TH1D("dth3","", nbin, 0., 480.);
  TH1D* dth4 = new TH1D("dth4","", nbin, 0., 480.);
  TH1D* dth23 = new TH1D("dth23","", nbin, 0., 480.);

  getFakeData( rbin, fakedata, ttstk, dth0, dth1, dth2, dth3, dth4 );

  // pre-fit background
  //TString theFileName = "mgtt_"+mName+".root" ;
  fprintf(logfile," === Chi2 Test for   === \n" );
  double bestMass = Chi2Test(fakedata, lowBound, upBound ,12 );
  // Fit the data
  TF1 *func1 = new TF1("func1",MassFitFunction::fitData1, lowBound, upBound,12);

  // get the extrapolated parameters
  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 12, -1 );

  func1->SetParLimits(0,  25., 100.);
  func1->FixParameter(1, fpar[1] );
  func1->FixParameter(2, fpar[2] );
  func1->SetParLimits(3, fpar[3]- 0.1*fpar[3], fpar[3]+0.1*fpar[3] );
  func1->FixParameter(4, fpar[4] );
  func1->FixParameter(5, fpar[5] );
  func1->FixParameter(6, fpar[6] );
  func1->FixParameter(7, fpar[7] );
  func1->FixParameter(8, fpar[8] );
  func1->FixParameter(9, fpar[9] );
  func1->SetParLimits(10,fpar[10]- 0.1*fpar[10], fpar[10]+0.1*fpar[10] );
  func1->FixParameter(11, fpar[11] );

  c1->cd(1);
  ttstk->Draw();
  func1->SetLineColor(4);
  func1->SetLineStyle(2);
  func1->SetLineWidth(3);
  fakedata->Fit( func1, "R","sames",110,330);
  c1->Update();

  // Draw the expected signal
  fprintf(logfile," ---------- Fitting Parameters --------- \n" );
  Double_t apars[12];
  for (int i=0; i< 12; i++) { 
      apars[i]   = func1->GetParameter(i);
      fprintf(logfile," p%d = %f \n", i, apars[i] );
  }
  fprintf(logfile," =========== Fitting Done ========= \n" );
  fprintf(logfile," \n" );

  TF1 *func2 = new TF1("func2",MassFitFunction::fitSG, 100, 360, 6 );
  for (int j=0; j<6; j++ ) {
      func2->FixParameter( j, apars[j] );
  }
  func2->SetLineColor(1);
  func2->SetLineStyle(5);
  func2->SetLineWidth(3);
  c1->cd(2);
  dth0->Draw();
  func2->Draw("sames");
  c1->Update();

  // Draw the expected background
  TF1 *func3 = new TF1("func3",MassFitFunction::fitLD, 100, 360, 3);
  func3->FixParameter(0, apars[0]*apars[10] );
  func3->FixParameter(1, apars[6] );
  func3->FixParameter(2, apars[7] );
  func3->SetLineColor(1);
  func3->SetLineStyle(2);
  func3->SetLineWidth(3);
  c1->cd(3);
  dth1->Draw();
  func3->Draw("sames");
  c1->Update();

  // Draw the expected background
  TF1 *func4 = new TF1("func4",MassFitFunction::fitLD, 100, 360, 3);
  func4->FixParameter(0, apars[0]*apars[11] );
  func4->FixParameter(1, apars[8] );
  func4->FixParameter(2, apars[9] );
  func4->SetLineColor(1);
  func4->SetLineStyle(2);
  func4->SetLineWidth(3);
  c1->cd(4);
  dth23->Add(dth2, 1);
  dth23->Add(dth3, 1);
  dth23->Add(dth4, 1);
  dth23->SetFillColor(2);
  dth23->Draw();
  func4->Draw("sames");
  c1->Update();

  c1->Print(plot1);

  fclose(logfile);

  delete func1;
  delete func2;
  delete func3;
  delete func4;
  delete ttstk;
  delete dth0;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth4;
  delete dth23;
  delete fakedata;
}

// combined all background : tt-wrong permutation, wjets and qcd
void TemplateMassFit::CombinedFitting( int rbin, int lowBound, int upBound, int NBTag ){

  logfile = fopen(hfolder+"/Outputf.log","a"); 

  gStyle->SetOptStat("i");
  gStyle->SetOptFit(111);
  gStyle->SetStatY(0.99);
  gStyle->SetStatX(0.99);
  gStyle->SetStatTextColor(1);
  c3 = new TCanvas("c3","", 900, 800);
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->Divide(2,2);

  int nbin = 480./rbin ;

  // Get fake data information
  THStack* ttstk = new THStack("ttstk", "Combined Fitting"); 
  TH1D* fakedata = new TH1D("fakedata","", nbin, 0., 480.);
  TH1D* dth0 = new TH1D("dth0","", nbin, 0., 480.);
  TH1D* dth1 = new TH1D("dth1","", nbin, 0., 480.);
  TH1D* dth2 = new TH1D("dth2","", nbin, 0., 480.);
  TH1D* dth3 = new TH1D("dth3","", nbin, 0., 480.);
  TH1D* dth4 = new TH1D("dth4","", nbin, 0., 480.);
  TH1D* dth123 = new TH1D("dt123","", nbin, 0., 480.);

  getFakeData( rbin, fakedata, ttstk, dth0, dth1, dth2, dth3, dth4 );

  double bestMass = Chi2Test(fakedata, lowBound, upBound ,9 , NBTag);
  // Fit the data
  TF1 *func1 = new TF1("func1",MassFitFunction::fitData, lowBound, upBound,9);

  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 9, NBTag );

  func1->SetParLimits(0,  5., 200.);
  func1->FixParameter(1, fpar[1] );
  func1->FixParameter(2, fpar[2] );
  func1->SetParLimits(3, fpar[3]- 0.1*fpar[3], fpar[3]+0.1*fpar[3] );
  func1->FixParameter(4, fpar[4] );
  func1->FixParameter(5, fpar[5] );
  func1->FixParameter(6, fpar[8] );
  func1->FixParameter(7, fpar[9] );
  func1->SetParLimits(8, fpar[11]- 0.1*fpar[11], fpar[11]+0.1*fpar[11] );

  c3->cd(1);
  ttstk->Draw();

  func1->SetLineColor(4);
  func1->SetLineStyle(2);
  func1->SetLineWidth(2);
  fakedata->Fit( func1, "R","sames",lowBound,upBound);
  func1->Draw("same");
  c3->Update();

  // Draw the expected signal
  Double_t apars[9];
  fprintf(logfile," ---------- Fitting Parameters --------- \n" );
  for (int i=0; i< 9; i++) { 
      apars[i]   = func1->GetParameter(i);
      fprintf(logfile," p%d = %f \n", i, apars[i] );
  }
  fprintf(logfile," ========== Fitting Finished =========== \n" );

  TF1 *func2 = new TF1("func2",MassFitFunction::fitSG, 100, 360, 6 );
  for (int j=0; j<6; j++ ) {
      func2->FixParameter( j, apars[j] );
  }
  func2->SetLineColor(1);
  func2->SetLineStyle(5);
  func2->SetLineWidth(3);

  c3->cd(2);
  dth0->Draw();
  func2->Draw("sames");
  c3->Update();

  // Draw the expected background
 
  TF1 *func3 = new TF1("func3",MassFitFunction::fitLD, 100, 360, 3);
  func3->FixParameter(0, apars[8]*apars[0] );
  func3->FixParameter(1, apars[6] );
  func3->FixParameter(2, apars[7] );
  func3->SetLineColor(1);
  func3->SetLineStyle(2);
  func3->SetLineWidth(3);

  c3->cd(3);
  dth123->Add(dth1, 1);
  dth123->Add(dth2, 1);
  dth123->Add(dth3, 1);
  dth123->Add(dth4, 1);
  dth123->SetFillColor(kOrange+7);
  dth123->Draw();
  func3->Draw("sames");
  c3->Update();

  c3->Print(plot2);

  delete func1;
  delete func2;
  delete func3;
  delete ttstk;
  delete dth0;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth4;
  delete dth123;
  delete fakedata;
 
  fclose(logfile);
}


// chi2 test + minimum chi2 fit, output best estimated mass
double TemplateMassFit::Chi2Test( TH1D* theData, int lowBound, int upBound, int nPar, int NBTag ) {  

   gStyle->SetOptFit(11);
   c9 = new TCanvas("c9","", 800, 600);
   c9->SetGrid();
   c9->SetFillColor(10);
   c9->SetFillColor(10);
   c9->cd();

   Double_t mAssumption[16] = {0.}; 
   Double_t nChi2[16] ={0.}; 
   Double_t fpar[12];
   double maxChi2 = 0. ;
   TF1 *func0 = new TF1("func0", MassFitFunction::fitData1, lowBound, upBound, 12);
   TF1 *fS    = new TF1("fS"   , MassFitFunction::fitSG, lowBound, upBound, 6);
   TF1 *fW    = new TF1("fW", MassFitFunction::fitLD, lowBound, upBound, 3);
   TF1 *fB    = new TF1("fB", MassFitFunction::fitLD, lowBound, upBound, 3);
   for (int i =0; i < 16; i++ ) {

       // set up the mass assumption
       mAssumption[i] = 156. + (i*2.) ;
       SetFitParameters( mAssumption[i], fpar, nPar, NBTag );

       // Fit the data
       func0->SetParLimits(0,  25., 100.);
       func0->FixParameter(1, fpar[1] );
       func0->FixParameter(2, fpar[2] );
       func0->SetParLimits(3, fpar[3]-0.1*fpar[3], fpar[3]+0.1*fpar[3] );
       func0->FixParameter(4, fpar[4] );
       func0->FixParameter(5, fpar[5] );
       if ( nPar == 12 ) {
          func0->FixParameter(6, fpar[6] );
          func0->FixParameter(7, fpar[7] );
          func0->FixParameter(8, fpar[8] );
          func0->FixParameter(9, fpar[9] );
          func0->SetParLimits(10,fpar[10]-0.1*fpar[10], fpar[10]+0.1*fpar[10] );
          func0->FixParameter(11,fpar[11]);
       }
       if ( nPar == 9 ) {
          func0->FixParameter(6, fpar[8] );
          func0->FixParameter(7, fpar[9] );
          func0->SetParLimits(8, fpar[11]-0.1*fpar[11], fpar[11]+0.1*fpar[11] );
       }

       theData->Fit( func0, "RQ0","", lowBound, upBound);

       nChi2[i] = getChi2(theData, func0, fS, fB, fW, 130, 210, nPar );
       //nChi2[i] = getChi2(theData, func0, 130,210);
       //nChi2[i] =  func0->GetChisquare() ;
       if ( nChi2[i] > maxChi2 ) maxChi2 = nChi2[i] ;
   }

   double mErr[16]={0.};
   double x2Err[16]={0.};
   TGraph* x2test = new TGraph(16, mAssumption, nChi2);
   /*
   x2test->SetTitle(" chi2 test ");
   x2test->SetMarkerColor(4);
   x2test->SetMarkerStyle(21);
   x2test->SetMaximum(maxChi2*1.25);
   x2test->SetMinimum(0.0);
   x2test->GetXaxis()->SetTitle(" input mass assumption ");
   x2test->GetYaxis()->SetTitle(" chi2  ");
   x2test->Draw("AP");
   */
   TF1 *func5 = new TF1("func5",MassFitFunction::fitParabola, 154, 188, 3 );
   func5->SetParLimits(0, 150, 190);
   func5->SetParameter(1, 3.5 );
   x2test->Fit( func5, "NR0", "",154,188 );

   std::vector<bool> rejList = fitFunc->DataRejection( func5, mAssumption, nChi2, 16);
   
   for (int i=0; i<16; i++) {
       if ( !rejList[i] ) x2Err[i] = 0.1 ;
       if ( rejList[i] ) x2Err[i] = 10. ;
   }
   TGraphErrors* x2testErr = new TGraphErrors(16, mAssumption, nChi2, mErr, x2Err );
   x2testErr->SetTitle(" chi2 test ");
   x2testErr->SetMarkerColor(4);
   x2testErr->SetMarkerStyle(21);
   x2testErr->SetMaximum(maxChi2*1.25);
   x2testErr->SetMinimum(0.0);
   x2testErr->GetXaxis()->SetTitle(" input mass assumption ");
   x2testErr->GetYaxis()->SetTitle(" chi2  ");
   x2testErr->Draw("AP");
   x2testErr->Fit( func5, "R", "sames",154,188 );
   
   c9->Update();
   c9->Print(plot9);

   double massCandidate =  func5->GetParameter(0);
   gStyle->SetOptFit(111);

   delete fS;
   delete fW;
   delete fB;
   delete func0;
   delete func5;
   delete x2test; 
   delete x2testErr;

   return massCandidate ;
}
 
double TemplateMassFit::getChi2( TH1D* theData1,  TF1* theFunc, TF1* fS, TF1* fB, TF1* fW, double lowBound, double upBound, int nPar ) {

  for (int s=0; s<6; s++) {
      fS->FixParameter( s, theFunc->GetParameter(s) );
  }

  fW->FixParameter( 0, theFunc->GetParameter(11) );
  fW->FixParameter( 1, theFunc->GetParameter(8) );
  fW->FixParameter( 2, theFunc->GetParameter(9) );

  fB->FixParameter( 0, theFunc->GetParameter(10) );
  fB->FixParameter( 1, theFunc->GetParameter(6) );
  fB->FixParameter( 2, theFunc->GetParameter(7) );

  int nbin = theData1->GetNbinsX() ;
  double bW = 480 / nbin ;
  Int_t b1 = lowBound / bW ;
  Int_t b2 =  upBound / bW ;
  double chi2 = 0;
  for (int i= b1 ; i< b2+1; i++) {
      double k = (i-1)*bW + bW/2.  ;
      double theF = theFunc->Eval(k);
      double theW = fW->Eval(k);
      double theB = fB->Eval(k);
      //double theS = fS->Eval(k);
      double theH = theData1->GetBinContent(i);
      //double x2 = ( theH - theF )*( theH - theF ) / (theF - theW) ; 
      //double x2 = ( theH - theF )*( theH - theF ) / theF ; 
      //double x2 = ( theH - theF )*( theH - theF ) / theS ; 
      double x2 = ( theH - theF )*( theH - theF ) / ( theF + theB + theW ); 
      if ( nPar == 9 ) x2 = ( theH - theF )*( theH - theF ) / ( theF + theW ); 
      chi2 += x2;
  }
  //chi2 = chi2 / ( b2+1 - b1 - 3 ) ;
  
  return chi2;
}

void TemplateMassFit::MultiTemplatesFitting( int rbin, int lowBound, int upBound ){

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatTextColor(1);

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  // Get fake data information
  THStack* ttstk = new THStack("ttstk", "Template Fitting"); 
  TH1D* dt   = new TH1D("dt","", nbin, 0., 480.);
  TH1D* dth0 = new TH1D("dth0","", nbin, 0., 480.);
  TH1D* dth1 = new TH1D("dth1","", nbin, 0., 480.);
  TH1D* dth2 = new TH1D("dth2","", nbin, 0., 480.);
  TH1D* dth3 = new TH1D("dth3","", nbin, 0., 480.);
  TH1D* dth23 = new TH1D("dth23","", nbin, 0., 480.);
  getFakeData( rbin, dt, ttstk, dth0, dth1, dth2, dth3 );

  // pre-fitting
  Double_t sPred[sz] ;
  Double_t tPred[sz] ;
  Double_t bPred[sz] ;
  double chi2[5] = {99999.};
  chi2[0] = TemplateTest("161", dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  chi2[1] = TemplateTest("166", dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  chi2[2] = TemplateTest("171", dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  chi2[3] = TemplateTest("176", dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  chi2[4] = TemplateTest("181", dt, rbin, lowBound, upBound, sPred, tPred, bPred );

  double minChi2 = 99999.;
  int bestMass = -1 ;
  for (int s1 = 0; s1 < 5; s1++) {
      cout<<" chi2("<<s1<<") : "<< chi2[s1] << endl;
      if ( chi2[s1] < minChi2 ) {
         minChi2 = chi2[s1] ;
         bestMass = s1 ;
      }
  }
  cout<<"  => "<< bestMass <<" win "<<endl;

  double fChi2 = -1; 
  if ( bestMass == 0 ) fChi2 = TemplateTest( "161" , dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  if ( bestMass == 1 ) fChi2 = TemplateTest( "166" , dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  if ( bestMass == 2 ) fChi2 = TemplateTest( "171" , dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  if ( bestMass == 3 ) fChi2 = TemplateTest( "176" , dt, rbin, lowBound, upBound, sPred, tPred, bPred );
  if ( bestMass == 4 ) fChi2 = TemplateTest( "181" , dt, rbin, lowBound, upBound, sPred, tPred, bPred );


    c6 = new TCanvas("c5","", 900, 800);
    c6->SetFillColor(10);
    c6->SetFillColor(10);
    c6->Divide(2,2);
    
    TH1D* hPred  = new TH1D("hPred" ,"", nbin, 0., 480.);
    TH1D* hPred0 = new TH1D("hPred0","", nbin, 0., 480.);
    TH1D* hPred1 = new TH1D("hPred1","", nbin, 0., 480.);
    TH1D* hPred2 = new TH1D("hPred2","", nbin, 0., 480.);
    for (int  j= b1; j<= b2 ; j++) {
        hPred0->SetBinContent(j, sPred[j-b1] );
        hPred1->SetBinContent(j, tPred[j-b1] );
        hPred2->SetBinContent(j, bPred[j-b1] );
        hPred->SetBinContent(j,  sPred[j-b1] + tPred[j-b1] + bPred[j-b1] );
    }

    // overall fitting
    c6->cd(1);
    ttstk->SetTitle(" Template Fitting Result ");
    ttstk->Draw();

    hPred->SetMarkerStyle(21);
    hPred->SetMarkerColor(1);
    hPred->SetMarkerSize(0.8);
    hPred->Draw("samesPE");

    c6->Update();
    
    // signal component
    c6->cd(2);
    dth0->SetFillColor(3);
    dth0->Draw();
    hPred0->SetMarkerStyle(21);
    hPred0->SetMarkerColor(1);
    hPred0->SetMarkerSize(0.8);
    hPred0->Draw("samesPE");

    TF1* func0 = new TF1("func0", MassFitFunction::fitSG , 100, 350, 6);
    func0->SetParLimits(0, 40.,120.);
    func0->SetParLimits(1, 140.,230.);
    func0->SetParLimits(2, 10.,100.);
    func0->SetParLimits(3, 10., 100.);
    func0->SetParLimits(4, 5., 10.);
    func0->SetParLimits(5, 0.3, 10.);
    func0->SetLineStyle(2);
    func0->SetLineWidth(2);
    hPred0->Fit( func0, "R","sames", 100, 350 );

    c6->Update();

    c6->cd(3);
    dth1->SetFillColor(7);
    dth1->Draw();
    hPred1->SetMarkerStyle(21);
    hPred1->SetMarkerColor(1);
    hPred1->SetMarkerSize(0.8);
    hPred1->Draw("samesEP");
    TF1* func1 = new TF1("func1", MassFitFunction::fitLD , 100, 380, 3);
    func1->SetParLimits(0, 100.,2000.);
    func1->SetParLimits(1, 50.,500.);
    func1->SetParLimits(2,  1., 100.);
    func1->SetLineStyle(2);
    func1->SetLineWidth(2);
    hPred1->Fit( func1, "R","sames", 100, 380 );

    c6->cd(4);
    dth23->Add(dth2, 1);
    dth23->Add(dth3, 1);
    dth23->SetFillColor(6);
    dth23->Draw();
    hPred2->SetMarkerStyle(21);
    hPred2->SetMarkerColor(1);
    hPred2->SetMarkerSize(0.8);
    hPred2->Draw("samesEP");
    TF1* func2 = new TF1("func2", MassFitFunction::fitLD , 100, 380, 3);
    func2->SetParLimits(0, 100.,2000.);
    func2->SetParLimits(1, 50.,500.);
    func2->SetParLimits(2,  1., 100.);
    func2->SetLineStyle(2);
    func2->SetLineWidth(2);
    hPred2->Fit( func2, "R","sames", 100, 380 );
    c6->Update();

    c6->Print(plot6);

  delete hPred;
  delete hPred0;
  delete hPred1;
  delete hPred2;
  delete dt;
  delete ttstk;
  delete dth0;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth23;
  delete func0;
  delete func1;
  delete func2;

}

void TemplateMassFit::TemplateFitting( int rbin, int lowBound, int upBound ){

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatTextColor(1);

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  // Get fake data information
  THStack* ttstk = new THStack("ttstk", "Template Fitting"); 
  TH1D* dt   = new TH1D("dt","", nbin, 0., 480.);
  TH1D* dth0 = new TH1D("dth0","", nbin, 0., 480.);
  TH1D* dth1 = new TH1D("dth1","", nbin, 0., 480.);
  TH1D* dth2 = new TH1D("dth2","", nbin, 0., 480.);
  TH1D* dth3 = new TH1D("dth3","", nbin, 0., 480.);
  TH1D* dth123 = new TH1D("dth123","", nbin, 0., 480.);
  getFakeData( rbin, dt, ttstk, dth0, dth1, dth2, dth3 );

  // pre-fitting
  Double_t sPred[sz] ;
  Double_t bPred[sz] ;
  double chi2[5] = {99999.};
  chi2[0] = TemplateTest("161", dt, rbin, lowBound, upBound, sPred, bPred );
  chi2[1] = TemplateTest("166", dt, rbin, lowBound, upBound, sPred, bPred );
  chi2[2] = TemplateTest("171", dt, rbin, lowBound, upBound, sPred, bPred );
  chi2[3] = TemplateTest("176", dt, rbin, lowBound, upBound, sPred, bPred );
  chi2[4] = TemplateTest("181", dt, rbin, lowBound, upBound, sPred, bPred );

  double minChi2 = 99999.;
  int bestMass = -1 ;
  for (int s1 = 0; s1 < 5; s1++) {
      cout<<" chi2("<<s1<<") : "<< chi2[s1] << endl;
      if ( chi2[s1] < minChi2 ) {
         minChi2 = chi2[s1] ;
         bestMass = s1 ;
      }
  }
  cout<<"  => "<< bestMass <<" win "<<endl;

  double fChi2 = -1; 
  if ( bestMass == 0 ) fChi2 = TemplateTest( "161" , dt, rbin, lowBound, upBound, sPred, bPred );
  if ( bestMass == 1 ) fChi2 = TemplateTest( "166" , dt, rbin, lowBound, upBound, sPred, bPred );
  if ( bestMass == 2 ) fChi2 = TemplateTest( "171" , dt, rbin, lowBound, upBound, sPred, bPred );
  if ( bestMass == 3 ) fChi2 = TemplateTest( "176" , dt, rbin, lowBound, upBound, sPred, bPred );
  if ( bestMass == 4 ) fChi2 = TemplateTest( "181" , dt, rbin, lowBound, upBound, sPred, bPred );


    c4 = new TCanvas("c4","", 900, 700);
    c4->SetFillColor(10);
    c4->SetFillColor(10);
    c4->Divide(2,2);
    
    TH1D* hPred  = new TH1D("hPred" ,"", nbin, 0., 480.);
    TH1D* hPred0 = new TH1D("hPred0","", nbin, 0., 480.);
    TH1D* hPred1 = new TH1D("hPred1","", nbin, 0., 480.);
    for (int  j= b1; j<= b2 ; j++) {
        hPred0->SetBinContent(j, sPred[j-b1] );
        hPred1->SetBinContent(j, bPred[j-b1] );
        hPred->SetBinContent(j,  sPred[j-b1]+bPred[j-b1] );
    }

    // overall fitting
    c4->cd(1);
    ttstk->SetTitle(" Template Fitting Result ");
    ttstk->Draw();

    hPred->SetMarkerStyle(21);
    hPred->SetMarkerColor(1);
    hPred->SetMarkerSize(0.8);
    hPred->Draw("samesPE");

    c4->Update();
    
    // signal component
    c4->cd(2);
    dth0->SetFillColor(3);
    dth0->Draw();
    hPred0->SetMarkerStyle(21);
    hPred0->SetMarkerColor(1);
    hPred0->SetMarkerSize(0.8);
    hPred0->Draw("samesPE");

    TF1* func0 = new TF1("func0", MassFitFunction::fitSG , 100, 350, 6);
    func0->SetParLimits(0, 40.,120.);
    func0->SetParLimits(1, 140.,230.);
    func0->SetParLimits(2, 10.,100.);
    func0->SetParLimits(3, 10., 100.);
    func0->SetParLimits(4, 5., 10.);
    func0->SetParLimits(5, 0.3, 10.);
    func0->SetLineStyle(2);
    func0->SetLineWidth(2);
    hPred0->Fit( func0, "R","sames", 100, 350 );

    c4->Update();

    c4->cd(3);
    dth123->Add(dth1, 1);
    dth123->Add(dth2, 1);
    dth123->Add(dth3, 1);
    dth123->SetFillColor(6);
    dth123->Draw();
    hPred1->SetMarkerStyle(21);
    hPred1->SetMarkerColor(1);
    hPred1->SetMarkerSize(0.8);
    hPred1->Draw("samesEP");
    TF1* func1 = new TF1("func1", MassFitFunction::fitLD , 100, 380, 3);
    func1->SetParLimits(0, 100.,2000.);
    func1->SetParLimits(1, 50.,500.);
    func1->SetParLimits(2,  1., 100.);
    func1->SetLineStyle(2);
    func1->SetLineWidth(2);
    hPred1->Fit( func1, "R","sames", 100, 380 );
    c4->Update();

    c4->Print(plot4);

  delete hPred;
  delete hPred0;
  delete hPred1;
  delete dt;
  delete ttstk;
  delete dth0;
  delete dth1;
  delete dth2;
  delete dth3;
  delete dth123;

}
// Multiple Templates
double TemplateMassFit::TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *tPred, Double_t *bPred ) {  

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
 
  // tmp1: tt wrong combinatorics
  TH1D* tt = new TH1D("tt","", nbin, 0., 480.);
  getBackground( tt, 1, rbin, fname );

  // tmp2: wjets and QCD 
  TH1D* wj = new TH1D("wj","", nbin, 0., 480.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 480.);
  TH1D* bg = new TH1D("bg","", nbin, 0., 480.);
  getBackground( wj, 2, rbin, fname );
  getBackground( qm, 3, rbin, fname );
  bg->Add(wj,qm, 1., 0.);

  // get template signal distribution
  TH1D* sg = new TH1D("sg","", nbin, 0., 480.);
  getSignal( sg, rbin , mName);

  // normalize signal and background shape
  ScaleTemplates( 5000., sg, b1, b2 );
  ScaleTemplates( 5000., tt, b1, b2 );
  ScaleTemplates( 5000., bg, b1, b2 );

  Double_t tmpar[6] ;
  TH1D* parasg = new TH1D("parasg","", nbin, 0., 480.);
  SmoothTemplate( 1 , sg, parasg, lowBound, upBound, tmpar );
  TH1D* paratt = new TH1D("paratt","", nbin, 0., 480.);
  SmoothTemplate( 0 , tt, paratt, lowBound, upBound, tmpar );
  TH1D* parabg = new TH1D("parabg","", nbin, 0., 480.);
  SmoothTemplate( 0 , bg, parabg, lowBound, upBound, tmpar );

  // start template fitting
  TObjArray *mctmp = new TObjArray(3);
  mctmp->Add(parasg);
  mctmp->Add(paratt);
  mctmp->Add(parabg);
 
  TFractionFitter* mfit = new TFractionFitter( theData, mctmp );
  mfit->Constrain( 0, 0.0, 0.15);
  mfit->Constrain( 1, 0.6, 0.8);
  mfit->Constrain( 2, 0.2, 0.4);
  mfit->SetRangeX(b1,b2);
  Int_t status = mfit->Fit();
  cout<<" fit status : "<< status << endl;
 
  double x2 = 99999.; 

  if ( status == 0 ) {

     //TH1D* result  = (TH1D*) mfit->GetPlot();
     double ratio[3];
     double err[3];
     mfit->GetResult(0,ratio[0],err[0]);
     mfit->GetResult(1,ratio[1],err[1]);
     mfit->GetResult(2,ratio[2],err[2]);

     double nData = theData->Integral(b1, b2);
     double nSig  = sg->Integral(b1, b2);
     double nBkgd = bg->Integral(b1, b2);
     double nBgtt = tt->Integral(b1, b2);
     cout<<"   nData= "<<nData<<"  nSignal= "<<nSig<<"  nBgtt="<< nBgtt <<"  nBkgd="<< nBkgd << endl;
     cout<<"   r0:"<<ratio[0] <<"  r1:"<<ratio[1]<<"  r2:"<<ratio[2]<<endl;

     for(int k=0; k<sz ; k++) {
        int i = b1 + k ;
        double sN = sg->GetBinContent(i); 
        double bN = bg->GetBinContent(i); 
        double tN = tt->GetBinContent(i); 
        sPred[k] = (nData/nSig)*ratio[0]*sN ;
        tPred[k] = (nData/nBgtt)*ratio[1]*tN ;
        bPred[k] = (nData/nBkgd)*ratio[2]*bN ;
     }
     x2 = mfit->GetChisquare() / mfit->GetNDF() ;
 
     //delete result ;
  }

  delete bg;
  delete sg;
  delete tt;
  delete qm;
  delete wj;
  delete parabg;
  delete parasg;
  delete paratt;
  delete mctmp;
  delete mfit;
  
  return x2;

}

// Single background Template
double TemplateMassFit::TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *bPred ) {  

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
 
  // combine all background channel for background template
  TH1D* bg = new TH1D("bg","", nbin, 0., 480.);
  combineBG( mName, bg, rbin );
  // get template signal distribution
  TH1D* sg = new TH1D("sg","", nbin, 0., 480.);
  getSignal( sg, rbin , mName);

  // normalize signal and background shape
  ScaleTemplates( 5000., sg, b1, b2 );
  ScaleTemplates( 5000., bg, b1, b2 );

  Double_t tmpar[6] ;
  TH1D* parasg = new TH1D("parasg","", nbin, 0., 480.);
  SmoothTemplate( 1 , sg, parasg, lowBound, upBound, tmpar );
  TH1D* parabg = new TH1D("parabg","", nbin, 0., 480.);
  SmoothTemplate( 0 , bg, parabg, lowBound, upBound, tmpar );

  // start template fitting
  TObjArray *mctmp = new TObjArray(2);
  mctmp->Add(parasg);
  mctmp->Add(parabg);
 
  TFractionFitter* mfit = new TFractionFitter( theData, mctmp );
  mfit->Constrain( 0, 0.0, 0.3);
  mfit->Constrain( 1, 0.7, 1.0);
  mfit->SetRangeX(b1,b2);
  Int_t status = mfit->Fit();
  cout<<" fit status : "<< status << endl;
 
  double x2 = 99999.; 

  if ( status == 0 ) {

     //TH1D* result  = (TH1D*) mfit->GetPlot();
     double ratio[2];
     double err[2];
     mfit->GetResult(0,ratio[0],err[0]);
     mfit->GetResult(1,ratio[1],err[1]);

     double nData = theData->Integral(b1, b2);
     double nSig  = sg->Integral(b1, b2);
     double nBkgd = bg->Integral(b1, b2);
     cout<<"   nData= "<<nData<<"  nSignal= "<<nSig<<"  nBgd="<<nBkgd<<endl;
     cout<<"   r0:"<<ratio[0] <<"  r1:"<<ratio[1]<<endl;

     for(int k=0; k<sz ; k++) {
        int i = b1 + k ;
        double sN = sg->GetBinContent(i); 
        double bN = bg->GetBinContent(i); 
        sPred[k] = (nData/nSig)*ratio[0]*sN ;
        bPred[k] = (nData/nBkgd)*ratio[1]*bN ;
     }
     x2 = mfit->GetChisquare() / mfit->GetNDF() ;
 
     //delete result ;
  }

  delete bg;
  delete sg;
  delete parabg;
  delete parasg;
  delete mctmp;
  delete mfit;
  
  return x2;

}

void TemplateMassFit::FitSignal( TString mName, int rbin, int lowBound, int upBound, Double_t* para, Double_t* perr ){

  FILE* testfile = fopen(hfolder+"/Sgpara.log","a"); 

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  plot7 = hfolder+"/FitTest_"+mName+"_"+cname+".gif";
  plot3 = hfolder+"/MFit_"+mName+"_"+cname+".gif";

  c7 = new TCanvas("c7","", 1000, 800);
  c7->SetFillColor(10);
  c7->SetFillColor(10);
  c7->Divide(2,2);

  double m1 = MassDigi(mName);
  lowBound = m1 - 50 ;
  upBound  = m1 + 170 ;

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  // set up the y-axis max
  double yMax = 120. *( rbin / 20. ) ;
  //double yMax = 220. *( rbin / 20. ) ;

  c7->cd(1);

  // tmp0: get template signal distribution
  TH1D* sg = new TH1D("sg","", nbin, 0., 480.);
  getSignal( sg, rbin , mName );

  sg->SetFillColor(7);
  //int maxbin = sg->GetMaximumBin() ;
  //double maxval = sg->GetBinContent(maxbin) ;
  //sg->SetMaximum( maxval*1.2 );
  sg->SetMaximum( yMax );

  sg->Draw();

  TF1* func0 = new TF1("func0", MassFitFunction::fitLG , 80, 450, 3);
  TF1* func1 = new TF1("func1", MassFitFunction::fitGS , lowBound, upBound, 3);
  TF1* func2 = new TF1("func2", MassFitFunction::fitSG , lowBound, upBound, 6);

  // pre-set the value
  double p0 = 50. ;
  double p1 = m1  ;
  double p2 = 20. ;
  double p3 = 33. ; 
  double p4 = log(m1) + sqrt(1./20.) ;
  double p5 =  5. ;

  // 1st Fit, Fix "mean" value for gaussisan & log-normal and allow normalization and width vary
  func2->FixParameter( 1, p1 );
  func2->SetParLimits(2, p2-0.3*p2, p2+1.0*p2);
  func2->SetParLimits(3, p3-0.1*p3, p3+0.1*p3);
  func2->FixParameter(4, p4 );
  if ( cname == "Had" ) func2->FixParameter(5, p5 );
  if ( cname == "Lep" ) func2->FixParameter(5, 3. );
  
  sg->Fit( func2, "RQ0","", lowBound, upBound );

  p0 = func2->GetParameter(0); 
  p1 = func2->GetParameter(1); 
  p2 = func2->GetParameter(2); 
  p3 = func2->GetParameter(3); 
  p5 = func2->GetParameter(5); 

  // Draw the gaussian
  func1->FixParameter(0, p0 );
  func1->FixParameter(1, p1 );
  func1->FixParameter(2, p2 );
  func1->Draw("sames");

  c7->Update();

  c7->cd(2);

  sg->SetFillColor(7);
  sg->SetMaximum( yMax );
  sg->Draw();

  // 2nd Fit , Allow gaussian change
  func2->SetParLimits(0, p0-0.1*p0, p0+0.1*p0);
  func2->SetParLimits(1, m1 - 1., m1 + 1. );
  func2->SetParLimits(2, p2-0.1*p2, p2+0.1*p2 );
  func2->SetParLimits(3, p3-0.05*p3,p3+0.05*p3);
  func2->FixParameter(4, p4 );
  func2->FixParameter(5, p5 );
  sg->Fit( func2, "R","sames", lowBound, upBound); 

  p0 = func2->GetParameter(0); 
  p1 = func2->GetParameter(1); 
  p2 = func2->GetParameter(2); 
  p3 = func2->GetParameter(3); 
  p4 = func2->GetParameter(4); 
  p5 = func2->GetParameter(5); 
 
  func0->FixParameter(0, p0*p3);
  func0->FixParameter(1, p4   );
  func0->FixParameter(2, p5   );

  func0->SetLineColor(2);
  func0->SetLineWidth(3);
  func0->SetLineStyle(2);
  func0->Draw("sames");

  c7->Update();

  c7->cd(3);
  // 3rd Fit , Final tunning
  func2->SetParLimits(0, p0- 0.1*p0 , p0+0.1*p0);
  func2->SetParLimits(1, m1 - 1.   , m1 + 1. );
  func2->SetParLimits(2, p2- 0.1*p2, p2+ 0.1*p2);
  func2->SetParLimits(3, p3-0.01*p3, p3+0.01*p3);
  func2->FixParameter(4, p4 );
  func2->SetParLimits(5, p5-0.01*p5, p5+0.01*p5 );
  sg->SetFillColor(7);
  sg->SetMaximum( yMax );
  sg->Draw();
  sg->Fit( func2, "R","sames", lowBound, upBound); 

  // print out covariant matrics
  /*
  double mtx[5][5];
  gMinuit->mnemat(&mtx[0][0],5);
  for (int j=0; j < 5 ; j++) {
      for (int i=0; i<5; i++) {
          fprintf(testfile,"  %.3f",  mtx[i][j] );
          if ( i == 4 ) fprintf(testfile," \n" );
      }
  }
  fprintf(testfile," ------ \n" );
  */

  p0 = func2->GetParameter(0); 
  p1 = func2->GetParameter(1); 
  p2 = func2->GetParameter(2); 
  p3 = func2->GetParameter(3); 
  p4 = func2->GetParameter(4); 
  p5 = func2->GetParameter(5); 

  fprintf(testfile," %.1f", m1 );
  for (int i=0; i<6; i++) {  
      fprintf(testfile,"  %.3f  %.3f",  func2->GetParameter(i), func2->GetParError(i) );
      if ( para != NULL ) para[i] = func2->GetParameter(i);
      if ( perr != NULL ) perr[i] = func2->GetParError(i) ;
  }
  fprintf(testfile," \n" );

  func0->FixParameter(0, p0*p3);
  func0->FixParameter(1, p4   );
  func0->FixParameter(2, p5   );

  func0->SetLineColor(2);
  func0->SetLineWidth(3);
  func0->SetLineStyle(2);
  func0->Draw("sames");

  c7->Update();
  c7->Print(plot7);
 
  c3 = new TCanvas("c3","", 800, 600);
  c3->SetFillColor(10);
  c3->SetFillColor(10);

  sg->Draw();
  func2->Draw("sames");
  func0->Draw("sames");
  c3->Update();
  c3->Print(plot3);

  delete func2;
  delete func1;
  delete func0;
  delete sg;

  fclose(testfile);
}

void TemplateMassFit::FitBackground( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp, Double_t *para, Double_t *perr ){

  FILE* bgpara = fopen(hfolder+"/Bgpara.log","a"); 

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  int nbin = 480./rbin ;

  double m1 = MassDigi(mName);

  plot8 = hfolder+"/"+channel+"_BG_"+mName+".gif";
  // tmp1: tt wrong combinatorics
  TH1D* tt = new TH1D("tt","", nbin, 0., 480.);
  // tmp2: wjets, single Top t, single Top tW, QCD
  TH1D* wj  = new TH1D("wj" ,"", nbin, 0., 480.);
  TH1D* stt = new TH1D("stt","", nbin, 0., 480.);
  TH1D* stw = new TH1D("stw","", nbin, 0., 480.);
  getBackground( tt, 1, rbin, mName );
  getBackground( wj, 2, rbin, mName );
  getBackground( stt, 3, rbin, mName );
  getBackground( stw, 4, rbin, mName );
  // mix the background
  TH1D* bg = new TH1D("bg","", nbin, 0., 480.);
  if ( comp[1] ) bg->Add(tt,1.);
  if ( comp[2] ) bg->Add(wj,1.);
  if ( comp[3] ) bg->Add(stt,1.);
  if ( comp[4] ) bg->Add(stw,1.);
  c8 = new TCanvas("c8","", 800, 600);
  c8->SetFillColor(10);
  c8->SetFillColor(10);
  c8->Divide(1,2);

  c8->cd(1);

  // Draw the major background
  tt->SetLineWidth(2);
  tt->SetLineColor(4);
  tt->Draw();

  // Fit tt-wrong permutation
  TF1* func0 = new TF1("func0", MassFitFunction::fitLD , 90, 380, 3);
  func0->SetParLimits(0, 100.,3000.);
  func0->SetParLimits(1, m1-10, m1+20.);
  func0->SetParLimits(2,  1., 100.);
  func0->SetLineColor(4);
  func0->SetLineStyle(2);
  tt->Fit( func0, "R","sames", lowBound, upBound );
  c8->Update();

  c8->cd(2);
  bg->SetLineWidth(2);
  bg->SetLineColor(1);
  bg->Draw();
  stt->SetLineColor(4);
  stt->SetLineWidth(3);
  stt->Draw("sames");
  stw->SetLineColor(6);
  stw->SetLineWidth(3);
  stw->Draw("sames");
  wj->SetLineColor(2);
  wj->SetLineWidth(3);
  wj->Draw("sames");
  // Fit the combined backgrounds
  TF1* func1 = new TF1("func1", MassFitFunction::fitLD , 90, 380, 3);
  func1->SetParLimits(0, 100.,3000.);
  func1->SetParLimits(1, m1-10, m1+20.);
  func1->SetParLimits(2,  1., 100.);
  func1->SetLineColor(1);
  func1->SetLineStyle(2);
  bg->Fit( func1, "R","sames", lowBound, upBound );

  // para[0~2] : tt wrong permutation 
  // para[3~5] : combined backgrounds 
  fprintf(bgpara," %.1f", m1 );
  for (int i=0; i<6; i++) {  
      if ( i < 3 ) {
         fprintf(bgpara,"  %.3f  %.3f",  func0->GetParameter(i), func0->GetParError(i) );
         if ( para != NULL ) para[i] = func0->GetParameter(i);
         if ( perr != NULL ) perr[i] = func0->GetParError(i) ;
      } else {
         fprintf(bgpara,"  %.3f  %.3f",  func1->GetParameter(i-3), func1->GetParError(i-3) );
         if ( para != NULL ) para[i] = func1->GetParameter(i-3);
         if ( perr != NULL ) perr[i] = func1->GetParError(i-3) ;
      }
  }
  fprintf(bgpara," \n" );

  c8->Update();
  c8->Print(plot8);
 
  delete bg;
  delete tt;
  delete wj;
  delete stt;
  delete stw;
  delete func0;
  delete func1;
  delete c8;

  fclose(bgpara);
}

void TemplateMassFit::GetAllCoeff( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp ) {

     Double_t sPar[6];
     Double_t sErr[6];
     Double_t bPar[6];
     Double_t bErr[6];
     FitSignal( mName, rbin, lowBound, upBound, sPar, sErr );
     FitBackground( mName, rbin, lowBound, upBound, comp, bPar, bErr );

     parfile = fopen(hfolder+"/paraf.log","a");
     errfile = fopen(hfolder+"/perrf.log","a");
     for ( int i =0; i< 6; i++) {
         fprintf(parfile,"  %.2f",  sPar[i] );
         fprintf(errfile,"  %.2f",  sErr[i] );
     }
     for ( int i =0; i< 6; i++) {
         if ( i != 0 && i != 3 ) {
            fprintf(parfile,"  %.2f",  bPar[i] );
            fprintf(errfile,"  %.2f",  bErr[i] );
         }
     }

     double err_bs0 = sqrt( pow(bErr[0]/sPar[0], 2) + pow( (bPar[0]*sErr[0]/(sPar[0]*sPar[0])),2 ) ) ;
     double err_bs3 = sqrt( pow(bErr[3]/sPar[0], 2) + pow( (bPar[3]*sErr[0]/(sPar[0]*sPar[0])),2 ) ) ;
     fprintf(parfile,"  %.2f",  bPar[0]/sPar[0] );
     fprintf(errfile,"  %.2f",  err_bs0 );
     fprintf(parfile,"  %.2f",  bPar[3]/sPar[0] );
     fprintf(errfile,"  %.2f",  err_bs3 );
     fprintf(parfile," \n" );
     fprintf(errfile," \n" );
     
     fclose(parfile);
     fclose(errfile);

}

void TemplateMassFit::TemplateDrawer( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp ){

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  plot8 = hfolder+"/"+channel+"_SBR_"+mName+".gif";

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  // tmp0: get template signal distribution
  TH1D* sg = new TH1D("sg","", nbin, 0., 480.);
  getSignal( sg, rbin , mName );
  // tmp1: tt wrong combinatorics
  TH1D* tt = new TH1D("tt","", nbin, 0., 480.);
  getBackground( tt, 1, rbin, mName );
  // tmp2: wjets and QCD 
  TH1D* wj  = new TH1D( "wj","", nbin, 0., 480.);
  TH1D* stt = new TH1D("stt","", nbin, 0., 480.);
  TH1D* stw = new TH1D("stw","", nbin, 0., 480.);
  TH1D* qm  = new TH1D( "qm","", nbin, 0., 480.);
  TH1D* bg  = new TH1D( "bg","", nbin, 0., 480.);
  getBackground(  wj, 2, rbin, mName );
  getBackground( stt, 3, rbin, mName );
  getBackground( stw, 4, rbin, mName );
  if ( comp[1] ) bg->Add(tt,  1.);
  if ( comp[2] ) bg->Add(wj,  1.);
  if ( comp[3] ) bg->Add(stt, 1.);
  if ( comp[4] ) bg->Add(stw, 1.);
  if ( comp[5] ) bg->Add(qm,  2.01);

  // get the fixed ratio btw wjets and ttbar signal
  Double_t nspar[6] ;
  Double_t ntpar[6] ;
  Double_t nbpar[6] ;
  Double_t nserr[6] ;
  Double_t nterr[6] ;
  Double_t nberr[6] ;
  TH1D* presg = new TH1D("presg","", nbin, 0., 480.);
  SmoothTemplate( 1 , sg, presg, lowBound, upBound, nspar, nserr );
  TH1D* prett = new TH1D("prett","", nbin, 0., 480.);
  SmoothTemplate( 0 , tt, prett, lowBound, upBound, ntpar, nterr );
  TH1D* prebg = new TH1D("prebg","", nbin, 0., 480.);
  SmoothTemplate( 0 , bg, prebg, lowBound, upBound, nbpar, nberr );

  double tWRatio = nbpar[0] / nspar[0] ;
  double tbRatio = ntpar[0] / nspar[0] ;
 
  c8 = new TCanvas("c8","", 1000, 800);
  c8->SetFillColor(10);
  c8->SetFillColor(10);
  c8->Divide(2,2);

  // Draw the Signal - tt-correct permutation
  c8->cd(1);
  sg->SetLineColor(1);
  sg->SetLineWidth(2);
  sg->SetTitle(" Signal");
  sg->Draw();
  c8->Update();

  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  presg->SetLineColor(2);
  presg->SetLineStyle(2);
  presg->SetLineWidth(3);
  presg->Draw("sames");
  c8->Update();
  gStyle->SetStatX(0.95);
  
  // Draw the combined Background - (tt wrong permutation) + (W+jet) + (QCD)
  c8->cd(2);
  gStyle->SetStatTextColor(1);
  bg->SetLineColor(1);
  bg->SetLineWidth(2);
  bg->SetTitle("combined backgrounds");
  bg->Draw(); 
  c8->Update();

  gStyle->SetStatX(0.32);
  gStyle->SetStatTextColor(2);
  prebg->SetLineColor(2);
  prebg->SetLineStyle(2);
  prebg->SetLineWidth(3);
  prebg->Draw("sames");
  c8->Update();
  gStyle->SetStatX(0.95);

  // Draw the background - W+Jets only 
  if ( !comp[1] ) { 
     c8->cd(3);
     gStyle->SetStatTextColor(1);
     tt->SetLineColor(1);
     tt->SetLineWidth(2);
     tt->SetTitle("Tt Wrong Combinatorics");
     tt->Draw(); 
     c8->Update();

     gStyle->SetStatX(0.32);
     gStyle->SetStatTextColor(2);
     prett->SetLineColor(2);
     prett->SetLineStyle(2);
     prett->SetLineWidth(3);
     prett->Draw("sames");
     c8->Update();
     gStyle->SetStatX(0.95);
  }

  gStyle->SetStatTextColor(1);
  c8->Print(plot8);

  // record parameters
  //  p0-p5 : signal  ;  p6-p7 : background(tt-wrong permuation) ; p8-p9 : background(combined background/w+jets) 
  //  p10: bg/sg      ;  p11: tt-wrong/sg
  Double_t tmpars[10] ;
  Double_t tmerrs[10] ;
  for ( int i =0; i< 10; i++) {
      if ( i   < 6  )  tmpars[i] = nspar[i];
      if ( i  == 6  )  tmpars[i] = ntpar[1];
      if ( i  == 7  )  tmpars[i] = ntpar[2];
      if ( i  == 8  )  tmpars[i] = nbpar[1];
      if ( i  == 9  )  tmpars[i] = nbpar[2];

      if ( i   < 6  )  tmerrs[i] = nserr[i];
      if ( i  == 6  )  tmerrs[i] = nterr[1];
      if ( i  == 7  )  tmerrs[i] = nterr[2];
      if ( i  == 8  )  tmerrs[i] = nberr[1];
      if ( i  == 9  )  tmerrs[i] = nberr[2];
  }

  parfile = fopen(hfolder+"/paraf.log","a"); 
  errfile = fopen(hfolder+"/perrf.log","a"); 
  for ( int i =0; i< 10; i++) {
      fprintf(parfile,"  %.2f",  tmpars[i] );
      fprintf(errfile,"  %.2f",  tmerrs[i] );
  }
  fprintf(parfile,"  %.2f",  tbRatio );
  fprintf(parfile,"  %.2f",  tWRatio );
  fprintf(parfile," \n" );
  fprintf(errfile,"  %.2f",  0. );
  fprintf(errfile,"  %.2f",  0. );
  fprintf(errfile," \n" );
  fclose(parfile);
  fclose(errfile);

  delete sg;
  delete bg;
  delete tt;
  delete wj;
  delete qm;
  delete presg;
  delete prebg;
  delete prett;

}



void TemplateMassFit::ScaleTemplates( double factor, TH1D* tmp, int B1, int B2 ){

  double tmpN = tmp->Integral(B1,B2);
  double Scal = factor / tmpN ;
  tmp->Scale(Scal);
  
}

// luminosity of data, nEvents of template
void TemplateMassFit::ScaleTemplates( double lumi, double nEvents, int channel, TH1D* tmp ){

  //**old for 100 /pb => ttbar signal = 8992 events ; wjets = 30393 ; QCD = 2841112
  //**new for 100 /pb => ttbar signal = 9067 events ; wjets = 34000 ; 
  //                     st_tw = 2730, st_t = 6360 
  //                     
  // channel #   => ttbar signal =  1 ; wjets = 2 ; Signle Top t_ch = 3 ; Single Top tW_ch = 4 ; QCD = 5
  double lumiScale = lumi / 100 ;
  double nBase = nEvents;
  if ( channel == 1 ) nBase = 9067 ;
  if ( channel == 2 ) nBase = 34000*1.402 ;
  if ( channel == 3 ) nBase = 6360 ;
  if ( channel == 4 ) nBase = 2730 ;
  if ( channel == 5 ) nBase = 2841112 ;
  double Scal = (nBase*lumiScale) / nEvents ;
  tmp->Scale(Scal);
  
}

void TemplateMassFit::SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, Double_t* pars, Double_t* perr ){

  TF1* func0;
  if (type == 0) {
     func0 = new TF1("func0", MassFitFunction::fitLD , 90, 380, 3);
     func0->SetParLimits(0, 100.,3000.);
     func0->SetParLimits(1, 50.,500.);
     func0->SetParLimits(2,  1., 100.);
     ds->Fit( func0, "RQ0","", 110, 380 );
  }
  if (type == 1) {

     // using my def gaus + logNormal
     func0 = new TF1("func0", MassFitFunction::fitSG , 100, 350, 6);
     
     func0->SetParLimits(0, 20.,1000.);
     func0->SetParLimits(1, 140.,230.);
     func0->SetParLimits(2, 22.,24.);
     func0->SetParLimits(3, 27., 45.);
     func0->FixParameter(4, log(175.) );
     func0->FixParameter(5, 5.);
     
     // convoluted BW-Gaus
     //func0 = new TF1("func0", MassFitFunction::ConvBWGS , 100, 300, 5);
     //func0->SetParameters( bw0, bw1, bw2, bw1, 20. );
     //func0->SetParLimits(0, bw0-20., bw0+20.);
     //func0->SetParLimits(1, bw1-20., bw1+20.);

     // Landau distribution
     ds->Fit( func0, "RQ0","", 100, 350 );
  }

  for (int h=0; h< 6; h++) {
      pars[h] = func0->GetParameter(h);
      if ( perr != NULL ) { 
         perr[h] = func0->GetParError(h);
      }
  }

  int nbin = ds->GetNbinsX() ;
  double bW = 480./nbin ;
  for (int i=1; i< nbin+1; i++) {
      double k = (i-1)*bW + bW/2.  ;
      double theN = func0->Eval(k);
      ds1->SetBinContent(i,theN);
  }

  double totalN = ds->Integral() ;
  ScaleTemplates( totalN, ds1, 1, nbin );

  delete func0;

}

void TemplateMassFit::combineBG( TString mName, TH1D* allbg, int rbin ) {

  int nbin = 480./rbin ;
  TH1D* tt = new TH1D("tt","", nbin, 0., 480.);
  TH1D* wj = new TH1D("wj","", nbin, 0., 480.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 480.);

  getBackground( tt, 1, rbin, mName );
  getBackground( wj, 2, rbin, mName );
  getBackground( qm, 3, rbin, mName );
  allbg->Add(tt, 1.);
  allbg->Add(wj, 1.);
  allbg->Add(qm, 0.00);

  delete tt;
  delete wj;
  delete qm;
}

void TemplateMassFit::getFakeData( int rbin, TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, TH1D* dth4, TH1D* dth5 ){

  int nbin = 480./rbin ;

  TFile *datafile1 = TFile::Open("pseud"+fname+"B0.root");
  if (dth5 != NULL ) getFakeData(dth5, "pseud_QCDB0", rbin , 2.01 );
  if (dth4 != NULL ) getFakeData(dth4, "pseud_STTWB0", rbin );
  if (dth3 != NULL ) getFakeData(dth3, "pseud_STTB0", rbin );
  if (dth2 != NULL ) getFakeData(dth2, "pseud_WJets_B0", rbin );

  // pseudo data from ttbar
  TH2D* hdata  = (TH2D *) datafile1->Get("Tops/"+hname);
  TH1D* hdata_pj = hdata->ProjectionX("hdata_pj",1,480,"");
  hdata_pj->Rebin(rbin);

  TH2D* hMC   = (TH2D *) datafile1->Get("MObjs/"+sname);
  TH1D* hMC_px = hMC->ProjectionX("hMC_px",1,480,"");
  TH1D* hMC_py = hMC->ProjectionY("hMC_py",1,480,"");
  hMC_px->Rebin(rbin);
  hMC_py->Rebin(rbin);

  TH1D* sg1 = new TH1D("sg1","", nbin, 0., 480.);
  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      if ( cname == "Had" ) theN = hMC_py->GetBinContent(i);
      if ( cname == "Lep" ) theN = hMC_px->GetBinContent(i);
      sg1->SetBinContent(i,theN);
  }
  TH1D* tt1 = new TH1D("tt1","", nbin, 0., 480.);
  tt1->Add(hdata_pj,sg1,1,-1);

  // sum all the ingredients 
  double allN[6] = {0.0};
  for (int i=1; i< nbin; i++) {
      allN[0] = sg1->GetBinContent(i);
      allN[1] = tt1->GetBinContent(i);
      if (dth2 != NULL ) allN[2] = dth2->GetBinContent(i);
      if (dth3 != NULL ) allN[3] = dth3->GetBinContent(i);
      if (dth4 != NULL ) allN[4] = dth4->GetBinContent(i);
      if (dth5 != NULL ) allN[5] = dth5->GetBinContent(i);
      double sumN = allN[0] + allN[1] + allN[2] + allN[3] + allN[4] + allN[5] ;
      ttadd->SetBinContent( i, sumN );
      dth0->SetBinContent( i, allN[0] );
      dth1->SetBinContent( i, allN[1] );
  }
  if (dth5 != NULL ) dth5->SetFillColor(5);  // QCD
  if (dth4 != NULL ) dth4->SetFillColor(6);  // single top tw channel
  if (dth3 != NULL ) dth3->SetFillColor(4);  // single top t channel
  if (dth2 != NULL ) dth2->SetFillColor(2);  // w+jets
  dth1->SetFillColor(7);
  dth0->SetFillColor(3);

  if (dth5 != NULL ) ttstk->Add( dth5 );
  if (dth4 != NULL ) ttstk->Add( dth4 );
  if (dth3 != NULL ) ttstk->Add( dth3 );
  if (dth2 != NULL ) ttstk->Add( dth2 );
  ttstk->Add( dth1 );
  ttstk->Add( dth0 );

  delete hMC_px;
  delete hMC_py;
  delete hMC;
  delete hdata_pj;
  delete hdata;
  delete tt1;
  delete sg1;
  
  datafile1->Close();
}

void TemplateMassFit::getFakeData( TH1D* ttadd, TString thefileName, int rbin, double theScale ){

  int nbin = 480./rbin ;
  TFile *datafile0  = TFile::Open(thefileName+".root");

  TH2D* fdata  = (TH2D *) datafile0->Get("Tops/"+hname);
  TH1D* fdata_pj = fdata->ProjectionX("fdata_pj",1,480,"");
  fdata_pj->Rebin(rbin);

  ttadd->Add(fdata_pj, theScale);
  ttadd->SetBinContent(nbin, 0.);
  /*
  for (int i=0; i< nbin; i++) {
      double binContain = fdata_pj->GetBinContent(i);
      ttadd->SetBinContent( i, binContain );
  }
  ttadd->Scale( theScale );
  */
  delete fdata_pj;
  delete fdata;

  datafile0->Close();

}

void TemplateMassFit::getData( TH1D* h_data, TString thefileName, int rbin ){

  int nbin = 480./rbin ;

  TFile *file0  = TFile::Open(thefileName+".root");
  TH2D* hdata  = (TH2D *) file0->Get("Tops/"+hname);

  TH1D* hdata_pj = hdata->ProjectionX("hdata_pj",1,200,"");
  hdata_pj->Rebin(rbin);

  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      theN  = hdata_pj->GetBinContent(i);
      h_data->SetBinContent(i,theN);
  }
  delete hdata_pj;
  delete hdata;
  file0->Close();

} 

void TemplateMassFit::getSignal(TH1D* h_Sg, int rbin, TString mName ) {

  int nbin = 480./rbin ;

  //TString theFileName = "tt_"+mName+"c.root" ;
  TString theFileName = "BTag_"+mName+"_TCHE5.root" ;

  thefile  = TFile::Open( theFileName );
  TH2D* ttMC  = (TH2D *) thefile->Get("MObjs/"+sname);

  TH1D* ttMC_px = ttMC->ProjectionX("ttMC_px",0,-1,"");
  TH1D* ttMC_py = ttMC->ProjectionY("ttMC_py",0,-1,"");
  ttMC_px->Rebin(rbin);
  ttMC_py->Rebin(rbin);
  //ttMC_px->Scale(0.428);
  //ttMC_py->Scale(0.428);
  

  double theN = 0. ;
  for (int i=1; i< nbin; i++) {
      if ( cname == "Had" ) theN = ttMC_py->GetBinContent(i);
      if ( cname == "Lep" ) theN = ttMC_px->GetBinContent(i);
      h_Sg->SetBinContent(i,theN);
  }

  //ScaleTemplates(100., 24000., 1., h_Sg );
  ScaleTemplates(100., 12000., 1., h_Sg );

  delete ttMC_px;
  delete ttMC_py;
  delete ttMC;
  thefile->Close();
}

void TemplateMassFit::getBackground(TH1D* h_Bg, int type, int rbin, TString mName ){

  int nbin = 480./rbin ;

  TString theFileName1 = "BTag_"+mName+"_TCHE5.root" ;
  TString theFileName2 = "BTag_WJets_TCHE5.root" ;
  TString theFileName3 = "BTag_STT_TCHE5.root" ;
  TString theFileName4 = "BTag_STTW_TCHE5.root" ;
  TString theFileName5 = "BTag_QCD_TCHE5.root" ;

  if (type == 1) thefile  = TFile::Open( theFileName1 );
  if (type == 2) file2  = TFile::Open( theFileName2 );
  if (type == 3) file2  = TFile::Open( theFileName3 );
  if (type == 4) file2  = TFile::Open( theFileName4 );

  if ( type == 1 ) ttmass = (TH2D *) thefile->Get("Tops/"+hname); 
  if ( type >= 2 ) ttmass = (TH2D *) file2->Get("Tops/"+hname); 
  
  TH1D* ttmass_pj = ttmass->ProjectionX("ttmass_pj",1,480,"");
  ttmass_pj->Rebin(rbin);

  if( type == 1 ) {

    TH1D* sgl = new TH1D("sgl","", nbin, 0., 480.);
    TH2D* sgnl   = (TH2D *) thefile->Get("MObjs/"+sname);
    TH1D* sgnl_px = sgnl->ProjectionX("sgnl_px",1,480,"");
    TH1D* sgnl_py = sgnl->ProjectionY("sgnl_py",1,480,"");
    sgnl_px->Rebin(rbin);
    sgnl_py->Rebin(rbin);
    double theN = 0. ;
    for (int i=1; i< nbin; i++) {
        if ( cname == "Had" ) theN = sgnl_py->GetBinContent(i);
        if ( cname == "Lep" ) theN = sgnl_px->GetBinContent(i);
        sgl->SetBinContent(i,theN);
    }

    ttmass_pj->Add(sgl, -1 );

    delete sgnl;
    delete sgnl_px;
    delete sgnl_py;
    delete sgl;
  }

  for (int i=1; i< nbin; i++) {
      double theN = ttmass_pj->GetBinContent(i);
      h_Bg->SetBinContent(i,theN);
  }

  // scale to 100 /pb
  if ( type == 1 ) ScaleTemplates( 100., 12000.,  1, h_Bg );  // ttbar
  if ( type == 2 ) ScaleTemplates( 100., 108980., 2, h_Bg );  // wjets
  if ( type == 3 ) ScaleTemplates( 100., 45018.,  3, h_Bg );  // single top t
  if ( type == 4 ) ScaleTemplates( 100., 25935.,  4, h_Bg );  // single top tW

  delete ttmass_pj;
  delete ttmass;

  if (type == 1 ) thefile->Close();
  if (type >= 2 ) file2->Close();

}

vector<TString> TemplateMassFit::FillMassAssumption( int npoints ){

TString marr5[]  = { "161", "166", "171", "176", "181" }; 
TString marr9[]  = { "161", "163", "166", "168", "171", "173", "176", "178", "181" }; 
TString marr10[] = { "161", "163", "166", "168", "171", "173", "176", "178", "181", "183" }; 

  vector<TString> mps;
  for (int k=0; k<npoints; k++) {
      if ( npoints ==  5 ) mps.push_back( marr5[k] );
      if ( npoints ==  9 ) mps.push_back( marr9[k] );
      if ( npoints == 10 ) mps.push_back( marr10[k] );
  }
  return mps ;
}

double TemplateMassFit::MassDigi( TString mString ) {

  double mDigi = 171.2;
  TString marr10[] = { "161", "163", "166", "168", "171", "173", "176", "178", "181", "183" }; 
  for (int k=0; k<10; k++) {
      if ( marr10[k] == mString ) mDigi = 161.2 + (2.5*k) ;
  }
  return mDigi ;
}

void TemplateMassFit::SetFitParameters( double mass, Double_t* para, int nPara, int NBTag ) {

     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau Wjets , 10: tt-ttwrong ratio , 11: tt-Wjets ratio 
     //                      p0      p1      p2     p3      p4     p5      p6     p7       p8      p9     p10      p11

     // measured by 9 parameters                      
     //Double_t b0[12] = { 26.736, 24.566, 23.000, 14.189, 5.550, 4.500, 91.761, 7.290, 116.716, 21.875, 25.346,  34.303 };
     //Double_t b1[12] = {  0.232,  0.861,  0.000,  0.098, 0.000, 0.000,  0.554, 0.214,   0.426,  0.152, -0.051,  -0.066 };

     // for 10 GeV Bin
     //Double_t a0[12] = { 26.438, 26.270, 23.000, 2.764, 5.550, 4.500, 101.981, 8.633, 197.690, 59.470, 25.159, 10.000 };
     //Double_t a1[12] = {  0.233,  0.851,  0.000, 0.166, 0.000, 0.000,   0.493, 0.206,   0.000,  0.000, -0.050, -0.021 };
     Double_t a0[12] = { 26.438, 26.270,  5.545, 2.764, 5.550, 4.500, 101.981, 7.290, 197.690, 59.470, 25.159, 10.000 };
     Double_t a1[12] = {  0.233,  0.851,  0.108, 0.166, 0.000, 0.000,   0.493, 0.214,   0.000,  0.000, -0.050, -0.021 };

     // for 5 GeV Bin
     //Double_t a0[12] = { 9.915, 25.160, 22.762, 17.651, 5.550, 4.500, 102.263, 9.471, 198.160, 57.9, 26.380, 10.579 };
     //Double_t a1[12] = { 0.138,  0.860,  0.001,  0.069, 0.000, 0.000,   0.491, 0.201,   0.000,  0.0, -0.058, -0.025 };


     // value for All BTagging
     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau/all Bg , 10: tt-ttwrong ratio , 11: tt-allBg ratio 
     //                      p0     p1      p2     p3     p4     p5      p6      p7      p8      p9     p10      p11
     //measured by CSV8
     //Double_t b0[12] = { -5.987, 5.768, -7.593, 51.243, 4.379, 4.803, 82.103, 23.808, 82.136, 24.728, 22.532,  22.613 };
     //Double_t b1[12] = {  0.316, 0.969,  0.172, -0.101, 0.006, 0.001,  0.591,  0.117,  0.592,  0.114, -0.047,  -0.048 };
     Double_t b0[12] = { -3.495, 7.495, -13.289, 67.106, 4.379,  5.162, 74.599, 14.40, 98.408, 20.172, 25.443,  38.488 };
     Double_t b1[12] = {  0.284, 0.960,   0.202, -0.186, 0.006, -0.001,  0.642,  0.17,  0.502,  0.144, -0.063,  -0.104 };
     // value for 1 BTagging
     Double_t c0[12] = {  8.264, 10.812, -13.0, 42.157, 4.379,  5.285, 103.493, 12.790, 97.420, 17.041, 25.294, 38.586 };
     Double_t c1[12] = {  0.128,  0.939,   0.2, -0.041, 0.006, -0.002,   0.480,  0.176,  0.526,  0.161, -0.046, -0.079 };
     // value for 2 BTagging
     Double_t d0[12] = { -14.213, -0.438, -10.402, 76.653, 4.379,  5.362, 30.953, 6.180, 62.611, 29.670, 22.394, 28.682 };
     Double_t d1[12] = {  -0.170,  1.005,   0.175, -0.258, 0.006, -0.002,  0.868, 0.206,  0.713,  0.084, -0.079, -0.101 };
     
     for ( int i =0; i < 12; i++ ) {
         if ( nPara == 12 && NBTag == -1 ) para[i] = a0[i] + a1[i]*mass ;
         if ( nPara == 9 && NBTag == 0 ) para[i] = b0[i] + b1[i]*mass ;
         if ( nPara == 9 && NBTag == 1 ) para[i] = c0[i] + c1[i]*mass ;
         if ( nPara == 9 && NBTag == 2 ) para[i] = d0[i] + d1[i]*mass ;
         if ( nPara != 9 && nPara != 12  ) para[i] = 0. ;
     }

}
