#include "TemplateMassFit.h"
#include "MassFitFunction.h"
TemplateMassFit::TemplateMassFit( TString aname, TString channel ){

  //file2  = TFile::Open("wjets_fall08_"+aname+".root");
  //file2  = TFile::Open("wj_fall08c.root");
  //file3  = TFile::Open("qcd_fall08_"+aname+".root");
  //file3  = TFile::Open("qcd_fall08c.root");

  fname = aname ;
  cname = channel ;
  hname = channel+"_selTmass_pt0";
  sname = "mcTtMass0";
  hfolder = "tt_fall08k10";


  plot1 = hfolder+"/"+channel+"_CombinedFit2_"+aname+".gif";
  plot2 = hfolder+"/"+channel+"_CombinedFit_"+aname+".gif";
  //plot3 = hfolder+"/"+channel+"_Norm_"+aname+".gif";
  plot4 = hfolder+"/"+channel+"_TemplateFit_"+aname+".gif";
  plot5 = hfolder+"/"+channel+"_Norm2_"+aname+".gif";
  plot6 = hfolder+"/"+channel+"_TemplateFit2_"+aname+".gif";
  plot8 = hfolder+"/"+channel+"_SBR_"+aname+".gif";
  plot9 = hfolder+"/"+channel+"_X2Test"+aname+".gif";

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

  
}

void TemplateMassFit::MoreCombinedFitting( int rbin, int lowBound, int upBound, int nMass ){

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
  TH1D* dth23 = new TH1D("dth23","", nbin, 0., 480.);

  getFakeData( fakedata, ttstk, dth0, dth1, dth2, dth3, rbin );

  // pre-fit background
  // smooth the background channels and get parameters for background shape
  Double_t bgpar[6] ;
  Double_t tbpar[6] ;
  Double_t sgpar[6] ;
  vector<TString> mvec = FillMassAssumption( nMass );
  vector<double> chi2(nMass, 99999.) ;
  
  for(int i = 0; i < nMass; i++) {
     chi2[i] = Chi2Test( mvec[i] , fakedata, rbin, lowBound, upBound, sgpar, tbpar, bgpar );
  }

  // record the process
  fprintf(logfile," === Chi2 Test for   === \n" );
  double minChi2 = 99999;
  int theBest = -1 ;
  for (int s1 = 0; s1 < nMass; s1++) {
      fprintf(logfile,"   %f \n", chi2[s1] );
      if ( chi2[s1] < minChi2 ) {
         minChi2 = chi2[s1] ;
         theBest = s1 ;
      }
  }
  fprintf(logfile," winner template is  %d \n", theBest );

  double bestMass = MinimumChi2( nMass, chi2 );

  //double fChi2 = -1; 
  //if ( theBest >= 0 ) fChi2 = Chi2Test( mvec[theBest], fakedata, rbin, lowBound, upBound, sgpar, tbpar, bgpar );
 
  // Fit the data
  TF1 *func1 = new TF1("func1",MassFitFunction::fitData1, lowBound, upBound,12);

  // get the extrapolated parameters
  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 12 );

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

  /*
  Double_t Rtw[5] = { 6.6,       6.4,       6.42,      6.5,       5.8 };
  Double_t Rtw[9] = { 6.6, 6.41, 6.4, 6.31, 6.42, 6.2, 6.5, 6.39, 5.8 };
  func1->SetParLimits(0,  25., 100.);
  func1->FixParameter(1, sgpar[1] );
  func1->FixParameter(2, sgpar[2] );
  func1->SetParLimits(3, 27., 45. );
  func1->FixParameter(4, sgpar[4] );
  func1->FixParameter(5, sgpar[5] );
  func1->FixParameter(6, tbpar[1] );
  func1->FixParameter(7, tbpar[2] );
  func1->FixParameter(8, bgpar[1] );
  func1->FixParameter(9, bgpar[2] );
  func1->SetParLimits(10, 12., 28.);
  func1->FixParameter(11,  Rtw[theBest] );
  */  

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
  dth23->Add(dth2,dth3, 1, 1);
  dth23->SetFillColor(6);
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
  delete dth23;
  delete fakedata;
}

void TemplateMassFit::CombinedFitting( int rbin, int lowBound, int upBound, int nMass ){

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
  TH1D* dth123 = new TH1D("dth23","", nbin, 0., 480.);

  getFakeData( fakedata, ttstk, dth0, dth1, dth2, dth3, rbin );

  // pre-fit background
  // smooth the background channels and get parameters for background shape
  Double_t bgpar[6] ;
  Double_t sgpar[6] ;
  vector<TString> mvec = FillMassAssumption( nMass );
  vector<double> chi2(nMass, 99999.) ;
  
  for(int i = 0; i < nMass; i++) {
     chi2[i] = Chi2Test( mvec[i] , fakedata, rbin, lowBound, upBound, sgpar, bgpar );
  }

  double bestMass = MinimumChi2( nMass, chi2 );

  // record the process
  fprintf(logfile," === Chi2 Test for   === \n" );
  double minChi2 = 99999;
  int theBest = -1 ;
  for (int s1 = 0; s1 < nMass; s1++) {
      fprintf(logfile,"   %f \n", chi2[s1] );
      if ( chi2[s1] < minChi2 ) {
         minChi2 = chi2[s1] ;
         theBest = s1 ;
      }
  }
  fprintf(logfile," winner template is  %d \n", theBest );

  // pre-fit background
  // smooth the background channels and get parameters for background shape
  /*
  Double_t bgpar[6] ;
  Double_t sgpar[6] ;
  double chi2[5] = {99999.};
  cout<<" *************************** "<<endl;
  chi2[0] = Chi2Test( "161" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );  
  chi2[1] = Chi2Test( "166" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );  
  chi2[2] = Chi2Test( "171" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );  
  chi2[3] = Chi2Test( "176" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );  
  chi2[4] = Chi2Test( "181" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );  
  cout<<" *************************** "<<endl;

  double minChi2 = 99999;
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
  if ( bestMass == 0 ) fChi2 = Chi2Test( "161" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );
  if ( bestMass == 1 ) fChi2 = Chi2Test( "166" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );
  if ( bestMass == 2 ) fChi2 = Chi2Test( "171" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );
  if ( bestMass == 3 ) fChi2 = Chi2Test( "176" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );
  if ( bestMass == 4 ) fChi2 = Chi2Test( "181" , fakedata, rbin, lowBound, upBound, sgpar, bgpar );

  cout<<"  final fitting chi2 = "<< fChi2 <<endl;
  */

  // Fit the data
  TF1 *func1 = new TF1("func1",MassFitFunction::fitData, lowBound, upBound,9);

  Double_t fpar[12];
  SetFitParameters( bestMass, fpar, 9 );

  func1->SetParLimits(0,  25., 100.);
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
  for (int i=0; i< 9; i++) { 
      apars[i]   = func1->GetParameter(i);
  }

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
  dth123->SetFillColor(6);
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
  delete dth123;
  delete fakedata;
 
  fclose(logfile);
}



double TemplateMassFit::MinimumChi2( int nMass, vector<double> chi2 ) {  

  c9 = new TCanvas("c9","", 800, 600);
  c9->SetGrid();
  c9->SetFillColor(10);
  c9->SetFillColor(10);
  c9->cd();

  const Int_t sz = nMass ;
  Double_t mInput[sz];
  Double_t nChi2[sz];
  for (int j=0; j< nMass; j++) {
      if ( nMass == 5 )  mInput[j] = 161.2 + j*5.  ;
      if ( nMass == 9 )  mInput[j] = 161.2 + j*2.5 ;
      nChi2[j]  = chi2[j];
  }

  TF1 *func5 = new TF1("func5",MassFitFunction::fitParabola, 160, 185, 3 );
  func5->SetParLimits(0, 150, 190);

  x2test = new TGraph(sz, mInput, nChi2);
  //x2test = new TGraph(5, mInput, chi2);
  //x2test = new TGraph(8, mInput, chi2);
  x2test->SetTitle(" chi2 test ");
  x2test->SetMarkerColor(4);
  x2test->SetMarkerStyle(21);
  x2test->SetMaximum(7.5);
  x2test->SetMinimum(0.0);
  x2test->GetXaxis()->SetTitle(" input template mass  ");
  x2test->GetYaxis()->SetTitle(" Normalized chi2  ");
  x2test->Draw("AP");
  x2test->Fit( func5, "R", "sames",155,185);
  c9->Update();
  c9->Print(plot9);

  double massCandidate =  func5->GetParameter(0);

  delete func5;
  delete x2test; 

  return massCandidate ;
}

// for multiple templates
double TemplateMassFit::Chi2Test( TString mName1, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sgpar, Double_t *tbpar, Double_t *bgpar ) {  

  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  int nbin = 480./rbin ;
 
  TH1D* bgadd = new TH1D("bgadd","", nbin, 0., 480.);
  TH1D* wj    = new TH1D("wj","", nbin, 0., 480.);
  TH1D* qm    = new TH1D("qm","", nbin, 0., 480.);
  getBackground( qm, 3, rbin, mName1 );
  bgadd->Add(qm, 0.0);
  getBackground( wj, 2, rbin, mName1 );
  bgadd->Add(wj);

  TH1D* parabg = new TH1D("parabg","", nbin, 0., 480.);
  SmoothTemplate( 0 , bgadd, parabg, lowBound, upBound, bgpar );

  TH1D* tt    = new TH1D("tt","", nbin, 0., 480.);
  getBackground( tt, 1, rbin, mName1 );

  TH1D* paratt = new TH1D("paratt","", nbin, 0., 480.);
  SmoothTemplate( 0 , tt, paratt, lowBound, upBound, tbpar );

  TH1D* sg = new TH1D("sg","", nbin, 0., 480.);
  getSignal( sg, rbin, mName1 );

  TH1D* parasg = new TH1D("parasg","", nbin, 0., 480.);
  SmoothTemplate( 1 , sg, parasg, lowBound, upBound, sgpar );

  Double_t Rtw = bgpar[0] / sgpar[0] ;
  Double_t Rtb = tbpar[0] / sgpar[0] ;

  // Fit the data

  TF1 *func1 = new TF1("func1",MassFitFunction::fitData1, lowBound, upBound, 12);
  func1->SetParLimits(0,  25., 100.);
  func1->FixParameter(1, sgpar[1] );
  func1->FixParameter(2, sgpar[2] );
  //func1->FixParameter(3, sgpar[3] );
  func1->SetParLimits(3, sgpar[3]-0.05*sgpar[3], sgpar[3]+0.05*sgpar[3] );
  func1->FixParameter(4, sgpar[4] );
  func1->FixParameter(5, sgpar[5] );
  func1->FixParameter(6, tbpar[1] );
  func1->FixParameter(7, tbpar[2] );
  func1->FixParameter(8, bgpar[1] );
  func1->FixParameter(9, bgpar[2] );
  //func1->SetParametr(10, 17.14);
  func1->SetParLimits(10, Rtb-0.05*Rtb , Rtb+0.05*Rtb );
  //func1->SetParLimits(11,  6., 7.5);
  func1->FixParameter(11,  Rtw);

  theData->Fit( func1, "RQ0","", lowBound, upBound);
 
  double x2 = func1->GetChisquare() / func1->GetNDF() ;
  //double x2 = func1->GetChisquare() ;
 
  cout<<" Fit X2 = "<< func1->GetChisquare() <<" / "<<func1->GetNDF() <<" = "<< x2 << endl;  
 
  delete bgadd;
  delete qm;
  delete wj;
  delete tt;
  delete sg;
  delete parabg;
  delete parasg;
  delete paratt;
  delete func1;

  return x2;
}

// for single background template
double TemplateMassFit::Chi2Test( TString mName1, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sgpar, Double_t *bgpar ) {  

  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  int nbin = 480./rbin ;
 
  TH1D* bgadd = new TH1D("bgadd","", nbin, 0., 480.);
  TH1D* tt    = new TH1D("tt","", nbin, 0., 480.);
  TH1D* wj    = new TH1D("wj","", nbin, 0., 480.);
  TH1D* qm    = new TH1D("qm","", nbin, 0., 480.);
  getBackground( qm, 3, rbin, mName1 );
  bgadd->Add(qm, 0.0);
  getBackground( wj, 2, rbin, mName1 );
  bgadd->Add(wj);
  getBackground( tt, 1, rbin, mName1 );
  bgadd->Add(tt);

  TH1D* parabg = new TH1D("parabg","", nbin, 0., 480.);
  SmoothTemplate( 0 , bgadd, parabg, lowBound, upBound, bgpar );

  TH1D* sg = new TH1D("sg","", nbin, 0., 480.);
  getSignal( sg, rbin, mName1 );

  TH1D* parasg = new TH1D("parasg","", nbin, 0., 480.);
  SmoothTemplate( 1 , sg, parasg, lowBound, upBound, sgpar );

  Double_t Rsb = bgpar[0] / sgpar[0] ;

  // Fit the data
  TF1 *func1 = new TF1("func1",MassFitFunction::fitData, lowBound, upBound, 9);
  func1->SetParLimits(0,  25., 100.);
  func1->FixParameter(1, sgpar[1] );
  func1->FixParameter(2, sgpar[2] );
  //func1->FixParameter(3, sgpar[3] );
  func1->SetParLimits(3, sgpar[3]-0.05*sgpar[3], sgpar[3]+0.05*sgpar[3] );
  func1->FixParameter(4, sgpar[4] );
  func1->FixParameter(5, sgpar[5] );
  func1->FixParameter(6, bgpar[1] );
  func1->FixParameter(7, bgpar[2] );
  func1->SetParLimits(8, Rsb-0.05*Rsb, Rsb+0.05*Rsb );

  theData->Fit( func1, "RQ0","", lowBound, upBound );
 
  double x2 = func1->GetChisquare() / func1->GetNDF() ;
 
  cout<<" Fit X2 = "<< func1->GetChisquare() <<" / "<<func1->GetNDF() <<" = "<< x2 << endl;  
 
  delete bgadd;
  delete qm;
  delete wj;
  delete tt;
  delete sg;
  delete parabg;
  delete parasg;
  delete func1;

  return x2;
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
  getFakeData( dt, ttstk, dth0, dth1, dth2, dth3, rbin );

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

    func0 = new TF1("func0", MassFitFunction::fitSG , 100, 350, 6);
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
    func1 = new TF1("func1", MassFitFunction::fitLD , 100, 380, 3);
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
    func2 = new TF1("func2", MassFitFunction::fitLD , 100, 380, 3);
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
  getFakeData( dt, ttstk, dth0, dth1, dth2, dth3, rbin );

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

    func0 = new TF1("func0", MassFitFunction::fitSG , 100, 350, 6);
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
    func1 = new TF1("func1", MassFitFunction::fitLD , 100, 380, 3);
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
  mfit->Constrain( 0, 0.0, 0.3);
  mfit->Constrain( 1, 0.3, 0.9);
  mfit->Constrain( 2, 0.1, 0.7);
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

void TemplateMassFit::FitTester( TString mName, int rbin, int lowBound, int upBound ){

  testfile = fopen(hfolder+"/testf.log","a"); 

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  plot7 = hfolder+"/_FitTest"+mName+".gif";

  c7 = new TCanvas("c7","", 1000, 800);
  c7->SetFillColor(10);
  c7->SetFillColor(10);
  c7->Divide(2,2);

  int nbin = 480./rbin ;
  Int_t b1 = lowBound / (rbin*1);
  Int_t b2 =  upBound / (rbin*1);
  const Int_t sz = b2 - b1 + 1 ;
  cout<<" b1:"<<b1<<" b2:"<<b2<<" size:"<< sz <<endl;
 
  // set up the y-axis max
  double yMax = 180. *( rbin / 20. ) ;

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
  func1->SetParLimits(1, 120.,220.);
  func1->SetParLimits(2, 10.,100.);
  sg->Fit( func1, "RQ0","", 140, 200 );

  double p0 = func1->GetParameter(0); 
  double p1 = func1->GetParameter(1); 
  double p2 = func1->GetParameter(2); 

  double r1 = p1 - 2*p2 ;
  double r2 = p1 + 2*p2 ;
  //double r1 = p1 - 30. ;
  //double r2 = p1 + 30. ;

  sg->Fit( func1, "R","sames", r1, r2 );
  p0 = func1->GetParameter(0); 
  p1 = func1->GetParameter(1); 
  p2 = func1->GetParameter(2); 

  c7->Update();

  c7->cd(2);

  TF1* func2 = new TF1("func2", MassFitFunction::fitSG , lowBound, upBound, 6);
  func2->FixParameter(0, p0);
  func2->FixParameter(1, p1);
  func2->FixParameter(2, p2);
  func2->SetParLimits(3, 10., 100.);
  //func2->SetParLimits(4, 3., 10.);
  func2->SetParLimits(5, 1., 10.);
  func2->FixParameter(4, 5.55 );

  sg->SetFillColor(7);
  sg->SetMaximum( yMax );
  sg->Draw();
  sg->Fit( func2, "RQ0","", lowBound, upBound); 

  double p3 = func2->GetParameter(3); 
  double p4 = func2->GetParameter(4); 
  double p5 = func2->GetParameter(5); 
  //func2->SetParLimits(4, p4 - sqrt(p4), p4 + sqrt(p4) );
  //func2->SetParLimits(5, p5 - sqrt(p5), p5 + sqrt(p5) );
  func2->FixParameter(4, 5.55 );
  //func2->FixParameter(5, 4.5 );
  func2->SetParLimits(5, 3. , 5. );
  sg->Fit( func2, "R","sames", lowBound, upBound); 
  p0 = func2->GetParameter(0); 
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

  TF1* func3 = new TF1("func3", MassFitFunction::fitSG , lowBound, upBound, 6);
  func3->SetParameter(0, p0);
  func3->SetParameter(1, p1);
  func3->SetParameter(2, p2);
  func3->SetParameter(3, p3);
  //func3->SetParameter(4, p4);
  func3->SetParameter(5, p5);
  /*
  func3->FixParameter(3, p3);
  func3->FixParameter(4, p4);
  func3->FixParameter(5, p5);
   */
  func3->SetParLimits(1, p1-0.1*p1, p1+0.1*p1);
  func3->SetParLimits(2, p2-0.1*p2, p2+0.1*p2);
  //func3->SetParLimits(3, p3-0.1*p3, p3+0.1*p3);
  func3->FixParameter(4, 5.55 );
  func3->SetParLimits(5, p5-0.1*p5, p5+0.1*p5);

  sg->SetFillColor(7);
  sg->SetMaximum( yMax );
  sg->Draw();
  sg->Fit( func3, "R","sames", lowBound, upBound); 

  // print out covariant matrics
  double mtx[5][5];
  gMinuit->mnemat(&mtx[0][0],5);
  for (int j=0; j < 5 ; j++) {
      for (int i=0; i<5; i++) {
          fprintf(testfile,"  %.3f",  mtx[i][j] );
          if ( i == 4 ) fprintf(testfile," \n" );
      }
  }
  fprintf(testfile," ------ \n" );

  p0 = func3->GetParameter(0); 
  p3 = func3->GetParameter(3); 
  p4 = func3->GetParameter(4); 
  p5 = func3->GetParameter(5); 

  func0->FixParameter(0, p0*p3);
  func0->FixParameter(1, p4   );
  func0->FixParameter(2, p5   );

  func0->SetLineColor(2);
  func0->SetLineWidth(3);
  func0->SetLineStyle(2);
  func0->Draw("sames");


  c7->Update();
  c7->Print(plot7);
 
  delete func3;
  delete func2;
  delete func1;
  delete sg;

  fclose(testfile);
}


void TemplateMassFit::TemplateDrawer( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp ){

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

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
  getBackground( tt, 1, rbin, fname );
  // tmp2: wjets and QCD 
  TH1D* wj = new TH1D("wj","", nbin, 0., 480.);
  TH1D* qm = new TH1D("qm","", nbin, 0., 480.);
  TH1D* bg = new TH1D("bg","", nbin, 0., 480.);
  getBackground( wj, 2, rbin, fname );
  getBackground( qm, 3, rbin, fname );
  if ( comp[1] ) bg->Add(tt,1.);
  if ( comp[2] ) bg->Add(wj,1.);
  if ( comp[3] ) bg->Add(qm,2.01);

  // get the fixed ratio btw wjets and ttbar signal
  Double_t nspar[6] ;
  Double_t ntpar[6] ;
  Double_t nbpar[6] ;
  TH1D* presg = new TH1D("presg","", nbin, 0., 480.);
  SmoothTemplate( 1 , sg, presg, lowBound, upBound, nspar );
  TH1D* prett = new TH1D("prett","", nbin, 0., 480.);
  SmoothTemplate( 0 , tt, prett, lowBound, upBound, ntpar );
  TH1D* prebg = new TH1D("prebg","", nbin, 0., 480.);
  SmoothTemplate( 0 , bg, prebg, lowBound, upBound, nbpar );

  double tWRatio = nbpar[0] / nspar[0] ;
  double tbRatio = ntpar[0] / nspar[0] ;
 
  c8 = new TCanvas("c8","", 1000, 800);
  c8->SetFillColor(10);
  c8->SetFillColor(10);
  c8->Divide(2,2);

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

  Double_t tmpars[10] ;
  for ( int i =0; i< 10; i++) {
      if ( i   < 6  )  tmpars[i] = nspar[i];
      if ( i  == 6  )  tmpars[i] = ntpar[1];
      if ( i  == 7  )  tmpars[i] = ntpar[2];
      if ( i  == 8  )  tmpars[i] = nbpar[1];
      if ( i  == 9  )  tmpars[i] = nbpar[2];
  }

  parfile = fopen(hfolder+"/paraf.log","a"); 
  for ( int i =0; i< 10; i++) {
      fprintf(parfile,"  %.2f",  tmpars[i] );
  }
  fprintf(parfile,"  %.2f",  tbRatio );
  fprintf(parfile,"  %.2f",  tWRatio );
  fprintf(parfile," \n" );
  fclose(parfile);

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

  // for 100 /pb => ttbar signal = 8992 events ; wjets = 30393 ; QCD = 2841112
  // channel #   => ttbar signal =  1 ; wjets = 2 ; QCD = 3 
  double lumiScale = lumi / 100 ;
  double nBase = nEvents;
  if ( channel == 1 ) nBase = 8992 ;
  if ( channel == 2 ) nBase = 30393*1.402 ;
  if ( channel == 3 ) nBase = 2841112 ;
  double Scal = (nBase*lumiScale) / nEvents ;
  tmp->Scale(Scal);
  
}

void TemplateMassFit::SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, Double_t* pars ){

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
     func0->SetParLimits(2, 22.,23.);
     //func0->FixParameter(2, 22.62);
     func0->SetParLimits(3, 27., 45.);
     //func0->SetParLimits(4, 5., 10.);
     //func0->SetParLimits(5, 0.3, 10.);
     func0->FixParameter(4, 5.55);
     func0->FixParameter(5, 4.5);
     
     // convoluted BW-Gaus
     //func0 = new TF1("func0", MassFitFunction::ConvBWGS , 100, 300, 5);
     //func0->SetParameters( bw0, bw1, bw2, bw1, 20. );
     //func0->SetParLimits(0, bw0-20., bw0+20.);
     //func0->SetParLimits(1, bw1-20., bw1+20.);

     // log-normal distribution
     /*
     ds->Fit( func0a, "R0","", 120, 230 );
     Double_t bw0 = func0a->GetParameter(0) ;
     Double_t bw1 = func0a->GetParameter(1) ;
     Double_t bw2 = func0a->GetParameter(2) ;
     func0 = new TF1("func0", MassFitFunction::fitSG , 100, 350, 6);
     func0->SetParLimits(0, bw0-20., bw0+20.);
     func0->SetParLimits(1, bw1-20., bw1+20.);
     func0->SetParLimits(2, bw2-100., bw2+100.);
     func0->SetParLimits(3, 10., 40000.);
     func0->SetParLimits(4, 0.01, 10000.);
     func0->SetParLimits(5, 0.1,  10000.);
     */

     // Landau distribution
     ds->Fit( func0, "RQ0","", 100, 350 );
  }

  for (int h=0; h< 6; h++) {
      pars[h] = func0->GetParameter(h);
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

void TemplateMassFit::getFakeData( TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, int rbin ){

  TFile *datafile1  = TFile::Open("pseud"+fname+"c.root");
  TFile *datafile2  = TFile::Open("pseudWjc.root");
  TFile *datafile3  = TFile::Open("pseudQcd.root");
  int nbin = 480./rbin ;

  // pseudo data from qcd
  
  qdata  = (TH2F *) datafile3->Get("Tops/"+hname);
  qdata->ProjectionX("qdata_pj",1,480,"");
  qdata_pj->Rebin(rbin);
  //qdata_pj->Scale(2.01);
  qdata_pj->Scale(0.0);
  
  // pseudo data from wjets
  wdata  = (TH2F *) datafile2->Get("Tops/"+hname);
  wdata->ProjectionX("wdata_pj",1,480,"");
  wdata_pj->Rebin(rbin);
  wdata_pj->Scale(1.402);

  // pseudo data from ttbar
  hdata  = (TH2F *) datafile1->Get("Tops/"+hname);
  hdata->ProjectionX("hdata_pj",1,480,"");
  hdata_pj->Rebin(rbin);

  hMC   = (TH2D *) datafile1->Get("MObjs/"+sname);
  hMC->ProjectionX("hMC_px",1,480,"");
  hMC->ProjectionY("hMC_py",1,480,"");
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
  double allN[4] = {0.0};
  for (int i=0; i< nbin; i++) {
      allN[0] = sg1->GetBinContent(i);
      allN[1] = tt1->GetBinContent(i);
      allN[2] = wdata_pj->GetBinContent(i);
      allN[3] = qdata_pj->GetBinContent(i);
      double sumN = allN[0] + allN[1] + allN[2] + allN[3] ;
      ttadd->SetBinContent( i, sumN );
      dth0->SetBinContent( i, allN[0] );
      dth1->SetBinContent( i, allN[1] );
      dth2->SetBinContent( i, allN[2] );
      dth3->SetBinContent( i, allN[3] );
  }
  dth3->SetFillColor(4);
  dth2->SetFillColor(2);
  dth1->SetFillColor(7);
  dth0->SetFillColor(3);
  ttstk->Add( dth3 );
  ttstk->Add( dth2 );
  ttstk->Add( dth1 );
  ttstk->Add( dth0 );

  delete hMC_px;
  delete hMC_py;
  delete hMC;
  delete hdata_pj;
  delete hdata;
  delete wdata_pj;
  delete wdata;
  delete qdata_pj;
  delete qdata;
  delete tt1;
  delete sg1;
  
  datafile1->Close();
  datafile2->Close();
  datafile3->Close();
}

void TemplateMassFit::getFakeData( TH1D* ttadd,  int rbin ){

  TFile *datafile1  = TFile::Open("pseud"+fname+"c.root");
  TFile *datafile2  = TFile::Open("pseudWjc.root");
  TFile *datafile3  = TFile::Open("pseudQcdc.root");

  // pseudo data from pythia tt
  int nbin = 480./rbin ;
  hdata  = (TH2F *) datafile1->Get("Tops/"+hname);
  hdata->ProjectionX("hdata_pj",1,480,"");
  hdata_pj->Rebin(rbin);

  // pseudo data from wjets
  wdata  = (TH2F *) datafile2->Get("Tops/"+hname);
  wdata->ProjectionX("wdata_pj",1,480,"");
  wdata_pj->Rebin(rbin);

  // pseudo data from qcd
  qdata  = (TH2F *) datafile3->Get("Tops/"+hname);
  qdata->ProjectionX("qdata_pj",1,480,"");
  qdata_pj->Rebin(rbin);
  qdata_pj->Scale(2.01);

  double allN = 0.0;
  for (int i=0; i< nbin; i++) {
      allN = hdata_pj->GetBinContent(i);
      allN += wdata_pj->GetBinContent(i);
      allN += qdata_pj->GetBinContent(i);
      ttadd->SetBinContent( i, allN );
  }
  
  //delete sg2;
  //delete tt2;
  //delete hMC_px;
  //delete hMC_py;
  //delete hMC;
  delete hdata_pj;
  delete hdata;
  delete wdata_pj;
  delete wdata;
  delete qdata_pj;
  delete qdata;

  datafile1->Close();
  datafile2->Close();
  datafile3->Close();

}
void TemplateMassFit::getData( TH1D* h_data, int rbin ){

  TFile *file0  = TFile::Open("pseud"+fname+".root");
  int nbin = 480./rbin ;
  hdata  = (TH2F *) file0->Get("Tops/"+hname);

  hdata->ProjectionX("hdata_pj",1,200,"");
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

  TString theFileName = "tt_"+mName+"c.root" ;

  thefile  = TFile::Open( theFileName );
  ttMC   = (TH2D *) thefile->Get("MObjs/"+sname);

  ttMC->ProjectionX("ttMC_px",1,480,"");
  ttMC->ProjectionY("ttMC_py",1,480,"");
  ttMC_px->Rebin(rbin);
  ttMC_py->Rebin(rbin);
  //ttMC_px->Scale(0.428);
  //ttMC_py->Scale(0.428);
  

  double theN = 0. ;
  for (int i=0; i< nbin; i++) {
      if ( cname == "Had" ) theN = ttMC_py->GetBinContent(i);
      if ( cname == "Lep" ) theN = ttMC_px->GetBinContent(i);
      h_Sg->SetBinContent(i,theN);
  }

  ScaleTemplates(100., 24000., 1., h_Sg );

  delete ttMC_px;
  delete ttMC_py;
  delete ttMC;
  thefile->Close();
}

void TemplateMassFit::getBackground(TH1D* h_Bg, int type, int rbin, TString mName ){

  TString theFileName1 = "tt_"+mName+"c.root" ;
  thefile  = TFile::Open( theFileName1 );
  if ( type < 2 ) ttmass = (TH2D *) thefile->Get("Tops/"+hname); 

  //file2  = TFile::Open("wjets_fall08_"+aname+".root");
  file2  = TFile::Open("wj_fall08c.root");
  //file3  = TFile::Open("qcd_fall08_"+aname+".root");
  file3  = TFile::Open("qcd_fall08c.root");

  if (type == 2) ttmass = (TH2D *) file2->Get("Tops/"+hname); 
  if (type == 3) ttmass = (TH2D *) file3->Get("Tops/"+hname);

  int nbin = 480./rbin ;

  ttmass->ProjectionX("ttmass_pj",1,480,"");
  ttmass_pj->Rebin(rbin);

  if( type == 1 ) {
    //ttmass_pj->Scale(0.428);

    TH1D* sgl = new TH1D("sgl","", nbin, 0., 480.);
    //getSignal( sgl, rbin, fname );

    TH2D* sgnl   = (TH2D *) thefile->Get("MObjs/"+sname);
    sgnl->ProjectionX("sgnl_px",1,480,"");
    sgnl->ProjectionY("sgnl_py",1,480,"");
    sgnl_px->Rebin(rbin);
    sgnl_py->Rebin(rbin);
    //sgnl_px->Scale(0.428);
    //sgnl_py->Scale(0.428);
    double theN = 0. ;
    for (int i=0; i< nbin; i++) {
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

  for (int i=0; i< nbin; i++) {
      double theN = ttmass_pj->GetBinContent(i);
      h_Bg->SetBinContent(i,theN);
  }

  if ( type == 1 ) ScaleTemplates( 100., 24000., 1, h_Bg );
  if ( type == 2 ) ScaleTemplates( 100., 86998., 2, h_Bg );

  delete ttmass_pj;
  delete ttmass;

  thefile->Close();
  file2->Close();
  file3->Close();
}

vector<TString> TemplateMassFit::FillMassAssumption( int npoints ){

  TString marr5[] = { "161", "166", "171", "176", "181" }; 
  TString marr9[] = { "161", "163", "166", "168", "171", "173", "176", "178", "181" }; 

  vector<TString> mps;
  for (int k=0; k<npoints; k++) {
      if ( npoints == 5 ) mps.push_back( marr5[k] );
      if ( npoints == 9 ) mps.push_back( marr9[k] );
  }
  return mps ;
}

void TemplateMassFit::SetFitParameters( double mass, Double_t* para, int nPara ) {

     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau Wjets , 10: tt-ttwrong ratio , 11: tt-Wjets ratio 
     //                      p0      p1      p2     p3      p4     p5      p6     p7       p8      p9     p10      p11
     Double_t a0[12] = { 822.587, 24.566, 23.000, 14.189, 5.550, 4.500, 91.761, 7.290, 197.690, 59.470, 25.346,  9.955 };
     Double_t a1[12] = {  -1.137,  0.861,  0.000,  0.098, 0.000, 0.000,  0.554, 0.214,   0.000,  0.000, -0.051, -0.021 };

     //  0~2: gaus , 3~5: log-normal ,  6,7:landau ttwrong ,  8,9:landau all Bg , 10: tt-ttwrong ratio , 11: tt-allBg ratio 
     //                      p0      p1      p2     p3      p4     p5      p6     p7       p8      p9     p10      p11
     Double_t b0[12] = { 26.736, 24.566, 23.000, 14.189, 5.550, 4.500, 91.761, 7.290, 116.716, 21.875, 25.346,  34.303 };
     Double_t b1[12] = {  0.232,  0.861,  0.000,  0.098, 0.000, 0.000,  0.554, 0.214,   0.426,  0.152, -0.051,  -0.066 };

     for ( int i =0; i < 12; i++ ) {
         if ( nPara == 12 ) para[i] = a0[i] + a1[i]*mass ;
         if ( nPara == 9  ) para[i] = b0[i] + b1[i]*mass ;
         if ( nPara != 9 && nPara != 12  ) para[i] = 0. ;
     }

}
