#include "MassAna.h"

MassAna::MassAna( TString channel, int NBTag, double massL, double massH  ){

  ptype  = ".gif";

  cname = channel;

  mL = massL ;
  mH = massH ;

  fitFunc  = new MassFitFunction();
  fitInput = new MassAnaInput( channel, NBTag, massL, massH );
  fitInput->Initialize( &hfolder  );

  gSystem->mkdir(hfolder);

}

MassAna::~MassAna(){

 delete fitFunc;
 delete fitInput;

}

// chi2 test + minimum chi2 fit, output best estimated mass
double MassAna::Chi2Test( TString mName, TH1D* theData, int lowBound, int upBound, int nPar, int NBTag, Double_t* statErr, int rbin, bool isWeight ) {  

   plot9 = hfolder+"Chi2Test_"+mName+".gif";

   gStyle->SetOptFit(111);
   c9 = new TCanvas("c9","", 800, 600);
   c9->SetGrid();
   c9->SetFillColor(10);
   c9->SetFillColor(10);
   c9->cd();

   Double_t mAssumption[16] = {0.}; 
   Double_t nChi2[16] ={0.}; 
   Double_t fpar[12];
   double maxChi2 = 0. ;
   TF1 *func0 = 0;
   if ( cname == "had" && !isWeight ) func0 = new TF1("func0", MassFitFunction::fitData1, lowBound, upBound, 12);
   if ( cname == "lep" ||  isWeight ) func0 = new TF1("func0", MassFitFunction::fitData2, lowBound, upBound, 12);
   TF1 *fS  =0;
   if ( cname == "had" && !isWeight ) fS    = new TF1("fS"   , MassFitFunction::fitSG,  lowBound, upBound, 6);
   if ( cname == "lep" ||  isWeight ) fS    = new TF1("fS"   , MassFitFunction::fitSG1, lowBound, upBound, 6);

   TF1 *fW    = new TF1("fW", MassFitFunction::fitLD, lowBound, upBound, 3);
   TF1 *fB    = new TF1("fB", MassFitFunction::fitLD, lowBound, upBound, 3);
   for (int i =0; i < 16; i++ ) {

       // set up the mass assumption
       mAssumption[i] = 156. + (i*2.) ;
       SetFitParameters( mAssumption[i], fpar, nPar, NBTag, rbin );

       // Fit the data
       func0->SetParLimits(0,  10., 100.);
       func0->FixParameter(1, fpar[1] );
       func0->FixParameter(2, fpar[2] );
       func0->SetParLimits(3, fpar[3]-0.1*fpar[3], fpar[3]+0.1*fpar[3] );
       func0->FixParameter(4, fpar[4] );
       func0->FixParameter(5, fpar[5] );
       func0->FixParameter(6, fpar[6] );
       func0->FixParameter(7, fpar[7] );
       func0->FixParameter(8, fpar[8] );
       func0->FixParameter(9, fpar[9] );
       func0->SetParLimits(10,fpar[10]-0.1*fpar[10], fpar[10]+0.1*fpar[10] );
       func0->FixParameter(11,fpar[11]);

       theData->Fit( func0, "RQ0","", lowBound, upBound);

       nChi2[i] = getChi2(theData, func0, fS, fB, fW, 130, 210, nPar );
       if ( nChi2[i] > maxChi2 ) maxChi2 = nChi2[i] ;
   }

   double mErr[16]={0.};
   double x2Err[16]={0.};
   TGraph* x2test = new TGraph(16, mAssumption, nChi2);
   
   x2test->SetTitle(" chi2 test ");
   x2test->SetMarkerColor(4);
   x2test->SetMarkerStyle(21);
   x2test->SetMaximum(maxChi2*1.25);
   x2test->SetMinimum(0.0);
   x2test->GetXaxis()->SetTitle(" input mass assumption ");
   x2test->GetYaxis()->SetTitle(" chi2  ");
   x2test->Draw("AP");


   TF1 *func5 = new TF1("func5",MassFitFunction::fitParabola, 154, 188, 3 );
   func5->SetParLimits(0, 150, 190);
   func5->SetParameter(1, 3.5 );
   x2test->Fit( func5, "NR0", "",154,188 );

   std::vector<bool> rejList = fitFunc->DataRejection( func5, mAssumption, nChi2, 16);
   
   for (int i=0; i<16; i++) {
       if ( !rejList[i] ) x2Err[i] = 0.1 ;
       if (  rejList[i] ) x2Err[i] = 10. ;
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

   // re-fit the chi2 curve
   if ( func5->GetChisquare() > 650.  ) {
      double p1 = func5->GetParameter(0);
      double L1 = p1 - (p1 - 158)/2 ;
      double H1 = p1 + (188 - p1)/2 ;
      if ( p1 > 188 ) {
         L1 = 171;
         H1 = 188;
      }
      if ( p1 < 158 ) {
         L1 = 158;
         H1 = 171;
      }
      x2testErr->Fit( func5, "R", "sames", L1, H1 );
   }


   c9->Update();
   c9->Print(plot9);

   double massCandidate =  func5->GetParameter(0);
   if (statErr != NULL ) *statErr = func5->GetParameter(1);

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
 
double MassAna::getChi2( TH1D* theData1,  TF1* theFunc, TF1* fS, TF1* fB, TF1* fW, double lowBound, double upBound, int nPar ) {

  for (int s=0; s<6; s++) {
      fS->FixParameter( s, theFunc->GetParameter(s) );
  }

  fW->FixParameter( 0, theFunc->GetParameter(11)*theFunc->GetParameter(0) );
  fW->FixParameter( 1, theFunc->GetParameter(8) );
  fW->FixParameter( 2, theFunc->GetParameter(9) );

  fB->FixParameter( 0, theFunc->GetParameter(10)*theFunc->GetParameter(0) );
  fB->FixParameter( 1, theFunc->GetParameter(6) );
  fB->FixParameter( 2, theFunc->GetParameter(7) );

  int nbin = theData1->GetNbinsX() ;
  double bW = (mH -mL) / nbin ;
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
      double x2 = ( theH - theF )*( theH - theF ) / ( theF + theB + theW ); 
      if ( nPar == 9 ) x2 = ( theH - theF )*( theH - theF ) / ( theF + theB ); 
      chi2 += x2;
  }
  //chi2 = chi2 / ( b2+1 - b1 - 3 ) ;
  
  return chi2;
}


void MassAna::GetAllCoeff( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp ) {

     Double_t sPar[6];
     Double_t sErr[6];
     Double_t bPar[6];
     Double_t bErr[6];
     if ( comp == NULL ) {
       FitTtbar( mName, rbin, sPar, sErr );
       FitBackground( mName, rbin, lowBound, upBound, bPar, bErr );
     } else {
       FitSignal( mName, rbin, sPar, sErr );
       FitBackground( mName, rbin, lowBound, upBound, comp, bPar, bErr );
     }
     parfile = fopen(hfolder+"/paraf.log","a");
     errfile = fopen(hfolder+"/perrf.log","a");
     for ( int i =0; i< 6; i++) {
         fprintf(parfile,"  %.3f",  sPar[i] );
         fprintf(errfile,"  %.3f",  sErr[i] );
     }
     for ( int i =0; i< 6; i++) {
         if ( i != 0 && i != 3 ) {
            fprintf(parfile,"  %.3f",  bPar[i] );
            fprintf(errfile,"  %.3f",  bErr[i] );
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

void MassAna::FitSignal(TString mName, int rbin, Double_t* para, Double_t* perr) {

  FILE* testfile = fopen(hfolder+"/Sgpara.log","a");

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  TString plot7 = hfolder+"FitTest_"+mName+".gif";

  c7 = new TCanvas("c7","", 1000, 800);
  c7->SetFillColor(10);
  c7->SetFillColor(10);
  c7->Divide(2,2);

  double m1 = MassDigi(mName);
  double lowBound = m1 - 50 ;
  double upBound  = m1 + 170 ;
  
  c7->cd(1);

  // tmp0: get template signal distribution
  int nbin = ( mH - mL ) / rbin ;
  TH1D* sg = new TH1D("sg","", nbin, mL, mH );
  fitInput->getSignal( sg, 1, mName );

  sg->SetFillColor(7);
  //sg->SetMaximum( yMax );
  sg->Draw();

  TF1* func0 = new TF1("func0", MassFitFunction::fitLG , 80, 450, 3);
  TF1* func7 = new TF1("func7", MassFitFunction::fitGS , lowBound, upBound, 3);
  TF1* func2 = new TF1("func2", MassFitFunction::fitSG , lowBound-20, upBound, 6);
  if ( cname == "Lep" ) {
     func0 = new TF1("func0", MassFitFunction::fitLD , 80, 450, 3);
     func2 = new TF1("func2", MassFitFunction::fitSG1 , lowBound-20, upBound, 6);
  }
  TF1* func3 = new TF1("func3", MassFitFunction::fitLD , 0, 480, 3);

  // pre-set the value
  // Hadronic top : gaus + log-normal
  // Leptonic top : gaus + landau
  double p0 = 50. ;
  double p1 = m1  ;
  double p2 = 20. ;
  double p3 = 33. ;
  double p4 = log(m1) + sqrt(1./20.) ;
  double p5 =  5. ;
  if ( cname == "Lep" ) {
     p3 = 1.8;
     p4 = m1 ;
     p5 = 30 ;
  }
  // 1st Fit, Fix "mean" value for gaussisan & log-normal and allow normalization and width vary
  func2->SetParLimits( 0, 10., p0+1.0*p0);
  func2->FixParameter( 1, p1 );
  func2->SetParLimits( 2, p2-0.3*p2, p2+1.0*p2);
  if ( cname == "Had" )  func2->SetParLimits(3, p3-0.1*p3, p3+0.1*p3);
  if ( cname == "Lep" )  func2->SetParLimits(3, p3-0.3*p3, p3+0.3*p3);

  func2->FixParameter(4, p4 );
  func2->FixParameter(5, p5 );

  sg->Fit( func2, "RQ0","", lowBound, upBound );

  p0 = func2->GetParameter(0);
  p1 = func2->GetParameter(1);
  p2 = func2->GetParameter(2);
  p3 = func2->GetParameter(3);
  p5 = func2->GetParameter(5);

  // Draw the gaussian
  func7->FixParameter(0, p0 );
  func7->FixParameter(1, p1 );
  func7->FixParameter(2, p2 );
  func7->SetLineColor(4);
  func7->Draw("sames");

  c7->Update();

  /*
  c7->cd(4);
  func3->FixParameter(0, p3*p0 );
  func3->FixParameter(1, p4 );
  func3->FixParameter(2, p5 );
  func3->Draw();
  c7->Update();
  */

  c7->cd(2);

  sg->SetFillColor(7);
  //sg->SetMaximum( yMax );
  sg->Draw();

  // 2nd Fit , Allow gaussian change
  func2->SetParLimits(0, p0-0.1*p0, p0+0.1*p0);
  func2->SetParLimits(1, m1 - 5., m1 + 5. );
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
  //sg->SetMaximum( yMax );
  sg->Draw();
  sg->Fit( func2, "R","sames", lowBound, upBound);
  func2->Draw("sames");

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

  delete c7;
  delete func3;
  delete func2;
  delete func7;
  delete func0;
  delete sg;

  fclose(testfile);

}

// for kinematic constrain
void MassAna::FitTtbar(TString mName, int rbin, Double_t* para, Double_t* perr) {

  FILE* testfile = fopen(hfolder+"Sgpara.log","a");

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  TString plot7 = hfolder+"FitTest_"+mName+".gif";

  c7 = new TCanvas("c7","", 1000, 800);
  c7->SetFillColor(10);
  c7->SetFillColor(10);
  c7->Divide(2,2);

  double m1 = MassDigi(mName);
  double lowBound = m1 - 50 ;
  double upBound  = m1 + 170 ;
  int nbin = ( mH - mL ) / rbin ;
  
  TF1* func7 = new TF1("func7", MassFitFunction::fitGS , lowBound, upBound, 3);
  TF1* func0 = new TF1("func0", MassFitFunction::fitLD , 80, 450, 3);
  TF1* func2 = new TF1("func2", MassFitFunction::fitSG1 , lowBound, upBound, 6);

  // A testing platform
  c7->cd(4);

  gStyle->SetOptStat("ei");
  gStyle->SetStatY(0.99);
  TString theBrName = cname+"TM" ;

  TH1D* mph  = new TH1D("mph","", nbin, mL, mH );
  fitInput->getMostProb( mName, theBrName, mph );
  mph->Draw();
  c7->Update();

  gStyle->SetStatY(0.84);
  gStyle->SetStatTextColor(4);
  TH1D* mch  = new TH1D("mch","", nbin, mL, mH );
  fitInput->getMcMatching( mName, theBrName, mch );
  mch->SetLineColor(4);
  func7->SetParLimits( 0, 1, 20 );
  func7->SetParameter( 1, m1 );
  func7->SetParLimits( 2, 10, 30 );
  func7->SetLineColor(4);
  mch->Fit( func7, "R","sames", lowBound-20, upBound );
  c7->Update();

  gStyle->SetStatY(0.55);
  gStyle->SetStatTextColor(2);
  TH1D* tbg  = new TH1D("tbg","", nbin, mL, mH );
  tbg->Add( mph );
  tbg->Add( mch, -1. );
  tbg->SetLineColor(2);
  func0->SetParLimits( 0, 1, 100 );
  func0->SetParLimits( 1, m1+10, m1+200 );
  func0->SetParLimits( 2, 5, 55 );
  func0->SetLineColor(2);
  tbg->Fit( func0, "R","sames", lowBound-20, upBound );
  c7->Update();

  c7->cd(1);

  gStyle->SetStatTextColor(1);
  gStyle->SetStatY(0.95);
  // tmp0: get template signal distribution
  TH1D* sg = new TH1D("sg","", nbin, mL, mH );
  fitInput->getSignal( sg, 0, mName );

  sg->SetFillColor(7);
  //sg->SetMaximum( yMax );
  sg->Draw();

  // pre-set the value
  double p0 = 50. ;
  double p1 = m1  ;
  double p2 = func7->GetParameter(2);
  double p3 = func0->GetParameter(0);
  double p4 = func0->GetParameter(1);
  double p5 = func0->GetParameter(2);

  // 1st Fit, Fix "mean" value for gaussisan & landau and allow normalization and width vary
  func2->SetParLimits( 0, p0-0.1*p0, p0+1.0*p0);
  func2->FixParameter( 1, p1 );
  func2->SetParLimits( 2, p2-0.2*p2, p2+0.2*p2);
  func2->SetParLimits( 3, 0.1, 8 );
  func2->SetParameter( 4, p4 );
  func2->SetParLimits( 5, p5-0.1*p5, p5+0.1*p5);

  sg->Fit( func2, "RQ0","", lowBound, upBound );
  p0 = func2->GetParameter(0);
  p1 = func2->GetParameter(1);
  p2 = func2->GetParameter(2);
  p3 = func2->GetParameter(3);
  p5 = func2->GetParameter(5);

  // Draw the gaussian
  func7->FixParameter(0, p0 );
  func7->FixParameter(1, p1 );
  func7->FixParameter(2, p2 );
  func7->SetLineColor(4);
  func7->Draw("sames");
  c7->Update();


  c7->cd(2);

  gStyle->SetOptStat("nirm");
  sg->SetFillColor(7);
  //sg->SetMaximum( yMax );
  sg->Draw();

  // 2nd Fit , Allow gaussian change
  func2->SetParLimits(0, p0-0.1*p0, p0+0.1*p0);
  func2->SetParLimits(1, m1 - 5.  , m1 + 5. );
  func2->SetParLimits(2, p2-0.1*p2, p2+0.1*p2 );
  func2->SetParLimits(3, p3-0.1*p3, p3+0.1*p3);
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
  func2->SetParLimits(4, p4-0.01*p4, p4+0.01*p4);
  func2->SetParLimits(5, p5-0.01*p5, p5+0.01*p5 );
  sg->SetFillColor(7);
  //sg->SetMaximum( yMax );
  sg->Draw();
  sg->Fit( func2, "R","sames", lowBound, upBound);
  func2->Draw("sames");

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

  delete c7;
  delete func0;
  delete func2;
  delete func7;
  delete sg;
  delete tbg;
  delete mch;
  delete mph;

  fclose(testfile);

}


void MassAna::FitBackground( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp, Double_t *para, Double_t *perr ){

  FILE* bgpara = fopen(hfolder+"Bgpara.log","a"); 

  gStyle->SetOptFit(111);
  gStyle->SetOptStat("nirm");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(1);

  double m1 = MassDigi(mName);
  plot8 = hfolder+cname+"_BG_"+mName+".gif";

  int nbin = ( mH - mL ) / rbin ;
  TH1D* tt  = new TH1D("tt",  "", nbin, mL, mH );   // tmp1: tt wrong combinatorics
  TH1D* wj  = new TH1D("wj",  "", nbin, mL, mH );   // tmp2: wjets, single Top t, single Top tW, QCD
  TH1D* stt = new TH1D("stt", "", nbin, mL, mH );
  TH1D* stw = new TH1D("stw", "", nbin, mL, mH );
  fitInput->getBackground( tt,  1, nbin, mName );
  fitInput->getBackground( wj,  2, nbin, mName );
  fitInput->getBackground( stt, 3, nbin, mName );
  fitInput->getBackground( stw, 4, nbin, mName );
  // mix the background
  TH1D* bg1  = new TH1D("bg1","Background group 1", nbin, mL, mH );
  if ( comp[1] ) bg1->Add(tt, 1.);
  if ( comp[2] ) bg1->Add(wj, 1.);
  if ( comp[3] ) bg1->Add(stt,1.);
  if ( comp[4] ) bg1->Add(stw,1.);

  TH1D* bg2 = new TH1D("bg2","Background group 2", nbin, mL, mH );
  if ( !comp[1] ) bg2->Add(tt, 1.);
  if ( !comp[2] ) bg2->Add(wj, 1.);
  if ( !comp[3] ) bg2->Add(stt,1.);
  if ( !comp[4] ) bg2->Add(stw,1.);

  c8 = new TCanvas("c8","", 800, 600);
  c8->SetFillColor(10);
  c8->SetFillColor(10);
  c8->Divide(2,2);

  // show all background compositions
  c8->cd(1);
  gStyle->SetOptStat("ni");
  THStack* allbgstk = new THStack("allbgstk", "Sum of backgrounds");

  stt->SetFillColor(4);
  allbgstk->Add(stt);
  stw->SetFillColor(6);
  allbgstk->Add(stw);
  wj->SetFillColor(2);
  allbgstk->Add(wj);
  tt->SetFillColor(7);
  allbgstk->Add(tt);
  allbgstk->Draw();

  gStyle->SetStatY(0.30);
  gStyle->SetStatTextColor(2);
  wj->SetName(" W + Jets ");
  wj->DrawCopy("sames");
  c8->Update();
  gStyle->SetStatY(0.50);
  gStyle->SetStatTextColor(4);
  stt->SetName("SingleTop_t");
  stt->DrawCopy("sames");
  c8->Update();
  gStyle->SetStatY(0.70);
  gStyle->SetStatTextColor(6);
  stw->SetName("SingleTop_tW");
  stw->DrawCopy("sames");
  c8->Update();
  gStyle->SetStatY(0.90);
  gStyle->SetStatTextColor(7);
  tt->SetName("Tt-Wrong");
  tt->DrawCopy("sames");
  c8->Update();

  allbgstk->Draw("same");
  c8->Update();
  gStyle->SetStatTextColor(1);
 
  // sum of backgrounds
  c8->cd(3);
  THStack* bg1stk = new THStack("bg1stk", "background group1 ");
  if ( comp[3] ) bg1stk->Add(stt);
  if ( comp[4] ) bg1stk->Add(stw);
  if ( comp[2] ) bg1stk->Add(wj);
  if ( comp[1] ) bg1stk->Add(tt);
  bg1stk->Draw();

  // Fit the combined backgrounds
  gStyle->SetStatTextColor(1);
  TF1* func7 = new TF1("func7", MassFitFunction::fitLD , 90, 380, 3);
  func7->SetParLimits(0, 100., bg1->Integral() );
  if (cname == "Had" ) func7->SetParLimits(1, m1-10, m1+20.);
  if (cname == "Lep" ) func7->SetParLimits(1, m1-20, m1+15.);
  func7->SetParLimits(2,  1., 100.);
  func7->SetLineColor(1);
  bg1->Fit( func7, "R","sames", lowBound, upBound );

  c8->Update();

  // sum of top related backgrounds
  c8->cd(4);
  gStyle->SetOptStat("ni");
  THStack* bg2stk = new THStack("bg2stk", "background group 2");
  if ( !comp[3] ) bg2stk->Add(stt);
  if ( !comp[4] ) bg2stk->Add(stw);
  if ( !comp[2] ) bg2stk->Add(wj);
  if ( !comp[1] ) bg2stk->Add(tt);

  bg2stk->Draw();

  gStyle->SetStatTextColor(1);
  // Fit the combined backgrounds
  TF1* func2 = new TF1("func2", MassFitFunction::fitLD , 90, 380, 3);
  func2->SetParLimits(0, 10., bg2->Integral() );
  func2->SetParLimits(1, m1-10, m1+20.);
  func2->SetParLimits(2,  1., 100.);
  bg2->Fit( func2, "N0R","", lowBound, upBound );
  double p0 = func2->GetParameter(0);
  double p1 = func2->GetParameter(1);
  double p2 = func2->GetParameter(2);
  func2->SetParLimits(0, p0-0.1*p0, p0+0.1*p0 );
  func2->SetParLimits(1, p1-0.1*p1, p1+0.1*p1 );
  func2->SetParLimits(2, p2-0.1*p2, p2+0.1*p2 );
  bg2->Fit( func2, "R","sames", lowBound, upBound );

  if ( comp[1] && comp[2] && comp[3] && comp[4] ) {
     bg2stk->Add(stt);
     bg2stk->Add(stw);
     bg2stk->Add(wj);
     bg2stk->Draw();
  }

  c8->Update();

  // Draw the tt-wrong combination background
  c8->cd(2);
  tt->SetLineWidth(2);
  tt->SetTitle("tt-wrong combinations");
  tt->Draw();

  // Fit tt-wrong permutation
  TF1* func0 = new TF1("func0", MassFitFunction::fitLD , 90, 380, 3);
  func0->SetParLimits(0, 100., tt->Integral() );
  if (cname == "Had" ) func7->SetParLimits(1, m1-10, m1+20.);
  if (cname == "Lep" ) func7->SetParLimits(1, m1-20, m1+15.);
  func0->SetParLimits(2,  1., 100.);
  func0->SetLineStyle(2);
  tt->Fit( func0, "R","sames", lowBound, upBound );
  c8->Update();

  c8->Print(plot8);


  // para[0~2] : background group 1 / tt wrong permutation 
  // para[3~5] : background group 2 / combined backgrounds 
  fprintf(bgpara," %.1f", m1 );
  for (int i=0; i<6; i++) {  
      if ( i < 3 ) {
         fprintf(bgpara,"  %.3f  %.3f",  func7->GetParameter(i), func7->GetParError(i) );
         if ( para != NULL ) para[i] = func7->GetParameter(i);
         if ( perr != NULL ) perr[i] = func7->GetParError(i) ;
      } else {
         fprintf(bgpara,"  %.3f  %.3f",  func2->GetParameter(i-3), func2->GetParError(i-3) );
         if ( para != NULL ) para[i] = func2->GetParameter(i-3);
         if ( perr != NULL ) perr[i] = func2->GetParError(i-3) ;
      }
  }
  fprintf(bgpara," \n" );

  delete c8;
  delete bg1;
  delete bg2;
  delete tt;
  delete wj;
  delete stt;
  delete stw;
  delete func0;
  delete func7;
  delete func2;
  delete bg1stk;
  delete bg2stk;
  delete allbgstk;

  fclose(bgpara);
}

// another method for kinematic constrain case
void MassAna::FitBackground( TString mName, int rbin, int lowBound, int upBound, Double_t *para, Double_t *perr ){

  FILE* bgpara = fopen(hfolder+"/Bgpara.log","a"); 

  double m1 = MassDigi(mName);
  plot8 = hfolder+"/"+cname+"_BG_"+mName+".gif";

  int nbin = ( mH - mL ) / rbin ;
  TH1D* wj  = new TH1D("wj",  "", nbin, mL, mH );   // tmp2: wjets, single Top t, single Top tW, QCD
  TH1D* stt = new TH1D("stt", "", nbin, mL, mH );
  TH1D* stw = new TH1D("stw", "", nbin, mL, mH );
  fitInput->getBackground( wj,  2, nbin, mName );
  fitInput->getBackground( stt, 3, nbin, mName );
  fitInput->getBackground( stw, 4, nbin, mName );

  TH1D* bg1  = new TH1D("bg1","Background group 1", nbin, mL, mH );
  bg1->Add(wj, 1.);
  bg1->Add(stt,1.);
  bg1->Add(stw,1.);

  c8 = new TCanvas("c8","", 800, 600);
  c8->SetFillColor(10);
  c8->SetFillColor(10);

  // show the backgrounds compositions
  c8->cd();
  gStyle->SetOptStat("ni");
  THStack* allbgstk = new THStack("allbgstk", "Sum of backgrounds");

  stt->SetFillColor(4);
  allbgstk->Add(stt);
  stw->SetFillColor(6);
  allbgstk->Add(stw);
  wj->SetFillColor(2);
  allbgstk->Add(wj);
  allbgstk->Draw();

  gStyle->SetStatY(0.30);
  gStyle->SetStatTextColor(2);
  wj->SetName(" W + Jets ");
  wj->DrawCopy("sames");
  c8->Update();
  gStyle->SetStatY(0.45);
  gStyle->SetStatTextColor(4);
  stt->SetName("SingleTop_t");
  stt->DrawCopy("sames");
  c8->Update();
  gStyle->SetStatY(0.60);
  gStyle->SetStatTextColor(6);
  stw->SetName("SingleTop_tW");
  stw->DrawCopy("sames");
  c8->Update();

  gStyle->SetStatY(0.95);
  allbgstk->Draw("same");
  c8->Update();
  gStyle->SetStatTextColor(1);

  // Fit the combined backgrounds
  gStyle->SetStatTextColor(1);
  TF1* func7 = new TF1("func7", MassFitFunction::fitLD , 90, 350, 3);
  func7->SetParLimits(0, 10., bg1->Integral() );
  func7->SetParLimits(1, 180, 220.);
  func7->SetParLimits(2,  10., 200.);
  func7->SetLineColor(1);
  bg1->Fit( func7, "R","sames", lowBound, upBound );

  c8->Update();
  c8->Print(plot8);

  fprintf(bgpara," %.1f", m1 );
  for (int i=0; i<3; i++) {  
      fprintf(bgpara,"  %.3f  %.3f",  func7->GetParameter(i), func7->GetParError(i) );
      if ( para != NULL ) para[i] = func7->GetParameter(i);
      if ( perr != NULL ) perr[i] = func7->GetParError(i) ;
  }
  fprintf(bgpara," \n" );

  delete c8;
  delete bg1;
  delete wj;
  delete stt;
  delete stw;
  delete func7;
  delete allbgstk;

  fclose(bgpara);
}

double MassAna::MassDigi( TString mString ) {

  double mDigi = 171.2;
  TString marr10[] = { "161", "163", "166", "168", "171", "173", "176", "178", "181", "183" };
  for (int k=0; k<10; k++) {
      if ( marr10[k] == mString ) mDigi = 161.2 + (2.5*k) ;
  }
  return mDigi ;
}

void MassAna::SetFitParameters( double mass, Double_t* para, int nPara, int NBTag, int rbin ) {

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
     // 10 GeV Bin - dM Constrain
     //Double_t kk0[12] = {  1.904, 12.051, 5.438,  4.202, 117.285, 21.724, 195.073, 39.88,   0.,    0.,  3.820,  0. };
     //Double_t kk1[12] = {  0.258,  0.933, 0.139, -0.010,   0.653,  0.115,   0.000,  0.00,   0.,    0., -0.011,  0. };
     // 10 GeV Bin - dM*hadW Constrain
     Double_t kk0[12] = {  62.167, 4.508, 2.097, -3.413, 161.463, 10.687, 183.757, 37.475,   0.,    0.,  0.919,  0. };
     Double_t kk1[12] = {  -0.099, 0.969, 0.142,  0.038,   0.205,  0.134,   0.000,  0.000,   0.,    0.,  0.007,  0. };
     for ( int i =0; i < 12; i++ ) {
         if ( cname == "had" ) {
            if ( nPara ==  9 && NBTag == -1 ) para[i] = kk0[i] + kk1[i]*mass ;
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

