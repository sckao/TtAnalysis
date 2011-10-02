#include "WAnalysis.h"
#include "WFormat.h"

WAnalysis::WAnalysis(){

  fitInput = new MassAnaInput();
  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "bThreshold", &bTh);
  fitInput->GetParameters( "n_btag", &n_btag);
  fitInput->GetParameters( "n_Jets", &n_Jets);
  fitInput->GetParameters( "PlotType", &plotType );

  fitInput->GetParameters( "PhaseSmear", &phaseSmear );
  smearing = ( phaseSmear == "ON" ) ? true : false ;

  fitInput->GetParameters( "InputMean", &inputMean );

  wmfitter  = new HadWMassFitter();
  ltmfitter = new LepTopMassFitter();
  pseudoExp = new PseudoExp();

  isFolderCreate = false;

}

WAnalysis::~WAnalysis(){

  delete fitInput ;
  delete wmfitter ;
  delete ltmfitter ; 	
  delete pseudoExp ;
}

void WAnalysis::CreateFolders(){

  if ( !isFolderCreate ) {
     TString theFolder = hfolder ;

     gSystem->mkdir( theFolder );
     gSystem->cd( theFolder );
     gSystem->mkdir( "M2M3" );
     gSystem->mkdir( "M3M3" );
     gSystem->mkdir( "M2M2t" );
     gSystem->mkdir( "EtaM2" );
     gSystem->mkdir( "EtaM3" );
     gSystem->mkdir( "YM2" );
     gSystem->mkdir( "YM3" );
     gSystem->mkdir( "M2MET" );
     gSystem->mkdir( "M2dF" );
     gSystem->mkdir( "M2M3BG" );
     
     gSystem->mkdir( "M2M3Lep" );
     gSystem->mkdir( "M3M3ReFit" );
     gSystem->mkdir( "M2tM3" );
     gSystem->mkdir( "PhiM3" );
     gSystem->mkdir( "M2M3ReFit" );

     gSystem->cd("../");
     isFolderCreate = true ;
  }

}
// For new SolTree
// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::HadTopFitter( string fileName, TString DrawOpt, bool doScale, bool isMCMatched  ){

  double scale = 1 ;
  string channelType ;
  if ( fileName.size() > 2) {
     vector<string>  channel;
     fitInput->GetParameters( "channel" , &channel );
     for (size_t i=0; i < channel.size(); i++ ) {
         if ( channel[i] == fileName.substr(0,2) )  channelType = channel[i] ;
         if ( fileName.substr(0,2) == "qc" )  channelType = channel[i] ;
         if ( fileName.substr(0,2) == "da"  )  channelType = "data" ;
     }
  }
  if ( channelType != "data" )   scale = fitInput->NormalizeComponents( channelType );

  // refit the solutions
  hadWBoson* wbh = new hadWBoson();
  // do the MC-Data normalization
  wmfitter->SetMCNormalization( doScale );

  if ( isMCMatched && channelType != "data" ) {
     wmfitter->MCSolution( fileName, wbh );
  } else {
     wmfitter->ReFitSolution( fileName, wbh, 4, scale, NULL, 0, smearing );
  }
  // retrieve the histograms
  //wbh->scale( scale );	
  vector<TH2D*> h2Ds = wbh->Output2D();

  // regular, more plots
  //M2M3Plotter( h2Ds, fileName, DrawOpt, isMCMatched );
  // For analysis note use
  AN_M2M3Plotter( h2Ds, fileName, DrawOpt, isMCMatched );

  delete wbh ;
  cout<<" DONE !! "<<endl;
}

void WAnalysis::LepTopFitter( string fileName, TString DrawOpt, bool isMCMatched ){
 
  string channelType ;
  if ( fileName.size() > 2) {
     vector<string>  channel;
     fitInput->GetParameters( "channel" , &channel );
     for (size_t i=0; i < channel.size(); i++ ) {
         if ( channel[i][0] != fileName[0] ) continue ;
         if ( channel[i][1] != fileName[1] ) continue ;
         channelType = channel[i] ;
     }
  }

  double scale = fitInput->NormalizeComponents( channelType );

  // refit the solutions
  lepWBoson* wbh = new lepWBoson();
  if ( isMCMatched ) {
     ltmfitter->LepTopMCSolution( fileName, wbh );
  } else {
     ltmfitter->ReFitLepTopSolution( fileName, wbh );
  }
  // retrieve the histograms
  wbh->scale( scale );	
  vector<TH2D*> h2Ds;
  wbh->Fill2DVec( h2Ds );

  LepTopPlotter( h2Ds, fileName, DrawOpt, isMCMatched );

  delete wbh ;
  cout<<" DONE !!! "<<endl;
}

void WAnalysis::M2M3_1DPlots( string dataName, vector<string>& mcFiles, bool doScale ){

  TString theFolder = hfolder ;
  TString theSubFolder = ( doScale == false ) ?  "hTopo/" : "hTopo_Scale/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  gStyle->SetOptStat("");
  gStyle->SetLabelSize( 0.04, "X");
  gStyle->SetLabelSize( 0.04, "Y");
  gStyle->SetLabelSize( 0.04, "Z");
  gStyle->SetOptTitle(0);


  wmfitter->SetMCNormalization( doScale );
  double scaleV[7] = {1} ;
  scaleV[0]  = fitInput->NormalizeComponents( "tt" );
  scaleV[1]  = fitInput->NormalizeComponents( "wj" );
  scaleV[2]  = fitInput->NormalizeComponents( "zj" );
  scaleV[3]  = fitInput->NormalizeComponents( "tq" );
  scaleV[4]  = fitInput->NormalizeComponents( "tw" );
  scaleV[5]  = fitInput->NormalizeComponents( "ww" );
  scaleV[6]  = fitInput->NormalizeComponents( "qcd" );

  hTopo* hdt = new hTopo();
  wmfitter->ReFitSolution( dataName, hdt, 4, 1, NULL, 0 );
  vector<TH1D*> dtV ;
  hdt->Fill1DVec( dtV );

  TString hNames[20] = { "Mt",   "hadM2",   "hadM3",   "lepM3",   "Mt_1", "hadM2_1", "hadM3_1", "lepM3_1", 
                         "Mt_2", "hadM2_2", "hadM3_2", "lepM3_2", "Mt_3", "hadM2_3", "hadM3_3", "lepM3_3",
                         "Mt_4", "hadM2_4", "hadM3_4", "lepM3_4"  };
  vector<TH1D*> mcV ;
  vector< vector<TH1D*> > mcVV(20) ;
  for ( int i =0 ; i < 7 ; i++) {
       hTopo* hmc = new hTopo( mcFiles[i].substr(0,2) );
       wmfitter->ReFitSolution( mcFiles[i], hmc, 4, scaleV[i], NULL, 0, smearing );
       hmc->Fill1DVec( mcV );
       for ( int k=0; k<20; k++) {
           mcVV[k].push_back( (TH1D*)mcV[k]->Clone( mcFiles[i].substr(0,2).c_str() ) ) ;
       }
       mcV.clear();
       delete hmc ;
  }
  for ( int k =0 ; k < 20; k++) {

      vector<TH1D*> m2hV = mcVV[k];
      THStack* hStk = new THStack("hStk", "M2M3" );
      m2hV[0]->SetFillColor( kRed+1 ) ;
      m2hV[1]->SetFillColor( kGreen ) ;
      m2hV[2]->SetFillColor( kAzure-2 ) ;
      m2hV[3]->SetFillColor( kMagenta+2 ) ;
      m2hV[4]->SetFillColor( kMagenta ) ;
      m2hV[5]->SetFillColor( kWhite ) ;
      m2hV[6]->SetFillColor( kYellow ) ;

      TH1D* All_MC  = (TH1D*) m2hV[0]->Clone("combined MC");
      All_MC->Add( m2hV[1] );
      All_MC->Add( m2hV[2] );
      All_MC->Add( m2hV[3] );
      All_MC->Add( m2hV[4] );
      All_MC->Add( m2hV[5] );
      All_MC->Add( m2hV[6] );

      TLegend *leg = new TLegend(.55, .58, .84, .9 );
      TString sName[7] = { "Ttbar", "WJets", "ZJets", "single Top (t)", "single Top (tW)", "WW", "QCD" } ;
      for ( int i = 6 ; i >= 0  ; i--) {
          int ii = 6 - i ;
          hStk->Add( m2hV[i] ) ; 
          leg->AddEntry( m2hV[ii],  sName[ii],  "F");
      }

      TCanvas* c1 = new TCanvas("c1","", 800, 700);
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd();

      double hMax = 1.45* dtV[k]->GetBinContent( dtV[k]->GetMaximumBin() ) ;
      dtV[k]->SetMaximum( hMax );
      dtV[k]->SetMarkerSize(1);
      dtV[k]->SetMarkerStyle(21);
      dtV[k]->SetMinimum( 0. );
      dtV[k]->Draw("PE");
      c1->Update();

      hStk->Draw("same");
      c1->Update();
  
      dtV[k]->Draw("PE SAME");
      c1->Update();

      leg->Draw("same");
      c1->Update();

      //KS Test
      /*
      double ksP = All_MC->KolmogorovTest( dtV[k]  ) ;
      double ksD = All_MC->KolmogorovTest( dtV[k], "M" ) ;
      cout<<" KSP = "<< ksP <<"  KS = "<<ksD <<endl ;
      ostringstream ksStr ;
      ksStr << "KS = ";
      ksStr << setprecision(4) << ksP ;
      //ksStr << " D = " ;
      //ksStr << ksD ;
      TString ksReport = ksStr.str() ;
      TPaveText *pvtxt = new TPaveText( 0.35, 0.82, .55, .9,"NDC" );
      //pvtxt->SetFillColor(0);
      pvtxt->SetTextSize(0.03) ;
      pvtxt->AddText( ksReport ) ;

      pvtxt->Draw("same");
      c1->Update();
      */

      TString plotname1 = hfolder + theSubFolder+ hNames[k] + "."+plotType ;
      c1->Print( plotname1 );

      delete c1 ;
      delete leg ;
      delete hStk ;
  }
  delete hdt ;

}

void WAnalysis::Had_SBRatio( vector<string>& mcFile ){

  TString filepath = hfolder+"SBRatio.log" ;
  FILE* logfile = fopen( filepath,"a");
  fprintf(logfile," M2M3   dM3   S/B   SS   N_tt  N_bg  N_vj  N_qcd  Eff_tt  Eff_bg  Eff_vj  Eff_qcd  \n" ) ;  

  vector<string> File4J ;
  fitInput->GetParameters( "FakeData", &File4J );

  vector<double> nEvts ;
  fitInput->GetParameters( "nEvents" , &nEvts );

  double scaleTt  = fitInput->NormalizeComponents( "tt" );
  double scaleWJ  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ  = fitInput->NormalizeComponents( "zj" );
  double scaleTq  = fitInput->NormalizeComponents( "tq" );
  double scaleTW  = fitInput->NormalizeComponents( "tw" );
  double scaleWW  = fitInput->NormalizeComponents( "ww" );
  double scaleQCD = fitInput->NormalizeComponents( "qcd" );

  // normalize the MC for the counter 
  string mcNormalization ;
  fitInput->GetParameters( "MCNormalization", &mcNormalization );
  bool mcNorm = ( mcNormalization == "YES" ) ? true : false ;
  wmfitter->SetMCNormalization( mcNorm ) ;

  TH2D* hSB = new TH2D("hSB", " S/B with cuts  dM3(X), M2M3(Y) ",     8, 2.5, 42.5,  8, 2.5, 42.5 );
  TH2D* hPR = new TH2D("hPR", " Purity with cuts  dM3(X), M2M3(Y) ",  8, 2.5, 42.5,  8, 2.5, 42.5 );
  TH2D* hCE = new TH2D("hCE", " Cut Eff for signal dM3(X), M2M3(Y) ", 8, 2.5, 42.5,  8, 2.5, 42.5 );
 
  double nAll[4] = {0.};
  double cEff[4] = {0.};
  cout<<" ---  Measuring SB Ratio( tt / wj+Zj+qcd ) --- "<<endl;
  for ( int j = 0; j< 8; j++ ) {
      double m3w[7] = {0.} ;
      double Effa[7] = {0.} ;
      double SBRa[7] = {0.} ;
      double PURa[7] = {0.} ;
      for ( int k= 0; k< 8; k++ ) {

         // reset the M2M3 windows
	 double window  = 5. + j*5 ;
	 double window2 = 5. + k*5 ;
         double dM3     = 999 ;
	 //double window2  = window*1.5 ;
         //double dM3     = 0 + k*5 ;

	 double m2L = 80. - window ;
	 double m2H = 80. + window ;
	 double m3L = 170. - window2 ;
	 double m3H = 170. + window2 ;
         
	 if ( j == 0 && k == 0 ) { wmfitter->ResetCuts(  0., 999.,  0., 999.,      0., 999. );
                                             window = 999. ;
                                             dM3    = 999. ;    }
         /*
         if ( j != 0 && k == 0 ) { wmfitter->ResetCuts( m2L,  m2H,  0,  999,      0., dM3 );
                                             window2  = 999. ;    }
         if ( j == 0 && k != 0 ) { wmfitter->ResetCuts(  0., 999., m3L, m3H,      0., dM3 );
                                             window   = 999. ;    }
         if ( j != 0 && k == 0 ) { wmfitter->ResetCuts( m2L,  m2H, m3L,  m3H,      0., 999. );
                                             dM3    = 999. ;    }
         if ( j == 0 && k != 0 ) { wmfitter->ResetCuts(  0., 999.,  0., 999.,      0.,  dM3 );
                                             window = 999. ;    }
         */
         if ( j == 0 && k != 0 ) continue;
         if ( j != 0 && k == 0 ) continue;
        
         if ( j != 0 && k != 0 ) wmfitter->ResetCuts( m2L, m2H, m3L, m3H, 0, dM3 );

         ACounter* sg4j = new ACounter();
	 wmfitter->ReFitSolution( mcFile[0], sg4j, 4, scaleTt,  NULL, 0, smearing );
	 vector<double> nSg4J = sg4j->Output();

	 ACounter* bg4j = new ACounter();
	 wmfitter->ReFitSolution( mcFile[1], bg4j, 4, scaleWJ,  NULL, 0, smearing );
	 wmfitter->ReFitSolution( mcFile[2], bg4j, 4, scaleZJ,  NULL, 0, smearing );
	 wmfitter->ReFitSolution( mcFile[3], bg4j, 4, scaleTq,  NULL, 0, smearing );
	 wmfitter->ReFitSolution( mcFile[4], bg4j, 4, scaleTW,  NULL, 0, smearing );
	 wmfitter->ReFitSolution( mcFile[5], bg4j, 4, scaleWW,  NULL, 0, smearing );
	 vector<double> nWj4J = bg4j->Output();
	 wmfitter->ReFitSolution( mcFile[6], bg4j, 4, scaleQCD, NULL, 0, smearing );
         vector<double> nBg4J = bg4j->Output();

         // get the denumerator and cut efficiency for each channel
         if ( j == 0 && k == 0 ) {
            nAll[0] = nSg4J[0] ;
	    nAll[1] = nWj4J[0] ;
	    nAll[2] = nBg4J[0] - nWj4J[0];
	    nAll[3] = nBg4J[0] ;
         }
	 cEff[0] = nSg4J[0]/nAll[0] ;                 // signal eff
	 cEff[1] = nWj4J[0]/nAll[1] ;                 // vjets  eff
	 cEff[2] = (nBg4J[0] - nWj4J[0]) /nAll[2] ;   // qcd    eff
	 cEff[3] = nBg4J[0]/nAll[3] ;                 // all bg eff

         double ssR = nSg4J[0] / (nSg4J[0] + nBg4J[0]) ; // purity
         double sbR = nSg4J[0] / nBg4J[0] ;              // S/B ratio

         cout<<" Tt: "<< nSg4J[0] <<" Bg: "<< nBg4J[0] <<" SB: "<< nSg4J[0] / nBg4J[0] ;
         cout<<" SS: "<< nSg4J[0] / sqrt( nBg4J[0] + nSg4J[0] ) <<" CutEff ="<< cEff[0]<<" / "<<cEff[1]<<endl;

         if ( j != 0 || k !=0 ) {
            hSB->Fill( window2, window, sbR );
            hPR->Fill( window2, window, ssR );
            hCE->Fill( window2, window, cEff[0] );
         }
         if ( k != 0 ) {
            Effa[k-1] = cEff[0] ;
            PURa[k-1] = ssR ;
            SBRa[k-1] = sbR ;
            m3w[k-1]  = window2 ;
         }
         if (  cEff[0] > 0. ) {
            fprintf(logfile," %.1f   %1.f  %.2f %.2f   %.1f    %.1f      %.1f       %.1f   %.2f   %.2f   %.2f   %.2f \n",  
                           window, window2, sbR, ssR, nSg4J[0], nBg4J[0], nWj4J[0], (nBg4J[0] - nWj4J[0]), cEff[0], cEff[3], cEff[1], cEff[2]);
         }

          delete sg4j;
          delete bg4j;
      }
  }
  cout<<" SB Ration Measured !! "<<endl;

  fclose(logfile);

  gStyle->SetOptStat("");
  //gStyle->SetNumberContours(10);
  //gStyle->SetStatX(0.95);
  //gStyle->SetStatY(0.99);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  gStyle->SetPaintTextFormat("2.2f");
  c1->cd();
  hSB->SetMarkerSize(2.2);
  hSB->Draw("TEXT");
  c1->Update();

  TString plotname1 = hfolder + "SBRatio."+plotType ;
  c1->Print( plotname1 );

  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  gStyle->SetPaintTextFormat("2.2f");
  c2->cd();
  hCE->SetMarkerSize(2.2);
  hCE->Draw("TEXT");
  c2->Update();

  TString plotname2 = hfolder + "CutEff."+plotType ;
  c2->Print( plotname2 );

  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  gStyle->SetPaintTextFormat("2.2f");
  c3->cd();
  hPR->SetMarkerSize(2.2);
  hPR->Draw("TEXT");
  c3->Update();

  TString plotname3 = hfolder + "Purity."+plotType ;
  c3->Print( plotname3 );

  delete c1;
  delete c2;
  delete c3;
  delete hSB ;
  delete hPR ;
  delete hCE ;
}

void WAnalysis::SBPlotter(){

  FILE* pfile ;
  TString fpath = hfolder+"SBRatio1.log" ;
  pfile = fopen( fpath,"r");

  const int sz = 7 ;

  gStyle->SetOptTitle(0);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleOffset(1.5, "Y") ;
  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  TLegend *leg1 = new TLegend(.75, .57, .9, .9 );

  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  TLegend *leg2 = new TLegend(.73, .12, .88, .45 );

  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  TLegend *leg3 = new TLegend(.75, .57, .9, .9 );

  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  TLegend *leg4 = new TLegend(.18, .57, .33, .9 );

  TCanvas* c5 = new TCanvas("c5","", 800, 600);
  c5->SetGrid();
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  TLegend *leg5 = new TLegend(.73, .12, .88, .45 );

  float f0 ;
  float m2w[sz] ;
  float m3w[sz] ;
  float sbR[sz] ;
  float eff[sz] ;
  float fkR[sz] ;
  float pur[sz] ;
  float Ntt[sz] ;
  float Nbg[sz] ;
  float sgR[sz] ;
  string m2strV[sz] = { "10", "15", "20", "25", "30", "35", "40"} ;
  int colorcode[sz] = { 2, 4, 3, 7, 8, 6, 1 } ;
  float m3w1[sz] = { 10, 15, 20, 25, 30, 35, 40 } ;

  for (int k=0 ; k< 7; k++) {               // M2 width
      for ( int i=0 ; i< 7; i++) {          // test value
          for ( int j=0 ; j< 12; j++) {     // reading entries
              int fo = fscanf(pfile, "%f", &f0 );
	      cout<<" f"<<j<<" = "<< f0;
	      if ( fo != 1 ) cout<<" Logfile Reading Error "<<endl;
	      if ( j == 0 && i == 0 ) m2w[k] = f0 ;
	      if ( j == 1 ) m3w[i] = f0 ;
	      if ( j == 2 ) sbR[i] = f0 ;
	      if ( j == 4 ) Ntt[i] = f0 ;
	      if ( j == 5 ) Nbg[i] = f0 ;
	      if ( j == 3 ) pur[i] = f0 ;
	      if ( j == 8 ) eff[i] = f0 ; // signal efficiency
	      if ( j == 9 ) fkR[i] = f0 ; // background efficiency
          }
          sgR[i] = Ntt[i] / sqrt( Ntt[i] + Nbg[i] ) ;
          cout<< endl ;
      }
      string m2str = m2strV[k]  ;
      TString xtitle = "Width of M3 window";
      //TString xtitle = " Width of #Delta M3";

      c1->cd();
      TGraph* gSBR = new TGraph(sz, m3w1, sbR );
      gSBR->SetMarkerStyle(20);
      gSBR->SetMarkerSize(1);
      gSBR->SetMarkerColor( colorcode[k] );
      gSBR->SetLineColor( colorcode[k] );
      //if ( k==0 ) gSBR->SetMaximum(3.0);
      if ( k==0 ) gSBR->SetMaximum(2.0);
      if ( k==0 ) gSBR->SetMinimum(1.0);
      gSBR->GetXaxis()->SetTitle( xtitle );
      gSBR->GetYaxis()->SetTitle(" S/B ");
      if ( k == 0 ) gSBR->Draw("ACP");
      if ( k  > 0 ) gSBR->Draw("CPSAME");
      leg1->AddEntry(gSBR,  m2str.c_str(), "P");
      c1->Update();

      c2->cd();
      TGraph* gEff = new TGraph(sz, m3w1, eff );
      gEff->SetMarkerStyle(21);
      gEff->SetMarkerSize(1);
      gEff->SetMarkerColor( colorcode[k] );
      gEff->SetLineColor( colorcode[k] );
      //if ( k==0 ) gEff->SetMaximum(0.5);
      //if ( k==0 ) gEff->SetMinimum(0.0);
      if ( k==0 ) gEff->SetMaximum(0.9);
      if ( k==0 ) gEff->SetMinimum(0.2);
      gEff->GetXaxis()->SetTitle( xtitle );
      gEff->GetYaxis()->SetTitle("  Efficiency  ");
      if ( k == 0 ) gEff->Draw("ACP");
      if ( k  > 0 ) gEff->Draw("CPSAME");
      leg2->AddEntry(gEff,  m2str.c_str(), "P");
      c2->Update();

      c3->cd();
      TGraph* gPur = new TGraph(sz, m3w1, pur );
      gPur->SetMarkerStyle(22);
      gPur->SetMarkerSize(1);
      gPur->SetMarkerColor( colorcode[k] );
      gPur->SetLineColor( colorcode[k] );
      //if ( k==0 ) gPur->SetMaximum(0.8);
      if ( k==0 ) gPur->SetMaximum(0.68);
      if ( k==0 ) gPur->SetMinimum(0.5);
      gPur->GetXaxis()->SetTitle( xtitle );
      gPur->GetYaxis()->SetTitle("  Purity  #frac{S}{S+B} ");
      if ( k == 0 ) gPur->Draw("ACP");
      if ( k  > 0 ) gPur->Draw("CPSAME");
      leg3->AddEntry(gPur,  m2str.c_str(), "P");
      c3->Update();

      c4->cd();
      TGraph* gFkr = new TGraph(sz, m3w1, fkR );
      gFkr->SetMarkerStyle(21);
      gFkr->SetMarkerSize(1);
      gFkr->SetMarkerColor( colorcode[k] );
      gFkr->SetLineColor( colorcode[k] );
      if ( k==0 ) gFkr->SetMaximum(0.9);
      if ( k==0 ) gFkr->SetMinimum(0.1);
      //if ( k==0 ) gFkr->SetMaximum(0.4);
      //if ( k==0 ) gFkr->SetMinimum(0.0);
      gFkr->GetXaxis()->SetTitle( xtitle );
      gFkr->GetYaxis()->SetTitle("  Background Efficiency  ");
      if ( k == 0 ) gFkr->Draw("ACP");
      if ( k  > 0 ) gFkr->Draw("CPSAME");
      leg4->AddEntry(gFkr,  m2str.c_str(), "P");
      c4->Update();

      c5->cd();
      TGraph* gSGR = new TGraph(sz, m3w1, sgR );
      gSGR->SetMarkerStyle(20);
      gSGR->SetMarkerSize(1);
      gSGR->SetMarkerColor( colorcode[k] );
      gSGR->SetLineColor( colorcode[k] );
      if ( k==0 ) gSGR->SetMaximum(11.0);
      if ( k==0 ) gSGR->SetMinimum(6.0);
      //if ( k==0 ) gSGR->SetMaximum(9.0);
      //if ( k==0 ) gSGR->SetMinimum(4.0);
      gSGR->GetXaxis()->SetTitle( xtitle );
      gSGR->GetYaxis()->SetTitle(" #frac{S}{ #sqrt{S+B} }");
      if ( k == 0 ) gSGR->Draw("ACP");
      if ( k  > 0 ) gSGR->Draw("CPSAME");
      leg5->AddEntry(gSGR,  m2str.c_str(), "P");
      c5->Update();

  }
  c1->cd();
  leg1->Draw("SAME");
  c1->Update();
  c2->cd();
  leg2->Draw("SAME");
  c2->Update();
  c3->cd();
  leg3->Draw("SAME");
  c3->Update();
  c4->cd();
  leg4->Draw("SAME");
  c4->Update();
  c5->cd();
  leg5->Draw("SAME");
  c5->Update();

  TString plotname1 = hfolder + "SBRatio."+plotType ;
  c1->Print( plotname1 );

  TString plotname2 = hfolder + "TtEff."+plotType ;
  c2->Print( plotname2 );

  TString plotname3 = hfolder + "Purity."+plotType ;
  c3->Print( plotname3 );

  TString plotname4 = hfolder + "BgEff."+plotType ;
  c4->Print( plotname4 );

  TString plotname5 = hfolder + "Significancy."+plotType ;
  c5->Print( plotname5 );

  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
  fclose(pfile);

}

void WAnalysis::SBCEPlotter(){

  TH1D* hCE0 = new TH1D("hCE0", " Cut Eff for M2M3 ", 10, 7.5, 57.5 );
  TH1D* hSB0 = new TH1D("hSB0", " S/B     for M2M3 ", 10, 7.5, 57.5 );
  TH1D* hCE1 = new TH1D("hCE1", " Cut Eff for M2M3 , dM = 30 ", 10, 2.5, 52.5 );
  TH1D* hSB1 = new TH1D("hSB1", " S/B     for M2M3 , dM = 30 ", 10, 2.5, 52.5 );
  TH1D* hCE2 = new TH1D("hCE2", " Cut Eff for M2M3 , dM = 20 ", 10, 2.5, 52.5 );
  TH1D* hSB2 = new TH1D("hSB2", " S/B     for M2M3 , dM = 20 ", 10, 2.5, 52.5 );

  // dM3 = 999
  double sbR[8]  = { 1.97, 1.65, 1.42, 1.31, 1.26, 1.20, 1.17, 1.15 };   
  double cEff[8] = { 0.46, 0.65, 0.76, 0.82, 0.86, 0.89, 0.91, 0.92 };   
  // dM3 = 30
  double sbR1[8]  = { 2.43, 2.02, 1.74, 1.54, 1.48, 1.40, 1.35, 1.32 };   
  double cEff1[8] = { 0.36, 0.51, 0.60, 0.66, 0.69, 0.72, 0.74, 0.75 };   
  // dM3 = 20
  double sbR2[8]  = { 2.83, 2.25, 1.88, 1.67, 1.57, 1.50, 1.43, 1.40 };   
  double cEff2[8] = { 0.31, 0.44, 0.53, 0.58, 0.62, 0.64, 0.66, 0.68 };   

  for ( int i=0; i< 8; i++){
      double window  = 15. + i*5 ;
      double sR0 = 1. / ( (1./sbR[i] ) + 1.) ;
      double sR1 = 1. / ( (1./sbR1[i]) + 1.) ;
      double sR2 = 1. / ( (1./sbR2[i]) + 1.) ;

      hSB0->Fill( window, sR0 );
      hCE0->Fill( window, cEff[i] );
      hSB1->Fill( window, sR1 );
      hCE1->Fill( window, cEff1[i] );
      hSB2->Fill( window, sR2 );
      hCE2->Fill( window, cEff2[i] );
  }

  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->cd();
  hSB0->SetAxisRange(0.1, 1., "Y");
  hSB0->SetMarkerSize(1.5);
  hSB0->SetMarkerStyle(20);
  hSB0->SetXTitle(" Width of M2 window, #times 1.5 = Width of M3 window");
  hSB0->SetYTitle(" Purity ( #frac{S}{S+B} ) ");
  hSB0->SetTitleOffset(1.5, "Y") ;
  hSB0->Draw("P");
  c1->Update();

  hSB1->SetMarkerSize(1.5);
  hSB1->SetMarkerStyle(21);
  hSB1->Draw("sameP");

  hSB2->SetMarkerSize(1.5);
  hSB2->SetMarkerStyle(22);
  hSB2->Draw("sameP");
  c1->Update();

  //float rMax   = 1.1 * hCE0->GetMaximum();
  float rMax   = 1. ;
  float rScale = gPad->GetUymax() / rMax ;
  hCE0->Scale(rScale)   ;
  hCE0->SetMarkerColor(2);
  hCE0->SetMarkerSize(2);
  hCE0->SetMarkerStyle(24);
  hCE0->Draw("sameP");
  c1->Update();

  hCE1->Scale(rScale)   ;
  hCE1->SetMarkerColor(2);
  hCE1->SetMarkerSize(2);
  hCE1->SetMarkerStyle(25);
  hCE1->Draw("sameP");
  c1->Update();

  hCE2->Scale(rScale)   ;
  hCE2->SetMarkerColor(2);
  hCE2->SetMarkerSize(2);
  hCE2->SetMarkerStyle(26);
  hCE2->Draw("sameP");
  c1->Update();

  TLegend *leg = new TLegend(.50, .12, .75, .40 );
  leg->AddEntry(hSB0,  "Purity( No dM cut)", "P");
  leg->AddEntry(hCE0,  "Efficiency( No dM cut)", "P");
  leg->AddEntry(hSB1,  "Purity( dM < 30)",   "P");
  leg->AddEntry(hCE1,  "Efficiency( dM < 30)",   "P");
  leg->AddEntry(hSB2,  "Purity( dM < 20)",   "P");
  leg->AddEntry(hCE2,  "Efficiency( dM < 20)",   "P");
  leg->Draw();

  TGaxis *axis = new TGaxis( gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(), 0.1, rMax, 510,"+L");
  axis->SetLineColor(2);
  axis->SetLabelColor(2);
  axis->SetTitle("Signal Efficiency") ;
  axis->SetTitleColor(2) ;
  axis->Draw();

  TString plotname3 = hfolder + "SB_Eff."+plotType ;
  c1->Print( plotname3 );

  delete c1;
  delete leg;
  delete hSB0;
  delete hCE0;
  delete hSB1;
  delete hCE1;
}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::MixAll( vector<string>& flist, TString DrawOpt ){
 
  cout<<" Mix All  "<<endl;

  // combined plotting
  /*
  gStyle->SetOptStat("neirom");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.90);
  gStyle->SetStatTextColor(1);
  */
  gStyle->SetPalette(1);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  //c1->Divide(2,2);

  hadWBoson* allh  = new hadWBoson();

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  wmfitter->ReFitSolution( flist[1], allh, 4, scale1, NULL, 0, smearing );
  wmfitter->ReFitSolution( flist[2], allh, 4, scale2, NULL, 0, smearing );
  wmfitter->ReFitSolution( flist[3], allh, 4, scale3, NULL, 0, smearing );
  wmfitter->ReFitSolution( flist[4], allh, 4, scale4, NULL, 0, smearing );
  wmfitter->ReFitSolution( flist[5], allh, 4, scale5, NULL, 0, smearing );
  wmfitter->ReFitSolution( flist[6], allh, 4, scale6, NULL, 0, smearing );

  vector<TH2D*> hbgs ;
  allh->Fill2DVec( hbgs );
  AN_M2M3Plotter( hbgs, "Bg", DrawOpt ) ;

  wmfitter->ReFitSolution( flist[0], allh, 4, scale0, NULL, 0, smearing );

  vector<TH2D*> h2Ds ;
  allh->Fill2DVec( h2Ds );
  AN_M2M3Plotter( h2Ds, "ALL", DrawOpt ) ;

  /*
  c1->cd(2);
  TH1D* bgM2M3_py = hbgs[0]->ProjectionY("bgM2M3_py");
  bgM2M3_py->SetFillColor(2);
  bgM2M3_py->SetMaximum(3.5);
  bgM2M3_py->Draw();
  c1->Update();
  c1->cd(3);
  TH1D* bgM2M3_px = hbgs[0]->ProjectionX("bgM2M3_px");
  bgM2M3_px->SetFillColor(2);
  bgM2M3_px->SetMaximum(3.5);
  bgM2M3_px->Draw();
  c1->Update();
  double scale0 = fitInput->NormalizeComponents( "tt" );
  wmfitter->ReFitSolution( flist[0], allh, 4, scale0, NULL, 0, smearing );
  */
 
  //M2M3Plotter( h2Ds, "ALL", DrawOpt ) ;

  /*
  c1->cd(1);
  gStyle->SetStatX(0.30);
  gStyle->SetNumberContours(5);
  h2Ds[0]->Draw( DrawOpt );
  c1->Update();

  c1->cd(2);
  gStyle->SetStatX(0.90);
  TH1D* hM2M3_py = h2Ds[0]->ProjectionY("hM2M3_py");
  hM2M3_py->Draw("sames");
  c1->Update();

  c1->cd(3);
  TH1D* hM2M3_px = h2Ds[0]->ProjectionX("hM2M3_px");
  hM2M3_px->Draw("sames");
  c1->Update();

  c1->cd(4);
  h2Ds[0]->Draw( "LEGO2Z" );
  c1->Update();

  TString plotname1 = hfolder + "M2M3/ALL_M2M3."+plotType ;
  c1->Print( plotname1 );
  */

  delete allh;

}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::EnsembleTest( int randomSeed, TString DrawOpt ){
 
  vector<string> flist;
  fitInput->GetParameters( "FakeData", &flist );

  TString treeName = "muJets";

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  hadWBoson* sg = new hadWBoson();

  // ttbar
  vector<int> sglist = pseudoExp->GetEnsemble( flist[0], treeName, inputMean[0], randomSeed );
  wmfitter->ReFitSolution( flist[0], sg, 4, 1, &sglist );
  cout<<" tt test done "<<endl;

  // Making signal plots
  vector<TH2D*> hs0 = sg->Output2D();
  M2M3Plotter( hs0, "pseudoSG", DrawOpt );

  hadWBoson* bg = new hadWBoson();
  // w+jets
  vector<int> bg1list = pseudoExp->GetEnsemble( flist[1], treeName, inputMean[1], randomSeed );
  wmfitter->ReFitSolution( flist[1], bg , 4, 1, &bg1list, 0, true );
  cout<<" wj test done "<<endl;

  // QCD
  vector<int> bg2list = pseudoExp->GetEnsemble( flist[2], treeName, inputMean[2], randomSeed );
  wmfitter->ReFitSolution( flist[2], bg , 4, 1, &bg2list, 0, true );
  cout<<" qcd test done "<<endl;

  // retrieve the histograms
  vector<TH2D*> hBgs = bg->Output2D();
  M2M3Plotter( hBgs, "pseudoBG", DrawOpt );

  for ( size_t i =0; i < hs0.size(); i++ ) {
      hs0[i]->Add( hBgs[i] );
  }
  M2M3Plotter( hs0, "pseudoExp", DrawOpt );

  delete sg ;
  delete bg ;
  cout<<" YA !!! "<<endl;
}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
/*
void WAnalysis::HadTEnsembleTest( int randomSeed, TString DrawOpt ){
 
  vector<string> flist;
  fitInput->GetParameters( "FakeData", &flist );

  TString treeName = "muJets";

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  // ttbar
  hadWBoson* sg = new hadWBoson();
  double inputMean_tt = ( n_Jets == 4 ) ? 9.28 : 8.701 ;
  vector<int> sglist = pseudoExp->GetEnsemble( flist[0], treeName, inputMean_tt, randomSeed );
  wmfitter->ReFitLepTopSolution( flist[0], sg, 1, &sglist );
  cout<<" tt test done "<<endl;

  // Making signal plots
  vector<TH2D*> h_sg;
  sg->Fill2DVec( h_sg );
  M2M3Plotter( h_sg, "pseudoSG", DrawOpt );

  // w+jets
  hadWBoson* bg = new hadWBoson();
  double inputMean_wj = ( n_Jets == 4 ) ? 2.458 : 92.0 ;
  vector<int> bg1list = pseudoExp->GetEnsemble( flist[1], treeName, inputMean_wj, randomSeed );
  wmfitter->ReFitLepTopSolution( flist[1], bg, 1, &bg1list );
  cout<<" wj test done "<<endl;

  // QCD
  double inputMean_qcd = ( n_Jets == 4 ) ?  3.31 : 131.352 ;
  vector<int> bg2list = pseudoExp->GetEnsemble( flist[2], treeName, inputMean_qcd, randomSeed );
  wmfitter->ReFitLepTopSolution( flist[2], bg, 1, &bg2list );
  cout<<" qcd test done "<<endl;

  // retrieve the histograms
  vector<TH2D*> h_bg;
  bg->Fill2DVec( h_bg );
  M2M3Plotter( h_bg, "pseudoBG", DrawOpt );

  for ( size_t i =0; i < h_sg.size(); i++ ) {
      h_sg[i]->Add( h_bg[i] );
  }
  M2M3Plotter( h_sg, "pseudoExp", DrawOpt );

  delete sg ;
  delete bg ;
  cout<<" YA !!! "<<endl;

}
*/

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::LepTEnsembleTest( int randomSeed, TString DrawOpt ){
 
  vector<string> flist;
  fitInput->GetParameters( "FakeData", &flist );

  TString treeName = "muJets";

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  // ttbar
  lepWBoson* sg = new lepWBoson();
  double inputMean_tt = ( n_Jets == 4 ) ? 9.28 : 8.701 ;
  vector<int> sglist = pseudoExp->GetEnsemble( flist[0], treeName, inputMean_tt, randomSeed );
  ltmfitter->ReFitLepTopSolution( flist[0], sg, 1, &sglist );
  cout<<" tt test done "<<endl;

  // Making signal plots
  vector<TH2D*> h_sg;
  sg->Fill2DVec( h_sg );
  LepTopPlotter( h_sg, "pseudoSG", DrawOpt );

  // w+jets
  lepWBoson* bg = new lepWBoson();
  double inputMean_wj = ( n_Jets == 4 ) ? 2.458 : 92.0 ;
  vector<int> bg1list = pseudoExp->GetEnsemble( flist[1], treeName, inputMean_wj, randomSeed );
  ltmfitter->ReFitLepTopSolution( flist[1], bg, 1, &bg1list );
  cout<<" wj test done "<<endl;

  // QCD
  double inputMean_qcd = ( n_Jets == 4 ) ?  3.31 : 131.352 ;
  vector<int> bg2list = pseudoExp->GetEnsemble( flist[2], treeName, inputMean_qcd, randomSeed );
  ltmfitter->ReFitLepTopSolution( flist[2], bg, 1, &bg2list );
  cout<<" qcd test done "<<endl;

  // retrieve the histograms
  vector<TH2D*> h_bg;
  bg->Fill2DVec( h_bg );
  LepTopPlotter( h_bg, "pseudoBG", DrawOpt );

  for ( size_t i =0; i < h_sg.size(); i++ ) {
      h_sg[i]->Add( h_bg[i] );
  }
  LepTopPlotter( h_sg, "pseudoExp", DrawOpt );

  delete sg ;
  delete bg ;
  cout<<" YA !!! "<<endl;

}

void WAnalysis::M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt, bool isMCMatched ){

  TString dNames[10] = { "M2M3",  "M3M3",  "EtaM2", "EtaM3", "YM2",
                          "YM3", "M2MET",   "M3Pt", "METMt", "M2M3BG" };

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );

  
  gStyle->SetOptStat("neirom");
  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(1);
  
  gStyle->SetPalette(1);
  gStyle->SetLabelSize( 0.06, "X");
  gStyle->SetLabelSize( 0.06, "Y");
  gStyle->SetLabelSize( 0.06, "Z");

  for (size_t i=0; i< h2Ds.size(); i++ ) {

      gSystem->cd( theFolder );
      gSystem->mkdir( dNames[i] );
      gSystem->cd("../");

      TCanvas* c1 = new TCanvas("c1","", 1000, 800);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->Divide(2,2);
      c1->cd(1);
      gStyle->SetStatX(0.30);
      gStyle->SetNumberContours(5);
      h2Ds[i]->Draw( DrawOpt );
      c1->Update();

      c1->cd(2);
      gStyle->SetStatX(0.90);
      TH1D* h2Di_py = h2Ds[i]->ProjectionY("h2Di_py");
      h2Di_py->Draw();
      c1->Update();

      c1->cd(3);
      TH1D* h2Di_px = h2Ds[i]->ProjectionX("h2Di_px");
      h2Di_px->Draw();
      c1->Update();

      c1->cd(4);
      h2Ds[i]->Draw( "LEGO2Z" );
      c1->Update();

      TString plotname1 = hfolder + dNames[i]+"/"+  fileName +"_"+ dNames[i] + "."+plotType ;
      if ( isMCMatched ) plotname1 = hfolder + dNames[i] + "/" + fileName + "_"+ dNames[i]+"_MC."+plotType ;
      c1->Print( plotname1 );

      delete c1;
      delete h2Di_py;
      delete h2Di_px;
  }

  cout<<" FINISHED  "<<endl;
}

void WAnalysis::AN_M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt, bool isMCMatched ){

  TString dNames[10] = { "M2M3",  "M3M3",  "EtaM2", "EtaM3", "YM2",
                          "YM3", "M2MET",   "M3Pt", "METMt", "M2M3BG" };

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );

  
  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);

  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPalette(1);
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");

  for (size_t i=0; i< h2Ds.size(); i++ ) {

      if ( i > 1 ) continue ;
      gSystem->cd( theFolder );
      gSystem->mkdir( dNames[i] );
      gSystem->cd("../");

      TCanvas* c1 = new TCanvas("c1","", 1000, 800);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd(1);
      gStyle->SetNumberContours(5);
      h2Ds[i]->SetTitleSize(0.05, "X") ;
      h2Ds[i]->SetTitleSize(0.05, "Y") ;
      h2Ds[i]->SetTitleOffset(1.5, "X") ;
      h2Ds[i]->SetTitleOffset(1.5, "Y") ;
      if ( i == 0 ) {
         h2Ds[i]->SetXTitle("M3h ( GeV/c^{2} )");
         h2Ds[i]->SetYTitle("M2h ( GeV/c^{2} )");
      }
      if ( i == 1 ) {
         h2Ds[i]->SetXTitle(" M3h ( GeV/c^{2} )");
         h2Ds[i]->SetYTitle(" M3L ( GeV/c^{2} )");
      }
      h2Ds[i]->Draw( DrawOpt );
      c1->Update();

      TString plotname1 = hfolder + dNames[i]+"/"+  fileName +"_"+ dNames[i] + "."+plotType ;
      if ( isMCMatched ) plotname1 = hfolder + dNames[i] + "/" + fileName + "_"+ dNames[i]+"_MC."+plotType ;
      c1->Print( plotname1 );

      delete c1;
  }

  cout<<" FINISHED  "<<endl;
}


void WAnalysis::LepTopPlotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt, bool isMCMatched ){
 
  CreateFolders(); 

  TH2D* hM2M3L = h2Ds[0] ;
  TH2D* hM2tM3 = h2Ds[1] ;
  TH2D* hM3M3  = h2Ds[2] ;
  TH2D* hPhiM3 = h2Ds[3] ;
  TH2D* hM2M3  = h2Ds[4] ;

  gStyle->SetOptStat("neirom");
  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(1);
  gStyle->SetPalette(1);

  // plot  M2 vs M3 
  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1);
  gStyle->SetNumberContours(5);
  hM2M3L->Draw( DrawOpt );
  c1->Update();

  c1->cd(2);
  TH1D* hM2M3L_py = hM2M3L->ProjectionY("hM2M3L_py");
  hM2M3L_py->Draw();
  c1->Update();

  c1->cd(3);
  TH1D* hM2M3L_px = hM2M3L->ProjectionX("hM2M3L_px");
  hM2M3L_px->Draw();
  c1->Update();

  c1->cd(4);
  hM2M3L->Draw( "LEGO2Z" );
  c1->Update();

  TString plotname1 = hfolder + "M2M3Lep/"+  fileName +  "_M2M3L."+plotType ;
  if ( isMCMatched ) plotname1 = hfolder + "M2M3Lep/" + fileName + "_M2M3L_MC."+plotType ;
  c1->Print( plotname1 );

  // plot  hadronic M2t M3
  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->Divide(2,2);

  c2->cd(1);
  gStyle->SetStatY(0.55);
  gStyle->SetNumberContours(5);
  hM2tM3->Draw( DrawOpt );
  c2->Update();

  c2->cd(2);
  gStyle->SetStatY(0.90);
  TH1D* hM2tM3_py = hM2tM3->ProjectionY("hM2tM3_py");
  hM2tM3_py->Draw();
  c2->Update();

  c2->cd(3);
  TH1D* hM2tM3_px = hM2tM3->ProjectionX("hM2tM3_px");
  hM2tM3_px->Draw();
  c2->Update();

  c2->cd(4);
  hM2tM3->Draw("LEGO2Z" );
  c2->Update();

  TString plotname2 = hfolder + "M2tM3/" + fileName + "_M2tM3."+plotType ;
  if ( isMCMatched ) plotname2 = hfolder + "M2tM3/" + fileName + "_M2tM3_MC."+plotType ;
  c2->Print( plotname2 );

  // plot  M3(lep) vs M3(had) without JES
  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->Divide(2,2);
  c3->cd(1);
  gStyle->SetStatX(0.90);
  gStyle->SetNumberContours(5);
  hM3M3->Draw( DrawOpt );
  c3->Update();

  c3->cd(2);
  TH1D* hM3M3_py = hM3M3->ProjectionY("hM3M3_py");
  hM3M3_py->Draw();
  c3->Update();

  c3->cd(3);
  TH1D* hM3M3_px = hM3M3->ProjectionX("hM3M3_px");
  hM3M3_px->Draw();
  c3->Update();

  c3->cd(4);
  hM3M3->Draw( "LEGO2Z" );
  c3->Update();

  TString plotname3 = hfolder + "M3M3ReFit/"+ fileName + "_M3M3."+plotType ;
  if ( isMCMatched ) plotname3 = hfolder + "M3M3ReFit/" + fileName + "_M3M3_MC."+plotType ;
  c3->Print( plotname3 );

  // plot  Phi M3 
  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->Divide(2,2);

  c4->cd(1);
  gStyle->SetNumberContours(5);
  hPhiM3->Draw( DrawOpt );
  c4->Update();

  c4->cd(2);
  TH1D* hPhiM3_py = hPhiM3->ProjectionY("hPhiM3_py");
  hPhiM3_py->Draw();
  c4->Update();

  c4->cd(3);
  TH1D* hPhiM3_px = hPhiM3->ProjectionX("hPhiM3_px");
  hPhiM3_px->Draw();
  c4->Update();

  c4->cd(4);
  hPhiM3->Draw("LEGO2Z" );
  c4->Update();

  TString plotname4 = hfolder + "PhiM3/" + fileName + "_PhiM3."+plotType ;
  if ( isMCMatched ) plotname4 = hfolder + "PhiM3/" + fileName + "_PhiM3_MC."+plotType ;
  c4->Print( plotname4 );

  // plot  had M2 vs had M3 
  TCanvas* c5 = new TCanvas("c5","", 800, 600);
  c5->SetGrid();
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->Divide(2,2);
  c5->cd(1);
  gStyle->SetNumberContours(5);
  hM2M3->Draw( DrawOpt );
  c5->Update();

  c5->cd(2);
  TH1D* hM2M3_py = hM2M3->ProjectionY("hM2M3_py");
  hM2M3_py->Draw();
  c5->Update();

  c5->cd(3);
  TH1D* hM2M3_px = hM2M3->ProjectionX("hM2M3_px");
  hM2M3_px->Draw();
  c5->Update();

  c5->cd(4);
  hM2M3->Draw( "LEGO2Z" );
  c5->Update();

  TString plotname5 = hfolder + "M2M3ReFit/"+  fileName +  "_M2M3."+plotType ;
  if ( isMCMatched ) plotname5 = hfolder + "M2M3ReFit/" + fileName + "_M2M3_MC."+plotType ;
  c5->Print( plotname5 );

  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
  cout<<" FINISHED !! "<<endl;

}


void WAnalysis::BJetEff( string fileName ) {

  TFile*  file = NULL ;
  TTree* tr = fitInput->GetTree( fileName, "selJet",  file );

  double pt, bDis ;
  int evtId,objId ;
  tr->SetBranchAddress("pt"    ,&pt);
  tr->SetBranchAddress("qCut"  ,&bDis );
  tr->SetBranchAddress("evtId" ,&evtId);
  tr->SetBranchAddress("objId" ,&objId);

  int NEvt = 0 ;
  int nbtagged = 0 ;
  int NbJ[3] = { 0. } ;
  int evtId0 = 0 ;
  for ( int j= 0 ; j< tr->GetEntries() ; j++ ) {
      tr->GetEntry(j);
      if ( evtId0 !=  evtId ) {
         NEvt ++;
         if ( nbtagged == 0 ) NbJ[0] += 1 ;
         if ( nbtagged == 1 ) NbJ[1] += 1 ;
         if ( nbtagged == 2 ) NbJ[2] += 1 ;
         nbtagged = 0;
         evtId0 = evtId ;
      }
      if ( pt < 30 ) continue;
      if ( bDis > bTh ) nbtagged++ ;

  }
  cout<<" **ALLEvt = "<< NEvt <<"  0B = "<< NbJ[0] <<"  1B = "<< NbJ[1] <<" 2B = "<< NbJ[2] <<endl;
  if ( file != NULL ) file->Close() ;

}


