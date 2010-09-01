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
void WAnalysis::HadTopFitter( string fileName, TString DrawOpt, bool isMCMatched  ){
 
  double scale = 1 ;
  string channelType ;
  if ( fileName.size() > 2) {
     vector<string>  channel;
     fitInput->GetParameters( "channel" , &channel );
     for (size_t i=0; i < channel.size(); i++ ) {
         if ( channel[i] == fileName.substr(0,2) )  channelType = channel[i] ;
         if ( fileName.substr(0,2) == "qc" )  channelType = channel[i] ;
         if ( fileName.substr(0,4) == "data" )  channelType = "data" ;
     }
  }
  if ( channelType != "data" )   scale = fitInput->NormalizeComponents( channelType );

  // refit the solutions
  hadWBoson* wbh = new hadWBoson();
  if ( isMCMatched && channelType != "data" ) {
     wmfitter->MCSolution( fileName, wbh );
  } else {
     wmfitter->ReFitSolution( fileName, wbh, n_Jets, scale, NULL, 0, smearing );
  }
  // retrieve the histograms
  //wbh->scale( scale );	
  vector<TH2D*> h2Ds = wbh->Output2D();

  M2M3Plotter( h2Ds, fileName, DrawOpt, isMCMatched );

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

void WAnalysis::Had_SBRatio(){

  TString filepath = hfolder+"SBRatio40.log" ;
  FILE* logfile = fopen( filepath,"a");
  fprintf(logfile," WMt    M2M3   dM3   S/B   SS   N_tt  N_bg  N_wj  N_qcd  Eff_tt  Eff_bg  Eff_wj  Eff_qcd  \n" ) ;  

  vector<string> File4J ;
  fitInput->GetParameters( "FakeData", &File4J );

  //vector<TTree*> TrNJ = fitInput->GetForest( "FakeData", "muJets" );

  vector<double> nEvts ;
  fitInput->GetParameters( "nEvents" , &nEvts );

  double scaleTt  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ  = fitInput->NormalizeComponents( "wj" );

  TH2D* hSB = new TH2D("hSB", " S/B with cuts  dM3(X), M2M3(Y) ", 7, 7.5, 42.5, 14, 2.5, 72.5 );
  TH2D* hCE = new TH2D("hCE", " Cut Eff for signal dM3(X), M2M3(Y) ", 7, 7.5, 42.5, 14, 2.5, 72.5 );

  double nAll[4] = {0.};
  double cEff[4] = {0.};
  cout<<" ---  Measuring SB Ratio( tt / wj+qcd ) --- "<<endl;
  for ( int i = 0; i< 1; i++ ) {
    for ( int j = 0; j< 6; j++ ) {
      for ( int k= 0; k< 5; k++ ) {

         // reset the M2M3 windows
         double lepm2tL = 40. + i*10 ;
	 double window  = 10. + j*10 ;
	 double m2L = 80. - window ;
	 double m2H = 80. + window ;
	 double m3L = 170. - window ;
	 double m3H = 170. + window ;
         double dM3 = 10. + k*5 ;

	 if ( i == 0 && j == 0 && k == 0 ) wmfitter->ResetCuts(  0., 300.,  0., 350.,      0., 999. );

         if ( i != 0 && j == 0 && k == 0 ) wmfitter->ResetCuts(  0., 300.,  0., 350., lepm2tL, 999. );
         if ( i == 0 && j != 0 && k == 0 ) wmfitter->ResetCuts( m2L,  m2H, m3L,  m3H,      0., 999. );
         if ( i == 0 && j == 0 && k != 0 ) wmfitter->ResetCuts(  50., 110.,140., 200.,     0.,  dM3 );

         //if ( i == 0 && j != 0 && k != 0 ) continue;         
         if ( i != 0 && j == 0 && k != 0 ) continue;         
         if ( i != 0 && j != 0 && k == 0 ) continue;         

         if ( i ==0 && j!=0 && k !=0 ) {
            wmfitter->ResetCuts( m2L, m2H, m3L, m3H, lepm2tL, dM3 );
            //cout<<" new LepM2t threshold = "<< lepm2tL<<" , ("<< m2L<<","<<m2H<<","<<m3L<<","<<m3H<<")"<< endl;
         }

         ACounter* sg4j = new ACounter();
	 wmfitter->ReFitSolution( File4J[0], sg4j, 4, scaleTt, NULL, 0, smearing );
	 vector<double> nSg4J = sg4j->Output();

	 ACounter* bg4j = new ACounter();
	 wmfitter->ReFitSolution( File4J[1], bg4j, 4, scaleWJ, NULL, 0, smearing );
	 vector<double> nWj4J = bg4j->Output();
	 wmfitter->ReFitSolution( File4J[2], bg4j, 4, scaleQCD, NULL, 0, smearing );
         vector<double> nBg4J = bg4j->Output();

         // get the denumerator and cut efficiency for each channel
         if ( i == 0 && j == 0 && k == 0 ) {
            nAll[0] = nSg4J[0] ;
	    nAll[1] = nWj4J[0] ;
	    nAll[2] = nBg4J[0] - nWj4J[0];
	    nAll[3] = nBg4J[0] ;
         }
	 cEff[0] = nSg4J[0]/nAll[0] ;
	 cEff[1] = nWj4J[0]/nAll[1] ;
	 cEff[2] = (nBg4J[0] - nWj4J[0]) /nAll[2] ;
	 cEff[3] = nBg4J[0]/nAll[3] ;

         double ssR = nSg4J[0] / sqrt ( nSg4J[0] + nBg4J[0] );
         //double sbR = nSg4J[0] / sqrt ( nBg4J[0] );
         double sbR = nSg4J[0] / nBg4J[0] ;
         //cout<<" M2M3Counter Signal = "<< nSg4J[0] <<" Background = "<< nBg4J[0]  <<endl;
	 //cout<<" The Cut Efficiency = "<< cEff[0] <<" WJ = "<< cEff[1] <<" QCD = "<< cEff[2] <<" BG = "<< cEff[3] <<endl;
	 //cout<<"  *** S/B = "<< nSg4J[0] / nBg4J[0] << endl;
	 //cout<<""<<endl; 

         if ( i == 0 && j != 0 && k !=0 ) {
            hSB->Fill( dM3, window, sbR );
            hCE->Fill( dM3, window, cEff[0] );
         }
         if (  cEff[0] > 0. ) {
            fprintf(logfile," %.1f   %.1f   %1.f   %.1f   %.1f   %.1f   %.1f   %.1f   %.1f   %.2f   %.2f   %.2f   %.2f \n",  
                         lepm2tL, window, dM3, sbR, ssR, nSg4J[0], nBg4J[0], nWj4J[0], (nBg4J[0] - nWj4J[0]), cEff[0], cEff[3], cEff[1], cEff[2]);
         }

          delete sg4j;
          delete bg4j;
      }
    }
  }
  cout<<" SB Ration Measured !! "<<endl;

  fclose(logfile);

  gStyle->SetOptStat("nuom");
  //gStyle->SetNumberContours(10);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.99);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  gStyle->SetPaintTextFormat("2.1f");
  c1->cd();
  hSB->SetMarkerSize(2.5);
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
  hCE->SetMarkerSize(2.5);
  hCE->Draw("TEXT");
  c2->Update();

  TString plotname2 = hfolder + "CutEff."+plotType ;
  c2->Print( plotname2 );

  delete c1;
  delete c2;
  delete hSB ;
  delete hCE ;

}

// for new soltree
void WAnalysis::MixBG( TString DrawOpt ){
 
  // refit the solutions
  cout<<" Mix Background  "<<endl;
  vector<string> flist;
  //fitInput->GetParameters( "FakeData", &flist );
  fitInput->GetParameters( "2JSamples", &flist );

  // define the histograms
  hadWBoson* bg = new hadWBoson();

  double scale1 = fitInput->NormalizeComponents( "wj" );
  cout<<" all bg wjet scale = "<< scale1 <<endl;
  wmfitter->ReFitSolution( flist[1], bg, 2, scale1, NULL, 0, smearing );

  double scale2 = fitInput->NormalizeComponents( "qcd" );
  cout<<" all bg qcd scale = "<< scale2 <<endl;
  wmfitter->ReFitSolution( flist[2], bg, 2, scale2, NULL, 0, smearing );

  vector<TH2D*> h2Ds ;
  bg->Fill2DVec( h2Ds );

  M2M3Plotter( h2Ds, "allBG", DrawOpt ) ;

  delete bg;

}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::MixAll( vector<string>& flist, TString DrawOpt ){
 
  cout<<" Mix All  "<<endl;

  // combined plotting
  gStyle->SetOptStat("neirom");
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.90);
  gStyle->SetStatTextColor(1);
  gStyle->SetPalette(1);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);

  hadWBoson* allh  = new hadWBoson();

  double scale1 = fitInput->NormalizeComponents( "wj" );
  wmfitter->ReFitSolution( flist[1], allh, 4, scale1, NULL, 0, smearing );

  double scale2 = fitInput->NormalizeComponents( "qcd" );
  wmfitter->ReFitSolution( flist[2], allh, 4, scale2, NULL, 0, smearing );

  double scale3 = fitInput->NormalizeComponents( "zj" );
  wmfitter->ReFitSolution( flist[3], allh, 4, scale3, NULL, 0, smearing );

  vector<TH2D*> hbgs ;
  allh->Fill2DVec( hbgs );
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
 
  vector<TH2D*> h2Ds ;
  allh->Fill2DVec( h2Ds );

  //M2M3Plotter( h2Ds, "ALL", DrawOpt ) ;

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

  TString dNames[10] = { "M2M3", "M3M3",  "EtaM2", "EtaM3", "YM2",
                       "YM3",  "M2MET", "M3Pt",  "METMt", "M2M3BG" };

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

      TCanvas* c1 = new TCanvas("c1","", 800, 600);
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


