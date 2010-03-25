#include "WAnalysis.h"
#include "WFormat.h"

WAnalysis::WAnalysis( double massL, double massH ){

  fitInput = new MassAnaInput( "had", massL, massH );
  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "bThreshold", &bTh);
  fitInput->GetParameters( "n_btag", &n_btag);

  wmfitter  = new HadWMassFitter( massL, massH );
  pseudoExp = new PseudoExp( massL, massH );

  mL = massL;
  mH = massH;

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
  gSystem->cd("../../");

}

WAnalysis::~WAnalysis(){

  delete fitInput ;

}


// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::TWFitter( string channelName, int type, TString DrawOpt, bool isMCMatched ){
 
  string fileName;

  string channelType = channelName.substr( 0, 2 )  ;
  vector<string> chlist;
  fitInput->GetParameters( "channel", &chlist );

  double scale = 1. ;
  if ( channelType == "tt" ) {
     vector<string> mlist;
     fitInput->GetParameters( "Signals", &mlist );
     for ( size_t j = 0; j< mlist.size(); j++){
         if ( channelName == mlist[j] ) fileName =  mlist[j];
     }
     scale = fitInput->NormalizeComponents( chlist[0] );
     
  } else { 
     vector<string> bglist;
     fitInput->GetParameters( "Backgrounds", &bglist );
     for ( size_t j = 0; j< bglist.size(); j++) {
         if ( channelName == bglist[j] ) {
            fileName = bglist[j] ;
            scale = fitInput->NormalizeComponents( chlist[j+1] );
         }
     }
  }

  // refit the solutions
  hadWBoson* wbh = new hadWBoson();
  wmfitter->ReFitSolution( fileName, wbh , type, isMCMatched );

  // retrieve the histograms
  wbh->scale( scale );	
  vector<TH2D*> h2Ds = wbh->Output2D();

  M2M3Plotter( h2Ds, fileName, DrawOpt, isMCMatched );

  delete wbh ;
  cout<<" DONE !! "<<endl;
}

// For new SolTree
// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::TWFitter1( int njets, string channelName, int type, TString DrawOpt, bool isMCMatched  ){
 
  string fileName;

  string channelType = channelName.substr( 0, 2 )  ;
  vector<string> chlist;
  fitInput->GetParameters( "channel", &chlist );

  double scale = 1. ;
  if ( channelType == "tt" ) {
     vector<string> mlist;
     fitInput->GetParameters( "Signals", &mlist );
     for ( size_t j = 0; j< mlist.size(); j++){
         if ( channelName == mlist[j] ) fileName =  mlist[j];
     }
     scale = fitInput->NormalizeComponents( chlist[0] );
     
  } else { 
     vector<string> bglist;
     fitInput->GetParameters( "Backgrounds", &bglist );
     for ( size_t j = 0; j< bglist.size(); j++) {
         if ( channelName == bglist[j] ) {
            fileName = bglist[j] ;
            scale = fitInput->NormalizeComponents( chlist[j+1] );
         }
     }
  }

  // refit the solutions
  hadWBoson* wbh = new hadWBoson();
  if ( isMCMatched ) {
     wmfitter->MCSolution( fileName, wbh, type );
  } else {
     if ( type == 0 ) wmfitter->ReFitSolution1( fileName, njets, wbh );
     if ( type  > 0 ) wmfitter->ReFitSolution1( fileName, njets, wbh, type );
  }
  // retrieve the histograms
  wbh->scale( scale );	
  vector<TH2D*> h2Ds = wbh->Output2D();

  M2M3Plotter( h2Ds, fileName, DrawOpt, isMCMatched );

  delete wbh ;
  cout<<" DONE !! "<<endl;
}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::EnsembleTest( int type, int randomSeed, TString DrawOpt ){
 
  vector<string> flist;
  fitInput->GetParameters( "FakeData", &flist );

  // ttbar
  hadWBoson* sg1 = new hadWBoson();
  vector< pair<int,int> > sglist = pseudoExp->GetEnsemble( flist[0], 15.7, randomSeed );
  wmfitter->ReFitSolution( flist[0], sg1 , sglist, type );
  vector<TH2D*> hs0 = sg1->Output2D();

  cout<<" tt test done "<<endl;

  // w+jets
  hadWBoson* bg1 = new hadWBoson();
  vector< pair<int,int> > bg1list = pseudoExp->GetEnsemble( flist[1], 1.15, randomSeed );
  wmfitter->ReFitSolution( flist[1], bg1 , bg1list, type );
  vector<TH2D*> hs1 = bg1->Output2D();

  cout<<" wj test done "<<endl;

  // QCD
  hadWBoson* bg2 = new hadWBoson();
  vector< pair<int,int> > bg2list = pseudoExp->GetEnsemble( flist[2], 7.31, randomSeed );
  wmfitter->ReFitSolution( flist[2], bg2 , bg2list, type );
  // retrieve the histograms
  vector<TH2D*> hs2 = bg2->Output2D();

  cout<<" qcd test done "<<endl;

  // Making plots
  M2M3Plotter( hs0, "pseudoSG", DrawOpt );

  TH2D* h0M2M3  = new TH2D("h0M2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
  h0M2M3->Add( hs1[0] ) ;
  h0M2M3->Add( hs2[0] ) ;
  TH2D* h0M3M3  = new TH2D("h0M3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
  h0M3M3->Add( hs1[1] ) ;
  h0M3M3->Add( hs2[1] ) ;
  TH2D* h0M2M2t = new TH2D("h0M2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310); 
  h0M2M2t->Add( hs1[2] );
  h0M2M2t->Add( hs2[2] );
  TH2D* h0EtaM2 = new TH2D("h0EtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM2->Add( hs1[3] );
  h0EtaM2->Add( hs2[3] );
  TH2D* h0EtaM3 = new TH2D("h0EtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM3->Add( hs1[4] );
  h0EtaM3->Add( hs2[4] );
  TH2D* h0YM2 = new TH2D("h0YM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM2->Add( hs1[5] );
  h0YM2->Add( hs2[5] );
  TH2D* h0YM3 = new TH2D("h0YM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM3->Add( hs1[6] );
  h0YM3->Add( hs2[6] );

  vector<TH2D*> h2Ds;
  h2Ds.push_back( h0M2M3 ); 
  h2Ds.push_back( h0M3M3 ); 
  h2Ds.push_back( h0M2M2t ); 
  h2Ds.push_back( h0EtaM2 ); 
  h2Ds.push_back( h0EtaM3 ); 
  h2Ds.push_back( h0YM2 ); 
  h2Ds.push_back( h0YM3 ); 

  M2M3Plotter( h2Ds, "pseudoBG", DrawOpt );

  h0M2M3->Add( hs0[0] ) ;
  h0M3M3->Add( hs0[1] ) ;
  h0M2M2t->Add( hs0[2] );
  h0EtaM2->Add( hs0[3] );
  h0EtaM3->Add( hs0[4] );
  h0YM2->Add( hs0[5] );
  h0YM3->Add( hs0[6] );

  /*
  vector<TH2D*> hAll;
  hAll.push_back( h0M2M3 ); 
  hAll.push_back( h0M3M3 ); 
  hAll.push_back( h0M2M2t ); 
  hAll.push_back( h0EtaM2 ); 
  hAll.push_back( h0EtaM3 ); 
  hAll.push_back( h0YM2 ); 
  hAll.push_back( h0YM3 ); 
  M2M3Plotter( hAll, "pseudoExp", DrawOpt );
  */
  M2M3Plotter( h2Ds, "pseudoExp", DrawOpt );

  delete sg1 ;
  delete bg1 ;
  delete bg2 ;
  cout<<" YA !!! "<<endl;
}


// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::MixBG( int type, TString DrawOpt ){
 
  // refit the solutions
  hadWBoson* bg1 = new hadWBoson();
  wmfitter->ReFitSolution( "wj_336", bg1 , type );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  bg1->scale( scale1 );	
  cout<<" all bg wjet scale = "<< scale1 <<endl;
  vector<TH2D*> hs1 = bg1->Output2D();

  hadWBoson* bg2 = new hadWBoson();
  wmfitter->ReFitSolution( "QCD+", bg2 , type );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  bg2->scale( scale2 );	
  cout<<" all bg qcd scale = "<< scale2 <<endl;
  // retrieve the histograms
  vector<TH2D*> hs2 = bg2->Output2D();

  TH2D* h0M2M3  = new TH2D("h0M2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
  h0M2M3->Add( hs1[0] ) ;
  h0M2M3->Add( hs2[0] ) ;
  TH2D* h0M3M3  = new TH2D("h0M3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
  h0M3M3->Add( hs1[1] ) ;
  h0M3M3->Add( hs2[1] ) ;
  TH2D* h0M2M2t = new TH2D("h0M2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310); 
  h0M2M2t->Add( hs1[2] );
  h0M2M2t->Add( hs2[2] );
  TH2D* h0EtaM2 = new TH2D("h0EtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM2->Add( hs1[3] );
  h0EtaM2->Add( hs2[3] );
  TH2D* h0EtaM3 = new TH2D("h0EtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM3->Add( hs1[4] );
  h0EtaM3->Add( hs2[4] );
  TH2D* h0YM2 = new TH2D("h0YM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM2->Add( hs1[5] );
  h0YM2->Add( hs2[5] );
  TH2D* h0YM3 = new TH2D("h0YM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM3->Add( hs1[6] );
  h0YM3->Add( hs2[6] );

 
  vector<TH2D*> h2Ds;
  h2Ds.push_back( h0M2M3 ); 
  h2Ds.push_back( h0M3M3 ); 
  h2Ds.push_back( h0M2M2t ); 
  h2Ds.push_back( h0EtaM2 ); 
  h2Ds.push_back( h0EtaM3 ); 
  h2Ds.push_back( h0YM2 ); 
  h2Ds.push_back( h0YM3 ); 
  M2M3Plotter( h2Ds, "allBG", DrawOpt ) ;

  delete bg1;
  delete bg2;

}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::MixAll( int type, TString DrawOpt ){
 
  cout<<" Mix All  "<<endl;
  hadWBoson* sg1 = new hadWBoson();
  wmfitter->ReFitSolution( "tt171_336", sg1 , type );
  double scale0 = fitInput->NormalizeComponents( "tt" );
  sg1->scale( scale0 );	
  vector<TH2D*> hs0 = sg1->Output2D();
 
  hadWBoson* bg1 = new hadWBoson();
  wmfitter->ReFitSolution( "wj_336", bg1 , type );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  bg1->scale( scale1 );	
  vector<TH2D*> hs1 = bg1->Output2D();

  hadWBoson* bg2 = new hadWBoson();
  wmfitter->ReFitSolution( "QCD+", bg2 , type );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  bg2->scale( scale2 );	
  vector<TH2D*> hs2 = bg2->Output2D();

  TH2D* h0M2M3  = new TH2D("h0M2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
  h0M2M3->Add( hs0[0] ) ;
  h0M2M3->Add( hs1[0] ) ;
  h0M2M3->Add( hs2[0] ) ;
  TH2D* h0M3M3  = new TH2D("h0M3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
  h0M3M3->Add( hs0[1] ) ;
  h0M3M3->Add( hs1[1] ) ;
  h0M3M3->Add( hs2[1] ) ;
  TH2D* h0M2M2t = new TH2D("h0M2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310); 
  h0M2M2t->Add( hs0[2] );
  h0M2M2t->Add( hs1[2] );
  h0M2M2t->Add( hs2[2] );
  TH2D* h0EtaM2 = new TH2D("h0EtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM2->Add( hs0[3] );
  h0EtaM2->Add( hs1[3] );
  h0EtaM2->Add( hs2[3] );
  TH2D* h0EtaM3 = new TH2D("h0EtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM3->Add( hs0[4] );
  h0EtaM3->Add( hs1[4] );
  h0EtaM3->Add( hs2[4] );
  TH2D* h0YM2 = new TH2D("h0YM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM2->Add( hs0[5] );
  h0YM2->Add( hs1[5] );
  h0YM2->Add( hs2[5] );
  TH2D* h0YM3 = new TH2D("h0YM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM3->Add( hs0[6] );
  h0YM3->Add( hs1[6] );
  h0YM3->Add( hs2[6] );
 
  vector<TH2D*> h2Ds;
  h2Ds.push_back( h0M2M3 ); 
  h2Ds.push_back( h0M3M3 ); 
  h2Ds.push_back( h0M2M2t ); 
  h2Ds.push_back( h0EtaM2 ); 
  h2Ds.push_back( h0EtaM3 ); 
  h2Ds.push_back( h0YM2 ); 
  h2Ds.push_back( h0YM3 ); 

  M2M3Plotter( h2Ds, "ALL", DrawOpt ) ;
  delete sg1;
  delete bg1;
  delete bg2;

}

void WAnalysis::M2M3Plotter( vector<TH2D*> h2Ds, string fileName, TString DrawOpt, bool isMCMatched ){
 
  TH2D* hM2M3  = h2Ds[0] ;
  TH2D* hM3M3  = h2Ds[1] ;
  TH2D* hM2M2t = h2Ds[2] ;
  TH2D* hEtaM2 = h2Ds[3] ;
  TH2D* hEtaM3 = h2Ds[4] ;
  TH2D* hYM2   = h2Ds[5] ;
  TH2D* hYM3   = h2Ds[6] ;

  gStyle->SetOptStat("neirom");
  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(1);
  gStyle->SetPalette(1);

  // plot  M2 vs M3 without JES
  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1);
  gStyle->SetStatX(0.30);
  gStyle->SetNumberContours(5);
  hM2M3->Draw( DrawOpt );
  c1->Update();

  c1->cd(2);
  gStyle->SetStatX(0.90);
  TH1D* hM2M3_py = hM2M3->ProjectionY("hM2M3_py");
  hM2M3_py->Draw();
  c1->Update();

  c1->cd(3);
  TH1D* hM2M3_px = hM2M3->ProjectionX("hM2M3_px");
  hM2M3_px->Draw();
  c1->Update();

  c1->cd(4);
  hM2M3->Draw( "LEGO2Z" );
  c1->Update();

  TString plotname1 = hfolder + "M2M3/"+  fileName +  "_M2M3.gif" ;
  if ( isMCMatched ) plotname1 = hfolder + "M2M3/" + fileName + "_M2M3_MC.gif" ;
  c1->Print( plotname1 );

  // plot  M3(lep) vs M3(had) without JES
  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->Divide(2,2);
  c2->cd(1);
  gStyle->SetStatX(0.90);
  gStyle->SetNumberContours(5);
  hM3M3->Draw( DrawOpt );
  c2->Update();

  c2->cd(2);
  TH1D* hM3M3_py = hM3M3->ProjectionY("hM3M3_py");
  hM3M3_py->Draw();
  c2->Update();

  c2->cd(3);
  TH1D* hM3M3_px = hM2M3->ProjectionX("hM3M3_px");
  hM3M3_px->Draw();
  c2->Update();

  c2->cd(4);
  hM3M3->Draw( "LEGO2Z" );
  c2->Update();

  TString plotname2 = hfolder + "M3M3/"+ fileName + "_M3M3.gif" ;
  if ( isMCMatched ) plotname2 = hfolder + "M3M3/" + fileName + "_M3M3_MC.gif" ;
  c2->Print( plotname2 );

  // plot  Eta M2 
  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->Divide(2,2);
  c3->cd(1);
  gStyle->SetNumberContours(5);
  hEtaM2->Draw( DrawOpt );
  c3->Update();

  c3->cd(2);
  TH1D* hEtaM2_py = hEtaM2->ProjectionY("hEtaM2_py");
  hEtaM2_py->Draw();
  c3->Update();

  c3->cd(3);
  TH1D* hEtaM2_px = hEtaM2->ProjectionX("hEtaM2_px");
  hEtaM2_px->Draw();
  c3->Update();

  c3->cd(4);
  hEtaM2->Draw("LEGO2Z" );
  c3->Update();

  TString plotname3 = hfolder + "EtaM2/" + fileName + "_EtaM2.gif" ;
  if ( isMCMatched ) plotname3 = hfolder + "EtaM2/" + fileName + "_EtaM2_MC.gif" ;
  c3->Print( plotname3 );

  // plot  Eta M3 
  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->Divide(2,2);

  c4->cd(1);
  gStyle->SetNumberContours(5);
  hEtaM3->Draw( DrawOpt );
  c4->Update();

  c4->cd(2);
  TH1D* hEtaM3_py = hEtaM3->ProjectionY("hEtaM3_py");
  hEtaM3_py->Draw();
  c4->Update();

  c4->cd(3);
  TH1D* hEtaM3_px = hEtaM3->ProjectionX("hEtaM3_px");
  hEtaM3_px->Draw();
  c4->Update();

  c4->cd(4);
  hEtaM3->Draw("LEGO2Z" );
  c4->Update();

  TString plotname4 = hfolder + "EtaM3/" + fileName + "_EtaM3.gif" ;
  if ( isMCMatched ) plotname4 = hfolder + "EtaM3/" + fileName + "_EtaM3_MC.gif" ;
  c4->Print( plotname4 );

  // plot  hadronic M2 M2t
  TCanvas* c5 = new TCanvas("c5","", 800, 600);
  c5->SetGrid();
  c5->SetFillColor(10);
  c5->Divide(2,2);

  c5->cd(1);
  gStyle->SetStatY(0.55);
  gStyle->SetNumberContours(5);
  hM2M2t->Draw( DrawOpt );
  c5->Update();

  c5->cd(2);
  gStyle->SetStatY(0.90);
  TH1D* hM2M2t_py = hM2M2t->ProjectionY("hM2M2t_py");
  hM2M2t_py->Draw();
  c5->Update();

  c5->cd(3);
  TH1D* hM2M2t_px = hM2M2t->ProjectionX("hM2M2t_px");
  hM2M2t_px->Draw();
  c5->Update();

  c5->cd(4);
  hM2M2t->Draw("LEGO2Z" );
  c5->Update();

  TString plotname5 = hfolder + "M2M2t/" + fileName + "_M2M2t.gif" ;
  if ( isMCMatched ) plotname5 = hfolder + "M2M2t/" + fileName + "_M2M2t_MC.gif" ;
  c5->Print( plotname5 );

  // plot  Y M2 
  TCanvas* c6 = new TCanvas("c6","", 800, 600);
  c6->SetGrid();
  c6->SetFillColor(10);
  c6->Divide(2,2);
  c6->cd(1);
  gStyle->SetNumberContours(5);
  hYM2->Draw( DrawOpt );
  c6->Update();

  c6->cd(2);
  TH1D* hYM2_py = hYM2->ProjectionY("hYM2_py");
  hYM2_py->Draw();
  c6->Update();

  c6->cd(3);
  TH1D* hYM2_px = hYM2->ProjectionX("hYM2_px");
  hYM2_px->Draw();
  c6->Update();

  c6->cd(4);
  hYM2->Draw("LEGO2Z" );
  c6->Update();

  TString plotname6 = hfolder + "YM2/" + fileName + "_YM2.gif" ;
  if ( isMCMatched ) plotname6 = hfolder + "YM2/" + fileName + "_YM2_MC.gif" ;
  c6->Print( plotname6 );

  // plot  Eta M3 
  TCanvas* c7 = new TCanvas("c7","", 800, 600);
  c7->SetGrid();
  c7->SetFillColor(10);
  c7->Divide(2,2);

  c7->cd(1);
  gStyle->SetNumberContours(5);
  hYM3->Draw( DrawOpt );
  c7->Update();

  c7->cd(2);
  TH1D* hYM3_py = hYM3->ProjectionY("hYM3_py");
  hYM3_py->Draw();
  c7->Update();

  c7->cd(3);
  TH1D* hYM3_px = hYM3->ProjectionX("hYM3_px");
  hYM3_px->Draw();
  c7->Update();

  c7->cd(4);
  hYM3->Draw("LEGO2Z" );
  c7->Update();

  TString plotname7 = hfolder + "YM3/" + fileName + "_YM3.gif" ;
  if ( isMCMatched ) plotname7 = hfolder + "YM3/" + fileName + "_YM3_MC.gif" ;
  c7->Print( plotname7 );

  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
  delete c6;
  delete c7;

  cout<<" FINISHED  "<<endl;
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


// for new soltree
void WAnalysis::MixBG1( int njets, int type, TString DrawOpt ){
 
  // refit the solutions
  hadWBoson* bg1 = new hadWBoson();
  if ( type==0 ) wmfitter->ReFitSolution1( "wj_336a", njets, bg1, type );
  if ( type >0 ) wmfitter->ReFitSolution1( "wj_336a", njets, bg1 );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  bg1->scale( scale1 );	
  cout<<" all bg wjet scale = "<< scale1 <<endl;
  vector<TH2D*> hs1 = bg1->Output2D();

  hadWBoson* bg2 = new hadWBoson();
  wmfitter->ReFitSolution1( "qcd1_336a", njets, bg2 );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  bg2->scale( scale2 );	
  cout<<" all bg qcd scale = "<< scale2 <<endl;

  // retrieve the histograms
  vector<TH2D*> hs2 = bg2->Output2D();

  TH2D* h0M2M3  = new TH2D("h0M2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
  h0M2M3->Add( hs1[0] ) ;
  h0M2M3->Add( hs2[0] ) ;
  TH2D* h0M3M3  = new TH2D("h0M3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
  h0M3M3->Add( hs1[1] ) ;
  h0M3M3->Add( hs2[1] ) ;
  TH2D* h0M2M2t = new TH2D("h0M2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310); 
  h0M2M2t->Add( hs1[2] );
  h0M2M2t->Add( hs2[2] );
  TH2D* h0EtaM2 = new TH2D("h0EtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM2->Add( hs1[3] );
  h0EtaM2->Add( hs2[3] );
  TH2D* h0EtaM3 = new TH2D("h0EtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM3->Add( hs1[4] );
  h0EtaM3->Add( hs2[4] );
  TH2D* h0YM2 = new TH2D("h0YM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM2->Add( hs1[5] );
  h0YM2->Add( hs2[5] );
  TH2D* h0YM3 = new TH2D("h0YM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM3->Add( hs1[6] );
  h0YM3->Add( hs2[6] );

 
  vector<TH2D*> h2Ds;
  h2Ds.push_back( h0M2M3 ); 
  h2Ds.push_back( h0M3M3 ); 
  h2Ds.push_back( h0M2M2t ); 
  h2Ds.push_back( h0EtaM2 ); 
  h2Ds.push_back( h0EtaM3 ); 
  h2Ds.push_back( h0YM2 ); 
  h2Ds.push_back( h0YM3 ); 
  M2M3Plotter( h2Ds, "allBG", DrawOpt ) ;

  delete bg1;
  delete bg2;

}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::MixAll1( int njets, TString DrawOpt ){
 
  cout<<" Mix All  "<<endl;
  hadWBoson* sg1 = new hadWBoson();
  wmfitter->ReFitSolution1( "tt171_336a", njets, sg1 );
  double scale0 = fitInput->NormalizeComponents( "tt" );
  sg1->scale( scale0 );	
  vector<TH2D*> hs0 = sg1->Output2D();
 
  hadWBoson* bg1 = new hadWBoson();
  wmfitter->ReFitSolution1( "wj_336a", njets, bg1 );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  bg1->scale( scale1 );	
  vector<TH2D*> hs1 = bg1->Output2D();

  hadWBoson* bg2 = new hadWBoson();
  wmfitter->ReFitSolution1( "qcd1_336a", njets, bg2 );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  bg2->scale( scale2 );	
  vector<TH2D*> hs2 = bg2->Output2D();

  TH2D* h0M2M3  = new TH2D("h0M2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
  h0M2M3->Add( hs0[0] ) ;
  h0M2M3->Add( hs1[0] ) ;
  h0M2M3->Add( hs2[0] ) ;
  TH2D* h0M3M3  = new TH2D("h0M3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
  h0M3M3->Add( hs0[1] ) ;
  h0M3M3->Add( hs1[1] ) ;
  h0M3M3->Add( hs2[1] ) ;
  TH2D* h0M2M2t = new TH2D("h0M2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310); 
  h0M2M2t->Add( hs0[2] );
  h0M2M2t->Add( hs1[2] );
  h0M2M2t->Add( hs2[2] );
  TH2D* h0EtaM2 = new TH2D("h0EtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM2->Add( hs0[3] );
  h0EtaM2->Add( hs1[3] );
  h0EtaM2->Add( hs2[3] );
  TH2D* h0EtaM3 = new TH2D("h0EtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM3->Add( hs0[4] );
  h0EtaM3->Add( hs1[4] );
  h0EtaM3->Add( hs2[4] );
  TH2D* h0YM2 = new TH2D("h0YM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM2->Add( hs0[5] );
  h0YM2->Add( hs1[5] );
  h0YM2->Add( hs2[5] );
  TH2D* h0YM3 = new TH2D("h0YM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM3->Add( hs0[6] );
  h0YM3->Add( hs1[6] );
  h0YM3->Add( hs2[6] );
 
  vector<TH2D*> h2Ds;
  h2Ds.push_back( h0M2M3 ); 
  h2Ds.push_back( h0M3M3 ); 
  h2Ds.push_back( h0M2M2t ); 
  h2Ds.push_back( h0EtaM2 ); 
  h2Ds.push_back( h0EtaM3 ); 
  h2Ds.push_back( h0YM2 ); 
  h2Ds.push_back( h0YM3 ); 

  M2M3Plotter( h2Ds, "ALL", DrawOpt ) ;
  delete sg1;
  delete bg1;
  delete bg2;

}


// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void WAnalysis::EnsembleTest1( int njets, int randomSeed, TString DrawOpt ){
 
  vector<string> flist;
  fitInput->GetParameters( "FakeData", &flist );

  TString treeName ;
  if ( njets == 3 ) treeName = "mu3Jets" ;
  if ( njets == 4 ) treeName = "mu4Jets" ;

  // ttbar
  hadWBoson* sg1 = new hadWBoson();
  vector<int> sglist = pseudoExp->GetEnsemble( flist[0], treeName, 15.7, randomSeed );
  wmfitter->ReFitSolution1( flist[0], njets, sg1 , &sglist );
  vector<TH2D*> hs0 = sg1->Output2D();

  cout<<" tt test done "<<endl;

  // w+jets
  hadWBoson* bg1 = new hadWBoson();
  vector<int> bg1list = pseudoExp->GetEnsemble( flist[1], treeName, 1.15, randomSeed );
  wmfitter->ReFitSolution1( flist[1], njets, bg1 , &bg1list );
  vector<TH2D*> hs1 = bg1->Output2D();

  cout<<" wj test done "<<endl;

  // QCD
  hadWBoson* bg2 = new hadWBoson();
  vector<int> bg2list = pseudoExp->GetEnsemble( flist[2], treeName, 7.31, randomSeed );
  wmfitter->ReFitSolution1( flist[2], njets, bg2 , &bg2list );
  // retrieve the histograms
  vector<TH2D*> hs2 = bg2->Output2D();

  cout<<" qcd test done "<<endl;

  // Making plots
  M2M3Plotter( hs0, "pseudoSG", DrawOpt );

  TH2D* h0M2M3  = new TH2D("h0M2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
  h0M2M3->Add( hs1[0] ) ;
  h0M2M3->Add( hs2[0] ) ;
  TH2D* h0M3M3  = new TH2D("h0M3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
  h0M3M3->Add( hs1[1] ) ;
  h0M3M3->Add( hs2[1] ) ;
  TH2D* h0M2M2t = new TH2D("h0M2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310); 
  h0M2M2t->Add( hs1[2] );
  h0M2M2t->Add( hs2[2] );
  TH2D* h0EtaM2 = new TH2D("h0EtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM2->Add( hs1[3] );
  h0EtaM2->Add( hs2[3] );
  TH2D* h0EtaM3 = new TH2D("h0EtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0EtaM3->Add( hs1[4] );
  h0EtaM3->Add( hs2[4] );
  TH2D* h0YM2 = new TH2D("h0YM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM2->Add( hs1[5] );
  h0YM2->Add( hs2[5] );
  TH2D* h0YM3 = new TH2D("h0YM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  h0YM3->Add( hs1[6] );
  h0YM3->Add( hs2[6] );

  vector<TH2D*> h2Ds;
  h2Ds.push_back( h0M2M3 ); 
  h2Ds.push_back( h0M3M3 ); 
  h2Ds.push_back( h0M2M2t ); 
  h2Ds.push_back( h0EtaM2 ); 
  h2Ds.push_back( h0EtaM3 ); 
  h2Ds.push_back( h0YM2 ); 
  h2Ds.push_back( h0YM3 ); 

  M2M3Plotter( h2Ds, "pseudoBG", DrawOpt );

  h0M2M3->Add( hs0[0] ) ;
  h0M3M3->Add( hs0[1] ) ;
  h0M2M2t->Add( hs0[2] );
  h0EtaM2->Add( hs0[3] );
  h0EtaM3->Add( hs0[4] );
  h0YM2->Add( hs0[5] );
  h0YM3->Add( hs0[6] );

  /*
  vector<TH2D*> hAll;
  hAll.push_back( h0M2M3 ); 
  hAll.push_back( h0M3M3 ); 
  hAll.push_back( h0M2M2t ); 
  hAll.push_back( h0EtaM2 ); 
  hAll.push_back( h0EtaM3 ); 
  hAll.push_back( h0YM2 ); 
  hAll.push_back( h0YM3 ); 
  M2M3Plotter( hAll, "pseudoExp", DrawOpt );
  */
  M2M3Plotter( h2Ds, "pseudoExp", DrawOpt );

  delete sg1 ;
  delete bg1 ;
  delete bg2 ;
  cout<<" YA !!! "<<endl;
}

