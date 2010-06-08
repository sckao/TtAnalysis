#include "JetSpectrum.h"

JetSpectrum::JetSpectrum() {

  fitInput = new MassAnaInput();
  pseudoExp = new PseudoExp();

  fitInput->GetParameters( "Path", &hfolder );

}

JetSpectrum::~JetSpectrum(){

  delete fitInput ;
  delete pseudoExp;
}


void JetSpectrum::EtSpectrum( string fileName, recoObj* histos, bool smearing ) {

  TFile* file ;
  TTree* tr = fitInput->GetTree( fileName, "muJets",  file );

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

  // this is only for MC Matching
  double jpx[5],jpy[5],jpz[5],jE[5], bDis[5];
  double mpx[2],mpy[2],mpz[2],mE[2] ;
  double npx[2],npy[2],npz[2],nE[2] ;
  int nNu, nJ ;
  tr->SetBranchAddress("jpx"    ,&jpx);
  tr->SetBranchAddress("jpy"    ,&jpy);
  tr->SetBranchAddress("jpz"    ,&jpz);
  tr->SetBranchAddress("jE"     ,&jE);
  tr->SetBranchAddress("nJ"     ,&nJ);
  tr->SetBranchAddress("mpx"    ,&mpx);
  tr->SetBranchAddress("mpy"    ,&mpy);
  tr->SetBranchAddress("mpz"    ,&mpz);
  tr->SetBranchAddress("mE"     ,&mE);
  tr->SetBranchAddress("npx"    ,&npx);
  tr->SetBranchAddress("npy"    ,&npy);
  tr->SetBranchAddress("npz"    ,&npz);
  tr->SetBranchAddress("nE"     ,&nE);
  tr->SetBranchAddress("bTh"    ,&bDis );
  tr->SetBranchAddress("nNu"    ,&nNu );

  double bwmt = 0 ;
  double awmt = 0 ;
  for (int k=0; k< tr->GetEntries() ; k++) {

      tr->GetEntry(k);
      if ( nJ < 2 ) continue ;
      
      vector<TLorentzVector> objlist;
      for ( int i = 0; i< nJ ; i ++) {
          TLorentzVector jn( jpx[i], jpy[i], jpz[i], jE[i] );
          objlist.push_back( jn );
      }
      TLorentzVector m0( mpx[0], mpy[0], mpz[0], mE[0] );
      TLorentzVector n0( npx[0], npy[0], npz[0], nE[0] );
      objlist.push_back( m0 );
      objlist.push_back( n0 );

      if ( smearing )  {
         pseudoExp->PhaseSmearing( objlist, 0 );
         pseudoExp->JetEtReSort( objlist );
      }

      int sz = objlist.size() ;
      TLorentzVector j0( 0., 0., 0., 0. );
      if ( nJ >= 4 ) {      
         histos->gethad( objlist[0], objlist[1], objlist[2] );   
         histos->getlep( objlist[3], objlist[sz-2], objlist[sz-1] );  
         histos->Fillh( 1., scale ) ;

         TLorentzVector vlepM2 = objlist[sz-2] + objlist[sz-1] ;
         double dphi = objlist[sz-2].DeltaPhi( objlist[sz-1] ) ;
         double lepMt2 = 2.*objlist[sz-2].Pt()*objlist[sz-1].Pt()*( 1. - cos(dphi) );
         if ( sqrt(lepMt2) < 40 ) {
            bwmt = bwmt + 1. ;
         } else {
            awmt = awmt + 1. ;
         }
      } 
      else if ( nJ == 3 ) {
         histos->gethad( objlist[0], objlist[1], objlist[2] );   
         histos->getlep( j0, objlist[sz-2], objlist[sz-1] );  
         histos->Fillh( 1., scale ) ;
      } else {
         histos->gethad( objlist[0], objlist[1], j0 );   
         histos->getlep( j0, objlist[sz-2], objlist[sz-1] );  
         histos->Fillh( 1., scale ) ;
      }
  }
  cout<<" ratio wmt = " <<  bwmt / (awmt + bwmt) <<endl;
}


void JetSpectrum::ObjHistoPlotter( string fileName, bool smearing ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hObjects" );
  gSystem->cd( "../" );

  hObjs* hs = new hObjs() ;
  EtSpectrum( fileName, hs, smearing );

  /*
  double scale = 1. ;
  if ( fileName.substr(0,2) == "tt" )  scale = fitInput->NormalizeComponents( "tt" );
  if ( fileName.substr(0,2) == "wj" )  scale = fitInput->NormalizeComponents( "wj" );
  if ( fileName.substr(0,2) == "qc" )  scale = fitInput->NormalizeComponents( "qcd" );
  cout<<" the scale = "<< scale << endl;
  */

  vector<TH1D*> h1Ds ;
  hs->Fill1DVec( h1Ds );

  gStyle->SetOptStat("ieruom");
  gStyle->SetLabelSize( 0.06, "X");
  gStyle->SetLabelSize( 0.06, "Y");
  gStyle->SetLabelSize( 0.06, "Z");

  // Jet Et Spectrum
  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetLogy();
  c1->cd();

  //gStyle->SetTextSize(0.6);
  //gStyle->SetStatFontSize(0.6);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(7);
  h1Ds[3]->SetLineColor(7);
  h1Ds[3]->SetLineWidth(2);
  h1Ds[3]->Draw();
  c1->Update();

  gStyle->SetStatX(0.75);
  gStyle->SetStatTextColor(6);
  h1Ds[2]->SetLineColor(6);
  h1Ds[2]->SetLineWidth(2);
  h1Ds[2]->DrawCopy("sames");
  c1->Update();

  gStyle->SetStatY(0.70);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(4);
  h1Ds[1]->SetLineColor(4);
  h1Ds[1]->SetLineWidth(2);
  h1Ds[1]->DrawCopy("sames");
  c1->Update();

  gStyle->SetStatX(0.75);
  gStyle->SetStatTextColor(2);
  h1Ds[0]->SetLineColor(2);
  h1Ds[0]->SetLineWidth(2);
  h1Ds[0]->DrawCopy("sames");
  c1->Update();

  TString plotname1 = hfolder + "hObjects/"+ fileName + "_JetEt.gif" ;
  c1->Print( plotname1 );

  // Muon and MET Spectrum
  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->cd();
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(2);
  h1Ds[4]->SetLineColor(2);
  h1Ds[4]->SetLineWidth(2);
  h1Ds[4]->Draw();
  c2->Update();

  gStyle->SetStatY(0.7);
  gStyle->SetStatTextColor(4);
  h1Ds[5]->SetLineColor(4);
  h1Ds[5]->SetLineWidth(2);
  h1Ds[5]->DrawCopy("sames");
  c2->Update();

  TString plotname2 = hfolder + "hObjects/"+ fileName + "_LepEt.gif" ;
  c2->Print( plotname2 );


  // MET Eta distribution
  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->cd();

  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(2);
  h1Ds[6]->SetLineColor(1);
  h1Ds[6]->SetLineWidth(2);
  h1Ds[6]->Draw();
  c3->Update();

  TString plotname3 = hfolder + "hObjects/"+ fileName + "_MuonEta.gif" ;
  c3->Print( plotname3 );

  // Leptonic W transverse mass 
  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  c4->cd();

  h1Ds[7]->SetLineWidth(2);
  h1Ds[7]->Draw();
  c4->Update();

  TString plotname4 = hfolder + "hObjects/"+ fileName + "_LepWMt.gif" ;
  c4->Print( plotname4 );

  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete hs;

}

