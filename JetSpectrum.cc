#include "JetSpectrum.h"

JetSpectrum::JetSpectrum( double massL, double massH ) {

  fitInput = new MassAnaInput( "had", massL, massH );

  fitInput->GetParameters( "Path", &hfolder );
  TString theFolder = hfolder ;

  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hObjects" );
  gSystem->cd( "../" );


}

JetSpectrum::~JetSpectrum(){

  delete fitInput ;

}


void JetSpectrum::EtSpectrum( string fileName, recoObj* histos ) {

  TFile* file ;
  TTree* tr = fitInput->GetTree( fileName, "muJets",  file );

  // this is only for MC Matching
  double jpx[5],jpy[5],jpz[5],jE[5], bDis[5];
  double mpx[2],mpy[2],mpz[2],mE[2] ;
  double npx[2],npy[2],npz[2],nE[2] ;
  int nNu ;
  tr->SetBranchAddress("jpx"    ,&jpx);
  tr->SetBranchAddress("jpy"    ,&jpy);
  tr->SetBranchAddress("jpz"    ,&jpz);
  tr->SetBranchAddress("jE"     ,&jE);
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

  for (int k=0; k< tr->GetEntries() ; k++) {
      tr->GetEntry(k);
      //if ( nNu < 2 ) continue ;
      TLorentzVector j0( jpx[0], jpy[0], jpz[0], jE[0] );
      TLorentzVector j1( jpx[1], jpy[1], jpz[1], jE[1] );
      TLorentzVector j2( jpx[2], jpy[2], jpz[2], jE[2] );
      TLorentzVector j3( jpx[3], jpy[3], jpz[3], jE[3] );
      TLorentzVector m0( mpx[0], mpy[0], mpz[0], mE[0] );
      TLorentzVector n0( npx[0], npy[0], npz[0], nE[0] );
      histos->gethad(j0,j1,j2);   
      histos->getlep(j3,m0,n0);  
      histos->Fillh( 1. ) ;
  }
  
}


void JetSpectrum::ObjHistoPlotter( string fileName ){

  hObjs* hs = new hObjs() ;
  EtSpectrum( fileName, hs );

  vector<TH1D*> h1Ds ;
  hs->Fill1DVec( h1Ds );


  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetLogy();
  c1->cd();

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

}

