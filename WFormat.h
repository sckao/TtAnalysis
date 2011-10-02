#ifndef WFormat_H
#define WFormat_H

#include <TObject.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

class recoObj {

  public:

  recoObj(){}
  recoObj( string fsuffix, double iniV ){}
  virtual ~recoObj(){}
  virtual void Fillh( double weight, double scale = 1. )  = 0;
  virtual void getErr( double nPass, double nFail, double scale = 1. ) = 0 ;
  virtual void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ) = 0 ;
  virtual void getlep( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ) = 0 ;
  virtual void getFloats( double fArr[] ) = 0;
  virtual void getIntegrals( int iArr[] ) = 0;
  virtual void scale( double scale ) = 0;
  //virtual vector<TH2D*> Output2D() = 0;

};

class hadWBoson : public recoObj {

  double hadM2 ;
  double hadM3 ;
  double lepM3 ;
  double lepM2_eta;
  double hadM2_eta;
  double lepM2_Y;
  double hadM2_Y;
  double lepM2_Mt ;
  double lepM3_eta;
  double hadM3_eta;
  double lepM3_Y;
  double hadM3_Y;
  double lepM3_Pt;
  double hadM3_Pt;

  double MET ;

  TLorentzVector allJP4;

  public:

  hadWBoson(){
    hM2M3    = new TH2D("hM2M3",  " M3(X) vs M2(Y) ", 20, 50, 350, 16, 0, 240);
    //hM2M3    = new TH2D("hM2M3",  " M3(X) vs M2(Y) ", 25, 100, 600, 20, 0, 400);
    hM2M3BG  = new TH2D("hM2M3BG"," LepM3(X) vs HadM2(Y) ", 18, 0, 360, 15, 0, 300);
    hM3M3    = new TH2D("hM3M3", " M3 had(X) vs M3 lep(Y) ", 20, 50, 350, 20, 50, 350);
    //hM3M3    = new TH2D("hM3M3", " M3 had(X) vs M3 lep(Y) ", 25, 100, 600, 25, 100, 600);
    hEtaM2   = new TH2D("hEtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hEtaM3   = new TH2D("hEtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hYM2     = new TH2D("hYM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hYM3     = new TH2D("hYM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hM2MET   = new TH2D("hM2MET", " M2 had(X) vs MET(Y) ", 15, 10, 310, 12, 0, 240);
    hM3Pt    = new TH2D("hM3Pt", " Pt of hadT(X) vs Pt of lepT(Y) ", 20, 0, 200, 20, 0, 200);
    hMETM2t  = new TH2D("hMETM2t", " MET(X) vs lep Mt(Y) ", 24, 0, 240, 15, 0, 150);
  }

  virtual ~hadWBoson(){
    delete hM2M3 ;
    delete hM2M3BG ;
    delete hM3M3 ;
    delete hEtaM2 ;
    delete hEtaM3 ;
    delete hYM2 ;
    delete hYM3 ;
    delete hM2MET ;
    delete hM3Pt ;
    delete hMETM2t ;
  }

  void Fillh( double weight, double scale = 1. ) {
       hM2M3->Fill( hadM3, hadM2, weight*scale );
       hM2M3BG->Fill( lepM3, hadM2, weight*scale );
       if (hadM2 > 70 && hadM2 < 90 ) hM3M3->Fill( hadM3, lepM3, weight*scale );
       hEtaM2->Fill( hadM2_eta, lepM2_eta, weight*scale );
       hEtaM3->Fill( hadM3_eta, lepM3_eta, weight*scale );
       hYM2->Fill( hadM2_Y, lepM2_Y, weight*scale );
       hYM3->Fill( hadM3_Y, lepM3_Y, weight*scale );
       hM2MET->Fill( hadM2, MET, weight*scale );
       hM3Pt->Fill( hadM3_Pt, lepM3_Pt , weight*scale );
       hMETM2t->Fill( MET, lepM2_Mt , weight*scale );
  }
  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       TLorentzVector vM2 = v0 + v1 ;
       TLorentzVector vM3 = v0 + v1 + v2;
       //double dphi = v0.DeltaPhi( v1 ) ;
       //double Mt2 = 2.*v0.Pt()*v1.Pt()*( 1. - cos(dphi) );
       hadM2 = vM2.M() ;
       hadM3 = vM3.M() ;
       hadM2_eta = vM2.Eta() ;
       hadM3_eta = vM3.Eta() ;
       hadM2_Y = vM2.Rapidity() ;
       hadM3_Y = vM3.Rapidity() ;
       hadM3_Pt = vM3.Pt() ;
       //cout<<" hadM3: "<< hadM3 <<"  hadM2: "<< hadM2 ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       TLorentzVector vM3 = v3 + v4 + v5;
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       lepM2_Mt = sqrt( Mt2 ) ;
       lepM3 = vM3.M() ;
       lepM2_Y = vM2.Rapidity() ;
       lepM3_Y = vM3.Rapidity() ;
       lepM2_eta = vM2.Eta() ;
       lepM3_eta = vM3.Eta() ;
       MET = v5.Pt() ;
       lepM3_Pt = vM3.Pt() ;
       //cout<<"   lepM3: "<< lepM3 <<"   Neutrino Pz: "<< v5.Pz() <<endl;
  }

  void getFloats( double fArr[] ) { }
  void getIntegrals( int iArr[] ) { }
  void getErr( double np, double nf, double scale ) { }

  void scale( double scale ) {
       hM2M3->Scale( scale );
       hM2M3BG->Scale( scale );
       hM3M3->Scale( scale );
       hEtaM2->Scale( scale );
       hEtaM3->Scale( scale );
       hYM2->Scale( scale );
       hYM3->Scale( scale );
       hM2MET->Scale( scale );
       hM3Pt->Scale( scale );
       hMETM2t->Scale( scale );
  }
  vector<TH2D*> Output2D(){
       h2Dlist.clear() ;
       h2Dlist.push_back(hM2M3);
       h2Dlist.push_back(hM3M3);
       h2Dlist.push_back(hEtaM2);
       h2Dlist.push_back(hEtaM3);
       h2Dlist.push_back(hYM2);
       h2Dlist.push_back(hYM3);
       h2Dlist.push_back(hM2MET);
       h2Dlist.push_back(hM3Pt);
       h2Dlist.push_back(hMETM2t);
       h2Dlist.push_back(hM2M3BG);
       return h2Dlist ;
  }
  void Fill2DVec( vector<TH2D*>& hList ){
       hList.push_back(hM2M3);
       hList.push_back(hM3M3);
       hList.push_back(hEtaM2);
       hList.push_back(hEtaM3);
       hList.push_back(hYM2);
       hList.push_back(hYM3);
       hList.push_back(hM2MET);
       hList.push_back(hM3Pt);
       hList.push_back(hMETM2t);
       hList.push_back(hM2M3BG);
  }

  TH2D* hM2M3 ;
  TH2D* hM3M3 ;
  TH2D* hEtaM2 ;
  TH2D* hEtaM3 ;
  TH2D* hYM2 ;
  TH2D* hYM3 ;
  TH2D* hM2MET ;
  TH2D* hM3Pt ;
  TH2D* hMETM2t ;
  TH2D* hM2M3BG ;
  vector<TH2D*>  h2Dlist ;

  //ClassDef(hadWBoson, 1);
};

class hTopo : public recoObj {

  double met ;
  double Mt ;
  double hadM2 ;
  double hadM3 ;
  double lepM3 ;

  public:

  hTopo( string fsfx = "" ){
    TString sfx = fsfx ;
    hMT      = new TH1D("hMT"+sfx,    " Mt          ", 15,  0, 150);
    M2had    = new TH1D("M2had"+sfx,  " hadronic M2 ", 16,  0, 240);
    M3had    = new TH1D("M3had"+sfx,  " hadronic M3 ", 20, 50, 350);
    M3lep    = new TH1D("M3lep"+sfx,  " leptonic M3 ", 20, 50, 350);

    hMT_1    = new TH1D("hMT_1"+sfx,    " Mt          ", 15,  0, 150);
    M2had_1  = new TH1D("M2had_1"+sfx,  " hadronic M2 ", 12,  0, 240);
    M3had_1  = new TH1D("M3had_1"+sfx,  " hadronic M3 ", 15, 50, 350);
    M3lep_1  = new TH1D("M3lep_1"+sfx,  " leptonic M3 ", 15, 50, 350);

    hMT_2    = new TH1D("hMT_2"+sfx,    " Mt          ", 15,  0, 150);
    M2had_2  = new TH1D("M2had_2"+sfx,  " hadronic M2 ", 12,  0, 240);
    M3had_2  = new TH1D("M3had_2"+sfx,  " hadronic M3 ", 15, 50, 350);
    M3lep_2  = new TH1D("M3lep_2"+sfx,  " leptonic M3 ", 15, 50, 350);

    hMT_3    = new TH1D("hMT_3"+sfx,    " Mt          ", 15,  0, 150);
    M2had_3  = new TH1D("M2had_3"+sfx,  " hadronic M2 ", 12,  0, 240);
    M3had_3  = new TH1D("M3had_3"+sfx,  " hadronic M3 ", 15, 50, 350);
    M3lep_3  = new TH1D("M3lep_3"+sfx,  " leptonic M3 ", 15, 50, 350);

    hMT_4    = new TH1D("hMT_4"+sfx,    " Mt          ", 15,  0, 150);
    M2had_4  = new TH1D("M2had_4"+sfx,  " hadronic M2 ", 12,  0, 240);
    M3had_4  = new TH1D("M3had_4"+sfx,  " hadronic M3 ", 15, 50, 350);
    M3lep_4  = new TH1D("M3lep_4"+sfx,  " leptonic M3 ", 15, 50, 350);
  }

  virtual ~hTopo(){
    delete hMT ;
    delete M2had ;
    delete M3had ;
    delete M3lep ;
    delete hMT_1 ;
    delete M2had_1 ;
    delete M3had_1 ;
    delete M3lep_1 ;
    delete hMT_2 ;
    delete M2had_2 ;
    delete M3had_2 ;
    delete M3lep_2 ;
    delete hMT_3 ;
    delete M2had_3 ;
    delete M3had_3 ;
    delete M3lep_3 ;
    delete hMT_4 ;
    delete M2had_4 ;
    delete M3had_4 ;
    delete M3lep_4 ;
  }

  void Fillh( double weight, double scale = 1. ) {
       M2had->Fill( hadM2, weight*scale );
       M3had->Fill( hadM3, weight*scale );
       M3lep->Fill( lepM3, weight*scale );
       hMT->Fill(      Mt, weight*scale );
       
       if ( hadM2 <  60  )                 hMT_1->Fill(      Mt, weight*scale );
       if ( hadM2 >=  60 && hadM2 <= 100 ) hMT_2->Fill(      Mt, weight*scale );
       if ( hadM2 > 100  )                 hMT_3->Fill(      Mt, weight*scale );
       if ( met   > 20   )                 hMT_4->Fill(      Mt, weight*scale );

       if ( hadM3 < 140  )                 M2had_1->Fill( hadM2, weight*scale );
       if ( hadM3 >= 140 && hadM3 <= 170 ) M2had_2->Fill( hadM2, weight*scale );
       if ( hadM3 > 170  )                 M2had_3->Fill( hadM2, weight*scale );
       if ( met   > 20   )                 M2had_4->Fill( hadM2, weight*scale );

       if ( hadM2 <  60  )                 M3had_1->Fill( hadM3, weight*scale );
       if ( hadM2 >=  60 && hadM2 <= 100 ) M3had_2->Fill( hadM3, weight*scale );
       if ( hadM2 > 100  )                 M3had_3->Fill( hadM3, weight*scale );
       if ( met   > 20   )                 M3had_4->Fill( hadM3, weight*scale );

       if ( hadM3 < 140                  && hadM2 >= 65 && hadM2 <= 95 ) M3lep_1->Fill( lepM3, weight*scale );
       if ( hadM3 >= 140 && hadM3 <= 170 && hadM2 >= 65 && hadM2 <= 95 ) M3lep_2->Fill( lepM3, weight*scale );
       if ( hadM3 > 170                  && hadM2 >= 65 && hadM2 <= 95 ) M3lep_3->Fill( lepM3, weight*scale );
       if ( met   > 20                   && hadM2 >= 65 && hadM2 <= 95 ) M3lep_4->Fill( lepM3, weight*scale );
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       TLorentzVector vM2 = v0 + v1 ;
       TLorentzVector vM3 = v0 + v1 + v2;
       hadM2 = ( vM2.M() < 240 ) ? vM2.M() : 239. ;
       hadM3 = ( vM3.M() < 350 ) ? vM3.M() : 349. ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       TLorentzVector vM3 = v3 + v4 + v5;
       met = v5.Pt() ;
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       Mt    = ( sqrt( Mt2 ) < 150 ) ? sqrt( Mt2 ) : 149. ;
       lepM3 = (     vM3.M() < 350 ) ?     vM3.M() : 349. ;
  }

  void getFloats( double fArr[] ) { }
  void getIntegrals( int iArr[] ) { }
  void getErr( double np, double nf, double scale ) { }

  void scale( double scale ) {
       hMT->Scale( scale );
       M2had->Scale( scale );
       M3had->Scale( scale );
       M3lep->Scale( scale );
       hMT_1->Scale( scale );
       M2had_1->Scale( scale );
       M3had_1->Scale( scale );
       M3lep_1->Scale( scale );
       hMT_2->Scale( scale );
       M2had_2->Scale( scale );
       M3had_2->Scale( scale );
       M3lep_2->Scale( scale );
       hMT_3->Scale( scale );
       M2had_3->Scale( scale );
       M3had_3->Scale( scale );
       M3lep_3->Scale( scale );
       hMT_4->Scale( scale );
       M2had_4->Scale( scale );
       M3had_4->Scale( scale );
       M3lep_4->Scale( scale );
  }

  void Fill1DVec( vector<TH1D*>& hList ){
       hList.push_back(hMT);
       hList.push_back(M2had);
       hList.push_back(M3had);
       hList.push_back(M3lep);
       hList.push_back(hMT_1);
       hList.push_back(M2had_1);
       hList.push_back(M3had_1);
       hList.push_back(M3lep_1);
       hList.push_back(hMT_2);
       hList.push_back(M2had_2);
       hList.push_back(M3had_2);
       hList.push_back(M3lep_2);
       hList.push_back(hMT_3);
       hList.push_back(M2had_3);
       hList.push_back(M3had_3);
       hList.push_back(M3lep_3);
       hList.push_back(hMT_4);
       hList.push_back(M2had_4);
       hList.push_back(M3had_4);
       hList.push_back(M3lep_4);
  }

  TH1D* hMT ;
  TH1D* M2had ;
  TH1D* M3had ;
  TH1D* M3lep ;
  TH1D* hMT_1 ;
  TH1D* M2had_1 ;
  TH1D* M3had_1 ;
  TH1D* M3lep_1 ;
  TH1D* hMT_2 ;
  TH1D* M2had_2 ;
  TH1D* M3had_2 ;
  TH1D* M3lep_2 ;
  TH1D* hMT_3 ;
  TH1D* M2had_3 ;
  TH1D* M3had_3 ;
  TH1D* M3lep_3 ;
  TH1D* hMT_4 ;
  TH1D* M2had_4 ;
  TH1D* M3had_4 ;
  TH1D* M3lep_4 ;

  //ClassDef(hTopo, 1);
};

class ACounter : public recoObj {

  public:

  ACounter( double muCut = 20 ){
     hadM2M3  = 0 ;
     shadM2M3 = 0 ;
     w2np     = 0 ;
     w2nf     = 0 ;
     wnp      = 0 ;
     wnf      = 0 ;
     wn       = 0 ;

     outV.clear();
  }

  virtual ~ACounter(){

  }

  void Fillh( double weight, double scale = 1. ) {
       hadM2M3  = hadM2M3  + (weight*scale) ; 
       shadM2M3 = shadM2M3 + (weight*scale*scale) ; 
  }

  void getErr( double np, double nf, double scale = 1. ) {
       wn  = wn + scale*( np + nf ) ;
       wnp = wnp + scale*np ;
       wnf = wnf + scale*nf ;
       w2np = w2np + scale*scale*np ;
       w2nf = w2nf + scale*scale*nf ;
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){ }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){ }

  void getFloats( double fArr[] ) { }
  void getIntegrals( int iArr[] ) { }

  void scale( double scale ) {
       hadM2M3 = hadM2M3 * scale ;
  }
  void Reset(){
     hadM2M3  = 0 ;
     shadM2M3 = 0 ;

     w2np = 0 ;
     w2nf = 0 ;
     wnp = 0 ;
     wnf = 0 ;
     wn  = 0 ;
     outV.clear();
  }

  vector<double> Output(){
      outV.clear();
      outV.push_back( hadM2M3 );
      outV.push_back( shadM2M3 );
      outV.push_back( wnp );
      outV.push_back( w2np );
      outV.push_back( wnf );
      outV.push_back( w2nf );
      outV.push_back( wn );
      
      return outV ;
  }

  double hadM2M3 ;
  double shadM2M3 ;
 
  double w2np ;
  double w2nf ;
  double wnp ;
  double wnf ;
  double wn ;

  vector<double> outV  ;

  //ClassDef(ACounter, 1);
};

class lepWBoson : public recoObj {

  double hadM2 ;
  double hadM3 ;
  double lepM3 ;
  double lepM2t;
  double lepM3_phi;
  double hadM3_phi;

  public:

  lepWBoson(){
    hM2M3    = new TH2D("hM2M3", " M3 had(X) vs M2 had(Y) ", 15, 50, 350, 15, 10, 310);
    hM2M3L   = new TH2D("hM2M3L", " M3 lep(X) vs M2 had(Y) ", 15, 50, 350, 15, 10, 310);
    hM2tM3   = new TH2D("hM2tM3", " M3(X) vs M2t(Y) ", 15, 50, 350, 20, 0, 100);
    hM3M3    = new TH2D("hM3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
    hPhiM3   = new TH2D("hPhiM3", " Phi of M3 had(X) vs lep(Y) ", 32, -3.2, 3.2, 32, -3.2, 3.2);
  }

  virtual ~lepWBoson(){
    delete hM2M3  ;
    delete hM2M3L  ;
    delete hM2tM3 ;
    delete hM3M3  ;
    delete hPhiM3 ;
  }

  void Fillh( double weight, double scale = 1. ) {
       hM2M3->Fill(  hadM3, hadM2, weight*scale );
       hM2M3L->Fill(  lepM3, hadM2, weight*scale );
       hM2tM3->Fill( lepM3, lepM2t, weight*scale );
       hM3M3->Fill( hadM3, lepM3, weight*scale );
       hPhiM3->Fill( hadM3_phi, lepM3_phi, weight*scale );
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       TLorentzVector vM2 = v0 + v1 ;
       TLorentzVector vM3 = v0 + v1 + v2;
       hadM2 = vM2.M() ;
       hadM3 = vM3.M() ;
       hadM3_phi = vM3.Phi() ;
  }

  //             v3 : bjet          v4 : muon        v5 : neutrino             
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       TLorentzVector vM3 = v3 + v4 + v5;
       //double Mt2 = (v4.Et()+v5.Et())*(v4.Et()+v5.Et()) - ( vM2.Pt()*vM2.Pt() );
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       lepM2t = sqrt( Mt2 ) ;
       lepM3 = vM3.M() ;
       lepM3_phi = vM3.Phi() ;
  }

  void getFloats( double fArr[] ) { }
  void getIntegrals( int iArr[] ) { }
  void getErr( double np, double nf, double scale ) { }

  void scale( double scale ) {
       hM2M3->Scale( scale );
       hM2M3L->Scale( scale );
       hM2tM3->Scale( scale );
       hM3M3->Scale( scale );
       hPhiM3->Scale( scale );
  }

  void Fill2DVec( vector<TH2D*>& hList ){
       hList.push_back(hM2M3L);
       hList.push_back(hM2tM3);
       hList.push_back(hM3M3);
       hList.push_back(hPhiM3);
       hList.push_back(hM2M3);
  }

  TH2D* hM2M3L ;
  TH2D* hM2tM3 ;
  TH2D* hM3M3 ;
  TH2D* hPhiM3 ;
  TH2D* hM2M3 ;
  vector<TH2D*>  h2Dlist ;

  //ClassDef(lepWBoson, 1);
};


class hObjs : public recoObj {

  double pt0 ;
  double pt1 ;
  double pt2 ;
  double pt3 ;
  double pt4 ;
  double pt5 ;
  double j0_eta ;
  double j1_eta ;
  double mu_eta ;
  double mu_iso ;
  double mu_nhits ;
  double mu_d0 ;
  double mu_x2 ;
  double dRmj ;
  double relPt ;
  double Htlep ;
  double Ht ;

  double lepW_pt ;
  double lepW_mt ;
  double lepW_mt1 ;
  double lepW_mt2 ;
  double lepW_mt3 ;
  double lepW_mt4 ;
  double NlepW_all ;
  double NlepW ;
  double WPtCut ;
  double dPhi_MuMet ;

  int njets ;

  public:

  hObjs( string fsfx = "", int bs = 2. ){
    TString sfx = fsfx ;
    J0Pt    = new TH1D("J0Pt"+sfx, " J0 Pt  ", 15*bs, 0, 150);
    J1Pt    = new TH1D("J1Pt"+sfx, " J1 Pt  ", 15*bs, 0, 150);
    J2Pt    = new TH1D("J2Pt"+sfx, " J2 Pt  ", 15*bs, 0, 150);
    J3Pt    = new TH1D("J3Pt"+sfx, " J3 Pt  ", 15*bs, 0, 150);
    J0Eta   = new TH1D("J0Eta"+sfx, " Jet 0 Eta ", 25, -2.6, 2.6 );
    J1Eta   = new TH1D("J1Eta"+sfx, " Jet 1 Eta ", 25, -2.6, 2.6 );

    muPt    = new TH1D("muPt"+sfx, " muon Pt ", 15*bs, 0, 150);
    metH    = new TH1D("metH"+sfx, " MET     ", 15*bs, 0, 150);
    muEta   = new TH1D("muEta"+sfx, " muon Eta ", 25, -2.6, 2.6 );
    lepWPt  = new TH1D("lepWPt"+sfx, "  Pt(Mu,MET)", 10*bs, 0, 200);
    lepWMt  = new TH1D("lepWMt"+sfx, " Mt(Mu,MET) ",  15*bs, 0, 150);
    lepWMt1  = new TH1D("lepWMt1"+sfx, " Mt(Mu,MET)", 15*bs, 0, 150);
    lepWMt2  = new TH1D("lepWMt2"+sfx, " Mt(Mu,MET)", 15*bs, 0, 150);
    lepWMt3  = new TH1D("lepWMt3"+sfx, " Mt(Mu,MET)", 15*bs, 0, 150);
    lepWMt4  = new TH1D("lepWMt4"+sfx, " Mt(Mu,MET)", 15, 0, 150);

    hNJ      = new TH1D("hNJ"+sfx,     " n jets", 7, -0.5, 6.5) ;
    //hNJ      = new TH1D("hNJ"+sfx,     " n jets", 6, -0.5, 5.5) ;
    hMuIso   = new TH1D("hMuIso"+sfx,  " muon RelIso ", 10*bs, 0, 0.15) ; 
    hMuNHits = new TH1D("hMuNHits"+sfx," muon N Hits ", 20*bs, 0, 40 ) ; 
    hMuD0    = new TH1D("hMuD0"+sfx,   " muon d0(Bsp)", 20*bs, -0.02, 0.02 ) ; 
    hMuX2    = new TH1D("hMuX2"+sfx,   " muon X2 "    , 110, 0, 11 ) ; 
    hdRmj    = new TH1D("hdRmj"+sfx,   " dR( mu, jet) "     ,  63, 0, 6.3 ) ; 
    hRelPt   = new TH1D("hRelPt"+sfx,  " RelPt( mu, jet)"  ,  25*bs, 0, 100 ) ; 
    hHtlep   = new TH1D("hHtlep"+sfx,  " Ht_lep"  ,      20*bs, 0, 600 ) ; 
    hHt      = new TH1D("hHt"+sfx,     " Ht_tot"  ,      20*bs, 0, 600 ) ; 

    hMETMt   = new TH2D("hMETMt"+sfx, " MET(X) vs lep Mt(Y) ",          15*bs, 0, 150, 15*bs, 0, 150);
    hWPtMt   = new TH2D("hWPtMt"+sfx, " Pt of W(X) vs lep Mt(Y) ",      10*bs, 0, 200, 15*bs, 0, 150);
    hWPtdPhi = new TH2D("hWPtdPhi"+sfx, " Pt of W(X) vs dPhi(mu, MET)", 10*bs, 0, 200, 16*bs, 0, 3.2);
    hMtdPhi  = new TH2D("hMtdPhi"+sfx, " Mt vs dPhi(mu, MET) ",         15*bs, 0, 150, 16*bs, 0, 3.2);
    hMuMET   = new TH2D("hMuMET"+sfx, " Mu(X) vs MET(Y) ",              15*bs, 0, 120, 15*bs, 0, 120);

    NlepW     = 0 ;
    NlepW_all = 0 ;
    WPtCut    = 20 ;
  }

  virtual ~hObjs(){
    delete J0Pt ;
    delete J1Pt ;
    delete J2Pt ;
    delete J3Pt ;
    delete J0Eta ;
    delete J1Eta ;
    delete muPt ;
    delete metH ;
    delete muEta ;
    delete lepWPt ;
    delete lepWMt ;
    delete lepWMt1 ;
    delete lepWMt2 ;
    delete lepWMt3 ;
    delete lepWMt4 ;
    delete hNJ ;
    delete hMuIso ;
    delete hMuNHits ;
    delete hMuD0 ;
    delete hMuX2 ;
    delete hdRmj ;
    delete hRelPt ;
    delete hHtlep;
    delete hHt;

    delete hMETMt ;
    delete hWPtMt ;
    delete hWPtdPhi ;
    delete hMtdPhi ;
    delete hMuMET ;
  }

  void Fillh( double weight =  1., double scale = 1. ) {
       J0Pt->Fill( pt0, weight*scale );
       J1Pt->Fill( pt1, weight*scale );
       J2Pt->Fill( pt2, weight*scale );
       J3Pt->Fill( pt3, weight*scale );
       J0Eta->Fill( j0_eta, weight*scale );
       J1Eta->Fill( j1_eta, weight*scale );

       muPt->Fill( pt4, weight*scale );
       metH->Fill( pt5, weight*scale );
       muEta->Fill( mu_eta, weight*scale );
       hNJ->Fill( njets, weight*scale );
       hMuIso->Fill( mu_iso, weight*scale );
       hMuNHits->Fill( mu_nhits, weight*scale );
       hMuD0->Fill( mu_d0, weight*scale );
       hMuX2->Fill( mu_x2, weight*scale );
       hdRmj->Fill( dRmj, weight*scale );
       hRelPt->Fill( relPt, weight*scale );

       lepWPt->Fill( lepW_pt, weight*scale );
       lepWMt->Fill( lepW_mt, weight*scale );
       hHtlep->Fill( Htlep, weight*scale );
       hHt->Fill( Ht, weight*scale );

       hMETMt->Fill(       pt5,    lepW_mt, weight*scale );
       hWPtMt->Fill(   lepW_pt,    lepW_mt, weight*scale );
       hWPtdPhi->Fill( lepW_pt, dPhi_MuMet, weight*scale );
       hMtdPhi->Fill(  lepW_mt, dPhi_MuMet, weight*scale );
       hMuMET->Fill(       pt4,        pt5, weight*scale );

       //if (lepW_pt <   WPtCut )                        lepWMt1->Fill( lepW_mt, weight*scale );
       //if (lepW_pt < 2*WPtCut && lepW_pt >=   WPtCut ) lepWMt2->Fill( lepW_mt, weight*scale );
       //if (lepW_pt < 3*WPtCut && lepW_pt >= 2*WPtCut ) lepWMt3->Fill( lepW_mt, weight*scale );
       //if (                      lepW_pt >= 3*WPtCut ) lepWMt4->Fill( lepW_mt, weight*scale );

       if ( njets == 1 ) lepWMt1->Fill( lepW_mt, weight*scale );
       if ( njets == 2 ) lepWMt2->Fill( lepW_mt, weight*scale );
       if ( njets == 3 ) lepWMt3->Fill( lepW_mt, weight*scale );
       if ( njets >= 4 ) lepWMt4->Fill( lepW_mt, weight*scale );

       if (lepW_pt < 3*WPtCut && lepW_mt >= 40 ) NlepW     = NlepW + (weight*scale) ;
       if (lepW_pt < 3*WPtCut )                  NlepW_all = NlepW_all + (weight*scale) ;
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       pt0 = v0.Pt() ;
       pt1 = ( v1.Pt() != 0 ) ? v1.Pt() : -1 ;
       pt2 = ( v2.Pt() != 0 ) ? v2.Pt() : -1 ;
       j0_eta = ( v0.Pt() != 0 ) ? v0.Eta() : 3. ;
       j1_eta = ( v1.Pt() != 0 ) ? v1.Eta() : 3. ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       //double Mt2 = (v4.Et()+v5.Et())*(v4.Et()+v5.Et()) - ( vM2.Pt()*vM2.Pt() );
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       pt3 = v3.Pt() ;
       pt4 = v4.Pt() ;
       pt5 = v5.Pt() ;
       mu_eta = v4.Eta() ;
       lepW_mt = sqrt( Mt2 );
       if ( lepW_mt > 150 ) lepW_mt = 149.9 ;
       lepW_pt = vM2.Pt() ;
       dPhi_MuMet = fabs( dphi );
  }


  void scale( double scale = 1. ) {
       J0Pt->Scale( scale );
       J1Pt->Scale( scale );
       J2Pt->Scale( scale );
       J3Pt->Scale( scale );
       J0Eta->Scale( scale );
       J1Eta->Scale( scale );
       muPt->Scale( scale );
       metH->Scale( scale );
       muEta->Scale( scale );
       lepWMt->Scale( scale );
       lepWMt1->Scale( scale );
       lepWMt2->Scale( scale );
       lepWMt3->Scale( scale );
       lepWMt4->Scale( scale );
       lepWPt->Scale( scale );
       hNJ->Scale(scale);
       hMuIso->Scale(scale);
       hMuNHits->Scale(scale);
       hMuD0->Scale(scale);
       hMuX2->Scale(scale);
       hdRmj->Scale(scale);
       hRelPt->Scale(scale);
       hHt->Scale(scale);
       hHtlep->Scale(scale);
  }

  void getFloats( double fArr[] ) { 
       mu_iso   = fArr[0] ;
       mu_d0    = fArr[1] ;
       mu_x2    = fArr[2] ;
       dRmj     = ( fArr[3] < 6.3 ) ? fArr[3] : 6.29 ;
       relPt    = fArr[4] ;
       Ht       = fArr[5] ;
       Htlep    = fArr[6] ;
  }
  void getIntegrals( int iArr[] ) { 
       njets = iArr[0] ;
       //if ( iArr[0] > 4 ) njets = 4 ;
       if ( iArr[0] > 6 ) njets = 6 ;
       mu_nhits = iArr[1] ;
  }
  void getErr( double np, double nf, double scale ) { }

  void Fill1DVec( vector<TH1D*>& hList ){
       hList.push_back(J0Pt);
       hList.push_back(J1Pt);
       hList.push_back(J2Pt);
       hList.push_back(J3Pt);
       hList.push_back(muPt);
       hList.push_back(metH);
       hList.push_back(muEta);
       hList.push_back(hNJ);
       hList.push_back(hMuIso);
       hList.push_back(hMuNHits);

       hList.push_back(hMuD0);
       hList.push_back(hMuX2);
       hList.push_back(lepWPt);
       hList.push_back(lepWMt);
       hList.push_back(lepWMt1);
       hList.push_back(lepWMt2);
       hList.push_back(lepWMt3);
       hList.push_back(lepWMt4);
       hList.push_back(hdRmj);
       hList.push_back(hRelPt);

       hList.push_back(hHtlep);
       hList.push_back(hHt);
       hList.push_back(J0Eta);
       hList.push_back(J1Eta);
  }

  void Fill2DVec( vector<TH2D*>& hList ){
       hList.push_back( hMETMt );
       hList.push_back( hWPtMt );
       hList.push_back( hWPtdPhi );
       hList.push_back( hMuMET );
       hList.push_back( hMtdPhi );
  }

  void CounterVec( vector<double>& cList ) {
      cList.push_back( NlepW );
      cList.push_back( NlepW_all );
  }

  TH1D* J0Pt ;
  TH1D* J1Pt ;
  TH1D* J2Pt ;
  TH1D* J3Pt ;
  TH1D* J0Eta ;
  TH1D* J1Eta ;
  TH1D* muPt ;
  TH1D* metH ;
  TH1D* muEta ;
  TH1D* lepWMt ;
  TH1D* lepWPt ;
  TH1D* hNJ ;
  TH1D* hMuIso ;
  TH1D* hMuNHits ;
  TH1D* hMuD0 ;
  TH1D* hMuX2 ;
  TH1D* hdRmj ;
  TH1D* hRelPt ;
  TH1D* lepWMt1 ;
  TH1D* lepWMt2 ;
  TH1D* lepWMt3 ;
  TH1D* lepWMt4 ;
  TH1D* hHtlep ;
  TH1D* hHt ;

  TH2D* hMETMt ;  
  TH2D* hWPtMt ;  
  TH2D* hWPtdPhi ;  
  TH2D* hMtdPhi ;  
  TH2D* hMuMET ;

  //ClassDef(hObjs, 1);
};

/*
class hNorm : public recoObj {

  double lepW_pt ;
  double lepW_mt ;
  double lepW_mt1 ;
  double lepW_mt2 ;
  double lepW_mt3 ;
  double lepW_mt4 ;
  double NlepW_all ;
  double NlepW ;
  double WPtCut ;

  public:

  hObjs( string fsfx = "", int nbin = 30 ){
    TString sfx = fsfx ;

    lepWMt   = new TH1D("lepWMt"+sfx, " Mt(Mu,MET) ",  nbin, 0, 150);
    lepWMt1  = new TH1D("lepWMt1"+sfx, " Mt(Mu,MET)", nbin, 0, 150);
    lepWMt2  = new TH1D("lepWMt2"+sfx, " Mt(Mu,MET)", nbin, 0, 150);
    lepWMt3  = new TH1D("lepWMt3"+sfx, " Mt(Mu,MET)", nbin, 0, 150);
    lepWMt4  = new TH1D("lepWMt4"+sfx, " Mt(Mu,MET)", 15, 0, 150);

    NlepW     = 0 ;
    NlepW_all = 0 ;
    WPtCut    = 20 ;
  }

  virtual ~hObjs(){
    delete lepWMt ;
    delete lepWMt1 ;
    delete lepWMt2 ;
    delete lepWMt3 ;
    delete lepWMt4 ;
  }

  void Fillh( double weight =  1., double scale = 1. ) {
       lepWMt->Fill( lepW_mt, weight*scale );

       if (lepW_pt <   WPtCut )                        lepWMt1->Fill( lepW_mt, weight*scale );
       if (lepW_pt < 2*WPtCut && lepW_pt >=   WPtCut ) lepWMt2->Fill( lepW_mt, weight*scale );
       if (lepW_pt < 3*WPtCut && lepW_pt >= 2*WPtCut ) lepWMt3->Fill( lepW_mt, weight*scale );
       if (                      lepW_pt >= 3*WPtCut ) lepWMt4->Fill( lepW_mt, weight*scale );

       if (lepW_pt < 3*WPtCut && lepW_mt >= 40 ) NlepW     = NlepW + (weight*scale) ;
       if (lepW_pt < 3*WPtCut )                  NlepW_all = NlepW_all + (weight*scale) ;
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       pt0 = v0.Et() ;
       pt1 = v1.Et() ;
       pt2 = v2.Et() ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       //double Mt2 = (v4.Et()+v5.Et())*(v4.Et()+v5.Et()) - ( vM2.Pt()*vM2.Pt() );
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       pt3 = v3.Et() ;
       pt4 = v4.Pt() ;
       pt5 = v5.Pt() ;
       mu_eta = v4.Eta() ;
       lepW_mt = sqrt( Mt2 );
       lepW_pt = vM2.Pt() ;
       dPhi_MuMet = fabs( dphi );
  }


  void scale( double scale = 1. ) {
       lepWMt->Scale( scale );
       lepWMt1->Scale( scale );
       lepWMt2->Scale( scale );
       lepWMt3->Scale( scale );
       lepWMt4->Scale( scale );
  }

  void getFloats( double fArr[] ) { 
  }
  void getIntegrals( int iArr[] ) { 
       njets = iArr[0] ;
       if ( iArr[0] > 6 ) njets = 6 ;
       mu_nhits = iArr[1] ;
  }
  void getErr( double np, double nf, double scale ) { }

  void Fill1DVec( vector<TH1D*>& hList ){
       hList.push_back(lepWMt);
       hList.push_back(lepWMt1);
       hList.push_back(lepWMt2);
       hList.push_back(lepWMt3);
       hList.push_back(lepWMt4);
  }

  void Fill2DVec( vector<TH2D*>& hList ){
  }


  TH1D* lepWMt ;
  TH1D* lepWMt1 ;
  TH1D* lepWMt2 ;
  TH1D* lepWMt3 ;
  TH1D* lepWMt4 ;

  //ClassDef(hObjs, 1);
};
*/

class hLepM2 : public recoObj {

  double sum_JPt ;
  double lepW_pt ;
  double lepW_mt ;

  double WPtCut ;

  double nW[5] ;
  double sW[5] ;

  TLorentzVector allJP4;

  public:

  hLepM2( string fsfx = "", int nbin = 30 ){
    TString sfx = fsfx ;
    WPtCut = 20. ;

    hWPt0 = new TH1D("hWPt0"+sfx, " lep M2 pt 0 ", 40, 0, 200);
    hWMt0 = new TH1D("hWMt0"+sfx, " lep M2t   0 ", nbin, 0, 150);
    sumJetPt0 = new TH1D("sumJetPt0"+sfx, " vector sum of jets pt 0", 40, 0, 200);

    hWMt1 = new TH1D("hWMt1"+sfx, " lep M2t   1 ", nbin, 0, 150);
    sumJetPt1 = new TH1D("sumJetPt1"+sfx, " vector sum of jets pt 1", 40, 0, 200);

    hWMt2 = new TH1D("hWMt2"+sfx, " lep M2t   2 ", nbin, 0, 150);
    sumJetPt2 = new TH1D("sumJetPt2"+sfx, " vector sum of jets pt 2", 40, 0, 200);

    hWMt3 = new TH1D("hWMt3"+sfx, " lep M2t   3 ", 15, 0, 150);
    sumJetPt3 = new TH1D("sumJetPt3"+sfx, " vector sum of jets pt 3", 40, 0, 200);

    hWMt4 = new TH1D("hWMt4"+sfx, " lep M2t   4 ", 15, 0, 150);
    sumJetPt4 = new TH1D("sumJetPt4"+sfx, " vector sum of jets pt 4", 40, 0, 200);


    for (int i=0; i<5; i++) {
        nW[i] = 0 ;
        sW[i] = 0 ;
    }
  }

  virtual ~hLepM2(){

    delete hWPt0 ;
    delete hWMt0 ;
    delete sumJetPt0 ;
    
    delete hWMt1 ;
    delete sumJetPt1 ;
    
    delete hWMt2 ;
    delete sumJetPt2 ;
    
    delete hWMt3 ;
    delete sumJetPt3 ;
    
    delete hWMt4 ;
    delete sumJetPt4 ;
    
  }

  void Fillh( double weight =  1., double scale = 1. ) {

       hWPt0->Fill( lepW_pt, weight*scale );
       hWMt0->Fill( lepW_mt, weight*scale );
       sumJetPt0->Fill( sum_JPt, weight*scale );
       nW[0] = nW[0] + (weight*scale) ;
       sW[0] = sW[0] + (weight*scale*scale) ;

       if ( lepW_pt < WPtCut ) {
          hWMt1->Fill( lepW_mt, weight*scale );
          sumJetPt1->Fill( sum_JPt, weight*scale );
          nW[1] = nW[1] + (weight*scale) ;
          sW[1] = sW[1] + (weight*scale*scale) ;
       }

       if ( lepW_pt >= 1*WPtCut && lepW_pt < WPtCut*2 ) {
          hWMt2->Fill( lepW_mt, weight*scale );
          sumJetPt2->Fill( sum_JPt, weight*scale );
          nW[2] = nW[2] + (weight*scale) ;
          sW[2] = sW[2] + (weight*scale*scale) ;
       }

       if ( lepW_pt >= 2*WPtCut && lepW_pt < WPtCut*3 ) {
          hWMt3->Fill( lepW_mt, weight*scale );
          sumJetPt3->Fill( sum_JPt, weight*scale );
          nW[3] = nW[3] + (weight*scale) ;
          sW[3] = sW[3] + (weight*scale*scale) ;
       }

       if ( lepW_pt >= 3*WPtCut  ) {
          hWMt4->Fill( lepW_mt, weight*scale );
          sumJetPt4->Fill( sum_JPt, weight*scale );
          nW[4] = nW[4] + (weight*scale) ;
          sW[4] = sW[4] + (weight*scale*scale) ;
       }
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       allJP4 = v0 + v1 + v2 ;
       sum_JPt = ( allJP4.Pt() < 200. ) ? allJP4.Pt() : 200. ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       lepW_mt = sqrt( Mt2 );
       lepW_pt = ( vM2.Pt() < 200. ) ? vM2.Pt() : 200. ;

       allJP4 = allJP4 + v3 ;
       sum_JPt = ( allJP4.Pt() < 200. ) ? allJP4.Pt() : 200. ;
  }

  void getFloats( double fArr[] ) { } 
  void getIntegrals( int iArr[] ) { }
  void getErr( double np, double nf, double scale ) { }

  void scale( double scale = 1. ) {

       hWPt0->Scale( scale );
       hWMt0->Scale( scale );
       sumJetPt0->Scale( scale );
       hWMt1->Scale( scale );
       sumJetPt1->Scale( scale );
       hWMt2->Scale( scale );
       sumJetPt2->Scale( scale );
       hWMt3->Scale( scale );
       sumJetPt3->Scale( scale );
       hWMt4->Scale( scale );
       sumJetPt4->Scale( scale );
  }

  void Fill1DVec( vector<TH1D*>& hList ){
       hList.push_back(hWPt0);
       hList.push_back(hWMt0);
       hList.push_back(sumJetPt0);
       hList.push_back(hWMt1);
       hList.push_back(sumJetPt1);
       hList.push_back(hWMt2);
       hList.push_back(sumJetPt2);
       hList.push_back(hWMt3);
       hList.push_back(sumJetPt3);
       hList.push_back(hWMt4);
       hList.push_back(sumJetPt4);
  }

  void CounterVec( vector<double>& cList ) {
      cList.push_back( nW[0] );
      cList.push_back( nW[1] );
      cList.push_back( nW[2] );
      cList.push_back( nW[3] );
      cList.push_back( nW[4] );
      cList.push_back( sqrt( sW[0] ) );
      cList.push_back( sqrt( sW[1] ) );
      cList.push_back( sqrt( sW[2] ) );
      cList.push_back( sqrt( sW[3] ) );
      cList.push_back( sqrt( sW[4] ) );
  }

  TH1D* hWMt0 ;
  TH1D* hWPt0 ;
  TH1D* sumJetPt0 ;
  TH1D* hWMt1 ;
  TH1D* sumJetPt1 ;
  TH1D* hWMt2 ;
  TH1D* sumJetPt2 ;
  TH1D* hWMt3 ;
  TH1D* sumJetPt3 ;
  TH1D* hWMt4 ;
  TH1D* sumJetPt4 ;

  //ClassDef(hLepM2, 1);
};

class bgCounter : public recoObj {

  double sum_JPt ;
  double lepW_pt ;
  double lepW_mt ;
  double WPtCut ;

  double nW[5] ;
  double sW[5] ;

  TLorentzVector allJP4;

  public:

  bgCounter( string fsfx = "", double iniPt = 20. ){
    WPtCut = iniPt ;

    for (int i=0; i<5; i++) {
        nW[i] = 0 ;
        sW[i] = 0 ;
    }
  }

  virtual ~bgCounter(){

  }

  void Fillh( double weight =  1., double scale = 1. ) {

       nW[0] = nW[0] + (weight*scale) ;
       sW[0] = sW[0] + (weight*scale*scale) ;
       if ( lepW_pt < WPtCut ) {
          nW[1] = nW[1] + (weight*scale) ;
          sW[1] = sW[1] + (weight*scale*scale) ;
       }

       if ( lepW_pt >= 1*WPtCut && lepW_pt < WPtCut*2 ) {
          nW[2] = nW[2] + (weight*scale) ;
          sW[2] = sW[2] + (weight*scale*scale) ;
       }

       if ( lepW_pt >= 2*WPtCut && lepW_pt < WPtCut*3 ) {
          nW[3] = nW[3] + (weight*scale) ;
          sW[3] = sW[3] + (weight*scale*scale) ;
       }

       if ( lepW_pt >= 3*WPtCut  ) {
          nW[4] = nW[4] + (weight*scale) ;
          sW[4] = sW[4] + (weight*scale*scale) ;
       }

  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       allJP4 = v0 + v1 + v2 ;
       sum_JPt = ( allJP4.Pt() < 200. ) ? allJP4.Pt() : 200. ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       lepW_mt = sqrt( Mt2 );
       lepW_pt = ( vM2.Pt() < 200. ) ? vM2.Pt() : 200. ;

       allJP4 = allJP4 + v3 ;
       sum_JPt = ( allJP4.Pt() < 200. ) ? allJP4.Pt() : 200. ;
  }

  void getFloats( double fArr[] ) { } 
  void getIntegrals( int iArr[] ) { }
  void getErr( double np, double nf, double scale ) { }
  void Fill1DVec( vector<TH1D*>& hList ){ }

  void Reset( double iniPt = 20. ){
    WPtCut = iniPt ;

    for (int i=0; i<5; i++) {
        nW[i] = 0 ;
        sW[i] = 0 ;
    }
  }

  void scale( double scale = 1. ) { 
       nW[0] = nW[0]*scale ;
       nW[1] = nW[1]*scale ;
       nW[2] = nW[2]*scale ;
       nW[3] = nW[3]*scale ;
       nW[4] = nW[4]*scale ;
  }

  void CounterVec( vector<double>& cList ) {
      cList.push_back( nW[0] );
      cList.push_back( nW[1] );
      cList.push_back( nW[2] );
      cList.push_back( nW[3] );
      cList.push_back( nW[4] );
      cList.push_back( sqrt( sW[0] ) );
      cList.push_back( sqrt( sW[1] ) );
      cList.push_back( sqrt( sW[2] ) );
      cList.push_back( sqrt( sW[3] ) );
      cList.push_back( sqrt( sW[4] ) );
  }

  //ClassDef(bgCounter, 1);
};

//#if !defined(__CINT__)
//    ClassImp(WFormat);
#endif

