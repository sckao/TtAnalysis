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
    hM2M3    = new TH2D("hM2M3",  " M3(X) vs M2(Y) ", 18, 0, 360, 15, 0, 300);
    hM2M3BG  = new TH2D("hM2M3BG"," LepM3(X) vs HadM2(Y) ", 18, 0, 360, 15, 0, 300);
    hM3M3    = new TH2D("hM3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
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
       hM3M3->Fill( hadM3, lepM3, weight*scale );
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
       //double Mt2 = (v0.Et()+v1.Et())*(v0.Et()+v1.Et()) - vM2.Pt()*vM2.Pt() ;
       double dphi = v0.DeltaPhi( v1 ) ;
       double Mt2 = 2.*v0.Pt()*v1.Pt()*( 1. - cos(dphi) );
       hadM2 = vM2.M() ;
       hadM3 = vM3.M() ;
       hadM2_eta = vM2.Eta() ;
       hadM3_eta = vM3.Eta() ;
       hadM2_Y = vM2.Rapidity() ;
       hadM3_Y = vM3.Rapidity() ;
       hadM3_Pt = vM3.Pt() ;
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
  }

  void getFloats( double fArr[] ) { }
  void getIntegrals( int iArr[] ) { }

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

class ACounter : public recoObj {

  public:

  ACounter(){
     hadM2M3  = 0 ;
     shadM2M3 = 0 ;
     M3M3     = 0 ;
     sM3M3    = 0 ;
     NlepW_Mt = 0 ;

     lepM3 = 0 ;
     hadM3 = 0 ;
     hadM2 = 0 ;
     lepW_mt = 0 ;
     lepW_pt = 0 ;
     outV.clear();
  }

  virtual ~ACounter(){

  }

  void Fillh( double weight, double scale = 1. ) {
       hadM2M3  = hadM2M3  + (weight*scale) ; 
       shadM2M3 = shadM2M3 + (weight*scale*scale) ; 
       M3M3   = M3M3   + (weight*scale) ; 
       sM3M3  = sM3M3  + (weight*scale*scale) ; 
       if ( lepW_mt > 40. ) NlepW_Mt = NlepW_Mt + (weight*scale) ;
  }

  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       TLorentzVector vM2 = v0 + v1 ;
       TLorentzVector vM3 = v0 + v1 + v2;
       hadM3 = vM3.M() ;
       hadM2 = vM2.M() ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM3 = v3 + v4 + v5;
       lepM3 = vM3.M() ;

       TLorentzVector vM2 = v4 + v5 ;
       double dphi = v4.DeltaPhi( v5 ) ;
       double Mt2 = 2.*v4.Pt()*v5.Pt()*( 1. - cos(dphi) );
       lepW_mt = sqrt( Mt2 );
       lepW_pt = vM2.Pt() ;
  }

  void getFloats( double fArr[] ) { }
  void getIntegrals( int iArr[] ) { }

  void scale( double scale ) {
       hadM2M3 = hadM2M3 * scale ;
       M3M3    = M3M3 * scale;
       NlepW_Mt = NlepW_Mt * scale ;
  }
  vector<double> Output(){
      outV.clear();
      outV.push_back( hadM2M3 );
      outV.push_back( shadM2M3 );
      outV.push_back( M3M3 );
      outV.push_back( sM3M3 );
      outV.push_back( NlepW_Mt );
      return outV ;
  }

  double hadM2M3 ;
  double shadM2M3 ;
  double M3M3 ;
  double sM3M3 ;
  double NlepW_Mt ;
 
  double hadM2 ;
  double hadM3 ;
  double lepM3 ;

  double lepW_mt ;
  double lepW_pt ;
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
  double mu_eta ;
  double mu_iso ;
  double mu_nhits ;
  double mu_d0 ;
  double mu_x2 ;

  double lepW_pt ;
  double lepW_mt ;
  double lepW_mt1 ;
  double lepW_mt2 ;
  double lepW_mt3 ;
  double lepW_mt4 ;
  double NlepW_all ;
  double NlepW ;
  double WPtCut ;

  int njets ;

  public:

  hObjs( string fsfx = "", int nbin = 30 ){
    TString sfx = fsfx ;
    J0Pt    = new TH1D("J0Pt"+sfx, " J0 Pt  ", 30, 0, 150);
    J1Pt    = new TH1D("J1Pt"+sfx, " J1 Pt  ", 30, 0, 150);
    J2Pt    = new TH1D("J2Pt"+sfx, " J2 Pt  ", 30, 0, 150);
    J3Pt    = new TH1D("J3Pt"+sfx, " J3 Pt  ", 30, 0, 150);
    muPt    = new TH1D("muPt"+sfx, " muon Pt ", 30, 0, 150);
    metH    = new TH1D("metH"+sfx, " MET     ", 30, 0, 150);

    muEta   = new TH1D("muEta"+sfx, " muon Eta ", 25, -2.6, 2.6 );
    lepWPt  = new TH1D("lepWPt"+sfx, "  Pt(Mu,MET)", 20, 0, 200);
    lepWMt  = new TH1D("lepWMt"+sfx, " Mt(Mu,MET) ",  nbin, 0, 150);
    lepWMt1  = new TH1D("lepWMt1"+sfx, " Mt(Mu,MET)", nbin, 0, 150);
    lepWMt2  = new TH1D("lepWMt2"+sfx, " Mt(Mu,MET)", nbin, 0, 150);
    lepWMt3  = new TH1D("lepWMt3"+sfx, " Mt(Mu,MET)", 15, 0, 150);
    lepWMt4  = new TH1D("lepWMt4"+sfx, " Mt(Mu,MET)", 15, 0, 150);

    hNJ      = new TH1D("hNJ"+sfx,     " n jets", 9, -0.5, 8.5) ;
    hMuIso   = new TH1D("hMuIso"+sfx,  " muon RelIso ", 20, 0, 0.1) ; 
    hMuNHits = new TH1D("hMuNHits"+sfx," muon N Hits ", 40, 0, 40 ) ; 
    hMuD0    = new TH1D("hMuD0"+sfx,   " muon d0(Bsp)", 40, -0.02, 0.02 ) ; 
    hMuX2    = new TH1D("hMuX2"+sfx,   " muon X2 "    , 110, 0, 11 ) ; 

    hMETMt   = new TH2D("hMETMt"+sfx, " MET(X) vs lep Mt(Y) ", 30, 0, 150, nbin, 0, 150);

    NlepW     = 0 ;
    NlepW_all = 0 ;
    WPtCut    = 20 ;
  }

  virtual ~hObjs(){
    delete J0Pt ;
    delete J1Pt ;
    delete J2Pt ;
    delete J3Pt ;
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
    delete hMETMt ;
    
  }

  void Fillh( double weight =  1., double scale = 1. ) {
       J0Pt->Fill( pt0, weight*scale );
       J1Pt->Fill( pt1, weight*scale );
       J2Pt->Fill( pt2, weight*scale );
       J3Pt->Fill( pt3, weight*scale );
       muPt->Fill( pt4, weight*scale );
       metH->Fill( pt5, weight*scale );

       muEta->Fill( mu_eta, weight*scale );
       hNJ->Fill( njets, weight*scale );
       hMuIso->Fill( mu_iso, weight*scale );
       hMuNHits->Fill( mu_nhits, weight*scale );
       hMuD0->Fill( mu_d0, weight*scale );
       hMuX2->Fill( mu_x2, weight*scale );
       lepWPt->Fill( lepW_pt, weight*scale );
       lepWMt->Fill( lepW_mt, weight*scale );
       hMETMt->Fill( pt5, lepW_mt, weight*scale );

       if (lepW_pt <   WPtCut )                        lepWMt1->Fill( lepW_mt, weight*scale );
       if (lepW_pt < 2*WPtCut && lepW_pt >=   WPtCut ) lepWMt2->Fill( lepW_mt, weight*scale );
       if (lepW_pt < 3*WPtCut && lepW_pt >= 2*WPtCut ) lepWMt3->Fill( lepW_mt, weight*scale );
       if (                      lepW_pt >= 3*WPtCut ) lepWMt4->Fill( lepW_mt, weight*scale );

       if (lepW_pt < 3*WPtCut && lepW_mt >=40 ) NlepW = NlepW + (weight*scale) ;
       if (lepW_pt < 3*WPtCut )                 NlepW_all = NlepW_all + (weight*scale) ;
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
  }


  void scale( double scale = 1. ) {
       J0Pt->Scale( scale );
       J1Pt->Scale( scale );
       J2Pt->Scale( scale );
       J3Pt->Scale( scale );
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
  }

  void getFloats( double fArr[] ) { 
       mu_iso   = fArr[0] ;
       mu_d0    = fArr[1] ;
       mu_x2    = fArr[2] ;
  }
  void getIntegrals( int iArr[] ) { 
       njets = iArr[0] ;
       if ( iArr[0] > 8 ) njets = 8 ;
       mu_nhits = iArr[1] ;
  }

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
  }

  void Fill2DVec( vector<TH2D*>& hList ){
       hList.push_back( hMETMt );
  }

  void CounterVec( vector<double>& cList ) {
      cList.push_back( NlepW );
      cList.push_back( NlepW_all );
  }

  TH1D* J0Pt ;
  TH1D* J1Pt ;
  TH1D* J2Pt ;
  TH1D* J3Pt ;
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
  TH1D* lepWMt1 ;
  TH1D* lepWMt2 ;
  TH1D* lepWMt3 ;
  TH1D* lepWMt4 ;
  TH2D* hMETMt ;  

  //ClassDef(hObjs, 1);
};

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

  void getFloats( double fArr[] ) {} 
  void getIntegrals( int iArr[] ) { }

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
    TString sfx = fsfx ;
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
  void Fill1DVec( vector<TH1D*>& hList ){ }

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

