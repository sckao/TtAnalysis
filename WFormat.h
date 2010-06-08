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
  virtual ~recoObj(){}
  virtual void Fillh( double weight, double scale = 1. )  = 0;
  virtual void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ) = 0 ;
  virtual void getlep( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ) = 0 ;
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
  double hadM2_Mt ;
  double lepM2_Mt ;
  double lepM3_eta;
  double hadM3_eta;
  double lepM3_Y;
  double hadM3_Y;

  double MET ;
  double hadF ;
  double lepF ;

  public:

  hadWBoson(){
    hM2M3    = new TH2D("hM2M3", " M3(X) vs M2(Y) ", 18, 0, 360, 15, 0, 300);
    hM3M3    = new TH2D("hM3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
    hM2M2t   = new TH2D("hM2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310);
    hEtaM2   = new TH2D("hEtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hEtaM3   = new TH2D("hEtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hYM2     = new TH2D("hYM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hYM3     = new TH2D("hYM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hM2MET   = new TH2D("hM2MET", " M2 had(X) vs MET(Y) ", 15, 10, 310, 12, 0, 240);
    hM2dF    = new TH2D("hM2dF", " M2 had(X) vs dPhi( M2 ) ", 15, 10, 310, 32, 0., 3.2);
    hMETM2t  = new TH2D("hMETM2t", " MET(X) vs lep M2t(Y) ", 24, 0, 240, 15, 0, 150);
  }

  virtual ~hadWBoson(){
    delete hM2M3 ;
    delete hM3M3 ;
    delete hM2M2t ;
    delete hEtaM2 ;
    delete hEtaM3 ;
    delete hYM2 ;
    delete hYM3 ;
    delete hM2MET ;
    delete hM2dF ;
    delete hMETM2t ;
  }

  void Fillh( double weight, double scale = 1. ) {
       hM2M3->Fill( hadM3, hadM2, weight*scale );
       hM3M3->Fill( hadM3, lepM3, weight*scale );
       hM2M2t->Fill( hadM2_Mt, hadM2, weight*scale );
       hEtaM2->Fill( hadM2_eta, lepM2_eta, weight*scale );
       hEtaM3->Fill( hadM3_eta, lepM3_eta, weight*scale );
       hYM2->Fill( hadM2_Y, lepM2_Y, weight*scale );
       hYM3->Fill( hadM3_Y, lepM3_Y, weight*scale );
       hM2MET->Fill( hadM2, MET, weight*scale );
       double dF = ( fabs(hadF-lepF) < 3.1416 ) ? fabs(hadF-lepF) : 6.2832 - fabs(hadF-lepF) ;
       hM2dF->Fill( hadM2, dF , weight*scale );
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
       hadM2_Mt  = ( Mt2 >= 0 ) ? sqrt(Mt2) : -1 ;
       hadM2_eta = vM2.Eta() ;
       hadM3_eta = vM3.Eta() ;
       hadM2_Y = vM2.Rapidity() ;
       hadM3_Y = vM3.Rapidity() ;
       hadF = vM2.Phi() ;
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
       lepF = vM2.Phi() ;
  }
  void scale( double scale ) {
       hM2M3->Scale( scale );
       hM3M3->Scale( scale );
       hM2M2t->Scale( scale );
       hEtaM2->Scale( scale );
       hEtaM3->Scale( scale );
       hYM2->Scale( scale );
       hYM3->Scale( scale );
       hM2MET->Scale( scale );
       hM2dF->Scale( scale );
       hMETM2t->Scale( scale );
  }
  vector<TH2D*> Output2D(){
       h2Dlist.clear() ;
       h2Dlist.push_back(hM2M3);
       h2Dlist.push_back(hM3M3);
       h2Dlist.push_back(hM2M2t);
       h2Dlist.push_back(hEtaM2);
       h2Dlist.push_back(hEtaM3);
       h2Dlist.push_back(hYM2);
       h2Dlist.push_back(hYM3);
       h2Dlist.push_back(hM2MET);
       h2Dlist.push_back(hM2dF);
       h2Dlist.push_back(hMETM2t);
       return h2Dlist ;
  }
  void Fill2DVec( vector<TH2D*>& hList ){
       hList.push_back(hM2M3);
       hList.push_back(hM3M3);
       hList.push_back(hM2M2t);
       hList.push_back(hEtaM2);
       hList.push_back(hEtaM3);
       hList.push_back(hYM2);
       hList.push_back(hYM3);
       hList.push_back(hM2MET);
       hList.push_back(hM2dF);
       hList.push_back(hMETM2t);
  }

  TH2D* hM2M3 ;
  TH2D* hM3M3 ;
  TH2D* hEtaM2 ;
  TH2D* hEtaM3 ;
  TH2D* hYM2 ;
  TH2D* hYM3 ;
  TH2D* hM2M2t ;
  TH2D* hM2MET ;
  TH2D* hM2dF ;
  TH2D* hMETM2t ;
  vector<TH2D*>  h2Dlist ;

  //ClassDef(hadWBoson, 1);
};

class bgCounter : public recoObj {

  public:

  bgCounter(){
     hadM2M3  = 0 ;
     shadM2M3 = 0 ;
     M3M3     = 0 ;
     sM3M3    = 0 ;

     lepM3 = 0 ;
     hadM3 = 0 ;
     hadM2 = 0 ;
     outV.clear();
  }

  virtual ~bgCounter(){

  }

  void Fillh( double weight, double scale = 1. ) {
       hadM2M3  = hadM2M3  + (weight*scale) ; 
       shadM2M3 = shadM2M3 + (weight*scale*scale) ; 
       M3M3   = M3M3   + (weight*scale) ; 
       sM3M3  = sM3M3  + (weight*scale*scale) ; 
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
  }
  void scale( double scale ) {
       hadM2M3 = hadM2M3 * scale ;
       M3M3    = M3M3 * scale;
  }
  vector<double> Output(){
      outV.clear();
      outV.push_back( hadM2M3 );
      outV.push_back( shadM2M3 );
      outV.push_back( M3M3 );
      outV.push_back( sM3M3 );
      return outV ;
  }

  double hadM2M3 ;
  double shadM2M3 ;
  double M3M3 ;
  double sM3M3 ;
 
  double hadM2 ;
  double hadM3 ;
  double lepM3 ;

  vector<double> outV  ;

  //ClassDef(bgCounter, 1);
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
  double lepW_mt ;

  public:

  hObjs(){
    J0Pt    = new TH1D("J0Pt", " J0 Pt  ", 70, 0, 350);
    J1Pt    = new TH1D("J1Pt", " J1 Pt  ", 70, 0, 350);
    J2Pt    = new TH1D("J2Pt", " J2 Pt  ", 70, 0, 350);
    J3Pt    = new TH1D("J3Pt", " J3 Pt  ", 70, 0, 350);
    muPt    = new TH1D("muPt", " muon Pt ", 70, 0, 350);
    metH    = new TH1D("metH", " MET     ", 70, 0, 350);

    muEta   = new TH1D("muEta", " muon Eta ", 55, -2.7, 2.7 );
    lepWMt  = new TH1D("lepWMt", " lep M2t ", 15, 0, 150);
  }

  virtual ~hObjs(){
    delete J0Pt ;
    delete J1Pt ;
    delete J2Pt ;
    delete J3Pt ;
    delete muPt ;
    delete metH ;
    delete muEta ;
    delete lepWMt ;
  }

  void Fillh( double weight =  1., double scale = 1. ) {
       J0Pt->Fill( pt0, weight*scale );
       J1Pt->Fill( pt1, weight*scale );
       J2Pt->Fill( pt2, weight*scale );
       J3Pt->Fill( pt3, weight*scale );
       muPt->Fill( pt4, weight*scale );
       metH->Fill( pt5, weight*scale );

       muEta->Fill( mu_eta, weight*scale );
       lepWMt->Fill( lepW_mt, weight*scale );
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
  }

  void Fill1DVec( vector<TH1D*>& hList ){
       hList.push_back(J0Pt);
       hList.push_back(J1Pt);
       hList.push_back(J2Pt);
       hList.push_back(J3Pt);
       hList.push_back(muPt);
       hList.push_back(metH);
       hList.push_back(muEta);
       hList.push_back(lepWMt);
  }

  TH1D* J0Pt ;
  TH1D* J1Pt ;
  TH1D* J2Pt ;
  TH1D* J3Pt ;
  TH1D* muPt ;
  TH1D* metH ;
  TH1D* muEta ;
  TH1D* lepWMt ;

  //ClassDef(hObjs, 1);
};

//#if !defined(__CINT__)
//    ClassImp(WFormat);
#endif

