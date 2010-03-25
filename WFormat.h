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
  virtual void Fillh( double weight )  = 0;
  virtual void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ) = 0 ;
  virtual void getlep( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ) = 0 ;
  virtual void scale( double scale ) = 0;
  virtual vector<TH2D*> Output2D() = 0;
  //virtual vector<TH1D*> Output1D() = 0;

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
  double lepM3_eta;
  double hadM3_eta;
  double lepM3_Y;
  double hadM3_Y;

  public:

  hadWBoson(){
    hM2M3    = new TH2D("hM2M3", " M3(X) vs M2(Y) ", 15, 50, 350, 15, 10, 310);
    hM3M3    = new TH2D("hM3M3", " M3 had(X) vs M3 lep(Y) ", 15,50,350, 15, 50, 350);
    hM2M2t   = new TH2D("hM2M2t", " M2t(X) vs M2(Y) ", 15, 10, 310, 15, 10, 310);
    hEtaM2   = new TH2D("hEtaM2", " Eta of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hEtaM3   = new TH2D("hEtaM3", " Eta of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hYM2     = new TH2D("hYM2", " Rapidity of M2 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
    hYM3     = new TH2D("hYM3", " Rapidity of M3 had(X) vs lep(Y) ", 25, -5, 5, 25, -5, 5);
  }
  virtual ~hadWBoson(){
    delete hM2M3 ;
    delete hM3M3 ;
    delete hM2M2t ;
    delete hEtaM2 ;
    delete hEtaM3 ;
    delete hYM2 ;
    delete hYM3 ;
  }

  void Fillh( double weight) {
       hM2M3->Fill( hadM3, hadM2, weight );
       hM3M3->Fill( hadM3, lepM3, weight );
       hM2M2t->Fill( hadM2_Mt, hadM2, weight );
       hEtaM2->Fill( hadM2_eta, lepM2_eta, weight );
       hEtaM3->Fill( hadM3_eta, lepM3_eta, weight );
       hYM2->Fill( hadM2_Y, lepM2_Y, weight );
       hYM3->Fill( hadM3_Y, lepM3_Y, weight );
  }
  void gethad( TLorentzVector v0, TLorentzVector v1, TLorentzVector v2 ){
       TLorentzVector vM2 = v0 + v1 ;
       TLorentzVector vM3 = v0 + v1 + v2;
       double Mt2 = vM2.Et()*vM2.Et() - vM2.Pt()*vM2.Pt() ;
       hadM2 = vM2.M() ;
       hadM3 = vM3.M() ;
       hadM2_Mt  = ( Mt2 >= 0 ) ? sqrt(Mt2) : -1 ;
       hadM2_eta = vM2.Eta() ;
       hadM3_eta = vM3.Eta() ;
       hadM2_Y = vM2.Rapidity() ;
       hadM3_Y = vM3.Rapidity() ;
  }
  void getlep( TLorentzVector v3, TLorentzVector v4, TLorentzVector v5 ){
       TLorentzVector vM2 = v4 + v5 ;
       TLorentzVector vM3 = v3 + v4 + v5;
       lepM3 = vM3.M() ;
       lepM2_Y = vM2.Rapidity() ;
       lepM3_Y = vM3.Rapidity() ;
       lepM2_eta = vM2.Eta() ;
       lepM3_eta = vM3.Eta() ;
  }
  void scale( double scale ) {
       hM2M3->Scale( scale );
       hM3M3->Scale( scale );
       hM2M2t->Scale( scale );
       hEtaM2->Scale( scale );
       hEtaM3->Scale( scale );
       hYM2->Scale( scale );
       hYM3->Scale( scale );
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
       return h2Dlist ;
  }
  /*
  vector<TH1D*> Output1D(){
       vector<TH1D*>  h1Dlist ;
       return h1Dlist ;
  }
  */
  TH2D* hM2M3 ;
  TH2D* hM3M3 ;
  TH2D* hEtaM2 ;
  TH2D* hEtaM3 ;
  TH2D* hYM2 ;
  TH2D* hYM3 ;
  TH2D* hM2M2t ;
  vector<TH2D*>  h2Dlist ;

  //ClassDef(hadWBoson, 1);
};


//#if !defined(__CINT__)
//    ClassImp(WFormat);
#endif

