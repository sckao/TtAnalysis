#ifndef TtAnalysisHisto_H
#define TtAnalysisHisto_H

/** \class TopAnalysisHistograms
 *  Collection of histograms for TopAnalysis test.
 *
 * Author: S.C. Kao  - UC Riverside
 */

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <iostream>

using namespace std;

class HTOP1 {
public:
 
 HTOP1(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hEt    = new TH1F("hEt", " Et distribution",500,0,500);
    hEovH  = new TH1F("hEovH"," em Energy / hadronic Energy ", 500, -2., 23.);
    hNCont = new TH1F("hNCont"," # of constituents of a jet", 50, 0, 50);
    hNTk   = new TH1F("hNTk"," # of charge track of a jet", 50, 0, 50);
    hR60   = new TH1F("hR60", "n60/ total # of constituents", 80, -0.5, 1.5);
    hArea_pt= new TH2F("hArea_pt"," tower area vs pt ", 500, 0., 500., 200, 0., 1.);
    hEovH_pt= new TH2F("hEovH_pt"," em Energy / hadronic Energy vs pt", 500, 0., 500., 500, -2., 23.);
    hEovH_r = new TH2F("hEovH_r"," E/H vs. nChargeTrk/nConstituent ", 12, -0.05, 1.15, 500, -2., 23.);
    hEovH_A = new TH2F("hEovH_A"," E/H vs tower area ", 100, 0., 1., 500, -2., 23.);
    hemF    = new TH2F("hemF", " emF(x) vs emFCalo(y)", 150, -0.2, 1.3, 150, -0.2, 1.3);

    Et_Pt        = new TH2F("Et_Pt", " PAT Et vs Pt (all jets)", 500, 0, 500, 500, 0., 500);
    
    hWmass = new TH1F("hWmass"," W mass from gen-reco matching jets",240,30.,150.);
    hWp    = new TH1F("hWp", " momentum of W ",500,0.,500.);

    EovH1    = new TH1F("EovH1"," em Energy / hadronic Energy ", 500, -2., 23.);
    EovH2    = new TH1F("EovH2"," em Energy / hadronic Energy ", 500, -2., 23.);
    Cont_Trk1= new TH2F("Cont_Trk1","N Constituent vs N Charged Track ", 50, -0.5, 49.5, 50, -0.5, 49.5);
    Cont_Trk2= new TH2F("Cont_Trk2","N Constituent vs N Charged Track ", 50, -0.5, 49.5, 50, -0.5, 49.5);

    hNJets = new TH1F("hNJets", " # of Jets",65,-0.5,64.5);
    hNjet20= new TH1F("hNjet20"," # of Jets, pt > 20 GeV ",65,-0.5,64.5);
    hNjet30= new TH1F("hNjet30"," # of Jets, pt > 30 GeV ",65,-0.5,64.5);
    hNjet40= new TH1F("hNjet40"," # of Jets, pt > 40 GeV ",65,-0.5,64.5);

    dRjj = new TH1F("dRjj"," dR( Jx, Jy ) ", 200,-0.025, 9.975);
    
    hEovH_A1 = new TH2F("hEovH_A1", " E/H vs towersArea", 100, 0., 1., 500, -2., 23.);
    hEovH_r1 = new TH2F("hEovH_r1", " E/H vs nConstituent/nChargeTrk ", 12, -0.05, 1.15, 500, -2., 23.);
    hEovH_N1 = new TH2F("hEovH_N1", " E/H vs nConstituent ", 50, 0., 50., 500, -2., 23.);
    hEovH_C1 = new TH2F("hEovH_C1", " E/H vs nCharge Trk ", 50, 0., 50., 500, -2., 23.);
    hR60_1  = new TH1F("hR60_1", "n60/ total # of constituents", 80, -0.5, 1.5);

    gEta_Pt  = new TH2F("gEta_Pt", " eta vs pt (all gen jets)", 71, -3.55, 3.55, 500, 0., 500);
    gEovH    = new TH1F("gEovH"," E/H for gen Jets", 500, -2., 23.);
    gNJets   = new TH1F("gNJets"," # of genJets",65,-0.5,64.5);
    gNjet20  = new TH1F("gNjet20"," # of genJets, pt > 20 GeV ",65,-0.5,64.5);
    gNjet40  = new TH1F("gNjet40"," # of genJets, pt > 40 GeV ",65,-0.5,64.5);
    gNjeth35 = new TH1F("gNjeth35"," # of genJets, |h| > 3.5 GeV ",65,-0.5,64.5);

    Jet1Et = new TH1F("Jet1Et", " ET of the 1st highest Et jet ", 500,0.,500. ); 
    Jet2Et = new TH1F("Jet2Et", " ET of the 2nd highest Et jet ", 500,0.,500. ); 
    Jet3Et = new TH1F("Jet3Et", " ET of the 3rd highest Et jet ", 500,0.,500. ); 
    Jet4Et = new TH1F("Jet4Et", " ET of the 4th highest Et jet ", 500,0.,500. ); 
    thirdCalEt = new TH1F("thirdCalEt", " ET of the 3rd highest Et jet ", 500,0.,500. ); 

    eta_njets  = new TH2F("eta_njets", " eta vs njets ", 91, -4.55, 4.55, 25, -0.5, 24.5);
    wEta_njets = new TH2F("wEta_njets", "w eta vs njets ", 91, -4.55, 4.55, 25, -0.5, 24.5);

    wEta_dPhi0 = new TH2F("wEta_dPhi0", " w eta vs dPhi(w,jet1)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    wEta_dPhi1 = new TH2F("wEta_dPhi1", " w eta vs dPhi(w,jet2)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    wEta_dPhi2 = new TH2F("wEta_dPhi2", " w eta vs dPhi(w,jet3)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    njets_dPhi0= new TH2F("njets_dPhi0", " njets vs dPhi(w,jet1)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    njets_dPhi1= new TH2F("njets_dPhi1", " njets vs dPhi(w,jet2)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    njets_dPhi2= new TH2F("njets_dPhi2", " njets vs dPhi(w,jet3)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    
    njets_dRMuJ  = new TH2F("njets_dRMuJ", " njets vs min_dR(mu,jet)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    dRMuJ_RelPt  = new TH2F("dRMuJ_RelPt", " min_dR(mu,jet), RelPt  ", 200, 0., 5., 100, 0., 200 );
    muNjets_dPhi0= new TH2F("muNjets_dPhi0"," njets vs dPhi(mu,jet1)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    muNjets_dPhi1= new TH2F("muNjets_dPhi1"," njets vs dPhi(mu,jet2)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    muNjets_dPhi2= new TH2F("muNjets_dPhi2"," njets vs dPhi(mu,jet3)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);

    njets_dR0= new TH2F("njets_dR0", " njets vs dR(w,jet1)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    njets_dR1= new TH2F("njets_dR1", " njets vs dR(w,jet2)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    njets_dR2= new TH2F("njets_dR2", " njets vs dR(w,jet3)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);

    njets_dR = new TH2F("njets_dR", " njets vs dR(w,jet1+2+3)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    njets_m3 = new TH2F("njets_m3", " njets vs m3 of Jet1,2,3", 20, -0.5, 19.5, 500,0.,500. );
    m3_dR    = new TH2F("m3_dR", " m3 vs dR(w,jet1+2+3)", 500, 0., 500., 100, -0.05, 9.95);
    m3_j3    = new TH2F("m3_j3", " m3 vs J3 ET ", 500, 0., 500., 500, 0., 500. );

    allJ_selJ = new TH2F("allJ_selJ", "N of All jets vs N of selected jets in a event", 21, -0.5, 20.5, 21, -0.5, 20.5);
 } 

 HTOP1(TString name_, TFile* file) {
    name=name_;

    hEt    = (TH1F *) file->Get("hEt");
    hEovH  = (TH1F *) file->Get("hEovH");
    hNCont = (TH1F *) file->Get("hNCont");
    hNTk   = (TH1F *) file->Get("hNTK");
    hR60   = (TH1F *) file->Get("hR60");
    hArea_pt = (TH2F *) file->Get("hArea_pt");
    hEovH_pt = (TH2F *) file->Get("hEovH_pt");
    hEovH_r  = (TH2F *) file->Get("hEovH_r");
    hEovH_A  = (TH2F *) file->Get("hEovH_A");
    Et_Pt  = (TH2F *) file->Get("Et_Pt");
    hemF     = (TH2F *) file->Get("hemF");

    hWmass = (TH1F *) file->Get("hWmass");
    hWp    = (TH1F *) file->Get("hWp");

    EovH1   = (TH1F *) file->Get("EovH1");
    EovH2   = (TH1F *) file->Get("EovH2");
    Cont_Trk1 = (TH2F *) file->Get("Cont_Trk1");
    Cont_Trk2 = (TH2F *) file->Get("Cont_Trk2");

    hNJets = (TH1F *) file->Get("hNJets");
    hNjet20= (TH1F *) file->Get("hNjet20");
    hNjet30= (TH1F *) file->Get("hNjet30");
    hNjet40= (TH1F *) file->Get("hNjet40");
    dRjj   = (TH1F *) file->Get("dRjj");

    hEovH_A1  = (TH2F *) file->Get("hEovH_A1");
    hEovH_N1  = (TH2F *) file->Get("hEovH_N1");
    hEovH_C1  = (TH2F *) file->Get("hEovH_C1");
    hEovH_r1  = (TH2F *) file->Get("hEovH_r1");
    hR60_1    = (TH1F *) file->Get("hR60_1");

    gEta_Pt = (TH2F *) file->Get("gEta_Pt");
    gEovH   = (TH1F *) file->Get("gEovH");
    gNJets  = (TH1F *) file->Get("gNJets");
    gNjet20 = (TH1F *) file->Get("gNjet20");
    gNjet40 = (TH1F *) file->Get("gNjet40");
    gNjeth35= (TH1F *) file->Get("gNjeth35");

    Jet1Et = (TH1F *) file->Get("Jet1Et");
    Jet2Et = (TH1F *) file->Get("Jet2Et");
    Jet3Et = (TH1F *) file->Get("Jet3Et");
    Jet4Et = (TH1F *) file->Get("Jet4Et");
    thirdCalEt = (TH1F *) file->Get("thirdCalEt");
    eta_njets  = (TH2F *) file->Get("eta_njets");
    wEta_njets  = (TH2F *) file->Get("wEta_njets");

    wEta_dPhi0 = (TH2F *) file->Get("wEta_dPhi0");
    wEta_dPhi1 = (TH2F *) file->Get("wEta_dPhi1");
    wEta_dPhi2 = (TH2F *) file->Get("wEta_dPhi2");
    njets_dPhi0= (TH2F *) file->Get("njets_dPhi0");
    njets_dPhi1= (TH2F *) file->Get("njets_dPhi1");
    njets_dPhi2= (TH2F *) file->Get("njets_dPhi2");
 
    njets_dRMuJ  = (TH2F *) file->Get("njets_dRMuJ");
    dRMuJ_RelPt  = (TH2F *) file->Get("dRMuJ_RelPt");
    muNjets_dPhi0= (TH2F *) file->Get("muNjets_dPhi0");
    muNjets_dPhi1= (TH2F *) file->Get("muNjets_dPhi1");
    muNjets_dPhi2= (TH2F *) file->Get("muNjets_dPhi2");

    njets_dR0= (TH2F *) file->Get("njets_dR0");
    njets_dR1= (TH2F *) file->Get("njets_dR1");
    njets_dR2= (TH2F *) file->Get("njets_dR2");

    njets_dR = (TH2F *) file->Get("njets_dR");
    m3_dR    = (TH2F *) file->Get("m3_dR");
    m3_j3    = (TH2F *) file->Get("m3_j3");
    njets_m3 = (TH2F *) file->Get("njets_m3");


    allJ_selJ = (TH2F *) file->Get("allJ_selJ");
 }

 /// Destructor
 virtual ~HTOP1() {
    //delete hcsc_q;
    delete hEt;
    delete hEovH;
    delete hNCont;
    delete hNTk;
    delete hR60;
    delete hArea_pt;
    delete hEovH_pt;
    delete hEovH_r;
    delete hEovH_A;
    delete Et_Pt;
    delete hemF;

    delete hWmass;
    delete hWp;

    delete EovH1;
    delete Cont_Trk1;
    delete EovH2;
    delete Cont_Trk2;

    delete hNJets;
    delete hNjet20;
    delete hNjet30;
    delete hNjet40;
    delete dRjj;

    delete hEovH_A1;
    delete hEovH_N1;
    delete hEovH_C1;
    delete hEovH_r1;
    delete hR60_1;

    delete gEta_Pt;
    delete gEovH;
    delete gNJets;
    delete gNjet20;
    delete gNjet40;
    delete gNjeth35;

    delete Jet1Et; 
    delete Jet2Et; 
    delete Jet3Et; 
    delete Jet4Et; 
    delete thirdCalEt;
    delete eta_njets;
    delete wEta_njets;

    delete wEta_dPhi0;
    delete wEta_dPhi1;
    delete wEta_dPhi2;
    delete njets_dPhi0;
    delete njets_dPhi1;
    delete njets_dPhi2;

    delete njets_dRMuJ;
    delete dRMuJ_RelPt;
    delete muNjets_dPhi0;
    delete muNjets_dPhi1;
    delete muNjets_dPhi2;

    delete njets_dR0;
    delete njets_dR1;
    delete njets_dR2;

    delete njets_dR;
    delete m3_dR;
    delete m3_j3;
    delete njets_m3;

    delete allJ_selJ;
 }

 void Fill1a( double pt, double et, double EovH, int nCont, double r60, double r , double area, int nTk, double emF, double emFCalo ){
    hEt->Fill(et);
    hEovH->Fill(EovH);
    hNCont->Fill(nCont);
    hNTk->Fill(nTk);
    hR60->Fill(r60);
    hArea_pt->Fill(pt,area);
    hEovH_pt->Fill(pt,EovH);
    hEovH_r->Fill(r,EovH);
    hEovH_A->Fill(area,EovH);
    hemF->Fill(emF, emFCalo);

    Et_Pt->Fill(et, pt);
 } 
 void Fill1b( double wmass, double wp ){
    hWmass->Fill(wmass);
    hWp->Fill(wp);
 }
 void Fill1c( int nCon, int nTrk, double EovH) {
    EovH1->Fill( EovH );
    Cont_Trk1->Fill(nCon, nTrk);
 }
 void Fill1d( int nCon, int nTrk, double EovH) {
    Cont_Trk2->Fill(nCon, nTrk);
 }
 void Fill1f( int njets, int nj20, int nj40, int njh35 ){
    hNJets->Fill(njets);  
    hNjet20->Fill(nj20);  
    hNjet30->Fill(nj40);  
    hNjet40->Fill(njh35);  
 }
 void Fill1g( double dR ) {
    dRjj->Fill( dR );
 }
 void Fill1h( double EovH, double area, int nCon, int nTrk, double r, double r60 ) {
    hEovH_A1->Fill(area,EovH);
    hEovH_N1->Fill(nCon,EovH);
    hEovH_C1->Fill(nTrk,EovH);
    hEovH_r1->Fill(r,EovH);
    hR60_1->Fill(r60);
 }
 void Fill1i( int njet, double dR, double relPt ) {
     njets_dRMuJ->Fill( njet, dR );
     dRMuJ_RelPt->Fill( dR, relPt );
 }
 void Fill1j( double eta, double pt, double EovH) {
    gEta_Pt->Fill( eta, pt);
    gEovH->Fill(EovH);
 }
 void Fill1k( int nj, int nj20, int nj40, int njh35 ) {
    gNJets ->Fill(nj);  
    gNjet20->Fill(nj20);  
    gNjet40->Fill(nj40);  
    gNjeth35->Fill(njh35);  
 }
 void Fill1l( double et1, double et2, double et3, double et4, double calEt3, double m3 ) {
    Jet1Et->Fill( et1 );
    Jet2Et->Fill( et2 );
    Jet3Et->Fill( et3 );
    Jet4Et->Fill( et4 );
    thirdCalEt->Fill( calEt3 );
    m3_j3->Fill( m3, et3 );
 }
 void Fill1m( double njets, double eta ) {
    eta_njets->Fill( eta, njets );
 }
 void Fill1n( int njets, double eta, double dphi0, double dphi1, double dphi2, double dR0, double dR1, double dR2, double dR, double dF, double m3 ) {
    wEta_dPhi0->Fill( eta, dphi0 );
    wEta_dPhi1->Fill( eta, dphi1 );
    wEta_dPhi2->Fill( eta, dphi2 );
    njets_dPhi0->Fill( njets, dphi0 );
    njets_dPhi1->Fill( njets, dphi1 );
    njets_dPhi2->Fill( njets, dphi2 );
    njets_dR0->Fill( njets, dR0 );
    njets_dR1->Fill( njets, dR1 );
    njets_dR2->Fill( njets, dR2 );
    njets_dR->Fill( njets, dR );
    njets_m3->Fill( njets, m3 );
    m3_dR->Fill( m3, dR );
 } 
 void Fill1o( double njets, double dphi0, double dphi1, double dphi2 ) {
    muNjets_dPhi0->Fill( njets, dphi0 );
    muNjets_dPhi1->Fill( njets, dphi1 );
    muNjets_dPhi2->Fill( njets, dphi2 );
 } 
 void Fill1p( double njets, double eta ) {
    wEta_njets->Fill( eta, njets );
 }

 void Fill1s( int allJ, int selJ ) {
    allJ_selJ->Fill( allJ, selJ );
 }

 void Write() {
    hEt->Write();
    hEovH->Write();
    hNCont->Write();
    hNTk->Write();
    hR60->Write();
    hArea_pt->Write();
    hEovH_pt->Write();
    hEovH_r->Write();
    hEovH_A->Write();
    Et_Pt->Write();
    hemF->Write();

    hWmass->Write();
    hWp->Write();

    EovH1->Write();
    Cont_Trk1->Write();
    EovH2->Write();
    Cont_Trk2->Write();

    hNJets->Write();
    hNjet20->Write();
    hNjet30->Write();
    hNjet40->Write();
    dRjj->Write();
 
    hEovH_A1->Write();
    hEovH_N1->Write();
    hEovH_C1->Write();
    hEovH_r1->Write();
    hR60_1->Write();

    gEta_Pt->Write();
    gEovH->Write();
    gNJets->Write();
    gNjet20->Write();
    gNjet40->Write();
    gNjeth35->Write();

    Jet1Et->Write();
    Jet2Et->Write();
    Jet3Et->Write();
    Jet4Et->Write();
    thirdCalEt->Write();
    eta_njets->Write();
    wEta_njets->Write();

    wEta_dPhi0->Write();
    wEta_dPhi1->Write();
    wEta_dPhi2->Write();
    njets_dPhi0->Write();
    njets_dPhi1->Write();
    njets_dPhi2->Write();
    
    njets_dRMuJ->Write();
    dRMuJ_RelPt->Write();
    muNjets_dPhi0->Write();
    muNjets_dPhi1->Write();
    muNjets_dPhi2->Write();

    njets_dR0->Write();
    njets_dR1->Write();
    njets_dR2->Write();

    njets_dR->Write();
    m3_dR->Write();
    m3_j3->Write();
    njets_m3->Write();


    allJ_selJ->Write();
 }

  TH1F *hEt;
  TH1F *hEovH;
  TH1F *hNCont;
  TH1F *hNTk;
  TH1F *hR60;
  TH2F *hArea_pt;
  TH2F *hEovH_pt;
  TH2F *hEovH_r;
  TH2F *hEovH_A;
  TH2F *Et_Pt;
  TH2F *hemF;

  TH1F *hWmass;
  TH1F *hWp;

  TH1F *EovH1;
  TH2F *Cont_Trk1;
  TH1F *EovH2;
  TH2F *Cont_Trk2;

  TH1F *hNJets;
  TH1F *hNjet20;
  TH1F *hNjet30;
  TH1F *hNjet40;
  TH1F *dRjj;

  TH2F *hEovH_A1;
  TH2F *hEovH_N1;
  TH2F *hEovH_C1;
  TH2F *hEovH_r1;
  TH1F *hR60_1;

  TH2F *gEta_Pt;
  TH1F *gEovH;
  TH1F *gNJets;
  TH1F *gNjet20;
  TH1F *gNjet40;
  TH1F *gNjeth35;

  TH1F *Jet1Et;
  TH1F *Jet2Et;
  TH1F *Jet3Et;
  TH1F *Jet4Et;
  TH1F *thirdCalEt;
  TH2F *eta_njets;
  TH2F *wEta_njets;

  TH2F *wEta_dPhi0;
  TH2F *wEta_dPhi1;
  TH2F *wEta_dPhi2;
  TH2F *njets_dPhi0;
  TH2F *njets_dPhi1;
  TH2F *njets_dPhi2;

  TH2F *njets_dRMuJ;
  TH2F *dRMuJ_RelPt;
  TH2F *muNjets_dPhi0;
  TH2F *muNjets_dPhi1;
  TH2F *muNjets_dPhi2;

  TH2F *njets_dR0;
  TH2F *njets_dR1;
  TH2F *njets_dR2;
  TH2F *njets_dR;
  TH2F *m3_dR;
  TH2F *m3_j3;
  TH2F *njets_m3;

  TH2F *allJ_selJ;

 TString name;

};

class HTOP2 {
public:
 
 HTOP2(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

   // reconstructed objects masses
    hemfMET   = new TH2F("hemfMET"," emfCalo(Y), emfPAT(X)", 150, -0.2, 1.3, 150, -0.2, 1.3);
    hSumEt_MET= new TH2F("hSumEt_MET","sumEt vs MET ", 1500, 0,1500, 500, 0,500);

    hResPAT = new TH1F("ResPAT","PAT MET Resolution", 100, -5.,5.);
    hResCal = new TH1F("ResCal","Cal MET Resolution", 100, -5.,5.);
    hResCor = new TH1F("ResCor","Cor MET Resolution", 100, -5.,5.);
    hPhiPAT = new TH1F("PhiPAT","PAT Phi Resolution", 100, -5.,5.);
    hPhiCal = new TH1F("PhiCal","Cal Phi Resolution", 100, -5.,5.);
    hPhiCor = new TH1F("PhiCor","Cor Phi Resolution", 100, -5.,5.);

    MET_dPhi0 = new TH2F("MET_dPhi0","MET vs dPhi(Mu,MET) > 3j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi1 = new TH2F("MET_dPhi1","MET vs dPhi(Mu,MET) 4j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi2 = new TH2F("MET_dPhi2","MET vs dPhi(Mu,MET) 5j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi3 = new TH2F("MET_dPhi3","MET vs dPhi(Mu,MET) 6j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi4 = new TH2F("MET_dPhi4","MET vs dPhi(Mu,MET) > 6j", 500, 0,500, 32, -0.05,3.15);

    MET_muPt = new TH2F("MET_muPt"," MET vs IsoMuon Pt ",  400, 0,400, 400, 0,400);
    dPhi_METJ1  = new TH2F("dPhi_METJ1"," MET vs 1st Jet ", 400, 0,400,32, -0.05,3.15);
    dPhi_METJ12 = new TH2F("dPhi_METJ12"," MET vs 1st + 2nd Jet ", 400, 0,400,32, -0.05,3.15);
 } 

 HTOP2(TString name_, TFile* file) {
    name=name_;

    hemfMET     = (TH2F *) file->Get("hemfMET");
    hSumEt_MET  = (TH2F *) file->Get("hSumEt_MET");
    hResPAT     = (TH1F *) file->Get("hResPAT");
    hResCal     = (TH1F *) file->Get("hResCal");
    hResCor     = (TH1F *) file->Get("hResCor");
    hPhiPAT     = (TH1F *) file->Get("hPhiPAT");
    hPhiCal     = (TH1F *) file->Get("hPhiCal");
    hPhiCor     = (TH1F *) file->Get("hPhiCor");

    MET_dPhi0   = (TH2F *) file->Get("MET_dPhi0");
    MET_dPhi1   = (TH2F *) file->Get("MET_dPhi1");
    MET_dPhi2   = (TH2F *) file->Get("MET_dPhi2");
    MET_dPhi3   = (TH2F *) file->Get("MET_dPhi3");
    MET_dPhi4   = (TH2F *) file->Get("MET_dPhi4");

    MET_muPt    = (TH2F *) file->Get("MET_muPt");
    dPhi_METJ1  = (TH2F *) file->Get("dPhi_METJ1");
    dPhi_METJ12 = (TH2F *) file->Get("dPhi_METJ12");
 }

 /// Destructor
 virtual ~HTOP2() {

    delete hemfMET;
    delete hSumEt_MET;
    delete hResPAT;
    delete hResCal;
    delete hResCor;
    delete hPhiPAT;
    delete hPhiCal;
    delete hPhiCor;

    delete MET_dPhi0;
    delete MET_dPhi1;
    delete MET_dPhi2;
    delete MET_dPhi3;
    delete MET_dPhi4;

    delete MET_muPt;
    delete dPhi_METJ1;
    delete dPhi_METJ12;
 }

 void Fill2a(double met, double emf, double emfCalo, double sumEt){

    hemfMET->Fill(emf, emfCalo);
    hSumEt_MET->Fill(sumEt, met);
 }

 void Fill2b(double respat, double rescal, double rescor, double phipat, double phical, double phicor){
    hResPAT->Fill(respat);
    hResCal->Fill(rescal);
    hResCor->Fill(rescor);
    hPhiPAT->Fill(phipat);
    hPhiCal->Fill(phical);
    hPhiCor->Fill(phicor);
 }

 void Fill2c0( double met, double dPhi) {
    MET_dPhi0->Fill(met, dPhi);
 }
 void Fill2c1( double met, double dPhi) {
    MET_dPhi1->Fill(met, dPhi);
 }
 void Fill2c2( double met, double dPhi) {
    MET_dPhi2->Fill(met, dPhi);
 }
 void Fill2c3( double met, double dPhi) {
    MET_dPhi3->Fill(met, dPhi);
 }
 void Fill2c4( double met, double dPhi) {
    MET_dPhi4->Fill(met, dPhi);
 }
 void Fill2d(double met, double muPt) {
    MET_muPt->Fill(met, muPt);
 }
 void Fill2e(double met, double df1, double df12) {
    dPhi_METJ1->Fill(met, df1);
    dPhi_METJ12->Fill(met, df12);
 }

 void Write() {

    hemfMET->Write();
    hSumEt_MET->Write();

    hResPAT->Write();
    hResCal->Write();
    hResCor->Write();
    hPhiPAT->Write();
    hPhiCal->Write();
    hPhiCor->Write();

    MET_dPhi0->Write();
    MET_dPhi1->Write();
    MET_dPhi2->Write();
    MET_dPhi3->Write();
    MET_dPhi4->Write();

    MET_muPt->Write();
    dPhi_METJ1->Write();
    dPhi_METJ12->Write();
 }

  TH2F *hemfMET;
  TH2F *hSumEt_MET;

  TH1F *hResPAT;
  TH1F *hResCal;
  TH1F *hResCor;
  TH1F *hPhiPAT;
  TH1F *hPhiCal;
  TH1F *hPhiCor;

  TH2F *MET_dPhi0;
  TH2F *MET_dPhi1;
  TH2F *MET_dPhi2;
  TH2F *MET_dPhi3;
  TH2F *MET_dPhi4;

  TH2F *MET_muPt;
  TH2F *dPhi_METJ1;
  TH2F *dPhi_METJ12;

 TString name;

};

class HTOP3 {
public:
 
 HTOP3(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

   // reconstructed objects masses
    hPt   = new TH1F("hPt", " Pt distribution",500,0,500);
    hEta  = new TH1F("hEta"," eta distribution",71,-3.55,3.55);
    iPt   = new TH1F("iPt", " iso muon Pt distribution",500,0,500);
    iEta  = new TH1F("iEta"," iso muon eta distribution",71,-3.55,3.55);

    // parton level
    hEta1a  = new TH1F("hEta1a"," eta distribution for muon of semi Tt(parton) ",71,-3.55,3.55);
    hEta1b  = new TH1F("hEta1b"," eta distribution for muon of semi Tt(mcMatched)",71,-3.55,3.55);
    hEta1c  = new TH1F("hEta1c"," eta distribution for muon of semi Tt(recoUsed)",71,-3.55,3.55);
    hPt1a   = new TH1F("hPt1a", " Pt distribution for muon of semi Tt(parton)",30,0,300 );
    hPt1b   = new TH2F("hPt1b", " Pt vs Pt_Resolution for muon of semi Tt(mcMatched)",30,0,300, 31,-1.55,1.55 );
    hPt1c   = new TH2F("hPt1c", " Pt vs Pt_Resolution for muon of semi Tt(recoUsed)",30,0,300, 31,-1.55,1.55 );

    hIMu_caloE = new TH2F("hIMu_caloE","IsoMu caloE vs muon p ",800,0.,200., 800,0.,200. );
    hJMu_caloE = new TH2F("hJMu_caloE","JetMu caloE vs muon p ",800,0.,200., 800,0.,200. );
    hIMu_caloE_jE = new TH2F("hIMu_caloE_jE","IsoMu caloE vs jet E ",800,0.,200., 800,0.,200. );
    hJMu_caloE_jE = new TH2F("hJMu_caloE_jE","JetMu caloE vs jet E ",800,0.,200., 800,0.,200. );

   // Isolation

    Pt_emIso1    = new TH2F("Pt_emIso1", "pt vs iso emE for all",200,0.,200., 100, 0., 100.);
    Pt_hdIso1    = new TH2F("Pt_hdIso1", "pt vs iso hdE for all",200,0.,200., 100, 0., 100.);
    Pt_tkIso1    = new TH2F("Pt_tkIso1", "pt vs iso tkP for all",200,0.,200., 100, 0., 100.);
    Pt_CaloIso1  = new TH2F("Pt_CaloIso1","pt vs calo iso for all",200,0.,200., 100, 0., 100.);
    Pt_Iso1      = new TH2F("Pt_Iso1"  , "pt vs Iso value for all",200,0.,200., 200, -0.5, 1.5 );
    Eta_Iso1     = new TH2F("Eta_Iso1", "eta vs Iso for all ",61, -3.05,3.05 , 200, -0.5, 1.5);  
    Eta_CaloIso1 = new TH2F("Eta_CaloIso1", "eta vs CaloIso for all ",61, -3.05,3.05 , 200, -0.5, 1.5);  

    Pt_emIso2    = new TH2F("Pt_emIso2", "pt vs iso emE after cut",200,0.,200., 100, 0., 100.);
    Pt_hdIso2    = new TH2F("Pt_hdIso2", "pt vs iso hdE after cut",200,0.,200., 100, 0., 100.);
    Pt_tkIso2    = new TH2F("Pt_tkIso2", "pt vs iso tkP after cut",200,0.,200., 100, 0., 100.);
    Pt_CaloIso2  = new TH2F("Pt_CaloIso2","pt vs calo iso after cut",200,0.,200., 100, 0., 100.);
    Pt_Iso2      = new TH2F("Pt_Iso2"  , "pt vs Iso value after cut",200,0.,200., 200, -0.5, 1.5 );
    Eta_Iso2     = new TH2F("Eta_Iso2", "eta vs Iso after cut ",61, -3.05,3.05 , 200, -0.5, 1.5);  
    Eta_CaloIso2 = new TH2F("Eta_CaloIso2", "eta vs CaloIso after cut ",61, -3.05,3.05 , 200, -0.5, 1.5);  

    Pt_emIso    = new TH2F("Pt_emIso", "pt vs iso emE  MC Matched",200,0.,200., 100, 0., 100.);
    Pt_hdIso    = new TH2F("Pt_hdIso", "pt vs iso hdE  MC Matched",200,0.,200., 100, 0., 100.);
    Pt_tkIso    = new TH2F("Pt_tkIso", "pt vs iso tkP  MC Matched",200,0.,200., 100, 0., 100.);
    Pt_CaloIso  = new TH2F("Pt_CaloIso","pt vs calo iso  MC Matched",200,0.,200., 100, 0., 100.);
    Pt_Iso      = new TH2F("Pt_Iso"  , "pt vs Iso value  MC Matched",200,0.,200., 200, -0.5, 1.5 );
    Eta_Iso     = new TH2F("Eta_Iso", "eta vs Iso  MC Matched ",61, -3.05,3.05 , 200, -0.5, 1.5);  
    Eta_CaloIso = new TH2F("Eta_CaloIso", "eta vs CaloIso  MC Matched ",61, -3.05,3.05 , 200, -0.5, 1.5);  

    allMu_isoMu = new TH2F("allMu_isoMu", " N of muons vs N IsoMuon in a event ", 21, -0.5, 20.5 , 21, -0.5, 20.5);

    glb_trk_Pt = new TH2F("glb_trk_Pt"," global Pt vs tracker Pt ", 200,0.,200., 200, 0., 200. );
    dPt_eta    = new TH2F("dPt_eta", " eta vs Glb Pt -trk Pt ", 61, -3.05,3.05 , 201, -1.005, 1.005 );

 } 

 HTOP3(TString name_, TFile* file) {
    name=name_;

    hPt    = (TH1F *) file->Get("hPt");
    hEta   = (TH1F *) file->Get("hEta");
    iPt    = (TH1F *) file->Get("iPt");
    iEta   = (TH1F *) file->Get("iEta");

    hEta1a = (TH1F *) file->Get("hEta1a");
    hEta1b = (TH1F *) file->Get("hEta1b");
    hEta1c = (TH1F *) file->Get("hEta1c");
    hPt1a  = (TH1F *) file->Get("hPt1a");
    hPt1b  = (TH2F *) file->Get("hPt1b");
    hPt1c  = (TH2F *) file->Get("hPt1c");

    Pt_emIso1    = (TH2F *) file->Get("Pt_emIso1");
    Pt_hdIso1    = (TH2F *) file->Get("Pt_hdIso1");
    Pt_tkIso1    = (TH2F *) file->Get("Pt_tkIso1");
    Pt_CaloIso1  = (TH2F *) file->Get("Pt_CaloIso1");
    Pt_Iso1      = (TH2F *) file->Get("Pt_Iso1");
    Eta_Iso1     = (TH2F *) file->Get("Eta_Iso1");
    Eta_CaloIso1 = (TH2F *) file->Get("Eta_CaloIso1");

    Pt_emIso2    = (TH2F *) file->Get("Pt_emIso2");
    Pt_hdIso2    = (TH2F *) file->Get("Pt_hdIso2");
    Pt_tkIso2    = (TH2F *) file->Get("Pt_tkIso2");
    Pt_CaloIso2  = (TH2F *) file->Get("Pt_CaloIso2");
    Pt_Iso2      = (TH2F *) file->Get("Pt_Iso2");
    Eta_Iso2     = (TH2F *) file->Get("Eta_Iso2");
    Eta_CaloIso2 = (TH2F *) file->Get("Eta_CaloIso2");

    Pt_emIso    = (TH2F *) file->Get("Pt_emIso");
    Pt_hdIso    = (TH2F *) file->Get("Pt_hdIso");
    Pt_tkIso    = (TH2F *) file->Get("Pt_tkIso");
    Pt_CaloIso  = (TH2F *) file->Get("Pt_CaloIso");
    Pt_Iso      = (TH2F *) file->Get("Pt_Iso");
    Eta_Iso     = (TH2F *) file->Get("Eta_Iso");
    Eta_CaloIso = (TH2F *) file->Get("Eta_CaloIso");

    hIMu_caloE = (TH2F *) file->Get("hIMu_caloE");
    hJMu_caloE = (TH2F *) file->Get("hJMu_caloE");
    hIMu_caloE_jE = (TH2F *) file->Get("hIMu_caloE_jE");
    hJMu_caloE_jE = (TH2F *) file->Get("hJMu_caloE_jE");

    allMu_isoMu = (TH2F *) file->Get("allMu_isoMu");

    glb_trk_Pt = (TH2F *) file->Get("glb_trk_Pt");
    dPt_eta    = (TH2F *) file->Get("dPt_eta");
 }

 /// Destructor
 virtual ~HTOP3() {

    delete hPt;
    delete hEta;
    delete iPt;
    delete iEta;

    delete hEta1a;
    delete hEta1b;
    delete hEta1c;
    delete hPt1a;
    delete hPt1b;
    delete hPt1c;
 
    delete Pt_emIso1;
    delete Pt_hdIso1;
    delete Pt_tkIso1;
    delete Pt_CaloIso1;
    delete Pt_Iso1;
    delete Eta_Iso1;
    delete Eta_CaloIso1;

    delete Pt_emIso2;
    delete Pt_hdIso2;
    delete Pt_tkIso2;
    delete Pt_CaloIso2;
    delete Pt_Iso2;
    delete Eta_Iso2;
    delete Eta_CaloIso2;

    delete Pt_emIso;
    delete Pt_hdIso;
    delete Pt_tkIso;
    delete Pt_CaloIso;
    delete Pt_Iso;
    delete Eta_Iso;
    delete Eta_CaloIso;

    delete hIMu_caloE;
    delete hJMu_caloE;
    delete hIMu_caloE_jE;
    delete hJMu_caloE_jE;

    delete allMu_isoMu;
    delete glb_trk_Pt;
    delete dPt_eta;
 }

 void Fill3a(double Pt,double Eta ){
    hPt->Fill(Pt);
    hEta->Fill(Eta);
 }
 
 void Fill3b(double pt, double eta, double emE, double hdE, double tkP, double sumE, double Iso ){

    Pt_emIso1->Fill(pt, emE);
    Pt_hdIso1->Fill(pt, hdE);
    Pt_tkIso1->Fill(pt, tkP);
    Pt_CaloIso1->Fill(pt, sumE);
    Pt_Iso1->Fill(pt, Iso);
    Eta_Iso1->Fill(eta, Iso);
    Eta_CaloIso1->Fill(eta, sumE);

 }
 
 void Fill3c(double pt, double eta, double emE, double hdE, double tkP, double sumE, double Iso ){

    Pt_emIso2->Fill(pt, emE);
    Pt_hdIso2->Fill(pt, hdE);
    Pt_tkIso2->Fill(pt, tkP);
    Pt_CaloIso2->Fill(pt, sumE);
    Pt_Iso2->Fill(pt, Iso);
    Eta_Iso2->Fill(eta, Iso);
    Eta_CaloIso2->Fill(eta, sumE);

 }

 void Fill3d(double pt, double eta, double emE, double hdE, double tkP, double sumE, double Iso ){

    Pt_emIso->Fill(pt, emE);
    Pt_hdIso->Fill(pt, hdE);
    Pt_tkIso->Fill(pt, tkP);
    Pt_CaloIso->Fill(pt, sumE);
    Pt_Iso->Fill(pt, Iso);
    Eta_Iso->Fill(eta, Iso);
    Eta_CaloIso->Fill(eta, sumE);

 }

 void Fill3e(int totalMu, int isoMu ) {
    allMu_isoMu->Fill( totalMu, isoMu );
 }

 void Fill3f(double eta, double pt) {
    hEta1a->Fill(eta);
    hPt1a->Fill(pt);
 }
 void Fill3g(double eta, double pt, double ptRes) {
    hEta1b->Fill(eta);
    hPt1b->Fill(pt, ptRes);
 }
 void Fill3h(double eta, double pt, double ptRes) {
    hEta1c->Fill(eta);
    hPt1c->Fill(pt, ptRes);
 }
 void Fill3i(double Pt,double Eta ){
    iPt->Fill(Pt);
    iEta->Fill(Eta);
 }
 void Fill3j(double isoMuE, double mu_p, double jet_E ) {
    hIMu_caloE->Fill(isoMuE, mu_p);
    hIMu_caloE_jE->Fill(isoMuE, jet_E);
 }
 void Fill3k(double jetMuE, double mu_p, double jet_E ) {
    hJMu_caloE->Fill(jetMuE, mu_p);
    hJMu_caloE_jE->Fill(jetMuE, jet_E);
 }
 void Fill3l(double glbPt, double trkPt, double dPt, double eta ) {
    glb_trk_Pt->Fill(glbPt, trkPt);
    dPt_eta->Fill( eta, dPt);
 }

 void Write() {

    hPt->Write();
    hEta->Write();
    iPt->Write();
    iEta->Write();
    hEta1a->Write();
    hEta1b->Write();
    hEta1c->Write();
    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
 
    Pt_emIso1->Write();
    Pt_hdIso1->Write();
    Pt_tkIso1->Write();
    Pt_CaloIso1->Write();
    Pt_Iso1->Write();
    Eta_Iso1->Write();
    Eta_CaloIso1->Write();

    Pt_emIso2->Write();
    Pt_hdIso2->Write();
    Pt_tkIso2->Write();
    Pt_CaloIso2->Write();
    Pt_Iso2->Write();
    Eta_Iso2->Write();
    Eta_CaloIso2->Write();

    Pt_emIso->Write();
    Pt_hdIso->Write();
    Pt_tkIso->Write();
    Pt_CaloIso->Write();
    Pt_Iso->Write();
    Eta_Iso->Write();
    Eta_CaloIso->Write();

    hIMu_caloE->Write();
    hJMu_caloE->Write();
    hIMu_caloE_jE->Write();
    hJMu_caloE_jE->Write();

    allMu_isoMu->Write();
    glb_trk_Pt->Write();
    dPt_eta->Write();
 }

  TH1F *hPt;
  TH1F *hEta;
  TH1F *iPt;
  TH1F *iEta;
  TH1F *hEta1a;
  TH1F *hEta1b;
  TH1F *hEta1c;
  TH1F *hPt1a;
  TH2F *hPt1b;
  TH2F *hPt1c;

  TH2F *Pt_emIso1;
  TH2F *Pt_hdIso1;
  TH2F *Pt_tkIso1;
  TH2F *Pt_CaloIso1;
  TH2F *Pt_Iso1;
  TH2F *Eta_Iso1;
  TH2F *Eta_CaloIso1;

  TH2F *Pt_emIso2;
  TH2F *Pt_hdIso2;
  TH2F *Pt_tkIso2;
  TH2F *Pt_CaloIso2;
  TH2F *Pt_Iso2;
  TH2F *Eta_Iso2;
  TH2F *Eta_CaloIso2;

  TH2F *Pt_emIso;
  TH2F *Pt_hdIso;
  TH2F *Pt_tkIso;
  TH2F *Pt_CaloIso;
  TH2F *Pt_Iso;
  TH2F *Eta_Iso;
  TH2F *Eta_CaloIso;

  TH2F *hIMu_caloE;
  TH2F *hJMu_caloE;
  TH2F *hIMu_caloE_jE;
  TH2F *hJMu_caloE_jE;

  TH2F *allMu_isoMu;

  TH2F *glb_trk_Pt;
  TH2F *dPt_eta;

 TString name;

};

class HTOP4 {
public:
 
 HTOP4(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

   // reconstructed objects masses
    hPt   = new TH1F("hPt", " Pt distribution",500,0,500);
    hEta  = new TH1F("hEta"," eta distribution",59,-2.95,2.95);

    // parton level
    hEta1a = new TH1F("hEta1a"," eta distribution for electron of semi Tt(parton) ",71,-3.55,3.55);
    hEta1b = new TH1F("hEta1b"," eta distribution for electron of semi Tt(mcMatched)",71,-3.55,3.55);
    hEta1c = new TH1F("hEta1c"," eta distribution for electron of semi Tt(recoUsed)",71,-3.55,3.55);
    hPt1a  = new TH1F("hPt1a", " Pt distribution for e of semi Tt(parton)",30,0,300 );
    hPt1b  = new TH2F("hPt1b", " Pt vs Pt_Resolution for e of semi Tt(mcMatched)",30,0,300, 31,-1.55,1.55 );
    hPt1c  = new TH2F("hPt1c", " Pt vs Pt_Resolution for e of semi Tt(recoUsed)",30,0,300, 31,-1.55,1.55 );
   // Isolation
    hEta_Iso1= new TH2F("hEta_Iso1", "eta vs Iso for all ",61, -3.05,3.05 , 200, -0.5, 1.5);  
    Pt_emIso1 = new TH2F("Pt_emIso1", "pt vs iso emE for all",200,0.,200., 100, 0., 100.);
    Pt_hdIso1 = new TH2F("Pt_hdIso1", "pt vs iso hdE for all",200,0.,200., 100, 0., 100.);
    Pt_tkIso1 = new TH2F("Pt_tkIso1", "pt vs iso tkP for all",200,0.,200., 100, 0., 100.);
    hPt_sum1 = new TH2F("hPt_sum1","pt vs iso caloE for all",200,0.,200., 100, 0., 100.);
    hPt_Cal1 = new TH2F("hPt_Cal1","pt vs # of CaloDeposit in isoR3 for all ",200,0.,200., 25, -0.5, 24.5 );
    hPt_Trk1 = new TH2F("hPt_Trk1","pt vs # of Track in isoR3 for all ",200,0.,200., 25, -0.5, 24.5 );
    hPt_Iso1 = new TH2F("hPt_Iso1","pt vs Iso value for all",200,0.,200., 200, -0.5, 1.5 );
    Iso_EovP1 = new TH2F("Iso_EovP1","RelIso vs emE/P for all",55,-0.05, 1.05, 150, -0.25, 1.25 );
    Iso_HovE1 = new TH2F("Iso_HovE1","RelIso vs H/E  for all ",55,-0.05, 1.05, 60, -0.05, 0.25 );

    Pt_emIso2 = new TH2F("Pt_emIso2", "pt vs iso emE after cut",200,0.,200., 100, 0., 100.);
    Pt_hdIso2 = new TH2F("Pt_hdIso2", "pt vs iso hdE after cut",200,0.,200., 100, 0., 100.);
    Pt_tkIso2 = new TH2F("Pt_tkIso2", "pt vs iso tkP for all",200,0.,200., 100, 0., 100.);
    hPt_sum2 = new TH2F("hPt_sum2","pt vs iso caloE after cut",200,0.,200., 100, 0., 100.);
    hPt_Cal2 = new TH2F("hPt_Cal2","pt vs # of CaloDeposit in isoR3 after cut",200,0.,200., 25, -0.5, 24.5 );
    hPt_Trk2 = new TH2F("hPt_Trk2","pt vs # of Track in isoR3 after cut",200,0.,200., 25, -0.5, 24.5 );
    hPt_Iso2 = new TH2F("hPt_Iso2","pt vs Iso value after cut",200,0.,200., 200, -0.5, 1.5 );
    Iso_EovP2 = new TH2F("Iso_EovP2","RelIso vs emE/P after cut",55,-0.05, 1.05, 150, -0.25, 1.25 );
    Iso_HovE2 = new TH2F("Iso_HovE2","RelIso vs H/E after cut ",55,-0.05, 1.05, 60, -0.05, 0.25 );

   // MC matching electrons
    hEta_Iso = new TH2F("hEta_Iso", "eta vs Iso - MC matched ",61, -3.05,3.05, 200, -0.5, 1.5);  
    hPt_nTrk = new TH2F("hPt_nTrk","pt vs # of track in isoR3 ",200,0.,200., 25, -0.5, 24.5 );
    hPt_nCal = new TH2F("hPt_nCal","pt vs # of CaloDeposit in isoR3 ",200,0.,200., 25, -0.5, 24.5 );
    hPt_Iso  = new TH2F("hPt_Iso","pt vs Iso value - MC matched ",200,0.,200., 200, -0.5, 1.5);
    Pt_tkIso = new TH2F("Pt_tkIso", "pt vs trkIso ", 200, 0., 200. , 100, 0., 100.);  
    Pt_emIso = new TH2F("Pt_emIso","pt vs ecalIso ", 200, 0., 200. , 100, 0., 100.);  
    Pt_hdIso = new TH2F("Pt_hdIso","pt vs hcalIso ", 200, 0., 200. , 100, 0., 100.);  
    Iso_EovP = new TH2F("Iso_EovP","RelIso vs emE/P ",55,-0.05, 1.05, 150, -0.25, 1.25 );
    Iso_HovE = new TH2F("Iso_HovE","RelIso vs H/E   ",55,-0.05, 1.05, 60, -0.05, 0.25 );

    allEl_isoEl = new TH2F("allEl_isoEl", " N of Electrons vs N IsoElectrons in a event ", 21, -0.5, 20.5 , 21, -0.5, 20.5);
    isoEleCut = new TH2F("isoEleCut","N of isoEle vs N of SelectedJets  after isoMu =1 cut ", 21, -0.5, 20.5, 21, -0.5, 20.5 );
 } 

 HTOP4(TString name_, TFile* file) {
    name=name_;

    hPt    = (TH1F *) file->Get("hPt");
    hEta   = (TH1F *) file->Get("hEta");
    hEta1a = (TH1F *) file->Get("hEta1a");
    hEta1b = (TH1F *) file->Get("hEta1b");
    hEta1c = (TH1F *) file->Get("hEta1c");
    hPt1a  = (TH1F *) file->Get("hPt1a");
    hPt1b  = (TH2F *) file->Get("hPt1b");
    hPt1c  = (TH2F *) file->Get("hPt1c");

    Pt_emIso1 = (TH2F *) file->Get("Pt_emIso1");
    Pt_hdIso1 = (TH2F *) file->Get("Pt_hdIso1");
    Pt_tkIso1 = (TH2F *) file->Get("Pt_tkIso1");
    hPt_sum1 = (TH2F *) file->Get("hPt_sum1");
    hPt_Cal1 = (TH2F *) file->Get("hPt_Cal1");
    hPt_Trk1 = (TH2F *) file->Get("hPt_Trk1");
    hPt_Iso1 = (TH2F *) file->Get("hPt_Iso1");
    hEta_Iso1 = (TH2F *) file->Get("hEta_Iso1");
    Iso_EovP1 = (TH2F *) file->Get("Iso_EovP1");
    Iso_HovE1 = (TH2F *) file->Get("Iso_HovE1");

    Pt_emIso2 = (TH2F *) file->Get("Pt_emIso2");
    Pt_hdIso2 = (TH2F *) file->Get("Pt_hdIso2");
    Pt_tkIso2 = (TH2F *) file->Get("Pt_tkIso2");
    hPt_sum2 = (TH2F *) file->Get("hPt_sum2");
    hPt_Cal2 = (TH2F *) file->Get("hPt_Cal2");
    hPt_Trk2 = (TH2F *) file->Get("hPt_Trk2");
    hPt_Iso2 = (TH2F *) file->Get("hPt_Iso2");
    Iso_EovP2 = (TH2F *) file->Get("Iso_EovP2");
    Iso_HovE2 = (TH2F *) file->Get("Iso_HovE2");

    Iso_EovP = (TH2F *) file->Get("Iso_EovP");
    Iso_HovE = (TH2F *) file->Get("Iso_HovE");
    hPt_nTrk = (TH2F *) file->Get("hPt_nTrk");
    hPt_nCal = (TH2F *) file->Get("hPt_nCal");
    hPt_Iso  = (TH2F *) file->Get("hPt_Iso");
    hEta_Iso = (TH2F *) file->Get("hEta_Iso");
    Pt_tkIso  = (TH2F *) file->Get("Pt_tkIso");
    Pt_emIso = (TH2F *) file->Get("Pt_emIso");
    Pt_hdIso = (TH2F *) file->Get("Pt_hdIso");

    allEl_isoEl =(TH2F *) file->Get("allEl_isoEl");
    isoEleCut   =(TH2F *) file->Get("isoEleCut");
 }

 /// Destructor
 virtual ~HTOP4() {

    delete hPt;
    delete hEta;
    delete hEta1a;
    delete hEta1b;
    delete hEta1c;
    delete hPt1a;
    delete hPt1b;
    delete hPt1c;
 
    delete Pt_emIso1;
    delete Pt_hdIso1;
    delete Pt_tkIso1;
    delete hPt_sum1;
    delete hPt_Cal1;
    delete hPt_Trk1;
    delete hPt_Iso1;
    delete hEta_Iso1;
    delete Iso_EovP1;
    delete Iso_HovE1;

    delete Pt_emIso2;
    delete Pt_hdIso2;
    delete Pt_tkIso2;
    delete hPt_sum2;
    delete hPt_Cal2;
    delete hPt_Trk2;
    delete hPt_Iso2;
    delete Iso_EovP2;
    delete Iso_HovE2;

    delete Iso_EovP;
    delete Iso_HovE;
    delete hPt_nTrk;
    delete hPt_nCal;
    delete hPt_Iso;
    delete hEta_Iso;
    delete Pt_tkIso;
    delete Pt_emIso;
    delete Pt_hdIso;

    delete allEl_isoEl;
    delete isoEleCut;
 }

 void Fill4a(double Pt,double Eta ){
    hPt->Fill(Pt);
    hEta->Fill(Eta);
 }
 
 void Fill4b(double pt, double eta, double emE, double hdE, double tkP, double sumE, int nCal3, int ntrack3, 
             double Iso, double HovE, double EovP ){

    Pt_emIso1->Fill(pt, emE);
    Pt_hdIso1->Fill(pt, hdE);
    Pt_tkIso1->Fill(pt, tkP);
    hPt_sum1->Fill(pt, sumE);
    hPt_Cal1->Fill(pt, nCal3);
    hPt_Trk1->Fill(pt, ntrack3);
    hPt_Iso1->Fill(pt, Iso);
    hEta_Iso1->Fill(eta, Iso);
    Iso_EovP1->Fill(Iso, EovP);
    Iso_HovE1->Fill(Iso, HovE);
 }
 
 void Fill4c(double pt, double emE, double hdE, double tkP, double sumE, int nCal3, int ntrack3, 
             double Iso, double HovE, double EovP){

    Pt_emIso2->Fill(pt, emE);
    Pt_hdIso2->Fill(pt, hdE);
    Pt_tkIso2->Fill(pt, tkP);
    hPt_sum2->Fill(pt, sumE);
    hPt_Cal2->Fill(pt, nCal3);
    hPt_Trk2->Fill(pt, ntrack3);
    hPt_Iso2->Fill(pt, Iso);
    Iso_EovP2->Fill(Iso, EovP);
    Iso_HovE2->Fill(Iso, HovE);
 }

 void Fill4d(double pt, double eta, double EovP, double HovE, int nTrk, int nCal, float trkIso, float ecalIso, float hcalIso, double Iso ) {
    hEta_Iso->Fill(eta, Iso);
    hPt_nTrk->Fill(pt, nTrk);
    hPt_nCal->Fill(pt, nCal);
    hPt_Iso->Fill(pt, Iso);
    Pt_tkIso->Fill(pt, trkIso);
    Pt_emIso->Fill(pt, ecalIso);
    Pt_hdIso->Fill(pt, hcalIso);
    Iso_EovP->Fill(Iso, EovP);
    Iso_HovE->Fill(Iso, HovE);
 }
 void Fill4e( int allEl, int isoEl ){
    allEl_isoEl->Fill(allEl, isoEl);
 }

 void Fill4f(double eta, double pt) {
    hEta1a->Fill(eta);
    hPt1a->Fill(pt);
 }
 void Fill4g(double eta, double pt, double ptRes){
    hEta1b->Fill(eta);
    hPt1b->Fill(pt, ptRes);
 }
 void Fill4h(double eta, double pt, double ptRes) {
    hEta1c->Fill(eta);
    hPt1c->Fill(pt, ptRes);
 }
 void Fill4i( int NIsoEle, int NSelJet ){
    isoEleCut->Fill( NIsoEle, NSelJet );
 }

 void Write() {

    hPt->Write();
    hEta->Write();
    hEta1a->Write();
    hEta1b->Write();
    hEta1c->Write();
    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
 
    Pt_emIso1->Write();
    Pt_hdIso1->Write();
    Pt_tkIso1->Write();
    hPt_sum1->Write();
    hPt_Cal1->Write();
    hPt_Trk1->Write();
    hPt_Iso1->Write();
    hEta_Iso1->Write();
    Iso_EovP1->Write();
    Iso_HovE1->Write();

    Pt_emIso2->Write();
    Pt_hdIso2->Write();
    Pt_tkIso2->Write();
    hPt_sum2->Write();
    hPt_Cal2->Write();
    hPt_Trk2->Write();
    hPt_Iso2->Write();
    Iso_EovP2->Write();
    Iso_HovE2->Write();

    Iso_EovP->Write();
    Iso_HovE->Write();
    hPt_nTrk->Write();
    hPt_nCal->Write();
    hPt_Iso->Write();
    hEta_Iso->Write();
    Pt_tkIso->Write();
    Pt_emIso->Write();
    Pt_hdIso->Write();

    allEl_isoEl->Write();
    isoEleCut->Write();
 }

  TH1F *hPt;
  TH1F *hEta;
  TH1F *hEta1a;
  TH1F *hEta1b;
  TH1F *hEta1c;
  TH1F *hPt1a;
  TH2F *hPt1b;
  TH2F *hPt1c;

  TH2F *Pt_emIso1;
  TH2F *Pt_hdIso1;
  TH2F *Pt_tkIso1;
  TH2F *hPt_sum1;
  TH2F *hPt_Cal1;
  TH2F *hPt_Trk1;
  TH2F *hPt_Iso1;
  TH2F *hEta_Iso1;
  TH2F *Iso_EovP1;
  TH2F *Iso_HovE1;

  TH2F *Pt_emIso2;
  TH2F *Pt_hdIso2;
  TH2F *Pt_tkIso2;
  TH2F *hPt_sum2;
  TH2F *hPt_Cal2;
  TH2F *hPt_Trk2;
  TH2F *hPt_Iso2;
  TH2F *Iso_EovP2;
  TH2F *Iso_HovE2;

  TH2F *Iso_EovP;
  TH2F *Iso_HovE;
  TH2F *hPt_nTrk;
  TH2F *hPt_nCal;
  TH2F *hPt_Iso;
  TH2F *hEta_Iso;
  TH2F *Pt_tkIso;
  TH2F *Pt_emIso;
  TH2F *Pt_hdIso;

  TH2F *allEl_isoEl;
  TH2F *isoEleCut;

 TString name;

};

class HTOP5 {
public:
 
 HTOP5(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hPt   = new TH1F("hPt", " Pt distribution",500,0,500);
    hEta  = new TH1F("hEta"," eta distribution",61,-3.05,3.05);

 } 

 HTOP5(TString name_, TFile* file) {
    name=name_;

    hPt    = (TH1F *) file->Get("hPt");
    hEta   = (TH1F *) file->Get("hEta");

 }

 /// Destructor
 virtual ~HTOP5() {
    delete hPt;
    delete hEta;

 }

 void Fill5a(double Pt,double Eta ){
    hPt->Fill(Pt);
    hEta->Fill(Eta);
 }


 void Write() {
    hPt->Write();
    hEta->Write();
 }

  TH1F *hPt;
  TH1F *hEta;

 TString name;

};

class HTOP6 {
public:
 
 HTOP6(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    dRb_dRw_lep = new TH2F("dRb_dRw_lep","leptonic T dR(mcB, recoB), dR(mcW, recoW) ", 200,-0.025,9.975, 200, -0.025, 9.975 );
    dRb_dRw_had = new TH2F("dRb_dRw_had","hadronic T dR(mcB, recoB), dR(mcW, recoW) ", 200,-0.025,9.975, 200, -0.025, 9.975 );
    dRtt        = new TH2F("dRtt"," dR(rec, mc)_lep, dR(rec, mc)_had", 200, -0.025, 9.975, 200, -0.025, 9.975 );

    dRbw_had = new TH2F("dRbw_had","hadronic T dR(mcb, mcW), dR( rcb, rcW )", 200, -0.025, 9.975, 200, -0.025, 9.975 );
    dRbw_lep = new TH2F("dRbw_lep","leptonic T dR(mcb, mcW), dR( rcb, rcW )", 200, -0.025, 9.975, 200, -0.025, 9.975 );

    dR_lepW    = new TH1F("dR_lepW"," dR(reco - gen) ",200,0.,10.);
    PtRes_lepW = new TH1F("PtRes_lepW"," Pt Res ",200,-1.005,0.995);
    MRes_lepW  = new TH1F("MRes_lepW"," Mass Res ",200,-1.005,0.995);

 } 

 HTOP6(TString name_, TFile* file) {
    name=name_;

    dRb_dRw_lep= (TH2F *) file->Get("dRb_dRw_lep");
    dRb_dRw_had= (TH2F *) file->Get("dRb_dRw_had");
    dRbw_had   = (TH2F *) file->Get("dRbw_had");
    dRbw_lep   = (TH2F *) file->Get("dRbw_lep");
    dRtt       = (TH2F *) file->Get("dRtt");

    dR_lepW    = (TH1F *) file->Get("dR_lepW");
    PtRes_lepW = (TH1F *) file->Get("PtRes_lepW");
    MRes_lepW  = (TH1F *) file->Get("MRes_lepW");

 }

 /// Destructor
 virtual ~HTOP6() {

    delete dRb_dRw_lep;
    delete dRb_dRw_had;
    delete dRbw_had;
    delete dRbw_lep;
    delete dRtt;

    delete dR_lepW;
    delete PtRes_lepW;
    delete MRes_lepW;

 }

 void Fill6a( double dRt0, double dRbj0, double dRw0, double dRt1, double dRbj1, double dRw1 ) {
    dRb_dRw_lep->Fill( dRbj0, dRw0);
    dRb_dRw_had->Fill( dRbj1, dRw1);
    dRtt->Fill(dRt0,dRt1);
 }
 void Fill6b( double dRbw_mcl, double dRbw_rcl, double dRbw_mch, double dRbw_rch ){
    dRbw_had->Fill( dRbw_mch, dRbw_rch );
    dRbw_lep->Fill( dRbw_mcl, dRbw_rcl );
 } 
 void Fill6d( double dR, double PtRes, double MRes ) {
    dR_lepW->Fill(dR);
    PtRes_lepW->Fill(PtRes);
    MRes_lepW->Fill(MRes);
 }


 void Write() {

    dRb_dRw_lep->Write();
    dRb_dRw_had->Write();
    dRbw_had->Write();
    dRbw_lep->Write();
    dRtt->Write();

    dR_lepW->Write();
    PtRes_lepW->Write();
    MRes_lepW->Write();
  
 }

  TH2F *dRb_dRw_lep;
  TH2F *dRb_dRw_had;
  TH2F *dRbw_had;
  TH2F *dRbw_lep;
  TH2F *dRtt;

  TH1F *dR_lepW;
  TH1F *PtRes_lepW;
  TH1F *MRes_lepW;

 TString name;

};

class HTOP7 {
public:
 
 HTOP7(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hRes_Pt   = new TH1F("hRes_Pt","Resolutio  of b jet pt ",400, -2., 2. ); 
    dR_RelPt_bmu = new TH2F("dR_RelPt_bmu","dRmin vs RelPt for b and  iso muon ",200, 0.,5., 100, 0.,200. );  
    hdRbWj    = new TH1F("hdRbWj","dR for b and  W matched jets ",300, 0.,15. );  
    hbDis_bCand  = new TH2F("hbDis_bCand","bCandidates jets JProb(x),TkCount(y)", 150,-0.1,1.4, 500,-20.,80.);
    hbDis_all   = new TH2F("hbDis_all","all jets bDis JProb(x),TkCount(y) ",150,-0.1,1.4, 500,-20.,80.);

    hNTk_mc  = new TH1F("hNTk_mc"," # of charge track of a bjet", 50, -1.5, 48.5);
    hemE     = new TH1F("hemE"," em Energy fraction ", 150, -0.2, 1.3);
    hR60_mc  = new TH2F("hR60_mc"," n60 vs nConstituents of a bjet", 20, -1.5, 18.5, 50, -0.5, 49.5);
    hR90_mc  = new TH2F("hR90_mc"," n90 vs nConstituents of a bjet", 20, -1.5, 18.5, 50, -0.5, 49.5);
    hArea_pt_mc= new TH2F("hArea_pt_mc","tower area vs pt of a bjet", 500, 0., 500., 200, 0., 1.);
    hEovH_p_mc = new TH2F("hEovH_p_mc","E/H of a bjet", 500, 0., 500., 500, -2., 23.);
    hEovH_n_mc = new TH2F("hEovH_n_mc","E/H vs. nConstituent of a bjet", 50, -0.5, 49.5, 500, -2., 23.);
    hEovH_A_mc = new TH2F("hEovH_A_mc"," E/H vs tower area of a  bjet", 200, 0., 1., 500, -2., 23.);
    hEovH_h_mc = new TH2F("hEovH_h_mc"," E/H vs eta of a bjet", 61, -3.05, 3.05, 500, -2., 23.);
    hPt_Eta_mc = new TH2F("hPt_Eta_mc", " eta vs pt (bjets)", 500, 0.,250, 71, -3.55, 3.55);

    gbJpt_h   = new TH2F("gbJpt_h"," gen b jet pt vs eta", 500,0,250, 71,-3.55,3.55);
    gEovH     = new TH1F("gEovH"," E/H for gen Jets", 500, -2., 23.);

    // parton level
    hPt1a   = new TH1F("hPt1a", " Pt distribution  for b of semi Tt(parton)",30,0,300 );
    hPt1b   = new TH2F("hPt1b", " Pt vs Pt_Resolution for b of semi Tt(mcMatched)",30,0,300, 31,-1.55,1.55 );
    hPt1c   = new TH2F("hPt1c", " Pt vs Pt_Resolution for b of semi Tt(recoUsed)",30,0,300, 31,-1.55,1.55 );
    // 
    dR_RelPtg = new TH2F("dR_RelPtg","dR vs RelPt for b and mu - gen ",200, 0.,5., 500, 0.,50. );  
    dR_RelPtr = new TH2F("dR_RelPtr","dR vs RelPt for b and mu - reco after matching ",200, 0.,5., 500, 0.,50. );  
 } 

 HTOP7(TString name_, TFile* file) {
    name=name_;

    hRes_Pt = (TH1F *) file->Get("hRes_Pt");
    dR_RelPt_bmu = (TH2F *) file->Get("dR_RelPt_bmu");
    hdRbWj = (TH1F *) file->Get("hdRbWj");
    hbDis_bCand =(TH2F *) file->Get("hbDis_bCand");
    hbDis_all   =(TH2F *) file->Get("hbDis_all");

    hNTk_mc = (TH1F *) file->Get("hNTK_mc");
    hemE    = (TH1F *) file->Get("hemE");
    hR60_mc = (TH2F *) file->Get("hR60_mc");
    hR90_mc = (TH2F *) file->Get("hR90_mc");
    hArea_pt_mc = (TH2F *) file->Get("hArea_pt_mc");
    hEovH_p_mc  = (TH2F *) file->Get("hEovH_p_mc");
    hEovH_n_mc  = (TH2F *) file->Get("hEovH_n_mc");
    hEovH_A_mc  = (TH2F *) file->Get("hEovH_A_mc");
    hEovH_h_mc  = (TH2F *) file->Get("hEovH_h_mc");
    hPt_Eta_mc  = (TH2F *) file->Get("hPt_Eta_mc");

    gbJpt_h = (TH2F *) file->Get("gbJpt_h");
    gEovH   = (TH1F *) file->Get("gEovH");

    hPt1a  = (TH1F *) file->Get("hPt1a");
    hPt1b  = (TH2F *) file->Get("hPt1b");
    hPt1c  = (TH2F *) file->Get("hPt1c");

    dR_RelPtg = (TH2F *) file->Get("dR_RelPtg");
    dR_RelPtr = (TH2F *) file->Get("dR_RelPtr");

 }

 /// Destructor
 virtual ~HTOP7() {

    delete hRes_Pt;
    delete dR_RelPt_bmu;
    delete hdRbWj;
    delete hbDis_bCand;
    delete hbDis_all;

    delete hNTk_mc;
    delete hR60_mc;
    delete hR90_mc;
    delete hArea_pt_mc;
    delete hEovH_p_mc;
    delete hEovH_n_mc;
    delete hEovH_A_mc;
    delete hEovH_h_mc;
    delete hemE;
    delete hPt_Eta_mc;

    delete gbJpt_h;
    delete gEovH;

    delete hPt1a;
    delete hPt1b;
    delete hPt1c;

    delete dR_RelPtg;
    delete dR_RelPtr;
 }

 void Fill7a(double pt, double eta, double EovH, int nCon, int n60, int n90, double area, size_t nTrk, double emE, double Res_pt, double bDisJProb, double bDisTkCount ){

    hNTk_mc->Fill(nTrk);
    hemE->Fill(emE);
    hR60_mc->Fill(n60, nCon);
    hR90_mc->Fill(n90, nCon);
    hArea_pt_mc->Fill(pt,area);
    hEovH_p_mc->Fill(pt,EovH);
    hEovH_n_mc->Fill(nCon,EovH);
    hEovH_A_mc->Fill(area,EovH);
    hEovH_h_mc->Fill(eta,EovH);
    hPt_Eta_mc->Fill(pt,eta);

    hRes_Pt->Fill(Res_pt);
    hbDis_bCand->Fill(bDisJProb,bDisTkCount );
 }
 void Fill7b( double dR_bmu , double RelPt ){
    dR_RelPt_bmu->Fill(dR_bmu, RelPt);
 }
 void Fill7c( double dR_bwj ){
    hdRbWj->Fill(dR_bwj);
 }
 void Fill7d( double bDisJProb, double bDisTkCount ){
    hbDis_all->Fill(bDisJProb,bDisTkCount );
 }

 void Fill7e( double dR, double RelPt, double gdR, double gRelPt) {
    dR_RelPtr->Fill(dR, RelPt);
    dR_RelPtg->Fill(gdR, gRelPt);
 }

 void Fill7f( double pt, double eta, double EovH) {
    gbJpt_h->Fill(pt, eta);
    gEovH->Fill(EovH);
 }

 void Fill7g( double pt ) {
    hPt1a->Fill(pt);
 }
 void Fill7h( double pt, double ptRes) {
    hPt1b->Fill(pt, ptRes);
 }
 void Fill7i( double pt, double ptRes) {
    hPt1c->Fill(pt, ptRes);
 }

 void Write() {

    hRes_Pt->Write();
    dR_RelPt_bmu->Write();
    hdRbWj->Write();
    hbDis_bCand->Write();
    hbDis_all->Write();

    hNTk_mc->Write();
    hemE->Write();
    hR60_mc->Write();
    hR90_mc->Write();
    hArea_pt_mc->Write();
    hEovH_p_mc->Write();
    hEovH_n_mc->Write();
    hEovH_A_mc->Write();
    hEovH_h_mc->Write();
    hPt_Eta_mc->Write();

    gbJpt_h->Write();
    gEovH->Write();

    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
  
    dR_RelPtg->Write();
    dR_RelPtr->Write();
 }

  TH1F *hRes_Pt;
  TH2F *dR_RelPt_bmu;
  TH1F *hdRbWj;
  TH2F *hbDis_bCand;
  TH2F *hbDis_all;

  TH1F *hNTk_mc;
  TH1F *hemE;
  TH2F *hR60_mc;
  TH2F *hR90_mc;
  TH2F *hArea_pt_mc;
  TH2F *hEovH_p_mc;
  TH2F *hEovH_n_mc;
  TH2F *hEovH_A_mc;
  TH2F *hEovH_h_mc;
  TH2F *hPt_Eta_mc;

  TH2F *gbJpt_h;
  TH1F *gEovH;

  TH1F *hPt1a;
  TH2F *hPt1b;
  TH2F *hPt1c;

  TH2F *dR_RelPtg;
  TH2F *dR_RelPtr;
 TString name;

};

class HTOP8 {
public:
 
 HTOP8(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hEta_Pt1 = new TH2F("hEta_Pt1", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hEta_Pt2 = new TH2F("hEta_Pt2", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hNJets   = new TH1F("hNJets"," # of selected W jets ",65,-0.5,64.5);
    hWp_mass = new TH2F("hWp_mass"," W mass vs W momentum from selected reco jets", 500,0.,500.,240,30.,150.);

    hRes_Pt    = new TH1F("hRes_Pt", " Pt(jet) - Pt(q) / Pt(q)   ",400,-2.,2.);

    hdRWjMu = new TH1F("hdRWjMu","dR for isoMuon and W matched jets ",200, -0.025,9.975 );  
    hdRWjj  = new TH2F("hdRWjj" ,"dR(Wq1, Wq2) , dR ( matched Wj1 , Wj2)  " ,200, -0.025,9.975 ,200, -0.025,9.975 );  

    hbDis_WCand  = new TH2F("hbDis_WCand","bCandidates jets JProb(x),TkCount(y)", 150,-0.1,1.4, 500,-20.,80.);

    gwJpt_h   = new TH2F("gwJpt_h"," gen w jet pt vs eta", 500,0,250, 71,-3.55,3.55);
    gEovH    = new TH1F("gEovH"," E/H for gen Jets", 500, -2., 23.);

    // parton level
    hPt1a   = new TH1F("hPt1a", "Pt distribution for Wjets of semi Tt(parton)",30,0,300 );
    hPt1b   = new TH2F("hPt1b", "Pt vs Pt_Resolution for Wjets of semi Tt(mcMatched)",30,0,300, 31,-1.55,1.55);
    hPt1c   = new TH2F("hPt1c", "Pt vs Pt_Resolution for Wjets of semi Tt(recoUsed)",30,0,300, 31,-1.55,1.55);
 } 

 HTOP8(TString name_, TFile* file) {
    name=name_;

    hEta_Pt1 = (TH2F *) file->Get("hEta_Pt1");
    hEta_Pt2 = (TH2F *) file->Get("hEta_Pt2");
    hNJets = (TH1F *) file->Get("hNJets");
    hWp_mass   = (TH2F *) file->Get("hWp_mass");

    hRes_Pt     = (TH1F *) file->Get("hRes_Pt");

    hdRWjMu = (TH1F *) file->Get("hdRWjMu");
    hdRWjj  = (TH2F *) file->Get("hdRWjj");

    hbDis_WCand =(TH2F *) file->Get("hbDis_WCand");

    gwJpt_h = (TH2F *) file->Get("gwJpt_h");
    gEovH   = (TH1F *) file->Get("gEovH");

    hPt1a  = (TH1F *) file->Get("hPt1a");
    hPt1b  = (TH2F *) file->Get("hPt1b");
    hPt1c  = (TH2F *) file->Get("hPt1c");
 }

 /// Destructor
 virtual ~HTOP8() {

    delete hEta_Pt1;
    delete hEta_Pt2;
    delete hNJets;
    delete hWp_mass;

    delete hRes_Pt;
    delete hdRWjMu;
    delete hdRWjj;

    delete hbDis_WCand;

    delete gwJpt_h;
    delete gEovH;

    delete hPt1a;
    delete hPt1b;
    delete hPt1c;
 }

 void Fill8a( double pt1, double pt2, double eta1, double eta2, int njets){
    hEta_Pt1->Fill(eta1, pt1);
    hEta_Pt2->Fill(eta2, pt2);
    hNJets->Fill(njets);
 }

 void Fill8b( double wmass, double wp ){
    hWp_mass->Fill(wp,wmass);
 }

 void Fill8c( double Res_Pt,  double bDisJProb, double bDisTkCount ){
    hRes_Pt->Fill(Res_Pt);
    hbDis_WCand->Fill(bDisJProb,bDisTkCount );
 }

 void Fill8d( double dR_WjMu ){
    hdRWjMu->Fill(dR_WjMu);
 }

 void Fill8f( double pt, double eta, double EovH) {
    gwJpt_h->Fill(pt, eta);
    gEovH->Fill(EovH);
 }

 void Fill8g( double pt ) {
    hPt1a->Fill(pt );
 }
 void Fill8h( double pt, double ptRes) {
    hPt1b->Fill(pt, ptRes);
 }
 void Fill8i( double pt, double ptRes) {
    hPt1c->Fill(pt, ptRes);
 }
 void Fill8j( double dRjj, double dRqq ) {
    hdRWjj->Fill(dRqq, dRjj);
 }

 void Write() {
    hEta_Pt1->Write();
    hEta_Pt2->Write();
    hNJets->Write();
    hWp_mass->Write();

    hRes_Pt->Write();
    hdRWjMu->Write();
    hdRWjj->Write();

    hbDis_WCand->Write();

    gwJpt_h->Write();
    gEovH->Write();

    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
 }

  TH2F *hEta_Pt1;
  TH2F *hEta_Pt2;
  TH1F *hNJets;
  TH2F *hWp_mass;

  TH1F *hRes_Pt;
  TH1F *hdRWjMu;
  TH2F *hdRWjj;

  TH2F *hbDis_WCand;
  TH2F *gwJpt_h;
  TH1F *gEovH;

  TH1F *hPt1a;
  TH2F *hPt1b;
  TH2F *hPt1c;

 TString name;

};


class HTOP9 {
public:
 
 HTOP9(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

   // reconstructed objects masses
    hMCtt     = new TH2F("hMCtt","MCMatched top mass ",200,0,400, 200,0,400);
    hMCtt1    = new TH2F("hMCtt1","MCMatched top mass passing evt selection",200,0,400, 200,0,400);
    hRecott0   = new TH2F("hRecott0","Reco top mass w/ > 3j ", 200,0,400, 200,0,400);
    hRecott1   = new TH2F("hRecott1","Reco top mass w/ 4j", 200,0,400, 200,0,400);
    hRecott2   = new TH2F("hRecott2","Reco top mass w/ 5j", 200,0,400, 200,0,400);
    hRecott3   = new TH2F("hRecott3","Reco top mass w/ 6j", 200,0,400, 200,0,400);
    hRecott4   = new TH2F("hRecott4","Reco top mass w/ >= 7j", 200,0,400, 200,0,400);

    PtMCtt      = new TH1F("PtMCtt",   "Pt of MCMatched tt", 200,0,400);
    PtRecott0   = new TH1F("PtRecott0","Pt of Reco top mass w/ > 3j", 200,0,400);
    PtRecott1   = new TH1F("PtRecott1","Pt of Reco top mass w/ 4j", 200,0,400);
    PtRecott2   = new TH1F("PtRecott2","Pt of Reco top mass w/ 5j", 200,0,400);
    PtRecott3   = new TH1F("PtRecott3","Pt of Reco top mass w/ 6j", 200,0,400);
    PtRecott4   = new TH1F("PtRecott4","Pt of Reco top mass w/ >= 7j", 200,0,400);

    WMCtt      = new TH2F("WMCtt",   "W vs Mt1+Mt2 of MCMatched tt", 500,0,1000, 400,0,800);
    WRecott0   = new TH2F("WRecott0","W vs Mt1+Mt2 of Reco top mass w/ > 3j", 500,0,1000, 400,0,800);
    WRecott1   = new TH2F("WRecott1","W vs Mt1+Mt2 of Reco top mass w/ 4j", 500,0,1000, 400,0,800);
    WRecott2   = new TH2F("WRecott2","W vs Mt1+Mt2 of Reco top mass w/ 5j", 500,0,1000, 400,0,800);
    WRecott3   = new TH2F("WRecott3","W vs Mt1+Mt2 of Reco top mass w/ 6j", 500,0,1000, 400,0,800);
    WRecott4   = new TH2F("WRecott4","W vs Mt1+Mt2of Reco top mass w/ >= 7j", 500,0,1000, 400,0,800);

    MClepW   = new TH1F("MClepW",  " mt of leptonic W w/ highest Pt"   ,320,0,160.);
    MChadW   = new TH1F("MChadW",  " MC hadronic W  w/ closest w mass"  ,320,0,160.);
    RecolepW = new TH1F("RecolepW"," mt of leptonic W w/ highest Pt"    ,320,0,160.);
    RecohadW = new TH1F("RecohadW"," Reco hadronic W w/ closest w mass" ,320,0,160.);
    allMClepW   = new TH1F("allMClepW",  " all mt of leptonic W ",320,0,160.);
    allMChadW   = new TH1F("allMChadW",  " all MC hadronic W "   ,320,0,160.);
    allRecolepW = new TH1F("allRecolepW"," all mt of leptonic W ",320,0,160.);
    allRecohadW = new TH1F("allRecohadW"," all Reco hadronic W " ,320,0,160.);

    selMCW   = new TH2F("selMCW",  " MC  : mt of chosen lepW , hadW ",320,0,160.,320,0,160.);
    selRecoW = new TH2F("selRecoW"," Reco: mt of chosen lepW , hadW ",320,0,160.,320,0,160.);
    lRelPtMC  = new TH2F("lRelPtMC", " MC  : RelPt of W w.r.t LepTop ",32,-0.06,1.22,320,0,160.);
    hRelPtMC  = new TH2F("hRelPtMC", " MC  : RelPt of W w.r.t HadTop ",32,-0.06,1.22,320,0,160.);
    lRelPtRC  = new TH2F("lRelPtRC", " Reco: RelPt of W w.r.t LepTop ",32,-0.06,1.22,320,0,160.);
    hRelPtRC  = new TH2F("hRelPtRC", " Reco: RelPt of W w.r.t HadTop ",32,-0.06,1.22,320,0,160.);

    hEvtEff   = new TH1F("hEvtEff","(pass,fail)=> had(0,0.5) Mu(1,1.5) dilep(2,2.5) Ele(3,3.5), Tau(4,4.5), other(5,5.5)", 12,-0.25,5.75 );
    hObjEff   = new TH1F("hObjEff"," object selection Eff", 10, -0.5,9.5 );
    hbJetEff   = new TH2F("hbJetEff"," # of correctly matched bjet ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hWJetEff   = new TH2F("hWJetEff"," # of correctly matched Wjet ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hWLepEff   = new TH2F("hWLepEff"," # of correctly matched lep. ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hHLTBits   = new TH2F("hHLTBits","HLT Trigger bits(names) for hadronic tt", 166, -0.5, 165.5, 6,-0.5,5.5 );
    hHLTSelect = new TH1F("hHLTSelect","HLT Trigger Selected result "  , 15, -7.5, 7.5 );

    lepRCTMass = new TH1F("lepRCTMass",  " mass of Reco leptonic T w/ highest Pt " ,200,0,400.);
    lepMCTMass = new TH1F("lepMCTMass",  " mass of MC   leptonic T w/ highest Pt " ,200,0,400.);
    hadRCTMass = new TH1F("hadRCTMass",  " mass of Reco hadronic T w/ highest Pt " ,200,0,400.);
    hadMCTMass = new TH1F("hadMCTMass",  " mass of MCi  hadronic T w/ highest Pt " ,200,0,400.);

    alepRCTMass = new TH1F("alepRCTMass",  " mass of all Reco leptonic T " ,200,0,400.);
    alepMCTMass = new TH1F("alepMCTMass",  " mass of all MC   leptonic T " ,200,0,400.);
    ahadRCTMass = new TH1F("ahadRCTMass",  " mass of all Reco hadronic T " ,200,0,400.);
    ahadMCTMass = new TH1F("ahadMCTMass",  " mass of all MC   hadronic T " ,200,0,400.);

 } 

 HTOP9(TString name_, TFile* file) {
    name=name_;

    hMCtt      = (TH2F *) file->Get("hMCtt");
    hMCtt1     = (TH2F *) file->Get("hMCtt1");
    hRecott0   = (TH2F *) file->Get("hRecott0");
    hRecott1   = (TH2F *) file->Get("hRecott1");
    hRecott2   = (TH2F *) file->Get("hRecott2");
    hRecott3   = (TH2F *) file->Get("hRecott3");
    hRecott4   = (TH2F *) file->Get("hRecott4");
    PtMCtt     = (TH1F *) file->Get("PtMCtt");
    PtRecott0  = (TH1F *) file->Get("PtRecott0");
    PtRecott1  = (TH1F *) file->Get("PtRecott1");
    PtRecott2  = (TH1F *) file->Get("PtRecott2");
    PtRecott3  = (TH1F *) file->Get("PtRecott3");
    PtRecott4  = (TH1F *) file->Get("PtRecott4");
    WMCtt      = (TH2F *) file->Get("WMCtt");
    WRecott0   = (TH2F *) file->Get("WRecott0");
    WRecott1   = (TH2F *) file->Get("WRecott1");
    WRecott2   = (TH2F *) file->Get("WRecott2");
    WRecott3   = (TH2F *) file->Get("WRecott3");
    WRecott4   = (TH2F *) file->Get("WRecott4");

    MClepW    = (TH1F *) file->Get("MClepW");
    MChadW    = (TH1F *) file->Get("MChadW");
    RecolepW  = (TH1F *) file->Get("RecolepW");
    RecohadW  = (TH1F *) file->Get("RecohadW");
    allMClepW    = (TH1F *) file->Get("allMClepW");
    allMChadW    = (TH1F *) file->Get("allMChadW");
    allRecolepW  = (TH1F *) file->Get("allRecolepW");
    allRecohadW  = (TH1F *) file->Get("allRecohadW");
    selMCW    = (TH2F *) file->Get("selMCW");
    selRecoW  = (TH2F *) file->Get("selRecoW");
    lRelPtMC   = (TH2F *) file->Get("lRelPtMC");
    hRelPtMC   = (TH2F *) file->Get("hRelPtMC");
    lRelPtRC   = (TH2F *) file->Get("lRelPtRC");
    hRelPtRC   = (TH2F *) file->Get("hRelPtRC");

    hEvtEff    = (TH1F *) file->Get("hEvtEff");
    hObjEff    = (TH1F *) file->Get("hObjEff");
    hbJetEff   = (TH2F *) file->Get("hbJetEff");
    hWJetEff   = (TH2F *) file->Get("hWJetEff");
    hWLepEff   = (TH2F *) file->Get("hWLepEff");
    hHLTBits   = (TH2F *) file->Get("hHLTBits");
    hHLTSelect = (TH1F *) file->Get("hHLTSelect");

    lepRCTMass  = (TH1F *) file->Get("lepRCTMass");
    lepMCTMass  = (TH1F *) file->Get("lepMCTMass");
    hadRCTMass  = (TH1F *) file->Get("hadRCTMass");
    hadMCTMass  = (TH1F *) file->Get("hadMCTMass");

    alepRCTMass  = (TH1F *) file->Get("alepRCTMass");
    alepMCTMass  = (TH1F *) file->Get("alepMCTMass");
    ahadRCTMass  = (TH1F *) file->Get("ahadRCTMass");
    ahadMCTMass  = (TH1F *) file->Get("ahadMCTMass");
 }

 /// Destructor
 virtual ~HTOP9() {

    delete hMCtt;
    delete hMCtt1;
    delete hRecott0;
    delete hRecott1;
    delete hRecott2;
    delete hRecott3;
    delete hRecott4;
    delete PtMCtt;
    delete PtRecott0;
    delete PtRecott1;
    delete PtRecott2;
    delete PtRecott3;
    delete PtRecott4;
    delete WMCtt;
    delete WRecott0;
    delete WRecott1;
    delete WRecott2;
    delete WRecott3;
    delete WRecott4;

    delete MClepW;
    delete MChadW;
    delete RecolepW;
    delete RecohadW;
    delete allMClepW;
    delete allMChadW;
    delete allRecolepW;
    delete allRecohadW;
    delete selMCW;
    delete selRecoW;
    delete lRelPtMC;
    delete hRelPtMC;
    delete lRelPtRC;
    delete hRelPtRC;

    delete hEvtEff;
    delete hObjEff;
    delete hbJetEff;
    delete hWJetEff;
    delete hWLepEff;
    delete hHLTBits;
    delete hHLTSelect;
   
    delete lepRCTMass;
    delete lepMCTMass;
    delete hadRCTMass;
    delete hadMCTMass;

    delete alepRCTMass;
    delete alepMCTMass;
    delete ahadRCTMass;
    delete ahadMCTMass;
 }

 void Fill9(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hMCtt->Fill( lep_mt, had_mt );
    PtMCtt->Fill( PtTt );
    WMCtt->Fill( Wtt, Wt2 );
 }
 void Fill9g(double lep_mt, double had_mt){
    hMCtt1->Fill( lep_mt, had_mt );
 }
 void Fill9a0(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott0->Fill( lep_mt, had_mt );
    PtRecott0->Fill( PtTt );
    WRecott0->Fill( Wtt, Wt2 );
 }
 void Fill9a1(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott1->Fill( lep_mt, had_mt );
    PtRecott1->Fill( PtTt );
    WRecott1->Fill( Wtt, Wt2 );
 }
 void Fill9a2(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott2->Fill( lep_mt, had_mt );
    PtRecott2->Fill( PtTt );
    WRecott2->Fill( Wtt, Wt2 );
 }
 void Fill9a3(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott3->Fill( lep_mt, had_mt );
    PtRecott3->Fill( PtTt );
    WRecott3->Fill( Wtt, Wt2 );
 }
 void Fill9a4(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott4->Fill( lep_mt, had_mt );
    PtRecott4->Fill( PtTt );
    WRecott4->Fill( Wtt, Wt2 );
 }
 void Fill9b(double mw) {
    MClepW->Fill(mw) ; 
 }
 void Fill9c(double mtw) {
    RecolepW->Fill(mtw) ; 
 }
 void Fill9d(double mw) {
    MChadW->Fill(mw) ; 
 }
 void Fill9e(double mw) {
    RecohadW->Fill(mw) ; 
 }
 void Fill9f(float count ) {
    hEvtEff->Fill( count );
 }
 void Fill9i( int count ) {
    hObjEff->Fill( count );
 }
 void Fill9j( double dm, int nb, int nw, int nl ) {
    hbJetEff->Fill( nb , dm );
    hWJetEff->Fill( nw , dm );
    hWLepEff->Fill( nl , dm );
 }

 void Fill9k( int hlt, int topo ) {
    hHLTBits->Fill( hlt, topo );
 }
 void Fill9l( int hlt ) {
    hHLTSelect->Fill( hlt );
 }
 void Fill9m( double wmt, double wmass, double lRelPt, double hRelPt, double ltbeta, double htbeta ) {
    selMCW->Fill(wmt, wmass);
    lRelPtMC->Fill(ltbeta, lRelPt);
    hRelPtMC->Fill(htbeta, hRelPt);
 }
 void Fill9n( double wmt, double wmass, double lRelPt, double hRelPt, double ltbeta, double htbeta ) {
    selRecoW->Fill(wmt, wmass);
    lRelPtRC->Fill(ltbeta, lRelPt);
    hRelPtRC->Fill(htbeta, hRelPt);
 }

 void Fill9p1(double mtw) {
    allMClepW->Fill(mtw) ; 
 }
 void Fill9p2(double mtw) {
    allRecolepW->Fill(mtw) ; 
 }
 void Fill9p3(double mw) {
    allMChadW->Fill(mw) ; 
 }
 void Fill9p4(double mw) {
    allRecohadW->Fill(mw) ; 
 }

 void Fill9q1( double mass ){
     lepRCTMass->Fill(mass) ; 
 }
 void Fill9q2( double mass ){
     lepMCTMass->Fill(mass) ; 
 }
 void Fill9q3( double mass ){
     hadRCTMass->Fill(mass) ; 
 }
 void Fill9q4( double mass ){
     hadMCTMass->Fill(mass) ; 
 }

 void Fill9r1( double mass ){
     alepRCTMass->Fill(mass) ; 
 }
 void Fill9r2( double mass ){
     alepMCTMass->Fill(mass) ; 
 }
 void Fill9r3( double mass ){
     ahadRCTMass->Fill(mass) ; 
 }
 void Fill9r4( double mass ){
     ahadMCTMass->Fill(mass) ; 
 }

 void Write() {

    hMCtt->Write();
    hMCtt1->Write();
    hRecott0->Write();
    hRecott1->Write();
    hRecott2->Write();
    hRecott3->Write();
    hRecott4->Write();
    PtMCtt->Write();
    PtRecott0->Write();
    PtRecott1->Write();
    PtRecott2->Write();
    PtRecott3->Write();
    PtRecott4->Write();
    WMCtt->Write();
    WRecott0->Write();
    WRecott1->Write();
    WRecott2->Write();
    WRecott3->Write();
    WRecott4->Write();

    MClepW->Write();
    MChadW->Write();
    RecolepW->Write();
    RecohadW->Write();
    allMClepW->Write();
    allMChadW->Write();
    allRecolepW->Write();
    allRecohadW->Write();
    selMCW->Write();
    selRecoW->Write();
    lRelPtMC->Write();
    hRelPtMC->Write();
    lRelPtRC->Write();
    hRelPtRC->Write();

    hEvtEff->Write();
    hObjEff->Write();
    hbJetEff->Write();
    hWJetEff->Write();
    hWLepEff->Write();
    hHLTBits->Write();
    hHLTSelect->Write();

    lepRCTMass->Write();
    lepMCTMass->Write();
    hadRCTMass->Write();
    hadMCTMass->Write();

    alepRCTMass->Write();
    alepMCTMass->Write();
    ahadRCTMass->Write();
    ahadMCTMass->Write();
 }

  TH2F *hMCtt;
  TH2F *hMCtt1;
  TH2F *hRecott0;
  TH2F *hRecott1;
  TH2F *hRecott2;
  TH2F *hRecott3;
  TH2F *hRecott4;

  TH1F *PtMCtt;
  TH1F *PtRecott0;
  TH1F *PtRecott1;
  TH1F *PtRecott2;
  TH1F *PtRecott3;
  TH1F *PtRecott4;
  TH2F *WMCtt;
  TH2F *WRecott0;
  TH2F *WRecott1;
  TH2F *WRecott2;
  TH2F *WRecott3;
  TH2F *WRecott4;

  TH1F *MClepW;
  TH1F *MChadW;
  TH1F *RecolepW;
  TH1F *RecohadW;
  TH1F *allMClepW;
  TH1F *allMChadW;
  TH1F *allRecolepW;
  TH1F *allRecohadW;
  TH2F *selMCW;
  TH2F *selRecoW;
  TH2F *lRelPtMC;
  TH2F *hRelPtMC;
  TH2F *lRelPtRC;
  TH2F *hRelPtRC;

  TH1F *hEvtEff; 
  TH1F *hObjEff; 
  TH2F *hbJetEff;
  TH2F *hWJetEff;
  TH2F *hWLepEff;
  TH2F *hHLTBits;
  TH1F *hHLTSelect;

  TH1F *lepRCTMass;
  TH1F *lepMCTMass;
  TH1F *hadRCTMass;
  TH1F *hadMCTMass;

  TH1F *alepRCTMass;
  TH1F *alepMCTMass;
  TH1F *ahadRCTMass;
  TH1F *ahadMCTMass;

 TString name;

};
 
#endif
