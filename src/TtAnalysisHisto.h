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
 
 HTOP1() {

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

    Et_Pt   = new TH2F("Et_Pt", " PAT Et vs Pt (all jets)", 500, 0, 500, 500, 0., 500);
    
    hWmass = new TH1F("hWmass"," W mass from gen-reco matching jets",240,30.,150.);
    hWp    = new TH1F("hWp", " momentum of W ",500,0.,500.);

    EovH1    = new TH1F("EovH1"," em Energy / hadronic Energy ", 500, -2., 23.);
    EovH2    = new TH1F("EovH2"," em Energy / hadronic Energy ", 500, -2., 23.);
    Cont_Trk1= new TH2F("Cont_Trk1","N Constituent vs N Charged Track ", 50, -0.5, 49.5, 50, -0.5, 49.5);
    Cont_Trk2= new TH2F("Cont_Trk2","N Constituent vs N Charged Track ", 50, -0.5, 49.5, 50, -0.5, 49.5);

    hNJets = new TH1F("hNJets", " # of Jets",65,-0.5,64.5);

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

    wEta_njets = new TH2F("wEta_njets", "w eta vs njets ", 91, -4.55, 4.55, 25, -0.5, 24.5);

    wEta_dPhi0 = new TH2F("wEta_dPhi0", " w eta vs dPhi(w,jet1)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    wEta_dPhi1 = new TH2F("wEta_dPhi1", " w eta vs dPhi(w,jet2)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    wEta_dPhi2 = new TH2F("wEta_dPhi2", " w eta vs dPhi(w,jet3)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    njets_dPhi0= new TH2F("njets_dPhi0", " njets vs dPhi(w,jet1)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    njets_dPhi1= new TH2F("njets_dPhi1", " njets vs dPhi(w,jet2)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    njets_dPhi2= new TH2F("njets_dPhi2", " njets vs dPhi(w,jet3)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    
    njets_dR0= new TH2F("njets_dR0", " njets vs dR(w,jet1)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    njets_dR1= new TH2F("njets_dR1", " njets vs dR(w,jet2)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    njets_dR2= new TH2F("njets_dR2", " njets vs dR(w,jet3)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);

    njets_dR = new TH2F("njets_dR", " njets vs dR(w,jet1+2+3)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    njets_m3 = new TH2F("njets_m3", " njets vs m3 of Jet1,2,3", 20, -0.5, 19.5, 500,0.,500. );
    m3_dR    = new TH2F("m3_dR", " m3 vs dR(w,jet1+2+3)", 500, 0., 500., 100, -0.05, 9.95);

 } 

 HTOP1( TString theFolder, TFile* file ) {

    hEt    = (TH1F *) file->Get(theFolder+"hEt");
    hEovH  = (TH1F *) file->Get(theFolder+"hEovH");
    hNCont = (TH1F *) file->Get(theFolder+"hNCont");
    hNTk   = (TH1F *) file->Get(theFolder+"hNTK");
    hR60   = (TH1F *) file->Get(theFolder+"hR60");
    hArea_pt = (TH2F *) file->Get(theFolder+"hArea_pt");
    hEovH_pt = (TH2F *) file->Get(theFolder+"hEovH_pt");
    hEovH_r  = (TH2F *) file->Get(theFolder+"hEovH_r");
    hEovH_A  = (TH2F *) file->Get(theFolder+"hEovH_A");
    hemF     = (TH2F *) file->Get(theFolder+"hemF");

    Et_Pt  = (TH2F *) file->Get(theFolder+"Et_Pt");

    hWmass = (TH1F *) file->Get(theFolder+"hWmass");
    hWp    = (TH1F *) file->Get(theFolder+"hWp");

    EovH1   = (TH1F *) file->Get(theFolder+"EovH1");
    EovH2   = (TH1F *) file->Get(theFolder+"EovH2");
    Cont_Trk1 = (TH2F *) file->Get(theFolder+"Cont_Trk1");
    Cont_Trk2 = (TH2F *) file->Get(theFolder+"Cont_Trk2");
    
    hNJets = (TH1F *) file->Get(theFolder+"hNJets");

    hEovH_A1  = (TH2F *) file->Get(theFolder+"hEovH_A1");
    hEovH_r1  = (TH2F *) file->Get(theFolder+"hEovH_r1");
    hEovH_N1  = (TH2F *) file->Get(theFolder+"hEovH_N1");
    hEovH_C1  = (TH2F *) file->Get(theFolder+"hEovH_C1");
    hR60_1    = (TH1F *) file->Get(theFolder+"hR60_1");

    gEta_Pt = (TH2F *) file->Get(theFolder+"gEta_Pt");
    gEovH   = (TH1F *) file->Get(theFolder+"gEovH");
    gNJets  = (TH1F *) file->Get(theFolder+"gNJets");
    gNjet20 = (TH1F *) file->Get(theFolder+"gNjet20");
    gNjet40 = (TH1F *) file->Get(theFolder+"gNjet40");
    gNjeth35= (TH1F *) file->Get(theFolder+"gNjeth35");

    wEta_njets  = (TH2F *) file->Get(theFolder+"wEta_njets");

    wEta_dPhi0 = (TH2F *) file->Get(theFolder+"wEta_dPhi0");
    wEta_dPhi1 = (TH2F *) file->Get(theFolder+"wEta_dPhi1");
    wEta_dPhi2 = (TH2F *) file->Get(theFolder+"wEta_dPhi2");
    njets_dPhi0= (TH2F *) file->Get(theFolder+"njets_dPhi0");
    njets_dPhi1= (TH2F *) file->Get(theFolder+"njets_dPhi1");
    njets_dPhi2= (TH2F *) file->Get(theFolder+"njets_dPhi2");
 
    njets_dR0= (TH2F *) file->Get(theFolder+"njets_dR0");
    njets_dR1= (TH2F *) file->Get(theFolder+"njets_dR1");
    njets_dR2= (TH2F *) file->Get(theFolder+"njets_dR2");

    njets_dR = (TH2F *) file->Get(theFolder+"njets_dR");
    m3_dR    = (TH2F *) file->Get(theFolder+"m3_dR");
    njets_m3 = (TH2F *) file->Get(theFolder+"njets_m3");

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

    delete wEta_njets;

    delete wEta_dPhi0;
    delete wEta_dPhi1;
    delete wEta_dPhi2;
    delete njets_dPhi0;
    delete njets_dPhi1;
    delete njets_dPhi2;

    delete njets_dR0;
    delete njets_dR1;
    delete njets_dR2;

    delete njets_dR;
    delete m3_dR;
    delete njets_m3;

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
    EovH2->Fill( EovH );
    Cont_Trk2->Fill(nCon, nTrk);
 }
 void Fill1f( int njets ){
    hNJets->Fill(njets);  
 }
 void Fill1h( double EovH, double area, int nCon, int nTrk, double r, double r60 ) {
    hEovH_A1->Fill(area,EovH);
    hEovH_N1->Fill(nCon,EovH);
    hEovH_C1->Fill(nTrk,EovH);
    hEovH_r1->Fill(r,EovH);
    hR60_1->Fill(r60);
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
 void Fill1p( double njets, double eta ) {
    wEta_njets->Fill( eta, njets );
 }


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

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

    wEta_njets->Write();

    wEta_dPhi0->Write();
    wEta_dPhi1->Write();
    wEta_dPhi2->Write();
    njets_dPhi0->Write();
    njets_dPhi1->Write();
    njets_dPhi2->Write();
    
    njets_dR0->Write();
    njets_dR1->Write();
    njets_dR2->Write();

    njets_dR->Write();
    m3_dR->Write();
    njets_m3->Write();

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

  TH2F *wEta_njets;

  TH2F *wEta_dPhi0;
  TH2F *wEta_dPhi1;
  TH2F *wEta_dPhi2;
  TH2F *njets_dPhi0;
  TH2F *njets_dPhi1;
  TH2F *njets_dPhi2;

  TH2F *njets_dR0;
  TH2F *njets_dR1;
  TH2F *njets_dR2;
  TH2F *njets_dR;
  TH2F *m3_dR;
  TH2F *njets_m3;

};

class HTOP2 {
public:
 
 HTOP2() {

   // reconstructed objects masses
    hemfMET   = new TH2F("hemfMET"," emfCalo(Y), emfPAT(X)", 150, -0.2, 1.3, 150, -0.2, 1.3);
    hSumEt_MET= new TH2F("hSumEt_MET","sumEt vs MET ", 1500, 0,1500, 500, 0,500);

    hResPAT = new TH1F("ResPAT","PAT MET Resolution", 100, -5.,5.);
    hResCal = new TH1F("ResCal","Cal MET Resolution", 100, -5.,5.);
    hResCor = new TH1F("ResCor","Cor MET Resolution", 100, -5.,5.);
    hResTC  = new TH1F("ResTC", " TC MET Resolution", 100, -5.,5.);
    hPhiPAT = new TH1F("PhiPAT","PAT Phi Resolution", 100, -5.,5.);
    hPhiCal = new TH1F("PhiCal","Cal Phi Resolution", 100, -5.,5.);
    hPhiCor = new TH1F("PhiCor","Cor Phi Resolution", 100, -5.,5.);
    hPhiTC  = new TH1F("PhiTC" ," TC Phi Resolution", 100, -5.,5.);

    tc_patMET = new TH2F("tc_patMET"," tcMET vs patMet ",  400, 0,400, 400, 0,400);
    tc_patPhi = new TH2F("tc_patPhi"," tcMET vs patMet ",  165, -3.15,3.15, 165, -3.15,3.15);

    MET_dPhi0 = new TH2F("MET_dPhi0","MET vs dPhi(Mu,MET) > 3j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi1 = new TH2F("MET_dPhi1","MET vs dPhi(Mu,MET) 4j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi2 = new TH2F("MET_dPhi2","MET vs dPhi(Mu,MET) 5j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi3 = new TH2F("MET_dPhi3","MET vs dPhi(Mu,MET) 6j", 500, 0,500, 32, -0.05,3.15);
    MET_dPhi4 = new TH2F("MET_dPhi4","MET vs dPhi(Mu,MET) > 6j", 500, 0,500, 32, -0.05,3.15);

    MET_muPt = new TH2F("MET_muPt"," MET vs IsoMuon Pt ",  400, 0,400, 400, 0,400);
    dPhi_METJ1  = new TH2F("dPhi_METJ1"," MET vs 1st Jet ", 400, 0,400,32, -0.05,3.15);
    dPhi_METJ12 = new TH2F("dPhi_METJ12"," MET vs 1st + 2nd Jet ", 400, 0,400,32, -0.05,3.15);
 } 

 HTOP2( TString theFolder, TFile* file ) {

    hemfMET     = (TH2F *) file->Get(theFolder+"hemfMET");
    hSumEt_MET  = (TH2F *) file->Get(theFolder+"hSumEt_MET");
    hResPAT     = (TH1F *) file->Get(theFolder+"hResPAT");
    hResCal     = (TH1F *) file->Get(theFolder+"hResCal");
    hResCor     = (TH1F *) file->Get(theFolder+"hResCor");
    hResTC      = (TH1F *) file->Get(theFolder+"hResTC");
    hPhiPAT     = (TH1F *) file->Get(theFolder+"hPhiPAT");
    hPhiCal     = (TH1F *) file->Get(theFolder+"hPhiCal");
    hPhiCor     = (TH1F *) file->Get(theFolder+"hPhiCor");
    hPhiTC      = (TH1F *) file->Get(theFolder+"hPhiTC");

    tc_patMET   = (TH2F *) file->Get(theFolder+"tc_patMET");
    tc_patPhi   = (TH2F *) file->Get(theFolder+"tc_patPhi");

    MET_dPhi0   = (TH2F *) file->Get(theFolder+"MET_dPhi0");
    MET_dPhi1   = (TH2F *) file->Get(theFolder+"MET_dPhi1");
    MET_dPhi2   = (TH2F *) file->Get(theFolder+"MET_dPhi2");
    MET_dPhi3   = (TH2F *) file->Get(theFolder+"MET_dPhi3");
    MET_dPhi4   = (TH2F *) file->Get(theFolder+"MET_dPhi4");

    MET_muPt    = (TH2F *) file->Get(theFolder+"MET_muPt");
    dPhi_METJ1  = (TH2F *) file->Get(theFolder+"dPhi_METJ1");
    dPhi_METJ12 = (TH2F *) file->Get(theFolder+"dPhi_METJ12");
 }

 /// Destructor
 virtual ~HTOP2() {

    delete hemfMET;
    delete hSumEt_MET;
    delete hResPAT;
    delete hResCal;
    delete hResCor;
    delete hResTC;
    delete hPhiPAT;
    delete hPhiCal;
    delete hPhiCor;
    delete hPhiTC;

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

 void Fill2b(double respat, double rescal, double rescor, double restc, double phipat, double phical, double phicor, double phitc ){
    hResPAT->Fill(respat);
    hResCal->Fill(rescal);
    hResCor->Fill(rescor);
    hResTC->Fill(restc);
    hPhiPAT->Fill(phipat);
    hPhiCal->Fill(phical);
    hPhiCor->Fill(phicor);
    hPhiTC->Fill(phitc);
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
 void Fill2f( double tcmet, double patmet, double tcphi, double patphi ){
    tc_patMET->Fill(tcmet, patmet);
    tc_patPhi->Fill(tcphi, patphi);
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    hemfMET->Write();
    hSumEt_MET->Write();

    hResPAT->Write();
    hResCal->Write();
    hResCor->Write();
    hResTC->Write();
    hPhiPAT->Write();
    hPhiCal->Write();
    hPhiCor->Write();
    hPhiTC->Write();

    tc_patMET->Write();
    tc_patPhi->Write();

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
  TH1F *hResTC;
  TH1F *hPhiPAT;
  TH1F *hPhiCal;
  TH1F *hPhiCor;
  TH1F *hPhiTC;

  TH2F *tc_patMET;
  TH2F *tc_patPhi;

  TH2F *MET_dPhi0;
  TH2F *MET_dPhi1;
  TH2F *MET_dPhi2;
  TH2F *MET_dPhi3;
  TH2F *MET_dPhi4;

  TH2F *MET_muPt;
  TH2F *dPhi_METJ1;
  TH2F *dPhi_METJ12;

};

class HTOP3 {
public:
 
 HTOP3() {

    // parton level
    hEta     = new TH2F("hEta"," MC(X) vs reco(Y) , eta distribution for muon of semiTt ",71,-3.55,3.55,71,-3.55,3.55);
    hPhi     = new TH2F("hPhi"," MC(X) vs reco(Y) , phi distribution for muon of semiTt ",63,-3.15,3.15,71,-3.15,3.15);
    Pt_Resol = new TH2F("hPt_Resol", " Pt MC(X) vs Resol_Pt(Y) distribution for muon of semiTt",60,0,300, 200,-0.995,1.005 );

    hIMu_caloE = new TH2F("hIMu_caloE","IsoMu caloE vs muon p ",800,0.,200., 800,0.,200. );
    hJMu_caloE = new TH2F("hJMu_caloE","JetMu caloE vs muon p ",800,0.,200., 800,0.,200. );
    hIMu_caloE_jE = new TH2F("hIMu_caloE_jE","IsoMu caloE vs jet E ",800,0.,200., 800,0.,200. );
    hJMu_caloE_jE = new TH2F("hJMu_caloE_jE","JetMu caloE vs jet E ",800,0.,200., 800,0.,200. );


 } 

 HTOP3( TString theFolder, TFile* file ) {

    hEta     = (TH2F *) file->Get(theFolder+"hEta");
    hPhi     = (TH2F *) file->Get(theFolder+"hPhi");
    Pt_Resol = (TH2F *) file->Get(theFolder+"Pt_Resol");

    hIMu_caloE = (TH2F *) file->Get(theFolder+"hIMu_caloE");
    hJMu_caloE = (TH2F *) file->Get(theFolder+"hJMu_caloE");
    hIMu_caloE_jE = (TH2F *) file->Get(theFolder+"hIMu_caloE_jE");
    hJMu_caloE_jE = (TH2F *) file->Get(theFolder+"hJMu_caloE_jE");

 }

 /// Destructor
 virtual ~HTOP3() {

    delete hEta;
    delete hPhi;
    delete Pt_Resol;
 
    delete hIMu_caloE;
    delete hJMu_caloE;
    delete hIMu_caloE_jE;
    delete hJMu_caloE_jE;

 }


 void Fill3b(double eta_mc, double phi_mc, double eta_rc, double phi_rc, double pt, double ptResol) {
    hEta->Fill(eta_mc, eta_rc);
    hPhi->Fill(phi_mc, phi_rc);
    Pt_Resol->Fill(pt, ptResol);
 }
 void Fill3j(double isoMuE, double mu_p, double jet_E ) {
    hIMu_caloE->Fill(isoMuE, mu_p);
    hIMu_caloE_jE->Fill(isoMuE, jet_E);
 }
 void Fill3k(double jetMuE, double mu_p, double jet_E ) {
    hJMu_caloE->Fill(jetMuE, mu_p);
    hJMu_caloE_jE->Fill(jetMuE, jet_E);
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    hEta->Write();
    hPhi->Write();
    Pt_Resol->Write();
 
    hIMu_caloE->Write();
    hJMu_caloE->Write();
    hIMu_caloE_jE->Write();
    hJMu_caloE_jE->Write();

 }

  TH2F *hEta;
  TH2F *hPhi;
  TH2F *Pt_Resol;

  TH2F *hIMu_caloE;
  TH2F *hJMu_caloE;
  TH2F *hIMu_caloE_jE;
  TH2F *hJMu_caloE_jE;

};

class HTOP4 {
public:
 
 HTOP4() {

   // reconstructed objects masses
    allPt   = new TH1F("allPt", " Pt distribution  for all e",500,0,500);
    allEta  = new TH1F("allEta"," eta distribution for all e",59,-2.95,2.95);

    // parton level
    hEta     = new TH2F("hEta"," MC(X) vs reco(Y) , eta distribution for electron of semiTt ",71,-3.55,3.55,71,-3.55,3.55);
    hPhi     = new TH2F("hPhi"," MC(X) vs reco(Y) , phi distribution for electron of semiTt ",63,-3.15,3.15,71,-3.15,3.15);
    Pt_Resol = new TH2F("hPt_Resol", " Pt MC(X) vs Resol_Pt(Y) distribution for electron of semiTt",60,0,300, 200,-0.995,1.005 );

    isoEleCut = new TH2F("isoEleCut","N of isoEle vs N of SelectedJets  after isoMu =1 cut ", 21, -0.5, 20.5, 21, -0.5, 20.5 );
 } 

 HTOP4( TString theFolder, TFile* file) {

    allPt    = (TH1F *) file->Get(theFolder+"allPt");
    allEta   = (TH1F *) file->Get(theFolder+"allEta");

    hEta     = (TH2F *) file->Get(theFolder+"hEta");
    hPhi     = (TH2F *) file->Get(theFolder+"hPhi");
    Pt_Resol = (TH2F *) file->Get(theFolder+"Pt_Resol");

    isoEleCut   =(TH2F *) file->Get(theFolder+"isoEleCut");

 }

 /// Destructor
 virtual ~HTOP4() {

    delete allPt;
    delete allEta;

    delete hEta;
    delete hPhi;
    delete Pt_Resol;
 
    delete isoEleCut;
 }

 void Fill4a(double Pt,double Eta ){
    allPt->Fill(Pt);
    allEta->Fill(Eta);
 }
 
 void Fill4g(double eta_mc, double phi_mc, double eta_rc, double phi_rc, double pt, double ptResol) {
    hEta->Fill(eta_mc, eta_rc);
    hPhi->Fill(phi_mc, phi_rc);
    Pt_Resol->Fill(pt, ptResol);
 }

 void Fill4i( int NIsoEle, int NSelJet ){
    isoEleCut->Fill( NIsoEle, NSelJet );
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    allPt->Write();
    allEta->Write();
  
    hEta->Write();
    hPhi->Write();
    Pt_Resol->Write();

    isoEleCut->Write();
 }

  TH1F *allPt;
  TH1F *allEta;

  TH2F *hEta;
  TH2F *hPhi;
  TH2F *Pt_Resol;

  TH2F *isoEleCut;

};

class HTOP5 {
public:
 
 HTOP5() {

    hPt   = new TH1F("hPt", " Pt distribution",500,0,500);
    hEta  = new TH1F("hEta"," eta distribution",61,-3.05,3.05);

 } 

 HTOP5( TString theFolder, TFile* file) {

    hPt    = (TH1F *) file->Get(theFolder+"hPt");
    hEta   = (TH1F *) file->Get(theFolder+"hEta");

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


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    hPt->Write();
    hEta->Write();
 }

  TH1F *hPt;
  TH1F *hEta;

};

class HTOP6 {
public:
 
 HTOP6() {

    dRb_dRw_lep = new TH2F("dRb_dRw_lep","leptonic T dR(mcB, recoB), dR(mcW, recoW) ", 200,-0.025,9.975, 200, -0.025, 9.975 );
    dRb_dRw_had = new TH2F("dRb_dRw_had","hadronic T dR(mcB, recoB), dR(mcW, recoW) ", 200,-0.025,9.975, 200, -0.025, 9.975 );
    dRtt        = new TH2F("dRtt"," dR(rec, mc)_lep, dR(rec, mc)_had", 400, -0.22, 9.78, 400, -0.22, 9.78 );

    dRbw_had = new TH2F("dRbw_had","hadronic T dR(mcb, mcW), dR( rcb, rcW )", 200, -0.025, 9.975, 200, -0.025, 9.975 );
    dRbw_lep = new TH2F("dRbw_lep","leptonic T dR(mcb, mcW), dR( rcb, rcW )", 200, -0.025, 9.975, 200, -0.025, 9.975 );

    dR_neu     = new TH1F("dR_neu"," dR(reco - gen) for neutrino ",200,0.,10.);
    PtRes_lepW = new TH1F("PtRes_lepW"," Pt Res for matched leptonic W ",200,-1.005,0.995);
    PzRes_lepW = new TH1F("PzRes_lepW"," Pz Res for matched leptonic W ",200,-1.005,0.995);
 
    mcTtMass   = new TH2F("mcTtMass"  ," mc Tt mass matched w/ reco cuts all", 480,0.,480., 480,0.,480.);
    mcTtMass0  = new TH2F("mcTtMass0" ," mc Tt mass matched w/ reco cuts  4j", 480,0.,480., 480,0.,480.);
    mcTtMass1  = new TH2F("mcTtMass1" ," mc Tt mass matched w/ reco cuts  5j", 480,0.,480., 480,0.,480.);
    mcTtMass2  = new TH2F("mcTtMass2" ," mc Tt mass matched w/ reco cuts  6j", 480,0.,480., 480,0.,480.);
    mcTtMass3  = new TH2F("mcTtMass3" ," mc Tt mass matched w/ reco cuts >7j", 480,0.,480., 480,0.,480.);

 } 

 HTOP6( TString theFolder, TFile* file ) {

    dRb_dRw_lep= (TH2F *) file->Get(theFolder+"dRb_dRw_lep");
    dRb_dRw_had= (TH2F *) file->Get(theFolder+"dRb_dRw_had");
    dRbw_had   = (TH2F *) file->Get(theFolder+"dRbw_had");
    dRbw_lep   = (TH2F *) file->Get(theFolder+"dRbw_lep");
    dRtt       = (TH2F *) file->Get(theFolder+"dRtt");

    dR_neu    = (TH1F *) file->Get(theFolder+"dR_neu");
    PtRes_lepW = (TH1F *) file->Get(theFolder+"PtRes_lepW");
    PzRes_lepW  = (TH1F *) file->Get(theFolder+"PzRes_lepW");

    mcTtMass    = (TH2F *) file->Get(theFolder+"mcTtMass");
    mcTtMass0   = (TH2F *) file->Get(theFolder+"mcTtMass0");
    mcTtMass1   = (TH2F *) file->Get(theFolder+"mcTtMass1");
    mcTtMass2   = (TH2F *) file->Get(theFolder+"mcTtMass2");
    mcTtMass3   = (TH2F *) file->Get(theFolder+"mcTtMass3");
 }

 /// Destructor
 virtual ~HTOP6() {

    delete dRb_dRw_lep;
    delete dRb_dRw_had;
    delete dRbw_had;
    delete dRbw_lep;
    delete dRtt;

    delete dR_neu;
    delete PtRes_lepW;
    delete PzRes_lepW;

    delete mcTtMass;
    delete mcTtMass0;
    delete mcTtMass1;
    delete mcTtMass2;
    delete mcTtMass3;
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
 void Fill6c( double mclepTMass, double mchadTMass ){
    mcTtMass->Fill( mclepTMass, mchadTMass );
 } 
 void Fill6c0( double mclepTMass, double mchadTMass ){
    mcTtMass0->Fill( mclepTMass, mchadTMass );
 } 
 void Fill6c1( double mclepTMass, double mchadTMass ){
    mcTtMass1->Fill( mclepTMass, mchadTMass );
 } 
 void Fill6c2( double mclepTMass, double mchadTMass ){
    mcTtMass2->Fill( mclepTMass, mchadTMass );
 } 
 void Fill6c3( double mclepTMass, double mchadTMass ){
    mcTtMass3->Fill( mclepTMass, mchadTMass );
 } 
 void Fill6d( double dR, double PtRes, double MRes ) {
    dR_neu->Fill(dR);
    PtRes_lepW->Fill(PtRes);
    PzRes_lepW->Fill(MRes);
 }


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    dRb_dRw_lep->Write();
    dRb_dRw_had->Write();
    dRbw_had->Write();
    dRbw_lep->Write();
    dRtt->Write();

    dR_neu->Write();
    PtRes_lepW->Write();
    PzRes_lepW->Write();
  
    mcTtMass->Write();
    mcTtMass0->Write();
    mcTtMass1->Write();
    mcTtMass2->Write();
    mcTtMass3->Write();
 }

  TH2F *dRb_dRw_lep;
  TH2F *dRb_dRw_had;
  TH2F *dRbw_had;
  TH2F *dRbw_lep;
  TH2F *dRtt;

  TH1F *dR_neu;
  TH1F *PtRes_lepW;
  TH1F *PzRes_lepW;

  TH2F *mcTtMass;
  TH2F *mcTtMass0;
  TH2F *mcTtMass1;
  TH2F *mcTtMass2;
  TH2F *mcTtMass3;
};

class HTOP7 {
public:
 
 HTOP7() {

    hRes_Pt   = new TH1F("hRes_Pt","Resolutio  of b jet pt ",400, -2., 2. ); 
    dR_RelPt_bmu = new TH2F("dR_RelPt_bmu","dRmin vs RelPt for b and  iso muon ",200, 0.,5., 100, 0.,200. );  
    hdRbWj = new TH1F("hdRbWj","dR for b and  W matched jets ",300, 0.,15. );  
    bDis_bj = new TH1F("bDis_bj","b Discriminator for matched b jets ", 200,-0.5,4.5 );
    bDis    = new TH1F("bDis","b Discriminator for  selected jets ", 200,-0.5,4.5 );

    bDis_dR_RelPt = new TH2F("bDis_dR_RelPt"," dR vs RelPt ",100, 0.,1., 300, 0.,60.);
    bDis_phi_md0  = new TH2F("bDis_phi_md0" ," phi vs md0  ",630, 0.,3.15, 200,0.,0.4);

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
    hEta     = new TH2F("hEta"," MC(X) vs reco(Y) , eta distribution for bjets of semiTt ",71,-3.55,3.55,71,-3.55,3.55);
    hPhi     = new TH2F("hPhi"," MC(X) vs reco(Y) , phi distribution for bjets of semiTt ",63,-3.15,3.15,71,-3.15,3.15);
    Pt_Resol = new TH2F("hPt_Resol", " Pt MC(X) vs Resol_Pt(Y) distribution for bjets of semiTt",60,0,300, 200,-0.995,1.005 );
    // 
    dR_RelPtg = new TH2F("dR_RelPtg","dR vs RelPt for b and mu - gen ",200, 0.,5., 500, 0.,50. );  
    dR_RelPtr = new TH2F("dR_RelPtr","dR vs RelPt for b and mu - reco after matching ",200, 0.,5., 500, 0.,50. );  
 } 

 HTOP7( TString theFolder, TFile* file) {

    hRes_Pt = (TH1F *) file->Get(theFolder+"hRes_Pt");
    dR_RelPt_bmu = (TH2F *) file->Get(theFolder+"dR_RelPt_bmu");
    hdRbWj = (TH1F *) file->Get(theFolder+"hdRbWj");
    bDis_bj   = (TH1F *) file->Get(theFolder+"bDis_bj");

    bDis   = (TH1F *) file->Get(theFolder+"bDis");
    bDis_dR_RelPt =(TH2F *) file->Get(theFolder+"bDis_dR_RelPt");
    bDis_phi_md0   =(TH2F *) file->Get(theFolder+"bDis_phi_md0");

    hNTk_mc = (TH1F *) file->Get(theFolder+"hNTK_mc");
    hemE    = (TH1F *) file->Get(theFolder+"hemE");
    hR60_mc = (TH2F *) file->Get(theFolder+"hR60_mc");
    hR90_mc = (TH2F *) file->Get(theFolder+"hR90_mc");
    hArea_pt_mc = (TH2F *) file->Get(theFolder+"hArea_pt_mc");
    hEovH_p_mc  = (TH2F *) file->Get(theFolder+"hEovH_p_mc");
    hEovH_n_mc  = (TH2F *) file->Get(theFolder+"hEovH_n_mc");
    hEovH_A_mc  = (TH2F *) file->Get(theFolder+"hEovH_A_mc");
    hEovH_h_mc  = (TH2F *) file->Get(theFolder+"hEovH_h_mc");
    hPt_Eta_mc  = (TH2F *) file->Get(theFolder+"hPt_Eta_mc");

    gbJpt_h = (TH2F *) file->Get(theFolder+"gbJpt_h");
    gEovH   = (TH1F *) file->Get(theFolder+"gEovH");

    hEta     = (TH2F *) file->Get(theFolder+"hEta");
    hPhi     = (TH2F *) file->Get(theFolder+"hPhi");
    Pt_Resol = (TH2F *) file->Get(theFolder+"Pt_Resol");

    dR_RelPtg = (TH2F *) file->Get(theFolder+"dR_RelPtg");
    dR_RelPtr = (TH2F *) file->Get(theFolder+"dR_RelPtr");
 }

 /// Destructor
 virtual ~HTOP7() {

    delete hRes_Pt;
    delete dR_RelPt_bmu;
    delete hdRbWj;
    delete bDis_bj;
    delete bDis;
    delete bDis_dR_RelPt;
    delete bDis_phi_md0;

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

    delete hEta;
    delete hPhi;
    delete Pt_Resol;

    delete dR_RelPtg;
    delete dR_RelPtr;
 }

 void Fill7a(double pt, double eta, double EovH, int nCon, int n60, int n90, double area, size_t nTrk, double emE, double Res_pt, double muTag ){

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
    bDis_bj->Fill( muTag );
 }
 void Fill7b( double dR_bmu , double RelPt ){
    dR_RelPt_bmu->Fill(dR_bmu, RelPt);
 }
 void Fill7c( double dR_bwj ){
    hdRbWj->Fill(dR_bwj);
 }
 void Fill7d( double dR, double RelPt, double phi, double md0 ){
    bDis_dR_RelPt->Fill( dR, RelPt );
    bDis_phi_md0->Fill( phi, md0 );
 }

 void Fill7e( double dR, double RelPt, double gdR, double gRelPt) {
    dR_RelPtr->Fill(dR, RelPt);
    dR_RelPtg->Fill(gdR, gRelPt);
 }

 void Fill7f( double pt, double eta, double EovH) {
    gbJpt_h->Fill(pt, eta);
    gEovH->Fill(EovH);
 }
 void Fill7g( double bDisValue ) {
     bDis->Fill( bDisValue );
 }
 void Fill7h(double eta_mc, double phi_mc, double eta_rc, double phi_rc, double pt, double ptResol) {
    hEta->Fill(eta_mc, eta_rc);
    hPhi->Fill(phi_mc, phi_rc);
    Pt_Resol->Fill(pt, ptResol);
 }


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    hRes_Pt->Write();
    dR_RelPt_bmu->Write();
    hdRbWj->Write();

    bDis_bj->Write();
    bDis->Write();
    bDis_dR_RelPt->Write();
    bDis_phi_md0->Write();

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

    hEta->Write();
    hPhi->Write();
    Pt_Resol->Write();
  
    dR_RelPtg->Write();
    dR_RelPtr->Write();
 }

  TH1F *hRes_Pt;
  TH2F *dR_RelPt_bmu;
  TH1F *hdRbWj;

  TH1F *bDis_bj;
  TH1F *bDis;
  TH2F *bDis_dR_RelPt;
  TH2F *bDis_phi_md0;

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

  TH2F *hEta;
  TH2F *hPhi;
  TH2F *Pt_Resol;

  TH2F *dR_RelPtg;
  TH2F *dR_RelPtr;

};

class HTOP8 {
public:
 
 HTOP8() {

    hEta_Pt1 = new TH2F("hEta_Pt1", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hEta_Pt2 = new TH2F("hEta_Pt2", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hNJets   = new TH1F("hNJets"," # of selected W jets ",65,-0.5,64.5);
    hWp_mass = new TH2F("hWp_mass"," W mass vs W momentum from selected reco jets", 500,0.,500.,240,30.,150.);

    hRes_Pt    = new TH1F("hRes_Pt", " Pt(jet) - Pt(q) / Pt(q)   ",400,-2.,2.);

    hdRWjj  = new TH2F("hdRWjj" ,"dR(Wq1, Wq2) , dR ( matched Wj1 , Wj2)  " ,200, -0.025,9.975 ,200, -0.025,9.975 );  

    bDis_dRWjMu  = new TH2F("bDis_dRWjMu","bCandidates jets JProb(x),TkCount(y)", 150,-0.1,1.4, 500,-20.,80.);

    gwJpt_h   = new TH2F("gwJpt_h"," gen w jet pt vs eta", 500,0,250, 71,-3.55,3.55);
    gEovH    = new TH1F("gEovH"," E/H for gen Jets", 500, -2., 23.);

    // parton level and best matched reco comparison
    hEta     = new TH2F("hEta"," MC(X) vs reco(Y) , eta distribution for wjets of semiTt ",71,-3.55,3.55,71,-3.55,3.55);
    hPhi     = new TH2F("hPhi"," MC(X) vs reco(Y) , phi distribution for wjets of semiTt ",63,-3.15,3.15,71,-3.15,3.15);
    Pt_Resol = new TH2F("hPt_Resol", " Pt MC(X) vs Resol_Pt(Y) distribution for wjets of semiTt",60,0,300, 200,-0.995,1.005 );
 } 

 HTOP8( TString theFolder, TFile* file ) {

    hEta_Pt1 = (TH2F *) file->Get(theFolder+"hEta_Pt1");
    hEta_Pt2 = (TH2F *) file->Get(theFolder+"hEta_Pt2");
    hNJets = (TH1F *) file->Get(theFolder+"hNJets");
    hWp_mass   = (TH2F *) file->Get(theFolder+"hWp_mass");

    hRes_Pt     = (TH1F *) file->Get(theFolder+"hRes_Pt");

    hdRWjj  = (TH2F *) file->Get(theFolder+"hdRWjj");

    bDis_dRWjMu =(TH2F *) file->Get(theFolder+"bDis_dRWjMu");

    gwJpt_h = (TH2F *) file->Get(theFolder+"gwJpt_h");
    gEovH   = (TH1F *) file->Get(theFolder+"gEovH");

    hEta     = (TH2F *) file->Get(theFolder+"hEta");
    hPhi     = (TH2F *) file->Get(theFolder+"hPhi");
    Pt_Resol = (TH2F *) file->Get(theFolder+"Pt_Resol");
 }


 /// Destructor
 virtual ~HTOP8() {

    delete hEta_Pt1;
    delete hEta_Pt2;
    delete hNJets;
    delete hWp_mass;

    delete hRes_Pt;
    delete hdRWjj;

    delete bDis_dRWjMu;

    delete gwJpt_h;
    delete gEovH;

    delete hEta;
    delete hPhi;
    delete Pt_Resol;
 }

 void Fill8a( double pt1, double pt2, double eta1, double eta2, int njets){
    hEta_Pt1->Fill(eta1, pt1);
    hEta_Pt2->Fill(eta2, pt2);
    hNJets->Fill(njets);
 }

 void Fill8b( double wmass, double wp ){
    hWp_mass->Fill(wp,wmass);
 }

 void Fill8c( double Res_Pt,  double bDis, double dR_WjMu ){
    hRes_Pt->Fill(Res_Pt);
    bDis_dRWjMu->Fill( bDis, dR_WjMu );
 }

 void Fill8f( double pt, double eta, double EovH) {
    gwJpt_h->Fill(pt, eta);
    gEovH->Fill(EovH);
 }

 void Fill8h(double eta_mc, double phi_mc, double eta_rc, double phi_rc, double pt, double ptResol) {
    hEta->Fill(eta_mc, eta_rc);
    hPhi->Fill(phi_mc, phi_rc);
    Pt_Resol->Fill(pt, ptResol);
 }

 void Fill8j( double dRjj, double dRqq ) {
    hdRWjj->Fill(dRqq, dRjj);
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    hEta_Pt1->Write();
    hEta_Pt2->Write();
    hNJets->Write();
    hWp_mass->Write();

    hRes_Pt->Write();
    hdRWjj->Write();

    bDis_dRWjMu->Write();

    gwJpt_h->Write();
    gEovH->Write();

    hEta->Write();
    hPhi->Write();
    Pt_Resol->Write();
 }

  TH2F *hEta_Pt1;
  TH2F *hEta_Pt2;
  TH1F *hNJets;
  TH2F *hWp_mass;

  TH1F *hRes_Pt;
  TH2F *hdRWjj;

  TH2F *bDis_dRWjMu;

  TH2F *gwJpt_h;
  TH1F *gEovH;

  TH2F *hEta;
  TH2F *hPhi;
  TH2F *Pt_Resol;
};


class HTOP9 {
public:
 
 HTOP9() {

   // reconstructed objects masses

    hEvtEff   = new TH1F("hEvtEff","(pass,fail)=> had(0,0.5) Mu(1,1.5) dilep(2,2.5) Ele(3,3.5), Tau(4,4.5), other(5,5.5)", 12,-0.25,5.75 );
    hObjEff   = new TH1F("hObjEff"," object selection Eff", 10, -0.5,9.5 );
    hbJetEff   = new TH2F("hbJetEff"," # of correctly matched bjet ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hWJetEff   = new TH2F("hWJetEff"," # of correctly matched Wjet ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hWLepEff   = new TH2F("hWLepEff"," # of correctly matched lep. ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hHLTBits   = new TH2F("hHLTBits","HLT Trigger bits(names) for hadronic tt", 166, -0.5, 165.5, 6,-0.5,5.5 );
    hHLTSelect = new TH1F("hHLTSelect","HLT Trigger Selected result "  , 15, -7.5, 7.5 );
    NJ_NBTags  = new TH2F("NJ_NBTags"," NJets vs NBTags", 10, 0.5, 10.5, 10, -0.5, 9.5 );

 } 

 HTOP9( TString theFolder, TFile* file ) {

    hEvtEff    = (TH1F *) file->Get(theFolder+"hEvtEff");
    hObjEff    = (TH1F *) file->Get(theFolder+"hObjEff");
    hbJetEff   = (TH2F *) file->Get(theFolder+"hbJetEff");
    hWJetEff   = (TH2F *) file->Get(theFolder+"hWJetEff");
    hWLepEff   = (TH2F *) file->Get(theFolder+"hWLepEff");
    hHLTBits   = (TH2F *) file->Get(theFolder+"hHLTBits");
    hHLTSelect = (TH1F *) file->Get(theFolder+"hHLTSelect");
    NJ_NBTags  = (TH2F *) file->Get(theFolder+"NJ_NBTags");

 }



 /// Destructor
 virtual ~HTOP9() {

    delete hEvtEff;
    delete hObjEff;
    delete hbJetEff;
    delete hWJetEff;
    delete hWLepEff;
    delete hHLTBits;
    delete hHLTSelect;
    delete NJ_NBTags;
 }

 void Fill9e( int NJ, int NBTags ) {
   NJ_NBTags->Fill( NJ, NBTags );
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


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    hEvtEff->Write();
    hObjEff->Write();
    hbJetEff->Write();
    hWJetEff->Write();
    hWLepEff->Write();
    hHLTBits->Write();
    hHLTSelect->Write();
    NJ_NBTags->Write();

 }

  TH1F *hEvtEff; 
  TH1F *hObjEff; 
  TH2F *hbJetEff;
  TH2F *hWJetEff;
  TH2F *hWLepEff;
  TH2F *hHLTBits;
  TH1F *hHLTSelect;
  TH2F *NJ_NBTags;

};
 
class HTOP10 {
public:

 HTOP10( int idx ) {

    TString N1 ;
    if ( idx == 0 ) N1 = "MC_Lep" ;
    if ( idx == 1 ) N1 = "Lep" ;
    if ( idx == 2 ) N1 = "MC_Had" ;
    if ( idx == 3 ) N1 = "Had" ;
 
    allWmass    = new TH1F(N1+"_allWmass", " mt of leptonic W or mass of hadronic W for all solutions"     ,320,0,160.);
    selWmass    = new TH1F(N1+"_selWmass", " mt of leptonic W or mass of hadronic W for selected solution" ,320,0,160.);
    beta_RelPt  = new TH2F(N1+"_RelPt", "RelPt of W w.r.t Top ",32, -0.06, 1.22, 320, 0, 160.);
    dRab_W      = new TH1F(N1+"_dRab_W","dR for (j1 , j2) or (u , v) ",200, -0.025,9.975 );  

 }

 HTOP10( TString theFolder, int idx, TFile* file ) {

    TString N1 ;
    if ( idx == 0 ) N1 = "/MC_Lep" ;
    if ( idx == 1 ) N1 = "/Lep" ;
    if ( idx == 2 ) N1 = "/MC_Had" ;
    if ( idx == 3 ) N1 = "/Had" ;

    allWmass   = (TH1F *) file->Get(theFolder+N1+"_allWmass");
    selWmass   = (TH1F *) file->Get(theFolder+N1+"_selWmass");
    beta_RelPt = (TH2F *) file->Get(theFolder+N1+"_beta_RelPt");
    dRab_W     = (TH1F *) file->Get(theFolder+N1+"_dRab_W");

 }

 /// Destructor
 virtual ~HTOP10() {

    delete allWmass;
    delete selWmass;
    delete beta_RelPt;
    delete dRab_W;
 }

 void Fill10a( double wmass ){
    allWmass->Fill( wmass );
 }
 void Fill10b( double wmass, double RelPt, double beta, double dRab ) {
    selWmass->Fill( wmass );
    beta_RelPt->Fill( beta, RelPt );
    dRab_W->Fill( dRab );
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    allWmass->Write();
    selWmass->Write();
    beta_RelPt->Write();
    dRab_W->Write();

 }

  TH1F *allWmass;
  TH1F *selWmass;
  TH2F *beta_RelPt;
  TH1F *dRab_W;

 };

class HTOP11 {
public:

 HTOP11( int idx ) {

    TString N1 ;
    if ( idx == 0 ) N1 = "MC_Lep" ;
    if ( idx == 1 ) N1 = "Lep" ;
    if ( idx == 2 ) N1 = "MC_Had" ;
    if ( idx == 3 ) N1 = "Had" ;
    TString N2 ;
    if ( idx == 0 ) N2 = "MC_" ;
    if ( idx == 1 ) N2 = "" ;

    allTmass_pt    = new TH2F(N1+"_allTmass_pt",   " Top Mass vs Pt for all solutions" ,480,0.,480., 480,0,480 );
    highPtTmass_pt = new TH2F(N1+"_highPtTmass_pt"," highest Pt Top Mass vs Pt"        ,480,0.,480., 480,0,480 );
    selTmass_Wmass = new TH2F(N1+"_selTmass_Wmass"," Top Mass vs W Mass for seleted solution" ,480,0.,480., 480,0,480 );
    selTmass_pt    = new TH2F(N1+"_selTmass_pt",   " Top Mass vs Pt for seleted solution >3J" ,480,0.,480., 480,0,480 );
    selTmass_pt0   = new TH2F(N1+"_selTmass_pt0",  " Top Mass vs Pt for seleted solution =4J" ,480,0.,480., 480,0,480 );
    selTmass_pt1   = new TH2F(N1+"_selTmass_pt1", " Top Mass vs Pt for seleted solution =5J"  ,480,0.,480., 480,0,480 );
    selTmass_pt2   = new TH2F(N1+"_selTmass_pt2", " Top Mass vs Pt for seleted solution =6J"  ,480,0.,480., 480,0,480 );
    selTmass_pt3   = new TH2F(N1+"_selTmass_pt3", " Top Mass vs Pt for seleted solution >6J"  ,480,0.,480., 480,0,480 );

    hRecott0   = new TH2F(N2+"hRecott0","Reco top mass w/ 4j", 480,0.,480., 480,0.,480.);
    hRecott1   = new TH2F(N2+"hRecott1","Reco top mass w/ 5j", 480,0.,480., 480,0.,480.);
    hRecott2   = new TH2F(N2+"hRecott2","Reco top mass w/ 6j", 480,0.,480., 480,0.,480.);
    hRecott3   = new TH2F(N2+"hRecott3","Reco top mass >= 7j", 480,0.,480., 480,0.,480.);
    hRecott    = new TH2F(N2+"hRecott" ,"Reco top mass -all ", 480,0.,480., 480,0.,480.);

    PtRecott0   = new TH1F(N2+"PtRecott0","Pt of Reco top mass w/ 4j", 480,0.,480.);
    PtRecott1   = new TH1F(N2+"PtRecott1","Pt of Reco top mass w/ 5j", 480,0.,480.);
    PtRecott2   = new TH1F(N2+"PtRecott2","Pt of Reco top mass w/ 6j", 480,0.,480.);
    PtRecott3   = new TH1F(N2+"PtRecott3","Pt of Reco top mass >= 7j", 480,0.,480.);
    PtRecott    = new TH1F(N2+"PtRecott" ,"Pt of Reco top mass -all ", 480,0.,480.);

    WRecott0   = new TH2F(N2+"WRecott0","W vs dM of Reco top mass w/ 4j", 500,0,1000, 400,-199.5,200.5);
    WRecott1   = new TH2F(N2+"WRecott1","W vs dM of Reco top mass w/ 5j", 500,0,1000, 400,-199.5,200.5);
    WRecott2   = new TH2F(N2+"WRecott2","W vs dM of Reco top mass w/ 6j", 500,0,1000, 400,-199.5,200.5);
    WRecott3   = new TH2F(N2+"WRecott3","W vs dM of Reco top mass >= 7j", 500,0,1000, 400,-199.5,200.5);
    WRecott    = new TH2F(N2+"WRecott" ,"W vs dM of Reco top mass  -all", 500,0,1000, 400,-199.5,200.5);

 }

 HTOP11( TString theFolder, int idx, TFile* file ) {

    TString N1 ;
    if ( idx == 0 ) N1 = "/MC_Lep" ;
    if ( idx == 1 ) N1 = "/Lep" ;
    if ( idx == 2 ) N1 = "/MC_Had" ;
    if ( idx == 3 ) N1 = "/Had" ;
    TString N2 ;
    if ( idx == 0 ) N2 = "/MC_" ;
    if ( idx == 1 ) N2 = "/" ;

    allTmass_pt    = (TH2F *) file->Get(theFolder+N1+"_allTmass_pt");
    highPtTmass_pt = (TH2F *) file->Get(theFolder+N1+"_highPtTmass_pt");
    selTmass_Wmass = (TH2F *) file->Get(theFolder+N1+"_selTmass_Wmass");
    selTmass_pt    = (TH2F *) file->Get(theFolder+N1+"_selTmass_pt");
    selTmass_pt0   = (TH2F *) file->Get(theFolder+N1+"_selTmass_pt0");
    selTmass_pt1   = (TH2F *) file->Get(theFolder+N1+"_selTmass_pt1");
    selTmass_pt2   = (TH2F *) file->Get(theFolder+N1+"_selTmass_pt2");
    selTmass_pt3   = (TH2F *) file->Get(theFolder+N1+"_selTmass_pt3");

    hRecott0   = (TH2F *) file->Get(theFolder+N2+"_hRecott0");
    hRecott1   = (TH2F *) file->Get(theFolder+N2+"_hRecott1");
    hRecott2   = (TH2F *) file->Get(theFolder+N2+"_hRecott2");
    hRecott3   = (TH2F *) file->Get(theFolder+N2+"_hRecott3");
    hRecott    = (TH2F *) file->Get(theFolder+N2+"_hRecott");
    PtRecott0  = (TH1F *) file->Get(theFolder+N2+"_PtRecott0");
    PtRecott1  = (TH1F *) file->Get(theFolder+N2+"_PtRecott1");
    PtRecott2  = (TH1F *) file->Get(theFolder+N2+"_PtRecott2");
    PtRecott3  = (TH1F *) file->Get(theFolder+N2+"_PtRecott3");
    PtRecott   = (TH1F *) file->Get(theFolder+N2+"_PtRecott");
    WRecott0   = (TH2F *) file->Get(theFolder+N2+"_WRecott0");
    WRecott1   = (TH2F *) file->Get(theFolder+N2+"_WRecott1");
    WRecott2   = (TH2F *) file->Get(theFolder+N2+"_WRecott2");
    WRecott3   = (TH2F *) file->Get(theFolder+N2+"_WRecott3");
    WRecott    = (TH2F *) file->Get(theFolder+N2+"_WRecott");
 }

 /// Destructor
 virtual ~HTOP11() {

    delete allTmass_pt;
    delete highPtTmass_pt;
    delete selTmass_Wmass;
    delete selTmass_pt;
    delete selTmass_pt0;
    delete selTmass_pt1;
    delete selTmass_pt2;
    delete selTmass_pt3;

    delete hRecott0;
    delete hRecott1;
    delete hRecott2;
    delete hRecott3;
    delete hRecott;
    delete PtRecott0;
    delete PtRecott1;
    delete PtRecott2;
    delete PtRecott3;
    delete PtRecott;
    delete WRecott0;
    delete WRecott1;
    delete WRecott2;
    delete WRecott3;
    delete WRecott;
 }

 void Fill11a( double tmass, double pt ){
    allTmass_pt->Fill( tmass, pt );
 }
 void Fill11b( double tmass, double pt ){
    highPtTmass_pt->Fill( tmass, pt );
 }
 void Fill11c( double tmass, double pt, double wmass ) {
    selTmass_pt->Fill( tmass, pt );
    selTmass_Wmass->Fill( tmass, wmass );
 }
 void Fill11c0( double tmass, double pt) {
    selTmass_pt0->Fill( tmass, pt );
 }
 void Fill11c1( double tmass, double pt) {
    selTmass_pt1->Fill( tmass, pt );
 }
 void Fill11c2( double tmass, double pt) {
    selTmass_pt2->Fill( tmass, pt );
 }
 void Fill11c3( double tmass, double pt) {
    selTmass_pt3->Fill( tmass, pt );
 }

 void Fill11d(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott->Fill( lep_mt, had_mt );
    PtRecott->Fill( PtTt );
    WRecott->Fill( Wtt, Wt2 );
 }
 void Fill11d0(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott0->Fill( lep_mt, had_mt );
    PtRecott0->Fill( PtTt );
    WRecott0->Fill( Wtt, Wt2 );
 }
 void Fill11d1(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott1->Fill( lep_mt, had_mt );
    PtRecott1->Fill( PtTt );
    WRecott1->Fill( Wtt, Wt2 );
 }
 void Fill11d2(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott2->Fill( lep_mt, had_mt );
    PtRecott2->Fill( PtTt );
    WRecott2->Fill( Wtt, Wt2 );
 }
 void Fill11d3(double lep_mt, double had_mt, double PtTt, double Wtt, double Wt2 ){
    hRecott3->Fill( lep_mt, had_mt );
    PtRecott3->Fill( PtTt );
    WRecott3->Fill( Wtt, Wt2 );
 }

 void Write( TString theFolder , TFile* file, int idx  ) {

    file->cd( theFolder );

    allTmass_pt->Write();
    highPtTmass_pt->Write();
    selTmass_Wmass->Write();
    selTmass_pt->Write();
    selTmass_pt0->Write();
    selTmass_pt1->Write();
    selTmass_pt2->Write();
    selTmass_pt3->Write();

    if ( idx == 0 || idx == 1 ) {
       hRecott0->Write();
       hRecott1->Write();
       hRecott2->Write();
       hRecott3->Write();
       hRecott->Write();
       PtRecott0->Write();
       PtRecott1->Write();
       PtRecott2->Write();
       PtRecott3->Write();
       PtRecott->Write();
       WRecott0->Write();
       WRecott1->Write();
       WRecott2->Write();
       WRecott3->Write();
       WRecott->Write();
    }
 }

  TH2F *allTmass_pt;
  TH2F *highPtTmass_pt;
  TH2F *selTmass_Wmass;
  TH2F *selTmass_pt;
  TH2F *selTmass_pt0;
  TH2F *selTmass_pt1;
  TH2F *selTmass_pt2;
  TH2F *selTmass_pt3;

  TH2F *hRecott0;
  TH2F *hRecott1;
  TH2F *hRecott2;
  TH2F *hRecott3;
  TH2F *hRecott;
  TH1F *PtRecott0;
  TH1F *PtRecott1;
  TH1F *PtRecott2;
  TH1F *PtRecott3;
  TH1F *PtRecott;
  TH2F *WRecott0;
  TH2F *WRecott1;
  TH2F *WRecott2;
  TH2F *WRecott3;
  TH2F *WRecott;

 };


#endif
