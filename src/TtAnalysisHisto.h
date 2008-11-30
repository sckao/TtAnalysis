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
    Uncorr_Et_Pt = new TH2F("Uncorr_Et_Pt", "Uncorr Et vs Pt (all jets)", 500, 0., 500., 500, 0., 500);
    
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

    thirdJetEt = new TH1F("thirdJetEt", " ET of the 3rd highest Et jet ", 500,0.,500. ); 
    thirdCalEt = new TH1F("thirdCalEt", " ET of the 3rd highest Et jet ", 500,0.,500. ); 
    j1Et_Area  = new TH2F("j1Et_Area", " 1st Jet Et vs towerArea ", 500, 0., 500., 100, 0., 1.);
    j3Et_Area  = new TH2F("j3Et_Area", " 3rd Jet Et vs towerArea ", 500, 0., 500., 100, 0., 1.);
    j3Cal_Area = new TH2F("j3Cal_Area", " 3rd Jet CalEt vs towerArea ", 500, 0., 500., 100, 0., 1.);
    m3cutEt    = new TH1F("m3cutEt", " 3rd highest Et jet w/ m3 < 150 ", 500,0.,500. );
    etcutM3    = new TH1F("etcutM3", " m3 w/ 3rd Jet Et < 20 ", 500, 0., 500. );

    eta_njets  = new TH2F("eta_njets", " eta vs njets ", 91, -4.55, 4.55, 25, -0.5, 24.5);
    wEta_njets = new TH2F("wEta_njets", "w eta vs njets ", 91, -4.55, 4.55, 25, -0.5, 24.5);

    wEta_dPhi0 = new TH2F("wEta_dPhi0", " w eta vs dPhi(w,jet1)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    wEta_dPhi1 = new TH2F("wEta_dPhi1", " w eta vs dPhi(w,jet2)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    wEta_dPhi2 = new TH2F("wEta_dPhi2", " w eta vs dPhi(w,jet3)  ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    njets_dPhi0= new TH2F("njets_dPhi0", " njets vs dPhi(w,jet1)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    njets_dPhi1= new TH2F("njets_dPhi1", " njets vs dPhi(w,jet2)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    njets_dPhi2= new TH2F("njets_dPhi2", " njets vs dPhi(w,jet3)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    
    muEta_dPhi0  = new TH2F("muEta_dPhi0",  " mu eta vs dPhi(mu,jet1) ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    muEta_dPhi1  = new TH2F("muEta_dPhi1",  " mu eta vs dPhi(mu,jet2) ", 91, -4.55, 4.55, 160, -0.05, 3.15);
    muEta_dPhi2  = new TH2F("muEta_dPhi2",  " mu eta vs dPhi(mu,jet3) ", 91, -4.55, 4.55, 160, -0.05, 3.15);
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

    dR_Ht2  = new TH2F("dR_Ht2", " dR(w,jet1+2+3) vs Ht for Nj<=2", 100, -0.05, 9.95, 500, 0., 500.);
    dR_Ht3  = new TH2F("dR_Ht3", " dR(w,jet1+2+3) vs Ht for Nj>=3", 100, -0.05, 9.95, 500, 0., 500.);

    dF_Ht2  = new TH2F("dF_Ht2", " dF(w,jet1+2+3) vs Ht for Nj<=2", 40, -0.05, 3.95, 500, 0., 500.);
    dF_Ht3  = new TH2F("dF_Ht3", " dF(w,jet1+2+3) vs Ht for Nj>=3", 40, -0.05, 3.95, 500, 0., 500.);

    dR0_et3 = new TH2F("dR0_et3", " dR(w,jet1) vs J3Et for Nj > 2 ", 100, -0.05, 9.95, 500, 0., 500.);
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
    Uncorr_Et_Pt  = (TH2F *) file->Get("Uncorr_Et_Pt");
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

    thirdJetEt = (TH1F *) file->Get("thirdJetEt");
    thirdCalEt = (TH1F *) file->Get("thirdCalEt");
    j3Et_Area  = (TH2F *) file->Get("j3Et_Area");
    j1Et_Area  = (TH2F *) file->Get("j1Et_Area");
    m3cutEt    = (TH1F *) file->Get("m3cutEt");
    etcutM3    = (TH1F *) file->Get("etcutM3");
    j3Cal_Area = (TH2F *) file->Get("j3Cal_Area");
    eta_njets  = (TH2F *) file->Get("eta_njets");
    wEta_njets  = (TH2F *) file->Get("wEta_njets");

    wEta_dPhi0 = (TH2F *) file->Get("wEta_dPhi0");
    wEta_dPhi1 = (TH2F *) file->Get("wEta_dPhi1");
    wEta_dPhi2 = (TH2F *) file->Get("wEta_dPhi2");
    njets_dPhi0= (TH2F *) file->Get("njets_dPhi0");
    njets_dPhi1= (TH2F *) file->Get("njets_dPhi1");
    njets_dPhi2= (TH2F *) file->Get("njets_dPhi2");
 
    muEta_dPhi0  = (TH2F *) file->Get("muEta_dPhi0");
    muNjets_dPhi0= (TH2F *) file->Get("muNjets_dPhi0");
    muEta_dPhi1  = (TH2F *) file->Get("muEta_dPhi1");
    muNjets_dPhi1= (TH2F *) file->Get("muNjets_dPhi1");
    muEta_dPhi2  = (TH2F *) file->Get("muEta_dPhi2");
    muNjets_dPhi2= (TH2F *) file->Get("muNjets_dPhi2");

    njets_dR0= (TH2F *) file->Get("njets_dR0");
    njets_dR1= (TH2F *) file->Get("njets_dR1");
    njets_dR2= (TH2F *) file->Get("njets_dR2");

    njets_dR = (TH2F *) file->Get("njets_dR");
    m3_dR    = (TH2F *) file->Get("m3_dR");
    m3_j3    = (TH2F *) file->Get("m3_j3");
    njets_m3 = (TH2F *) file->Get("njets_m3");

    dR_Ht2 = (TH2F *) file->Get("dR_Ht2");
    dR_Ht3 = (TH2F *) file->Get("dR_Ht3");
    dF_Ht2 = (TH2F *) file->Get("dF_Ht2");
    dF_Ht3 = (TH2F *) file->Get("dF_Ht3");

    dR0_et3 = (TH2F *) file->Get("dR0_et3");
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
    delete Uncorr_Et_Pt;
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

    delete thirdJetEt; 
    delete thirdCalEt;
    delete m3cutEt;
    delete etcutM3;
    delete j3Et_Area;
    delete j1Et_Area;
    delete j3Cal_Area;
    delete eta_njets;
    delete wEta_njets;

    delete wEta_dPhi0;
    delete wEta_dPhi1;
    delete wEta_dPhi2;
    delete njets_dPhi0;
    delete njets_dPhi1;
    delete njets_dPhi2;

    delete muEta_dPhi0;
    delete muNjets_dPhi0;
    delete muEta_dPhi1;
    delete muNjets_dPhi1;
    delete muEta_dPhi2;
    delete muNjets_dPhi2;

    delete njets_dR0;
    delete njets_dR1;
    delete njets_dR2;

    delete njets_dR;
    delete m3_dR;
    delete m3_j3;
    delete njets_m3;

    delete dR_Ht2;
    delete dR_Ht3;
    delete dF_Ht2;
    delete dF_Ht3;

    delete dR0_et3;

 }

 void Fill1a( double pt, double et, double EovH, int nCont, double r60, double uncorrPt, double uncorrEt, double r , double area, int nTk, double emF, double emFCalo ){
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
    Uncorr_Et_Pt->Fill(uncorrEt, uncorrPt);
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
 void Fill1e( double m3, double et ) {
    etcutM3->Fill( m3 );
    m3cutEt->Fill( et );
 }
 void Fill1f( int njets, int nj20, int nj40, int njh35 ){
    hNJets->Fill(njets);  
    hNjet20->Fill(nj20);  
    hNjet30->Fill(nj40);  
    hNjet40->Fill(njh35);  
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
 void Fill1l( double et3 , double calEt, double towerArea3, double et1, double towerArea1, double m3 ) {
    thirdJetEt->Fill( et3 );
    thirdCalEt->Fill( calEt );
    j3Cal_Area->Fill( calEt, towerArea3);
    j3Et_Area->Fill( et3, towerArea3);
    j1Et_Area->Fill( et1, towerArea1);
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
 void Fill1o( double njets, double eta, double dphi0, double dphi1, double dphi2 ) {
    muEta_dPhi0->Fill( eta, dphi0 );
    muNjets_dPhi0->Fill( njets, dphi0 );
    muEta_dPhi1->Fill( eta, dphi1 );
    muNjets_dPhi1->Fill( njets, dphi1 );
    muEta_dPhi2->Fill( eta, dphi2 );
    muNjets_dPhi2->Fill( njets, dphi2 );
 } 
 void Fill1p( double njets, double eta ) {
    wEta_njets->Fill( eta, njets );
 }

 void Fill1q2( double dR, double HT, double dF ) {
    dR_Ht2->Fill(dR, HT);
    dF_Ht2->Fill(dF, HT);
 }
 void Fill1q3( double dR, double HT, double dF ) {
    dR_Ht3->Fill(dR, HT);
    dF_Ht3->Fill(dF, HT);
 }
 void Fill1r( double dR0 , double et3 ) {
    dR0_et3->Fill(dR0, et3) ;
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
    Uncorr_Et_Pt->Write();
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

    thirdJetEt->Write();
    thirdCalEt->Write();
    m3cutEt->Write();
    etcutM3->Write();
    j3Et_Area->Write();
    j1Et_Area->Write();
    j3Cal_Area->Write();
    eta_njets->Write();
    wEta_njets->Write();

    wEta_dPhi0->Write();
    wEta_dPhi1->Write();
    wEta_dPhi2->Write();
    njets_dPhi0->Write();
    njets_dPhi1->Write();
    njets_dPhi2->Write();
    
    muEta_dPhi0->Write();
    muNjets_dPhi0->Write();
    muEta_dPhi1->Write();
    muNjets_dPhi1->Write();
    muEta_dPhi2->Write();
    muNjets_dPhi2->Write();

    njets_dR0->Write();
    njets_dR1->Write();
    njets_dR2->Write();

    njets_dR->Write();
    m3_dR->Write();
    m3_j3->Write();
    njets_m3->Write();

    dR_Ht2->Write();
    dR_Ht3->Write();
    dF_Ht2->Write();
    dF_Ht3->Write();

    dR0_et3->Write();
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
  TH2F *Uncorr_Et_Pt;
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

  TH1F *thirdJetEt;
  TH1F *thirdCalEt;
  TH1F *m3cutEt;
  TH1F *etcutM3;
  TH2F *j3Et_Area;
  TH2F *j1Et_Area;
  TH2F *j3Cal_Area;
  TH2F *eta_njets;
  TH2F *wEta_njets;

  TH2F *wEta_dPhi0;
  TH2F *wEta_dPhi1;
  TH2F *wEta_dPhi2;
  TH2F *njets_dPhi0;
  TH2F *njets_dPhi1;
  TH2F *njets_dPhi2;

  TH2F *muEta_dPhi0;
  TH2F *muNjets_dPhi0;
  TH2F *muEta_dPhi1;
  TH2F *muNjets_dPhi1;
  TH2F *muEta_dPhi2;
  TH2F *muNjets_dPhi2;

  TH2F *njets_dR0;
  TH2F *njets_dR1;
  TH2F *njets_dR2;
  TH2F *njets_dR;
  TH2F *m3_dR;
  TH2F *m3_j3;
  TH2F *njets_m3;

  TH2F *dR_Ht2;
  TH2F *dR_Ht3;
  TH2F *dF_Ht2;
  TH2F *dF_Ht3;
 
  TH2F *dR0_et3;

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

 void Write() {

    hemfMET->Write();
    hSumEt_MET->Write();

    hResPAT->Write();
    hResCal->Write();
    hResCor->Write();
    hPhiPAT->Write();
    hPhiCal->Write();
    hPhiCor->Write();

 }

  TH2F *hemfMET;
  TH2F *hSumEt_MET;

  TH1F *hResPAT;
  TH1F *hResCal;
  TH1F *hResCor;
  TH1F *hPhiPAT;
  TH1F *hPhiCal;
  TH1F *hPhiCor;

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

   // Isolation
    hemEt   = new TH2F("hemEt","emEt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hhadEt  = new TH2F("hhadEt","hadEt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hNCalo  = new TH2F("hNCalo","hNCalo: (R3x, R5y)",40,0.,40.,40,0.,40 );
    hnTrack = new TH2F("hnTrack","hnTrack: (R3x, R5y)",40,0.,40.,40,0.,40 );
    hsumPt  = new TH2F("hsumPt","hsumPt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hsumEt  = new TH2F("hsumEt","hsumEt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hcaloIso= new TH1F("hcaloIso"," caloIso value", 75,-0.5,14.5 );
    hcaloIso2D= new TH2F("hcaloIso2D"," ECaloIso(x) , HCaloIso(y)", 75,-0.5,14.5, 75,-0.5,14.5 );
    hIso_Pt = new TH2F("hIso_Pt"," IsoValue(x) , Pt(y)", 120,-0.1,1.1, 500,0,500);

    hemEt2   = new TH2F("hemEt2","IsoSelected emEt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hhadEt2  = new TH2F("hhadEt2","IsoSelected hadEt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hNCalo2   = new TH2F("hNCalo2","IsoSelected hNCalo: (R3x, R5y)",40,0.,40.,40,0.,40 );
    hnTrack2 = new TH2F("hnTrack2","IsoSelected hnTrack: (R3x, R5y)",40,0.,40.,40,0.,40 );
    hsumPt2  = new TH2F("hsumPt2","IsoSelected hsumPt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hsumEt2  = new TH2F("hsumEt2","IsoSelected hsumEt: (R3x, R5y)",400,0.,20.,400,0.,20 );
    hcaloIso_2 = new TH1F("hcaloIso_2","IsoSelected caloIso value", 75,-0.5,14.5 );
    hcaloIso2D_2= new TH2F("hcaloIso2D_2","IsoSelected ECalIso(x) , HCalIso(y)", 75,-0.5,14.5, 75,-0.5,14.5 );
    hIso_Pt2 = new TH2F("hIso_Pt2","IsoSelected IsoValue(x) , Pt(y)", 120,-0.1,1.1, 500,0,500);

    hIMu_caloE = new TH2F("hIMu_caloE","IsoMu caloE vs muon p ",800,0.,200., 800,0.,200. );
    hJMu_caloE = new TH2F("hJMu_caloE","JetMu caloE vs muon p ",800,0.,200., 800,0.,200. );
    hIMu_caloE_jE = new TH2F("hIMu_caloE_jE","IsoMu caloE vs jet E ",800,0.,200., 800,0.,200. );
    hJMu_caloE_jE = new TH2F("hJMu_caloE_jE","JetMu caloE vs jet E ",800,0.,200., 800,0.,200. );

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

    hemEt  = (TH2F *) file->Get("hemEt"); 
    hhadEt = (TH2F *) file->Get("hhadEt"); 
    hNCalo  = (TH2F *) file->Get("hNCalo"); 
    hnTrack = (TH2F *) file->Get("hnTrack"); 
    hsumPt  = (TH2F *) file->Get("hsumPt"); 
    hsumEt  = (TH2F *) file->Get("hsumEt"); 
    hcaloIso = (TH1F *) file->Get("hcaloIso");
    hcaloIso2D = (TH2F *) file->Get("hcaloIso2D");
    hIso_Pt = (TH2F *) file->Get("hIso_Pt");
    
    hemEt2  = (TH2F *) file->Get("hemEt2"); 
    hhadEt2 = (TH2F *) file->Get("hhadEt2"); 
    hNCalo2  = (TH2F *) file->Get("hNCalo2"); 
    hnTrack2 = (TH2F *) file->Get("hnTrack2"); 
    hsumPt2  = (TH2F *) file->Get("hsumPt2"); 
    hsumEt2  = (TH2F *) file->Get("hsumEt2"); 
    hcaloIso_2 = (TH1F *) file->Get("hcaloIso_2");
    hcaloIso2D_2 = (TH2F *) file->Get("hcaloIso2D_2");
    hIso_Pt2 = (TH2F *) file->Get("hIso_Pt2");

    hIMu_caloE = (TH2F *) file->Get("hIMu_caloE");
    hJMu_caloE = (TH2F *) file->Get("hJMu_caloE");
    hIMu_caloE_jE = (TH2F *) file->Get("hIMu_caloE_jE");
    hJMu_caloE_jE = (TH2F *) file->Get("hJMu_caloE_jE");
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
 
    delete hemEt;
    delete hhadEt;
    delete hNCalo;
    delete hnTrack;
    delete hsumPt;
    delete hsumEt;
    delete hcaloIso; 
    delete hcaloIso2D; 
    delete hIso_Pt; 

    delete hemEt2;
    delete hhadEt2;
    delete hNCalo2;
    delete hnTrack2;
    delete hsumPt2;
    delete hsumEt2;
    delete hcaloIso_2; 
    delete hcaloIso2D_2; 
    delete hIso_Pt2; 
 
    delete hIMu_caloE;
    delete hJMu_caloE;
    delete hIMu_caloE_jE;
    delete hJMu_caloE_jE;
 }

 void Fill3a(double Pt,double Eta ){
    hPt->Fill(Pt);
    hEta->Fill(Eta);
 }
 
 void Fill3b(double emEt3, double emEt5, double hdEt3, double hdEt5, int nCal3, int nCal5, 
             int ntrack3, int ntrack5, double sumPt3, double sumPt5, double sumEt3, double sumEt5,
             double caloIso, double ecalIso, double hcalIso, double Iso, double pt ) {

    hemEt->Fill(emEt3,emEt5);
    hhadEt->Fill(hdEt3,hdEt5);
    hNCalo->Fill(nCal3,nCal5);
    hnTrack->Fill(ntrack3,ntrack5);
    hsumPt->Fill(sumPt3,sumPt5);
    hsumEt->Fill(sumEt3,sumEt5);
    hcaloIso->Fill(caloIso);
    hcaloIso2D->Fill(ecalIso,hcalIso);
    hIso_Pt->Fill(Iso, pt);
 }
 
 void Fill3c(double emEt3, double emEt5, double hdEt3, double hdEt5, int nCal3, int nCal5,
             int ntrack3, int ntrack5, double sumPt3, double sumPt5, double sumEt3, double sumEt5,
             double caloIso, double ecalIso, double hcalIso, double Iso, double pt ) {

    hemEt2->Fill(emEt3,emEt5);
    hhadEt2->Fill(hdEt3,hdEt5);
    hNCalo2->Fill(nCal3,nCal5);
    hnTrack2->Fill(ntrack3,ntrack5);
    hsumPt2->Fill(sumPt3,sumPt5);
    hsumEt2->Fill(sumEt3,sumEt5);
    hcaloIso_2->Fill(caloIso);
    hcaloIso2D_2->Fill(ecalIso,hcalIso);
    hIso_Pt2->Fill(Iso, pt);
 }

 void Fill3d(double isoMuE, double mu_p, double jet_E ) {
    hIMu_caloE->Fill(isoMuE, mu_p);
    hIMu_caloE_jE->Fill(isoMuE, jet_E);
 }
 void Fill3e(double jetMuE, double mu_p, double jet_E ) {
    hJMu_caloE->Fill(jetMuE, mu_p);
    hJMu_caloE_jE->Fill(jetMuE, jet_E);
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
 
    hemEt->Write();
    hhadEt->Write();
    hNCalo->Write();
    hnTrack->Write();
    hsumPt->Write();
    hsumEt->Write();
    hcaloIso->Write();
    hcaloIso2D->Write();
    hIso_Pt->Write();

    hemEt2->Write();
    hhadEt2->Write();
    hNCalo2->Write();
    hnTrack2->Write();
    hsumPt2->Write();
    hsumEt2->Write();
    hcaloIso_2->Write();
    hcaloIso2D_2->Write();
    hIso_Pt2->Write();

    hIMu_caloE->Write();
    hJMu_caloE->Write();
    hIMu_caloE_jE->Write();
    hJMu_caloE_jE->Write();

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

  TH2F *hemEt;
  TH2F *hhadEt;
  TH2F *hNCalo;
  TH2F *hnTrack;
  TH2F *hsumPt;
  TH2F *hsumEt;
  TH1F *hcaloIso;
  TH2F *hcaloIso2D;
  TH2F *hIso_Pt;

  TH2F *hemEt2;
  TH2F *hhadEt2;
  TH2F *hNCalo2;
  TH2F *hnTrack2;
  TH2F *hsumPt2;
  TH2F *hsumEt2;
  TH1F *hcaloIso_2;
  TH2F *hcaloIso2D_2;
  TH2F *hIso_Pt2;

  TH2F *hIMu_caloE;
  TH2F *hJMu_caloE;
  TH2F *hIMu_caloE_jE;
  TH2F *hJMu_caloE_jE;

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
    hEta1a  = new TH1F("hEta1a"," eta distribution for electron of semi Tt(parton) ",71,-3.55,3.55);
    hEta1b  = new TH1F("hEta1b"," eta distribution for electron of semi Tt(mcMatched)",71,-3.55,3.55);
    hEta1c  = new TH1F("hEta1c"," eta distribution for electron of semi Tt(recoUsed)",71,-3.55,3.55);
    hPt1a   = new TH1F("hPt1a", " Pt distribution for e of semi Tt(parton)",30,0,300 );
    hPt1b   = new TH2F("hPt1b", " Pt vs Pt_Resolution for e of semi Tt(mcMatched)",30,0,300, 31,-1.55,1.55 );
    hPt1c   = new TH2F("hPt1c", " Pt vs Pt_Resolution for e of semi Tt(recoUsed)",30,0,300, 31,-1.55,1.55 );
   // Isolation
    hPt_em1  = new TH2F("hPt_em1", "pt vs iso emE for all",800,0.,200., 400, 0., 40.);
    hPt_hd1  = new TH2F("hPt_hd1", "pt vs iso hdE for all",800,0.,200., 400, 0., 40.);
    hPt_tk1  = new TH2F("hPt_tk1", "pt vs iso tkP for all",800,0.,200., 400, 0., 40.);
    hPt_sum1 = new TH2F("hPt_sum1","pt vs iso caloE for all",800,0.,200., 400, 0., 40.);
    hPt_Cal1 = new TH2F("hPt_Cal1","pt vs # of CaloDeposit in isoR3 for all ",800,0.,200., 11, -0.5, 10.5 );
    hPt_Trk1 = new TH2F("hPt_Trk1","pt vs # of Track in isoR3 for all ",800,0.,200., 11, -0.5, 10.5 );
    hPt_Iso1 = new TH2F("hPt_Iso1","pt vs Iso value for all",800,0.,200., 200, -0.5, 1.5 );
    hEta_Iso1= new TH2F("hEta_Iso1", "eta vs Iso for all ",61, -3.05,3.05 , 200, -0.5, 1.5);  
    hEta_trkIso1  = new TH2F("hEta_trkIso1", "eta vs trkIso for all ",61, -3.05,3.05 , 30, -0.25,14.75);  
    hEta_ecalIso1 = new TH2F("hEta_ecalIso1","eta vs ecalIso for all",61, -3.05,3.05 , 30, -0.25,14.75);  
    hEta_hcalIso1 = new TH2F("hEta_hcalIso1","eta vs hcalIso for all",61, -3.05,3.05 , 30, -0.25,14.75);  

    hPt_em2  = new TH2F("hPt_em2", "pt vs iso emE after cut",800,0.,200., 400, 0., 40.);
    hPt_hd2  = new TH2F("hPt_hd2", "pt vs iso hdE after cut",800,0.,200., 400, 0., 40.);
    hPt_tk2  = new TH2F("hPt_tk2", "pt vs iso tkP for all",800,0.,200., 400, 0., 40.);
    hPt_sum2 = new TH2F("hPt_sum2","pt vs iso caloE after cut",800,0.,200., 400, 0., 40.);
    hPt_Cal2 = new TH2F("hPt_Cal2","pt vs # of CaloDeposit in isoR3 after cut",800,0.,200., 11, -0.5, 10.5 );
    hPt_Trk2 = new TH2F("hPt_Trk2","pt vs # of Track in isoR3 after cut",800,0.,200., 11, -0.5, 10.5 );
    hPt_Iso2 = new TH2F("hPt_Iso2","pt vs Iso value after cut",800,0.,200., 200, -0.5, 1.5 );
    hEta_trkIso2  = new TH2F("hEta_trkIso2", "eta vs trkIso after cut",61, -3.05,3.05 , 30, -0.25,14.75);  
    hEta_ecalIso2 = new TH2F("hEta_ecalIso2","eta vs ecalIso after cut",61, -3.05,3.05 , 30, -0.25,14.75);  
    hEta_hcalIso2 = new TH2F("hEta_hcalIso2","eta vs hcalIso after cut",61, -3.05,3.05 , 30, -0.25,14.75);  

   // MC matching electrons
    hPt_EovP = new TH2F("hPt_EovP","pt vs emE/P  ",800,0.,200., 500, 0.25, 1.25 );
    hPt_HovE = new TH2F("hPt_HovE","pt vs H/E ",800,0.,200., 500, -0.25, 0.75 );
    hPt_nTrk = new TH2F("hPt_nTrk","pt vs # of track in isoR3 ",800,0.,200., 11, -0.5, 10.5 );
    hPt_nCal = new TH2F("hPt_nCal","pt vs # of CaloDeposit in isoR3 ",800,0.,200., 11, -0.5, 10.5 );
    hPt_Iso  = new TH2F("hPt_Iso","pt vs Iso value - MC matched ",800,0.,200., 200, -0.5, 1.5);
    hEta_Iso = new TH2F("hEta_Iso", "eta vs Iso - MC matched ",61, -3.05,3.05, 200, -0.5, 1.5);  
    hEta_trkIso  = new TH2F("hEta_trkIso", "eta vs trkIso  ",61, -3.05,3.05 , 30, -0.25,14.75);  
    hEta_ecalIso = new TH2F("hEta_ecalIso","eta vs ecalIso ",61, -3.05,3.05 , 30, -0.25,14.75);  
    hEta_hcalIso = new TH2F("hEta_hcalIso","eta vs hcalIso ",61, -3.05,3.05 , 30, -0.25,14.75);  
 

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

    hPt_em1 = (TH2F *) file->Get("hPt_em1");
    hPt_hd1 = (TH2F *) file->Get("hPt_hd1");
    hPt_tk1 = (TH2F *) file->Get("hPt_tk1");
    hPt_sum1 = (TH2F *) file->Get("hPt_sum1");
    hPt_Cal1 = (TH2F *) file->Get("hPt_Cal1");
    hPt_Trk1 = (TH2F *) file->Get("hPt_Trk1");
    hPt_Iso1 = (TH2F *) file->Get("hPt_Iso1");
    hEta_Iso1 = (TH2F *) file->Get("hEta_Iso1");
    hEta_trkIso1  = (TH2F *) file->Get("hEta_trkIso1");
    hEta_ecalIso1 = (TH2F *) file->Get("hEta_ecalIso1");
    hEta_hcalIso1 = (TH2F *) file->Get("hEta_hcalIso1");

    hPt_em2 = (TH2F *) file->Get("hPt_em2");
    hPt_hd2 = (TH2F *) file->Get("hPt_hd2");
    hPt_tk2 = (TH2F *) file->Get("hPt_tk2");
    hPt_sum2 = (TH2F *) file->Get("hPt_sum2");
    hPt_Cal2 = (TH2F *) file->Get("hPt_Cal2");
    hPt_Trk2 = (TH2F *) file->Get("hPt_Trk2");
    hPt_Iso2 = (TH2F *) file->Get("hPt_Iso2");
    hEta_trkIso2  = (TH2F *) file->Get("hEta_trkIso2");
    hEta_ecalIso2 = (TH2F *) file->Get("hEta_ecalIso2");
    hEta_hcalIso2 = (TH2F *) file->Get("hEta_hcalIso2");

    hPt_EovP = (TH2F *) file->Get("hPt_EovP");
    hPt_HovE = (TH2F *) file->Get("hPt_HovE");
    hPt_nTrk = (TH2F *) file->Get("hPt_nTrk");
    hPt_nCal = (TH2F *) file->Get("hPt_nCal");
    hPt_Iso  = (TH2F *) file->Get("hPt_Iso");
    hEta_Iso = (TH2F *) file->Get("hEta_Iso");
    hEta_trkIso  = (TH2F *) file->Get("hEta_trkIso");
    hEta_ecalIso = (TH2F *) file->Get("hEta_ecalIso");
    hEta_hcalIso = (TH2F *) file->Get("hEta_hcalIso");

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
 
    delete hPt_em1;
    delete hPt_hd1;
    delete hPt_tk1;
    delete hPt_sum1;
    delete hPt_Cal1;
    delete hPt_Trk1;
    delete hPt_Iso1;
    delete hEta_Iso1;
    delete hEta_trkIso1;
    delete hEta_ecalIso1;
    delete hEta_hcalIso1;

    delete hPt_em2;
    delete hPt_hd2;
    delete hPt_tk2;
    delete hPt_sum2;
    delete hPt_Cal2;
    delete hPt_Trk2;
    delete hEta_trkIso2;
    delete hEta_ecalIso2;
    delete hEta_hcalIso2;
    delete hPt_Iso2;

    delete hPt_EovP;
    delete hPt_HovE;
    delete hPt_nTrk;
    delete hPt_nCal;
    delete hPt_Iso;
    delete hEta_Iso;
    delete hEta_trkIso;
    delete hEta_ecalIso;
    delete hEta_hcalIso;

 }

 void Fill4a(double Pt,double Eta ){
    hPt->Fill(Pt);
    hEta->Fill(Eta);
 }
 
 void Fill4b(double pt, double eta, double emE, double hdE, double tkP, double sumE, int nCal3, int ntrack3, 
             double trackIso, double ecalIso, double hcalIso, double Iso){

    hPt_em1->Fill(pt, emE);
    hPt_hd1->Fill(pt, hdE);
    hPt_tk1->Fill(pt, tkP);
    hPt_sum1->Fill(pt, sumE);
    hPt_Cal1->Fill(pt, nCal3);
    hPt_Trk1->Fill(pt, ntrack3);
    hEta_trkIso1->Fill(eta, trackIso);
    hEta_ecalIso1->Fill(eta, ecalIso);
    hEta_hcalIso1->Fill(eta, hcalIso);
    hPt_Iso1->Fill(pt, Iso);
    hEta_Iso1->Fill(eta, Iso);
 }
 
 void Fill4c(double pt, double eta, double emE, double hdE, double tkP, double sumE, int nCal3, int ntrack3, 
             double trackIso, double ecalIso, double hcalIso, double Iso){

    hPt_em2->Fill(pt, emE);
    hPt_hd2->Fill(pt, hdE);
    hPt_tk2->Fill(pt, tkP);
    hPt_sum2->Fill(pt, sumE);
    hPt_Cal2->Fill(pt, nCal3);
    hPt_Trk2->Fill(pt, ntrack3);
    hEta_trkIso2->Fill(eta, trackIso);
    hEta_ecalIso2->Fill(eta, ecalIso);
    hEta_hcalIso2->Fill(eta, hcalIso);
    hPt_Iso2->Fill(pt, Iso);
 }

 void Fill4d(double pt, double eta, double EovP, double HovE, int nTrk, int nCal, float trkIso, float ecalIso, float hcalIso, double Iso ) {
    hPt_EovP->Fill(pt, EovP);
    hPt_HovE->Fill(pt, HovE);
    hPt_nTrk->Fill(pt, nTrk);
    hPt_nCal->Fill(pt, nCal);
    hPt_Iso->Fill(pt, Iso);
    hEta_Iso->Fill(eta, Iso);
    hEta_trkIso->Fill(eta, trkIso);
    hEta_ecalIso->Fill(eta, ecalIso);
    hEta_hcalIso->Fill(eta, hcalIso);
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

 void Write() {

    hPt->Write();
    hEta->Write();
    hEta1a->Write();
    hEta1b->Write();
    hEta1c->Write();
    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
 
    hPt_em1->Write();
    hPt_hd1->Write();
    hPt_tk1->Write();
    hPt_sum1->Write();
    hPt_Cal1->Write();
    hPt_Trk1->Write();
    hPt_Iso1->Write();
    hEta_Iso1->Write();
    hEta_trkIso1->Write();
    hEta_ecalIso1->Write();
    hEta_hcalIso1->Write();

    hPt_em2->Write();
    hPt_hd2->Write();
    hPt_tk2->Write();
    hPt_sum2->Write();
    hPt_Cal2->Write();
    hPt_Trk2->Write();
    hPt_Iso2->Write();
    hEta_trkIso2->Write();
    hEta_ecalIso2->Write();
    hEta_hcalIso2->Write();

    hPt_EovP->Write();
    hPt_HovE->Write();
    hPt_nTrk->Write();
    hPt_nCal->Write();
    hPt_Iso->Write();
    hEta_Iso->Write();
    hEta_trkIso->Write();
    hEta_ecalIso->Write();
    hEta_hcalIso->Write();

 }

  TH1F *hPt;
  TH1F *hEta;
  TH1F *hEta1a;
  TH1F *hEta1b;
  TH1F *hEta1c;
  TH1F *hPt1a;
  TH2F *hPt1b;
  TH2F *hPt1c;

  TH2F *hPt_em1;
  TH2F *hPt_hd1;
  TH2F *hPt_tk1;
  TH2F *hPt_sum1;
  TH2F *hPt_Cal1;
  TH2F *hPt_Trk1;
  TH2F *hPt_Iso1;
  TH2F *hEta_Iso1;
  TH2F *hEta_trkIso1;
  TH2F *hEta_ecalIso1;
  TH2F *hEta_hcalIso1;

  TH2F *hPt_em2;
  TH2F *hPt_hd2;
  TH2F *hPt_tk2;
  TH2F *hPt_sum2;
  TH2F *hPt_Cal2;
  TH2F *hPt_Trk2;
  TH2F *hPt_Iso2;
  TH2F *hEta_trkIso2;
  TH2F *hEta_ecalIso2;
  TH2F *hEta_hcalIso2;

  TH2F *hPt_EovP;
  TH2F *hPt_HovE;
  TH2F *hPt_nTrk;
  TH2F *hPt_nCal;
  TH2F *hPt_Iso;
  TH2F *hEta_Iso;
  TH2F *hEta_trkIso;
  TH2F *hEta_ecalIso;
  TH2F *hEta_hcalIso;

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

    hNJets = new TH1F("hNJets"," # of selected W jets ",65,-0.5,64.5);
    hWp_mass = new TH2F("hWp_mass"," W mass vs W momentum from selected reco jets", 500,0.,500.,240,30.,150.);
    hdRWjj   = new TH1F("hdRWjj","dR for 2 W matched jets ",200, 0.,5. );  
    hRes_dR  = new TH1F("hRes_dR","Resolutio  of dR of 2 W matched jets ",400, -2., 2. );  

 } 

 HTOP6(TString name_, TFile* file) {
    name=name_;

    hNJets = (TH1F *) file->Get("hNJets");
    hWp_mass   = (TH2F *) file->Get("hWp_mass");
    hdRWjj  = (TH1F *) file->Get("hdRWjj");
    hRes_dR = (TH1F *) file->Get("hRes_dR");

 }

 /// Destructor
 virtual ~HTOP6() {

    delete hNJets;
    delete hWp_mass;
    delete hdRWjj;
    delete hRes_dR;

 }

 void Fill6b( double wmass, double wp, double dRWjj, double Res_dR ){
    hWp_mass->Fill(wp,wmass);
    hdRWjj->Fill(dRWjj);
    hRes_dR->Fill(Res_dR);
 }

 void Fill6c( int njets){
    hNJets->Fill(njets);
 }

 void Write() {

    hNJets->Write();
    hWp_mass->Write();
    hdRWjj->Write();
    hRes_dR->Write();

 }

  TH1F *hNJets;
  TH2F *hWp_mass;
  TH1F *hdRWjj;
  TH1F *hRes_dR;


 TString name;

};

class HTOP7 {
public:
 
 HTOP7(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hRes_Pt   = new TH1F("hRes_Pt","Resolutio  of b jet pt ",400, -2., 2. ); 
    hdRbMu    = new TH1F("hdRbMu","dR for b and  iso muon ",300, 0.,15. );  
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

    hPt_Eta_NoCalo = new TH2F("hPt_Eta_NoCalo", " not-caloJet, pt vs eta (b jets)", 500, 0., 250, 71, -3.55, 3.55);
    hRes_Pt_NoCalo = new TH1F("hRes_Pt_NoCalo", " not-caloJet, Pt(jet) - Pt(q) / Pt(q)   ",400,-2.,2.);

    gbJpt_h   = new TH2F("gbJpt_h"," gen b jet pt vs eta", 500,0,250, 71,-3.55,3.55);
    gEovH     = new TH1F("gEovH"," E/H for gen Jets", 500, -2., 23.);

    // parton level
    hEta1a  = new TH1F("hEta1a"," eta distribution for bjet of semi Tt(parton) ",71,-3.55,3.55);
    hEta1b  = new TH1F("hEta1b"," eta distribution for bjet of semi Tt(mcMatched)",71,-3.55,3.55);
    hEta1c  = new TH1F("hEta1c"," eta distribution for bjet of semi Tt(recoUsed)",71,-3.55,3.55);
    hPt1a   = new TH1F("hPt1a", " Pt distribution  for b of semi Tt(parton)",30,0,300 );
    hPt1b   = new TH2F("hPt1b", " Pt vs Pt_Resolution for b of semi Tt(mcMatched)",30,0,300, 31,-1.55,1.55 );
    hPt1c   = new TH2F("hPt1c", " Pt vs Pt_Resolution for b of semi Tt(recoUsed)",30,0,300, 31,-1.55,1.55 );
 } 

 HTOP7(TString name_, TFile* file) {
    name=name_;

    hRes_Pt = (TH1F *) file->Get("hRes_Pt");
    hdRbMu = (TH1F *) file->Get("hdRbMu");
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

    hPt_Eta_NoCalo = (TH2F *) file->Get("hPt_Eta_NoCalo");
    hRes_Pt_NoCalo = (TH1F *) file->Get("hRes_Pt_NoCalo");

    gbJpt_h = (TH2F *) file->Get("gbJpt_h");
    gEovH   = (TH1F *) file->Get("gEovH");

    hEta1a = (TH1F *) file->Get("hEta1a");
    hEta1b = (TH1F *) file->Get("hEta1b");
    hEta1c = (TH1F *) file->Get("hEta1c");
    hPt1a  = (TH1F *) file->Get("hPt1a");
    hPt1b  = (TH2F *) file->Get("hPt1b");
    hPt1c  = (TH2F *) file->Get("hPt1c");
 }

 /// Destructor
 virtual ~HTOP7() {

    delete hRes_Pt;
    delete hdRbMu;
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

    delete hPt_Eta_NoCalo;
    delete hRes_Pt_NoCalo;

    delete gbJpt_h;
    delete gEovH;

    delete hEta1a;
    delete hEta1b;
    delete hEta1c;
    delete hPt1a;
    delete hPt1b;
    delete hPt1c;
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
 void Fill7b( double dR_bmu ){
    hdRbMu->Fill(dR_bmu);
 }
 void Fill7c( double dR_bwj ){
    hdRbWj->Fill(dR_bwj);
 }
 void Fill7d( double bDisJProb, double bDisTkCount ){
    hbDis_all->Fill(bDisJProb,bDisTkCount );
 }

 void Fill7e( double pt, double eta, double Res_Pt) { 
    hPt_Eta_NoCalo->Fill(pt,eta);
    hRes_Pt_NoCalo->Fill(Res_Pt);
 }
 void Fill7f( double pt, double eta, double EovH) {
    gbJpt_h->Fill(pt, eta);
    gEovH->Fill(EovH);
 }

 void Fill7g(double eta, double pt ) {
    hEta1a->Fill(eta);
    hPt1a->Fill(pt);
 }
 void Fill7h(double eta, double pt, double ptRes) {
    hEta1b->Fill(eta);
    hPt1b->Fill(pt, ptRes);
 }
 void Fill7i(double eta, double pt, double ptRes) {
    hEta1c->Fill(eta);
    hPt1c->Fill(pt, ptRes);
 }

 void Write() {

    hRes_Pt->Write();
    hdRbMu->Write();
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

    hPt_Eta_NoCalo->Write();
    hRes_Pt_NoCalo->Write();

    gbJpt_h->Write();
    gEovH->Write();

    hEta1a->Write();
    hEta1b->Write();
    hEta1c->Write();
    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
  
 }

  TH1F *hRes_Pt;
  TH1F *hdRbMu;
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

  TH2F *hPt_Eta_NoCalo;
  TH1F *hRes_Pt_NoCalo;

  TH2F *gbJpt_h;
  TH1F *gEovH;

  TH1F *hEta1a;
  TH1F *hEta1b;
  TH1F *hEta1c;
  TH1F *hPt1a;
  TH2F *hPt1b;
  TH2F *hPt1c;

 TString name;

};

class HTOP8 {
public:
 
 HTOP8(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hEta_Pt1 = new TH2F("hEta_Pt1", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hEta_Pt2 = new TH2F("hEta_Pt2", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hNJets = new TH1F("hNJets"," # of selected W jets ",65,-0.5,64.5);
    hWp_mass = new TH2F("hWp_mass"," W mass vs W momentum from selected reco jets", 500,0.,500.,240,30.,150.);

    hNTk_mc  = new TH1F("hNTk_mc"," # of charge track of a jet", 50, -1.5, 48.5);
    hemE     = new TH1F("hemE"," em Energy fraction ", 150, -0.2, 1.3);
    hR60_mc  = new TH2F("hR60_mc"," n60 vs nConstituents of a W jet", 20, -1.5, 18.5, 50, -0.5, 49.5);
    hR90_mc  = new TH2F("hR90_mc"," n90 vs nConstituents of a W jet", 20, -1.5, 18.5, 50, -0.5, 49.5);
    hArea_pt_mc= new TH2F("hArea_pt_mc","tower area vs pt of a W jet", 500, 0., 500., 200, 0., 1.);
    hEovH_p_mc = new TH2F("hEovH_p_mc","E/H of a W jet", 500, 0., 500., 500, -2., 23.);
    hEovH_n_mc = new TH2F("hEovH_n_mc","E/H vs. nConstituent of a W jet", 50, -0.5, 49.5, 500, -2., 23.);
    hEovH_A_mc = new TH2F("hEovH_A_mc"," E/H vs tower area of a W jet", 200, 0., 1., 500, -2., 23.);
    hEovH_h_mc = new TH2F("hEovH_h_mc"," E/H vs eta of a W jet", 61, -3.05, 3.05, 500, -2., 23.);
    hPt_Eta_mc = new TH2F("hPt_Eta_mc", " eta vs pt (W jets)", 500, 0., 250, 71, -3.55, 3.55);

    hRes_Pt    = new TH1F("hRes_Pt", " Pt(jet) - Pt(q) / Pt(q)   ",400,-2.,2.);
    hRes_Eta   = new TH1F("hRes_Eta"," Eta(jet) - Eta(q) / Eta(q)",500,-5.,5.);

    hdRWjMu  = new TH1F("hdRWjMu","dR for isoMuon and W matched jets ",300, 0.,15. );  

    hPt_Eta_NoCalo = new TH2F("hPt_Eta_NoCalo", " not-caloJet, eta vs pt (W jets)", 500, 0., 500, 59, -2.95, 2.95);
    hRes_Pt_NoCalo = new TH1F("hRes_Pt_NoCalo", " not-caloJet, Pt(jet) - Pt(q) / Pt(q)   ",400,-2.,2.);

    hbDis_WCand  = new TH2F("hbDis_WCand","bCandidates jets JProb(x),TkCount(y)", 150,-0.1,1.4, 500,-20.,80.);

    gwJpt_h   = new TH2F("gwJpt_h"," gen w jet pt vs eta", 500,0,250, 71,-3.55,3.55);
    gEovH    = new TH1F("gEovH"," E/H for gen Jets", 500, -2., 23.);

    // parton level
    hEta1a  = new TH1F("hEta1a"," eta distribution for Wjet of semi Tt(parton) ",71,-3.55,3.55);
    hEta1b  = new TH1F("hEta1b"," eta distribution for Wjet of semi Tt(mcMatched)",71,-3.55,3.55);
    hEta1c  = new TH1F("hEta1c"," eta distribution for Wjet of semi Tt(recoUsed)",71,-3.55,3.55);
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

    hRes_Pt     = (TH1F *) file->Get("hRes_Pt");
    hRes_Eta    = (TH1F *) file->Get("hRes_Eta");

    hdRWjMu  = (TH1F *) file->Get("hdRWjMu");

    hPt_Eta_NoCalo = (TH2F *) file->Get("hPt_Eta_NoCalo");
    hRes_Pt_NoCalo = (TH1F *) file->Get("hRes_Pt_NoCalo");

    hbDis_WCand =(TH2F *) file->Get("hbDis_WCand");

    gwJpt_h = (TH2F *) file->Get("gwJpt_h");
    gEovH   = (TH1F *) file->Get("gEovH");

    hEta1a = (TH1F *) file->Get("hEta1a");
    hEta1b = (TH1F *) file->Get("hEta1b");
    hEta1c = (TH1F *) file->Get("hEta1c");
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

    delete hRes_Pt;
    delete hRes_Eta;
    delete hdRWjMu;

    delete hPt_Eta_NoCalo;
    delete hRes_Pt_NoCalo;

    delete hbDis_WCand;

    delete gwJpt_h;
    delete gEovH;

    delete hEta1a;
    delete hEta1b;
    delete hEta1c;
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

 void Fill8c( double pt, double eta, double EovH, int nCont, double n60, double n90, double area, int nTk, double emE, double Res_Pt, double Res_Eta , double bDisJProb, double bDisTkCount ){

    hNTk_mc->Fill(nTk);
    hemE->Fill(emE);
    hR60_mc->Fill(n60, nCont);
    hR90_mc->Fill(n90, nCont);
    hArea_pt_mc->Fill(pt,area);
    hEovH_p_mc->Fill(pt,EovH);
    hEovH_n_mc->Fill(nCont,EovH);
    hEovH_A_mc->Fill(area,EovH);
    hEovH_h_mc->Fill(eta,EovH);
    hPt_Eta_mc->Fill(pt,eta);
    hRes_Pt->Fill(Res_Pt);
    hRes_Eta->Fill(Res_Eta);

    hbDis_WCand->Fill(bDisJProb,bDisTkCount );
 }

 void Fill8d( double dR_WjMu ){
    hdRWjMu->Fill(dR_WjMu);
 }

 void Fill8e( double pt, double eta, double Res_Pt) { 
    hPt_Eta_NoCalo->Fill(pt,eta);
    hRes_Pt_NoCalo->Fill(Res_Pt);
 }
 void Fill8f( double pt, double eta, double EovH) {
    gwJpt_h->Fill(pt, eta);
    gEovH->Fill(EovH);
 }

 void Fill8g(double eta, double pt ) {
    hEta1a->Fill(eta);
    hPt1a->Fill(pt );
 }
 void Fill8h(double eta, double pt, double ptRes) {
    hEta1b->Fill(eta);
    hPt1b->Fill(pt, ptRes);
 }
 void Fill8i(double eta, double pt, double ptRes) {
    hEta1c->Fill(eta);
    hPt1c->Fill(pt, ptRes);
 }
 
 void Write() {
    hEta_Pt1->Write();
    hEta_Pt2->Write();
    hNJets->Write();
    hWp_mass->Write();

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

    hRes_Pt->Write();
    hRes_Eta->Write();
    hdRWjMu->Write();

    hPt_Eta_NoCalo->Write();
    hRes_Pt_NoCalo->Write();

    hbDis_WCand->Write();

    gwJpt_h->Write();
    gEovH->Write();

    hEta1a->Write();
    hEta1b->Write();
    hEta1c->Write();
    hPt1a->Write();
    hPt1b->Write();
    hPt1c->Write();
 }

  TH2F *hEta_Pt1;
  TH2F *hEta_Pt2;
  TH1F *hNJets;
  TH2F *hWp_mass;

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

  TH1F *hRes_Pt;
  TH1F *hRes_Eta;
  TH1F *hdRWjMu;

  TH2F *hPt_Eta_NoCalo;
  TH1F *hRes_Pt_NoCalo;

  TH2F *hbDis_WCand;
  TH2F *gwJpt_h;
  TH1F *gEovH;

  TH1F *hEta1a;
  TH1F *hEta1b;
  TH1F *hEta1c;
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
    hMCtt     = new TH1F("hMCtt","MCMatched top mass ",200,0,400);
    hMCtt1    = new TH1F("hMCtt1","MCMatched top mass passing evt selection",200,0,400);
    hRecott   = new TH1F("hRecott","Reco top mass ",200,0,400);
    hMClepW   = new TH1F("hMClepW",  " mt of leptonic W ",320,0,160.);
    hRecolepW = new TH1F("hRecolepW"," mt of leptonic W ",320,0,160.);
    hMChadW   = new TH1F("hMChadW","MC Hadronic W",320,0,160.);
    hRecohadW = new TH1F("hRecohadW","Reco Hadronic W ",320,0,160.);
    hEvtEff   = new TH1F("hEvtEff"," event pre-selection Eff", 10,-1.25,3.75 );
    hObjEff   = new TH1F("hObjEff"," object selection Eff", 10, -0.5,9.5 );
    hEvtShape1 = new TH2F("hEvtShape1"," event topo vs # of bjets ", 32,-1.5,30.5,  6,-0.5,5.5 );
    hEvtShape2 = new TH2F("hEvtShape2"," event topo vs # of wjets ", 32,-1.5,30.5, 10,-0.5,9.5 );
    hbJetEff   = new TH2F("hbJetEff"," # of correctly matched bjet ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hWJetEff   = new TH2F("hWJetEff"," # of correctly matched Wjet ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hWLepEff   = new TH2F("hWLepEff"," # of correctly matched lep. ", 4,-0.5,3.5, 120,-3.5,116.5 );
    hHLTBits   = new TH2F("hHLTBits","HLT Trigger bits(names) for hadronic tt", 166, -0.5, 165.5, 6,-0.5,5.5 );
    hHLTSelect = new TH1F("hHLTSelect","HLT Trigger Selected result "  , 15, -7.5, 7.5 );

    dR_lepW    = new TH1F("dR_lepW"," dR(reco - gen) ",200,0.,10.);
    PtRes_lepW = new TH1F("PtRes_lepW"," Pt Res ",200,-1.005,0.995);
    MRes_lepW  = new TH1F("MRes_lepW"," Mass Res ",200,-1.005,0.995);

 } 

 HTOP9(TString name_, TFile* file) {
    name=name_;

    hMCtt      = (TH1F *) file->Get("hMCtt");
    hMCtt1     = (TH1F *) file->Get("hMCtt1");
    hRecott    = (TH1F *) file->Get("hRecott");
    hMClepW    = (TH1F *) file->Get("hMClepW");
    hRecolepW  = (TH1F *) file->Get("hRecolepW");
    hMChadW    = (TH1F *) file->Get("hMChadW");
    hRecohadW  = (TH1F *) file->Get("hRecohadW");
    hEvtEff    = (TH1F *) file->Get("hEvtEff");
    hObjEff    = (TH1F *) file->Get("hObjEff");
    hEvtShape1 = (TH2F *) file->Get("hEvtShape1");
    hEvtShape2 = (TH2F *) file->Get("hEvtShape2");
    hbJetEff   = (TH2F *) file->Get("hbJetEff");
    hWJetEff   = (TH2F *) file->Get("hWJetEff");
    hWLepEff   = (TH2F *) file->Get("hWLepEff");
    hHLTBits   = (TH2F *) file->Get("hHLTBits");
    hHLTSelect = (TH1F *) file->Get("hHLTSelect");

    dR_lepW    = (TH1F *) file->Get("dR_lepW");
    PtRes_lepW = (TH1F *) file->Get("PtRes_lepW");
    MRes_lepW  = (TH1F *) file->Get("MRes_lepW");
 }

 /// Destructor
 virtual ~HTOP9() {

    delete hMCtt;
    delete hMCtt1;
    delete hRecott;
    delete hMClepW;
    delete hRecolepW;
    delete hMChadW;
    delete hRecohadW;
    delete hEvtEff;
    delete hObjEff;
    delete hEvtShape1;
    delete hEvtShape2;
    delete hbJetEff;
    delete hWJetEff;
    delete hWLepEff;
    delete hHLTBits;
    delete hHLTSelect;
   
    delete dR_lepW;
    delete PtRes_lepW;
    delete MRes_lepW;
 }

 void Fill9(double mt){
    hMCtt->Fill(mt);
 }
 void Fill9h(double mt){
    hMCtt1->Fill(mt);
 }
 void Fill9a(double mt){
    hRecott->Fill(mt);
 }
 void Fill9b(double mw) {
    hMClepW->Fill(mw) ; 
 }
 void Fill9c(double mtw) {
    hRecolepW->Fill(mtw) ; 
 }
 void Fill9d(double mw) {
    hMChadW->Fill(mw) ; 
 }
 void Fill9e(double mw) {
    hRecohadW->Fill(mw) ; 
 }
 void Fill9f(float count ) {
    hEvtEff->Fill( count );
 }
 void Fill9g(int id, size_t nBJ, size_t nWJ ) {
    hEvtShape1->Fill( id, nBJ );
    hEvtShape2->Fill( id, nWJ );
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
 void Fill9m( double dR, double PtRes, double MRes ) {
    dR_lepW->Fill(dR);
    PtRes_lepW->Fill(PtRes);
    MRes_lepW->Fill(MRes);
 }

 void Write() {

    hMCtt->Write();
    hMCtt1->Write();
    hRecott->Write();
    hMClepW->Write();
    hMChadW->Write();
    hRecolepW->Write();
    hRecohadW->Write();
    hEvtEff->Write();
    hObjEff->Write();
    hEvtShape1->Write();
    hEvtShape2->Write();
    hbJetEff->Write();
    hWJetEff->Write();
    hWLepEff->Write();
    hHLTBits->Write();
    hHLTSelect->Write();

    dR_lepW->Write();
    PtRes_lepW->Write();
    MRes_lepW->Write();

 }

  TH1F *hMCtt;
  TH1F *hMCtt1;
  TH1F *hRecott;
  TH1F *hMClepW;
  TH1F *hMChadW;
  TH1F *hRecolepW;
  TH1F *hRecohadW;
  TH1F *hEvtEff; 
  TH1F *hObjEff; 
  TH2F *hEvtShape1; 
  TH2F *hEvtShape2; 
  TH2F *hbJetEff;
  TH2F *hWJetEff;
  TH2F *hWLepEff;
  TH2F *hHLTBits;
  TH1F *hHLTSelect;

  TH1F *dR_lepW;
  TH1F *PtRes_lepW;
  TH1F *MRes_lepW;

 TString name;

};

 
#endif
