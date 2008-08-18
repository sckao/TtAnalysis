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
    hEovH  = new TH1F("hEovH"," em Energy / hadronic Energy ", 1050, -2., 103.);
    hNCont = new TH1F("hNCont"," # of constituents of a jet", 50, 0, 50);
    hNTk   = new TH1F("hNTk"," # of charge track of a jet", 50, 0, 50);
    hR60   = new TH1F("hR60", "n60/ total # of constituents", 80, -0.5, 1.5);
    hR90   = new TH1F("hR90", "n90/ total # of constituents", 80, -0.5, 1.5);
    hArea_pt= new TH2F("hArea_pt"," tower area vs pt ", 1000, 0., 1000., 200, 0., 1.);
    hEovH_pt = new TH2F("hEovH_pt"," em Energy / hadronic Energy vs pt", 1000, 0., 1000., 1050, -2., 103.);
    hEovH_r = new TH2F("hEovH_r"," E/H vs. nConstituent/nChargeTrk ", 80, -1.45, 6.55, 1050, -2., 103.);
    hEovH_A = new TH2F("hEovH_A"," E/H vs tower area ", 200, 0., 1., 1050, -2., 103.);
    hEovH_h = new TH2F("hEovH_h"," E/H vs eta ", 59, -2.95, 2.95, 1050, -2., 103.);
    hEta_Pt = new TH2F("hEta_Pt", " eta vs pt (all jets)", 59, -2.95, 2.95, 500, 0., 500);
    hEta_uncorrPt = new TH2F("hEta_uncorrPt", "eta vs uncorreted pt (all jets)", 59, -2.95, 2.95, 500, 0., 500);
    hemF    = new TH2F("hemF", " emF(x) vs emFCalo(y)", 150, -0.2, 1.3, 150, -0.2, 1.3);
    
    hPtCorr_UnCorr = new TH2F("hPtCorr_UnCOrr","pt: corr(x), uncorr(y)", 500, 0., 500, 500, 0., 500);
    hECorr_UnCorr  = new TH2F("hECorr_UnCOrr","E: corr(x), uncorr(y)", 500, 0., 500, 500, 0., 500);
    hPtCorr_UnCorr2 = new TH2F("hPtCorr_UnCOrr2","pt: corr(x), uncorr(y)", 500, 0., 500, 500, 0., 500);
    hECorr_UnCorr2  = new TH2F("hECorr_UnCOrr2","E: corr(x), uncorr(y)", 500, 0., 500, 500, 0., 500);
    hPtCorr_UnCorr3 = new TH2F("hPtCorr_UnCOrr3","pt: corr(x), uncorr(y)", 500, 0., 500, 500, 0., 500);
    hECorr_UnCorr3  = new TH2F("hECorr_UnCOrr3","E: corr(x), uncorr(y)", 500, 0., 500, 500, 0., 500);

    hWmass = new TH1F("hWmass"," W mass from gen-reco matching jets",240,30.,150.);
    hWp    = new TH1F("hWp", " momentum of W ",500,0.,500.);

    hNJets = new TH1F("hNJets"," # of Jets",65,-0.5,64.5);
    hNjet20= new TH1F("hNjet20"," # of Jets, pt > 20 GeV ",65,-0.5,64.5);
    hNjet30= new TH1F("hNjet30"," # of Jets, pt > 30 GeV ",65,-0.5,64.5);
    hNjet40= new TH1F("hNjet40"," # of Jets, pt > 40 GeV ",65,-0.5,64.5);
    
    hWm2   = new TH1F("hWm2"," W mass from q qbar ",240,30.,150.);
    hWp2   = new TH1F("hWp2", " momentum of W from q qbar ",500,0.,500.);

    hEovH_h1 = new TH2F("hEovH_h1"," high E/H vs eta ", 59, -2.95, 2.95, 1050, -2., 103.);
    hEovH_A1 = new TH2F("hEovH_A1"," high E/H vs towersArea", 200, 0., 1., 1050, -2., 103.);
    hEovH_r1 = new TH2F("hEovH_r1"," high E/H vs nConstituent/nChargeTrk ", 80, -1.45, 6.55, 1050, -2., 103.);
    hEovH_N1 = new TH2F("hEovH_N1"," high E/H vs nConstituent ", 50, 0., 50., 1050, -2., 103.);
    hEovH_C1 = new TH2F("hEovH_C1"," high E/H vs nCharge Trk ", 50, 0., 50., 1050, -2., 103.);
    hR60_1  = new TH1F("hR60_1", "n60/ total # of constituents", 80, -0.5, 1.5);

    hEta_2   = new TH1F("hEta_2"," eta distribution",59,-2.95,2.95);
    hEovH_r_2 = new TH2F("hEovH_r_2"," E/H vs. nConstituent/nChargeTrk ", 80, -1.45, 6.55, 1050, -2., 103.);
    hNCont_2 = new TH1F("hNCont_2"," # of constituents of crapy jets", 50, 0, 50);
    hNTk_2   = new TH1F("hNTk_2"," # of charge track of crapy jets", 50, 0, 50);
    hR60_2   = new TH1F("hR60_2", "n60/ total # of constituents of crapy jets", 80, -0.5, 1.5);
    hR90_2   = new TH1F("hR90_2", "n90/ total # of constituents of crapy jets", 80, -0.5, 1.5);

 } 

 HTOP1(TString name_, TFile* file) {
    name=name_;

    hEt    = (TH1F *) file->Get("hEt");
    hEovH  = (TH1F *) file->Get("hEovH");
    hNCont = (TH1F *) file->Get("hNCont");
    hNTk   = (TH1F *) file->Get("hNTK");
    hR60   = (TH1F *) file->Get("hR60");
    hR90   = (TH1F *) file->Get("hR90");
    hArea_pt = (TH2F *) file->Get("hArea_pt");
    hEovH_pt  = (TH2F *) file->Get("hEovH_pt");
    hEovH_r  = (TH2F *) file->Get("hEovH_r");
    hEovH_A  = (TH2F *) file->Get("hEovH_A");
    hEovH_h  = (TH2F *) file->Get("hEovH_h");
    hEta_Pt  = (TH2F *) file->Get("hEta_Pt");
    hEta_uncorrPt  = (TH2F *) file->Get("hEta_uncorrPt");
    hemF     = (TH2F *) file->Get("hemF");

    hPtCorr_UnCorr = (TH2F *) file->Get("hPtCorr_UnCorr");
    hECorr_UnCorr  = (TH2F *) file->Get("hECorr_UnCorr");
    hPtCorr_UnCorr2 = (TH2F *) file->Get("hPtCorr_UnCorr2");
    hECorr_UnCorr2  = (TH2F *) file->Get("hECorr_UnCorr2");
    hPtCorr_UnCorr3 = (TH2F *) file->Get("hPtCorr_UnCorr3");
    hECorr_UnCorr3  = (TH2F *) file->Get("hECorr_UnCorr3");

    hWmass = (TH1F *) file->Get("hWmass");
    hWp    = (TH1F *) file->Get("hWp");

    hNJets = (TH1F *) file->Get("hNJets");
    hNjet20= (TH1F *) file->Get("hNjet20");
    hNjet30= (TH1F *) file->Get("hNjet30");
    hNjet40= (TH1F *) file->Get("hNjet40");

    hWm2   = (TH1F *) file->Get("hWm2");
    hWp2   = (TH1F *) file->Get("hWp2");

    hEovH_A1  = (TH2F *) file->Get("hEovH_A1");
    hEovH_h1  = (TH2F *) file->Get("hEovH_h1");
    hEovH_N1  = (TH2F *) file->Get("hEovH_N1");
    hEovH_C1  = (TH2F *) file->Get("hEovH_C1");
    hEovH_r1  = (TH2F *) file->Get("hEovH_r1");
    hR60_1    = (TH1F *) file->Get("hR60_1");

    hEta_2   = (TH1F *) file->Get("hEta_2");
    hEovH_r_2= (TH2F *) file->Get("hEovH_r_2");
    hNCont_2 = (TH1F *) file->Get("hNCont_2");
    hNTk_2   = (TH1F *) file->Get("hNTK_2");
    hR60_2   = (TH1F *) file->Get("hR60_2");
    hR90_2   = (TH1F *) file->Get("hR90_2");

 }

 /// Destructor
 virtual ~HTOP1() {
    //delete hcsc_q;
    delete hEt;
    delete hEovH;
    delete hNCont;
    delete hNTk;
    delete hR60;
    delete hR90;
    delete hArea_pt;
    delete hEovH_pt;
    delete hEovH_r;
    delete hEovH_A;
    delete hEovH_h;
    delete hEta_Pt;
    delete hEta_uncorrPt;
    delete hemF;

    delete hPtCorr_UnCorr; 
    delete hECorr_UnCorr; 
    delete hPtCorr_UnCorr2; 
    delete hECorr_UnCorr2; 
    delete hPtCorr_UnCorr3; 
    delete hECorr_UnCorr3; 

    delete hWmass;
    delete hWp;

    delete hNJets;
    delete hNjet20;
    delete hNjet30;
    delete hNjet40;

    delete hWm2;
    delete hWp2;

    delete hEovH_A1;
    delete hEovH_h1;
    delete hEovH_N1;
    delete hEovH_C1;
    delete hEovH_r1;
    delete hR60_1;

    delete hEta_2;
    delete hEovH_r_2;
    delete hNCont_2;
    delete hNTk_2;
    delete hR60_2;
    delete hR90_2;

 }

 void Fill1a(double eta, double pt, double et, double EovH, int nCont, double r60, double r90, double uncorrPt, double r , double area, int nTk, double emF, double emFCalo ){
    hEt->Fill(et);
    hEovH->Fill(EovH);
    hNCont->Fill(nCont);
    hNTk->Fill(nTk);
    hR60->Fill(r60);
    hR90->Fill(r90);
    hArea_pt->Fill(pt,area);
    hEovH_pt->Fill(pt,EovH);
    hEovH_r->Fill(r,EovH);
    hEovH_A->Fill(area,EovH);
    hEovH_h->Fill(eta,EovH);
    hEta_Pt->Fill(eta, pt);
    hEta_uncorrPt->Fill(eta, uncorrPt);
    hemF->Fill(emF, emFCalo);
 } 
 void Fill1b( double wmass, double wp ){
    hWmass->Fill(wmass);
    hWp->Fill(wp);
 }
 void Fill1c(double Pt, double uncorrPt, double energy, double caloE ) {
    hPtCorr_UnCorr->Fill(Pt,uncorrPt);   
    hECorr_UnCorr->Fill(energy,caloE);   
 } 
 void Fill1c2(double Pt, double uncorrPt, double energy, double caloE ) {
    hPtCorr_UnCorr2->Fill(Pt,uncorrPt);   
    hECorr_UnCorr2->Fill(energy,caloE);   
 } 
 void Fill1c3(double Pt, double uncorrPt, double energy, double caloE ) {
    hPtCorr_UnCorr3->Fill(Pt,uncorrPt);   
    hECorr_UnCorr3->Fill(energy,caloE);   
 } 
 void Fill1f( int njets, int nj20, int nj30, int nj40 ){
    hNJets->Fill(njets);  
    hNjet20->Fill(nj20);  
    hNjet30->Fill(nj30);  
    hNjet40->Fill(nj40);  
 } 
 void Fill1g( double wmass, double wp ){
    hWm2->Fill(wmass);
    hWp2->Fill(wp);
 }
 void Fill1h( double EovH, double eta, double area, int nCon, int nTrk, double r, double r60 ) {
    hEovH_A1->Fill(area,EovH);
    hEovH_h1->Fill(eta,EovH);
    hEovH_N1->Fill(nCon,EovH);
    hEovH_C1->Fill(nTrk,EovH);
    hEovH_r1->Fill(r,EovH);
    hR60_1->Fill(r60);
 }
 void Fill1i( double eta, double EovH, int nCont, int nTk, double r60, double r90, double r){
    hEta_2->Fill(eta);
    hEovH_r_2->Fill(r,EovH);
    hNCont_2->Fill(nCont);
    hNTk_2->Fill(nTk);
    hR60_2->Fill(r60);
    hR90_2->Fill(r90);
 }

 void Write() {
    hEt->Write();
    hEovH->Write();
    hNCont->Write();
    hNTk->Write();
    hR60->Write();
    hR90->Write();
    hArea_pt->Write();
    hEovH_pt->Write();
    hEovH_r->Write();
    hEovH_A->Write();
    hEovH_h->Write();
    hEta_Pt->Write();
    hEta_uncorrPt->Write();
    hemF->Write();

    hPtCorr_UnCorr->Write();
    hECorr_UnCorr->Write();
    hPtCorr_UnCorr2->Write();
    hECorr_UnCorr2->Write();
    hPtCorr_UnCorr3->Write();
    hECorr_UnCorr3->Write();

    hWmass->Write();
    hWp->Write();

    hNJets->Write();
    hNjet20->Write();
    hNjet30->Write();
    hNjet40->Write();

    hWm2->Write();
    hWp2->Write();

    hEovH_A1->Write();
    hEovH_h1->Write();
    hEovH_N1->Write();
    hEovH_C1->Write();
    hEovH_r1->Write();
    hR60_1->Write();

    hEta_2->Write();
    hEovH_r_2->Write();
    hNCont_2->Write();
    hNTk_2->Write();
    hR60_2->Write();
    hR90_2->Write();

 }

  TH1F *hEt;
  TH1F *hEovH;
  TH1F *hNCont;
  TH1F *hNTk;
  TH1F *hR60;
  TH1F *hR90;
  TH2F *hArea_pt;
  TH2F *hEovH_pt;
  TH2F *hEovH_r;
  TH2F *hEovH_A;
  TH2F *hEovH_h;
  TH2F *hEta_Pt;
  TH2F *hEta_uncorrPt;
  TH2F *hemF;

  TH2F *hPtCorr_UnCorr;
  TH2F *hECorr_UnCorr;
  TH2F *hPtCorr_UnCorr2;
  TH2F *hECorr_UnCorr2;
  TH2F *hPtCorr_UnCorr3;
  TH2F *hECorr_UnCorr3;

  TH1F *hWmass;
  TH1F *hWp;

  TH1F *hNJets;
  TH1F *hNjet20;
  TH1F *hNjet30;
  TH1F *hNjet40;

  TH1F *hWm2;
  TH1F *hWp2;

  TH2F *hEovH_A1;
  TH2F *hEovH_h1;
  TH2F *hEovH_N1;
  TH2F *hEovH_C1;
  TH2F *hEovH_r1;
  TH1F *hR60_1;

  TH1F *hEta_2;
  TH2F *hEovH_r_2;
  TH1F *hNCont_2;
  TH1F *hNTk_2;
  TH1F *hR60_2;
  TH1F *hR90_2;


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
    hEta  = new TH1F("hEta"," eta distribution",59,-2.95,2.95);

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


 void Write() {

    hPt->Write();
    hEta->Write();
 
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

 HTOP4(TString name_, TFile* file) {
    name=name_;

    hPt    = (TH1F *) file->Get("hPt");
    hEta   = (TH1F *) file->Get("hEta");

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
 virtual ~HTOP4() {

    delete hPt;
    delete hEta;
 
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

 void Fill4a(double Pt,double Eta ){
    hPt->Fill(Pt);
    hEta->Fill(Eta);
 }
 
 void Fill4b(double emEt3, double emEt5, double hdEt3, double hdEt5, int nCal3, int nCal5, 
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
 
 void Fill4c(double emEt3, double emEt5, double hdEt3, double hdEt5, int nCal3, int nCal5,
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

 void Fill4d(double isoMuE, double mu_p, double jet_E ) {
    hIMu_caloE->Fill(isoMuE, mu_p);
    hIMu_caloE_jE->Fill(isoMuE, jet_E);
 }
 void Fill4e(double jetMuE, double mu_p, double jet_E ) {
    hJMu_caloE->Fill(jetMuE, mu_p);
    hJMu_caloE_jE->Fill(jetMuE, jet_E);
 }


 void Write() {

    hPt->Write();
    hEta->Write();
 
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

class HTOP5 {
public:
 
 HTOP5(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hEta_Pt1 = new TH2F("hEta_Pt1", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );
    hEta_Pt2 = new TH2F("hEta_Pt2", " eta vs pt (selected W jets)", 59, -2.95, 2.95, 500, 0., 500 );

 } 

 HTOP5(TString name_, TFile* file) {
    name=name_;

    hEta_Pt1 = (TH2F *) file->Get("hEta_Pt1");
    hEta_Pt2 = (TH2F *) file->Get("hEta_Pt2");

 }

 /// Destructor
 virtual ~HTOP5() {
    delete hEta_Pt1;
    delete hEta_Pt2;

 }

 void Fill5a( double pt1, double pt2, double eta1, double eta2){
    hEta_Pt1->Fill(eta1, pt1);
    hEta_Pt2->Fill(eta2, pt2);
 }


 void Write() {
    hEta_Pt1->Write();
    hEta_Pt2->Write();
 }

  TH2F *hEta_Pt1;
  TH2F *hEta_Pt2;

 TString name;

};

class HTOP6 {
public:
 
 HTOP6(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

    hPt_j   = new TH1F("hPt_j", " Pt distribution of matched jets ",500,0,500);
    hEta_j  = new TH1F("hEta_j"," eta distribution of matched jets  ",59,-2.95,2.95);
    hEovH_j = new TH1F("hEovH_j"," em Energy / hadronic Energy ", 1050, -2., 103.);
    hNCont_j= new TH1F("hNCont_j"," # of constituents of a W jet", 50, 0, 50);
    hNTk_j  = new TH1F("hNTk_j"," # of charge track of a jet", 50, 0, 50);
    hR60_j  = new TH1F("hR60_j", "n60/ total # of constituents of a W jet", 80, -0.5, 1.5);
    hR90_j  = new TH1F("hR90_j", "n90/ total # of constituents of a W jet", 80, -0.5, 1.5);
    hArea_pt_j= new TH2F("hArea_pt_j","tower area vs pt of a W jet", 1000, 0., 1000., 200, 0., 1.);
    hEovH_p_j = new TH2F("hEovH_p_j","E/H of a W jet", 1000, 0., 1000., 1050, -2., 103.);
    hEovH_r_j = new TH2F("hEovH_r_j","E/H vs. nConstituent/nChargeTrk of a W jet", 100, -1.45, 6.55, 1050, -2., 103.);
    hEovH_A_j = new TH2F("hEovH_A_j"," E/H vs tower area of a W jet", 200, 0., 1., 1050, -2., 103.);
    hEovH_h_j = new TH2F("hEovH_h_j"," E/H vs eta of a W jet", 59, -2.95, 2.95, 1050, -2., 103.);
    hemE      = new TH1F("hemE"," em Energy fraction ", 150, -0.2, 1.3);
    hEta_Pt_j = new TH2F("hEta_Pt_j", " eta vs pt (W jets)", 500, 0., 500, 59, -2.95, 2.95);

    hNJets = new TH1F("hNJets"," # of selected W jets ",65,-0.5,64.5);
    hWp_mass = new TH2F("hWp_mass"," W mass vs W momentum from selected reco jets", 500,0.,500.,240,30.,150.);
    hdRWjj   = new TH1F("hdRWjj","dR for 2 W matched jets ",200, 0.,5. );  
    hRes_dR  = new TH1F("hRes_dR","Resolutio  of dR of 2 W matched jets ",400, -2., 2. );  

    hPt_q  = new TH1F("hPt_q", " Pt distribution of q from W ",500,0,500);
    hEta_q = new TH1F("hEta_q"," eta distribution of q from W ",59,-2.95,2.95);
    hjdPt    = new TH1F("hjdPt", " Pt(jet) - Pt(q) ",500,-5.,5.);
    hdRWjMu  = new TH1F("hdRWjMu","dR for isoMuon and W matched jets ",300, 0.,15. );  

 } 

 HTOP6(TString name_, TFile* file) {
    name=name_;

    hPt_j   = (TH1F *) file->Get("hPt_j");
    hEta_j  = (TH1F *) file->Get("hEta_j");
    hEovH_j = (TH1F *) file->Get("hEovH_j");
    hNCont_j = (TH1F *) file->Get("hNCont_j");
    hNTk_j   = (TH1F *) file->Get("hNTK_j");
    hR60_j = (TH1F *) file->Get("hR60_j");
    hR90_j = (TH1F *) file->Get("hR90_j");
    hArea_pt_j = (TH2F *) file->Get("hArea_pt_j");
    hEovH_p_j  = (TH2F *) file->Get("hEovH_p_j");
    hEovH_r_j  = (TH2F *) file->Get("hEovH_r_j");
    hEovH_A_j  = (TH2F *) file->Get("hEovH_A_j");
    hEovH_h_j  = (TH2F *) file->Get("hEovH_h_j");
    hemE       = (TH1F *) file->Get("hemE");
    hEta_Pt_j  = (TH2F *) file->Get("hEta_Pt_j");

    hNJets = (TH1F *) file->Get("hNJets");
    hWp_mass   = (TH2F *) file->Get("hWp_mass");
    hdRWjj  = (TH1F *) file->Get("hdRWjj");
    hRes_dR = (TH1F *) file->Get("hRes_dR");

    hPt_q  = (TH1F *) file->Get("hPt_q");
    hEta_q = (TH1F *) file->Get("hEta_q");
    hjdPt    = (TH1F *) file->Get("hjdPt");
    hdRWjMu  = (TH1F *) file->Get("hdRWjMu");

 }

 /// Destructor
 virtual ~HTOP6() {

    delete hPt_j;
    delete hEta_j;
    delete hEovH_j;
    delete hNCont_j;
    delete hNTk_j;
    delete hR60_j;
    delete hR90_j;
    delete hArea_pt_j;
    delete hEovH_p_j;
    delete hEovH_r_j;
    delete hEovH_A_j;
    delete hEovH_h_j;
    delete hemE;
    delete hEta_Pt_j;

    delete hNJets;
    delete hWp_mass;
    delete hdRWjj;
    delete hRes_dR;

    delete hPt_q;
    delete hEta_q;
    delete hjdPt;
    delete hdRWjMu;

 }

 void Fill6a( double pt, double eta_j, double EovH, int nCont, double r60, double r90, double r, double area, int nTk, double emE){
    hPt_j->Fill(pt);
    hEta_j->Fill(eta_j);
    hEovH_j->Fill(EovH);
    hNCont_j->Fill(nCont);
    hNTk_j->Fill(nTk);
    hR60_j->Fill(r60);
    hR90_j->Fill(r90);
    hArea_pt_j->Fill(pt,area);
    hEovH_p_j->Fill(pt,EovH);
    hEovH_r_j->Fill(r,EovH);
    hEovH_A_j->Fill(area,EovH);
    hEovH_h_j->Fill(eta_j,EovH);
    hemE->Fill(emE);
    hEta_Pt_j->Fill(pt,eta_j);
 }

 void Fill6b( double wmass, double wp, double dRWjj, double Res_dR ){
    hWp_mass->Fill(wp,wmass);
    hdRWjj->Fill(dRWjj);
    hRes_dR->Fill(Res_dR);
 }

 void Fill6c( int njets){
    hNJets->Fill(njets);
 }

 void Fill6d( double pt_q, double eta_q, double dPt, double dR_WjMu ){
    hPt_q->Fill(pt_q);
    hEta_q->Fill(eta_q);
    hjdPt->Fill(dPt);
    hdRWjMu->Fill(dR_WjMu);
 }

 void Write() {
    hPt_j->Write();
    hEta_j->Write();
    hEovH_j->Write();
    hNCont_j->Write();
    hNTk_j->Write();
    hR60_j->Write();
    hR90_j->Write();
    hArea_pt_j->Write();
    hEovH_p_j->Write();
    hEovH_r_j->Write();
    hEovH_A_j->Write();
    hEovH_h_j->Write();
    hemE->Write();
    hEta_Pt_j->Write();
 
    hNJets->Write();
    hWp_mass->Write();
    hdRWjj->Write();
    hRes_dR->Write();

    hPt_q->Write();
    hEta_q->Write();
    hjdPt->Write();
    hdRWjMu->Write();

 }

  TH1F *hPt_j;
  TH1F *hEta_j;
  TH1F *hEovH_j;
  TH1F *hNCont_j;
  TH1F *hNTk_j;
  TH1F *hR60_j;
  TH1F *hR90_j;
  TH2F *hArea_pt_j;
  TH2F *hEovH_p_j;
  TH2F *hEovH_r_j;
  TH2F *hEovH_A_j;
  TH2F *hEovH_h_j;
  TH1F *hemE;
  TH2F *hEta_Pt_j;

  TH1F *hNJets;
  TH2F *hWp_mass;
  TH1F *hdRWjj;
  TH1F *hRes_dR;

  TH1F *hPt_q;
  TH1F *hEta_q;
  TH1F *hjdPt;
  TH1F *hdRWjMu;

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

 } 

 HTOP7(TString name_, TFile* file) {
    name=name_;

    hRes_Pt = (TH1F *) file->Get("hRes_Pt");
    hdRbMu = (TH1F *) file->Get("hdRbMu");
    hdRbWj = (TH1F *) file->Get("hdRbWj");
    hbDis_bCand =(TH2F *) file->Get("hbDis_bCand");
    hbDis_all   =(TH2F *) file->Get("hbDis_all");

 }

 /// Destructor
 virtual ~HTOP7() {

    delete hRes_Pt;
    delete hdRbMu;
    delete hdRbWj;
    delete hbDis_bCand;
    delete hbDis_all;

 }

 void Fill7a( double Res_pt, double bDisJProb, double bDisTkCount ){
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


 void Write() {

    hRes_Pt->Write();
    hdRbMu->Write();
    hdRbWj->Write();
    hbDis_bCand->Write();
    hbDis_all->Write();

 }

  TH1F *hRes_Pt;
  TH1F *hdRbMu;
  TH1F *hdRbWj;
  TH2F *hbDis_bCand;
  TH2F *hbDis_all;

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

 } 

 HTOP8(TString name_, TFile* file) {
    name=name_;

    hEta_Pt1 = (TH2F *) file->Get("hEta_Pt1");
    hEta_Pt2 = (TH2F *) file->Get("hEta_Pt2");
    hNJets = (TH1F *) file->Get("hNJets");
    hWp_mass   = (TH2F *) file->Get("hWp_mass");

 }

 /// Destructor
 virtual ~HTOP8() {
    delete hEta_Pt1;
    delete hEta_Pt2;
    delete hNJets;
    delete hWp_mass;

 }

 void Fill8a( double pt1, double pt2, double eta1, double eta2, int njets){
    hEta_Pt1->Fill(eta1, pt1);
    hEta_Pt2->Fill(eta2, pt2);
    hNJets->Fill(njets);
 }

 void Fill8b( double wmass, double wp ){
    hWp_mass->Fill(wp,wmass);
 }

 void Write() {
    hEta_Pt1->Write();
    hEta_Pt2->Write();
    hNJets->Write();
    hWp_mass->Write();
 }

  TH2F *hEta_Pt1;
  TH2F *hEta_Pt2;
  TH1F *hNJets;
  TH2F *hWp_mass;

 TString name;

};


class HTOP9 {
public:
 
 HTOP9(std::string name_) {
    TString N1 = name_.c_str();
    name=N1;

   // reconstructed objects masses
    hMCtt    = new TH1F("hMCtt","",200,0,400);
    hRecott  = new TH1F("hRecott","",200,0,400);

 } 

 HTOP9(TString name_, TFile* file) {
    name=name_;

    hMCtt      = (TH1F *) file->Get("hMCtt");
    hRecott    = (TH1F *) file->Get("hRecott");
 }

 /// Destructor
 virtual ~HTOP9() {

    delete hMCtt;
    delete hRecott;
 }

 void Fill9(double mt){
    hMCtt->Fill(mt);
 }
 void Fill9a(double mt){
    hRecott->Fill(mt);
 }


 void Write() {

    hMCtt->Write();
    hRecott->Write();

 }

  TH1F *hMCtt;
  TH1F *hRecott;

 TString name;

};

 
#endif
