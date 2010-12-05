#ifndef TtAnalysisHisto_H
#define TtAnalysisHisto_H

/** \class TopAnalysisHistograms
 *  Collection of histograms for TopAnalysis test.
 *
 * Author: S.C. Kao  - UC Riverside
 */

#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
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

    hNJets = new TH1F("hNJets", " # of Jets",65,-0.5,64.5);

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

    delete hNJets;

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
 void Fill1f( int njets ){
    hNJets->Fill(njets);  
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

    hNJets->Write();
 
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

  TH1F *hNJets;

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
    hResTC  = new TH1F("ResTC", " TC MET Resolution", 100, -5.,5.);
    hPhiPAT = new TH1F("PhiPAT","PAT Phi Resolution", 100, -5.,5.);
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

 /// Destructor
 virtual ~HTOP2() {

    delete hemfMET;
    delete hSumEt_MET;
    delete hResPAT;
    delete hResTC;
    delete hPhiPAT;
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

 void Fill2b(double respat, double restc, double phipat, double phitc ){
    hResPAT->Fill(respat);
    hResTC->Fill(restc);
    hPhiPAT->Fill(phipat);
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
    hResTC->Write();
    hPhiPAT->Write();
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
  TH1F *hResTC;
  TH1F *hPhiPAT;
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

    dR_neu     = new TH1F("dR_neu"," dR(reco - gen) for neutrino ",200,0.,10.);
    PtRes_lepW = new TH1F("PtRes_lepW"," Pt Res for matched leptonic W ",200,-1.005,0.995);
    PzRes_lepW = new TH1F("PzRes_lepW"," Pz Res for matched leptonic W ",200,-1.005,0.995);
    NeuMET_Pz  = new TH2F("NeuMET_Pz"  ," Neutrio Pz vs MET Pz from W solution", 400,-200.,200., 400, -200.,200.);
    METPzSol   = new TH2F("METPzSol"  ," 2 MET Pz from W solution", 400,-200.,200., 400, -200.,200.);
 
 } 

 /// Destructor
 virtual ~HTOP6() {

    delete dR_neu;
    delete PtRes_lepW;
    delete PzRes_lepW;
    delete NeuMET_Pz;
    delete METPzSol;

 }

 void Fill6d( double dR, double PtRes, double MRes, double neuPz, double metPz ) {
    dR_neu->Fill(dR);
    PtRes_lepW->Fill(PtRes);
    PzRes_lepW->Fill(MRes);
    NeuMET_Pz->Fill( neuPz, metPz);
 }
 void Fill6e( double pz1, double pz2) {
    METPzSol->Fill(pz1, pz2);
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    dR_neu->Write();
    PtRes_lepW->Write();
    PzRes_lepW->Write();
    NeuMET_Pz->Write();  
    METPzSol->Write();  

 }

  TH1F *dR_neu;
  TH1F *PtRes_lepW;
  TH1F *PzRes_lepW;
  TH2F *NeuMET_Pz;
  TH2F *METPzSol;

};

class HTOP7 {
public:
 
 HTOP7() {

    bDis_phi_d0  = new TH2F("bDis_phi_d0" ," phi vs d0  ",63, -3.15,3.15, 200,-0.5,0.5);
    bDis_dR_RelPt = new TH2F("bDis_dR_RelPt"," dR vs RelPt ",100, 0.,1., 300, 0.,60.);

    // parton level
    hEta     = new TH2F("hEta"," MC(X) vs reco(Y) , eta distribution for bjets of semiTt ",71,-3.55,3.55,71,-3.55,3.55);
    hPhi     = new TH2F("hPhi"," MC(X) vs reco(Y) , phi distribution for bjets of semiTt ",63,-3.15,3.15,71,-3.15,3.15);
    Pt_Resol = new TH2F("hPt_Resol", " Pt MC(X) vs Resol_Pt(Y) distribution for bjets of semiTt",60,0,300, 200,-0.995,1.005 );
 } 

 /// Destructor
 virtual ~HTOP7() {

    delete bDis_dR_RelPt;
    delete bDis_phi_d0;

    delete hEta;
    delete hPhi;
    delete Pt_Resol;

 }

 void Fill7d( double phi, double d0 ){
    bDis_phi_d0->Fill( phi, d0 );
 }

 void Fill7h(double eta_mc, double phi_mc, double eta_rc, double phi_rc, double pt, double ptResol) {
    hEta->Fill(eta_mc, eta_rc);
    hPhi->Fill(phi_mc, phi_rc);
    Pt_Resol->Fill(pt, ptResol);
 }
 void Fill7i( double dR, double RelPt ) {
    bDis_dR_RelPt->Fill( dR, RelPt );
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );


    bDis_dR_RelPt->Write();
    bDis_phi_d0->Write();

    hEta->Write();
    hPhi->Write();
    Pt_Resol->Write();
  
 }

  TH2F *bDis_dR_RelPt;
  TH2F *bDis_phi_d0;

  TH2F *hEta;
  TH2F *hPhi;
  TH2F *Pt_Resol;

};

class HTOP8 {
public:
 
 HTOP8() {

    hRes_Pt    = new TH1F("hRes_Pt", " Pt(jet) - Pt(q) / Pt(q)   ",400,-2.,2.);

    hdRWjj  = new TH2F("hdRWjj" ,"dR(Wq1, Wq2) , dR ( matched Wj1 , Wj2)  " ,200, -0.025,9.975 ,200, -0.025,9.975 );  

    bDis_dRWjMu  = new TH2F("bDis_dRWjMu","bCandidates jets JProb(x),TkCount(y)", 150,-0.1,1.4, 500,-20.,80.);

    // parton level and best matched reco comparison
    hEta     = new TH2F("hEta"," MC(X) vs reco(Y) , eta distribution for wjets of semiTt ",71,-3.55,3.55,71,-3.55,3.55);
    hPhi     = new TH2F("hPhi"," MC(X) vs reco(Y) , phi distribution for wjets of semiTt ",63,-3.15,3.15,71,-3.15,3.15);
    Pt_Resol = new TH2F("hPt_Resol", " Pt MC(X) vs Resol_Pt(Y) distribution for wjets of semiTt",60,0,300, 200,-0.995,1.005 );
 } 

 /// Destructor
 virtual ~HTOP8() {


    delete hRes_Pt;
    delete hdRWjj;

    delete bDis_dRWjMu;

    delete hEta;
    delete hPhi;
    delete Pt_Resol;
 }

 void Fill8c( double Res_Pt,  double bDis, double dR_WjMu ){
    hRes_Pt->Fill(Res_Pt);
    bDis_dRWjMu->Fill( bDis, dR_WjMu );
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

    hRes_Pt->Write();
    hdRWjj->Write();

    bDis_dRWjMu->Write();

    hEta->Write();
    hPhi->Write();
    Pt_Resol->Write();
 }


  TH1F *hRes_Pt;
  TH2F *hdRWjj;

  TH2F *bDis_dRWjMu;

  TH2F *hEta;
  TH2F *hPhi;
  TH2F *Pt_Resol;
};


class HTOP9 {
public:
 
 HTOP9() {

    EvtProb   = new TH1F("EvtProb"," Highest prob for permutation selection ", 200, 0., 1. );
    hEvtEff   = new TH1F("hEvtEff","(pass,fail)=> had(0,0.5) Mu(1,1.5) dilep(2,2.5) Ele(3,3.5), Tau(4,4.5), other(5,5.5)", 12,-0.25,5.75 );
    hObjEff   = new TH1F("hObjEff"," object selection Eff", 10, -0.5,9.5 );
    hbJetEff   = new TH1F("hbJetEff"," # of correctly matched bjet ", 4,-0.5,3.5 );
    hWJetEff   = new TH1F("hWJetEff"," # of correctly matched Wjet ", 4,-0.5,3.5 );
    hWLepEff   = new TH1F("hWLepEff"," # of correctly matched lep. ", 4,-0.5,3.5 );
    hHLTBits   = new TH2F("hHLTBits","HLT Trigger bits(names) for hadronic tt", 166, -0.5, 165.5, 6,-0.5,5.5 );
    NJ_NBTags  = new TH2F("NJ_NBTags"," NJets vs NBTags", 10, 0.5, 10.5, 10, -0.5, 9.5 );

 } 

 /// Destructor
 virtual ~HTOP9() {

    delete EvtProb;
    delete hEvtEff;
    delete hObjEff;
    delete hbJetEff;
    delete hWJetEff;
    delete hWLepEff;
    delete hHLTBits;
    delete NJ_NBTags;
 }

 void Fill9a( double prob ){
    EvtProb->Fill(prob);
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
 void Fill9j( int nb, int nw, int nl ) {
    hbJetEff->Fill( nb );
    hWJetEff->Fill( nw );
    hWLepEff->Fill( nl );
 }

 void Fill9k( int hlt, int topo ) {
    hHLTBits->Fill( hlt, topo );
 }


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    EvtProb->Write();
    hEvtEff->Write();
    hObjEff->Write();
    hbJetEff->Write();
    hWJetEff->Write();
    hWLepEff->Write();
    hHLTBits->Write();
    NJ_NBTags->Write();

 }

  TH1F *EvtProb; 
  TH1F *hEvtEff; 
  TH1F *hObjEff; 
  TH1F *hbJetEff;
  TH1F *hWJetEff;
  TH1F *hWLepEff;
  TH2F *hHLTBits;
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

    allTmass_pt    = new TH2D(N1+"_allTmass_pt",   " Top Mass vs Pt for all solutions" ,480,0.,480., 480,0,480 );
    highPtTmass_pt = new TH2D(N1+"_highPtTmass_pt"," highest Pt Top Mass vs Pt"        ,480,0.,480., 480,0,480 );
    selTmass_Wmass = new TH2D(N1+"_selTmass_Wmass"," Top Mass vs W Mass for seleted solution" ,480,0.,480., 480,0,480 );
    selTmass_pt    = new TH2D(N1+"_selTmass_pt",   " Top Mass vs Pt for seleted solution >3J" ,480,0.,480., 480,0,480 );
    selTmass_pt0   = new TH2D(N1+"_selTmass_pt0",  " Top Mass vs Pt for seleted solution =4J" ,480,0.,480., 480,0,480 );
    selTmass_pt1   = new TH2D(N1+"_selTmass_pt1", " Top Mass vs Pt for seleted solution =5J"  ,480,0.,480., 480,0,480 );
    selTmass_pt2   = new TH2D(N1+"_selTmass_pt2", " Top Mass vs Pt for seleted solution =6J"  ,480,0.,480., 480,0,480 );
    selTmass_pt3   = new TH2D(N1+"_selTmass_pt3", " Top Mass vs Pt for seleted solution >6J"  ,480,0.,480., 480,0,480 );

    hRecott0   = new TH2D(N2+"hRecott0","Reco top mass w/ 4j", 480,0.,480., 480,0.,480.);
    hRecott1   = new TH2D(N2+"hRecott1","Reco top mass w/ 5j", 480,0.,480., 480,0.,480.);
    hRecott2   = new TH2D(N2+"hRecott2","Reco top mass w/ 6j", 480,0.,480., 480,0.,480.);
    hRecott3   = new TH2D(N2+"hRecott3","Reco top mass >= 7j", 480,0.,480., 480,0.,480.);
    hRecott    = new TH2D(N2+"hRecott" ,"Reco top mass -all ", 480,0.,480., 480,0.,480.);

    PtRecott0   = new TH1F(N2+"PtRecott0","Pt of Reco top mass w/ 4j", 480,0.,480.);
    PtRecott1   = new TH1F(N2+"PtRecott1","Pt of Reco top mass w/ 5j", 480,0.,480.);
    PtRecott2   = new TH1F(N2+"PtRecott2","Pt of Reco top mass w/ 6j", 480,0.,480.);
    PtRecott3   = new TH1F(N2+"PtRecott3","Pt of Reco top mass >= 7j", 480,0.,480.);
    PtRecott    = new TH1F(N2+"PtRecott" ,"Pt of Reco top mass -all ", 480,0.,480.);

    WRecott0   = new TH2D(N2+"WRecott0","W vs dM of Reco top mass w/ 4j", 500,0,1000, 400,-199.5,200.5);
    WRecott1   = new TH2D(N2+"WRecott1","W vs dM of Reco top mass w/ 5j", 500,0,1000, 400,-199.5,200.5);
    WRecott2   = new TH2D(N2+"WRecott2","W vs dM of Reco top mass w/ 6j", 500,0,1000, 400,-199.5,200.5);
    WRecott3   = new TH2D(N2+"WRecott3","W vs dM of Reco top mass >= 7j", 500,0,1000, 400,-199.5,200.5);
    WRecott    = new TH2D(N2+"WRecott" ,"W vs dM of Reco top mass  -all", 500,0,1000, 400,-199.5,200.5);

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

    allTmass_pt    = (TH2D *) file->Get(theFolder+N1+"_allTmass_pt");
    highPtTmass_pt = (TH2D *) file->Get(theFolder+N1+"_highPtTmass_pt");
    selTmass_Wmass = (TH2D *) file->Get(theFolder+N1+"_selTmass_Wmass");
    selTmass_pt    = (TH2D *) file->Get(theFolder+N1+"_selTmass_pt");
    selTmass_pt0   = (TH2D *) file->Get(theFolder+N1+"_selTmass_pt0");
    selTmass_pt1   = (TH2D *) file->Get(theFolder+N1+"_selTmass_pt1");
    selTmass_pt2   = (TH2D *) file->Get(theFolder+N1+"_selTmass_pt2");
    selTmass_pt3   = (TH2D *) file->Get(theFolder+N1+"_selTmass_pt3");


    hRecott0   = (TH2D *) file->Get(theFolder+N2+"_hRecott0");
    hRecott1   = (TH2D *) file->Get(theFolder+N2+"_hRecott1");
    hRecott2   = (TH2D *) file->Get(theFolder+N2+"_hRecott2");
    hRecott3   = (TH2D *) file->Get(theFolder+N2+"_hRecott3");
    hRecott    = (TH2D *) file->Get(theFolder+N2+"_hRecott");
    PtRecott0  = (TH1F *) file->Get(theFolder+N2+"_PtRecott0");
    PtRecott1  = (TH1F *) file->Get(theFolder+N2+"_PtRecott1");
    PtRecott2  = (TH1F *) file->Get(theFolder+N2+"_PtRecott2");
    PtRecott3  = (TH1F *) file->Get(theFolder+N2+"_PtRecott3");
    PtRecott   = (TH1F *) file->Get(theFolder+N2+"_PtRecott");
    WRecott0   = (TH2D *) file->Get(theFolder+N2+"_WRecott0");
    WRecott1   = (TH2D *) file->Get(theFolder+N2+"_WRecott1");
    WRecott2   = (TH2D *) file->Get(theFolder+N2+"_WRecott2");
    WRecott3   = (TH2D *) file->Get(theFolder+N2+"_WRecott3");
    WRecott    = (TH2D *) file->Get(theFolder+N2+"_WRecott");

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
 void Fill11c( double tmass, double pt, double wmass, double weight ) {
    selTmass_pt->Fill( tmass, pt, weight );
    selTmass_Wmass->Fill( tmass, wmass, weight );
 }
 void Fill11c0( double tmass, double pt, double weight ) {
    selTmass_pt0->Fill( tmass, pt, weight );
 }
 void Fill11c1( double tmass, double pt, double weight ) {
    selTmass_pt1->Fill( tmass, pt, weight );
 }
 void Fill11c2( double tmass, double pt, double weight ) {
    selTmass_pt2->Fill( tmass, pt, weight );
 }
 void Fill11c3( double tmass, double pt, double weight ) {
    selTmass_pt3->Fill( tmass, pt, weight );
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

  TH2D *allTmass_pt;
  TH2D *highPtTmass_pt;
  TH2D *selTmass_Wmass;
  TH2D *selTmass_pt;
  TH2D *selTmass_pt0;
  TH2D *selTmass_pt1;
  TH2D *selTmass_pt2;
  TH2D *selTmass_pt3;


  TH2D *hRecott0;
  TH2D *hRecott1;
  TH2D *hRecott2;
  TH2D *hRecott3;
  TH2D *hRecott;
  TH1F *PtRecott0;
  TH1F *PtRecott1;
  TH1F *PtRecott2;
  TH1F *PtRecott3;
  TH1F *PtRecott;
  TH2D *WRecott0;
  TH2D *WRecott1;
  TH2D *WRecott2;
  TH2D *WRecott3;
  TH2D *WRecott;

 };


#endif
