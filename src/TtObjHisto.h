#ifndef TtObjHisto_H
#define TtObjHisto_H

/** \class TopObjectHistograms
 *  Collection of histograms for TopObjects Analysis test.
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

class HOBJ1 {
public:
 
 HOBJ1() {

    dRjj = new TH1F("dRjj"," dR( Jx, Jy ) ", 200,-0.025, 9.975);
    
    NJets  = new TH1F("NJets", " # of Jets",65,-0.5,64.5);
    Jet1Et = new TH1F("Jet1Et", " ET of the 1st highest Et jet ", 500,0.,500. ); 
    Jet2Et = new TH1F("Jet2Et", " ET of the 2nd highest Et jet ", 500,0.,500. ); 
    Jet3Et = new TH1F("Jet3Et", " ET of the 3rd highest Et jet ", 500,0.,500. ); 
    Jet4Et = new TH1F("Jet4Et", " ET of the 4th highest Et jet ", 500,0.,500. ); 

    eta_njets  = new TH2F("eta_njets", " eta vs njets ", 91, -4.55, 4.55, 25, -0.5, 24.5);

    njets_dRMuJ  = new TH2F("njets_dRMuJ", " njets vs min_dR(mu,jet)  ", 20, -0.5, 19.5, 100, -0.05, 9.95);
    dRMuJ_RelPt  = new TH2F("dRMuJ_RelPt", " min_dR(mu,jet), RelPt  ", 200, 0., 5., 100, 0., 200 );
    muNjets_dPhi0= new TH2F("muNjets_dPhi0"," njets vs dPhi(mu,jet1)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    muNjets_dPhi1= new TH2F("muNjets_dPhi1"," njets vs dPhi(mu,jet2)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);
    muNjets_dPhi2= new TH2F("muNjets_dPhi2"," njets vs dPhi(mu,jet3)  ", 20, -0.5, 19.5, 160, -0.05, 3.15);

    m3_j3    = new TH2F("m3_j3", " m3 vs J3 ET ", 500, 0., 500., 500, 0., 500. );

    allJ_selJ = new TH2F("allJ_selJ", "N of All jets vs N of selected jets in a event", 21, -0.5, 20.5, 21, -0.5, 20.5);
    nBJets    = new TH1F("nBJets", "N of BJets from Selected Jets ", 11, -0.5, 10.5);

    isoEleCut = new TH2F("isoEleCut","N of isoEle vs N of SelectedJets  after isoMu =1 cut ", 21, -0.5, 20.5, 21, -0.5, 20.5 );

    NloserJets = new TH1F("NloserJets", " # of Jets",35,-0.5,34.5);
    loserJetEt = new TH1F("loserJetEt", " ET of the highest Et jet not be selected", 100,0.,50. ); 
 } 

 /// Destructor
 virtual ~HOBJ1() {

    delete dRjj;
    delete NJets;

    delete Jet1Et; 
    delete Jet2Et; 
    delete Jet3Et; 
    delete Jet4Et; 
    delete eta_njets;

    delete njets_dRMuJ;
    delete dRMuJ_RelPt;

    delete m3_j3;

    delete allJ_selJ;
    delete nBJets;

    delete isoEleCut;

    delete NloserJets;
    delete loserJetEt;
 }

 void Fill_1a( double dR ) {
    dRjj->Fill( dR );
 }
 void Fill_1b( int njet, double dR, double relPt ) {
     njets_dRMuJ->Fill( njet, dR );
     dRMuJ_RelPt->Fill( dR, relPt );
 }
 void Fill_1c( double et1, double et2, double et3, double et4, double m3 ) {
    Jet1Et->Fill( et1 );
    Jet2Et->Fill( et2 );
    Jet3Et->Fill( et3 );
    Jet4Et->Fill( et4 );
    m3_j3->Fill( m3, et3 );
 }
 void Fill_1d( double njets, double eta ) {
    eta_njets->Fill( eta, njets );
 }
 void Fill_1e( double njets, double dphi0, double dphi1, double dphi2 ) {
    muNjets_dPhi0->Fill( njets, dphi0 );
    muNjets_dPhi1->Fill( njets, dphi1 );
    muNjets_dPhi2->Fill( njets, dphi2 );
 } 
 void Fill_1f( double njets ) {
    NJets->Fill(njets);
 }

 void Fill_1g( int allJ, int selJ, int nBJ ) {
    allJ_selJ->Fill( allJ, selJ );
    nBJets->Fill( nBJ );
 }

 void Fill_1h( int NIsoEle, int NSelJet ){
    isoEleCut->Fill( NIsoEle, NSelJet );
 }

 void Fill_1i( int NJet, double jetEt ){
    NloserJets->Fill( NJet  );
    loserJetEt->Fill( jetEt  );
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    NJets->Write();
    dRjj->Write();
 
    Jet1Et->Write();
    Jet2Et->Write();
    Jet3Et->Write();
    Jet4Et->Write();
    eta_njets->Write();
    
    njets_dRMuJ->Write();
    dRMuJ_RelPt->Write();
    muNjets_dPhi0->Write();
    muNjets_dPhi1->Write();
    muNjets_dPhi2->Write();

    m3_j3->Write();
    allJ_selJ->Write();
    nBJets->Write();
    isoEleCut->Write();

    NloserJets->Write();
    loserJetEt->Write();

    file->cd("../");

 }

  TH1F *dRjj;

  TH1F *NJets;

  TH1F *Jet1Et;
  TH1F *Jet2Et;
  TH1F *Jet3Et;
  TH1F *Jet4Et;
  TH2F *eta_njets;

  TH2F *njets_dRMuJ;
  TH2F *dRMuJ_RelPt;
  TH2F *muNjets_dPhi0;
  TH2F *muNjets_dPhi1;
  TH2F *muNjets_dPhi2;

  TH2F *m3_j3;

  TH2F *allJ_selJ;
  TH1F *nBJets;

  TH2F *isoEleCut;
  TH1F *NloserJets;
  TH1F *loserJetEt;

};

class HOBJ2 {
public:
 
 HOBJ2() {
    PtResol_pat = new TH1F("PtResol_pat","Pt Resol, patMET ", 150, -1.01, 1.99 );
    PtResol_tc  = new TH1F("PtResol_tc", "Pt Resol, tcMET  ", 150, -1.01, 1.99 );
    PtResol_evt = new TH1F("PtResol_evt","Pt Resol, evtMET ", 150, -1.01, 1.99 );
    PtResol_gen = new TH1F("PtResol_gen","Pt Resol, genMET ", 150, -1.01, 1.99 );

    gen_pat_Resol = new TH1F("gen_pat_Resol","Pat Pt Resol w.r.t genMET ", 150, -1.01, 1.99 );
    gen_evt_Resol = new TH1F("gen_evt_Resol","Event Pt Resol w.r.t genMET ", 150, -1.01, 1.99 );

    dPhiResol_pat = new TH1F("dPhiResol_pat", "dPhi(neu, MET),  pat ",315, -3.15, 3.15 );
    dPhiResol_tc  = new TH1F("dPhiResol_tc", "dPhi(neu, MET),   tc  ",315, -3.15, 3.15 );
    dPhiResol_evt = new TH1F("dPhiResol_evt", "dPhi(neu, MET),  evt ",315, -3.15, 3.15 );
    dPhiResol_gen = new TH1F("dPhiResol_gen", "dPhi(neu, MET),  gen ",315, -3.15, 3.15 );

    gen_pat_dPhi = new TH1F("gen_pat_dPhi", " dPhi(pat, gen) ", 315, -3.15, 3.15 );
    gen_evt_dPhi = new TH1F("gen_evt_dPhi", " dPhi(evt, gen) ", 315, -3.15, 3.15 );
 
 }

 /// Destructor
 virtual ~HOBJ2() {

    delete PtResol_pat;
    delete PtResol_tc;
    delete PtResol_evt;
    delete PtResol_gen;
    delete gen_pat_Resol;
    delete gen_evt_Resol;

    delete dPhiResol_pat;
    delete dPhiResol_tc;
    delete dPhiResol_evt;
    delete dPhiResol_gen;
    delete gen_pat_dPhi;
    delete gen_evt_dPhi;

 }

 void Fill_2a( double patResol,  double evtResol, double tcResol, double genResol, double genpat_Res, double genevt_Res,
       double dPhi_neu_pat, double dPhi_neu_evt, double dPhi_neu_tc, double dPhi_neu_gen, double genpat_dPhi, double genevt_dPhi ){

    PtResol_pat->Fill( patResol );
    PtResol_tc->Fill( tcResol );
    PtResol_evt->Fill( evtResol );
    PtResol_gen->Fill( genResol );
    gen_pat_Resol->Fill( genpat_Res);
    gen_evt_Resol->Fill( genevt_Res);

    dPhiResol_pat->Fill( dPhi_neu_pat );
    dPhiResol_tc->Fill( dPhi_neu_tc );
    dPhiResol_evt->Fill( dPhi_neu_evt );
    dPhiResol_gen->Fill( dPhi_neu_gen );
    gen_pat_dPhi->Fill( genpat_dPhi );
    gen_evt_dPhi->Fill( genevt_dPhi );

 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    PtResol_pat->Write();
    PtResol_tc->Write();
    PtResol_evt->Write();
    PtResol_gen->Write();
    gen_pat_Resol->Write();
    gen_evt_Resol->Write();

    dPhiResol_pat->Write();
    dPhiResol_tc->Write();
    dPhiResol_evt->Write();
    dPhiResol_gen->Write();
    gen_pat_dPhi->Write();
    gen_evt_dPhi->Write();

    file->cd("../");
 }

  TH1F *PtResol_tc;
  TH1F *PtResol_pat;
  TH1F *PtResol_evt;
  TH1F *PtResol_gen;
  TH1F *gen_pat_Resol;
  TH1F *gen_evt_Resol;
  TH1F *dPhiResol_tc;
  TH1F *dPhiResol_pat;
  TH1F *dPhiResol_evt;
  TH1F *dPhiResol_gen;
  TH1F *gen_pat_dPhi;
  TH1F *gen_evt_dPhi;


};

class HOBJ3 {
public:
 
 HOBJ3() {

   // Isolation

    Pt_emIso    = new TH2F("Pt_emIso", "pt vs iso emE ",200,0.,200., 100, 0., 100.);
    Pt_hdIso    = new TH2F("Pt_hdIso", "pt vs iso hdE ",200,0.,200., 100, 0., 100.);
    Pt_tkIso    = new TH2F("Pt_tkIso", "pt vs iso tkP ",200,0.,200., 100, 0., 100.);
    Pt_CaloIso  = new TH2F("Pt_CaloIso","pt vs calo iso ",200,0.,200., 100, 0., 100.);
    Pt_Iso      = new TH2F("Pt_Iso"  , "pt vs Iso value ",200,0.,200., 150, -0.2, 1.3 );
    Eta_Iso     = new TH2F("Eta_Iso", "eta vs Iso ",61, -3.05,3.05 , 150, -0.2, 1.3);
    Eta_CaloIso = new TH2F("Eta_CaloIso", "eta vs CaloIso ",61, -3.05,3.05 , 150, -0.2, 1.3);

    glb_trk_Pt = new TH2F("glb_trk_Pt"," global Pt vs tracker Pt ", 200,0.,200., 200, 0., 200. );
    dPt_eta    = new TH2F("dPt_eta", " eta vs Glb Pt -trk Pt ", 61, -3.05,3.05 , 201, -1.005, 1.005 );

    allMu_selMu = new TH2F("allMu_selMu", " N of muons vs N IsoMuon in a event ", 21, -0.5, 20.5 , 21, -0.5, 20.5);
 }

 /// Destructor
 virtual ~HOBJ3() {

    delete Pt_emIso;
    delete Pt_hdIso;
    delete Pt_tkIso;
    delete Pt_CaloIso;
    delete Pt_Iso;
    delete Eta_Iso;
    delete Eta_CaloIso;

    delete glb_trk_Pt;
    delete dPt_eta;

    delete allMu_selMu;
 }

 void Fill_3a(double pt, double eta, double emE, double hdE, double tkP, double sumE, double Iso ){

    Pt_emIso->Fill(pt, emE);
    Pt_hdIso->Fill(pt, hdE);
    Pt_tkIso->Fill(pt, tkP);
    Pt_CaloIso->Fill(pt, sumE);
    Pt_Iso->Fill(pt, Iso);
    Eta_Iso->Fill(eta, Iso);
    Eta_CaloIso->Fill(eta, sumE);

 }

 void Fill_3b(double glbPt, double trkPt, double dPt, double eta ) {
    glb_trk_Pt->Fill(glbPt, trkPt);
    dPt_eta->Fill( eta, dPt);
 }

 void Fill_3c(int allMu, int isoMu ) {
    allMu_selMu->Fill( allMu, isoMu );
 }  


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    Pt_emIso->Write();
    Pt_hdIso->Write();
    Pt_tkIso->Write();
    Pt_CaloIso->Write();
    Pt_Iso->Write();
    Eta_Iso->Write();
    Eta_CaloIso->Write();

    glb_trk_Pt->Write();
    dPt_eta->Write();

    allMu_selMu->Write();
 }

  TH2F *Pt_emIso;
  TH2F *Pt_hdIso;
  TH2F *Pt_tkIso;
  TH2F *Pt_CaloIso;
  TH2F *Pt_Iso;
  TH2F *Eta_Iso;
  TH2F *Eta_CaloIso;

  TH2F *glb_trk_Pt;
  TH2F *dPt_eta;

  TH2F * allMu_selMu;

};

class HOBJ4 {
public:
 
 HOBJ4() {

   // Isolation
   hEta_Iso = new TH2F("hEta_Iso", "eta vs Iso for all ",61, -3.05,3.05 , 150, -0.2, 1.3);
   Pt_emIso = new TH2F("Pt_emIso", "pt vs iso emE for all",200,0.,200., 100, 0., 100.);
   Pt_hdIso = new TH2F("Pt_hdIso", "pt vs iso hdE for all",200,0.,200., 100, 0., 100.);
   Pt_tkIso = new TH2F("Pt_tkIso", "pt vs iso tkP for all",200,0.,200., 100, 0., 100.);
   hPt_sum  = new TH2F("hPt_sum",  "pt vs iso caloE for all",200,0.,200., 100, 0., 100.);
   hPt_Cal  = new TH2F("hPt_Cal",  "pt vs N of CaloDeposit in isoR3 for all ",200,0.,200., 50, -0.5, 49.5 );
   hPt_Trk  = new TH2F("hPt_Trk",  "pt vs N of Track in isoR3 for all ",200,0.,200., 50, -0.5, 49.5 );
   hPt_Iso  = new TH2F("hPt_Iso",  "pt vs Iso value for all",200,0.,200., 150, -0.2, 1.3 );
   Iso_EovP = new TH2F("Iso_EovP", "RelIso vs emE/P for all",55,-0.05, 1.05, 150, -0.25, 1.25 );
   Iso_HovE = new TH2F("Iso_HovE", "RelIso vs H/E  for all ",55,-0.05, 1.05, 60, -0.05, 0.25 );

   allEl_selEl = new TH2F("allEl_selEl", " N of Electrons vs N IsoElectrons in a event ", 21, -0.5, 20.5 , 21, -0.5, 20.5);

 }

 /// Destructor
 virtual ~HOBJ4() {

   delete Pt_emIso;
   delete Pt_hdIso;
   delete Pt_tkIso;
   delete hPt_sum;
   delete hPt_Cal;
   delete hPt_Trk;
   delete hPt_Iso;
   delete hEta_Iso;
   delete Iso_EovP;
   delete Iso_HovE;

   delete allEl_selEl;

 }

 void Fill_4a(double pt, double eta, double emE, double hdE, double tkP, double sumE, int nCal3, int ntrack3,
             double Iso, double HovE, double EovP ){

    Pt_emIso->Fill(pt, emE);
    Pt_hdIso->Fill(pt, hdE);
    Pt_tkIso->Fill(pt, tkP);
    hPt_sum->Fill(pt, sumE);
    hPt_Cal->Fill(pt, nCal3);
    hPt_Trk->Fill(pt, ntrack3);
    hPt_Iso->Fill(pt, Iso);
    hEta_Iso->Fill(eta, Iso);
    Iso_EovP->Fill(Iso, EovP);
    Iso_HovE->Fill(Iso, HovE);
 }

 void Fill_4b( int allEl, int isoEl ){
    allEl_selEl->Fill(allEl, isoEl);
 }  


 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    Pt_emIso->Write();
    Pt_hdIso->Write();
    Pt_tkIso->Write();
    hPt_sum->Write();
    hPt_Cal->Write();
    hPt_Trk->Write();
    hPt_Iso->Write();
    hEta_Iso->Write();
    Iso_EovP->Write();
    Iso_HovE->Write();

    allEl_selEl->Write();

 }

  TH2F *Pt_emIso;
  TH2F *Pt_hdIso;
  TH2F *Pt_tkIso;
  TH2F *hPt_sum;
  TH2F *hPt_Cal;
  TH2F *hPt_Trk;
  TH2F *hPt_Iso;
  TH2F *hEta_Iso;
  TH2F *Iso_EovP;
  TH2F *Iso_HovE;

  TH2F *allEl_selEl;

};

#endif
