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

    
    NJets  = new TH1F("NJets", " # of Jets",65,-0.5,64.5);
    Jet1Et = new TH1F("Jet1Et", " ET of the 1st highest Et jet ", 500,0.,500. ); 
    Jet2Et = new TH1F("Jet2Et", " ET of the 2nd highest Et jet ", 500,0.,500. ); 
    Jet3Et = new TH1F("Jet3Et", " ET of the 3rd highest Et jet ", 500,0.,500. ); 
    Jet4Et = new TH1F("Jet4Et", " ET of the 4th highest Et jet ", 500,0.,500. ); 

    LeadJetEt = new TH1F("LeadJetEt", " ET of the 1st highest Et jet ", 500,0.,500. ); 

    dRjj = new TH1F("dRjj"," dR( Jx, Jy ) ", 200,-0.025, 9.975);

    dRMuJ_ERatio = new TH2F("dRMuJ_ERatio", " min_dR(mu,jet), muE/jetE ", 200, 0., 5., 100, -0.05, 9.95);
    dRMuJ_RelPt  = new TH2F("dRMuJ_RelPt", " min_dR(mu,jet), RelPt ", 200, 0., 5., 100, 0., 200 );

    dRElJ_ERatio = new TH2F("dRElJ_ERatio", " min_dR(Ele,jet), EleE/jetE ", 200, 0., 5., 100, -0.05, 9.95);
    dRElJ_RelPt  = new TH2F("dRElJ_RelPt", " min_dR(Ele,jet), RelPt ", 200, 0., 5., 100, 0., 200 );

    m3_j3    = new TH2F("m3_j3", " m3 vs J3 ET ", 500, 0., 500., 500, 0., 500. );

    allJ_selJ   = new TH2F("allJ_selJ", "N of All jets vs N of selected jets in a event", 21, -0.5, 20.5, 21, -0.5, 20.5);
    nBJets_selJ = new TH2F("nBJets_selJ", "N of BJets from Selected Jets ", 11, -0.5, 10.5, 11, -0.5, 10.5 );

    NloserJets = new TH1F("NloserJets", " # of Jets",35,-0.5,34.5);
    loserJetEt = new TH1F("loserJetEt", " ET of the highest Et jet not be selected", 100,0.,50. ); 

 } 

 /// Destructor
 virtual ~HOBJ1() {

    delete dRjj;
    delete NJets;

    delete LeadJetEt; 

    delete Jet1Et; 
    delete Jet2Et; 
    delete Jet3Et; 
    delete Jet4Et; 

    delete dRMuJ_ERatio;
    delete dRMuJ_RelPt;

    delete dRElJ_ERatio;
    delete dRElJ_RelPt;

    delete m3_j3;

    delete allJ_selJ;
    delete nBJets_selJ;


    delete NloserJets;
    delete loserJetEt;
 }

 void Fill_1a( double dR ) {
    dRjj->Fill( dR );
 }
 void Fill_1b( double dR, double relPt, double muE, double jetE ) {
     dRMuJ_ERatio->Fill( dR, muE/jetE );
     dRMuJ_RelPt->Fill( dR, relPt );
 }
 void Fill_1c( double et1, double et2, double et3, double et4, double m3 ) {
    Jet1Et->Fill( et1 );
    Jet2Et->Fill( et2 );
    Jet3Et->Fill( et3 );
    Jet4Et->Fill( et4 );
    m3_j3->Fill( m3, et3 );
 }
 void Fill_1d( double et1 ) {
    LeadJetEt->Fill(et1);
 }
 void Fill_1f( double njets ) {
    NJets->Fill(njets);
 }

 void Fill_1g( int allJ, int selJ, int nBJ ) {
    allJ_selJ->Fill( allJ, selJ );
    nBJets_selJ->Fill( nBJ, selJ );
 }

 void Fill_1h( double dR, double relPt, double eleE, double jetE ) {
     dRElJ_ERatio->Fill( dR, eleE/jetE );
     dRElJ_RelPt->Fill( dR, relPt );
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
 
    LeadJetEt->Write();   
 
    dRMuJ_ERatio->Write();
    dRMuJ_RelPt->Write();

    dRElJ_ERatio->Write();
    dRElJ_RelPt->Write();

    m3_j3->Write();
    allJ_selJ->Write();
    nBJets_selJ->Write();

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

  TH1F *LeadJetEt;

  TH2F *dRMuJ_ERatio;
  TH2F *dRMuJ_RelPt;

  TH2F *dRElJ_ERatio;
  TH2F *dRElJ_RelPt;

  TH2F *m3_j3;

  TH2F *allJ_selJ;
  TH2F *nBJets_selJ;

  TH1F *NloserJets;
  TH1F *loserJetEt;

};

class HOBJ2 {
public:
 
 HOBJ2() {
    PtResol_pat = new TH1F("PtResol_pat","Pt Resol, patMET ", 315, -3.15, 3.15 );
    PtResol_evt = new TH1F("PtResol_evt","Pt Resol, evtMET ", 315, -3.15, 3.15 );
    PtResol_gen = new TH1F("PtResol_gen","Pt Resol, genMET ", 315, -3.15, 3.15 );

    gen_pat_Resol = new TH1F("gen_pat_Resol","Pat Pt Resol w.r.t genMET ", 315, -3.15, 3.15 );
    gen_evt_Resol = new TH1F("gen_evt_Resol","Event Pt Resol w.r.t genMET ", 315, -3.15, 3.15 );

    dPhiResol_pat = new TH1F("dPhiResol_pat", "dPhi(neu, MET),  pat ",315, -3.15, 3.15 );
    dPhiResol_evt = new TH1F("dPhiResol_evt", "dPhi(neu, MET),  evt ",315, -3.15, 3.15 );
    dPhiResol_gen = new TH1F("dPhiResol_gen", "dPhi(neu, MET),  gen ",315, -3.15, 3.15 );

    gen_pat_dPhi = new TH1F("gen_pat_dPhi", " dPhi(pat, gen) ", 315, -3.15, 3.15 );
    gen_evt_dPhi = new TH1F("gen_evt_dPhi", " dPhi(evt, gen) ", 315, -3.15, 3.15 );
 
 }

 /// Destructor
 virtual ~HOBJ2() {

    delete PtResol_pat;
    delete PtResol_evt;
    delete PtResol_gen;
    delete gen_pat_Resol;
    delete gen_evt_Resol;

    delete dPhiResol_pat;
    delete dPhiResol_evt;
    delete dPhiResol_gen;
    delete gen_pat_dPhi;
    delete gen_evt_dPhi;

 }

 void Fill_2a( double patResol,  double evtResol, double genResol, double genpat_Res, double genevt_Res,
       double dPhi_neu_pat, double dPhi_neu_evt, double dPhi_neu_gen, double genpat_dPhi, double genevt_dPhi ){

    PtResol_pat->Fill( patResol );
    PtResol_evt->Fill( evtResol );
    PtResol_gen->Fill( genResol );
    gen_pat_Resol->Fill( genpat_Res);
    gen_evt_Resol->Fill( genevt_Res);

    dPhiResol_pat->Fill( dPhi_neu_pat );
    dPhiResol_evt->Fill( dPhi_neu_evt );
    dPhiResol_gen->Fill( dPhi_neu_gen );
    gen_pat_dPhi->Fill( genpat_dPhi );
    gen_evt_dPhi->Fill( genevt_dPhi );

 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    PtResol_pat->Write();
    PtResol_evt->Write();
    PtResol_gen->Write();
    gen_pat_Resol->Write();
    gen_evt_Resol->Write();

    dPhiResol_pat->Write();
    dPhiResol_evt->Write();
    dPhiResol_gen->Write();
    gen_pat_dPhi->Write();
    gen_evt_dPhi->Write();

    file->cd("../");
 }

  TH1F *PtResol_pat;
  TH1F *PtResol_evt;
  TH1F *PtResol_gen;
  TH1F *gen_pat_Resol;
  TH1F *gen_evt_Resol;
  TH1F *dPhiResol_pat;
  TH1F *dPhiResol_evt;
  TH1F *dPhiResol_gen;
  TH1F *gen_pat_dPhi;
  TH1F *gen_evt_dPhi;


};

class HOBJ3 {
public:
 
 HOBJ3() {

    allMuPt     = new TH1F("allMuPt",  " pt " ,200,0.,200.);
    allMuEta    = new TH1F("allMuEta", " Eta " ,200,-5., 5.);
    allMuRelIso = new TH1F("allMuRelIso", " RelIso " ,150, -0.2, 1.3);

    nHitX2_sta  =  new TH2F("nHitX2_sta"," nSeg vs X2 ", 50,0,50, 500, 0, 50);
    nHitX2_trk  =  new TH2F("nHitX2_trk"," nHit vs X2 ", 50,0,50, 500, 0, 50);
    nHitX2_glb  =  new TH2F("nHitX2_glb"," n glbHit vs X2 ", 50,0,50, 500, 0, 50);

    IpSig_trk   = new TH1F("IpSig_trk",  "IP Significance " ,125, -30.,30.);
    glbMuPt     = new TH1F("glbMuPt",  " pt " ,200,0.,200.);
    glbMuEta    = new TH1F("glbMuEta", " Eta " ,200, -5., 5.);
    glbMuRelIso = new TH1F("glbMuRelIso", " RelIso " ,150, -0.2, 1.3);

   // Isolation
    Pt_emIso    = new TH2F("Pt_emIso", "pt vs iso emE ",200,0.,200., 100, 0., 100.);
    Pt_hdIso    = new TH2F("Pt_hdIso", "pt vs iso hdE ",200,0.,200., 100, 0., 100.);
    Pt_tkIso    = new TH2F("Pt_tkIso", "pt vs iso tkP ",200,0.,200., 100, 0., 100.);
    Eta_Iso     = new TH2F("Eta_Iso", "eta vs Iso ",61, -3.05,3.05 , 150, -0.2, 1.3);
    //Pt_CaloIso  = new TH2F("Pt_CaloIso","pt vs calo iso ",200,0.,200., 100, 0., 100.);
    //Pt_Iso      = new TH2F("Pt_Iso"  , "pt vs Iso value ",200,0.,200., 150, -0.2, 1.3 );
    //Eta_CaloIso = new TH2F("Eta_CaloIso", "eta vs CaloIso ",61, -3.05,3.05 , 150, -0.2, 1.3);

    glb_trk_Pt = new TH2F("glb_trk_Pt"," global Pt vs tracker Pt ", 200,0.,200., 200, 0., 200. );
    dPt_eta    = new TH2F("dPt_eta", " eta vs Glb Pt -trk Pt ", 61, -3.05,3.05 , 201, -1.005, 1.005 );

    allMu_selMu = new TH2F("allMu_selMu", " N of muons vs N SelMuon in a event ", 21, -0.5, 20.5 , 21, -0.5, 20.5);
 }

 /// Destructor
 virtual ~HOBJ3() {

    delete allMuPt;
    delete allMuEta;
    delete allMuRelIso;

    delete nHitX2_sta;
    delete nHitX2_trk;
    delete nHitX2_glb;

    delete IpSig_trk;
    delete glbMuPt;
    delete glbMuEta;
    delete glbMuRelIso;

    delete Pt_emIso;
    delete Pt_hdIso;
    delete Pt_tkIso;
    //delete Pt_CaloIso;
    //delete Pt_Iso;
    delete Eta_Iso;
    //delete Eta_CaloIso;

    delete glb_trk_Pt;
    delete dPt_eta;

    delete allMu_selMu;
 }

 /*
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
 */ 

 void Fill_3c(int allMu, int selMu ) {
    allMu_selMu->Fill( allMu, selMu );
 }  

 void Fill_3d( double pt, double eta, double relIso, double emIso, double hadIso, double trkIso ) {
    allMuPt->Fill(pt);
    allMuEta->Fill(eta);
    allMuRelIso->Fill(relIso);  
    Pt_emIso->Fill(pt, emIso);
    Pt_hdIso->Fill(pt, hadIso);
    Pt_tkIso->Fill(pt, trkIso);
    Eta_Iso->Fill(eta, relIso);
 }

 void Fill_3e( int nSeg, double X2 ) {
    nHitX2_sta->Fill(nSeg, X2);
 }
 void Fill_3f( int nHit, double X2, double ipSig ) {
    nHitX2_trk->Fill(nHit, X2);
    IpSig_trk->Fill( ipSig );
 }
 void Fill_3g( int nHit, double X2, double glbPt, double eta, double relIso, double trkPt, double pt_Res ) {
    nHitX2_glb->Fill(nHit, X2);
    glbMuPt->Fill(glbPt);
    glbMuEta->Fill(eta);
    glbMuRelIso->Fill(relIso);  
    glb_trk_Pt->Fill(glbPt, trkPt);
    dPt_eta->Fill( eta, pt_Res );
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    allMuPt->Write();
    allMuEta->Write();
    allMuRelIso->Write();

    nHitX2_sta->Write();
    nHitX2_trk->Write();
    nHitX2_glb->Write();
    IpSig_trk->Write();

    glbMuPt->Write();
    glbMuEta->Write();
    glbMuRelIso->Write();

    Pt_emIso->Write();
    Pt_hdIso->Write();
    Pt_tkIso->Write();
    Eta_Iso->Write();
    //Pt_CaloIso->Write();
    //Pt_Iso->Write();
    //Eta_CaloIso->Write();

    glb_trk_Pt->Write();
    dPt_eta->Write();

    allMu_selMu->Write();
 }

  TH1F *allMuPt;
  TH1F *allMuEta;
  TH1F *allMuRelIso;

  TH2F *nHitX2_sta;
  TH2F *nHitX2_trk;
  TH2F *nHitX2_glb;

  TH1F *IpSig_trk;

  TH1F *glbMuPt;
  TH1F *glbMuEta;
  TH1F *glbMuRelIso;

  TH2F *Pt_emIso;
  TH2F *Pt_hdIso;
  TH2F *Pt_tkIso;
  //TH2F *Pt_CaloIso;
  //TH2F *Pt_Iso;
  TH2F *Eta_Iso;
  //TH2F *Eta_CaloIso;

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


class HBJet {
public:
 
 HBJet() {

    //dRbMuJ_RelPt  = new TH2F("dRbMuJ_RelPt", " min_dR( MuInJet,bjet), RelPt ", 200, 0., 5., 100, 0., 200 );
    //dRiMuJ_RelPt  = new TH2F("dRiMuJ_RelPt", " min_dR( Mu ,jet), RelPt ", 200, 0., 5., 100, 0., 200 );
    //dR_relIso     = new TH2F("dR_relIso", " dR( MuInJet,bjet), RelIso " , 200, 0., 5., 25, -10.5, 14.5 );
    
    bDis_pt   = new TH2F("bDis_pt","b Discriminator for all jets ", 60,-10.5, 49.5, 500, 0., 500. );
    d0_pt     = new TH2F("d0_pt", " ave. d0 (X) vs jet pt (Y)", 500, -30., 30., 500, 0., 500. );
    bDis_pt1  = new TH2F("bDis_pt1","b Discriminator for all jets ", 60,-10.5, 49.5, 500, 0., 500. );
    d0_pt1    = new TH2F("d0_pt1", " ave. d0 (X) vs jet pt (Y)", 500, -30., 30., 500, 0., 500. );

    nJ_nbJ    = new TH2F("nJ_nbJ", "N of b jets vs N of jets", 10, -0.5, 9.5, 25, -0.5, 24.5);
    nJ_nbJ2   = new TH2F("nJ_nbJ2", "N of cleaned b jets vs N of jets", 10, -0.5, 9.5, 25, -0.5, 24.5);
 } 

 /// Destructor
 virtual ~HBJet() {

    //delete dRbMuJ_RelPt;
    //delete dRiMuJ_RelPt;
    //delete dR_relIso;

    delete bDis_pt;
    delete d0_pt;
    delete bDis_pt1;
    delete d0_pt1;

    delete nJ_nbJ;
    delete nJ_nbJ2;
 }

 /*
 void Fill_b1( double dR, double relPt, double relIso ) {
     dRbMuJ_RelPt->Fill( dR, relPt );
     dR_relIso->Fill( dR, relIso );
 }

 void Fill_b2( double dR, double relPt ) {
     dRiMuJ_RelPt->Fill( dR, relPt );
 }
 */
 void Fill_b3( double bdis, double d0Ave, double pt ) {
    bDis_pt->Fill( bdis, pt );
    d0_pt->Fill( bdis, pt );
 }

 void Fill_b4( double bdis, double d0Ave, double pt ) {
    bDis_pt1->Fill( bdis, pt );
    d0_pt1->Fill( bdis, pt );
 }

 void Fill_b5( int nBJet, int nJet, int nBJet2 ) {
    nJ_nbJ->Fill(  nBJet, nJet );
    nJ_nbJ2->Fill(  nBJet2, nJet );
 }

 void Write( TString theFolder , TFile* file  ) {

    file->cd( theFolder );

    //dRbMuJ_RelPt->Write();
    //dRiMuJ_RelPt->Write();
    //dR_relIso->Write();

    bDis_pt->Write();
    bDis_pt1->Write();
 
    d0_pt->Write();
    d0_pt1->Write();

    nJ_nbJ->Write(); 
    nJ_nbJ2->Write(); 

    file->cd("../");

 }

  //TH2F *dRbMuJ_RelPt;
  //TH2F *dRiMuJ_RelPt;
  //TH2F *dR_relIso;

  TH2F *bDis_pt;
  TH2F *bDis_pt1;

  TH2F *d0_pt;
  TH2F *d0_pt1;

  TH2F *nJ_nbJ;
  TH2F *nJ_nbJ2;

};

#endif
