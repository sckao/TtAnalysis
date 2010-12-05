#ifndef TtAnalysisNtuple_H
#define TtAnalysisNtuple_H

/** \class TtAnalysisNtuple
 *  Collection of Ntuples 
 *
 * Author: S.C. Kao  - UC Riverside
 */

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TtFormat.h"
#include <string>
#include <iostream>

typedef math::XYZTLorentzVector LorentzVector;

using namespace std;

// new solution tree
class SolNtp2 {
public:

 SolNtp2( TString treeName ) {

    solTree2 = new TTree( treeName ,"");
 
    solTree2->Branch("evtId"     ,&evtId,    "evtId/I");

    solTree2->Branch("nJ"    ,&nJ,     "nJ/I");
    solTree2->Branch("jpx"   ,&jpx,    "jpx[nJ]/D");
    solTree2->Branch("jpy"   ,&jpy,    "jpy[nJ]/D");
    solTree2->Branch("jpz"   ,&jpz,    "jpz[nJ]/D");
    solTree2->Branch("jE"    ,&jE,     "jE[nJ]/D");
    solTree2->Branch("bTh"   ,&bTh,    "bTh[nJ]/D");
    solTree2->Branch("emF"   ,&emF,    "emF[nJ]/D");
    solTree2->Branch("EvH"   ,&EvH,    "EvH[nJ]/D");
    solTree2->Branch("n90"   ,&n90,    "n90[nJ]/I");

    solTree2->Branch("nNu"    ,&nNu,   "nNu/I");
    solTree2->Branch("npx"    ,&npx,   "npx[nNu]/D");
    solTree2->Branch("npy"    ,&npy,   "npy[nNu]/D");
    solTree2->Branch("npz"    ,&npz,   "npz[nNu]/D");
    solTree2->Branch("nE"     ,&nE,    "nE[nNu]/D");

    solTree2->Branch("nMu"    ,&nMu,   "nMu/I");
    solTree2->Branch("mpx"    ,&mpx,   "mpx[nMu]/D");
    solTree2->Branch("mpy"    ,&mpy,   "mpy[nMu]/D");
    solTree2->Branch("mpz"    ,&mpz,   "mpz[nMu]/D");
    solTree2->Branch("mE"     ,&mE,    "mE[nMu]/D");
    solTree2->Branch("d0"     ,&d0,    "d0[nMu]/D");
    solTree2->Branch("X2"     ,&X2,    "X2[nMu]/D");
    solTree2->Branch("calE"   ,&calE,  "calE[nMu]/D");
    solTree2->Branch("mIso"   ,&mIso,  "mIso[nMu]/D");
    solTree2->Branch("nHits"  ,&nHits, "nHits[nMu]/I");

    solTree2->Branch("pvZ"    ,&pvZ,   "pvZ/D");
    solTree2->Branch("pvRho"  ,&pvRho, "pvRho/D");
    solTree2->Branch("pvNDF"  ,&pvNDF, "pvNDF/D");

    solTree2->Branch("ePt"    ,&ePt,   "ePt/D");
    solTree2->Branch("eEta"   ,&eEta,  "eEta/D");
    solTree2->Branch("eIso"   ,&eIso,  "eIso/D");
    solTree2->Branch("eHovE"  ,&eHovE, "eHovE/D");
    solTree2->Branch("eEovP"  ,&eEovP, "eEovP/D");

 } 

 /// Destructor
 virtual ~SolNtp2() {
    delete solTree2;
 }

 void FillB1( int eventId, vector<ttCandidate> jp4V, vector<LorentzVector> np4V, vector<ttCandidate> mp4V, vector<double> pv, vector<ttCandidate> vetoInfo ) {
      evtId   = eventId ;

      nJ      = jp4V.size() ;
      for (int i = 0; i < nJ; i++) {
          n90[i] = jp4V[i].nHits ;
          bTh[i] = jp4V[i].cuts[0] ;
          emF[i] = jp4V[i].cuts[1] ;
          EvH[i] = jp4V[i].cuts[2] ;
          jpx[i] = jp4V[i].p4.Px() ;
          jpy[i] = jp4V[i].p4.Py() ;
          jpz[i] = jp4V[i].p4.Pz() ;
          jE[i]  = jp4V[i].p4.E()  ;
      }

      nNu      = np4V.size() ;
      for (int i = 0; i < nNu; i++) {
          npx[i] = np4V[i].Px() ;
          npy[i] = np4V[i].Py() ;
          npz[i] = np4V[i].Pz() ;
          nE[i]  = np4V[i].E()  ;
      }

      int nGoodMu  = mp4V.size() ;
      for (int i = 0; i < nGoodMu; i++) {
          mpx[i] = mp4V[i].p4.Px() ;
          mpy[i] = mp4V[i].p4.Py() ;
          mpz[i] = mp4V[i].p4.Pz() ;
          mE[i]  = mp4V[i].p4.E()  ;
          mIso[i]  = mp4V[i].iso   ;
          nHits[i] = mp4V[i].nHits ;
          d0[i]    = mp4V[i].cuts[0] ;
          X2[i]    = mp4V[i].cuts[1] ;
          calE[i]  = mp4V[i].cuts[2] ;
      }

      pvZ   = pv[0] ;
      pvRho = pv[1] ;
      pvNDF = pv[2] ;

      int m = 0 ;
      for (size_t i =0; i < vetoInfo.size() ; i++ ) {
          if ( vetoInfo[i].pdgId == 11 ) { 
             ePt   = vetoInfo[i].p4.Pt() ;
	     eEta  = vetoInfo[i].eta ;
	     eIso  = vetoInfo[i].iso ;
	     eEovP = vetoInfo[i].cuts[0] ;
	     eHovE = vetoInfo[i].cuts[1] ;
          }
          if ( vetoInfo[i].pdgId == 13 ) {
             if ( m > 3 ) continue ;
             mpx[ nGoodMu + m ] = vetoInfo[i].p4.Px() ;
	     mpy[ nGoodMu + m ] = vetoInfo[i].p4.Py() ;
	     mpz[ nGoodMu + m ] = vetoInfo[i].p4.Pz() ;
	     mE[ nGoodMu + m ]  = vetoInfo[i].p4.E()  ;
	     mIso[ nGoodMu + m ]  = vetoInfo[i].iso   ;
	     nHits[ nGoodMu + m ] = vetoInfo[i].nHits ;
	     d0[ nGoodMu + m ]    = vetoInfo[i].cuts[0] ;
	     X2[ nGoodMu + m ]    = vetoInfo[i].cuts[1] ;
	     calE[ nGoodMu + m ]  = vetoInfo[i].cuts[2] ;
             m++ ;
          }
      }
      nMu = nGoodMu + m ;

      solTree2->Fill();
 }


 void Write() {
      solTree2->Write();
 }

 TTree *solTree2;

private:

     static const int maxCount  = 15 ;
     static const int maxCount1 = 5 ;

     int evtId;
     int nJ ;
     double jpx[maxCount];
     double jpy[maxCount];
     double jpz[maxCount];
     double jE[maxCount];
     double bTh[maxCount];
     double EvH[maxCount];
     double emF[maxCount];
     int    n90[maxCount];

     int nNu;
     double npx[2];
     double npy[2];
     double npz[2];
     double nE[2];

     int nMu;
     double mpx[maxCount1];
     double mpy[maxCount1];
     double mpz[maxCount1];
     double mE[maxCount1];
     double d0[maxCount1];
     double X2[maxCount1];
     double calE[maxCount1];
     double mIso[maxCount1];
     int    nHits[maxCount1];

     double pvZ ;
     double pvRho ;
     double pvNDF ;

     double ePt  ;
     double eEta ;
     double eIso ;
     double eEovP ;
     double eHovE ;

};
 
// new solution tree
class mcNtp {
public:

 mcNtp( TString treeName ) {

    mcTree = new TTree( treeName ,"");
 
    mcTree->Branch("evtId"   ,&evtId,   "evtId/I");
    mcTree->Branch("nParton" ,&nParton, "nParton/I");
    mcTree->Branch("pId"     ,&pId,     "pId[nParton]/I"); // particle ID or matched ID
    mcTree->Branch("px"      ,&px,      "px[nParton]/D");
    mcTree->Branch("py"      ,&py,      "py[nParton]/D");
    mcTree->Branch("pz"      ,&pz,      "pz[nParton]/D");
    mcTree->Branch("E"       ,&E,       "E[nParton]/D");

 } 

 /// Destructor
 virtual ~mcNtp() {
    delete mcTree;
 }

 
 void FillB( int eventId, int matchId[], vector<ttCandidate> jp4V, LorentzVector np4V, LorentzVector mp4V ) {

      evtId   = eventId ;
      nParton = jp4V.size() + 2;
      for (int i = 0; i < 4; i++) {
          pId[i] = matchId[i] ;
          px[i] = jp4V[ matchId[i] ].p4.Px() ;
          py[i] = jp4V[ matchId[i] ].p4.Py() ;
          pz[i] = jp4V[ matchId[i] ].p4.Pz() ;
          E[i]  = jp4V[ matchId[i] ].p4.E()  ;
      }
      pId[4]= matchId[4];
      px[4] = np4V.Px() ;
      py[4] = np4V.Py() ;
      pz[4] = np4V.Pz() ;
      E[4]  = np4V.E() ; 

      pId[5]= matchId[5];
      px[5] = mp4V.Px() ;
      py[5] = mp4V.Py() ;
      pz[5] = mp4V.Pz() ; 
      E[5]  = mp4V.E()  ;

      mcTree->Fill();
 }

 void FillB1( int eventId, vector<const reco::Candidate*> jp4V ) {

      evtId   = eventId ;
      nParton = jp4V.size() ;
      for (size_t i = 0; i < jp4V.size() ; i++) {
          pId[i] = jp4V[i]->pdgId() ;
          px[i]  = jp4V[i]->px() ;
          py[i]  = jp4V[i]->py() ;
          pz[i]  = jp4V[i]->pz() ;
          E[i]   = jp4V[i]->energy()  ;
      }

      mcTree->Fill();
 }


 void Write() {
      mcTree->Write();
 }

 TTree *mcTree;

private:

     static const int maxCount = 12;

     int evtId;
     int nParton;
     int pId[maxCount];
     double px[maxCount];
     double py[maxCount];
     double pz[maxCount];
     double E[maxCount];

};




class ObjNtp {
public:

 ObjNtp( TString treeName ) {

    objTree = new TTree( treeName ,"");
 
    objTree->Branch("evtId" ,&evtId,  "evtId/I");
    objTree->Branch("objId" ,&objId,  "objId/I");
    objTree->Branch("qCut"  ,&qCut,   "qCut/D");
    objTree->Branch("px"    ,&px,     "px/D");
    objTree->Branch("py"    ,&py,     "py/D");
    objTree->Branch("pz"    ,&pz,     "pz/D");
    objTree->Branch("E"     ,&E,      "E/D");
    objTree->Branch("pt"    ,&pt,     "pt/D");

 } 

 /// Destructor
 virtual ~ObjNtp() {
    delete objTree;
 }

 
 void FillB( int eventId, int objectId, double qCutVal, double _px, double _py, double _pz, double _E, double _pt ) {
      evtId = eventId ;
      objId = objectId;
      qCut  = qCutVal ;
      px    = _px     ;
      py    = _py     ;
      pz    = _pz     ;
      E     = _E      ;
      pt    = _pt     ;
      objTree->Fill();
 }


 void Write() {
      objTree->Write();
 }

 TTree *objTree;

private:

     int evtId;
     int objId;
     double qCut;
     double px;      
     double py;      
     double pz;      
     double E;      
     double pt;      

};


class SolNtp {
public:

 SolNtp( TString treeName ) {

    solTree = new TTree( treeName ,"");
 
    solTree->Branch("evtId"   ,&evtId,    "evtId/I");
    solTree->Branch("JId"     ,&JId,      "JId[5]/I");
    solTree->Branch("iniId"   ,&iniId,    "iniId/I");
    solTree->Branch("entSz"   ,&entSz,    "entSz/I");
    solTree->Branch("entId"   ,&entId,    "entId/I");
    solTree->Branch("nJ"      ,&nJ,       "nJ/I");
    solTree->Branch("probW"   ,&probW,    "probW/D");
    solTree->Branch("probM"   ,&probM,    "probM/D");
    solTree->Branch("probTt"  ,&probTt,   "probTt/D");
    solTree->Branch("hadTM"   ,&hadTM,    "hadTM/D");
    solTree->Branch("lepTM"   ,&lepTM,    "lepTM/D");
    solTree->Branch("hadWM"   ,&hadWM,    "hadWM/D");
    solTree->Branch("lepWM"   ,&lepWM,    "lepWM/D");

 } 

 /// Destructor
 virtual ~SolNtp() {
    delete solTree;
 }

 
 void FillB( int eventId, int wj1, int wj2, int bjh, int bjl, int ni, int _iniId, int _entSz, int _entId, int _nJ, double _probW, double _probM, double _probTt, double _hadTM, double _lepTM, double _hadWM, double _lepWM ) {
      evtId   = eventId ;
      JId[0]  = wj1 ;
      JId[1]  = wj2 ;
      JId[2]  = bjh ;
      JId[3]  = bjl ;
      JId[4]  = ni  ;
      iniId   = _iniId ;
      entSz   = _entSz ;
      entId   = _entId ;
      nJ      = _nJ ;
      probW   = _probW  ;
      probM   = _probM  ;
      probTt  = _probTt ;
      hadTM   = _hadTM ;
      lepTM   = _lepTM ;
      hadWM   = _hadWM ;
      lepWM   = _lepWM ;
      solTree->Fill();
 }


 void Write() {
      solTree->Write();
 }

 TTree *solTree;

private:

     int evtId;
     int JId[5];
     int iniId;
     int entSz;
     int entId;
     int nJ;
     double probW;
     double probM;
     double probTt;
     double hadTM;      
     double lepTM;      
     double hadWM;      
     double lepWM;      

};

struct tNtuple {

    SolNtp2  *muJets ; // Tree hold the solutions
    mcNtp    *mcmTree ;
    mcNtp    *genTree ;
};


#endif
