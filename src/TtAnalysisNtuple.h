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
#include <string>
#include <iostream>
//#include "TtFormat.h"
typedef math::XYZTLorentzVector LorentzVector;


using namespace std;

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
    solTree2->Branch("jpt"   ,&jpt,    "jpt[nJ]/D");
    solTree2->Branch("bTh"   ,&bTh,    "bTh[nJ]/D");

    solTree2->Branch("nNu"    ,&nNu,     "nNu/I");
    solTree2->Branch("npx"    ,&npx,     "npx[nNu]/D");
    solTree2->Branch("npy"    ,&npy,     "npy[nNu]/D");
    solTree2->Branch("npz"    ,&npz,     "npz[nNu]/D");
    solTree2->Branch("nE"     ,&nE,      "nE[nNu]/D");
    solTree2->Branch("npt"    ,&npt,     "npt[nNu]/D");

    solTree2->Branch("nMu"    ,&nMu,     "nMu/I");
    solTree2->Branch("mpx"    ,&mpx,     "mpx[nMu]/D");
    solTree2->Branch("mpy"    ,&mpy,     "mpy[nMu]/D");
    solTree2->Branch("mpz"    ,&mpz,     "mpz[nMu]/D");
    solTree2->Branch("mE"     ,&mE,      "mE[nMu]/D");
    solTree2->Branch("mpt"    ,&mpt,     "mpt[nMu]/D");
 } 

 /// Destructor
 virtual ~SolNtp2() {
    delete solTree2;
 }

 
 void FillB( int eventId, vector<double> bThV, vector<const reco::Candidate*> jp4V, vector<LorentzVector> np4V,
                          vector<const reco::Candidate*> mp4V ) {
      evtId   = eventId ;

      nJ      = jp4V.size() ;
      for (int i = 0; i < nJ; i++) {
          bTh[i] = bThV[i] ;
          jpx[i] = jp4V[i]->px() ;
          jpy[i] = jp4V[i]->py() ;
          jpz[i] = jp4V[i]->pz() ;
          jE[i]  = jp4V[i]->energy()  ;
          jpt[i] = jp4V[i]->pt() ;
      }

      nNu      = np4V.size() ;
      for (int i = 0; i < nNu; i++) {
          npx[i] = np4V[i].Px() ;
          npy[i] = np4V[i].Py() ;
          npz[i] = np4V[i].Pz() ;
          nE[i]  = np4V[i].E()  ;
          npt[i] = np4V[i].Pt() ;
      }

      nMu      = mp4V.size() ;
      for (int i = 0; i < nMu; i++) {
          mpx[i] = mp4V[i]->px() ;
          mpy[i] = mp4V[i]->py() ;
          mpz[i] = mp4V[i]->pz() ;
          mE[i]  = mp4V[i]->energy()  ;
          mpt[i] = mp4V[i]->pt() ;
      }

      solTree2->Fill();
 }


 void Write() {
      solTree2->Write();
 }

 TTree *solTree2;

private:

     static const int maxCount = 10 ;

     int evtId;
     int nJ;
     double bTh[maxCount];
     double jpx[maxCount];
     double jpy[maxCount];
     double jpz[maxCount];
     double jE[maxCount];
     double jpt[maxCount];

     int nNu;
     double npx[maxCount];
     double npy[maxCount];
     double npz[maxCount];
     double nE[maxCount];
     double npt[maxCount];

     int nMu;
     double mpx[maxCount];
     double mpy[maxCount];
     double mpz[maxCount];
     double mE[maxCount];
     double mpt[maxCount];
};

// deceased version
class TtNtp {
public:
 

 TtNtp() {

    topTree = new TTree("topTree"," Ttbar infomaton tree");
    genB     = topTree->Branch("gen"    ,&gen.evtId,   "evtId/I:objId/I:qCut/D:px:py:pz:E:pt"); 
    selJetB  = topTree->Branch("selJet" ,&selJet.evtId,"evtId/I:objId/I:qCut/D:px:py:pz:E:pt"); 
    selMuB   = topTree->Branch("selMu"  ,&selMu.evtId, "evtId/I:objId/I:qCut/D:px:py:pz:E:pt"); 
    solNeuB  = topTree->Branch("solNeu" ,&solNeu.evtId,"evtId/I:objId/I:qCut/D:px:py:pz:E:pt"); 
    solTtB   = topTree->Branch("solTt"  ,&solTt.evtId, "evtId/I:wj1:wj2:bjh:bjl:prob/D:hadTM:lepTM:hadWM:lepWM"); 
    mcmTtB   = topTree->Branch("mcmTt"  ,&mcmTt.evtId, "evtId/I:wj1:wj2:bjh:bjl:prob/D:hadTM:lepTM:hadWM:lepWM"); 

 } 

 /// Destructor
 virtual ~TtNtp() {
    delete topTree;
 }

 // Fill topTree
 void FillB_gen( int evtId, int objId, double qCut, double px, double py, double pz, double E, double pt) {
      gen.evtId  = evtId;
      gen.objId  = objId;
      gen.qCut   = qCut;
      gen.px     = px;
      gen.py     = py;
      gen.pz     = pz;
      gen.E      = E ;
      gen.pt     = pt;
      genB->Fill();
 }

 void FillB_selJet( int evtId, int objId, double qCut, double px, double py, double pz, double E, double pt) {
      selJet.evtId = evtId ;
      selJet.objId = objId ;
      selJet.qCut  = qCut  ;
      selJet.px    = px ;
      selJet.py    = py ;
      selJet.pz    = pz ;
      selJet.E     = E ;
      selJet.pt    = pt ;
      selJetB->Fill();
 }

 void FillB_selMu( int evtId, int objId, double qCut, double px, double py, double pz, double E, double pt) {
      selMu.evtId = evtId ;
      selMu.objId = objId ;
      selMu.qCut  = qCut  ;
      selMu.px    = px ;
      selMu.py    = py ;
      selMu.pz    = pz ;
      selMu.E     = E ;
      selMu.pt    = pt ;
      selMuB->Fill();
 }

 void FillB_solNeu( int evtId, int objId, double qCut, double px, double py, double pz, double E, double pt) {
      solNeu.evtId = evtId ;
      solNeu.objId = objId ;
      solNeu.qCut  = qCut  ;
      solNeu.px    = px ;
      solNeu.py    = py ;
      solNeu.pz    = pz ;
      solNeu.E     = E ;
      solNeu.pt    = pt ;
      solNeuB->Fill();
 }

 void FillB_solTt( int evtId, int wj1, int wj2, int bjh, int bjl, double prob, double hadTM, double lepTM, double hadWM, double lepWM ) {
      solTt.evtId = evtId ;
      solTt.wj1   = wj1 ;
      solTt.wj2   = wj2 ;
      solTt.bjh   = bjh ;
      solTt.bjl   = bjl ;
      solTt.prob  = prob  ;
      solTt.hadTM = hadTM ;
      solTt.lepTM = lepTM ;
      solTt.hadWM = hadWM ;
      solTt.lepWM = lepWM ;
      //solTtB->Fill();
 }

 void FillB_mcmTt( int evtId, int wj1, int wj2, int bjh, int bjl, double prob, double hadTM, double lepTM, double hadWM, double lepWM ) {
      mcmTt.evtId = evtId ;
      mcmTt.wj1   = wj1 ;
      mcmTt.wj2   = wj2 ;
      mcmTt.bjh   = bjh ;
      mcmTt.bjl   = bjl ;
      mcmTt.prob  = prob  ;
      mcmTt.hadTM = hadTM ;
      mcmTt.lepTM = lepTM ;
      mcmTt.hadWM = hadWM ;
      mcmTt.lepWM = lepWM ;
      mcmTtB->Fill();
 }
 
 void Fill() {
      topTree->Fill();
 }

 void Write() {
      topTree->Write();
 }

 TTree *topTree;
 TBranch *genB;
 TBranch *selJetB;
 TBranch *selMuB;
 TBranch *solNeuB;
 TBranch *solTtB;
 TBranch *mcmTtB;

private:

 struct ObjPar {
        int evtId;
        int objId;
        double qCut;
        double px;      
        double py;      
        double pz;      
        double E;      
        double pt;      
 } selJet,selMu,solNeu,gen;

 struct EvtSol {
        int evtId;
        int wj1;
        int wj2;
        int bjh;
        int bjl;
        double prob;
        double hadTM;
        double lepTM;
        double hadWM;
        double lepWM;
 } solTt, mcmTt;

};
#endif
