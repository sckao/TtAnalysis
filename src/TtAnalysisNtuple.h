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
#include <string>
#include <iostream>
//#include "TtFormat.h"

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
