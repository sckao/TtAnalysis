#ifndef TtAnalysisNtuple_H
#define TtAnalysisNtuple_H

/** \class TtAnalysisNtuple
 *  Collection of Ntuples 
 *
 * Author: S.C. Kao  - UC Riverside
 */

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include <string>
#include <iostream>

using namespace std;

class NJet {
public:
 

 NJet() {
    jetT = new TTree("jetT"," jet infomaton tree");
    BpatJ  = jetT->Branch("patJ" ,&patJ.eventId,"eventId/I:eta/D:phi:caloE:caloH:p:pt"); 
    BpatMu = jetT->Branch("patMu",&patMu.eventId,"eventId/I:eta/D:phi:caloE:caloH:p:pt"); 
    BpatGa = jetT->Branch("patGa",&patGa.eventId,"eventId/I:eta/D:phi:caloE:caloH:p:pt"); 
    BpatNu = jetT->Branch("patNu",&patNu.eventId,"eventId/I:eta/D:phi:caloE:caloH:p:pt"); 
    BpatE  = jetT->Branch("patE" ,&patE.eventId, "eventId/I:eta/D:phi:caloE:caloH:p:pt"); 
    Bgen   = jetT->Branch("gen" ,&gen.eventId,"eventId/I:pdgId/I:eta/D:phi:energy:pt"); 
 } 

 /// Destructor
 virtual ~NJet() {
    delete jetT;
 }


 void FillBpatJ(int eventId, double eta, double phi, double caloE, double caloH, double p,double pt ) 
 {
      patJ.eventId = eventId;
      patJ.eta   = eta;
      patJ.phi   = phi;
      patJ.caloE = caloE;
      patJ.caloH = caloH;
      patJ.p     = p;
      patJ.pt    = pt;
      BpatJ->Fill();
 }
 void FillBpatE(int eventId, double eta, double phi, double caloE, double caloH, double p,double pt ) 
 {
      patE.eventId = eventId;
      patE.eta   = eta;
      patE.phi   = phi;
      patE.caloE = caloE;
      patE.caloH = caloH;
      patE.p     = p;
      patE.pt    = pt;
      BpatE->Fill();
 }
 void FillBpatGa(int eventId, double eta, double phi, double caloE, double caloH, double p,double pt ) 
 {
      patGa.eventId = eventId;
      patGa.eta   = eta;
      patGa.phi   = phi;
      patGa.caloE = caloE;
      patGa.caloH = caloH;
      patGa.p     = p;
      patGa.pt    = pt;
      BpatGa->Fill();
 }
 void FillBpatNu(int eventId, double eta, double phi, double caloE, double caloH, double p,double pt ) 
 {
      patNu.eventId = eventId;
      patNu.eta   = eta;
      patNu.phi   = phi;
      patNu.caloE = caloE;
      patNu.caloH = caloH;
      patNu.p     = p;
      patNu.pt    = pt;
      BpatNu->Fill();
 }
 void FillBpatMu(int eventId, double eta, double phi, double caloE, double caloH, double p,double pt ) 
 {
      patMu.eventId = eventId;
      patMu.eta   = eta;
      patMu.phi   = phi;
      patMu.caloE = caloE;
      patMu.caloH = caloH;
      patMu.p     = p;
      patMu.pt    = pt;
      BpatMu->Fill();
 }
 
 void FillBgen(int eventId, int pdgId, double eta, double phi, double energy, double pt ) 
 {
      gen.eventId = eventId;
      gen.pdgId  = pdgId;
      gen.eta    = eta;
      gen.phi    = phi;
      gen.energy = energy;
      gen.pt     = pt;
      Bgen->Fill();
 }

 void Write() {
      jetT->Write();
 }

 TTree *jetT;
 TBranch *BpatJ;
 TBranch *Bgen;
 TBranch *BpatMu;
 TBranch *BpatGa;
 TBranch *BpatNu;
 TBranch *BpatE;

private:

 struct TreeObj {
        int eventId;
        double eta;
        double phi;
        double caloE;
        double caloH;
        double p;
        double pt;
 } patJ,patMu,patE,patGa,patNu;

 struct TreeMC {
        int eventId;
        int pdgId;
        double eta;
        double phi;
        double energy;
        double pt;
 } gen;

};


#endif
