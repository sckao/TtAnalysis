#ifndef TtEvtSelector_H
#define TtEvtSelector_H
// -*- C++ -*-
//
// Package:    TtEvtSelector
// Class:      TtEvtSelector
// 
/**\class TtEvtSelector TtEvtSelector.cc PhysicsTools/TtAnalysis/src/TtEvtSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtEvtSelector.h,v 1.10 2009/03/07 14:22:25 sckao Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"


#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtAnalysisHisto.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"
#include "TtMCMatching.h"


#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
/*
class TtMuon;
class TtMET;
class TtElectron;
*/
class TtJet;

class TtEvtSelector {
   public:
    /// Constructor
    explicit TtEvtSelector(const edm::ParameterSet& );
    //explicit TtEvtSelector();
    /// Destructor
    ~TtEvtSelector();

    /// Perform the real analysis
    int eventSelection(edm::Handle<std::vector<pat::Muon> > rMu, edm::Handle<std::vector<pat::Electron> > rE,
                       edm::Handle<std::vector<pat::Jet> > rJet, double jetEtThreshold );

    int eventSelection(edm::Handle<std::vector<pat::Muon> > rMu, edm::Handle<std::vector<pat::Electron> > rE,
                       edm::Handle<std::vector<reco::CaloJet> > rJet, double jetEtThreshold );

    int eventSelection(int topo, double JetEtCut, std::vector<const reco::Candidate*>& isoLep,  std::vector<const reco::Candidate*>& selectedJets, LorentzVector& metp4, const edm::Event& iEvent, string MetType, string JetType, std::vector<bool>* bTags = NULL );

    int eventSelection( int topo, double JetEtCut, std::vector<const reco::Candidate*>& isoLep,  std::vector<const reco::Candidate*>& selectedWJets, std::vector<const reco::Candidate*>& selectedbJets, LorentzVector& metp4, const edm::Event& iEvent, string MetType );

    int MCEvtSelection( edm::Handle<std::vector<reco::GenParticle> > genParticles );

    bool HLTSemiSelection( edm::Handle <edm::TriggerResults> triggers, int setup );
    void TriggerStudy( edm::Handle <edm::TriggerResults> triggers, int topo, int setup, HTOP9* histo9 );

   private:

   TtMuon*       ttMuon;
   TtElectron*   ttEle;
   TtJet*        ttJet;
   TtMET*        ttMET;
   TtMCMatching* mcMatch;

   edm::InputTag muonSrc;
   edm::InputTag electronSrc;
   edm::InputTag metSrc;
   edm::InputTag tcmetSrc;
   edm::InputTag jetSrc;
   edm::InputTag jptSrc;
   string bTagAlgo;

};

#endif
