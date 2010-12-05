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
// $Id: TtEvtSelector.h,v 1.13 2009/10/15 06:47:16 sckao Exp $
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
#include <FWCore/Utilities/interface/InputTag.h>


#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"
#include "TtMCMatching.h"
#include "TtFormat.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>

//
// class decleration

class TtEvtSelector {
   public:
    /// Constructor
    explicit TtEvtSelector(const edm::ParameterSet& );
    //explicit TtEvtSelector();
    /// Destructor
    ~TtEvtSelector();

    /// Perform the real analysis

    int eventSelection(int topo, std::vector<ttCandidate>& isoLep,  std::vector<ttCandidate>& selectedJets, std::vector<LorentzVector>& metp4, std::vector<ttCandidate>& vetoInfo, const edm::Event& iEvent, string MetType );

    int eventSelection( int topo, double JetEtCut, const edm::Event& iEvent, string MetType );

    int MCEvtSelection( edm::Handle<std::vector<reco::GenParticle> > genParticles );

    bool VertexSelection( const edm::Event& iEvent, std::vector<double>& pvInfo ) ;

    bool TriggerSelection( const edm::Event& iEvent, string trigPathName );

   private:

   TtMuon*       ttMuon;
   TtElectron*   ttEle;
   TtJet*        ttJet;
   TtMET*        ttMET;
   TtMCMatching* mcMatch;

   double pvZ ;
   double pvRho ;
   double pvNDF ;
   std::vector<double> jetSetup;
   edm::InputTag pvSrc;
   edm::InputTag beamSpotSrc;
   edm::InputTag muonSrc;
   edm::InputTag electronSrc;
   edm::InputTag metSrc;
   edm::InputTag jetSrc;
   edm::InputTag trigSrc;
   //edm::InputTag ecalRecHitSrc;
  
   edm::InputTag recoMetSrc;
   string bTagAlgo;

};

#endif
