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
// $Id$
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

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//

class TtEvtSelector {
   public:
    /// Constructor
    explicit TtEvtSelector();
    /// Destructor
    ~TtEvtSelector();

    /// Perform the real analysis
    bool eventSelection(edm::Handle<std::vector<pat::Muon> > rMu, edm::Handle<std::vector<pat::Electron> > rE,
                        edm::Handle<std::vector<pat::Jet> > rJet );
    int MCEvtSelection( edm::Handle<std::vector<reco::GenParticle> > genParticles );

    std::vector< const pat::Jet* > WJetSelection( edm::Handle<std::vector<pat::Jet> >  Jets );
    std::vector< const pat::Jet* > bJetSelection( edm::Handle<std::vector<pat::Jet> >  Jets );

    bool HLTSemiSelection( edm::Handle <edm::TriggerResults> triggers, int setup );
    void TriggerStudy( edm::Handle <edm::TriggerResults> triggers, int topo, int setup, HTOP9* histo9 );

    //std::vector<pat::Muon> MuonSelection( edm::Handle<std::vector<pat::Muon> >  Muons );
   private:

};

#endif
