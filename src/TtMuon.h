#ifndef TtMuon_H
#define TtMuon_H
// -*- C++ -*-
//
// Package:    TtMuon
// Class:      TtMuon
// 
/**\class TtMuon TtMuon.cc PhysicsTools/TtAnalysis/src/TtMuon.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtMuon.h,v 1.9 2009/03/07 14:22:25 sckao Exp $
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
#include "TtObjHisto.h"
#include "TtAnalysisNtuple.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//

class TtMuon {
   public:
    /// Constructor
    explicit TtMuon();
    /// Destructor
    ~TtMuon();

    /// Perform the real analysis
    //void MuonTreeFeeder(edm::Handle<std::vector<pat::Muon> > patMu, TtNtp* jtree, int eventId );

    std::vector<const reco::Candidate*> IsoMuonSelection( edm::Handle<std::vector<pat::Muon> > patMu, HOBJ3* histo1,  HOBJ3* histo2 );
    std::vector<const reco::Candidate*> IsoMuonSelection( edm::Handle<std::vector<pat::Muon> > patMu );
    std::vector<const reco::Candidate*> nonIsoMuonSelection( edm::Handle<std::vector<pat::Muon> > patMu );

    bool IsoMuonID( pat::Muon Mu, double isoCut );

    std::vector<double> MuonEtCorrection( edm::Handle<std::vector<pat::Muon> > mu );

    void MuonTrigger( edm::Handle<std::vector<pat::Muon> >patMu, edm::Handle <edm::TriggerResults> triggers);

    void matchedMuonAnalysis( std::vector<const reco::Candidate*>  matchedMuon, HOBJ3* histo  ); 

   private:

};

#endif
