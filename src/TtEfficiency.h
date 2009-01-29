#ifndef TtEfficiency_H
#define TtEfficiency_H
// -*- C++ -*-
//
// Package:    TtEfficiency
// Class:      TtEfficiency
// 
/**\class TtEfficiency TtEfficiency.cc PhysicsTools/TtAnalysis/src/TtEfficiency.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtEfficiency.h,v 1.5 2009/01/15 14:58:02 sckao Exp $
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

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>

#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"

//
// class decleration
//
typedef math::XYZTLorentzVector LorentzVector;

class TtEfficiency {
   public:
    /// Constructor
    explicit TtEfficiency();
    /// Destructor
    ~TtEfficiency();

    /// Perform the real analysis
    void EventEfficiency(int topo, bool pass, HTOP9* histo9 );
    void JetEfficiency(std::vector<const pat::Jet*> recojets, std::vector<const pat::Jet*> mcjets, HTOP9* histo9 );
    void IsoLeptonEfficiency(std::vector<const reco::Candidate*> isolep, std::vector<const reco::Candidate*> mclep, HTOP9* histo9 );

   private:

};

#endif
