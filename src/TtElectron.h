#ifndef TtElectron_H
#define TtElectron_H
// -*- C++ -*-
//
// Package:    TtElectron
// Class:      TtElectron
// 
/**\class TtElectron TtElectron.cc PhysicsTools/TtAnalysis/src/TtElectron.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtElectron.h,v 1.8 2009/03/07 14:22:25 sckao Exp $
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
typedef math::XYZTLorentzVector LorentzVector;

class TtElectron {
   public:
    /// Constructor
    explicit TtElectron();
    /// Destructor
    ~TtElectron();

    /// Perform the real analysis
    void ElectronTreeFeeder(edm::Handle<std::vector<pat::Electron> > patMu, ObjNtp* jtree, int eventId );
    void ElectronAnalysis(edm::Handle<std::vector<pat::Electron> > patMu, HTOP4* histo4 );

    void matchedElectronAnalysis( std::vector<const reco::Candidate*> patMu, HOBJ4* histo4 );

    std::vector<const reco::Candidate*> IsoEleSelection( edm::Handle<std::vector<pat::Electron> > patEle, 
                                                         HOBJ4* histo1, HOBJ4* histo2 );

    std::vector<const reco::Candidate*> IsoEleSelection( edm::Handle<std::vector<pat::Electron> > patEle );

    

   private:

};

#endif
