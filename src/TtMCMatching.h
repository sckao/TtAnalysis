#ifndef TtMCMatching_H
#define TtMCMatching_H
// -*- C++ -*-
//
// Package:    TtMCMatching
// Class:      TtMCMatching
// 
/**\class TtMCMatching TtMCMatching.cc PhysicsTools/TtAnalysis/src/TtMCMatching.cc

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

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtAnalysisNtuple.h"
#include "TtAnalysisHisto.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//
typedef math::XYZTLorentzVector LorentzVector;

struct jmatch {
       int MomIdx ;
       double res_P;
       LorentzVector sumP4 ;
       std::vector<pat::Jet> assJets ;
       pat::Jet leadingJet ;
       pat::Jet trueJet ;
       reco::Particle mom ;
       bool hasMatched; 
};

struct iTt {
    int pdgId;
    int momId;
    LorentzVector p4;
    GlobalVector  gv; 
    double pt;
};

class TtMCMatching {
   public:
    /// Constructor
    explicit TtMCMatching();
    /// Destructor
    ~TtMCMatching();

    /// Perform the real analysis
    void MCTreeFeeder(edm::Handle<std::vector<reco::GenParticle> > genParticles, NJet* jtree, int eventId);

    std::vector<jmatch> matchWJets(edm::Handle<std::vector<reco::GenParticle> > genParticles,
                                   edm::Handle<std::vector<pat::Jet> > jets, std::vector<pat::Jet> selectedWJets,
                                   HTOP8* histo8, bool fillhisto );
    std::vector<jmatch> matchbJets(edm::Handle<std::vector<reco::GenParticle> > genParticles,
                                   edm::Handle<std::vector<pat::Jet> > jets, std::vector<pat::Jet> selectedbJets,
                                   HTOP7* histo7, bool fillhisto );
    std::vector<const reco::Candidate*> matchMuon(edm::Handle<std::vector<reco::GenParticle> > genParticles,
                                        edm::Handle<std::vector<pat::Muon> > muons, 
                                        std::vector<const reco::Candidate*> isoMuons, HTOP3* histo3, bool fillhisto);

    std::vector<const reco::Candidate*> matchElectron(edm::Handle<std::vector<reco::GenParticle> > genParticles,
                                        edm::Handle<std::vector<pat::Electron> > electrons,
                                        std::vector<const reco::Candidate*> isoEle, HTOP4* histo4, bool fillhisto);

    bool matchingGeneral( LorentzVector p4_1, iTt ttObject, double& dR0, double& ptRes0 );

    bool matchingGeneral( LorentzVector aP4, LorentzVector bP4, double& dR0, double& ptRes0 );

    std::vector<iTt> TtObjects( edm::Handle<std::vector<reco::GenParticle> > genParticles ); 

    std::vector<reco::Particle> ttPartons( edm::Handle<std::vector<reco::GenParticle> > genParticles, int targetId ) ;


   private:


};

#endif
