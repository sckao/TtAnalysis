#ifndef TtMET_H
#define TtMET_H
// -*- C++ -*-
//
// Package:    TtMET
// Class:      TtMET
// 
/**\class TtMET TtMET.cc PhysicsTools/TtAnalysis/src/TtMET.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtMET.h,v 1.13 2009/10/15 06:47:16 sckao Exp $
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
#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtMuon.h"
#include "TtFormat.h"
#include "TtTools.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//
class TtMuon;

class TtMET {
   public:
    /// Constructor
    explicit TtMET(const edm::ParameterSet&);
    /// Destructor
    ~TtMET();

    /// Perform the real analysis
    //void MetTreeFeeder(edm::Handle<std::vector<pat::MET> > patMet, TtNtp* jtree, int eventId ); 

    void metAnalysis(const edm::Event & iEvent, HTOP2* histo2);

    LorentzVector METfromNeutrino( edm::Handle<std::vector<reco::GenParticle> > genParticles );
    LorentzVector METfromNeutrino( const edm::Event & iEvent );

    LorentzVector CaloMET( edm::Handle<CaloTowerCollection> calotowers );
    LorentzVector CaloMET( const edm::Event & iEvent );

    LorentzVector METfromObjects( std::vector<const reco::Candidate*> theLep, std::vector<const reco::Candidate*> theJets );

    void MetAndMuon( edm::Handle<std::vector<pat::MET> > patMet, std::vector<const reco::Candidate*> isoMuons, HTOP2* histo2, int njets  );
    void MetAndJets( edm::Handle<std::vector<pat::MET> > met, std::vector<const reco::Candidate*> theJets, HTOP2* histo2 );

    void METandNeutrino( std::vector<const reco::Candidate*> theLep, std::vector<const reco::Candidate*> theJets,
                         edm::Handle<std::vector<pat::MET> > met,
                         edm::Handle<std::vector<reco::GenParticle> > genParticles, HOBJ2* histo );

   private:

    TtMuon*   ttMuon;
    TtTools*  tools;

    edm::InputTag caloSrc;
    edm::InputTag muonSrc;
    edm::InputTag recoMetSrc;
    edm::InputTag metSrc;
    edm::InputTag genSrc;


};

#endif
