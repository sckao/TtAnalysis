#ifndef TtSemiEventSolution_H
#define TtSemiEventSolution_H
// -*- C++ -*-
//
// Package:    TtSemiEventSolution
// Class:      TtSemiEventSolution
// 
/**\class TtSemiEventSolution TtSemiEventSolution.cc PhysicsTools/TtAnalysis/src/TtSemiEventSolution.cc

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
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/JetReco/interface/GenJet.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "DataFormates/TrackReco/interface/Track.h"
//#include "DataFormates/TrackReco/interface/TrackFwd.h"

//#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
//#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

#include "TtAnalysisHisto.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"
#include "TtEfficiency.h"
#include "TtFormat.h"

#include "TFile.h"
#include "TVector3.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//
class TtEvtSelector;
class TtMCMatching;
class TtMuon;
class TtElectron;
class TtMET;
class TtJet;
class TtEfficiency; 

/// Lorentz vector
/*
typedef math::XYZTLorentzVector LorentzVector;
typedef std::pair<int, LorentzVector> iParton;

//
struct iReco{
    LorentzVector p4;
    std::pair<int,int> from;
    std::pair<const reco::Candidate*, const reco::Candidate*> ptr;
    //std::pair<LorentzVector, LorentzVector> q4;
    std::vector<iParton> q4v ;
    double dm;
    double mt; // only filled for leptonic W
};
*/

class TtSemiEventSolution {
   public:
    /// Constructor
  

    explicit TtSemiEventSolution(const edm::ParameterSet&);
    /// Destructor
    ~TtSemiEventSolution();

    /// Perform the real analysis
    void BuildSemiTt(const edm::Event & iEvent, const edm::EventSetup& iSetup, HTOP3* histo3, HTOP4* histo4,  HTOP7* histo7, HTOP8* histo8, HTOP9* histo9 );

    bool recoW( std::vector<const pat::Jet*> wjets, std::vector<iReco>& wCandidate  );
    bool recoW( std::vector<const reco::Candidate*> lepton, edm::Handle<std::vector<pat::MET> > met,
                std::vector<iReco>& wCandidate );
    bool recoTop( std::vector<iReco> wCandidate, std::vector<const pat::Jet*> bjets, std::vector<iReco>& tCandidate );

    std::vector<iReco> recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theWJets,
                  std::vector<const pat::Jet*> thebJets, std::vector<const reco::Candidate*> mcMuons,
                  std::vector<const reco::Candidate*> mcElectrons, edm::Handle<std::vector<pat::MET> > met, HTOP9* histo9  );
 
    double getInvMass( LorentzVector lv );

    void dmSortRecoObjects( std::vector<iReco>& objCand );

    void accuracySemiTt( std::vector<iReco> ttMC, std::vector<iReco> ttReco, HTOP9* histo9 ) ;

   private:
      // ----------member data ---------------------------

    TtEvtSelector* evtSelected;
    TtMCMatching*  MCMatching;
    TtMuon*        ttMuon;
    TtElectron*    ttEle;
    TtMET*         ttMET;
    TtJet*         ttJet;
    TtEfficiency*  ttEff;

    bool exclude;

    // Switch for debug output
    bool debug;

    edm::InputTag muonSrc;
    std::string recoMuon;
    edm::InputTag electronSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
    edm::InputTag genJetSrc;
    edm::InputTag jetObj;
    edm::InputTag genSrc;
    edm::InputTag caloSrc;


};

#endif
