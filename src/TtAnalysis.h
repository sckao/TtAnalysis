#ifndef TtAnalysis_H
#define TtAnalysis_H
// -*- C++ -*-
//
// Package:    TtAnalysis
// Class:      TtAnalysis
// 
/**\class TtAnalysis TtAnalysis.cc PhysicsTools/TtAnalysis/src/TtAnalysis.cc

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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"

#include "DataFormats/MuonReco/interface/Muon.h" 
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "DataFormates/TrackReco/interface/Track.h"
//#include "DataFormates/TrackReco/interface/TrackFwd.h"

//#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
//#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"

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

// hfPos[0]:eta, hfPos[1]:phi, hfPos[2]:pt 
typedef std::vector<double> hfPos ;
/// Lorentz vector
typedef math::XYZTLorentzVector LorentzVector;

class TtAnalysis : public edm::EDAnalyzer {
   public:
    /// Constructor
  

    explicit TtAnalysis(const edm::ParameterSet&);
    /// Destructor
    ~TtAnalysis();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);

    std::vector<int> findGrandMa(int pdgId, double eta, double phi,
                                 edm::Handle<std::vector<reco::GenParticle> > genParticles);
    hfPos findDaughter(int dauId, int momId, double eta, double phi,
                                 edm::Handle<std::vector<reco::GenParticle> > genParticles);
    //return Top 4-momentum
    LorentzVector findTop(LorentzVector qm1, LorentzVector qm2);
 
    bool recoW( std::vector<pat::Jet> wjets, std::vector<LorentzVector>& wCandidate );
    bool recoW( std::vector<const reco::Candidate*> lepton, edm::Handle<std::vector<pat::MET> > met,
                std::vector<LorentzVector>& wCandidate );
    bool recoTop( std::vector<LorentzVector> wCandidate, std::vector<pat::Jet> bjets, std::vector<LorentzVector>& tCandidate );

    std::vector<LorentzVector> recoSemiLeptonicTtEvent(int topo, std::vector<pat::Jet> theWJets,
                  std::vector<pat::Jet> thebJets, std::vector<const reco::Candidate*> mcMuons,
                  std::vector<const reco::Candidate*> mcElectrons, edm::Handle<std::vector<pat::MET> > met );
 
    double getInvMass( LorentzVector lv );


   private:
      // ----------member data ---------------------------

    TtEvtSelector* evtSelected;
    TtMCMatching*  MCMatching;
    TtMuon*        ttMuon;
    TtElectron*    ttEle;
    TtMET*         ttMET;
    TtJet*         ttJet;

    // Histograms
    HTOP1 *h_Jet;
    HTOP2 *h_MET;
    HTOP3 *h_Muon;
    HTOP4 *h_Ele;
    HTOP5 *h_Gam;
    HTOP6 *h_MJet;
    HTOP7 *h_BJet;
    HTOP8 *h_WJet;
    HTOP9 *h_Top;
    NJet  *t_Jet;

    // The file which will store the histos
    TFile *theFile;

    bool exclude;

    // Switch for debug output
    bool debug;
    int evtIt;

    std::string rootFileName;
    std::string leptonFlavour;
    edm::InputTag muonSrc;
    std::string recoMuon;
    edm::InputTag electronSrc;
    edm::InputTag photonSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
    edm::InputTag jetObj;
    edm::InputTag genSrc;
    edm::InputTag caloSrc;


};

#endif
