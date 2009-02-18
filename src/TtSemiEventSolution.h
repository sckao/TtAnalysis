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
// $Id: TtSemiEventSolution.h,v 1.4 2009/01/23 16:08:16 sckao Exp $
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

#include "TtAnalysisHisto.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"
#include "TtEfficiency.h"
#include "TtFormat.h"
#include "TtTools.h"

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

class TtSemiEventSolution {
   public:

    /// Constructor
    explicit TtSemiEventSolution( const edm::ParameterSet& );
    /// Destructor
    ~TtSemiEventSolution();

    /// Perform the real analysis
    void BuildSemiTt(const edm::Event & iEvent, const edm::EventSetup& iSetup, int topo, tHisto histos );
    void MCBuildSemiTt(const edm::Event & iEvent, const edm::EventSetup& iSetup, int topo, tHisto histos );

    bool recoW( std::vector<const pat::Jet*> wjets, std::vector<iReco>& wCandidate  );
    bool recoW( std::vector<const pat::Jet*> wjets, std::vector<iReco>& wCandidate, HTOP10* histo10 );

    bool recoW( std::vector<const reco::Candidate*> lepton, edm::Handle<std::vector<pat::MET> > met,
                std::vector<iReco>& wCandidate );
    bool recoW( std::vector<const reco::Candidate*> lepton, edm::Handle<std::vector<pat::MET> > met,
                std::vector<iReco>& wCandidates, bool FoundWSolution );
    bool recoW( std::vector<const reco::Candidate*> lepton, edm::Handle<std::vector<pat::MET> > met,
                std::vector<iReco>& wCandidate, HTOP10* hitso10 );

    bool recoTop( std::vector<iReco> wCandidate, std::vector<const pat::Jet*> bjets, std::vector<iReco>& tCandidate, bool btagging );
    bool recoTop( std::vector<iReco> wCandidate, std::vector<const pat::Jet*> bjets, std::vector<iReco>& tCandidate, bool btagging, HTOP11* histo11 );

    // with 2b tagging
    std::vector<iReco> recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theWJets,
                  std::vector<const pat::Jet*> thebJets, std::vector<const reco::Candidate*> theLep,
                  edm::Handle<std::vector<pat::MET> > met, tHisto histos  );
 
    // no b-tagging
    std::vector<iReco> recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theJets,
                  std::vector<const reco::Candidate*> theLep, edm::Handle<std::vector<pat::MET> > met, tHisto histos  );
 
    void dmSortRecoObjects( std::vector<iReco>& objCand );

    void accuracySemiTt( std::vector<iReco> ttMC, std::vector<iReco> ttReco, HTOP9* histo9 ) ;

    //void MCTruthCheckB( std::vector<iReco> mcTt, std::vector<iReco> rcTt, std::vector<const pat::Jet*> mcBJets,
    //                   std::vector<const pat::Jet*> rcBJets, HTOP9* histo9 ) ;
    void MCTruthCheck( std::vector<iReco> mcTt, std::vector<iReco> rcTt, std::vector<const pat::Jet*> mcWJets,
                       std::vector<const pat::Jet*> rcWJets, std::vector<const pat::Jet*> mcBJets,
                       std::vector<const pat::Jet*> rcBJets, HTOP6* histo6 ); 

    void McRecoCompare( int topo, int r,  tHisto histos );

    void KeepBuildInfo( bool isData );
  
   private:
      // ----------member data ---------------------------

    TtEvtSelector* evtSelected;
    TtMCMatching*  MCMatching;
    TtMuon*        ttMuon;
    TtElectron*    ttEle;
    TtMET*         ttMET;
    TtJet*         ttJet;
    TtEfficiency*  ttEff;
    TtTools*       tools;

    std::vector<iReco> semiTt;
    std::vector<const reco::Candidate*> isoLep;
    std::vector<const pat::Jet*> selectedJets;
    std::vector<const pat::Jet*> selectedbJets;
    std::vector<const pat::Jet*> selectedWJets;

    std::vector<iReco> semiMCTt; 
    std::vector<const reco::Candidate*> mcLep;
    std::vector<const pat::Jet*> mcWJets ;
    std::vector<const pat::Jet*> mcbJets ;
   
    std::vector<TtResult> AllTt ;

    bool exclude;

    // Switch for debug output
    bool debug;
    bool btag;

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
