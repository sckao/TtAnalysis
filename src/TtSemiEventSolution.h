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
// $Id: TtSemiEventSolution.h,v 1.15 2009/10/15 06:47:17 sckao Exp $
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
#include <FWCore/Utilities/interface/InputTag.h>


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
#include "TtAnalysisNtuple.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"
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
typedef std::vector<int> Idx;

class TtSemiEventSolution {
   public:

    /// Constructor
    explicit TtSemiEventSolution( const edm::ParameterSet& );
    /// Destructor
    ~TtSemiEventSolution();

    /// Perform the real analysis
    void RecordSolutions(const edm::Event & iEvent,   int topo, int evtId, int njets, SolNtp2* solTree );

    void MCBuildSemiTt(const edm::Event & iEvent, int topo, int evtId, tNtuple* ntuples );

    // 
    bool recoW( std::vector<ttCandidate>& wjets, std::vector<iReco>& wCandidate, HTOP10* histo10 = NULL );

    bool recoW( ttCandidate& lepton, LorentzVector metP4, std::vector<iReco>& wCandidate, HTOP10* hitso10 = NULL );
    bool recoW( ttCandidate& lepton, LorentzVector metP4, std::vector<iReco>& wCandidates, bool FoundWSolution, HTOP10* hitso10 = NULL );

   
   private:
      // ----------member data ---------------------------

    TtEvtSelector* evtSelected;
    TtMCMatching*  MCMatching;
    TtMuon*        ttMuon;
    TtElectron*    ttEle;
    TtMET*         ttMET;
    TtJet*         ttJet;
    TtTools*       tools;

    std::vector<ttCandidate> isoLep;
    std::vector<ttCandidate> selectedJets;
    std::vector<LorentzVector> solvedMetP4;
    std::vector<double> pvInfo ;
    std::vector<ttCandidate> vetoInfo ;

    std::vector<ttCandidate> mcLep;
    std::vector<ttCandidate> mcWJets ;
    std::vector<LorentzVector> mcbJets ;
   
    int ini_Id;
    int ent_Id;
    int ent_sz;
    int counter[9] ;
    int nGoodJet ;

    // variables for py config
    bool debug;
    std::vector<double> jetSetup;

    //std::string recoMuon;
    std::string trigTag;
    edm::InputTag muonSrc;
    edm::InputTag electronSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
    edm::InputTag genSrc;

};

#endif
