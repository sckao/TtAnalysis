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
// $Id: TtAnalysis.h,v 1.14 2009/10/15 06:47:16 sckao Exp $
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
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "DataFormates/TrackReco/interface/Track.h"
//#include "DataFormates/TrackReco/interface/TrackFwd.h"

//#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
//#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

#include "TtObjHisto.h"
#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtPhoton.h"
#include "TtMET.h"
#include "TtJet.h"
#include "TtEfficiency.h"
#include "TtSemiEventSolution.h"
#include "TtFormat.h"

#include "TFile.h"
#include "TTree.h"
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
class TtPhoton;
class TtMET;
class TtJet;
class TtEfficiency; 
class TtSemiEventSolution;

class TtAnalysis : public edm::EDAnalyzer {
   public:
    /// Constructor

    explicit TtAnalysis(const edm::ParameterSet&);
    /// Destructor
    ~TtAnalysis();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);

   private:
      // ----------member data ---------------------------

    TtEvtSelector* evtSelected;
    TtMCMatching*  MCMatching;
    TtMuon*        ttMuon;
    TtElectron*    ttEle;
    TtPhoton*      ttGam;
    TtMET*         ttMET;
    TtJet*         ttJet;
    TtEfficiency*  ttEff;
    TtSemiEventSolution* semiSol;

    // Histograms & Trees
    tHisto histos ;

    // The file which will store the histos
    TFile *theFile;

    bool exclude;

    // Switch for debug output
    bool debug;
    bool trigOn;
    int evtIt;
    double JEScale;

    std::string rootFileName;
    std::string leptonFlavour;
    edm::InputTag muonSrc;
    std::string recoMuon;
    edm::InputTag electronSrc;
    edm::InputTag photonSrc;
    edm::InputTag metSrc;
    edm::InputTag jetSrc;
    edm::InputTag genJetSrc;
    edm::InputTag jetObj;
    edm::InputTag genSrc;
    edm::InputTag caloSrc;
    edm::InputTag triggerSrc;

};

#endif
