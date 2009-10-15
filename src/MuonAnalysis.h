#ifndef MuonAnalysis_H
#define MuonAnalysis_H
// -*- C++ -*-
//
// Package:    MuonAnalysis
// Class:      MuonAnalysis
// 
/**\class MuonAnalysis MuonAnalysis.cc PhysicsTools/TtAnalysis/src/MuonAnalysis.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: MuonAnalysis.h,v 1.2 2009/03/07 14:22:25 sckao Exp $
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
#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

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
#include "TtObjHisto.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
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
//
class TtEvtSelector;
class TtElectron;
class TtMuon;
class TtJet;
class TtMCMatching;

class MuonAnalysis : public edm::EDAnalyzer {
   public:
    /// Constructor
    explicit MuonAnalysis(const edm::ParameterSet&);
    /// Destructor
    ~MuonAnalysis();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);



   private:
      // ----------member data ---------------------------

    TtEvtSelector*       evtSelected;
    TtMCMatching*        MCMatching;
    TtMuon*              ttMuon;
    TtJet*               ttJet;
    TtElectron*          ttEle;

    // Histograms
    HOBJ3 *hMu_Iso;
    HOBJ3 *hMu_All;
    HOBJ3 *hMu_MC;

    HOBJ4 *hEl_Iso;
    HOBJ4 *hEl_All;
    HOBJ4 *hEl_MC;

    // The file which will store the histos
    TFile *theFile;


    // Switch for debug output
    bool debug;
    int evtIt;
    double JEScale;

    std::string rootFileName;
    std::string recoMuon;
    edm::InputTag muonSrc;
    edm::InputTag electronSrc;
    edm::InputTag jetSrc;
    edm::InputTag metSrc;
    edm::InputTag genJetSrc;
    edm::InputTag jetObj;
    edm::InputTag genSrc;
    edm::InputTag caloSrc;


};

#endif
