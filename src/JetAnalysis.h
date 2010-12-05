#ifndef JetAnalysis_H
#define JetAnalysis_H
// -*- C++ -*-
//
// Package:    JetAnalysis
// Class:      JetAnalysis
// 
/**\class JetAnalysis JetAnalysis.cc PhysicsTools/TtAnalysis/src/JetAnalysis.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: JetAnalysis.h,v 1.6 2009/10/15 06:47:16 sckao Exp $
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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

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
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "DataFormates/TrackReco/interface/Track.h"
//#include "DataFormates/TrackReco/interface/TrackFwd.h"

//#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
//#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"

#include "TtObjHisto.h"
#include "TtEvtSelector.h"
#include "TtMuon.h"
#include "TtJet.h"
#include "TtMET.h"
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
/*
class TtEvtSelector;
class TtMuon;
class TtMET;
class TtJet;
*/
class JetAnalysis : public edm::EDAnalyzer {
   public:
    /// Constructor
    explicit JetAnalysis(const edm::ParameterSet&);
    /// Destructor
    ~JetAnalysis();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);



   private:
      // ----------member data ---------------------------

    TtEvtSelector*       evtSelected;
    TtMuon*              ttMuon;
    TtJet*               ttJet;
    TtMET*               ttMET;
    TtElectron*          ttEle;
 
    // Histograms
    HOBJ1 *hJ_Et20;
    HOBJ1 *hJ_Et25;
    HOBJ1 *hJ_Et30;
    HBJet *hb_Et20;
    HBJet *hb_Et30;
    HOBJ2 *hMET_J20;
    HOBJ2 *hMET_J25;
    HOBJ2 *hMET_J30;

    // The file which will store the histos
    TFile *theFile;


    // Switch for debug output
    bool debug;
    bool isData;
    int evtIt;
    std::vector<double> jetSetup;

    std::string rootFileName;
    edm::InputTag muonSrc;
    std::string recoMuon;
    edm::InputTag electronSrc;
    edm::InputTag jetSrc;
    edm::InputTag metSrc;
    edm::InputTag genJetSrc;
    edm::InputTag genSrc;
    //edm::InputTag caloSrc;
    string bTagAlgo;

};

#endif
