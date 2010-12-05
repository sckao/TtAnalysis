#ifndef TtNtupleProd_H
#define TtNtupleProd_H
// -*- C++ -*-
//
// Package:    TtNtupleProd
// Class:      TtNtupleProd
// 
/**\class TtNtupleProd TtNtupleProd.cc PhysicsTools/TtNtupleProd/src/TtNtupleProd.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtNtupleProd.h,v 1.1 2009/10/15 06:47:17 sckao Exp $
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
//#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "FWCore/Framework/interface/TriggerNames.h"

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

#include "TtMCMatching.h"
#include "TtSemiEventSolution.h"
#include "TtFormat.h"
#include "TtAnalysisNtuple.h"

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
class TtMCMatching;
class TtSemiEventSolution;

class TtNtupleProd : public edm::EDAnalyzer {
   public:
    /// Constructor

    explicit TtNtupleProd(const edm::ParameterSet&);
    /// Destructor
    ~TtNtupleProd();

    /// Perform the real analysis
    void analyze(const edm::Event & iEvent, const edm::EventSetup& iSetup);

   private:
      // ----------member data ---------------------------

    TtMCMatching*         MCMatching;
    TtSemiEventSolution*  semiSol;

    // Histograms & Trees
    tNtuple ntuples ;
    // The file which will store the histos
    TFile *theFile;

    bool exclude;

    // Switch for debug output
    bool debug;
    bool isData;
    int  evtIt;
    int  nJets ;

    string rootFileName;
    edm::InputTag genSrc;
    //edm::InputTag triggerSrc;

};

#endif
