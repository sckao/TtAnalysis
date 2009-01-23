#ifndef TtJet_H
#define TtJet_H
// -*- C++ -*-
//
// Package:    TtJet
// Class:      TtJet
// 
/**\class TtJet TtJet.cc PhysicsTools/TtAnalysis/src/TtJet.cc

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

// For propagation
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Track Calo Mapping!
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"

//#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
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
#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"


#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtMuon.h"
#include "TtMCMatching.h"
#include "TtEvtSelector.h"
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
//typedef math::XYZTLorentzVector LorentzVector;

class TtMuon;
class TtMCMatching;
class TtEvtSelector;

class TtJet {
   public:
    /// Constructor
    explicit TtJet(const edm::ParameterSet&);
    /// Destructor
    ~TtJet();

    /// Perform the real analysis
    void JetTreeFeeder(edm::Handle<std::vector<pat::Jet> > patJet, NJet* jtree, int eventId );

    void jetAnalysis(edm::Handle<std::vector<pat::Jet> > patJet, HTOP1* histo1);
    void JetdRAnalysis(edm::Handle<std::vector<pat::Jet> > patJet, HTOP1* histo1);

    void thirdETJetSpectrum( edm::Handle<std::vector<reco::GenJet> > genJet, HTOP1* histo1);
    void thirdETJetSpectrum( std::vector<pat::Jet> patJet, HTOP1* histo1);

    void MuonAndJet( std::vector<pat::Jet> patJet, LorentzVector p1, HTOP1* histo1 );
    void JetAndLepW( std::vector<pat::Jet> patJet,  LorentzVector p1, HTOP1* histo1 );

    void genJetInfo(edm::Handle<std::vector<reco::GenJet> > genJet,
                    edm::Handle<std::vector<reco::GenParticle> > genParticles,
                    HTOP1* histo1, HTOP7* histo7, HTOP8* histo8 );

    void matchedWJetsAnalysis(std::vector<jmatch> mcwjets, std::vector<const reco::Candidate*> isoMuons, HTOP8* histo8);

    void matchedbJetsAnalysis(std::vector<jmatch> mcbjets, std::vector<jmatch> mcwjets, 
                              std::vector<const reco::Candidate*> isoMuons, HTOP7* histo7);

    void selectedWJetsAnalysis(edm::Handle<std::vector<pat::Jet> > patJet, HTOP8* histo8 );

    std::vector< const pat::Jet* > WJetSelection( edm::Handle<std::vector<pat::Jet> >  Jets );
    std::vector< const pat::Jet* > bJetSelection( edm::Handle<std::vector<pat::Jet> >  Jets );    

    double NofJetConstituents( pat::Jet theJet );
    double EoverH( pat::Jet theJet );
    //return W 4-momentum
    LorentzVector findW(LorentzVector qm1, LorentzVector qm2);

    void JetMatchedMuon( edm::Handle<std::vector<pat::Jet> > patJet , edm::Handle<std::vector<pat::Muon> > patMuon
             , const edm::Event& iEvent, const edm::EventSetup& iSetup, HTOP3* histo3 );
    void JetMatchedMuon( edm::Handle<std::vector<pat::Jet> > patJet , edm::Handle<std::vector<pat::Muon> > patMuon
             , const edm::Event& iEvent, const edm::EventSetup& iSetup, HTOP3* histo3, bool Done  );

    FreeTrajectoryState getFTS(GlobalPoint GP, GlobalVector GV, int charge,
                               const AlgebraicSymMatrix66& cov, const MagneticField* field);

    double getEta(double vx, double vy, double vz );
    double getdPhi(  LorentzVector v1, LorentzVector v2 );
    double getdR(  LorentzVector v1, LorentzVector v2 );
    double getY( LorentzVector v1 );
    double getInvMass( std::vector<LorentzVector> vlist );

    //std::vector<pat::Jet> JetSelection( edm::Handle<std::vector<pat::Jet> > patJet, LorentzVector muP4 );
    std::vector<pat::Jet> JetSelection( edm::Handle<std::vector<pat::Jet> > patJet, std::vector<const reco::Candidate*> IsoMuons );
    std::vector< const pat::Jet* > JetSelection( edm::Handle<std::vector<pat::Jet> > patJet, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold );

    void bTagAnalysis( edm::Handle<std::vector<pat::Jet> > patJet, HTOP7* histo7 );

    void JetTrigger( edm::Handle<std::vector<pat::Jet> > jets, edm::Handle <edm::TriggerResults> triggers);

   private:

    TtMuon*        fromTtMuon;
    TtMCMatching*  JetMatching; 

    edm::InputTag caloSrc;
    edm::InputTag muonSrc;
    edm::InputTag genSrc;
    edm::ParameterSet tkParas; 
    TrackAssociatorParameters theParameters; 

};

#endif
