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
// $Id: TtJet.h,v 1.14 2009/10/15 06:47:16 sckao Exp $
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
//#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
//#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"

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
#include <FWCore/Utilities/interface/InputTag.h>

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

// for the fucking JEC uncertainty
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TtTools.h"
#include "TtObjHisto.h"
#include "TtAnalysisHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtMuon.h"
#include "TtMCMatching.h"
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

//class TtMuon;
//class TtMCMatching;

template<typename T> 
struct EtSorting{
       bool operator()( const T & t1, const T & t2 ) const {
		return t1->et() > t2->et();
	}

};

class TtJet {
   public:
    /// Constructor
    explicit TtJet(const edm::ParameterSet&);
    /// Destructor
    ~TtJet();

    /// Perform the real analysis
    //void JetTreeFeeder(edm::Handle<std::vector<pat::Jet> > patJet, TtNtp* jtree, int eventId );

    void jetAnalysis(edm::Handle<std::vector<pat::Jet> > patJet, HTOP1* histo1);

    void MuonAndJet( std::vector<const reco::Candidate*> patJet, const reco::Candidate* isoMu, HOBJ1* histo1 );
    void ElectronAndJet( std::vector<const reco::Candidate*> patJet, const reco::Candidate* isoEl, HOBJ1* histo1 );

    double NofJetConstituents( pat::Jet theJet );
    double EoverH( pat::Jet theJet );
    //return W 4-momentum
    LorentzVector findW(LorentzVector qm1, LorentzVector qm2);

    void JetMatchedMuon( edm::Handle<std::vector<pat::Jet> > patJet , edm::Handle<std::vector<pat::Muon> > patMuon
             , const edm::Event& iEvent, const edm::EventSetup& iSetup, HTOP3* histo3, bool Done  );

    FreeTrajectoryState getFTS(GlobalPoint GP, GlobalVector GV, int charge,
                               const AlgebraicSymMatrix66& cov, const MagneticField* field);

    void bTagAnalysis( const edm::Event& iEvent, edm::Handle<std::vector<pat::Jet> > patJet, double etaCut, double EtCut, HBJet* histo );

    //void JetTrigger( edm::Handle<std::vector<pat::Jet> > jets, edm::Handle <edm::TriggerResults> triggers);

    std::vector<const reco::Candidate*> JetSelection( edm::Handle<std::vector<pat::Jet> > patJet, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold, double fScale = 1., HOBJ1* histo1 = NULL, string bTagAlgo = "", std::vector<double>* bDisList = NULL );
    std::vector<const reco::Candidate*> JetSelection( edm::Handle<std::vector<reco::CaloJet> > caloJets, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold, double fScale = 1., HOBJ1* histo1 = NULL );

    std::vector< const reco::Candidate* > SoftJetSelection( edm::Handle<std::vector<pat::Jet> > recoJet, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold, double fScale = 1., string bTagAlgo="", HOBJ1* histo1 = NULL, std::vector<double>* bDisList = NULL );

    std::vector<ttCandidate> JetSelection1( edm::Handle<std::vector<pat::Jet> > patJet, std::vector<ttCandidate> IsoMuons, double EtThreshold, double fScale = 1., HOBJ1* histo1 = NULL, string bTagAlgo = "" );
    std::vector<ttCandidate> SoftJetSelection1( edm::Handle<std::vector<pat::Jet> > patJet, std::vector<ttCandidate> IsoMuons, double EtThreshold, double fScale = 1., HOBJ1* histo1 = NULL, string bTagAlgo = "" );

    double JES_Uncertainty( double pt, double eta ) ;

    template<typename jetT>
    void JetEtSpectrum( std::vector<jetT> theJet, HOBJ1* histo1);
    template<typename jetT>
    void JetdRAnalysis( std::vector<jetT> theJet, HOBJ1* histo1);

   private:

    TtTools*       tools;

    //edm::InputTag caloSrc;
    edm::InputTag muonSrc;
    edm::InputTag genSrc;
    double bCut;
    string bTagAlgo;
    std::vector<double> jetSetup;
    bool isData ;
    //edm::ParameterSet tkParas; 
    //TrackAssociatorParameters theParameters; 

};

// Template Functions
template<typename jetT>
void TtJet::JetEtSpectrum( std::vector<jetT> jet_temp, HOBJ1* histo1){

   histo1->Fill_1f( jet_temp.size() );
   if ( jet_temp.size() > 0) histo1->Fill_1d(jet_temp[0]->et() ) ;
   if ( jet_temp.size() > 3) {
      // get m3
      std::vector<LorentzVector> vlist;
      vlist.push_back( jet_temp[0]->p4() );
      vlist.push_back( jet_temp[1]->p4() );
      vlist.push_back( jet_temp[2]->p4() );
      double m3 = tools->getInvMass( vlist );
   
      histo1->Fill_1c(jet_temp[0]->et(), jet_temp[1]->et(), jet_temp[2]->et(), jet_temp[3]->et(), m3 ) ;
   }
}

// dR(j,j)
template<typename jetT>
void TtJet::JetdRAnalysis( std::vector<jetT>  theJet, HOBJ1* histo1 ){

    int jetSize = static_cast<int>( theJet.size() ) ;
    double dR = 99. ;
    for(int i=0; i< jetSize; i++) {

       if ( i+1 >= jetSize  ) continue;
       for (int j=i+1; j< jetSize; j++) {

           double dR1 = tools->getdR( theJet[i]->p4(), theJet[j]->p4() );
           if ( dR1 < dR ) dR = dR1;
       }
    }
    if ( dR < 98 ) histo1->Fill_1a( dR );
}


#endif
