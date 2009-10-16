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
#include "FWCore/ParameterSet/interface/InputTag.h"


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
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtTools.h"
#include "TtObjHisto.h"
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
    void JetdRAnalysis(std::vector<const reco::Candidate*>  patJet, HOBJ1* histo1);

    void MuonAndJet( std::vector<const reco::Candidate*> patJet, const reco::Candidate* isoMu, HOBJ1* histo1 );

    void JetAndLepW( std::vector<const reco::Candidate*> patJet,  LorentzVector p1, HTOP1* histo1 );

    void matchedWJetsAnalysis(std::vector<jmatch> mcwjets, std::vector<const reco::Candidate*> isoMuons,
                              edm::Handle<std::vector<pat::Jet> > patJets, HTOP8* histo8);

    void matchedbJetsAnalysis(std::vector<jmatch> mcjets, std::vector<const reco::Candidate*> isoMuons,
                              edm::Handle<std::vector<pat::Jet> > patJets, HTOP7* histo7);

    std::vector< const reco::Candidate* > WJetSelection( edm::Handle<std::vector<pat::Jet> > jets, std::vector<const reco::Candidate*> IsoMuons, double JetEtCut );
    std::vector< const reco::Candidate* > bJetSelection( edm::Handle<std::vector<pat::Jet> > jets, std::vector<const reco::Candidate*> IsoMuons, double JetEtCut, string bTagAlgo );    

    double NofJetConstituents( pat::Jet theJet );
    double EoverH( pat::Jet theJet );
    //return W 4-momentum
    LorentzVector findW(LorentzVector qm1, LorentzVector qm2);

    void JetMatchedMuon( edm::Handle<std::vector<pat::Jet> > patJet , edm::Handle<std::vector<pat::Muon> > patMuon
             , const edm::Event& iEvent, const edm::EventSetup& iSetup, HTOP3* histo3, bool Done  );

    FreeTrajectoryState getFTS(GlobalPoint GP, GlobalVector GV, int charge,
                               const AlgebraicSymMatrix66& cov, const MagneticField* field);

    void bTagAnalysis( edm::Handle<std::vector<pat::Jet> > patJet, edm::Handle<std::vector<pat::Muon> > patMuon, HTOP7* histo7 );
    //void JetTrigger( edm::Handle<std::vector<pat::Jet> > jets, edm::Handle <edm::TriggerResults> triggers);

    template<typename jetT>
    std::vector<const reco::Candidate*> JetSelection( edm::Handle<std::vector<jetT> > patJet, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold, double fScale = 1., HOBJ1* histo1 = NULL, std::vector<bool>* bTags = NULL, string bTagAlgo = "", std::vector<double>* bDisList = NULL );
    std::vector<const reco::Candidate*> JetSelection( edm::Handle<std::vector<reco::CaloJet> > Jpts, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold, double fScale = 1., HOBJ1* histo1 = NULL, std::vector<bool>* bTags = NULL );

    template<typename jetT>
    std::vector< const reco::Candidate* > SoftJetSelection( edm::Handle<std::vector<jetT> > recoJet, std::vector<const reco::Candidate*> IsoMuons, double EtThreshold, double fScale = 1., std::vector<bool>* bTags = NULL, string bTagAlgo="", HOBJ1* histo1 = NULL, std::vector<double>* bDisList = NULL );

    template<typename jetT>
    void JetEtSpectrum( std::vector<jetT> patJet, HOBJ1* histo1);

   private:

    TtTools*       tools;
    TtMuon*        ttMuon;
    TtMCMatching*  JetMatching; 

    edm::InputTag caloSrc;
    edm::InputTag muonSrc;
    edm::InputTag genSrc;
    edm::InputTag jptSrc;
    double bCut;
    string bTagAlgo;
    //edm::ParameterSet tkParas; 
    //TrackAssociatorParameters theParameters; 

};

// Template Functions 
template<typename jetT>
std::vector<const reco::Candidate* > TtJet::JetSelection( edm::Handle<std::vector<jetT> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold, double fScale, HOBJ1* histo, std::vector<bool>* bTags, string bTagAlgo1, std::vector<double>* bDisList ) {

   int nBJets = 0;
   std::vector<const reco::Candidate* > jet_temp ;
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( fabs(j1->eta()) > 2.7 ) continue;  
       if ( ( fScale * j1->et() ) < EtThreshold ) continue;

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+
                     j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       double bDis = j1->bDiscriminator( bTagAlgo1 ) ;
       if ( bDis >= bCut ) nBJets++ ;
       if ( bTags != NULL ) {
          if ( bDis >= bCut ) (*bTags).push_back( true );
          if ( bDis <  bCut ) (*bTags).push_back( false );
          if ( bDisList != NULL ) (*bDisList).push_back( bDis );
       }

       reco::Candidate* scaleJ = j1->clone() ;
       LorentzVector jP4 = scaleJ->p4() ;
       scaleJ->setP4( jP4*fScale ) ;
       jet_temp.push_back( scaleJ );
       //jet_temp.push_back( &*j1 );

   }

   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtSorting<const reco::Candidate*>() );

   if ( histo != NULL ) {
      int jetSize = static_cast<int>( jets->size() );
      if (jetSize > 20) jetSize = 20 ;
      histo->Fill_1g( jetSize, jet_temp.size(), nBJets );
   }
   return jet_temp;

}

template<typename jetT>
std::vector<const reco::Candidate* > TtJet::SoftJetSelection( edm::Handle<std::vector<jetT> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold, double fScale, std::vector<bool>* bTags, string bTagAlgo, HOBJ1* histo, std::vector<double>* bDisList) {

   std::vector<const reco::Candidate* > jet_temp ;
   for (typename std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( ( fScale * j1->et() ) >= EtThreshold ) continue; 
       if ( fabs(j1->eta()) > 2.7 ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ 
                     j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;
       // bTags has to be "NULL" if use non-PAT Jets
       if ( bTags != NULL ) {
          double bDis = j1->bDiscriminator( bTagAlgo ) ;
          if ( bDis >= bCut ) (*bTags).push_back( true );
          if ( bDis <  bCut ) (*bTags).push_back( false );
          if ( bDisList != NULL ) (*bDisList).push_back( bDis );
       }

       reco::Candidate* scaleJ = j1->clone() ;
       LorentzVector jP4 = scaleJ->p4() ;
       scaleJ->setP4( jP4*fScale ) ;

       jet_temp.push_back( scaleJ );
   }

   //if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtDecreasing2 );
   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtSorting<const reco::Candidate*>() );
   if ( histo != NULL ) {
      double leadingJetEt = ( jet_temp.size() > 0 ) ? jet_temp[0]->et() : 0.  ;
      histo->Fill_1i( jet_temp.size(), leadingJetEt );
   }
   return jet_temp;
}

template<typename jetT>
void TtJet::JetEtSpectrum( std::vector<jetT> jet_temp, HOBJ1* histo1){

   histo1->Fill_1f( jet_temp.size() );
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

#endif
