#ifndef TtMuon_H
#define TtMuon_H
// -*- C++ -*-
//
// Package:    TtMuon
// Class:      TtMuon
// 
/**\class TtMuon TtMuon.cc PhysicsTools/TtAnalysis/src/TtMuon.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtMuon.h,v 1.10 2009/10/15 06:47:17 sckao Exp $
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
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h" 
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
//#include "DataFormats/MuonReco/interface/MuonEnergy.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TtAnalysisHisto.h"
#include "TtObjHisto.h"
#include "TtAnalysisNtuple.h"
#include "TtFormat.h"
#include "TtTools.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//

class TtMuon {
   public:
    /// Constructor
    explicit TtMuon(const edm::ParameterSet&);
    //explicit TtMuon();

    /// Destructor
    ~TtMuon();

    /// Perform the real analysis
    std::vector<const reco::Candidate*> IsoMuonSelection( edm::Handle<std::vector<pat::Muon> > patMu );

    std::vector<ttCandidate> IsoMuonSelection1( edm::Handle<std::vector<pat::Muon> > patMu, edm::Handle<std::vector<pat::Jet> > patJet, edm::Handle<reco::BeamSpot> bSpot_, double PVz, std::vector<ttCandidate>& vetoInfo );

    bool IsoMuonID( pat::Muon Mu, double isoCut );
    double IsoMuonID( pat::Muon Mu );

    std::vector<double> MuonEtCorrection( edm::Handle<std::vector<pat::Muon> > mu );

    //void MuonTrigger( edm::Handle<std::vector<pat::Muon> >patMu, edm::Handle <edm::TriggerResults> triggers);

    void matchedMuonAnalysis( std::vector<const reco::Candidate*>  matchedMuon, HOBJ3* histo  ); 

    void PatMuonScope( edm::Handle<std::vector<pat::Muon> > mu, const edm::Event& iEvent, HOBJ3* histo);
    void PromptMuonScope( edm::Handle<std::vector<pat::Muon> > patMu, HOBJ3* histo ) ;

    double PatMuonRelIso( const pat::Muon& patMu, double ConeSize, double& ecalIso, double& hcalIso, double& trkIso );

    template<typename muonT>
    void MuonScope( std::vector<muonT*>& theMu, const edm::Event& iEvent, HOBJ3* histo);

   private:
   std::vector<double> muSetup ;
   //edm::InputTag muonSrc;
   TtTools*       tools;

};

template<typename muonT>
void TtMuon::MuonScope( std::vector<muonT*>&  muons, const edm::Event& iEvent, HOBJ3* histo ){
//void TtMuon::MuonScope( edm::Handle<std::vector<muonT> > muons, const edm::Event& iEvent, HOBJ3* histo ){

     // get the beam spot information
     //edm::Handle<reco::BeamSpot> theBeamSpot;
     //iEvent.getByLabel("offlineBeamSpot", theBeamSpot );

     //math::XYZPoint pos_beamspot( 0, 0, 0 );
     //if ( theBeamSpot->isValid() )  pos_beamspot = pos_beamspot( beamSpot.x0(), beamSpot.y0(), beamSpot.z0() );
     //math::XYZPoint pos_beamspot( theBeamSpot->x0(), theBeamSpot->y0(), theBeamSpot->z0() );

     int N_glbMu ;
     //for ( typename std::vector<muonT*>::iterator m1 = muons.begin(); m1 != muons.end(); m1++ ) {
     for ( int i=0; i < muons.size(); i++ ) {

         //1. Muon Isolation
         reco::MuonIsolation isoR3 = muons[i]->isolationR03();
         float emIso  = isoR3.emVetoEt ;
         float hadIso = isoR3.hadVetoEt + isoR3.hoVetoEt;
         float trkIso = isoR3.trackerVetoPt ;
         double RelIso = ( emIso + hadIso + trkIso ) / muons[i]->pt() ; 

         histo->Fill_3d( muons[i]->pt(), muons[i]->eta(), RelIso, emIso, hadIso, trkIso );

         // track quality
         int nSeg = 0 ;
         double stdX2 = -1 ;
         if ( muons[i]->isStandAloneMuon() ) {
            reco::TrackRef outTrack = muons[i]->outerTrack();
            nSeg  = outTrack->numberOfValidHits();
            stdX2 = outTrack->normalizedChi2();
            histo->Fill_3e( nSeg, stdX2 );
         }

         int nHit = 0 ;
         double trkX2 = -1 ;
         if ( muons[i]->isTrackerMuon() ) {
            reco::TrackRef inTrack = muons[i]->innerTrack();
            nHit  = inTrack->numberOfValidHits();
            trkX2 = inTrack->normalizedChi2();
            // beam spot study , according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideFindingBeamSpot
            //double d0    = -1.* muons[i]->innerTrack()->dxy( pos_beamspot );
            //double beamWdthX = theBeamSpot->BeamWidthX() ;
            //double beamWdthY = theBeamSpot->BeamWidthY() ;
            //double d0Err = sqrt( inTrack->d0Error()*inTrack->d0Error() + 0.5*beamWdthX*beamWdthX + 0.5* beamWdthY*beamWdthY );
            double d0    = -1.* muons[i]->innerTrack()->d0();
            double d0Err = inTrack->d0Error()  ;
            double ipSig = d0 / d0Err ;
            histo->Fill_3f( nHit, trkX2, ipSig );
         }

         int aHit = 0 ;
         double glbX2 = -1 ;
         if ( muons[i]->isGlobalMuon() ) {
            reco::TrackRef inTrack = muons[i]->innerTrack();
            reco::TrackRef glbTrack = muons[i]->globalTrack();
            aHit  = glbTrack->numberOfValidHits();
            glbX2 = glbTrack->normalizedChi2();

            double Pt_Res = (glbTrack->pt()/inTrack->pt()) - 1  ;
            histo->Fill_3g( aHit, glbX2, glbTrack->pt(), muons[i]->eta(), RelIso, inTrack->pt(), Pt_Res );
            N_glbMu++;
         }

     }
     int muSize = static_cast<int>( muons.size() );
     histo->Fill_3c( muSize, N_glbMu );

}

#endif
