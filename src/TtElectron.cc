// -*- C++ -*-
//
// Package:    TtElectron
// Class:      TtElectron
// 
/**\class TtElectron TtElectron.cc PhysicsTools/TtAnalysis/src/TtElectron.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
//
//


// system include files
#include <memory>

// user include files
#include "TtElectron.h"

#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"

//#include "FWCore/Framework/interface/Event.h"

//#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//


// constructors and destructor
using namespace edm;
using namespace std;
TtElectron::TtElectron( const edm::ParameterSet& iConfig )
{
  eleSetup     = iConfig.getParameter<std::vector<double> >("eleSetup");
  //electronSrc  = iConfig.getParameter<edm::InputTag> ("electronSource");

}


TtElectron::~TtElectron()
{
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------

void TtElectron::ElectronAnalysis(Handle<std::vector<pat::Electron> > patEle, HTOP4* histo4  ) {

 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {
     if ( it->pt() < 20. ) continue;
     histo4->Fill4a( it->pt(), it->eta() );
 }

}

void TtElectron::matchedElectronAnalysis( std::vector<const reco::Candidate*>  matchedEle, HOBJ4* histo4  ) {

 for (std::vector<const reco::Candidate*>::const_iterator it = matchedEle.begin(); it!= matchedEle.end(); it++) {
     const reco::Candidate* it0 = *it ; 
     const pat::Electron* it1 = dynamic_cast<const pat::Electron*>( it0 );

     const reco::IsoDeposit* ecalIso  = it1->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it1->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it1->trackIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3); 
     
     double caloE = it1->caloEnergy() ;
     double HovE = it1->hadronicOverEm() ;
     //double emE   = caloE / ( 1. + HovE ) ;
     //double hdE   = emE * HovE  ;
     double EovP = caloE / (*it)->p() ;
 
     double sumE = emR.first + hdR.first ;

     double emCompensation = ecalIso->depositWithin(0.055);
     //double sumIso = emR.first + hdR.first + tkR.first - emCompensation ;
     double sumIso = emR.first + hdR.first - emCompensation ;
     double IsoValue = it1->et() / (it1->et() + sumIso );

     histo4->Fill_4a( (*it)->pt(), (*it)->eta(), emR.first, hdR.first, tkR.first, sumE, IsoValue, HovE, EovP );

 }

}

/*
std::vector<const reco::Candidate*> TtElectron::IsoEleSelection( Handle<std::vector<pat::Electron> > patEle, HOBJ4* histo1, HOBJ4* histo2 ) {

 int elSize = static_cast<int>(  patEle->size() );
 if ( elSize > 20 ) elSize = 20;

 std::vector<const reco::Candidate*> isoEle;
 isoEle.clear();
 int N_goodEle = 0;
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3); 

     double sumE = emR.first + hdR.first ;
     int calR_count = emR.second + hdR.second ;

     double emCompensation = ecalIso->depositWithin(0.055);
     double sumIso = emR.first + hdR.first + tkR.first - emCompensation ;
     double IsoValue = it->et() / (it->et() + sumIso );
      
     double EovP = it->caloEnergy() / it->p() ;
     double HovE = it->hadronicOverEm() ;

     histo1->Fill_4a(it->pt(), it->eta(), emR.first, hdR.first, tkR.first, sumE, calR_count, tkR.second,
                    IsoValue, HovE, EovP );
     // quality cut  
     if ( it->pt() < eleSetup[0] )  continue ;
     if ( HovE > eleSetup[3] ) continue;
     if ( EovP < eleSetup[4] ) continue;
     N_goodEle++ ;

     // Isolation Cut
     if ( fabs( it->eta() ) > eleSetup[1] ) continue;    
     if ( IsoValue < eleSetup[2] ) continue;

     histo2->Fill_4a(it->pt(), it->eta(), emR.first, hdR.first, tkR.first, sumE, calR_count, tkR.second,
                    IsoValue, HovE, EovP );

     isoEle.push_back( &*it );

 }

 histo1->Fill_4b( elSize, N_goodEle ); 
 histo2->Fill_4b( elSize, isoEle.size() ); 

 return isoEle ;

}
*/
std::vector<const reco::Candidate*> TtElectron::IsoEleSelection( Handle<std::vector<pat::Electron> > patEle) {

 //std::vector<pat::Electron> isoEle;
 std::vector<const reco::Candidate*> isoEle;
 isoEle.clear();
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
     /*
     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3); 
     double emCompensation = ecalIso->depositWithin(0.055);
     double sumIso = emR.first + hdR.first + tkR.first - emCompensation;
     double IsoValue = it->et() / (it->et() + sumIso );
     */

     double emIso = it->dr04EcalRecHitSumEt() ;
     double hdIso = it->dr04HcalTowerSumEt();
     double tkIso = it->dr04TkSumPt() ;
     double IsoValue = (emIso + hdIso+ tkIso) / it->et() ;
     //double IsoValue =  it->et() / ( it->et() + emIso + hdIso );
      
     double EovP = it->caloEnergy() / it->p() ;
     double HovE = it->hadronicOverEm() ;

     // Isolation Cut
     if ( it->pt() < eleSetup[0] || fabs( it->eta() ) > eleSetup[1] )  continue ;
     if ( IsoValue > eleSetup[2] ) continue;     
     //if ( IsoValue < eleSetup[2] ) continue;     
     if ( HovE > eleSetup[3] ) continue;
     if ( EovP < eleSetup[4] || EovP > 1.2 ) continue;
     //if ( EovP < eleSetup[4]  ) continue;

     isoEle.push_back( &*it );

 }
 return isoEle ;

}

void TtElectron::PatEleScope( Handle<std::vector<pat::Electron> > patEle, HOBJ4* histo ) {

 int elSize = static_cast<int>(  patEle->size() );
 if ( elSize > 20 ) elSize = 20;

 int nIsoEl = 0 ;
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
     /*
     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3); 
     double emCompensation = ecalIso->depositWithin(0.055);
     double sumIso = emR.first + hdR.first + tkR.first - emCompensation;
     double IsoValue = it->et() / (it->et() + sumIso );
     */

     double emIso = it->dr04EcalRecHitSumEt() ;
     double hdIso = it->dr04HcalTowerSumEt();
     double tkIso = it->dr04TkSumPt() ;
     double IsoValue = (emIso + hdIso+ tkIso) / it->et() ;
     double sumE = emIso + hdIso ;
     double EovP = it->caloEnergy() / it->p() ;
     double HovE = it->hadronicOverEm() ;

     histo->Fill_4a( it->pt(), it->eta(), emIso, hdIso, tkIso, sumE, IsoValue, HovE, EovP );

     // Isolation Cut
     if ( it->pt() < eleSetup[0] || fabs( it->eta() ) > eleSetup[1] )  continue ;
     if ( IsoValue > eleSetup[2] ) continue;     
     if ( HovE > eleSetup[3] ) continue;
     if ( EovP < eleSetup[4] || EovP > 1.2 ) continue;
     nIsoEl++ ;
 }

 histo->Fill_4b( elSize, nIsoEl ); 

}

//std::vector<ttCandidate> TtElectron::IsoEleSelection1( Handle<std::vector<pat::Electron> > patEle, Handle<reco::BeamSpot> bSpot_, std::vector<ttCandidate>& vetoInfo, Handle<EcalRecHitCollection> recHits  ) {
std::vector<ttCandidate> TtElectron::IsoEleSelection1( Handle<std::vector<pat::Electron> > patEle, Handle<reco::BeamSpot> bSpot_, std::vector<ttCandidate>& vetoInfo  ) {

 //const  EcalRecHitCollection *myRecHits = recHits.product();

 std::vector<ttCandidate> isoEle;
 isoEle.clear();
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     // stupid ECal cleaning
     /*
     const reco::CaloClusterPtr    seed =   it->superCluster()->seed();
     const DetId seedId = seed->seed();
     EcalSeverityLevelAlgo severity;
     double myswissCross =  severity.swissCross(seedId, *myRecHits) ;
     if ( myswissCross > 0.95 ) continue ;
     */

     // Isolation Value
     double emIso = it->dr03EcalRecHitSumEt() ;
     double hdIso = it->dr03HcalTowerSumEt();
     double tkIso = it->dr03TkSumPt() ;
     double IsoValue = (emIso + hdIso+ tkIso) / it->et() ;
      
     double EovP = it->caloEnergy() / it->p() ;
     double HovE = it->hadronicOverEm() ;

     // beam Spot information
     reco::BeamSpot bSpot = *bSpot_ ;
     //double beamWdthX = bSpot.BeamWidthX() ;
     //double beamWdthY = bSpot.BeamWidthY() ;

     // d0 Cut , set pass = false in order to use the loose veto 
     bool pass = true ;
     reco::GsfTrackRef gsfTrk = it->gsfTrack() ;
     int nHit = gsfTrk->numberOfValidHits();
     //double d0 = it->dB() ;
     //double d0 = -1.*gsfTrk->dxy( bSpot.position() );
     //double d0Err = sqrt( gsfTrk->d0Error()*gsfTrk->d0Error() + 0.5*beamWdthX*beamWdthX + 0.5* beamWdthY*beamWdthY );

     if ( it->pt()          < eleSetup[0]  ) pass = false ;
     if ( fabs( it->eta() ) > eleSetup[1]  ) pass = false ;
     if ( IsoValue          > eleSetup[2]  ) pass = false ;     
     //if ( HovE > eleSetup[3]               ) pass = false ;
     //if ( EovP < eleSetup[4] || EovP > 1.2 ) pass = false ;
     //if ( d0   > 0.02                      ) pass = false ;
     /*
     if ( pass ) {
        ttCandidate ttEle ;
	ttEle.p4 = it->p4() ;
	ttEle.eta = it->eta() ;
	ttEle.iso =  IsoValue ;
	ttEle.cuts[0] = EovP ;
	ttEle.cuts[1] = HovE ;
	ttEle.cuts[2] = it->caloEnergy()  ;
	ttEle.charge = it->charge() ;
	ttEle.nHits = nHit ;
	ttEle.pdgId = it->pdgId() ;
	isoEle.push_back( ttEle );
     } else {
     */
     if ( pass ) {
        //if ( it->pt()          < 15. )  continue ;
	//if ( fabs( it->eta() ) > 2.5 )  continue ;
	//if ( IsoValue          > 0.2 )  continue ;
        ttCandidate vetoEle ;
	vetoEle.p4  = it->p4() ;
	vetoEle.eta = it->eta() ;
	vetoEle.iso =  IsoValue ;
	vetoEle.cuts[0] = EovP ;
	vetoEle.cuts[1] = HovE ;
	vetoEle.cuts[2] = it->caloEnergy()  ;
	vetoEle.charge = it->charge() ;
	vetoEle.nHits = nHit ;
	vetoEle.pdgId = 11 ;
	vetoInfo.push_back( vetoEle  ); 
     }

 }

 return isoEle ;

}

