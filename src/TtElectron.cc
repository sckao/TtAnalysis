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
//TtElectron::TtElectron(const edm::ParameterSet& iConfig)
TtElectron::TtElectron( const edm::ParameterSet& iConfig )
{
  eleSetup = iConfig.getParameter<std::vector<double> >("eleSetup");

}


TtElectron::~TtElectron()
{
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------
void TtElectron::ElectronTreeFeeder(Handle<std::vector<pat::Electron> > patEle, ObjNtp* eTree, int eventId ) {

 int i =0 ;
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {
     i++;
     /*
     double caloE = it->caloEnergy() ;
     double HovE = it->hadronicOverEm() ;
     double emE   = caloE / ( 1. + HovE ) ;
     double hdE   = emE * HovE  ;
     */
     eTree->FillB( eventId, i, -1, it->px(), it->py(), it->pz(), it->energy(),  it->pt() );
 }

}

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

std::vector<ttCandidate> TtElectron::IsoEleSelection1( Handle<std::vector<pat::Electron> > patEle) {

 std::vector<ttCandidate> isoEle;
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

     ttCandidate ttEle ;
     ttEle.p4 = it->p4() ;
     ttEle.eta = it->eta() ;
     ttEle.iso =  IsoValue ;
     ttEle.cuts[0] = EovP ;
     ttEle.cuts[1] = HovE ;
     ttEle.cuts[2] = it->caloEnergy()  ;
     ttEle.charge = it->charge() ;
     ttEle.nHits = 0 ;
     ttEle.pdgId = it->pdgId() ;
     isoEle.push_back( ttEle );

 }

 return isoEle ;

}

