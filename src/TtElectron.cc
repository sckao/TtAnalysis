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
TtElectron::TtElectron()
{

}


TtElectron::~TtElectron()
{
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------
void TtElectron::ElectronTreeFeeder(Handle<std::vector<pat::Electron> > patEle, NJet* jtree, int eventId ) {

 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     double caloE = it->caloEnergy() ;
     double HovE = it->hadronicOverEm() ;
     double emE   = caloE / ( 1. + HovE ) ;
     double hdE   = emE * HovE  ;
     jtree->FillBpatE( eventId, it->eta(), it->phi(), emE, hdE, it->p(), it->pt() );
 }

}

void TtElectron::ElectronAnalysis(Handle<std::vector<pat::Electron> > patEle, HTOP4* histo4  ) {

 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {
     histo4->Fill4a( it->pt(), it->eta() );
 }

}

void TtElectron::matchedElectronAnalysis( std::vector<const reco::Candidate*>  matchedEle, HTOP4* histo4  ) {

 for (std::vector<const reco::Candidate*>::const_iterator it = matchedEle.begin(); it!= matchedEle.end(); it++) {
     const reco::Candidate* it0 = *it ; 
     const pat::Electron* it1 = dynamic_cast<const pat::Electron*>( it0 );

     const reco::IsoDeposit* ecalIso  = it1->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it1->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it1->trackerIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3); 
     
     double caloE = it1->caloEnergy() ;
     double HovE = it1->hadronicOverEm() ;
     double emE   = caloE / ( 1. + HovE ) ;
     double hdE   = emE * HovE  ;
     double EovP = caloE / (*it)->p() ;
 
     double IsoValue = emR.first / (*it)->p() ;

     //cout<<" Pt:"<< (*it)->pt() <<" emR5= "<<emR5.first <<"  candE= "<<ecalIso->candEnergy() ;
     //cout<<" caloE: "<< it1->caloEnergy()<<endl;
     //cout<<" Pt:"<< (*it)->pt() <<" isoE:"<<emR.first<<"/"<<emR.second;
     //cout<< "  isoH:"<<hdR.first<<"/"<<hdR.second <<endl;

     histo4->Fill4d( (*it)->pt(), (*it)->eta(), EovP, HovE, tkR.second, emR.second,
                      it1->trackIso(), it1->ecalIso(), it1->hcalIso(), IsoValue );

 }

}

std::vector<const reco::Candidate*> TtElectron::IsoEleSelection( Handle<std::vector<pat::Electron> > patEle, HTOP4* histo4 ) {

 //std::vector<pat::Electron> isoEle;
 std::vector<const reco::Candidate*> isoEle;
 isoEle.clear();
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackerIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3); 
     //double thetaCal = (ecalIso->direction()).theta();
     //double thetaTrk = (trackIso->direction()).theta();
     double sumE = emR.first + hdR.first ;
     int calR_count = emR.second + hdR.second ;

     double IsoValue = emR.first / it->p() ;
      
     histo4->Fill4b(it->pt(), it->eta(), emR.first, hdR.first, sumE, calR_count, tkR.second,
                    it->trackIso(), it->ecalIso(), it->hcalIso(), IsoValue );
  
     // Isolation Cut
     if ( tkR.second   > 3 ) continue;
     if ( it->trackIso() > 3 ) continue;
     if ( IsoValue > 0.12 && fabs(it->eta()) < 1.5 ) continue;
     if ( it->ecalIso()  > 3 && fabs(it->eta()) >= 1.5 ) continue;
     //if ( it->ecalIso()  > 0.5 && fabs(it->eta()) <= 1.4 ) continue;
     //if ( it->ecalIso()  > 5   && fabs(it->eta())  > 1.4 ) continue;
     //if ( it->trackIso() > 3   && fabs(it->eta())  > 1.4 ) continue;
     histo4->Fill4c(it->pt(), it->eta(), emR.first, hdR.first, sumE, calR_count, tkR.second,
                    it->trackIso(), it->ecalIso(), it->hcalIso(), IsoValue );


     isoEle.push_back( &*it );

 }
 return isoEle ;

}


