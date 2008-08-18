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
void TtElectron::ElectronAnalysis(Handle<std::vector<pat::Electron> > patEle, HTOP4* histo4, NJet* jtree, int eventId ) {

 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     const reco::IsoDeposit* caloE = it->ecalIsoDeposit(); 
     const reco::IsoDeposit* caloH = it->hcalIsoDeposit(); 
     double emE = caloE->candEnergy() ;
     double hdE = caloH->candEnergy() ;
     jtree->FillBpatE( eventId, it->eta(), it->phi(), emE, hdE, it->p(), it->pt() );
     histo4->Fill4a( it->pt(), it->eta() );
 }
}


//std::vector<pat::Electron> TtElectron::IsoEleSelection( Handle<std::vector<pat::Electron> > patEle, HTOP4* histo4 ) {
std::vector<const reco::Candidate*> TtElectron::IsoEleSelection( Handle<std::vector<pat::Electron> > patEle, HTOP4* histo4 ) {

 //std::vector<pat::Electron> isoEle;
 std::vector<const reco::Candidate*> isoEle;
 isoEle.clear();
 for (std::vector<pat::Electron>::const_iterator it = patEle->begin(); it!= patEle->end(); it++) {

     //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackerIsoDeposit();
     std::pair<double, int> emR3 = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> emR5 = ecalIso->depositAndCountWithin(0.5); 
     std::pair<double, int> hdR3 = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR5 = hcalIso->depositAndCountWithin(0.5); 
     std::pair<double, int> tkR3 = trackIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR5 = trackIso->depositAndCountWithin(0.5);
     double thetaCal = (ecalIso->direction()).theta();
     double thetaTrk = (trackIso->direction()).theta();
     double sumEtR3 = (emR3.first + hdR3.first)* sin(thetaCal) ;
     double sumEtR5 = (emR5.first + hdR5.first)* sin(thetaCal) ;
     double sumPtR3 = (tkR3.first )* sin(thetaTrk) ;
     double sumPtR5 = (tkR5.first )* sin(thetaTrk) ;
     int calR3_count = emR3.second + hdR3.second ;
     int calR5_count = emR5.second + hdR5.second ;
     double IsoValue = it->pt() / sumEtR3  ;
      
     histo4->Fill4b(emR3.first,  emR5.first,  hdR3.first,   hdR5.first,  
                    calR3_count, calR5_count, tkR3.second, tkR5.second, 
                    sumPtR3, sumPtR5, sumEtR3, sumEtR5,
                    it->caloIso(), it->ecalIso(), it->hcalIso(), IsoValue, it->pt() );
  
     // Isolation Cut
     //if ( sumPtR3 > 3. || sumEtR3 > 5. ) continue; 
     if (IsoValue < 0.92 ) continue;
     //if ( tkR3.second   > 1 ) continue;
     //if ( it->ecalIso() > 3 ) continue;
     //if ( it->hcalIso() > 1 ) continue;

     histo4->Fill4c(emR3.first,  emR5.first,  hdR3.first,   hdR5.first,  
                    calR3_count, calR5_count, tkR3.second, tkR5.second,  
                    sumPtR3, sumPtR5, sumEtR3, sumEtR5,
                    it->caloIso(), it->ecalIso(), it->hcalIso(),  IsoValue, it->pt() );

     isoEle.push_back( &*it );

 }
 return isoEle ;

}

