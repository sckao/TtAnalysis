// -*- C++ -*-
//
// Package:    TtPhoton
// Class:      TtPhoton
// 
/**\class TtPhoton TtPhoton.cc PhysicsTools/TtAnalysis/src/TtPhoton.cc

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
#include "TtPhoton.h"

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
//TtPhoton::TtPhoton(const edm::ParameterSet& iConfig)
TtPhoton::TtPhoton()
{

}


TtPhoton::~TtPhoton()
{
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------
void TtPhoton::PhotonTreeFeeder(Handle<std::vector<pat::Photon> > patGamma, ObjNtp* gTree, int eventId ) {

 int i=0 ;
 for (std::vector<pat::Photon>::const_iterator it = patGamma->begin(); it!= patGamma->end(); it++) {
     i++;
     /*
     const reco::IsoDeposit* caloE = it->ecalIsoDeposit(); 
     const reco::IsoDeposit* caloH = it->hcalIsoDeposit(); 
     double emE = caloE->candEnergy() ;
     double hdE = caloH->candEnergy() ;
     */
     gTree->FillB( eventId, i, -1., it->px(), it->py(), it->pz(), it->energy(), it->pt() );
 }

}

void TtPhoton::PhotonAnalysis(Handle<std::vector<pat::Photon> > patGamma, HTOP5* histo5  ) {

 for (std::vector<pat::Photon>::const_iterator it = patGamma->begin(); it!= patGamma->end(); it++) {
     if ( it->pt() < 20. ) continue;
     histo5->Fill5a( it->pt(), it->eta() );
 }

}

