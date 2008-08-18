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
//
//


// system include files
#include <memory>

// user include files
#include "TtMuon.h"

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
//TtMuon::TtMuon(const edm::ParameterSet& iConfig)
TtMuon::TtMuon()
{
   //now do what ever initialization is needed
  /*
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 

  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

   */

}


TtMuon::~TtMuon()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtMuon Analysis] Destructor called" << endl;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------
void TtMuon::muonAnalysis(Handle<std::vector<pat::Muon> > patMu, HTOP3* histo3, NJet* jtree, int eventId ) {

 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     const reco::IsoDeposit* caloE = it->ecalIsoDeposit(); 
     const reco::IsoDeposit* caloH = it->hcalIsoDeposit(); 
     double emE = caloE->candEnergy() ;
     double hdE = caloH->candEnergy() ;
     jtree->FillBpatMu( eventId, it->eta(), it->phi(), emE, hdE, it->p(), it->pt() );
     histo3->Fill3a( it->pt(), it->eta() );
 }

}

std::vector<const reco::Candidate*> TtMuon::IsoMuonSelection( Handle<std::vector<pat::Muon> > patMu, HTOP3* histo3 ) {

 //std::vector<pat::Muon> isoMuons;
 std::vector<const reco::Candidate*> isoMuons;
 isoMuons.clear();
 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

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
     double IsoValue = it->pt() / ( it->pt() + sumEtR3 + sumPtR3 ) ;
      
     histo3->Fill3b(emR3.first,  emR5.first,  hdR3.first,   hdR5.first,  
                    calR3_count, calR5_count, tkR3.second, tkR5.second, 
                    sumPtR3, sumPtR5, sumEtR3, sumEtR5,
                    it->caloIso(), it->ecalIso(), it->hcalIso(), IsoValue, it->pt() );
  
     // Isolation Cut
     //if ( sumPtR3 > 3. || sumEtR3 > 5. ) continue; 
     if (IsoValue < 0.92 ) continue;
     //if ( tkR3.second   > 1 ) continue;
     //if ( it->ecalIso() > 3 ) continue;
     //if ( it->hcalIso() > 1 ) continue;

     histo3->Fill3c(emR3.first,  emR5.first,  hdR3.first,   hdR5.first,  
                    calR3_count, calR5_count, tkR3.second, tkR5.second,  
                    sumPtR3, sumPtR5, sumEtR3, sumEtR5,
                    it->caloIso(), it->ecalIso(), it->hcalIso(),  IsoValue, it->pt() );

     isoMuons.push_back( &*it );
 }
 
 return isoMuons ;
}

std::vector<double> TtMuon::MuonEtCorrection( Handle<std::vector<pat::Muon> > mu ) {

     std::vector<double> ptcorr;
     ptcorr.clear();
     double pxy[4] ={0.0};
     for (std::vector<pat::Muon>::const_iterator u1 = mu->begin(); u1 != mu->end(); u1++)
     {
         if ( !(*u1).isGlobalMuon() ) continue;
	 const reco::IsoDeposit* caloE = u1->ecalIsoDeposit(); 
	 const reco::IsoDeposit* caloH = u1->hcalIsoDeposit(); 
	 double Uem = caloE->candEnergy() ;
	 double Uhd = caloH->candEnergy() ;

         double ex = (Uem+Uhd) * sin((*u1).theta()) * cos((*u1).phi()) ;
         double ey = (Uem+Uhd) * sin((*u1).theta()) * sin((*u1).phi()) ;
         pxy[0] += ((*u1).px() - ex) ;
         pxy[1] += ((*u1).py() - ey) ;
         pxy[2] += (*u1).px()  ;
         pxy[3] += (*u1).py()  ;
     }
     ptcorr.push_back(pxy[0]);
     ptcorr.push_back(pxy[1]);
     ptcorr.push_back(pxy[2]);
     ptcorr.push_back(pxy[3]);
     return ptcorr ;
}

