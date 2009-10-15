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
static bool PtDecreasing(const reco::Candidate* s1, const reco::Candidate* s2) { return ( s1->pt() > s2->pt() ); }

// ------------ method called to for each event  ------------
/*
void TtMuon::MuonTreeFeeder(Handle<std::vector<pat::Muon> > patMu, TtNtp* jtree, int eventId ) {

 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     double emE = (it->calEnergy()).em ;
     double hdE = (it->calEnergy()).had ;

     jtree->FillBpatMu( eventId, it->eta(), it->phi(), emE, hdE, it->p(), it->pt() );
 }

}
*/

std::vector<const reco::Candidate*> TtMuon::IsoMuonSelection( Handle<std::vector<pat::Muon> > patMu, HOBJ3* histo1, HOBJ3* histo2 ) {

 int muSize = static_cast<int>( patMu->size() );
 if ( muSize > 20 ) muSize = 20;

 //std::vector<pat::Muon> isoMuons;
 std::vector<const reco::Candidate*> isoMuons;
 isoMuons.clear();
 int N_glbMu = 0;
 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackerIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3); 
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3);
     //double thetaCal = (ecalIso->direction()).theta();

     double sumE = emR.first + hdR.first ;
     double RelIso = it->pt() / ( it->pt() + emR.first + hdR.first + tkR.first ) ;
      
     // trust global muon only
     bool global =  it->isGlobalMuon() ;
     if ( !global ) continue;    
     N_glbMu++;

     reco::TrackRef glbTrack = it->globalTrack();
     reco::TrackRef inTrack = it->innerTrack();
     double dPt = (glbTrack->pt()/inTrack->pt()) - 1  ;

     histo1->Fill_3a( it->pt(), it->eta(), emR.first, hdR.first, tkR.first, sumE, RelIso );
     histo1->Fill_3b( glbTrack->pt(), inTrack->pt(), dPt, it->eta() );

     // Isolation Cut
     if ( RelIso < 0.9 ) continue;
     if ( fabs(it->eta()) > 2.1 ) continue;
     if ( it->pt() < 20. ) continue;

     // check the consistency of global pt and tracker pt
     histo2->Fill_3a( it->pt(), it->eta(), emR.first, hdR.first, tkR.first, sumE, RelIso );
     histo2->Fill_3b( glbTrack->pt(), inTrack->pt(), dPt, it->eta() );
     isoMuons.push_back( &*it );
 }

 if ( isoMuons.size() > 1 ) sort( isoMuons.begin(), isoMuons.end(), PtDecreasing );
 histo1->Fill_3c( muSize, N_glbMu );
 histo2->Fill_3c( muSize, isoMuons.size() );

 return isoMuons ;
}

std::vector<const reco::Candidate*> TtMuon::IsoMuonSelection( Handle<std::vector<pat::Muon> > patMu) {

 //std::vector<pat::Muon> isoMuons;
 std::vector<const reco::Candidate*> isoMuons;
 isoMuons.clear();
 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     bool isolated = IsoMuonID( *it, 0.9 );    

     if ( isolated && it->pt() > 20. && fabs(it->eta()) < 2.1 && it->isGlobalMuon() ) isoMuons.push_back( &*it );
     //if ( isolated && it->pt() > 20. && fabs(it->eta()) < 2.1 ) isoMuons.push_back( &*it );
 }
 if ( isoMuons.size() > 1 ) sort( isoMuons.begin(), isoMuons.end(), PtDecreasing );
 return isoMuons ;

}

std::vector<const reco::Candidate*> TtMuon::nonIsoMuonSelection( Handle<std::vector<pat::Muon> > patMu) {

 //std::vector<pat::Muon> isoMuons;
 std::vector<const reco::Candidate*> nonIsoMuons;
 nonIsoMuons.clear();
 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     bool isolated = IsoMuonID( *it, 0.9 );    

     //if ( isolated && it->pt() > 20. && fabs(it->eta()) < 2.1 && it->isGlobalMuon() ) isoMuons.push_back( &*it );
     if ( !isolated &&  fabs(it->eta()) < 2.1 ) nonIsoMuons.push_back( &*it );
 }
 if ( nonIsoMuons.size() > 1 ) sort( nonIsoMuons.begin(), nonIsoMuons.end(), PtDecreasing );
 return nonIsoMuons ;

}

bool TtMuon::IsoMuonID( pat::Muon muon, double isoCut ) {

    bool isolation = false;
    //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
    const reco::IsoDeposit* ecalIso  = muon.ecalIsoDeposit();
    const reco::IsoDeposit* hcalIso  = muon.hcalIsoDeposit();
    const reco::IsoDeposit* trackIso = muon.trackerIsoDeposit();
    std::pair<double, int> emR3 = ecalIso->depositAndCountWithin(0.3);
    std::pair<double, int> hdR3 = hcalIso->depositAndCountWithin(0.3);
    std::pair<double, int> tkR3 = trackIso->depositAndCountWithin(0.3);

    double RelIso = muon.pt() / ( muon.pt() + emR3.first + hdR3.first + tkR3.first ) ;

    // Isolation Cut
    if ( RelIso >= isoCut ) isolation = true;

    return isolation;

}

std::vector<double> TtMuon::MuonEtCorrection( Handle<std::vector<pat::Muon> > mu ) {

     std::vector<double> ptcorr;
     ptcorr.clear();
     double pxy[2] ={0.0};
     for (std::vector<pat::Muon>::const_iterator u1 = mu->begin(); u1 != mu->end(); u1++)
     {
         if ( !(*u1).isGlobalMuon() ) continue;
	 const reco::IsoDeposit* caloE = u1->ecalIsoDeposit(); 
	 const reco::IsoDeposit* caloH = u1->hcalIsoDeposit(); 
         // energy deposite in ECal and HCal
	 double Uem = caloE->candEnergy() ;
	 double Uhd = caloH->candEnergy() ;
         double ex = (Uem+Uhd) * sin((*u1).theta()) * cos((*u1).phi()) ;
         double ey = (Uem+Uhd) * sin((*u1).theta()) * sin((*u1).phi()) ;
         double et = sqrt( ex*ex + ey*ey );

         // over-correction
         bool badCorr =  ( et > (*u1).pt() ? true:false ) ;

         if ( ex > 5. || ey > 5. || badCorr ) { 
            pxy[0] += ((*u1).px() - ex) ;
            pxy[1] += ((*u1).py() - ey) ;
         }
     }
     ptcorr.push_back(pxy[0]);
     ptcorr.push_back(pxy[1]);
     return ptcorr ;
}

void TtMuon::MuonTrigger( Handle<std::vector<pat::Muon> >patMu, Handle <edm::TriggerResults> triggers ) {

   for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {
       std::vector<pat::TriggerPrimitive> trigInfo = it->triggerMatches() ;
       cout<<" === Muon trigger ==== "<<endl;
       for(size_t i=0; i< trigInfo.size(); i++) {
          cout<<"    ObjId:"<< trigInfo[i].triggerObjectId() ;
          cout<<"  ObjType:"<< trigInfo[i].triggerObjectType() <<endl;
          cout<<" filter name:"<<trigInfo[i].filterName() <<endl; 
       }
       /*
       if ( trigInfo.size() < 1) continue;
       edm::TriggerNames trigNames( *triggers );
       for (size_t i=0; i< triggers->size(); i++ ) {
           string triggered = triggers->accept(i) ? "Yes" : "No" ;
           cout<<" path("<<i<<") accepted ? "<< triggered ;
           cout<<" trigName: "<< trigNames.triggerName(i)<<endl;
       }
       */
   }

}

void TtMuon::matchedMuonAnalysis( std::vector<const reco::Candidate*>  matchedMuon, HOBJ3* histo ) {

 for (std::vector<const reco::Candidate*>::const_iterator it = matchedMuon.begin(); it!= matchedMuon.end(); it++) {

     const reco::Candidate* it0 = *it ;
     const pat::Muon* it1 = dynamic_cast<const pat::Muon*>( it0 );

     const reco::IsoDeposit* ecalIso  = it1->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it1->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it1->trackerIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3);
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3);
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3);

     double sumE = emR.first + hdR.first ;
     double RelIso = it1->pt() / ( it1->pt() + emR.first + hdR.first + tkR.first ) ;

     histo->Fill_3a( it1->pt(), it1->eta(), emR.first, hdR.first, tkR.first, sumE, RelIso );

 }

}
