// -*- C++ -*-
//
// Package:    TtEvtSelector
// Class:      TtEvtSelector
// 
/**\class TtEvtSelector TtEvtSelector.cc PhysicsTools/TtAnalysis/src/TtEvtSelector.cc

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
#include "TtEvtSelector.h"
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
//TtEvtSelector::TtEvtSelector(const edm::ParameterSet& iConfig)
TtEvtSelector::TtEvtSelector()
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


TtEvtSelector::~TtEvtSelector()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtEvtSelector Analysis] Destructor called" << endl;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;
static bool PtDecreasing(const pat::Jet s1, const pat::Jet s2) { return ( s1.pt() > s2.pt() ); }

// ------------ method called to for each event  ------------
bool TtEvtSelector::eventSelection(Handle<std::vector<pat::Muon> > rMu, Handle<std::vector<pat::Electron> > rE,
                                Handle<std::vector<pat::Jet> >   rJet) {

 bool pass = false; 

 int nMu = 0;
 for (std::vector<pat::Muon>::const_iterator it = rMu->begin(); it != rMu->end(); it++ ) {
     if ( (*it).pt() > 15.0 )  nMu++ ; 
 }
  
 int nE = 0;
 for (std::vector<pat::Electron>::const_iterator it = rE->begin(); it != rE->end(); it++ ) {
     if ( (*it).pt() > 15.0 )  nE++ ; 
 }
  
 int nJet = 0;
 for (std::vector<pat::Jet>::const_iterator it = rJet->begin(); it != rJet->end(); it++ ) {
     if ( (*it).et() > 20.0 )  nJet++ ;
 }
 
 int nLep = nMu + nE ;

 if ( nLep > 0 && nJet > 3 ) pass = true;

 return pass;

}

int TtEvtSelector::MCEvtSelection( Handle<std::vector<reco::GenParticle> > genParticles ) {

   // type 2: di-lep, 1:semi-lep, 0:hadron , -1:Non-Tt event
   int type = -1;

   int tt[2] = {0} ;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       if ( (*it).pdgId() ==  6 && (*it).status() == 3 ) tt[0]++;
       if ( (*it).pdgId() == -6 && (*it).status() == 3 ) tt[1]++;
       /*  if ( abs((*it).pdgId()) == 6 ) {
          cout<<" The Top status:" << (*it).status() << endl;
       }*/ 
   }
   //cout<<" - - - - - - - - - - - - - "<<endl;

   if ( tt[0] == 1 && tt[1] == 1 ) {

      int nlep = 0;
      for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
          // find the W from t
          bool Wfromt = false ;
          if ( abs((*it).pdgId()) != 24) continue;
          for (size_t q=0; q< (*it).numberOfMothers(); ++q) {
              const reco::Candidate *mom = (*it).mother(q) ;
              if ( abs(mom->pdgId()) == 6 ) Wfromt = true;
              //cout <<" mom("<< mom->pdgId()<<")     status:"<<mom->status()<<endl;
              //if ( abs(mom->pdgId()) == 6 && mom->status() == 1) continue;
          }
    
          // tag the decay mode of W
          if ( !Wfromt ) continue;
          for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
              const reco::Candidate *dau = (*it).daughter(q) ;
              if ( abs(dau->pdgId()) == 11 || abs(dau->pdgId()) == 13 ) nlep++;
	      //cout<<" daughter("<< dau->pdgId()<<") status:"<< dau->status()<<endl;
          }
      }
      if (nlep > 2 )  type = 3;
      if (nlep == 2 ) type = 2;
      if (nlep == 1 ) type = 1;
      if (nlep == 0 ) type = 0;
   }

   return type;
}

std::vector<pat::Jet> TtEvtSelector::WJetSelection( Handle<std::vector<pat::Jet> >  Jets ) {

   std::vector<pat::Jet> jCollection;
   jCollection.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = Jets->begin(); j1 != Jets->end(); j1++)
   {
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;
       if (assTk.size() == 0) continue;
       if ((*j1).pt() < 20. ) continue;

       double bDis_TkCount = j1->bDiscriminator("trackCountingHighEffBJetTags") ;
       if (bDis_TkCount < 2.  ) continue;
       
       //double EovH1 = EoverH(*j1) ;
       //if ((*j1).nConstituents() < 5 || EovH1 > 20 || EovH1 < 0.01) continue;
       //if ((*j1).towersArea()/(*j1).pt() > 0.005 ) continue;

       jCollection.push_back( *j1 );
   }
   
   // sort the seeds by # of own segments
   sort(jCollection.begin(), jCollection.end(), PtDecreasing ) ;

   return jCollection;

}


std::vector<pat::Jet> TtEvtSelector::bJetSelection( Handle<std::vector<pat::Jet> >  Jets ) {

   std::vector<pat::Jet> jCollection;
   jCollection.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = Jets->begin(); j1 != Jets->end(); j1++)
   {
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;
       if (assTk.size() == 0) continue;
       if ((*j1).pt() < 20. ) continue;
       double bDis_TkCount = j1->bDiscriminator("trackCountingHighEffBJetTags") ;
       if (bDis_TkCount < 2.  ) continue;

       //double EovH1 = EoverH(*j1) ;
       //if ((*j1).nConstituents() < 5 || EovH1 > 20 || EovH1 < 0.01) continue;
       //if ((*j1).towersArea()/(*j1).pt() > 0.005 ) continue;

       jCollection.push_back( *j1 );
   }

   // sort the seeds by # of own segments
   sort(jCollection.begin(), jCollection.end(), PtDecreasing ) ;

   return jCollection;
}



std::vector<pat::Muon> TtEvtSelector::MuonSelection( Handle<std::vector<pat::Muon> >  Muons ) {

   std::vector<pat::Muon> uCollection;
   uCollection.clear();
   for (std::vector<pat::Muon>::const_iterator u1 = Muons->begin(); u1 != Muons->end(); u1++)
   {
       if( ! u1->isIsolationValid()  ) continue;
       if( ! u1->isTrackerMuon() ) continue;
   
       uCollection.push_back( *u1 );
       //cout<<" collected muons  "<< uCollection.size()<<endl;
   }
   return uCollection;
}


//define this as a plug-in
