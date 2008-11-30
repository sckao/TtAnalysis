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
static bool PtDecreasing(const pat::Jet* s1, const pat::Jet* s2) { return ( s1->pt() > s2->pt() ); }
static bool PtDecreasing1(const pat::Jet s1, const pat::Jet s2) { return ( s1.pt() > s2.pt() ); }

// ------------ method called to for each event  ------------
bool TtEvtSelector::eventSelection(Handle<std::vector<pat::Muon> > rMu, Handle<std::vector<pat::Electron> > rE,
                                Handle<std::vector<pat::Jet> >   rJet) {

 bool pass = false; 

 int nMu = 0;
 for (std::vector<pat::Muon>::const_iterator it = rMu->begin(); it != rMu->end(); it++ ) {
     if ( (*it).pt() < 30.0 || fabs( it->eta() ) > 2.1 )  continue; 

     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackerIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3);
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3);
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3);

     double sumIso = emR.first + hdR.first + tkR.first ;
     double RelIso = it->pt() / (it->pt() + sumIso );

     if ( RelIso < 0.9 ) continue;

     nMu++;
 }
  
 int nE = 0;
 for (std::vector<pat::Electron>::const_iterator it = rE->begin(); it != rE->end(); it++ ) {
     if ( (*it).pt() < 30.0 || fabs( it->eta() ) > 2.4 )  continue ; 

     const reco::IsoDeposit* ecalIso  = it->ecalIsoDeposit();
     const reco::IsoDeposit* hcalIso  = it->hcalIsoDeposit();
     const reco::IsoDeposit* trackIso = it->trackerIsoDeposit();
     std::pair<double, int> emR = ecalIso->depositAndCountWithin(0.3);
     std::pair<double, int> hdR = hcalIso->depositAndCountWithin(0.3);
     std::pair<double, int> tkR = trackIso->depositAndCountWithin(0.3);

     double sumIso = emR.first + hdR.first + tkR.first ;
     double RelIso = it->pt() / (it->pt() + sumIso );

     if ( RelIso < 0.3 ) continue;

     nE++;
 }
  
 bool jetPreSelect = false;
 std::vector<pat::Jet> jet_temp ;
 for (std::vector<pat::Jet>::const_iterator it = rJet->begin(); it != rJet->end(); it++ ) {
     if ( (*it).pt()  < 40. ) continue;
     jet_temp.push_back( *it );
 }
 if ( jet_temp.size() > 0) { 
    sort(jet_temp.begin(), jet_temp.end(), PtDecreasing1 );
    if ( jet_temp[0].pt() > 60. && jet_temp.size() > 3  ) jetPreSelect = true ; 
 }
 int nLep = nMu + nE ;

 if ( nLep > 0 && jetPreSelect ) pass = true;

 return pass;

}

int TtEvtSelector::MCEvtSelection( Handle<std::vector<reco::GenParticle> > genParticles ) {

   // type 4:semi-electron 3:semi-Muon 2: di-lep, 1:semi-lep, 0:hadron , -1:Non-Tt event
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
      int nMu  = 0;
      int nEle = 0;
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
              if ( abs(dau->pdgId()) == 11 ) nEle++ ;
              if ( abs(dau->pdgId()) == 13 ) nMu++  ;
	      //cout<<" daughter("<< dau->pdgId()<<") status:"<< dau->status()<<endl;
          }
      }
      if (nlep == 1 ) type = 1;
      if (nEle == 1 && nlep ==1 ) type = 4;
      if (nMu  == 1 && nlep ==1 ) type = 3;
      if (nlep == 2 ) type = 2;
      if (nlep == 0 ) type = 0;
   }

   return type;
}

std::vector<const pat::Jet*> TtEvtSelector::WJetSelection( Handle<std::vector<pat::Jet> >  Jets ) {

   std::vector<const pat::Jet* > jCollection;
   jCollection.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = Jets->begin(); j1 != Jets->end(); j1++)
   {
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;
       if (assTk.size() == 0) continue;
       if ((*j1).pt() < 30. ) continue;

       double bDis_TkCount = j1->bDiscriminator("trackCountingHighEffBJetTags") ;
       double jProb        = j1->bDiscriminator("jetProbabilityBJetTags") ;
       if (bDis_TkCount >=  2. && jProb >=0.2  ) continue;

       if ( j1->towersArea() < 0.03 ) continue;
       
       //double EovH1 = EoverH(*j1) ;
       //if ((*j1).nConstituents() < 5 || EovH1 > 20 || EovH1 < 0.01) continue;
       //if ((*j1).towersArea()/(*j1).pt() > 0.005 ) continue;

       jCollection.push_back( &*j1 );
   }
   
   // sort the seeds by # of own segments
   sort(jCollection.begin(), jCollection.end(), PtDecreasing ) ;

   return jCollection;

}


std::vector<const pat::Jet*> TtEvtSelector::bJetSelection( Handle<std::vector<pat::Jet> >  Jets ) {

   std::vector<const pat::Jet* > jCollection;
   jCollection.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = Jets->begin(); j1 != Jets->end(); j1++)
   {
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;
       if (assTk.size() == 0) continue;
       if ((*j1).pt() < 30. ) continue;
       double bDis_TkCount = j1->bDiscriminator("trackCountingHighEffBJetTags") ;
       double jProb        = j1->bDiscriminator("jetProbabilityBJetTags") ;
       if (bDis_TkCount < 2. || jProb < 0.2 ) continue;

       //double EovH1 = EoverH(*j1) ;
       //if ((*j1).nConstituents() < 5 || EovH1 > 20 || EovH1 < 0.01) continue;
       //if ((*j1).towersArea()/(*j1).pt() > 0.005 ) continue;

       jCollection.push_back( &*j1 );
   }

   // sort the seeds by # of own segments
   sort(jCollection.begin(), jCollection.end(), PtDecreasing ) ;

   return jCollection;
}

bool  TtEvtSelector::HLTSemiSelection( Handle <edm::TriggerResults> triggers, int setup ) {

     bool final = false;
     bool hltJ   = false ;
     bool hltMET = false ;

     switch( setup ) {
    
       case 1:
	  bool hltMu  = triggers->accept(76) ; 
          hltJ   = triggers->accept(25) ;
	  hltMET = triggers->accept(29) ;
	  final = (hltJ && hltMET && hltMu) ? true : false ;
	  break;

       case 2:
	  bool hltEle = triggers->accept(48) ;
          hltJ   = triggers->accept(25) ;
	  hltMET = triggers->accept(29) ;
	  final = (hltJ && hltMET && hltEle) ? true : false ;
	  break;

       case 3:
          bool CandMu  = triggers->accept(139);
	  final =  CandMu  ? true : false ;
          break;

       case 4:
	  bool CandEle = triggers->accept(133);
	  final =  CandEle  ? true : false ;
     }

     return final ;
}

void TtEvtSelector::TriggerStudy( Handle <edm::TriggerResults> triggers, int topo, int setup ,HTOP9* histo9 ) {

   // trigger summeray for 3 different topo 
   edm::TriggerNames trigNames( *triggers );
   for (size_t i=0; i< triggers->size(); i++ ) {

       /*
       string triggered = triggers->accept(i) ? "Yes" : "No" ;
       cout<<" path("<<i<<") accepted ? "<< triggered ;
       cout<<" trigName: "<< trigNames.triggerName(i)<<endl;
       if ( triggers->accept(i) && topo == 0 ) histo9->Fill9k0( i );
       if ( triggers->accept(i) && topo == 1 ) histo9->Fill9k1( i );
       if ( triggers->accept(i) && topo == 2 ) histo9->Fill9k2( i );
       */

       if ( triggers->accept(i) ) histo9->Fill9k( i, topo );
       if ( triggers->accept(i) && (topo ==3 || topo ==4 )) histo9->Fill9k( i, 1 );

   }

   // study for the HLT selection
   bool hltMuSel  = HLTSemiSelection(triggers, 3) ;
   bool hltEleSel = HLTSemiSelection(triggers, 4) ;

   bool hltSel = HLTSemiSelection(triggers, setup) ; // temporary useless
   hltSel = ( hltMuSel || hltEleSel ) ? true:false ;

   int hlt = 0 ;
   if ( topo==3 || topo ==4 ) { 
      hlt = hltSel ? 1: -1 ;
      histo9->Fill9l(hlt) ;
   }

   if (topo==2) hlt = hltSel ? 2: -2 ;
   if (topo==3) hlt = hltMuSel ? 3: -3 ;
   if (topo==4) hlt = hltEleSel ? 4: -4 ;
   if (topo==0) hlt = hltSel ? 6: -6 ;
   histo9->Fill9l(hlt) ;

}


//define this as a plug-in
