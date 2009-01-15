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
TtEvtSelector::TtEvtSelector(const edm::ParameterSet& iConfig)
{
//TtEvtSelector::TtEvtSelector()
   //now do what ever initialization is needed
  /*
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
  ttMuon    = new TtMuon();
  ttEle     = new TtElectron();
  ttJet     = new TtJet( iConfig );
  ttMET     = new TtMET( iConfig );

}


TtEvtSelector::~TtEvtSelector()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete ttMuon;
   delete ttEle;
   delete ttJet;
   delete ttMET;
   
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------
int TtEvtSelector::eventSelection(Handle<std::vector<pat::Muon> > rMu, Handle<std::vector<pat::Electron> > rE,
                                Handle<std::vector<pat::Jet> >   rJet, double jetEtThreshold ) {

 int pass = -1 ; 
 std::vector<const reco::Candidate*> IsoElecs = ttEle->IsoEleSelection( rE ) ;
 int nEle = IsoElecs.size();
 std::vector<const reco::Candidate*> IsoMuons = ttMuon->IsoMuonSelection( rMu ) ;
 int nMu = IsoMuons.size();

 std::vector<pat::Jet> goodJets = ttJet->JetSelection( rJet, IsoMuons );
 
 //if ( goodJets.size() > 3 && goodJets[0].pt() > 60. && goodJets[3].pt() > 40. ) jetTightSelect = true ; 
 /*
 int njets = static_cast<int>(goodJets.size());
 int jetSelect = -1 ; 
 if ( njets > -1 ) jetSelect = 0 ; 
 if ( njets > 0 && goodJets[0].pt() > jetEtThreshold ) jetSelect = 1 ; 
 if ( njets > 1 && goodJets[1].pt() > jetEtThreshold ) jetSelect = 2 ; 
 if ( njets > 2 && goodJets[2].pt() > jetEtThreshold ) jetSelect = 3 ; 
 if ( njets > 3 && goodJets[3].pt() > jetEtThreshold ) jetSelect = 4 ; 
 if ( njets > 4 && goodJets[4].pt() > jetEtThreshold ) jetSelect = 5 ; 
 if ( njets > 5 && goodJets[5].pt() > jetEtThreshold ) jetSelect = 6 ; 
 if ( nMu == 1 && nEle == 0 && jetSelect > -1 ) pass = jetSelect ;
 */

 int nJetEvt = 0;
 for (size_t i=0; i< goodJets.size(); i++) {
     if ( goodJets[i].pt() > jetEtThreshold ) nJetEvt++ ;
 }

 if ( nMu == 1 && nEle == 0 ) pass = nJetEvt ;

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


