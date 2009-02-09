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
  mcMatch   = new TtMCMatching();

}


TtEvtSelector::~TtEvtSelector()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete ttMuon;
   delete ttEle;
   delete ttJet;
   delete ttMET;
   delete mcMatch;
   
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

 std::vector<pat::Jet> goodJets = ttJet->JetSelection( rJet, IsoMuons, jetEtThreshold );
 
 //if ( goodJets.size() > 3 && goodJets[0].pt() > 60. && goodJets[3].pt() > 40. ) jetTightSelect = true ; 

 int nJetEvt = static_cast<int>( goodJets.size() );

 if ( nMu == 1 && nEle == 0 ) pass = nJetEvt ;

 return pass;

}

int TtEvtSelector::MCEvtSelection( Handle<std::vector<reco::GenParticle> > genParticles ) {

   // type 3:semi-electron 2: di-lep, 1:semi-muon, 0:hadron , -1:Non-Tt event
   int type = -1;
   std::vector<reco::Particle> tauColl = mcMatch->ttDecay(genParticles, 15);
   std::vector<reco::Particle> muColl  = mcMatch->ttDecay(genParticles, 13);
   std::vector<reco::Particle> eColl   = mcMatch->ttDecay(genParticles, 11);

   
   int nTau = tauColl.size();
   int nMu  = muColl.size();
   int nEle = eColl.size();
   int nlep = nMu + nEle + nTau ;

   if (nTau == 1 && nlep ==1 ) type = 4;
   if (nEle == 1 && nlep ==1 ) type = 3;
   if (nMu  == 1 && nlep ==1 ) type = 1;
   if (nlep == 2 ) type = 2;
   if (nlep == 0 ) type = 0;
   if (nlep > 2 ) type = -1;

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


