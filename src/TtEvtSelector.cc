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
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");
   */
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  tcmetSrc          = iConfig.getParameter<edm::InputTag> ("tcMetSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  jptSrc            = iConfig.getParameter<edm::InputTag> ("jptSource");
  bTagAlgo          = iConfig.getUntrackedParameter<string> ("bTagAlgo");

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

// Event selection only

int TtEvtSelector::eventSelection( Handle<std::vector<pat::Muon> > rMu,  Handle<std::vector<pat::Electron> > rE,
                                   Handle<std::vector<pat::Jet> >  rJet, double jetEtThreshold ) {

 int pass = -1 ; 
 std::vector<const reco::Candidate*> IsoElecs = ttEle->IsoEleSelection( rE ) ;
 int nEle = IsoElecs.size();
 std::vector<const reco::Candidate*> IsoMuons = ttMuon->IsoMuonSelection( rMu ) ;
 int nMu = IsoMuons.size();

 std::vector<const reco::Candidate*> goodJets = ttJet->JetSelection( rJet, IsoMuons, jetEtThreshold );

 std::vector<const reco::Candidate*> jet20 = ttJet->JetSelection( rJet, IsoMuons, 20. );
 bool pure4Jet = (jet20.size() == 4 ) ? true : false ; 
 pure4Jet = true ; 

 //if ( goodJets.size() > 3 && goodJets[0].pt() > 60. && goodJets[3].pt() > 40. ) jetTightSelect = true ; 

 int nJetEvt = static_cast<int>( goodJets.size() );

 if ( nMu == 1 && nEle == 0 && pure4Jet ) pass = nJetEvt ;

 return pass;

}

// Event selection only
int TtEvtSelector::eventSelection( Handle<std::vector<pat::Muon> > rMu,  Handle<std::vector<pat::Electron> > rE,
                                   Handle<std::vector<reco::CaloJet> > rJet, double jetEtThreshold ) {
  
 int pass = -1 ;
 std::vector<const reco::Candidate*> IsoElecs = ttEle->IsoEleSelection( rE ) ;
 int nEle = IsoElecs.size();
 std::vector<const reco::Candidate*> IsoMuons = ttMuon->IsoMuonSelection( rMu ) ;
 int nMu = IsoMuons.size();

 std::vector<const reco::Candidate*> goodJets = ttJet->JetSelection( rJet, IsoMuons, jetEtThreshold );

 std::vector<const reco::Candidate*> jet20 = ttJet->JetSelection( rJet, IsoMuons, 20. );
 bool pure4Jet = (jet20.size() == 4 ) ? true : false ;
 pure4Jet = true ;

 //if ( goodJets.size() > 3 && goodJets[0].pt() > 60. && goodJets[3].pt() > 40. ) jetTightSelect = true ; 

 int nJetEvt = static_cast<int>( goodJets.size() );

 if ( nMu == 1 && nEle == 0 && pure4Jet ) pass = nJetEvt ;

 return pass;

}

// Event selection plus object selection 
// for btagging and no btagging
int TtEvtSelector::eventSelection( int topo, double JetEtCut, std::vector<const reco::Candidate*>& isoLep,  std::vector<const reco::Candidate*>& selectedJets, LorentzVector& metp4, const edm::Event& iEvent, string MetType, string JetType, std::vector<bool>* bTags ){

   // retrieve the reco-objects
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);
   
   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<reco::MET> > tcmet;
   iEvent.getByLabel(tcmetSrc, tcmet);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::CaloJet> > jpts;
   iEvent.getByLabel(jptSrc, jpts);

   std::vector<const reco::Candidate*> isoMu = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> isoEl = ttEle->IsoEleSelection( electrons );
   if ( JetType == "patJet" ) selectedJets = ttJet->JetSelection( jets, isoMu, JetEtCut, NULL ,bTags, bTagAlgo );
   if ( JetType == "jptJet" ) selectedJets = ttJet->JetSelection( jpts, isoMu, JetEtCut );
   
   int pass = -1;

   int nEle = isoEl.size();
   int nMu  = isoMu.size();
   int nJet = selectedJets.size();
   int nLep = nEle + nMu ;
   std::vector<bool> softJbTags;
   std::vector<const reco::Candidate*> additionalJets = ttJet->SoftJetSelection( jets, isoMu, JetEtCut, &softJbTags, bTagAlgo );
   if ( additionalJets.size() > 0  && additionalJets[0]->et() > 20. ) {
      selectedJets.push_back( additionalJets[0] );
      if ( bTags != NULL ) bTags->push_back( softJbTags[0] );
   }


   // hadronic 
   if ( topo == 0 && nMu == 0 && nEle == 0 ) {
      pass = nJet ;
   }

   // muon + jets
   if ( topo == 1 && nMu == 1 && nEle == 0 ) {
      pass = nJet ;
      isoLep = isoMu;
   }
   // dilepton e mu
   if ( topo == 2 && nLep == 2 ) {
      pass = nJet ;
      if (isoMu.size() ==2 ) isoLep = isoMu;
      if (isoEl.size() ==2 ) isoLep = isoEl;
      if (isoMu.size() ==1 ) {
         isoLep = isoMu;
         isoLep.push_back( isoEl[0] );
      }
   }
   // e + jets
   if ( topo == 3 && nMu == 0 && nEle == 1 ) {
      pass = nJet ;
      isoLep = isoEl;
   }

   if ( MetType == "patMet" && mets->size()  > 0 ) metp4 = (*mets)[0].p4() ;
   if ( MetType == "tcMet"  && tcmet->size() > 0 ) metp4 = (*tcmet)[0].p4() ;
   if ( MetType == "evtMet" ) metp4 = ttMET->METfromObjects( isoLep, selectedJets );
   if ( MetType == "neuMet" ) metp4 = ttMET->METfromNeutrino( iEvent );

   return pass;

}

// for B-tagging
int TtEvtSelector::eventSelection( int topo, double JetEtCut, std::vector<const reco::Candidate*>& isoLep,  std::vector<const reco::Candidate*>& selectedWJets, std::vector<const reco::Candidate*>& selectedbJets, LorentzVector& metp4, const edm::Event& iEvent, string MetType ){

   // retrieve the reco-objects
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);
   
   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<reco::MET> > tcmet;
   iEvent.getByLabel(tcmetSrc, tcmet);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::CaloJet> > jpts;
   iEvent.getByLabel(jptSrc, jpts);

   std::vector<const reco::Candidate*> isoMu    = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> isoEl    = ttEle->IsoEleSelection( electrons );
   //std::vector<const reco::Candidate*> goodJets = ttJet->JetSelection( jets, isoMu, JetEtCut );

   selectedWJets = ttJet->WJetSelection( jets, isoLep, JetEtCut );
   selectedbJets = ttJet->bJetSelection( jets, isoLep, JetEtCut, bTagAlgo );

   std::vector<bool> bTags;
   std::vector<const reco::Candidate*> additionalJets = ttJet->SoftJetSelection( jets, isoMu, JetEtCut, &bTags );
   if ( additionalJets.size() > 0  && additionalJets[0]->et() > 20. && !bTags[0] ) selectedWJets.push_back( additionalJets[0] );
   if ( additionalJets.size() > 0  && additionalJets[0]->et() > 20. &&  bTags[0] ) selectedbJets.push_back( additionalJets[0] );

  

   int pass = -1;

   int nEle = isoEl.size();
   int nMu  = isoMu.size();
   int nJet = selectedbJets.size() + selectedWJets.size();
   int nLep = nEle + nMu ;

   if ( topo == 0 && nMu == 0 && nEle == 0 ) {
      pass = nJet ;
   }
   if ( topo == 1 && nMu == 1 && nEle == 0 ) {
      pass = nJet ;
      isoLep = isoMu;
   }
   if ( topo == 2 && nLep == 2 ) {
      pass = nJet ;
      if (isoMu.size() ==2 ) isoLep = isoMu;
      if (isoEl.size() ==2 ) isoLep = isoEl;
      if (isoMu.size() ==1 ) {
         isoLep = isoMu;
         isoLep.push_back( isoEl[0] );
      }
   }
   if ( topo == 3 && nMu == 0 && nEle == 1 ) {
      pass = nJet ;
      isoLep = isoEl;
   }

   if ( MetType == "patMet" && mets->size()  > 0 ) metp4 = (*mets)[0].p4() ;
   if ( MetType ==  "tcMet" && tcmet->size() > 0 ) metp4 = (*tcmet)[0].p4() ;
   if ( MetType == "neuMet" ) metp4 = ttMET->METfromNeutrino( iEvent );

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


