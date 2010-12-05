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
  pvSrc             = iConfig.getParameter<edm::InputTag> ("pvSource");
  pvNDF             = iConfig.getUntrackedParameter<double> ("pvNDOF");
  pvZ               = iConfig.getUntrackedParameter<double> ("pvMaxZ");
  pvRho             = iConfig.getUntrackedParameter<double> ("pvMaxRho");
  beamSpotSrc       = iConfig.getParameter<edm::InputTag> ("beamSpotSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  trigSrc           = iConfig.getParameter<edm::InputTag> ("trigSource");
  recoMetSrc        = iConfig.getParameter<edm::InputTag> ("recoMetSource");
  bTagAlgo          = iConfig.getUntrackedParameter<string> ("bTagAlgo");
  jetSetup          = iConfig.getParameter<std::vector<double> >("jetSetup");
  //ecalRecHitSrc     = iConfig.getParameter<edm::InputTag>("ecalRecHitSource");
  //JEScale           = iConfig.getUntrackedParameter<double> ("JEScale");

  ttEle     = new TtElectron( iConfig );
  ttMuon    = new TtMuon( iConfig );
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

// Event selection plus object selection , new general method using ttCandidate
int TtEvtSelector::eventSelection( int topo, std::vector<ttCandidate>& isoLep,  std::vector<ttCandidate>& selectedJets, std::vector<LorentzVector>& metp4, std::vector<ttCandidate>& vetoInfo, const edm::Event& iEvent, string MetType ){

   // retrieve the reco-objects

   Handle<reco::BeamSpot> beamSpot;
   iEvent.getByLabel(beamSpotSrc, beamSpot );

   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);
   
   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::MET> > recomet;
   iEvent.getByLabel(recoMetSrc, recomet);

   Handle<std::vector<reco::Vertex> > primVtx;
   iEvent.getByLabel(pvSrc, primVtx);

   double PVz = 0 ;
   if (  primVtx->size()  >= 1 ) { 
      const reco::Vertex pv = primVtx->at(0);
      if (  !pv.isFake() )  PVz = pv.z() ;
   }
   //Handle<EcalRecHitCollection> recHits;
   //iEvent.getByLabel( ecalRecHitSrc, recHits);

   int pass = -1;

   std::vector<ttCandidate> isoMu = ttMuon->IsoMuonSelection1( muons, jets, beamSpot, PVz, vetoInfo );
   std::vector<ttCandidate> isoEl = ttEle->IsoEleSelection1( electrons, beamSpot, vetoInfo );
   selectedJets = ttJet->JetSelection1( jets, isoMu, jetSetup[0], jetSetup[2], NULL,  bTagAlgo );
  
   int nEle = isoEl.size();
   int nMu  = isoMu.size();
   int nJet = selectedJets.size();
   int nLep = nEle + nMu ;

   // this probably will exhaust memory fast
   //std::vector<ttCandidate> additionalJets = ttJet->SoftJetSelection1( jets, isoMu, jetSetup[0], jetSetup[2], NULL, bTagAlgo );
   //if ( additionalJets.size() > 0  && additionalJets[0].p4.Et() > 20. ) {
   //   selectedJets.push_back( additionalJets[0] );
   //}

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

   double metCut = 0. ;
   LorentzVector theMetP4 ;
   if ( MetType == "patMet" && mets->size()  > 0 ) theMetP4 = jetSetup[2] * (*mets)[0].p4() ;
   if ( MetType == "recoMet" && recomet->size() > 0 ) theMetP4 = jetSetup[2] * (*recomet)[0].p4() ;
   if ( MetType == "evtMet" ) theMetP4 = ttMET->METfromObjects( isoLep, selectedJets );
   if ( MetType == "neuMet" ) theMetP4 = ttMET->METfromNeutrino( iEvent );
   if ( theMetP4.Et() < metCut ) pass = -1 ;

   metp4.push_back( theMetP4 );

   return pass;

}


// Event selection only , not keeping objects
int TtEvtSelector::eventSelection( int topo, double JetEtCut, const edm::Event& iEvent, string MetType ){

   // retrieve the reco-objects
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);
   
   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<reco::MET> > recomet;
   iEvent.getByLabel(recoMetSrc, recomet);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);


   std::vector<const reco::Candidate*> isoMu = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> isoEl = ttEle->IsoEleSelection( electrons );
   std::vector<const reco::Candidate*> isoLep = isoMu ;
   for (size_t i=0; i< isoEl.size(); i++ ) isoLep.push_back( isoEl[i] ) ;

   std::vector<const reco::Candidate*> theJets;
   theJets = ttJet->JetSelection( jets, isoMu, jetSetup[0], jetSetup[2], NULL , bTagAlgo );
   
   int pass = -1;
   int nEle = isoEl.size();
   int nMu  = isoMu.size();
   int nJet = theJets.size();
   int nLep = nEle + nMu ;

   // hadronic 
   if ( topo == 0 && nMu == 0 && nEle == 0 )  pass = nJet ;
   
   // muon + jets
   if ( topo == 1 && nMu == 1 && nEle == 0 )  pass = nJet ;

   // dilepton e mu
   if ( topo == 2 && nLep == 2 )              pass = nJet ;
   // e + jets
   if ( topo == 3 && nMu == 0 && nEle == 1 )  pass = nJet ;

   //double metCut = 0. ;
   LorentzVector metp4;
   if ( MetType == "patMet" && mets->size()  > 0 ) metp4 = jetSetup[2] * (*mets)[0].p4() ;
   if ( MetType == "recoMet"  && recomet->size() > 0 ) metp4 = jetSetup[2] * (*recomet)[0].p4() ;
   if ( MetType == "evtMet" ) metp4 = ttMET->METfromObjects( isoLep, theJets );
   if ( MetType == "neuMet" ) metp4 = ttMET->METfromNeutrino( iEvent );
   //if ( metp4.Et() < metCut ) pass = -1 ;

   return pass;

}


int TtEvtSelector::MCEvtSelection( Handle<std::vector<reco::GenParticle> > genParticles ) {

   // type 4: tau+jets type 3:electron+jest 2: di-lep, 1:muon+jets, 0:hadron , -1:Non-Tt event
   int type = -1;
   std::vector<const reco::Candidate*> tauColl = mcMatch->ttDecay(genParticles, 15);
   std::vector<const reco::Candidate*> muColl  = mcMatch->ttDecay(genParticles, 13);
   std::vector<const reco::Candidate*> eColl   = mcMatch->ttDecay(genParticles, 11);
   
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

bool TtEvtSelector::VertexSelection( const edm::Event& iEvent, std::vector<double>& pvInfo ){

   // retrieve the reco-objects
   Handle<std::vector<reco::Vertex> > primVtx;
   iEvent.getByLabel(pvSrc, primVtx);

   // primVertex cut
   bool goodVtx = true ;
   if (  primVtx->size()  < 1 ) { 
      goodVtx = false ;  
   }  else {
      const reco::Vertex pv = primVtx->at(0);
      if (  pv.isFake()  )                  goodVtx = false ;
      if (  fabs( pv.z()) >= pvZ  )         goodVtx = false ;
      if (  pv.position().Rho() >= pvRho  ) goodVtx = false ;
      if (  pv.ndof() <= pvNDF  )           goodVtx = false ;
      pvInfo.push_back( pv.z() );
      pvInfo.push_back( pv.position().Rho() );
      pvInfo.push_back( pv.ndof() );
   }

   return goodVtx ;

}


bool TtEvtSelector::TriggerSelection( const edm::Event& iEvent, string trigPathName ) {

   Handle<pat::TriggerEvent> triggerEvent;
   iEvent.getByLabel(trigSrc, triggerEvent);

   /*
   Handle<std::vector<pat::TriggerPath> > trigPaths;
   iEvent.getByLabel("patTrigger", trigPaths);
   for(std::vector<pat::TriggerPath>::const_iterator t1 = trigPaths->begin(); t1!= trigPaths->end(); t1++){
      cout<<" trig path name : "<< t1->name() <<endl;
   } 
   */

   bool passTrig = false ;
   pat::TriggerEvent const * trig = &*triggerEvent;
   if ( trig->wasRun() && trig->wasAccept() ) {

       pat::TriggerPath const * thePath = trig->path( trigPathName );
       if ( thePath != 0 && thePath->wasAccept() ) passTrig = true;    
   }
  
   return passTrig ;
}


