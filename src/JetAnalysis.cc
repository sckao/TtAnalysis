// -*- C++ -*-
//
// Package:    JetAnalysis
// Class:      JetAnalysis
// 
/**\class JetAnalysis JetAnalysis.cc PhysicsTools/TtAnalysis/src/JetAnalysis.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Mon Feb 9 2009
//
//


// system include files
#include <memory>

// user include files
#include "JetAnalysis.h"
#include "FWCore/Framework/interface/MakerMacros.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//
// constructors and destructor
using namespace edm;
using namespace std;
JetAnalysis::JetAnalysis(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

  evtSelected = new TtEvtSelector( iConfig );
  ttMuon      = new TtMuon();
  ttJet       = new TtJet( iConfig );
  ttMET       = new TtMET( iConfig );
  ttEle       = new TtElectron();

  evtIt = 0;
  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");

  theFile->mkdir("Jet_Et20");
  theFile->cd();
  theFile->mkdir("Jet_Et25");
  theFile->cd();
  theFile->mkdir("Jet_Et30");
  theFile->cd();

  hJ_Et20   = new HOBJ1();
  hJ_Et25   = new HOBJ1();
  hJ_Et30   = new HOBJ1();
  hMET_J20  = new HOBJ2();
  hMET_J25  = new HOBJ2();
  hMET_J30  = new HOBJ2();

}


JetAnalysis::~JetAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (debug) cout << "[JetAnalysis Analysis] Destructor called" << endl;
   if (debug) cout << "  Number of Events"<< evtIt << endl;

   delete evtSelected;
   delete ttMuon;
   delete ttEle;
   delete ttJet;
   delete ttMET;

   theFile->cd();
   hJ_Et20->Write("Jet_Et20", theFile);
   hMET_J20->Write("Jet_Et20", theFile);

   //theFile->cd();
   hJ_Et25->Write("Jet_Et25", theFile);
   hMET_J25->Write("Jet_Et25", theFile);

   //theFile->cd();
   hJ_Et30->Write("Jet_Et30", theFile);
   hMET_J30->Write("Jet_Et30", theFile);


   //Release the memory
   delete hJ_Et20;
   delete hJ_Et25;
   delete hJ_Et30;
   delete hMET_J20;
   delete hMET_J25;
   delete hMET_J30;
  

   //Close the Root file
   theFile->Close();
   if (debug) cout << "************* Finished writing histograms to file" << endl;

}

//
// member functions
//

// ------------ method called to for each event  ------------
void JetAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // retrieve the reco-objects
   
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<reco::Muon> > recomuons;
   if ( recoMuon == "paramMuons" ) {
      iEvent.getByLabel(recoMuon,"ParamGlobalMuons",recomuons);
   }
   if ( recoMuon == "muons" ) {
      iEvent.getByLabel(recoMuon,"",recomuons);
   }
   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<pat::MET> > met;
   iEvent.getByLabel(metSrc, met);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   Handle<std::vector<reco::GenJet> > genJets;
   iEvent.getByLabel(genJetSrc, genJets);

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

   // Handle<std::vector<pat::TriggerPrimitive> >  triggers;
   // Handle<edm::TriggerResults>  triggers;
   // iEvent.getByLabel(triggerSrc, triggers);

   //cout<<" ***** new Event start ***** "<<endl;
   evtIt++;

   // 0. select the semi-lep events and objects
   int pass = evtSelected->eventSelection(muons, electrons, jets, 20.);
   //int topo = evtSelected->MCEvtSelection(genParticles);

   //1. Jet Et threshold analysis
   std::vector<const reco::Candidate*> isoMu = ttMuon->IsoMuonSelection( muons );
   std::vector<const pat::Jet*> theJets1 = ttJet->JetSelection( jets, isoMu, 20, hJ_Et20 ) ;
   std::vector<const pat::Jet*> theJets2 = ttJet->JetSelection( jets, isoMu, 25, hJ_Et25 ) ;
   std::vector<const pat::Jet*> theJets3 = ttJet->JetSelection( jets, isoMu, 30, hJ_Et30 ) ;

   if ( pass > -1 ) {
      ttJet->MuonAndJet( theJets1, isoMu[0] , hJ_Et20 );
      ttJet->JetEtSpectrum( theJets1, hJ_Et20 );
      ttJet->JetdRAnalysis( theJets1, hJ_Et20 );
      if ( theJets1.size() == 4 ) ttMET->METandNeutrino( isoMu, theJets1, met, genParticles, hMET_J20 );

      ttJet->MuonAndJet( theJets2, isoMu[0] , hJ_Et25 );
      ttJet->JetEtSpectrum( theJets2, hJ_Et25 );
      ttJet->JetdRAnalysis( theJets2, hJ_Et25 );
      if ( theJets2.size() == 4 ) ttMET->METandNeutrino( isoMu, theJets2, met, genParticles, hMET_J25 );

      ttJet->MuonAndJet( theJets3, isoMu[0] , hJ_Et30 );
      ttJet->JetEtSpectrum( theJets3, hJ_Et30 );
      ttJet->JetdRAnalysis( theJets3, hJ_Et30 );
      if ( theJets3.size() == 4 ) ttMET->METandNeutrino( isoMu, theJets3, met, genParticles, hMET_J30 );
   }

   std::vector<const reco::Candidate*> isoEle = ttEle->IsoEleSelection( electrons );
   if ( isoMu.size() == 1 ) {
      hJ_Et20->Fill_1h( isoEle.size(), theJets1.size() );
      hJ_Et25->Fill_1h( isoEle.size(), theJets2.size() );
      hJ_Et30->Fill_1h( isoEle.size(), theJets3.size() );
   }
  
   
}


//define this as a plug-in
//DEFINE_FWK_MODULE(JetAnalysis);
