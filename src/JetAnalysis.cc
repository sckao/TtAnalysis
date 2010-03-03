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
  isData            = iConfig.getUntrackedParameter<bool>   ("isData");
  jetSetup          = iConfig.getParameter<std::vector<double> > ("jetSetup");
  bTagAlgo          = iConfig.getUntrackedParameter<string> ( "bTagAlgo" );
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 

  evtSelected = new TtEvtSelector( iConfig );
  ttMuon      = new TtMuon( iConfig );
  ttJet       = new TtJet( iConfig );
  ttMET       = new TtMET( iConfig );
  ttEle       = new TtElectron( iConfig );

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

  hb_Et30   = new HBJet();

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
   hb_Et30->Write("Jet_Et30", theFile);
   hMET_J30->Write("Jet_Et30", theFile);


   //Release the memory
   delete hJ_Et20;
   delete hJ_Et25;
   delete hJ_Et30;

   delete hb_Et30;

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

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<pat::MET> > met;
   iEvent.getByLabel(metSrc, met);

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

   // Handle<std::vector<pat::TriggerPrimitive> >  triggers;
   // Handle<edm::TriggerResults>  triggers;
   // iEvent.getByLabel(triggerSrc, triggers);

   //cout<<" ***** new Event start ***** "<<endl;
   evtIt++;

   // 0. select the semi-lep events and objects
   int pass = evtSelected->eventSelection( 1, 5, iEvent, "patMET" );

   //1. Jet Et threshold analysis
   std::vector<bool> bTags1;
   std::vector<bool> bTags2;
   std::vector<bool> bTags3;
   std::vector<const reco::Candidate*> isoMu = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> theJets1 = ttJet->JetSelection( jets, isoMu, 20, jetSetup[2], hJ_Et20, &bTags1, bTagAlgo );
   std::vector<const reco::Candidate*> theJets2 = ttJet->JetSelection( jets, isoMu, 25, jetSetup[2], hJ_Et25, &bTags2, bTagAlgo );
   std::vector<const reco::Candidate*> theJets3 = ttJet->JetSelection( jets, isoMu, 30, jetSetup[2], hJ_Et30, &bTags3, bTagAlgo );

   ttJet->bTagAnalysis( iEvent, jets, hb_Et30 );

   if ( isoMu.size() > 0 ) {
      ttJet->MuonAndJet( theJets1, isoMu[0] , hJ_Et20 );
      ttJet->MuonAndJet( theJets2, isoMu[0] , hJ_Et25 );
      ttJet->MuonAndJet( theJets3, isoMu[0] , hJ_Et30 );
   }
   if ( pass > -2) {
      ttJet->JetEtSpectrum( theJets1, hJ_Et20 );
      ttJet->JetdRAnalysis( theJets1, hJ_Et20 );
      if ( theJets1.size() == 4 ) {
         std::vector<const reco::Candidate*> outJets1 = ttJet->SoftJetSelection( jets, isoMu, 20, jetSetup[2], &bTags1, bTagAlgo, hJ_Et20 ) ;
      }
      ttJet->JetEtSpectrum( theJets2, hJ_Et25 );
      ttJet->JetdRAnalysis( theJets2, hJ_Et25 );
      if ( theJets2.size() == 4 ) { 
         std::vector<const reco::Candidate*> outJets2 = ttJet->SoftJetSelection( jets, isoMu, 25, jetSetup[2], &bTags2, bTagAlgo, hJ_Et25 ) ;
      }
      ttJet->JetEtSpectrum( theJets3, hJ_Et30 );
      ttJet->JetdRAnalysis( theJets3, hJ_Et30 );
      if ( theJets3.size() == 4 ) { 
         std::vector<const reco::Candidate*> outJets3 = ttJet->SoftJetSelection( jets, isoMu, 30, jetSetup[2], &bTags3, bTagAlgo, hJ_Et30 ) ;
      }
   }

   std::vector<const reco::Candidate*> isoEle = ttEle->IsoEleSelection( electrons );
   if ( isoEle.size() > 0 ) {
      ttJet->ElectronAndJet( theJets1, isoEle[0] , hJ_Et20 );
      ttJet->ElectronAndJet( theJets2, isoEle[0] , hJ_Et25 );
      ttJet->ElectronAndJet( theJets3, isoEle[0] , hJ_Et30 );
   }
   
   if ( !isData ) {
      Handle<std::vector<reco::GenParticle> > genParticles;
      iEvent.getByLabel(genSrc, genParticles);

      Handle<std::vector<reco::GenJet> > genJets;
      iEvent.getByLabel(genJetSrc, genJets);

      int topo = evtSelected->MCEvtSelection(genParticles);
      if ( topo == 1 ) {
         ttMET->METandNeutrino( isoMu, theJets1, met, genParticles, hMET_J20 );
         ttMET->METandNeutrino( isoMu, theJets2, met, genParticles, hMET_J25 );
         ttMET->METandNeutrino( isoMu, theJets3, met, genParticles, hMET_J30 );
      }
   }

}


//define this as a plug-in
//DEFINE_FWK_MODULE(JetAnalysis);
