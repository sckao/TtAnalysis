// -*- C++ -*-
//
// Package:    TtAnalysis
// Class:      TtAnalysis
// 
/**\class TtAnalysis TtAnalysis.cc PhysicsTools/TtAnalysis/src/TtAnalysis.cc

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
#include "TtAnalysis.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// constants, enums and typedefs
//
// static data member definitions
//
// constructors and destructor
using namespace edm;
using namespace std;
TtAnalysis::TtAnalysis(const edm::ParameterSet& iConfig) 
: evtIt(0)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  photonSrc         = iConfig.getParameter<edm::InputTag> ("photonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  //caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 
  triggerSrc        = iConfig.getParameter<edm::InputTag> ("trigSource");
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");

  MCMatching  = new TtMCMatching();
  evtSelected = new TtEvtSelector( iConfig );
  ttMuon      = new TtMuon( iConfig );
  ttEle       = new TtElectron( iConfig );
  ttGam       = new TtPhoton();
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );
  semiSol     = new TtSemiEventSolution( iConfig );
  ttEff       = new TtEfficiency();

  theFile->cd();
  theFile->mkdir("Jets");
  theFile->cd();
  theFile->mkdir("METs");
  theFile->cd();
  theFile->mkdir("Muons");
  theFile->cd();
  theFile->mkdir("Ele");
  theFile->cd();
  theFile->mkdir("Gam");
  theFile->cd();
  theFile->mkdir("MObjs");
  theFile->cd();
  theFile->mkdir("BJets");
  theFile->cd();
  theFile->mkdir("WJets");
  theFile->cd();
  theFile->mkdir("Tops");
  theFile->cd();
  theFile->mkdir("Ws");
  theFile->cd();
  
  histos.hJet   = new HTOP1();
  histos.hMET   = new HTOP2();
  histos.hMuon  = new HTOP3();
  histos.hEle   = new HTOP4();
  histos.hGam   = new HTOP5();
  histos.hMObj  = new HTOP6();
  histos.hBJet  = new HTOP7();
  histos.hWJet  = new HTOP8();
  histos.hTop   = new HTOP9();
  for (int i=0; i< 4; i++) {
      histos.hWs[i]   = new HTOP10(i) ;
      histos.hTops[i] = new HTOP11(i) ;
  }

}


TtAnalysis::~TtAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (debug) cout << "[TtAnalysis Analysis] Destructor called" << endl;

   delete MCMatching;
   delete evtSelected;
   delete ttMuon;
   delete ttEle;
   delete ttGam;
   delete ttMET;
   delete ttJet;
   delete ttEff; 
   delete semiSol;
   
   theFile->cd();
   histos.hJet->Write("Jets",theFile);
   theFile->cd();
   histos.hMET->Write("METs",theFile);
   theFile->cd();
   histos.hMuon->Write("Muons",theFile);
   theFile->cd();
   histos.hEle->Write("Ele",theFile);
   theFile->cd();
   histos.hGam->Write("Gam",theFile);
   theFile->cd();
   histos.hMObj->Write("MObjs",theFile);
   theFile->cd();
   histos.hBJet->Write("BJets",theFile);
   theFile->cd();
   histos.hWJet->Write("WJets",theFile);
   theFile->cd();
   histos.hTop->Write("Tops",theFile);
   theFile->cd();
   for (int i=0; i<4; i++) {
     histos.hWs[i]->Write("Ws",theFile);
     theFile->cd();
     histos.hTops[i]->Write("Tops",theFile, i );
     theFile->cd();
   } 
   if (debug) cout << "[Wrote the histograms]" << endl;

   //Release the memory
   //delete &histos;

   if (debug) cout << "[Release the memory of historgrams]" << endl;

   //Close the Root file
   theFile->Close();
   if (debug) cout << "************* Finished writing histograms to file" << endl;

}

//
// member functions
//

// ------------ method called to for each event  ------------
void TtAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // retrieve the reco-objects
   cout<<"*** New Event Start **** "<<endl;
   
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<reco::Muon> > recomuons;
   iEvent.getByLabel(recoMuon, recomuons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<pat::Photon> > photons;
   iEvent.getByLabel(photonSrc, photons);

   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   Handle<std::vector<reco::GenJet> > genJets;
   iEvent.getByLabel(genJetSrc, genJets);

   //Handle<CaloTowerCollection>  caloTowers;
   //iEvent.getByLabel(caloSrc, caloTowers);

   //cout<<" ***** new Event start ***** "<<endl;
   evtIt++;
   //int eventId = evtIt + (iEvent.id().run()*100000) ;

   // The Feed the gen Tree 

   // 0. select the semi-lep events and objects
   int jpass = evtSelected->eventSelection( 1, 30., iEvent, "patMET");
   int topo  = evtSelected->MCEvtSelection(genParticles);

   cout<<" jpass : "<<jpass<<" topo :"<<topo <<endl;
   //  And calculate the selection efficiency
   bool passSelect = ( jpass > 3 ) ? true : false ;
   ttEff->EventEfficiency( topo, passSelect, histos.hTop );

   // check generator process 
   //MCMatching->CheckGenParticle( genParticles );

   // 1. Build semi-mu tt events and compare the result
 
   //3. Muon and Leptonic W analysis
   std::vector<const reco::Candidate*> isoMu   = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> theJets = ttJet->JetSelection( jets, isoMu, 30 ) ;

   //4. Jet Studies
   ttJet->jetAnalysis(jets, histos.hJet);
   //ttJet->JetTrigger( jets, triggers );

   std::vector<const reco::Candidate*> genCollects = MCMatching->GenTtCollection( genParticles ) ;

   ttJet->JetMatchedMuon( jets, muons, iEvent, iSetup, histos.hMuon, true );

   //5. MET from PAT
   ttMET->metAnalysis( iEvent, histos.hMET);
   if ( jpass > 3 ) {
      ttMET->MetAndMuon(mets, isoMu, histos.hMET, jpass );
      ttMET->MetAndJets(mets, theJets, histos.hMET );
   }
   
   //6. Electron Studies
   ttEle->ElectronAnalysis(electrons, histos.hEle);

   //7. photon studies
   ttGam->PhotonAnalysis(photons, histos.hGam);
   cout<<" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "<<endl;


}


//define this as a plug-in => In SealModule.C
//DEFINE_FWK_MODULE(TtAnalysis);
