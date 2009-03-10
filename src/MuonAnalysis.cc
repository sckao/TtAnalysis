// -*- C++ -*-
//
// Package:    MuonAnalysis
// Class:      MuonAnalysis
// 
/**\class MuonAnalysis MuonAnalysis.cc PhysicsTools/TtAnalysis/src/MuonAnalysis.cc

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
#include "MuonAnalysis.h"
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
MuonAnalysis::MuonAnalysis(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

  evtSelected = new TtEvtSelector( iConfig );
  ttMuon      = new TtMuon();
  ttEle       = new TtElectron();
  ttJet       = new TtJet( iConfig );
  MCMatching  = new TtMCMatching();

  evtIt = 0;
  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");

  theFile->mkdir("Muon_Iso");
  theFile->cd();
  theFile->mkdir("Muon_All");
  theFile->cd();
  theFile->mkdir("Muon_MC");
  theFile->cd();
  theFile->mkdir("Ele_Iso");
  theFile->cd();
  theFile->mkdir("Ele_All");
  theFile->cd();
  theFile->mkdir("Ele_MC");
  theFile->cd();


  hMu_All  = new HOBJ3();
  hMu_Iso  = new HOBJ3();
  hMu_MC   = new HOBJ3();

  hEl_All  = new HOBJ4();
  hEl_Iso  = new HOBJ4();
  hEl_MC   = new HOBJ4();

}


MuonAnalysis::~MuonAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (debug) cout << "[MuonAnalysis Analysis] Destructor called" << endl;

   delete evtSelected;
   delete ttMuon;
   delete ttEle;
   delete ttJet;
   delete MCMatching;

   theFile->cd();
   hMu_Iso->Write("Muon_Iso",theFile);
   theFile->cd();
   hMu_All->Write("Muon_All", theFile);
   theFile->cd();
   hMu_MC->Write("Muon_MC", theFile);
   theFile->cd();
   hEl_Iso->Write("Ele_Iso",theFile);
   theFile->cd();
   hEl_All->Write("Ele_All", theFile);
   theFile->cd();
   hEl_MC->Write("Ele_MC", theFile);
   theFile->cd();

   if (debug) cout << "  Number of Events"<< evtIt << endl;
   //Release the memory

   delete hMu_All;
   delete hMu_Iso;
   delete hMu_MC;

   delete hEl_All;
   delete hEl_Iso;
   delete hEl_MC;

   //Close the Root file
   theFile->Close();
   if (debug) cout << "************* Finished writing histograms to file" << endl;
}

//
// member functions
//

// ------------ method called to for each event  ------------
void MuonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   Handle<std::vector<reco::GenJet> > genJets;
   iEvent.getByLabel(genJetSrc, genJets);

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   //Handle<std::vector<pat::TriggerPrimitive> >  triggers;

   // Handle<edm::TriggerResults>  triggers;
   // iEvent.getByLabel(triggerSrc, triggers);

   // Initial the histograms
   HTOP3  *histo3 = 0;
   HTOP4  *histo4 = 0;

   //cout<<" ***** new Event start ***** "<<endl;
   evtIt++;

   // 0. select the semi-lep events and objects
   //int pass = evtSelected->eventSelection(muons, electrons, jets, 20.);
   //int topo = evtSelected->MCEvtSelection(genParticles);


   // 2. Muon Isolation analysis
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection( muons, hMu_All, hMu_Iso );
   std::vector<const reco::Candidate*> mcMuons = MCMatching->matchMuon(genParticles, isoMuons, histo3, false);
   ttMuon->matchedMuonAnalysis( mcMuons, hMu_MC );
     
   //ttMuon->MuonTrigger( muons, triggers );

   //3. Electron Studies
   std::vector<const reco::Candidate*> isoEle      = ttEle->IsoEleSelection(electrons, hEl_All, hEl_Iso );
   std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron(genParticles, isoEle, histo4, false);
   ttEle->matchedElectronAnalysis( mcElectrons, hEl_MC );

}


//define this as a plug-in
//DEFINE_FWK_MODULE(MuonAnalysis);
