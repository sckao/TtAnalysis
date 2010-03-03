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
  isData            = iConfig.getUntrackedParameter<bool>   ("isData");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  muSetup           = iConfig.getParameter<std::vector<double> >("muSetup");

  //evtSelected = new TtEvtSelector( iConfig );
  ttMuon      = new TtMuon( iConfig );
  MCMatching  = new TtMCMatching();
  ttEle       = new TtElectron( iConfig );

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

   //delete evtSelected;
   delete ttMuon;
   delete ttEle;
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
   iEvent.getByLabel(recoMuon,"",recomuons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   //Handle<std::vector<pat::TriggerPrimitive> >  triggers;

   // Handle<edm::TriggerResults>  triggers;
   // iEvent.getByLabel(triggerSrc, triggers);

   //cout<<" ***** new Event start ***** "<<endl;
   evtIt++;

   // 0. select the semi-lep events and objects
   //int pass = evtSelected->eventSelection( 1, 30., iEvent, "patMET" );
   //int topo = evtSelected->MCEvtSelection(genParticles);

   // 2. Muon Isolation analysis

   double newsets[] = { 5, 2.4, 10 };
   std::vector<double> muSetup0( newsets, newsets+3 ) ;
   std::vector<const pat::Muon*> mostMuons = ttMuon->GeneralMuonSelection( muons, iEvent, muSetup0 );
   ttMuon->PatMuonScope( mostMuons, iEvent, hMu_All);

   std::vector<double> muSetup1=  muSetup ;
   std::vector<const pat::Muon*> isoMuons = ttMuon->GeneralMuonSelection( muons, iEvent, muSetup1 );
   ttMuon->PatMuonScope( isoMuons, iEvent, hMu_Iso);

   //ttMuon->MuonTrigger( muons, triggers );
   //3. Electron Studies
   std::vector<const reco::Candidate*> isoEle   = ttEle->IsoEleSelection(electrons, hEl_All, hEl_Iso );

   if ( !isData ) {
      /*
      std::vector<const reco::Candidate*> isoMuReco ;
      for (size_t i=0; i< isoMuons.size(); i++) {
          isoMuReco.push_back( isoMuons[i] );
      }*/
      //std::vector<const reco::Candidate*> mcMuons = MCMatching->matchMuon(genParticles, isoMuReco );
      //ttMuon->matchedMuonAnalysis( mcMuons, hMu_MC );
      std::vector<const pat::Muon*> mcMuons = MCMatching->matchMuon(genParticles, isoMuons );
      ttMuon->PatMuonScope( mcMuons, iEvent, hMu_MC );

      std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron(genParticles, isoEle );
      ttEle->matchedElectronAnalysis( mcElectrons, hEl_MC );
   }

}

//define this as a plug-in
//DEFINE_FWK_MODULE(MuonAnalysis);
