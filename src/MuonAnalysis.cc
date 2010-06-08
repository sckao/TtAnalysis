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

  theFile->mkdir("Muon_All");
  theFile->cd();
  theFile->mkdir("Ele_Iso");
  theFile->cd();
  theFile->mkdir("Ele_All");
  theFile->cd();


  hMu_All  = new HOBJ3();
  hEl_All  = new HOBJ4();
  hEl_Iso  = new HOBJ4();

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
   hMu_All->Write("Muon_All", theFile);
   theFile->cd();
   hEl_Iso->Write("Ele_Iso",theFile);
   theFile->cd();
   hEl_All->Write("Ele_All", theFile);
   theFile->cd();

   if (debug) cout << "  Number of Events"<< evtIt << endl;
   //Release the memory

   delete hMu_All;

   delete hEl_All;
   delete hEl_Iso;

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
   ttMuon->PatMuonScope( muons, iEvent, hMu_All);

   //ttMuon->MuonTrigger( muons, triggers );
   //3. Electron Studies
   ttEle->PatEleScope( electrons, hEl_All );

   /*
   if ( !isData ) {
      //std::vector<const reco::Candidate*> mcMuons = MCMatching->matchMuon(genParticles, isoMuReco );
      //ttMuon->matchedMuonAnalysis( mcMuons, hMu_MC );
      std::vector<const reco::Candidate*> genCollects = MCMatching->GenTtCollection( genParticles ) ;

      std::vector<const pat::Muon*> mcMuons = MCMatching->matchMuon( genCollects, isoMuons );
      ttMuon->PatMuonScope( mcMuons, iEvent, hMu_MC );

      std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron( genCollects, isoEle );
      ttEle->matchedElectronAnalysis( mcElectrons, hEl_MC );
   }
   */
}

//define this as a plug-in
//DEFINE_FWK_MODULE(MuonAnalysis);
