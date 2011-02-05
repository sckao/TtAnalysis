// -*- C++ -*-
//
// Package:    TtNtupleProd
// Class:      TtNtupleProd
// 
/**\class TtNtupleProd TtNtupleProd.cc PhysicsTools/TtNtupleProd/src/TtNtupleProd.cc

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
#include "TtNtupleProd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// constants, enums and typedefs
//
// static data member definitions
//
// constructors and destructor
using namespace edm;
using namespace std;
TtNtupleProd::TtNtupleProd(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  //isData            = iConfig.getUntrackedParameter<bool>   ("isData");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  genSrc            = iConfig.getParameter<edm::InputTag>   ("genParticles"); 
  evtIt             = iConfig.getUntrackedParameter<int> ("eventId");
  nJets             = iConfig.getUntrackedParameter<int> ("numberOfJets");

  // Create the root file
  theFile     = new TFile(rootFileName.c_str(), "RECREATE");
  //MCMatching  = new TtMCMatching();
  semiSol     = new TtSemiEventSolution( iConfig );

  theFile->cd();

  ntuples.muJets = new SolNtp2("muJets" );
  //if ( !isData ) ntuples.mcmTree = new mcNtp("mcmTt");
  //if ( !isData ) ntuples.genTree = new mcNtp("genTt");

}

TtNtupleProd::~TtNtupleProd()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (debug) cout << "[TtNtupleProd Analysis] Destructor called" << endl;

   theFile->cd();
   if (debug) cout << "[Wrote the ntuples]" << endl;
   
   ntuples.muJets->Write();
   //if ( !isData ) ntuples.mcmTree->Write();
   //if ( !isData ) ntuples.genTree->Write();
   

   //Close the Root file
   //theFile->Write();
   theFile->Close();
   if (debug) cout << "************* Finished writing ntuples to file" << endl;

   //Release the memory
   //delete MCMatching;
   delete semiSol;
   
   if (debug) cout << "[Release the memory of ntuples]" << endl;

}

//
// member functions
//

// ------------ method called to for each event  ------------
void TtNtupleProd::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // retrieve the reco-objects
   //cout<<"*** New Event Start **** "<<endl;
   
   // Handle<std::vector<pat::TriggerPrimitive> >  triggers;
   // Handle<edm::TriggerResults>  triggers;
   // iEvent.getByLabel(triggerSrc, triggers);

   //cout<<" ***** new Event start ***** "<<endl;
   evtIt++;

   //cout<<" 1. building "<<endl;
   //semiSol->BuildSemiTt(iEvent, 1, evtIt, &ntuples );
   //cout<<" find N jets event "<<endl;
   semiSol->RecordSolutions(iEvent, 1, evtIt, nJets, ntuples.muJets );

   // The Feed the gen Tree 
   /*
   if ( !isData ) {
      //cout<<" the mc matched event "<<endl;
      Handle<std::vector<reco::GenParticle> > genParticles;
      iEvent.getByLabel(genSrc, genParticles);
      //cout<<" 2. GEN "<<endl;
      //cout<<" 3. MC Matching "<<endl;
      semiSol->MCBuildSemiTt(iEvent, 1, evtIt, &ntuples );
   }
   */
   //cout<<" Done !!! "<<endl;

}

//define this as a plug-in => In SealModule.C
//DEFINE_FWK_MODULE(TtNtupleProd);
