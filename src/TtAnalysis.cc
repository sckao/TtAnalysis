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

//
// static data member definitions
//
// constructors and destructor
using namespace edm;
using namespace std;
TtAnalysis::TtAnalysis(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  needTree          = iConfig.getUntrackedParameter<bool>   ("needTree");
  trigOn            = iConfig.getUntrackedParameter<bool>   ("trigOn");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  photonSrc         = iConfig.getParameter<edm::InputTag> ("photonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 
  triggerSrc        = iConfig.getParameter<edm::InputTag> ("triggerSource");

  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

  evtSelected = new TtEvtSelector();
  MCMatching  = new TtMCMatching();
  ttMuon      = new TtMuon();
  ttEle       = new TtElectron();
  ttGam       = new TtPhoton();
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );
  ttEff       = new TtEfficiency();
  semiSol     = new TtSemiEventSolution( iConfig );

  evtIt = 0;
  // Create the root file
  theFile = new TFile(rootFileName.c_str(), "RECREATE");

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
  theFile->mkdir("MJets");
  theFile->cd();
  theFile->mkdir("BJets");
  theFile->cd();
  theFile->mkdir("WJets");
  theFile->cd();
  theFile->mkdir("Tops");
  theFile->cd();

  h_Jet   = new HTOP1("Jets_");
  h_MET   = new HTOP2("MET_");
  h_Muon  = new HTOP3("Muons_");
  h_Ele   = new HTOP4("Ele_");
  h_Gam   = new HTOP5("Gam_");
  h_MJet  = new HTOP6("MJets_");
  h_BJet  = new HTOP7("BJets_");
  h_WJet  = new HTOP8("WJets_");
  h_Top   = new HTOP9("Tops_");

  if (needTree) t_Jet  = new NJet();
}


TtAnalysis::~TtAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (debug) cout << "[TtAnalysis Analysis] Destructor called" << endl;

   delete evtSelected;
   delete MCMatching;
   delete ttMuon;
   delete ttEle;
   delete ttGam;
   delete ttMET;
   delete ttJet;
   delete ttEff; 
   delete semiSol;

   theFile->cd();
   theFile->cd("Jets");
   h_Jet->Write();

   theFile->cd();
   theFile->cd("METs");
   h_MET->Write();

   theFile->cd();
   theFile->cd("Muons");
   h_Muon->Write();

   theFile->cd();
   theFile->cd("Ele");
   h_Ele->Write();

   theFile->cd();
   theFile->cd("Gam");
   h_Gam->Write();

   theFile->cd();
   theFile->cd("MJets");
   h_MJet->Write();

   theFile->cd();
   theFile->cd("BJets");
   h_BJet->Write();

   theFile->cd();
   theFile->cd("WJets");
   h_WJet->Write();

   theFile->cd();
   theFile->cd("Tops");
   h_Top->Write();

   //Release the memory
   delete h_Jet;
   delete h_MET;
   delete h_Muon;
   delete h_Ele;
   delete h_Gam;
   delete h_MJet;
   delete h_BJet;
   delete h_WJet;
   delete h_Top;
  
   if (needTree) {
      theFile->cd();
      t_Jet->Write();
      delete t_Jet;
   }

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

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

   //Handle<std::vector<pat::TriggerPrimitive> >  triggers;

   // Handle<edm::TriggerResults>  triggers;
   // iEvent.getByLabel(triggerSrc, triggers);

   // Initial the histograms
   HTOP1 *histo1 = 0;
   HTOP2 *histo2 = 0;
   HTOP3 *histo3 = 0;
   HTOP4 *histo4 = 0;
   HTOP5 *histo5 = 0;
   HTOP6 *histo6 = 0;
   HTOP7 *histo7 = 0;
   HTOP8 *histo8 = 0;
   HTOP9 *histo9 = 0;

   cout<<" ***** new Event start ***** "<<endl;

   evtIt++;
   int eventId = evtIt + (iEvent.id().run()*100000) ;

   // The Tree Feeder 
   if ( needTree ) {
      NJet *jtree = 0;
      jtree = t_Jet;
      MCMatching->MCTreeFeeder( genParticles, jtree, eventId);
      ttMuon->MuonTreeFeeder(muons, jtree, eventId);
      ttEle->ElectronTreeFeeder(electrons, jtree, eventId);
      ttGam->PhotonTreeFeeder( photons, jtree, eventId);
      ttJet->JetTreeFeeder(jets,jtree,eventId);
      ttMET->MetTreeFeeder(mets,jtree,eventId);
   }

   // 0. select the semi-lep events and objects
   bool pass = evtSelected->eventSelection(muons, electrons, jets);
   int  topo = evtSelected->MCEvtSelection(genParticles);


   histo1 = h_Jet;
   histo3 = h_Muon;
   histo4 = h_Ele;
   histo6 = h_MJet;
   histo7 = h_BJet;
   histo8 = h_WJet;
   histo9 = h_Top;

   //1. Isolation Leptons Selections
   ttMuon->muonAnalysis( muons, histo3 );
   //ttMuon->MuonTrigger( muons, triggers );

   ttEle->ElectronAnalysis(electrons, histo4);

   //2. build semi-Tt events

   if ( trigOn ) {
      Handle<edm::TriggerResults>  triggers;
      iEvent.getByLabel(triggerSrc, triggers);
      evtSelected->TriggerStudy( triggers, topo, 3 ,histo9 );
   }
   semiSol->BuildSemiTt(iEvent, iSetup, histo3, histo4, histo7, histo8, histo9 );
 
   //3. W+Jets analysis
   std::vector<const reco::Candidate*> isoMu = ttMuon->IsoMuonSelection( muons );
   //double minMET = ((*mets)[0]).et() ;

   if ( isoMu.size() == 1  ) {

      std::vector<pat::Jet> theJets = ttJet->JetSelection( jets, isoMu[0]->p4() ) ;

      ttJet->thirdETJetSpectrum( theJets, histo1 );
      //ttJet->thirdETJetSpectrum( genJets, histo1 );
      ttJet->dPhiMuJet( theJets, isoMu[0]->p4() , histo1 );
      
      std::vector<iReco> lepW;
      semiSol->recoW( isoMu, mets, lepW );
      if ( lepW.size() > 0 ) {
         int wl = MCMatching->matchLeptonicW( genParticles, lepW, histo9 );
         if ( wl != -1 ) ttJet->JetAndLepW( theJets, lepW[wl].p4, histo1 );
         //ttJet->JetAndLepW( theJets, lepW[wl].p4, histo1 );
      }
   }

   // looking for muon decay W
   /*
   LorentzVector wP4(0.,0.,0.,0.);
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       if ( abs((*it).pdgId()) != 24 ) continue;
       for (size_t q=0; q< (*it).numberOfDaughters(); q++) {
           const reco::Candidate *dau = (*it).daughter(q) ;
           if ( abs(dau->pdgId()) ==  13 && it->status() == 3 ) wP4 = it->p4() ;
       }
       LorentzVector zero(0.,0.,0.,0.);
       if ( wP4 != zero ) ttJet->JetAndLepW( jets, wP4, histo1 );
   }        
   */

   //4. Jet Studies
   ttJet->jetAnalysis(jets, histo1);
   //ttJet->JetTrigger( jets, triggers );

   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection(muons, histo3);
   std::vector<const pat::Jet*> tmpJets;
   std::vector<jmatch> mcwjets  = MCMatching->matchWJets(genParticles, jets, tmpJets, histo8, false);
   std::vector<jmatch> mcbjets  = MCMatching->matchbJets(genParticles, jets, tmpJets, histo7, false);
   ttJet->matchedWJetsAnalysis( mcwjets, isoMuons, histo8 );
   ttJet->matchedbJetsAnalysis( mcbjets, mcwjets, isoMuons, histo7 );
   ttJet->JetMatchedMuon( jets, muons, iEvent, iSetup, histo3, true );
   ttJet->bTagAnalysis(jets, histo7);
   if ( topo == 1) {
      ttJet->genJetInfo(genJets,genParticles, histo1, histo7, histo8);
   } 
   if ( pass && topo == 1) {
      ttJet->selectedWJetsAnalysis(jets,histo8);
   }

   //5. MET from PAT
   histo2 = h_MET;
   ttMET->metAnalysis(mets, iEvent, histo2);

   //6. Electron Studies
   std::vector<const reco::Candidate*> tempEle; 
   std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron(genParticles, electrons, tempEle, histo4, false);  
   ttEle->matchedElectronAnalysis( mcElectrons, histo4 );

   //7. photon studies
   histo5 = h_Gam;
   ttGam->PhotonAnalysis(photons, histo5);

}

//*************************************
//*        Utility Function           *
//*************************************

std::vector<int> TtAnalysis::findGrandMa(int momId, double eta, double phi, 
                                         Handle<std::vector<reco::GenParticle> > genParticles) {

 double dR = 99.;
 std::vector<int> sources; 
 for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
     if ( (*it).pdgId()!= momId ) continue;
     double dEta = (*it).eta() - eta ;
     double dPhi = (*it).phi() - phi ;
     double dR1 = sqrt( dEta*dEta + dPhi*dPhi );

     if ( dR1 < dR ) {
        dR = dR1;
        sources.clear();
        for (size_t i=0; i< (*it).numberOfMothers(); i++) {
            const reco::Candidate *grandma = (*it).mother(i);
            sources.push_back(grandma->pdgId()); 
            //cout<<" W/Z : "<<(*it).pdgId()<<" from "<<grandma->pdgId()<<" h:"<<(*it).eta()<<" f:"<<(*it).phi()<<endl;
        }
     } 
  }
  return sources;
}

// return the daughters 
hfPos TtAnalysis::findDaughter(int dauId, int momId, double mom_eta, double mom_phi,
                                         Handle<std::vector<reco::GenParticle> > genParticles) {
 double dR = 99.;
 int mom_status = -99;
 int dau_status = -99;
 hfPos decayPos(4);
 for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
     // find the daughter candidate
     if ( (*it).pdgId()!= dauId && (*it).status() == 1 ) continue;
    
     // match their mom
     for (size_t i=0; i< (*it).numberOfMothers(); i++) {
         const reco::Candidate *mom = (*it).mother(i) ;
         if ( mom->pdgId() != momId ) continue;
	 double dEta = mom->eta() - mom_eta ;
	 double dPhi = mom->phi() - mom_phi ;
	 double dR1 = sqrt( dEta*dEta + dPhi*dPhi );
         if ( dR1 < dR ) {
            dR  = dR1;
	    decayPos[0] = (*it).eta();
	    decayPos[1] = (*it).phi();
	    decayPos[2] = (*it).pt();
	    decayPos[3] = -1;
	    mom_status = mom->status();
	    dau_status = (*it).status();
         }
     }
    
  }

  //cout <<" The "<<momId<<" (h:"<<mom_eta<<" f:"<<mom_phi <<" st:"<<mom_status<<")";
  //cout <<" decays to status("<<dau_status<<") "<<dauId<<" h:"<<decayPos[0]<<" f:"<<decayPos[1]<<endl;

  return decayPos;

} 


//define this as a plug-in
DEFINE_FWK_MODULE(TtAnalysis);
