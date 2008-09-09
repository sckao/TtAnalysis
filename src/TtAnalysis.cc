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
static bool dmIncrease(const iReco t1, const iReco t2) { return ( t1.dm < t2.dm ); }

// constructors and destructor
using namespace edm;
using namespace std;
TtAnalysis::TtAnalysis(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  needTree          = iConfig.getUntrackedParameter<bool>   ("needTree");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  photonSrc         = iConfig.getParameter<edm::InputTag> ("photonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 

  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

  evtSelected = new TtEvtSelector();
  MCMatching  = new TtMCMatching();
  ttMuon      = new TtMuon();
  ttEle       = new TtElectron();
  ttGam       = new TtPhoton();
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );
  ttEff       = new TtEfficiency();

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

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

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
   /// calculate the selection efficiency
   histo9 = h_Top;
   ttEff->EventEfficiency( topo, pass, histo9 );

   //1. Isolation Leptons Selections
   histo3 = h_Muon;
   ttMuon->muonAnalysis( muons, histo3 );
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection(muons, histo3);
   histo4 = h_Ele;
   ttEle->ElectronAnalysis(electrons, histo4);
   std::vector<const reco::Candidate*> isoEle = ttEle->IsoEleSelection(electrons, histo4);
   //2. W and b jets selections
   std::vector<pat::Jet> selectedWJets = evtSelected->WJetSelection(jets);
   std::vector<pat::Jet> selectedbJets = evtSelected->bJetSelection(jets);
   
   //3. reconstruct semi-leptonic Tt events 
   std::vector<iReco> semiTt;
   if ( pass ) {
      semiTt = recoSemiLeptonicTtEvent(-1, selectedWJets, selectedbJets, isoMuons, isoEle, mets);
      ttEff->EventShape(topo, isoMuons.size(), isoEle.size(), selectedbJets.size(), selectedWJets.size(), histo9);
   }
   //4. gen-reco matching selection for leptons and jets 
   std::vector<const reco::Candidate*> mcMuons     = MCMatching->matchMuon(genParticles, muons);  
   std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron(genParticles, electrons);  
   std::vector<jmatch> mcwjets  = MCMatching->matchWJets(genParticles, jets);
   std::vector<jmatch> mcbjets  = MCMatching->matchbJets(genParticles, jets);

   //5. reconstruct semi-leptonic Tt events
   std::vector<pat::Jet> theWJets ; 
   for (size_t i =0; i< mcwjets.size(); i++ ){
       if ( mcwjets[i].hasMatched ) theWJets.push_back( mcwjets[i].trueJet ) ;
   }  
   std::vector<pat::Jet> thebJets ; 
   for (size_t i =0; i< mcbjets.size(); i++ ){
       if ( mcbjets[i].hasMatched ) thebJets.push_back( mcbjets[i].trueJet ) ;
   }
   std::vector<iReco> semiMCTt = recoSemiLeptonicTtEvent( topo, theWJets, thebJets, mcMuons, mcElectrons, mets);

   //6. look at top mass distribution 
   histo9 = h_Top;
   for(size_t i=0; i< semiMCTt.size(); i++ ){
      double tMass = getInvMass( semiMCTt[i].p4 ); 
      histo9->Fill9(tMass);
   }
   for(size_t i=0; i< semiTt.size(); i++ ){
      double tMass = getInvMass( semiTt[i].p4 ); 
      histo9->Fill9a(tMass);
   }

   //7. W studies 
   /// From Reco
   if ( pass && topo == 1) {
      bool isMuon     = ( isoMuons.size() ==1 && isoEle.size() ==0 ) ? true : false ;
      bool isElectron = ( isoMuons.size() ==0 && isoEle.size() ==1 ) ? true : false ;
      std::vector<const reco::Candidate*> leptons = ( isMuon && !isElectron  ) ? isoMuons : isoEle ;
      std::vector<iReco> LepWs;
      bool leptonic = recoW( leptons, mets ,LepWs );
      if (leptonic) dmSortRecoObjects( LepWs );
      for (size_t i =0; i< LepWs.size(); i++) {
          histo9->Fill9c( getInvMass( LepWs[i].p4 ) ); 
      }
      std::vector<iReco> HadWs;
      bool hadronic = recoW( selectedWJets , HadWs );
      if (hadronic) dmSortRecoObjects( HadWs );
      for (size_t i =0; i<HadWs.size(); i++) {
          histo9->Fill9e( getInvMass( HadWs[i].p4 ) ); 
      }
   }
   ///  From MC 
   if ( pass && topo == 1 ) {
      bool isMuon     = ( mcMuons.size() ==1 && mcElectrons.size() ==0 ) ? true : false ;
      bool isElectron = ( mcMuons.size() ==0 && mcElectrons.size() ==1 ) ? true : false ;
      std::vector<const reco::Candidate*> mclepton = ( isMuon && !isElectron  ) ? mcMuons : mcElectrons ;
      std::vector<iReco> mcLepWs;
      bool leptonic = recoW( mclepton, mets ,mcLepWs );
      if (leptonic) dmSortRecoObjects( mcLepWs );
      for (size_t i =0; i<mcLepWs.size(); i++) {
          histo9->Fill9b( getInvMass( mcLepWs[i].p4 ) ); 
      }
      std::vector<iReco> mcHadWs;
      bool hadronic = recoW( theWJets, mcHadWs );
      if (hadronic) dmSortRecoObjects( mcHadWs );
      for (size_t i =0; i<mcHadWs.size(); i++) {
          if (i > 0 ) break;
          histo9->Fill9d( getInvMass( mcHadWs[i].p4 )  ); 
      }
   }
 
   // 8. Jet Studies
   histo1 = h_Jet;
   ttJet->jetAnalysis(jets, iEvent, histo1);

   histo6 = h_MJet;
   histo7 = h_BJet;
   histo8 = h_WJet;
   ttJet->matchedWJetsAnalysis( mcwjets, isoMuons , histo8);
   ttJet->matchedbJetsAnalysis( mcbjets, mcwjets , isoMuons , histo7);
   ttJet->JetMatchedMuon( jets, muons, iEvent, iSetup, histo3, true);
   ttJet->bTagAnalysis(jets, histo7);

   if ( pass && topo == 1) {
      ttJet->selectedWJetsAnalysis(jets,histo8);
   }

   //9. MET from PAT
   histo2 = h_MET;
   ttMET->metAnalysis(mets, iEvent, histo2);
   histo5 = h_Gam;
   ttGam->PhotonAnalysis(photons, histo5);

   //10. Electron Studies
   ttEle->matchedElectronAnalysis( mcElectrons, histo4 );

   // Tt objects selection efficiency
   ttEff->JetEfficiency( selectedbJets, thebJets, histo9 );
   ttEff->JetEfficiency( selectedWJets, theWJets, histo9 );
   ttEff->IsoLeptonEfficiency( isoMuons, mcMuons, histo9 );
   ttEff->IsoLeptonEfficiency( isoEle, mcElectrons, histo9 );

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


// hadronic channel 
bool TtAnalysis::recoW( std::vector<pat::Jet> wjets, std::vector<iReco>& wCandidate ) {

    bool findW = false ; 
    if ( wjets.size() < 2 ) return findW ;

    for (size_t i=0; i< wjets.size(); i++ ) {
       for (size_t j =i+1; j < wjets.size(); j++ ) {
           double xW = wjets[i].p4().Px() + wjets[j].p4().Px() ;
           double yW = wjets[i].p4().Py() + wjets[j].p4().Py() ;
           double zW = wjets[i].p4().Pz() + wjets[j].p4().Pz() ;
           double EW = wjets[i].p4().E()  + wjets[j].p4().E() ;
	   double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;
	   if (massTest < 0 ) continue;
           iReco candW;
           LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
           std::pair<int,int> jetPair(i,j);
           candW.p4   = mW;
           candW.from = jetPair ;
           candW.dm   = fabs( sqrt(massTest) - 80.4 ) ;
           wCandidate.push_back( candW );
           findW = true;
       }
    }
    dmSortRecoObjects( wCandidate );
    return findW;
}

// semi-leptonic channel 
bool TtAnalysis::recoW( std::vector<const reco::Candidate*> lepton, Handle<std::vector<pat::MET> >  patMet, 
                        std::vector<iReco>& wCandidate ) {
    
    std::vector<pat::MET> met ;
    for (std::vector<pat::MET>::const_iterator m1 = patMet->begin(); m1 != patMet->end(); m1++){
        met.push_back( *m1 );     
    }

    bool findW = false ; 
    if ( lepton.size() != 1 || met.size() != 1) return findW ;

    // use w mass constrain to solve 4 momentum of W
    double xW = lepton[0]->p4().Px() + met[0].p4().Px() ;
    double yW = lepton[0]->p4().Py() + met[0].p4().Py() ;
    double nuPt2 = (met[0].p4().Px()*met[0].p4().Px()) + (met[0].p4().Py()*met[0].p4().Py()); 
 
    double zl = lepton[0]->p4().Pz() ;
    double El = lepton[0]->p4().E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    if ( (B*B) < (4.*A*C) ) return findW ;

    double nz1 = (-1.*B + sqrt(B*B - 4.*A*C )) / (2.*A) ;
    double nz2 = (-1.*B - sqrt(B*B - 4.*A*C )) / (2.*A) ;

    for (int i=1; i<3; i++ ) {
       double nz = 0.0;  
       if (i==1) nz =nz1;
       if (i==2) nz =nz2;
       double ENu2 = (met[0].p4().Px()*met[0].p4().Px()) + (met[0].p4().Py()*met[0].p4().Py()) + (nz*nz);
       double zW = lepton[0]->p4().Pz() + nz ;
       double EW = lepton[0]->p4().E() + sqrt(ENu2) ;
       double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

       if (massTest > 0. ) {
          iReco candW; 
          LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
          std::pair<int,int> lvPair(0,0);
          candW.p4   = mW;
          candW.from = lvPair;      
          candW.dm   = fabs( sqrt(massTest) - 80.4 );
          wCandidate.push_back( candW );
          findW = true;
       }
    }
    
    dmSortRecoObjects( wCandidate );
    return findW;
}

bool TtAnalysis::recoTop( std::vector<iReco> wCandidate, std::vector<pat::Jet> bjets, std::vector<iReco>& tCandidates ) {

  bool findt = false ;
  tCandidates.clear();
  if ( wCandidate.size() < 2 || bjets.size() < 2 ) return findt;

  int widx = 0;
  std::vector<iReco> bW1;
  std::vector<iReco> bW2;

  for(std::vector<iReco>::const_iterator i1 = wCandidate.begin(); i1 != wCandidate.end(); i1++ ){
     widx++;
     std::vector<iReco> bWpairList;  
     int bj = -1;
     for (std::vector<pat::Jet>::const_iterator i2 = bjets.begin(); i2 != bjets.end(); i2++ ){
         bj++;    
         double xT = i1->p4.Px() + (i2->p4()).Px() ;   
         double yT = i1->p4.Py() + (i2->p4()).Py() ;
         double zT = i1->p4.Pz() + (i2->p4()).Pz() ;
         double ET = i1->p4.E()  + (i2->p4()).E()  ;
	 LorentzVector mT = LorentzVector(xT,yT,zT,ET);
         double massTest =  (ET*ET) - (xT*xT) - (yT*yT) - (zT*zT) ;

	 if ( massTest <= 0. ) continue;
         iReco bWpair;
         bWpair.p4   = mT ;
         bWpair.from = std::pair<int,int>( widx-1 , bj );
         bWpair.dm   = fabs(sqrt(massTest) - 174.); 
	 bWpairList.push_back( bWpair );
	 findt = true;    
     }
     dmSortRecoObjects( bWpairList );
     if (widx == 1) bW1 = bWpairList;
     if (widx == 2) bW2 = bWpairList;
  }

  int sz = ( bW1.size() > bW2.size() ) ? bW1.size() : bW2.size() ;

  for(int i=0; i < sz; i++) {

     if (bW1[i].from.second != bW2[i].from.second ) {
        tCandidates.push_back(bW1[i]);
        tCandidates.push_back(bW2[i]);
        break;
     } else {
        if ( bW1[i].dm < bW2[i].dm ) {
           tCandidates.push_back(bW1[i]);
           tCandidates.push_back(bW2[i+1]);
           break;
        }
        if ( bW1[i].dm > bW2[i].dm ) {
           tCandidates.push_back(bW1[i+1]);
           tCandidates.push_back(bW2[i]);
           break;
        }
     }
     if ( tCandidates.size() > 1 ) break;

  }

  return findt;
}
 
std::vector<iReco> TtAnalysis::recoSemiLeptonicTtEvent(int topo, std::vector<pat::Jet> theWJets,
                            std::vector<pat::Jet> thebJets, std::vector<const reco::Candidate*> theMuons, 
                            std::vector<const reco::Candidate*> theElectrons, 
                            Handle<std::vector<pat::MET> >  mets ) {

   std::vector<iReco> TtCollection;
   if (topo != 1 && topo != -1) return TtCollection;

   // get the leptons 
   bool isMuon     = ( theMuons.size() ==1 && theElectrons.size() ==0 ) ? true : false ;
   bool isElectron = ( theMuons.size() ==0 && theElectrons.size() ==1 ) ? true : false ;
   std::vector<const reco::Candidate*> mclepton = ( isMuon && !isElectron  ) ? theMuons : theElectrons ;
   
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = recoW( theWJets, hadWs );
   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = recoW( mclepton, mets, lepWs );

   bool findtops = false ;
   if ( hadronic && leptonic ) {
      std::vector<iReco> wCandidates;
      std::vector<iReco> Tt_temp;
      for (size_t i=0; i< lepWs.size(); i++) {
         double dm_had = 999.;
         for (size_t j=0; j< hadWs.size(); j++) {
            wCandidates.push_back(lepWs[i]);
            wCandidates.push_back(hadWs[j]);
            findtops = recoTop( wCandidates ,thebJets, Tt_temp );
	    wCandidates.clear();
            if ( !findtops ) continue;
            // pick the best hadronic top
            if (Tt_temp[1].dm < dm_had) TtCollection = Tt_temp ;
         } 
      }
   }
   return TtCollection ;

}

double TtAnalysis::getInvMass( LorentzVector lv ) {

     double mom2  = (lv.Px()*lv.Px()) + (lv.Py()*lv.Py()) + (lv.Pz()*lv.Pz()) ;
     double mass2 = lv.E()*lv.E() - mom2;
     double mass = 0. ;
     if (mass2 < 0. ) mass2 = 0;
     mass = sqrt(mass2) ;
     return mass;

}


void TtAnalysis::dmSortRecoObjects( std::vector<iReco>& objCandidates ) {

   sort(objCandidates.begin(), objCandidates.end(), dmIncrease );

}
//define this as a plug-in
DEFINE_FWK_MODULE(TtAnalysis);
