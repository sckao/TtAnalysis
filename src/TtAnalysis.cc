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
typedef std::pair<double, LorentzVector> mtpair ;
static bool dmIncreasing(const mtpair t1, const mtpair t2) { return ( t1.first < t2.first ); }

// constructors and destructor
using namespace edm;
using namespace std;
TtAnalysis::TtAnalysis(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
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
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );

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

  t_Jet  = new NJet();
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
   delete ttMET;
   delete ttJet;

   //theFile->cd();
   t_Jet->Write();

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
  
   delete t_Jet;

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
   NJet  *jtree = 0;

   cout<<" ***** new Event start ***** "<<endl;

   // 0. select the semi-lep events and objects
   bool pass = evtSelected->eventSelection(muons, electrons, jets);
   int  topo = evtSelected->MCEvtSelection(genParticles);

   jtree = t_Jet;
   evtIt++;
   int eventId = evtIt + (iEvent.id().run()*100000) ;

   //2. Isolation Leptons 
   histo3 = h_Muon;
   ttMuon->muonAnalysis(muons, histo3, jtree, eventId);
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection(muons, histo3);
   histo4 = h_Ele;
   std::vector<const reco::Candidate*> isoEle = ttEle->IsoEleSelection(electrons, histo4);
   std::vector<pat::Jet> selectedWJets = evtSelected->bJetSelection(jets);
   std::vector<pat::Jet> selectedbJets = evtSelected->WJetSelection(jets);
   std::vector<LorentzVector> semiTt;
   if ( pass ) {
      semiTt = recoSemiLeptonicTtEvent(-1, selectedWJets, selectedbJets, isoMuons, isoEle, mets);
   }

   // 1. Jet Sector
   histo1 = h_Jet;
   ttJet->jetAnalysis(jets, iEvent, histo1, jtree, eventId);

   /// 1b) gen-reco matching 
   histo6 = h_MJet;
   histo7 = h_BJet;
   std::vector<const reco::Candidate*>  mcMuons     = MCMatching->matchMuon(genParticles, muons);  
   std::vector<const reco::Candidate*>  mcElectrons = MCMatching->matchElectron(genParticles, electrons);  
   std::vector<jmatch> mcwjets    = MCMatching->matchWJets(genParticles, jets);
   std::vector<jmatch> mcbjets    = MCMatching->matchbJets(genParticles, jets);

   std::vector<pat::Jet> theWJets ; 
   for (size_t i =0; i< mcwjets.size(); i++ ){
       theWJets.push_back( mcwjets[i].trueJet ) ;
   }  
   std::vector<pat::Jet> thebJets ; 
   for (size_t i =0; i< mcbjets.size(); i++ ){
       thebJets.push_back( mcbjets[i].trueJet ) ;
   }  
   std::vector<LorentzVector> semiMCTt = recoSemiLeptonicTtEvent(topo, theWJets, thebJets, mcMuons, mcElectrons, mets);
   histo9 = h_Top;
   for(size_t i=0; i< semiMCTt.size(); i++ ){
      double tMass = getInvMass( semiMCTt[i] ); 
      histo9->Fill9(tMass);
   }
   for(size_t i=0; i< semiTt.size(); i++ ){
      double tMass = getInvMass( semiTt[i] ); 
      histo9->Fill9a(tMass);
   }

   ttJet->matchedWJetsAnalysis( mcwjets, isoMuons , histo6);
   ttJet->matchedbJetsAnalysis( mcbjets, mcwjets , isoMuons , histo7);
   cout<<" ======== test another matching ======= "<<endl;
   ttJet->JetMatchedMuon( jets, muons, iEvent, iSetup, histo3, true);

   //  1c) jets selection
   ttJet->bTagAnalysis(jets, histo7);
   histo8 = h_WJet;
   if ( pass && topo == 1) {
      ttJet->selectedWJetsAnalysis(jets,histo8);
   }

   // look at the W for qq 
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       if ( abs((*it).pdgId()) != 24) continue;
       bool WfromT = false;
       for (size_t q=0; q< (*it).numberOfMothers(); q++) {
           const reco::Candidate *mom = (*it).mother(q) ;
           if ( abs(mom->pdgId()) != 6 ) continue;
           WfromT = true ;
       }
       if ( !WfromT ) continue;

       std::vector<LorentzVector> qm ;
       for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
           const reco::Candidate *dau = (*it).daughter(q) ;
           if( abs(dau->pdgId()) > 6 ) continue;
           qm.push_back( dau->p4() );
	   jtree->FillBgen( eventId, dau->pdgId(), dau->eta(), dau->phi(), dau->energy(), dau->pt() );
           //cout<<"quarks:"<<dau->pdgId() <<" h: "<<dau->eta()<<" f:"<<dau->phi()<<" It:"<<eventId;
           //cout<<" E:"<<dau->energy()<< endl;
       }
       if ( qm.size() == 2 ) {
          LorentzVector pW = ttJet->findW( qm[0], qm[1]);
	  double momW = sqrt( pW.Px()*pW.Px() + pW.Py()*pW.Py() + pW.Pz()*pW.Pz() );
	  double massW = sqrt ( pW.E()*pW.E() - pW.Px()*pW.Px() - pW.Py()*pW.Py() - pW.Pz()*pW.Pz() );
	  histo1->Fill1g( massW, momW );
       }
   }
   //cout<<" "<<endl;
   //cout<<"  ^^^^^^^ end of w reco ^^^^^^^^ "<<endl;

   // MET from PAT
   histo2 = h_MET;
   ttMET->metAnalysis(mets, iEvent, histo2, jtree, eventId);
   ttMuon->muonAnalysis(muons, histo3, jtree, eventId);

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


LorentzVector TtAnalysis::findTop(LorentzVector qm1, LorentzVector qm2 ) {
     
     //cout<<" ------------------ reco Top ------------------- "<<endl;
     double p1 = qm1.P();
     double p2 = qm2.P();
     //double p1p2 = qm1.Px()*qm2.Px() + qm1.Py()*qm2.Py() + qm1.Pz()*qm2.Pz() ;
     //double m1 = sqrt( (qm1.E()*qm1.E()) - (qm1.P()*qm1.P()) ) ;
     //double m2 = sqrt( (qm2.E()*qm2.E()) - (qm2.P()*qm2.P()) ) ;
     //double massT = sqrt( (m1*m1) + (m2*m2) + (2.0*p1*p2) - (2.0*p1p2) ) ;
     //double angleqq = acos( p1p2/(p1*p2) ) ;
     //cout<<" q1("<<p1<<") q2("<<p2<<") dP= "<< fabs(p1-p2) <<",open angle of Wb = "<< angleqq << endl;
     //cout<<" mass of Top = "<< massT <<endl;

     double xT = qm1.Px() + qm2.Px() ;
     double yT = qm1.Py() + qm2.Py() ;
     double zT = qm1.Pz() + qm2.Pz() ;
     double ET = p1 + p2 ;   
     LorentzVector mT = LorentzVector(xT,yT,zT,ET) ;
     //cout<<" pT: "<<xT<<","<<yT<<","<<zT<<","<<ET<<endl;
     //cout<<" mT: "<<mT.Px()<<","<<mT.Py()<<","<<mT.Pz()<<","<<mT.E()<<endl;

     return mT;
}

// hadronic channel 
bool TtAnalysis::recoW( std::vector<pat::Jet> wjets, std::vector<LorentzVector>& wCandidate ) {

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
           LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
           wCandidate.push_back( mW );
           findW = true;
       }
    }
    return findW;
}

// semi-leptonic channel 
bool TtAnalysis::recoW( std::vector<const reco::Candidate*> lepton, Handle<std::vector<pat::MET> >  patMet, 
                        std::vector<LorentzVector>& wCandidate ) {
    
    std::vector<pat::MET> met ;
    for (std::vector<pat::MET>::const_iterator m1 = patMet->begin(); m1 != patMet->end(); m1++){
        met.push_back( *m1 );     
    }

    bool findW = false ; 
    if ( lepton.size() != 1 || met.size() != 1) return findW ;

    double xW = lepton[0]->p4().Px() + met[0].p4().Px() ;
    double yW = lepton[0]->p4().Py() + met[0].p4().Py() ;
    double zW = lepton[0]->p4().Pz() + met[0].p4().Pz() ;
    double EW = lepton[0]->p4().E()  + met[0].p4().E() ;
    double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;
    if (massTest > 0. ) { 
       LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
       wCandidate.push_back( mW );
       findW = true;
    }
    return findW;
}

bool TtAnalysis::recoTop( std::vector<LorentzVector> wCandidate, std::vector<pat::Jet> bjets, std::vector<LorentzVector>& tCandidates ) {

  bool findt = false ;
  for(std::vector<LorentzVector>::const_iterator i1 = wCandidate.begin(); i1 != wCandidate.end(); i1++ ){ 
     for (std::vector<pat::Jet>::const_iterator i2 = bjets.begin(); i2 != bjets.end(); i2++ ){
         double xT = i1->Px() + (i2->p4()).Px() ;   
         double yT = i1->Py() + (i2->p4()).Py() ;
         double zT = i1->Pz() + (i2->p4()).Pz() ;
         double ET = i1->E()  + (i2->p4()).E()  ;
         double massTest =  (ET*ET) - (xT*xT) - (yT*yT) - (zT*zT) ;
         if (massTest > 0. ) {
            LorentzVector mT = LorentzVector(xT,yT,zT,ET);
            tCandidates.push_back( mT );
            findt = true;        
         }
     }
  }
  return findt;
}
 
std::vector<LorentzVector> TtAnalysis::recoSemiLeptonicTtEvent(int topo, std::vector<pat::Jet> theWJets,
                            std::vector<pat::Jet> thebJets, std::vector<const reco::Candidate*> mcMuons, 
                            std::vector<const reco::Candidate*> mcElectrons, 
                            Handle<std::vector<pat::MET> >  mets ) {

   std::vector<LorentzVector> TtCollection;
   if (topo != 1 && topo != -1) return TtCollection;

   std::vector<LorentzVector> wCadidates;
   bool isMuon     = ( mcMuons.size() > 0 && mcElectrons.size() ==0 ) ? true : false ;
   bool isElectron = ( mcMuons.size() ==0 && mcElectrons.size() > 0 ) ? true : false ;
   std::vector<const reco::Candidate*> mclepton = ( isMuon && !isElectron  ) ? mcMuons : mcElectrons ;
   bool hadronic = recoW( theWJets, wCadidates );
   bool leptonic = leptonic = recoW( mclepton, mets ,wCadidates );
   bool findtops = false ;
   if (hadronic && leptonic ) {
      findtops = recoTop( wCadidates,thebJets, TtCollection );
   }

   std::vector<mtpair> tlist;
   for (size_t i=0; i<TtCollection.size(); i++) {

       double tmass = getInvMass(TtCollection[i]) ;  
       double dm = fabs(tmass - 174.) ;
       mtpair tpair(dm , TtCollection[i]);  
       tlist.push_back( tpair );
   }
   sort(tlist.begin(), tlist.end(), dmIncreasing );
    
   TtCollection.clear();
   for(size_t i=0; i<tlist.size() ; i++) {
      if (i > 2) continue;
      TtCollection.push_back(tlist[i].second); 
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
//define this as a plug-in
DEFINE_FWK_MODULE(TtAnalysis);
