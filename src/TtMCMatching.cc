// -*- C++ -*-
//
// Package:    TtMCMatching
// Class:      TtMCMatching
// 
/**\class TtMCMatching TtMCMatching.cc PhysicsTools/TtAnalysis/src/TtMCMatching.cc

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
#include "TtMCMatching.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"

//#include "FWCore/Framework/interface/Event.h"

//#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//


// constructors and destructor
using namespace edm;
using namespace std;
//TtMCMatching::TtMCMatching(const edm::ParameterSet& iConfig)
TtMCMatching::TtMCMatching()
{
   //now do what ever initialization is needed
  /*
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 

  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");

   */
  tools           = new TtTools();
}


TtMCMatching::~TtMCMatching()
{
   delete tools;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   // if (debug) cout << "[TtMCMatching Analysis] Destructor called" << endl;
}

//
// member functions
//
static bool PtDecreasing(const reco::Particle s1, const reco::Particle s2) { return ( s1.pt() > s2.pt() ); }
static bool PidDecreasing( jmatch s1, jmatch s2) { return ( s1.MomIdx > s2.MomIdx ); }
static bool dRIncreasing( IDPair s1, IDPair s2) { return ( s1.second < s2.second ); }

void TtMCMatching::MCTreeFeeder(Handle<std::vector<reco::GenParticle> > genParticles , NJet* jtree, int eventId ){

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
       }
   } 

} 

std::vector<jmatch> TtMCMatching::matchWJets( Handle<std::vector<reco::GenParticle> > genParticles,
                               std::vector<const pat::Jet*> selectedWJets, HTOP8* histo8, bool fillhisto ) {

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> jetMom ;
   for (int i=1; i<5; i++) {
       std::vector<reco::Particle> tmpMom = ttDecay(genParticles, i) ;
       for (std::vector<reco::Particle>::iterator it= tmpMom.begin(); it!= tmpMom.end(); it++){
           jetMom.push_back( *it );
       }
   }
 
   // Loop all pat jets
   // get 4 collection w.r.t. possible 4 jet mom(quarks from W)
   // each collection contains candidate jets
   std::vector<int> q1; // q+ from w+
   std::vector<int> q2; // q- from w+
   std::vector<int> q3; // q+ from w-
   std::vector<int> q4; // q- from w-
   for (size_t i=0; i < selectedWJets.size(); i++){ 
       double dR0    = 99. ;
       double ptRes0 = 99. ;
       int wj = -1 ;
       for (size_t j=0; j < jetMom.size(); j++ ) {
           bool matched = matchingGeneral( jetMom[j].p4() , selectedWJets[i]->p4() , dR0 , ptRes0 ) ;
           if ( matched )  wj = static_cast<int>(j) ;
       }
       if ( wj!= -1 && (jetMom[wj].pdgId()== 2 || jetMom[wj].pdgId()== 4) ) q1.push_back( static_cast<int>(i) ) ;
       if ( wj!= -1 && (jetMom[wj].pdgId()==-1 || jetMom[wj].pdgId()==-3) ) q2.push_back( static_cast<int>(i) ) ;
       if ( wj!= -1 && (jetMom[wj].pdgId()== 1 || jetMom[wj].pdgId()== 3) ) q3.push_back( static_cast<int>(i) ) ;
       if ( wj!= -1 && (jetMom[wj].pdgId()==-2 || jetMom[wj].pdgId()==-4) ) q4.push_back( static_cast<int>(i) ) ;
   }
  
   // loop 2or4 collections to find the best matched jet 
   std::vector<jmatch> matchedJets ;
   std::vector<int> wtemp;   
   for (size_t j=0; j < jetMom.size(); j++  ) {
       
       if (jetMom[j].pdgId() ==  2 || jetMom[j].pdgId()==  4) wtemp = q1 ;
       if (jetMom[j].pdgId() == -1 || jetMom[j].pdgId()== -3) wtemp = q2 ;
       if (jetMom[j].pdgId() ==  1 || jetMom[j].pdgId()==  3) wtemp = q3 ;
       if (jetMom[j].pdgId() == -2 || jetMom[j].pdgId()== -4) wtemp = q4 ;

       double dR0 = 99.;
       double ptRes0 =99.;
       int wj = -1;
       for(std::vector<int>::iterator it=wtemp.begin(); it != wtemp.end(); it++ ){
          bool matched = matchingGeneral( jetMom[j].p4() , selectedWJets[*it]->p4() , dR0 , ptRes0 ); 
          if (matched) wj = *it ;
       }

       if (wj != -1) {
          jmatch theWjet;
	  theWjet.trueJet= selectedWJets[wj] ;
	  theWjet.mom    = jetMom[j] ;
	  theWjet.MomIdx = jetMom[j].pdgId() ;
	  theWjet.res_P  = (selectedWJets[wj]->pt() /jetMom[j].pt()) - 1. ;
	  theWjet.sumP4  = selectedWJets[wj]->p4() ;
	  theWjet.hasMatched = true ;
          matchedJets.push_back( theWjet );
          if ( fillhisto ) histo8->Fill8h( jetMom[j].eta(), jetMom[j].phi(), selectedWJets[wj]->eta(), selectedWJets[wj]->phi()
                                         , jetMom[j].pt(), theWjet.res_P  );
       }

   }
   return matchedJets ;
}

// general matching for w jets and bjets 
std::vector<jmatch> TtMCMatching::matchJets( Handle<std::vector<reco::GenParticle> > genParticles,
                                           std::vector<const pat::Jet*> selectedJets,
                                           HTOP7* histo7, HTOP8* histo8, bool fillhisto ) {

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> jetMom ;
   for (int i=1; i<6; i++) {
       std::vector<reco::Particle> tmpMom = ttDecay(genParticles, i) ;
       for (std::vector<reco::Particle>::iterator it= tmpMom.begin(); it!= tmpMom.end(); it++){
           jetMom.push_back( *it );
           //if (fillhisto) cout<<"tmpMom sz"<<tmpMom.size() << " jmom:"<< it->pdgId() <<endl;
       }
   }
   //if (fillhisto) cout<<"  === jetMom size == "<<jetMom.size() <<endl;

   // make a map for jet index and dR with all parton
   std::vector< std::vector<IDPair> > Qs;
   for( size_t i=0; i< jetMom.size(); i++ ) {
      std::vector<IDPair> js ;
      for( size_t j=0; j<selectedJets.size(); j++ ){
         double dR = tools->getdR( jetMom[i].p4(), selectedJets[j]->p4() ) ;
         IDPair jId_dR = make_pair(j,dR);
         js.push_back( jId_dR );
      }
      sort( js.begin(), js.end(), dRIncreasing );
      Qs.push_back(js);
   }

   // check the best matching for each parton
   std::vector<int>  qMatch( jetMom.size() , -1 );
   std::vector<bool> usedJet( selectedJets.size(), false );

   // loop index m => the order of match-sorting
   for( size_t m=0; m<selectedJets.size(); m++ ) {

      for( size_t j=0; j<selectedJets.size(); j++ ) {

         if (usedJet[j]) continue;
         // look at the best matched jet for each parton
         std::vector<int> qCand;
         for( size_t i=0; i< Qs.size(); i++) {
            if ( qMatch[i] != -1 ) continue;
            std::vector<IDPair> JdR = Qs[i];
            if( JdR[m].first == static_cast<int>(j) ) qCand.push_back( i ) ;
         }

         double dR0    = 9. ;
         double ptRes0 = 9. ;
         int wj = -1 ;
         for (size_t k=0; k < qCand.size(); k++ ) {
             bool matched = matchingGeneral( jetMom[ qCand[k] ].p4() , selectedJets[j]->p4() , dR0 , ptRes0 ) ;
             if ( matched )  wj = qCand[k] ;
         }
         if ( wj !=  -1 ) {
            qMatch[wj] = j ;
            usedJet[j] = true;

         }

      }

   }

   cout<<" ===== new matching ====== "<<endl;
   std::vector<jmatch> matchedJets ;
   for(size_t k=0; k < qMatch.size(); k++ ) {
      cout<<" q:"<< jetMom[k].pdgId() <<"  matchedJ:"<< qMatch[k] <<endl;
      int s = qMatch[k] ;
      if ( s == -1 )continue;
      jmatch theTjet;
      theTjet.trueJet= selectedJets[s] ;
      theTjet.mom    = jetMom[k] ;
      theTjet.MomIdx = jetMom[k].pdgId() ;
      theTjet.res_P  = ( selectedJets[s]->pt() /jetMom[k].pt()) - 1. ;
      theTjet.sumP4  = selectedJets[s]->p4() ;
      theTjet.hasMatched = true ;
      matchedJets.push_back( theTjet );
      if ( fillhisto && abs(theTjet.MomIdx) == 5  ) histo7->Fill7h( jetMom[k].eta(), jetMom[k].phi(),
		      selectedJets[s]->eta(), selectedJets[s]->phi(), jetMom[k].pt(), theTjet.res_P  );
      if ( fillhisto && abs(theTjet.MomIdx)  < 5  ) histo8->Fill8h( jetMom[k].eta(), jetMom[k].phi(),
		      selectedJets[s]->eta(), selectedJets[s]->phi(), jetMom[k].pt(), theTjet.res_P  );
   }
   cout<<" ======================== "<<endl;

   sort( matchedJets.begin(), matchedJets.end(), PidDecreasing );
   return matchedJets ;

}


std::vector<jmatch> TtMCMatching::matchbJets( Handle<std::vector<reco::GenParticle> > genParticles,
                                              std::vector<const pat::Jet*> selectedbJets ,HTOP7* histo7, bool fillhisto ) {
   // Accumulate the b quark from t
   std::vector<reco::Particle> jetMom = ttDecay(genParticles, 5) ;
 
   // Loop all selected reco bjets => rely on b-tagging
   // get 2 collections w.r.t. possible 2 bjet mom( b and bbar )
   // each collection contains candidate jets to avoid double using

   // 1. classifying every reco jet -> bp or bn are the list of candidates for b or bbar
   std::vector<int> bp;
   std::vector<int> bn;
   for (size_t i=0; i < selectedbJets.size(); i++){ 
       // exclude the jet with no charged track
       //edm::RefVector<reco::TrackCollection>  assTk = (*jets)[i].associatedTracks() ;
       //if ( assTk.size()== 0  ) continue;
       
       double dR0 = 99.;
       double ptRes0 =99.;
       int bj = -1;
       for (size_t j=0; j < jetMom.size(); j++ ) {
           bool matched = matchingGeneral( jetMom[j].p4() , selectedbJets[i]->p4() , dR0 , ptRes0 );
           if ( matched )  bj = static_cast<int>(j) ;
       }

       if ( bj!= -1 && jetMom[bj].pdgId() ==  5 ) bp.push_back(  static_cast<int>(i) ) ;
       if ( bj!= -1 && jetMom[bj].pdgId() == -5 ) bn.push_back(  static_cast<int>(i) ) ;
      
   }

   // loop bp/bn to find the best matching ; 
   std::vector<jmatch> matchedJets ;
   std::vector<int> btemp;   
   for (size_t j=0; j < jetMom.size(); j++  ) {
       
       if (jetMom[j].pdgId() ==  5) btemp = bp ;
       if (jetMom[j].pdgId() == -5) btemp = bn ;

       double dR0 = 99.;
       double ptRes0 =99.;
       int bj = -1;
       for(std::vector<int>::iterator it=btemp.begin(); it != btemp.end(); it++  ){
          bool matched = matchingGeneral( jetMom[j].p4() , selectedbJets[*it]->p4() , dR0 , ptRes0 ); 
          if (matched) bj = *it ;
       }

       if (bj != -1) {
          jmatch thebjet;
	  thebjet.trueJet= selectedbJets[bj] ;
	  thebjet.mom    = jetMom[j] ;
	  thebjet.MomIdx = jetMom[j].pdgId() ;
	  thebjet.res_P  = ( selectedbJets[bj]->pt() /jetMom[j].pt()) - 1. ;
	  thebjet.sumP4  = selectedbJets[bj]->p4() ;
	  thebjet.hasMatched = true ;
          matchedJets.push_back( thebjet );
          if ( fillhisto ) histo7->Fill7h( jetMom[j].eta(), jetMom[j].phi(), selectedbJets[bj]->eta(), selectedbJets[bj]->phi()
                                         , jetMom[j].pt(), thebjet.res_P );
       }


   }
   return matchedJets ;

}

/*
std::vector<const reco::Candidate*> TtMCMatching::matchGenJet(Handle<std::vector<reco::GenParticle> > genParticles, Handle<std::vector<reco::GenJet> > genJets, HTOP6* histo6 ) {

   // Accumulate the b quark from t
   std::vector<iTt> ttobjs = TtObjects(genParticles);
   std::vector<reco::Particle> jetMom = ttPartons(genParticles, 5) ;

}
*/

std::vector<const reco::Candidate*> TtMCMatching::matchMuon( Handle<std::vector<reco::GenParticle> > genParticles,
                          std::vector<const reco::Candidate*> isoMuons ,HTOP3* histo3, bool fillhisto ){

   // Accumulate the leptonic dauaghters from W
   std::vector<reco::Particle> mcMuon = ttDecay(genParticles, 13) ;

   // find the matched muon 
   std::vector<int> matchList;
   for (size_t i=0; i< mcMuon.size(); i++) {

      // matching selected isoMuons with parton
      double dR1    = 9.;
      double ptRes1 = 9.;
      int usedMuon = -1;
      for (size_t j=0; j < isoMuons.size(); j++ ) {
          bool matched = matchingGeneral( mcMuon[i].p4(), isoMuons[j]->p4(), dR1, ptRes1 );
          if ( matched && isoMuons[j]->charge() == mcMuon[i].charge() ) usedMuon = static_cast<int>(j);
      }
      if (usedMuon != -1) {
         matchList.push_back( usedMuon ); 
         double ptRes2 = (isoMuons[usedMuon]->pt()/mcMuon[i].pt()) - 1. ;
         if (fillhisto) histo3->Fill3b( mcMuon[i].eta(),  mcMuon[i].phi(), isoMuons[usedMuon]->eta(), isoMuons[usedMuon]->phi()
                                       ,mcMuon[i].pt(), ptRes2 );
      }

   }
   // using one more loops to check whether double counting 
   std::vector<const reco::Candidate*> matchedMuon ;

   std::vector<bool> used(isoMuons.size(), false  );

   for (std::vector<int>::const_iterator it= matchList.begin(); it != matchList.end(); it++) { 
       if ( used[*it] ) continue;
       matchedMuon.push_back( isoMuons[*it] );  
       used[*it] = true;
   }

   return matchedMuon ;
}

std::vector<const reco::Candidate*> TtMCMatching::matchMuon( Handle<std::vector<reco::GenParticle> > genParticles,
           Handle<std::vector<pat::Muon> > muons, HTOP3* histo3, bool fillhisto ){

   // Accumulate the leptonic dauaghters from W
   std::vector<reco::Particle> mcMuon = ttDecay(genParticles, 13) ;

   // find the matched muon 
   std::vector<int> matchList;
   for (size_t i=0; i< mcMuon.size(); i++) {

      // matching pat::muons with parton
      double dR0    = 9.;
      double ptRes0 = 9.;
      int theMuon = -1;
      for (size_t j=0; j < muons->size(); j++ ) {
          bool matched = matchingGeneral( (*muons)[j].p4() , mcMuon[i].p4(), dR0, ptRes0 );
          if ( matched && (*muons)[j].charge() == mcMuon[i].charge() ) theMuon = static_cast<int>(j);
      }
      if (theMuon != -1) {
         matchList.push_back( theMuon ); 
         double ptRes2 = ((*muons)[theMuon].pt()/mcMuon[i].pt()) - 1. ;
         if (fillhisto) histo3->Fill3b( mcMuon[i].eta(),  mcMuon[i].phi(), (*muons)[theMuon].eta(), (*muons)[theMuon].phi()
                                       ,mcMuon[i].pt(), ptRes2 );
      }
    
   }
   // using one more loops to check whether double counting 
   std::vector<const reco::Candidate*> matchedMuon ;

   std::vector<bool> used;
   for (size_t k=0; k < muons->size(); k++ ) {
       used.push_back(false);
   }
   for (std::vector<int>::const_iterator it= matchList.begin(); it != matchList.end(); it++) { 
       if ( used[*it] ) continue;
       matchedMuon.push_back( &(*muons)[*it] );  
       used[*it] = true;
   }

   if ( matchedMuon.size() > 2 ) cout<<" WRONG MATCHING "<<endl;

   return matchedMuon ;
}

std::vector<const reco::Candidate*> TtMCMatching::matchElectron(Handle<std::vector<reco::GenParticle> > genParticles,
                                        Handle<std::vector<pat::Electron> > electrons, HTOP4* histo4, bool fillhisto ){

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> mcElectron = ttDecay(genParticles, 11) ;

   // find the matched electron 
   std::vector<int> matchList;
   for (size_t i=0; i< mcElectron.size(); i++) {

      // matching pat::electrons with parton
      double dR0 = 99.;
      double ptRes0 =99.;
      int theElectron = -1;
      for (size_t j=0; j < electrons->size(); j++ ) {
          bool matched = matchingGeneral( (*electrons)[j].p4() , mcElectron[i].p4(), dR0, ptRes0 );
          if ( matched && (*electrons)[j].charge() == mcElectron[i].charge() ) theElectron = static_cast<int>(j);
      }
      if (theElectron != -1) {
         matchList.push_back( theElectron ); 
         double ptRes2 = ((*electrons)[theElectron].pt()/mcElectron[i].pt()) - 1. ;
         if (fillhisto) histo4->Fill4g( mcElectron[i].eta(), mcElectron[i].phi(), (*electrons)[theElectron].eta(), 
                                       (*electrons)[theElectron].phi() ,mcElectron[i].pt(), ptRes2 );
      }
   }

   // using one more loops to check whether double counting 
   std::vector<bool> used;
   std::vector<const reco::Candidate*> matchedElectron ;
   for (size_t k=0; k < electrons->size(); k++ ) {
       used.push_back(false);
   }
   for (std::vector<int>::const_iterator it= matchList.begin(); it != matchList.end(); it++) { 
       if ( used[*it] ) continue;
       matchedElectron.push_back( &(*electrons)[*it] );  
       used[*it] = true;
   }

   //std::vector<pat::Electron> matchedElectron ;
   if ( matchedElectron.size() > 2 ) cout<<" WRONG MATCHING "<<endl;

   return matchedElectron ;
}

std::vector<const reco::Candidate*> TtMCMatching::matchElectron(Handle<std::vector<reco::GenParticle> > genParticles,
                                           std::vector<const reco::Candidate*> isoEle, HTOP4* histo4, bool fillhisto ){

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> mcElectron = ttDecay(genParticles, 11) ;

   // find the matched electron 
   std::vector<int> matchList;
   for (size_t i=0; i< mcElectron.size(); i++) {

      // matching pat::electrons with parton
      double dR0 = 99.;
      double ptRes0 =99.;
      int theElectron = -1;
      for (size_t j=0; j < isoEle.size(); j++ ) {
          bool matched = matchingGeneral( isoEle[j]->p4() , mcElectron[i].p4(), dR0, ptRes0 );
          if ( matched &&  isoEle[j]->charge() == mcElectron[i].charge() ) theElectron = static_cast<int>(j);
      }
      if (theElectron != -1) {
         matchList.push_back( theElectron ); 
         double ptRes2 = (isoEle[theElectron]->pt()/mcElectron[i].pt()) - 1. ;
         if (fillhisto) histo4->Fill4g( mcElectron[i].eta(), mcElectron[i].phi(), isoEle[theElectron]->eta(), 
                                       isoEle[theElectron]->phi() ,mcElectron[i].pt(), ptRes2 );
      }
    
   }
   // using one more loops to check whether double counting 
   std::vector<bool> used;
   std::vector<const reco::Candidate*> matchedElectron ;
   for (size_t k=0; k < isoEle.size(); k++ ) {
       used.push_back(false);
   }
   for (std::vector<int>::const_iterator it= matchList.begin(); it != matchList.end(); it++) { 
       if ( used[*it] ) continue;
       matchedElectron.push_back( isoEle[*it] );  
       used[*it] = true;
   }

   //std::vector<pat::Electron> matchedElectron ;
   if ( matchedElectron.size() > 2 ) cout<<" WRONG MATCHING "<<endl;

   return matchedElectron ;
}


int TtMCMatching::matchLeptonicW( Handle<std::vector<reco::GenParticle> > genParticles, std::vector<iReco> wSolutions ){

   int wl = -1;
   std::vector<reco::Particle> nu_m = ttDecay( genParticles, 14 ) ;
   std::vector<reco::Particle> nu_e = ttDecay( genParticles, 12 ) ;
   std::vector<reco::Particle> neu;
   if ( nu_m.size() == 1 && nu_e.size() == 0 ) neu = nu_m ;
   if ( nu_m.size() == 0 && nu_e.size() == 1 ) neu = nu_e ;
   if ( neu.size() == 0 ) return wl;
   cout<<" nu_m:"<< nu_m.size() <<" nu_e:"<<nu_e.size()<<" p4:"<< neu[0].p4()<<endl;

   double dR0 = 9.;
   double ptRes0 = 9.;
   for (size_t j=0; j < wSolutions.size(); j++ ) {
       iParton reco_neu = wSolutions[j].q4v[1] ;
       bool matched = matchingGeneral( neu[0].p4() , reco_neu.second , dR0 , ptRes0 );
       if ( matched )  wl = static_cast<int>(j) ;
   }
   cout<<" selected w: "<< wl <<endl;
   return wl;

}


int TtMCMatching::matchLeptonicW( Handle<std::vector<reco::GenParticle> > genParticles, std::vector<iReco> wSolutions, HTOP6* histo6 ){

   int wl = -1;
   std::vector<reco::Particle> nu_m = ttDecay( genParticles, 14 ) ;
   std::vector<reco::Particle> nu_e = ttDecay( genParticles, 12 ) ;
   std::vector<reco::Particle> neu;
   if ( nu_m.size() == 1 && nu_e.size() == 0 ) neu = nu_m ;
   if ( nu_m.size() == 0 && nu_e.size() == 1 ) neu = nu_e ;
   if ( neu.size() == 0 ) return wl;

   double dR0 = 9.;
   double ptRes0 = 9.;
   for (size_t j=0; j < wSolutions.size(); j++ ) {
       iParton reco_neu = wSolutions[j].q4v[1] ;
       bool matched = matchingGeneral( neu[0].p4() , reco_neu.second , dR0 , ptRes0 );
       if ( matched )  wl = static_cast<int>(j) ;
   }

   if ( wl != -1 ) {
      LorentzVector lv = wSolutions[wl].p4 ;

      double ptRes = (lv.Pt()/neu[0].pt()) - 1 ;
      double pzRes = (lv.Pz()/neu[0].p4().Pz()) - 1 ;

      histo6->Fill6d( dR0, ptRes, pzRes );
   }
   return wl;
}

bool TtMCMatching::matchingGeneral( LorentzVector theP4 , iTt ttObj, double& dR0, double& ptRes0 ){

     bool betterMatch = false;

     double pt1 = sqrt( (theP4.Px()*theP4.Px())  + (theP4.Py()*theP4.Py()) );
   
     double dR = tools->getdR( theP4, ttObj.p4 );
     double ptRes = fabs( (ttObj.pt/pt1 ) - 1.) ;

     if ( (dR - dR0) < 0.1 && (dR < 0.5) && ( ptRes < ptRes0 ) ) {
        betterMatch = true ;
        dR0    = dR ;
        ptRes0 = ptRes;
     }
     if ( (dR - dR0) < 0.  && dR >= 0.5 ) {
        betterMatch = true ;
        dR0    = dR ;
        ptRes0 = ptRes;
     }

     return betterMatch ;
}

bool TtMCMatching::matchingGeneral( LorentzVector aP4 ,  LorentzVector bP4, double& dR0, double& ptRes0 ){

     bool betterMatch = false;

     double pt1 = sqrt( (aP4.Px()*aP4.Px())  + (aP4.Py()*aP4.Py()) );
     double pt2 = sqrt( (bP4.Px()*bP4.Px())  + (bP4.Py()*bP4.Py()) );

     double dR = tools->getdR( aP4, bP4 );
     // 1/pt resolution
     double ptRes = fabs( (pt2/pt1 ) - 1.) ;

     if ( (dR - dR0) < 0.05 && (dR < 0.5) && ( ptRes < ptRes0 ) ) {
        betterMatch = true ;
        dR0    = dR ;
        ptRes0 = ptRes;
     }
     if ( (dR - dR0) < 0.  && dR >= 0.5 ) {
        betterMatch = true ;
        dR0    = dR ;
        ptRes0 = ptRes;
     }
     return betterMatch ;
}


std::vector<iTt> TtMCMatching::TtObjects( Handle<std::vector<reco::GenParticle> > genParticles ) {

   std::vector<iTt> ttobjs;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for W from top quark
       if ( abs((*it).pdgId()) == 24) {
  
          bool WfromT = false;
          for (size_t q=0; q< (*it).numberOfMothers(); q++) {
              const reco::Candidate *mom = (*it).mother(q) ;
              if ( abs(mom->pdgId()) != 6 ) continue;
              WfromT = true ;
          }
          if ( !WfromT ) continue;  

          // looking for the electron from W 
          for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
              const reco::Candidate *dau = (*it).daughter(q) ;
              if( abs(dau->pdgId()) > 18 ) continue;

              GlobalVector gv( dau->p4().Px(), dau->p4().Py(), dau->p4().Pz() );
	      iTt wdau;
	      wdau.pdgId = dau->pdgId() ;  
	      wdau.momId = (*it).pdgId();
	      wdau.p4    = dau->p4();
	      wdau.pt    = dau->pt();
	      wdau.gv    = gv;
              ttobjs.push_back( wdau ) ;
          }
       } 
       if ( abs((*it).pdgId()) == 5) {
          for (size_t q=0; q< (*it).numberOfMothers(); q++) {
              const reco::Candidate *mom = (*it).mother(q) ;
              if ( abs(mom->pdgId()) != 6 ) continue;
              GlobalVector gv( (it->p4()).Px(), (it->p4()).Py(), (it->p4()).Pz() );
	      iTt bfromT;
	      bfromT.pdgId = it->pdgId() ;  
	      bfromT.momId = mom->pdgId();
	      bfromT.p4    = it->p4();
	      bfromT.pt    = it->pt();
	      bfromT.gv    = gv;
              ttobjs.push_back( bfromT ) ;
          }
       } 
   }
   return ttobjs; 
}

// looking for objects from T decay => b, lep or quarks from W 
std::vector<reco::Particle> TtMCMatching::ttPartons( Handle<std::vector<reco::GenParticle> > genParticles, int targetId  ){

   // Accumulate the leptonic dauaghters from W
   std::vector<reco::Particle> targets;
   bool fromT =  ( abs(targetId) == 5 || abs(targetId) == 24 ) ? true:false ;
   bool fromW =  ( abs(targetId)  < 5 || (abs(targetId) > 10 && abs(targetId) < 19) )  ? true:false ;

   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for objects from W 
       if ( abs((*it).pdgId()) == 24 && fromW ) {
  
         /// make sure the W from Top
         bool WfromT = false;
         for (size_t q=0; q< (*it).numberOfMothers(); q++) {
             const reco::Candidate *mom = (*it).mother(q) ;
             // allow w+jets event pass
             if ( abs(mom->pdgId()) != 6 ) continue;
             //if ( abs(mom->pdgId()) > 6 ) continue;
             WfromT = true ;
         }
         if ( !WfromT ) continue;  

         // looking for objects from W 
         for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
             const reco::Candidate *dau = (*it).daughter(q) ;
	     if( abs(dau->pdgId()) == abs(targetId) ) targets.push_back(*dau);
         }
       }
       // looking for b or W from top quark
       if ( fromT ) {
  
          if ( abs((*it).pdgId()) != 6) continue;
          // looking for the jet daughter of b quark 
          for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
              const reco::Candidate *dau = (*it).daughter(q) ;
              if (dau->status() != 3 ) continue;
	      if( abs(dau->pdgId()) == abs(targetId) ) targets.push_back(*dau);
          }

       }    
   }
   return targets ;
}

std::vector<const reco::Candidate*> TtMCMatching::matchMuonfromB( Handle<std::vector<reco::GenParticle> > genParticles,
                std::vector<const pat::Jet*> jets, std::vector<const reco::Candidate*> theMuons ,HTOP7* histo7, bool fillhisto ){

   // Accumulate the candidates of leptonic dauaghters from B
   std::vector<reco::Particle> theB = ttDecay( genParticles, 5 ) ;

   std::vector<reco::Particle> lepB ;
   std::vector<reco::Particle> bmuon ;
   for( size_t i=0; i< theB.size(); i++ ) {
      std::vector<reco::Particle> bmu = genMuonFromB( genParticles, theB[i] );  
      if( bmu.size() > 0 ) {
        bmuon.push_back( bmu[0] );
        lepB.push_back( theB[i] );
      }
   }
  
   // find the matched muon passed isolation cut
   std::vector<const reco::Candidate*> matchedMuon ;
   for (size_t i=0; i< bmuon.size(); i++) {

      int usedMuon = -1;
      int usedB    = -1;

      // matching selected muons with parton
      double dR1 = 9.;
      double ptRes1 = 1.;
      for (size_t j=0; j < theMuons.size(); j++ ) {
          bool matched = matchingGeneral( bmuon[i].p4(), theMuons[j]->p4(), dR1, ptRes1 );
          if ( matched ) usedMuon = static_cast<int>(j);
      }
      // matched b parton with jets
      double dR2 = 9.;
      double ptRes2 = 1.;
      for (size_t j=0; j < jets.size(); j++ ) {
          bool matched = matchingGeneral( theB[i].p4(), jets[j]->p4(), dR2, ptRes2 );
          if ( matched ) usedB = static_cast<int>(j);
      }
      // 
      if ( usedMuon != -1 && usedB != -1 ) {
         double dR     = tools->getdR( theMuons[usedMuon]->p4(), jets[usedB]->p4() );
         double RelPt  = tools->getRelPt( theMuons[usedMuon]->p4(), jets[usedB]->p4() );
         double gdR    = tools->getdR( bmuon[i].p4(), lepB[i].p4() );
         double gRelPt = tools->getRelPt( bmuon[i].p4(), lepB[i].p4() );
         if (fillhisto) histo7->Fill7e( dR, RelPt, gdR, gRelPt );
         matchedMuon.push_back( theMuons[usedMuon] );
      }
   }

   return matchedMuon ;
}

// B quarks have hadronic daughters => no directory link btw muon & b  
std::vector<reco::Particle> TtMCMatching::genMuonFromB( Handle<std::vector<reco::GenParticle> > genParticles, reco::Particle theB){

   std::vector<reco::Particle> muonfromB ;
   std::vector<reco::Particle> bestMatched ;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for muon in final status
       if ( abs(it->pdgId()) != 13 || it->status() != 1 ) continue;

       /// make sure "it" NOT from Top/W
       bool fromW = false;
       for (size_t q=0; q< (*it).numberOfMothers(); q++) {
           const reco::Candidate *mom = (*it).mother(q) ;
           if ( abs(mom->pdgId()) == 6 || abs(mom->pdgId()) == 24 || abs(mom->pdgId()) == 13 ) fromW = true ;
           if ( abs(mom->pdgId()) == 15 ) fromW = true ;
       }
       if ( fromW ) continue;

       /// match B and muon within a certain dR cone ;
       double dR0 = tools->getdR( it->p4(), theB.p4() );
       if ( dR0 < 0.3 ) { 
          bestMatched.push_back( *it ) ;
       }
   }
   if ( bestMatched.size() > 0 ) {
      sort(bestMatched.begin(), bestMatched.end(), PtDecreasing );
      muonfromB.push_back( bestMatched[0] );
   }
   return muonfromB ;

}

// objects from Top => b,leptons, W(not available all the time, some events has no W record )
std::vector<reco::Particle> TtMCMatching::ttDecay( Handle<std::vector<reco::GenParticle> > genParticles, int targetID ){

   // Accumulate the leptonic dauaghters from W
   std::vector<reco::Particle> targets;

   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){

       // looking for objects from W  or T
       if ( abs(it->pdgId()) != targetID || it->status() != 3 ) continue;

       bool fromT = false;
       bool fromW = false;
       for (size_t q=0; q< (*it).numberOfMothers(); q++) {
           const reco::Candidate *mom = (*it).mother(q) ;
           if ( abs(mom->pdgId()) ==  6 ) fromT = true ;
           if ( abs(mom->pdgId()) == 24 ) fromW = true ;
       }
       if ( fromT || fromW ) targets.push_back( *it );
   }
   return targets ;
}

void TtMCMatching::CheckGenParticle(  Handle<std::vector<reco::GenParticle> > genParticles ) {
   

   cout<<" === new EVENT ========= "<<endl;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       if ( abs(it->pdgId()) ==6 ) {
          cout<<" top: "<< it->pdgId()<<" status: "<<it->status() <<endl;
          for (size_t q=0; q< (*it).numberOfMothers(); q++) {
              const reco::Candidate *mom = (*it).mother(q) ;
              cout<<"  => its mom:"<< mom->pdgId()<<" @"<<mom->status() <<endl;
          }
          for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
              const reco::Candidate *dau = (*it).daughter(q) ;
              cout<<"  => its dau:"<< dau->pdgId()<<" @"<<dau->status() <<endl;
          }
       }
       if ( abs(it->pdgId()) ==24 ) {
          cout<<" W  : "<< it->pdgId()<<" status: "<<it->status() <<endl;
          for (size_t q=0; q< (*it).numberOfMothers(); q++) {
              const reco::Candidate *mom = (*it).mother(q) ;
              cout<<"  => its mom:"<< mom->pdgId()<<" @"<<mom->status() <<endl;
          }
          for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
              const reco::Candidate *dau = (*it).daughter(q) ;
              cout<<"  => its dau:"<< dau->pdgId()<<" @"<<dau->status() <<endl;
          }
       }
       if ( abs(it->pdgId()) ==5 ) {
          cout<<" B  : "<< it->pdgId()<<" status: "<<it->status() <<endl;
          for (size_t q=0; q< (*it).numberOfMothers(); q++) {
              const reco::Candidate *mom = (*it).mother(q) ;
              cout<<"  => its mom:"<< mom->pdgId()<<" @"<<mom->status() <<endl;
          }
          for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
              const reco::Candidate *dau = (*it).daughter(q) ;
              cout<<"  => its dau:"<< dau->pdgId()<<" @"<<dau->status() <<endl;
          }
       }
       if ( abs(it->pdgId()) ==13 ) {
          cout<<" mu : "<< it->pdgId()<<" status: "<<it->status() <<endl;
          for (size_t q=0; q< (*it).numberOfMothers(); q++) {
              const reco::Candidate *mom = (*it).mother(q) ;
              cout<<"  => its mom:"<< mom->pdgId()<<" @"<<mom->status() <<endl;
          }
       }
   }
   cout<<" "<<endl;

}
