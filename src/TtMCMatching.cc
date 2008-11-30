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

}


TtMCMatching::~TtMCMatching()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtMCMatching Analysis] Destructor called" << endl;
}

//
// member functions
//
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
                               Handle<std::vector<pat::Jet> > jets, std::vector<const pat::Jet*> selectedWJets,
                               HTOP8* histo8, bool fillhisto ) {


   // Accumulate the hadronic dauaghters from W
   
   std::vector<reco::Particle> jetMom ;
   for (int i=1; i<5; i++) {
       std::vector<reco::Particle> tmpMom = ttPartons(genParticles, i) ;
       for (std::vector<reco::Particle>::iterator it= tmpMom.begin(); it!= tmpMom.end(); it++){
           jetMom.push_back( *it );
           if ( fillhisto ) histo8->Fill8g( it->eta(), it->pt() );
       }
   }
 
   // Loop all pat jets
   // get 4 collection w.r.t. possible 4 jet mom(quarks from W)
   // each collection contains candidate jets
   std::vector<int> q1; // q+ from w+
   std::vector<int> q2; // q- from w+
   std::vector<int> q3; // q+ from w-
   std::vector<int> q4; // q- from w-
   for (size_t i=0; i < jets->size(); i++){ 
       // exclude the jet with no charged track
       //edm::RefVector<reco::TrackCollection>  assTk = (*jets)[i].associatedTracks() ;
       //if ( assTk.size()== 0  ) continue;
       
       double dR0 = 99.;
       double ptRes0 =99.;
       int wj = -1;
       for (size_t j=0; j < jetMom.size(); j++ ) {
           bool matched = matchingGeneral( jetMom[j].p4() , (*jets)[i].p4() , dR0 , ptRes0 );
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
          bool matched = matchingGeneral( jetMom[j].p4() , (*jets)[*it].p4() , dR0 , ptRes0 ); 
          if (matched) wj = *it ;
       }


       if (wj != -1) {
          jmatch theWjet;
	  theWjet.trueJet= &(*jets)[wj] ;
	  theWjet.mom    = jetMom[j] ;
	  theWjet.MomIdx = jetMom[j].pdgId() ;
	  theWjet.res_P  = ((*jets)[wj].pt() /jetMom[j].pt()) - 1. ;
	  theWjet.sumP4  = (*jets)[wj].p4() ;
	  theWjet.hasMatched = true ;
          matchedJets.push_back( theWjet );
          if ( fillhisto ) histo8->Fill8h( jetMom[j].eta(), jetMom[j].pt(), theWjet.res_P  );
       }

       // matched what have been selected => to look up the efficiency
       double dR1 = 99.;
       double ptRes1 =99.;
       int sj = -1;
       std::vector<bool> used( selectedWJets.size(), false );
       for(size_t k =0; k < selectedWJets.size(); k++ ){
          if ( used[k] ) continue;
          bool matched = matchingGeneral( jetMom[j].p4() , selectedWJets[k]->p4() , dR1 , ptRes1 ); 
          if (matched) sj = static_cast<int>(k) ;
       }
       if (sj != -1) {
          used[sj] = true;
          double ptRes2 = (selectedWJets[sj]->pt()/jetMom[j].pt()) - 1. ;
          if ( fillhisto ) histo8->Fill8i( jetMom[j].eta(), jetMom[j].pt(), ptRes2 );
       }
   }
   return matchedJets ;
}

std::vector<jmatch> TtMCMatching::matchbJets( Handle<std::vector<reco::GenParticle> > genParticles,
                                              Handle<std::vector<pat::Jet> > jets, std::vector<const pat::Jet*> selectedbJets
                                             ,HTOP7* histo7, bool fillhisto ) {
   // Accumulate the b quark from t
   std::vector<reco::Particle> jetMom = ttPartons(genParticles, 5) ;
   if (fillhisto) {
      for (std::vector<reco::Particle>::iterator it = jetMom.begin(); it!=jetMom.end(); it++ ) {
          histo7->Fill7g( it->eta(), it->pt() );
       }
   }
 
   // Loop all pat jets
   // get 2 collections w.r.t. possible 2 bjet mom( b and bbar )
   // each collection contains candidate jets to avoid double using
   std::vector<int> bp;
   std::vector<int> bn;
   for (size_t i=0; i < jets->size(); i++){ 
       // exclude the jet with no charged track
       //edm::RefVector<reco::TrackCollection>  assTk = (*jets)[i].associatedTracks() ;
       //if ( assTk.size()== 0  ) continue;
       
       double dR0 = 99.;
       double ptRes0 =99.;
       int bj = -1;
       for (size_t j=0; j < jetMom.size(); j++ ) {
           bool matched = matchingGeneral( jetMom[j].p4() , (*jets)[i].p4() , dR0 , ptRes0 );
           if ( matched )  bj = static_cast<int>(j) ;
       }

       if ( bj!= -1 && jetMom[bj].pdgId() ==  5 ) bp.push_back(  static_cast<int>(i) ) ;
       if ( bj!= -1 && jetMom[bj].pdgId() == -5 ) bn.push_back(  static_cast<int>(i) ) ;
      
   }

   // loop 2 collections to find the best matched jet 
   std::vector<jmatch> matchedJets ;
   std::vector<int> btemp;   
   for (size_t j=0; j < jetMom.size(); j++  ) {
       
       if (jetMom[j].pdgId() ==  5) btemp = bp ;
       if (jetMom[j].pdgId() == -5) btemp = bn ;

       double dR0 = 99.;
       double ptRes0 =99.;
       int bj = -1;
       for(std::vector<int>::iterator it=btemp.begin(); it != btemp.end(); it++  ){
          bool matched = matchingGeneral( jetMom[j].p4() , (*jets)[*it].p4() , dR0 , ptRes0 ); 
          if (matched) bj = *it ;
       }

       if (bj != -1) {
          jmatch thebjet;
	  thebjet.trueJet= &(*jets)[bj] ;
	  thebjet.mom    = jetMom[j] ;
	  thebjet.MomIdx = jetMom[j].pdgId() ;
	  thebjet.res_P  = ((*jets)[bj].pt() /jetMom[j].pt()) - 1. ;
	  thebjet.sumP4  = (*jets)[bj].p4() ;
	  thebjet.hasMatched = true ;
          matchedJets.push_back( thebjet );
          if ( fillhisto ) histo7->Fill7h( jetMom[j].eta(), jetMom[j].pt(), thebjet.res_P );
       }

       // matched what have been selected => to look up the efficiency
       double dR1 = 99.;
       double ptRes1 =99.;
       int sj = -1;
       std::vector<bool> used( selectedbJets.size(), false );
       for(size_t k =0; k < selectedbJets.size(); k++ ){
          if ( used[k] ) continue;
          bool matched = matchingGeneral( jetMom[j].p4() , selectedbJets[k]->p4() , dR1 , ptRes1 ); 
          if (matched) sj = static_cast<int>(k) ;
       }
       if (sj != -1) {
          used[sj] = true;
          double ptRes2 = (selectedbJets[sj]->pt()/jetMom[j].pt()) - 1. ;
          if ( fillhisto ) histo7->Fill7i( jetMom[j].eta(), jetMom[j].pt(), ptRes2 );
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
           Handle<std::vector<pat::Muon> > muons, std::vector<const reco::Candidate*> isoMuons ,HTOP3* histo3, bool fillhisto ){

   // Accumulate the leptonic dauaghters from W
   std::vector<reco::Particle> mcMuon = ttPartons(genParticles, 13) ;
   if (fillhisto) {
      for (std::vector<reco::Particle>::iterator it = mcMuon.begin(); it!=mcMuon.end(); it++ ) {
          histo3->Fill3f( it->eta(), it->pt() );
       }
   }

   // find the matched muon 
   std::vector<int> matchList;

   for (size_t i=0; i< mcMuon.size(); i++) {

      // matching pat::muons with parton
      double dR0 = 99.;
      double ptRes0 =99.;
      int theMuon = -1;
      for (size_t j=0; j < muons->size(); j++ ) {
          bool matched = matchingGeneral( (*muons)[j].p4() , mcMuon[i].p4(), dR0, ptRes0 );
          if ( matched ) theMuon = static_cast<int>(j);
      }
      if (theMuon != -1) {
         matchList.push_back( theMuon ); 
         double ptRes2 = ((*muons)[theMuon].pt()/mcMuon[i].pt()) - 1. ;
         if (fillhisto) histo3->Fill3g( mcMuon[i].eta(), mcMuon[i].pt(), ptRes2 );
      }
    
      // matching selected isoMuons with parton
      double dR1 = 99.;
      double ptRes1 =99.;
      int usedMuon = -1;
      for (size_t j=0; j < isoMuons.size(); j++ ) {
          bool matched = matchingGeneral( isoMuons[j]->p4() , mcMuon[i].p4(), dR1, ptRes1 );
          if ( matched ) usedMuon = static_cast<int>(j);
      }
      if (usedMuon != -1) {
         double ptRes2 = (isoMuons[usedMuon]->pt()/mcMuon[i].pt()) - 1. ;
         if (fillhisto) histo3->Fill3h( mcMuon[i].eta(), mcMuon[i].pt(), ptRes2 );
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
      Handle<std::vector<pat::Electron> > electrons, std::vector<const reco::Candidate*> isoEle, HTOP4* histo4, bool fillhisto ){

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> mcElectron = ttPartons(genParticles, 11) ;
   if (fillhisto) {
      for (std::vector<reco::Particle>::iterator it = mcElectron.begin(); it!=mcElectron.end(); it++ ) {
          histo4->Fill4f( it->eta() , it->pt() );
       }
   }

   // find the matched electron 
   std::vector<int> matchList;
   for (size_t i=0; i< mcElectron.size(); i++) {

      // matching pat::electrons with parton
      double dR0 = 99.;
      double ptRes0 =99.;
      int theElectron = -1;
      for (size_t j=0; j < electrons->size(); j++ ) {
          bool matched = matchingGeneral( (*electrons)[j].p4() , mcElectron[i].p4(), dR0, ptRes0 );
          if ( matched ) theElectron = static_cast<int>(j);
      }
      if (theElectron != -1) {
         matchList.push_back( theElectron ); 
         double ptRes2 = ((*electrons)[theElectron].pt()/mcElectron[i].pt()) - 1. ;
         if (fillhisto) histo4->Fill4g( mcElectron[i].eta(), mcElectron[i].pt(), ptRes2 );
      }
    
      // matching selected isoMuons with parton
      double dR1 = 99.;
      double ptRes1 =99.;
      int usedEle = -1;
      for (size_t j=0; j < isoEle.size(); j++ ) {
          bool matched = matchingGeneral( isoEle[j]->p4() , mcElectron[i].p4(), dR1, ptRes1 );
          if ( matched ) usedEle = static_cast<int>(j);
      }
      if (usedEle != -1) {
         double ptRes2 = (isoEle[usedEle]->pt()/mcElectron[i].pt()) - 1. ;
         if (fillhisto) histo4->Fill4h( mcElectron[i].eta(), mcElectron[i].pt(), ptRes2 );
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

int TtMCMatching::matchLeptonicW( Handle<std::vector<reco::GenParticle> > genParticles, std::vector<iReco> wSolutions ){

   int wl = -1;
   LorentzVector wP4(0.,0.,0.,0.);
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for object from W 
       if ( abs((*it).pdgId()) != 24 ) continue;
       for (size_t q=0; q< (*it).numberOfDaughters(); q++) {
           const reco::Candidate *dau = (*it).daughter(q) ;
           if ( abs(dau->pdgId()) ==  13 ) wP4 = it->p4() ;
       }
   }

   LorentzVector zero(0.,0.,0.,0.);
   if (wP4 == zero) return wl;
   
   double dR0 = 99.;
   double ptRes0 = 99.;
   for (size_t j=0; j < wSolutions.size(); j++ ) {
       bool matched = matchingGeneral( wP4 , wSolutions[j].p4 , dR0 , ptRes0 );
       if ( matched )  wl = static_cast<int>(j) ;
   }
   return wl;

}

int TtMCMatching::matchLeptonicW( Handle<std::vector<reco::GenParticle> > genParticles, std::vector<iReco> wSolutions, HTOP9* histo9 ){

   int wl = -1;
   LorentzVector wP4(0.,0.,0.,0.);
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for object from W
       if ( abs((*it).pdgId()) != 24 ) continue;
       for (size_t q=0; q< (*it).numberOfDaughters(); q++) {
           const reco::Candidate *dau = (*it).daughter(q) ;
           if ( abs(dau->pdgId()) ==  13 ) wP4 = it->p4() ;
       }
   }

   LorentzVector zero(0.,0.,0.,0.);
   if (wP4 == zero) return wl;

   double dR0 = 99.;
   double ptRes0 = 99.;
   for (size_t j=0; j < wSolutions.size(); j++ ) {
       bool matched = matchingGeneral( wP4 , wSolutions[j].p4 , dR0 , ptRes0 );
       if ( matched )  wl = static_cast<int>(j) ;
   }

   if ( wl != -1 ) {
      LorentzVector lv = wSolutions[wl].p4 ;
      double mom = (lv.Px()*lv.Px()) + (lv.Py()*lv.Py()) + (lv.Pz()*lv.Pz()) ;
      double mass = sqrt( lv.E()*lv.E() - mom );
      double gmom = (wP4.Px()*wP4.Px()) + (wP4.Py()*wP4.Py()) + (wP4.Pz()*wP4.Pz()) ;
      double gmass = sqrt( wP4.E()*wP4.E() - gmom );
      double dM = (mass/gmass) - 1. ;
      histo9->Fill9m( dR0, ptRes0, dM);
   }
   return wl;
}


bool TtMCMatching::matchingGeneral( LorentzVector theP4 , iTt ttObj, double& dR0, double& ptRes0 ){

     bool betterMatch = false;

     double pt1 = sqrt( (theP4.Px()*theP4.Px())  + (theP4.Py()*theP4.Py()) );
     double pa1 = sqrt( (pt1*pt1) + (theP4.Pz()*theP4.Pz()) );
     GlobalVector gp1( theP4.Px()/pa1, theP4.Py()/pa1, theP4.Pz()/pa1 );
   
     double dh = gp1.eta() - ttObj.gv.eta() ;
     double df = gp1.phi() - ttObj.gv.phi() ; 
     double dR = sqrt( (dh*dh) + (df*df) );
     // 1/pt resolution
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
     double pa1 = sqrt( (pt1*pt1) + (aP4.Pz()*aP4.Pz()) );
     GlobalVector gp1( aP4.Px()/pa1, aP4.Py()/pa1, aP4.Pz()/pa1 );
   
     double pt2 = sqrt( (bP4.Px()*bP4.Px())  + (bP4.Py()*bP4.Py()) );
     double pa2 = sqrt( (pt2*pt2) + (bP4.Pz()*bP4.Pz()) );
     GlobalVector gp2( bP4.Px()/pa2, bP4.Py()/pa2, bP4.Pz()/pa2 );
   
     double dh = gp1.eta() - gp2.eta() ;
     double df = gp1.phi() - gp2.phi() ; 
     double dR = sqrt( (dh*dh) + (df*df) );
     // 1/pt resolution
     double ptRes = fabs( (pt2/pt1 ) - 1.) ;

     if ( (dR - dR0) < 0.05 && (dR < 0.5) && ( ptRes < ptRes0 ) ) {
        betterMatch = true ;
        dR0    = dR ;
        ptRes0 = ptRes;
     }
     if ( (dR - dR0) < 0.  && dR >= 0.5 && dR < 0.7) {
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

std::vector<reco::Particle> TtMCMatching::ttPartons( Handle<std::vector<reco::GenParticle> > genParticles, int targetId  ){

   // Accumulate the leptonic dauaghters from W
   std::vector<reco::Particle> targets;
   bool fromT =  ( abs(targetId) == 5 || abs(targetId) == 24 ) ? true:false ;
   bool fromW =  ( abs(targetId)  < 5 || (abs(targetId) > 10 && abs(targetId) < 19) )  ? true:false ;

   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for object from W 
       if ( abs((*it).pdgId()) == 24 && fromW ) {
  
         /// make sure the W from proton
         bool WfromT = false;
         for (size_t q=0; q< (*it).numberOfMothers(); q++) {
             const reco::Candidate *mom = (*it).mother(q) ;
             // allow w+jets event pass
             //if ( abs(mom->pdgId()) != 6 ) continue;
             if ( abs(mom->pdgId()) > 6 ) continue;
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
	      if( abs(dau->pdgId()) == abs(targetId) ) targets.push_back(*dau);
          }

       }    
   }
   return targets ;
}

