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

std::vector<jmatch> TtMCMatching::matchWJets( Handle<std::vector<reco::GenParticle> > genParticles,
                               Handle<std::vector<pat::Jet> > jets ) {


   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> jetMom;
   jetMom.clear();
   bool WfromT = false;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for W from top quark
       if ( abs((*it).pdgId()) != 24) continue;

  
       for (size_t q=0; q< (*it).numberOfMothers(); q++) {
           const reco::Candidate *mom = (*it).mother(q) ;
           if ( abs(mom->pdgId()) != 6 ) continue;
           WfromT = true ;
       }
       if ( !WfromT ) continue;  

       // looking for the jet daughter of W 
       for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
           const reco::Candidate *dau = (*it).daughter(q) ;
	   if( abs(dau->pdgId()) > 6 ) continue;
	   jetMom.push_back( *dau );
       }
   }

   // Loop all pat jets
   // get 4 collection w.r.t. possible 4 jet mom(quarks from W)
   // each collection contains candidate jets
   int jid = -1;
   std::vector<int> q1;
   std::vector<int> q2;
   std::vector<int> q3;
   std::vector<int> q4;
   q1.clear();
   q2.clear();
   q3.clear();
   q4.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {
       jid++;

       // exclude the jet with no charged track
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;
       if ( assTk.size()== 0  ) continue;

       // Loop all q qbar from W to match with jets
       int jmflag=0;
       double dR1 = 999. ;
       int jm=0;
       for (std::vector<reco::Particle>::iterator it = jetMom.begin(); it != jetMom.end(); it++ ){

           jm++;

           double dh = (*j1).eta() - (*it).eta() ;
           double df = (*j1).phi() - (*it).phi() ;
           double dR = sqrt( dh*dh + df*df ) ;
           if (dR < dR1) {
              dR1 = dR;
              jmflag = jm ;
           }
       }
       if ( jmflag == 1 ) q1.push_back(jid) ;
       if ( jmflag == 2 ) q2.push_back(jid) ;
       if ( jmflag == 3 ) q3.push_back(jid) ;
       if ( jmflag == 4 ) q4.push_back(jid) ;
   }

   // best match
   // for each q collection, find the best matched patJet
   std::vector<jmatch> matchedJets ;
   std::vector<int> qq;
   int jm1 =0;
   for (std::vector<reco::Particle>::iterator it = jetMom.begin(); it != jetMom.end(); it++ ){

       qq.clear();
       jm1++;
       if ( jm1 == 1) qq = q1;
       if ( jm1 == 2) qq = q2;
       if ( jm1 == 3) qq = q3;
       if ( jm1 == 4) qq = q4;
       if ( qq.size() == 0 ) continue;

       std::vector<pat::Jet> jv = *jets ;

       std::vector<pat::Jet> bestJet ;
       LorentzVector sumP4 ;

       double jetERatio = 0.0;
       double biggestP = 0.0 ;
       double dR1 = 0.5;

       pat::Jet leader ;
       pat::Jet truth ;
       for (size_t i=0; i< qq.size(); i++) {

           pat::Jet jt = jv[ qq[i] ];
           double dh = jt.eta() - (*it).eta() ;
           double df = jt.phi() - (*it).phi() ;
           double dR = sqrt( dh*dh + df*df ) ;
           double ERatio = jt.p()/(*it).p() ;
           double theP = jt.p() ;

           if ( dR < 0.5 ) {
              if ( (jetERatio + ERatio) > 1.05 && jetERatio > 0.5 ) continue; 
              jetERatio += ERatio ;
              bestJet.push_back( jt ) ;
              sumP4 += jt.p4() ;
              if ( dR < dR1 && abs(jt.partonFlavour()) < 5 ){
                 dR1 = dR ;
                 truth = jt;
              }
              if (theP > biggestP ) {
                 biggestP = theP;
                 leader = jt ; 
              }
           }

       }

       if (bestJet.size() < 1) continue;

       jmatch wjets;
       wjets.assJets =  bestJet ;
       wjets.leadingJet = leader ;
       wjets.trueJet = truth ;
       wjets.mom = (*it) ;
       wjets.MomIdx = jm1;
       wjets.res_P = 1. - jetERatio ;
       wjets.sumP4 = sumP4 ;
      
       matchedJets.push_back(wjets);
   }
   return matchedJets ;
}

std::vector<jmatch> TtMCMatching::matchbJets( Handle<std::vector<reco::GenParticle> > genParticles,
                                              Handle<std::vector<pat::Jet> > jets ) {

   // Accumulate the b quark from t
   std::vector<reco::Particle> jetMom;
   jetMom.clear();
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for b from top quark
       if ( abs((*it).pdgId()) != 6) continue;
  
       // looking for the jet daughter of b quark 
       for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
           const reco::Candidate *dau = (*it).daughter(q) ;
	   if( abs(dau->pdgId()) != 5 ) continue;
	   jetMom.push_back( *dau );
       }
   }
 
   // Loop all pat jets
   // get 4 collection w.r.t. possible 4 jet mom(quarks from W)
   // each collection contains candidate jets
   int jid = -1;
   std::vector<int> q1;
   std::vector<int> q2;
   q1.clear();
   q2.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {
       jid++;
       // exclude the jet with no charged track
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;
       if ( assTk.size()== 0  ) continue;

       // Loop all q qbar from W to match with jets
       int jmflag=0;
       double dR1 = 999. ;
       int jm=0;
       for (std::vector<reco::Particle>::iterator it = jetMom.begin(); it != jetMom.end(); it++ ){

           jm++;
           double dh = (*j1).eta() - (*it).eta() ;
           double df = (*j1).phi() - (*it).phi() ;
           double dR = sqrt( dh*dh + df*df ) ;
           if (dR < dR1) {
              dR1 = dR;
              jmflag = jm ;
           }
       }
       if ( jmflag == 1 ) q1.push_back(jid) ;
       if ( jmflag == 2 ) q2.push_back(jid) ;
   }
   
   // best match
   // for each q collection, find the best matched patJet
   std::vector<jmatch> matchedJets ;
   std::vector<int> qq;
   int jm1 =0;
   for (std::vector<reco::Particle>::iterator it = jetMom.begin(); it != jetMom.end(); it++ ){

       qq.clear();
       jm1++;
       if ( jm1 == 1) qq = q1;
       if ( jm1 == 2) qq = q2;
       if ( qq.size() == 0 ) continue;

       std::vector<pat::Jet> jv = *jets ;
       std::vector<pat::Jet> bestJet ;

       LorentzVector sumP4 ;
       double jetERatio = 0.0;
       double biggestP = 0.0 ;
       double dR1 = 0.5;

       pat::Jet leader ;
       pat::Jet truth ;
       for (size_t i=0; i< qq.size(); i++) {

           pat::Jet jt = jv[ qq[i] ];
           double dh = jt.eta() - (*it).eta() ;
           double df = jt.phi() - (*it).phi() ;
           double dR = sqrt( dh*dh + df*df ) ;
           double ERatio = jt.p()/(*it).p() ;
           double theP = jt.p() ;

           if ( dR < 0.5 ) {
              if ( (jetERatio + ERatio) > 1.05 && jetERatio > 0.5 ) continue; 
              jetERatio += ERatio ;
              bestJet.push_back( jt ) ;
              sumP4 += jt.p4() ; 
              if ( dR < dR1 && abs(jt.partonFlavour()) == 5 ){
                 dR1 = dR ;
                 truth = jt;
              }
              if (theP > biggestP ) {
                 biggestP = theP;
                 leader = jt ; 
              }
           }
       }

       if (bestJet.size() < 1) continue;

       jmatch bjets;
       bjets.assJets =  bestJet ;
       bjets.leadingJet = leader ;
       bjets.trueJet = truth ;
       bjets.mom = (*it) ;
       bjets.MomIdx = jm1;
       bjets.res_P = 1. - jetERatio ;
       bjets.sumP4 = sumP4 ;
      
       matchedJets.push_back(bjets);
   }
   return matchedJets ;


}

//std::vector<pat::Muon> TtMCMatching::matchMuon( Handle<std::vector<reco::GenParticle> > genParticles,
std::vector<const reco::Candidate*> TtMCMatching::matchMuon( Handle<std::vector<reco::GenParticle> > genParticles,
                                               Handle<std::vector<pat::Muon> > muons ) {

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> mcMuon;
   mcMuon.clear();
   bool WfromT = false;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for W from top quark
       if ( abs((*it).pdgId()) != 24) continue;
  
       for (size_t q=0; q< (*it).numberOfMothers(); q++) {
           const reco::Candidate *mom = (*it).mother(q) ;
           if ( abs(mom->pdgId()) != 6 ) continue;
           WfromT = true ;
       }
       if ( !WfromT ) continue;  

       // looking for the muon from W 
       for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
           const reco::Candidate *dau = (*it).daughter(q) ;
	   if( abs(dau->pdgId()) != 13 ) continue;
           mcMuon.push_back(*dau);
       }
   }

   // find the matched muon 
   std::vector<int> matchList;
   for (size_t i=0; i< mcMuon.size(); i++) {
       
      double dRmin = 0.5;
      double dPmin = 9999.;
      int theMuon = -1;
      int idx = -1;
      for (std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end(); it++ ) {
          idx++;
          double dh = mcMuon[i].eta() - (*it).eta() ;
          double df = mcMuon[i].phi() - (*it).phi() ;
          double dR = sqrt((dh*dh) + (df*df)) ;
          double dP = fabs(it->pt() - mcMuon[i].pt());
          bool close = ( fabs(dR - dRmin) < 0.01 && dRmin != 0.5 ) ? true:false ; 
          if ( dR > dRmin && !close ) continue; 
          if ( fabs(dR - dRmin) < 0.01 ) {
             if ( dP > dPmin ) continue;
             dPmin = dP;
             dRmin = dR;
             theMuon = idx ;
          } else {  
             dPmin = dP;
             dRmin = dR;
             theMuon = idx ;
          }
      }
      if (theMuon != -1) matchList.push_back( theMuon ); 
   }
   // using one more loops to check whether double counting 
   int id = -1;
   //std::vector<pat::Muon> matchedMuon ;
   std::vector<const reco::Candidate*> matchedMuon ;
   for (std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end(); it++ ) {
       id++;
       int count =0 ;
       for (std::vector<int>::const_iterator it1= matchList.begin(); it1 != matchList.end(); it1++) { 
          if (id == (*it1) && count < 1 ) matchedMuon.push_back( &*it );
          count++; 
       }  
   }
   if ( matchedMuon.size() > 2 ) cout<<" WRONG MATCHING "<<endl;

   return matchedMuon ;
}

std::vector<const reco::Candidate*> TtMCMatching::matchElectron( Handle<std::vector<reco::GenParticle> > genParticles,
                                             Handle<std::vector<pat::Electron> > electrons ) {

   // Accumulate the hadronic dauaghters from W
   std::vector<reco::Particle> mcElectron;
   mcElectron.clear();
   bool WfromT = false;
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       // looking for W from top quark
       if ( abs((*it).pdgId()) != 24) continue;
  
       for (size_t q=0; q< (*it).numberOfMothers(); q++) {
           const reco::Candidate *mom = (*it).mother(q) ;
           if ( abs(mom->pdgId()) != 6 ) continue;
           WfromT = true ;
       }
       if ( !WfromT ) continue;  

       // looking for the electron from W 
       for (size_t q=0; q< (*it).numberOfDaughters(); ++q) {
           const reco::Candidate *dau = (*it).daughter(q) ;
	   if( abs(dau->pdgId()) != 11 ) continue;
           mcElectron.push_back(*dau);
       }
   }

   // find the matched electron 
   std::vector<int> matchList;
   for (size_t i=0; i< mcElectron.size(); i++) {
       
      double dRmin = 0.5;
      double dPmin = 9999.;
      int theElectron = -1;
      int idx = -1;
      for (std::vector<pat::Electron>::const_iterator it = electrons->begin(); it != electrons->end(); it++ ) {
          idx++;
          double dh = mcElectron[i].eta() - (*it).eta() ;
          double df = mcElectron[i].phi() - (*it).phi() ;
          double dR = sqrt((dh*dh) + (df*df)) ;
          double dP = fabs(it->pt() - mcElectron[i].pt());
          bool close = ( fabs(dR - dRmin) < 0.01 && dRmin != 0.5 ) ? true:false ; 
          if ( dR > dRmin && !close ) continue; 
          if ( fabs(dR - dRmin) < 0.01 ) {
             if ( dP > dPmin ) continue;
             dPmin = dP;
             dRmin = dR;
             theElectron = idx ;
          } else {  
             dPmin = dP;
             dRmin = dR;
             theElectron = idx ;
          }
      }
      if (theElectron != -1) matchList.push_back( theElectron ); 
   }
   // using one more loops to check whether double counting 
   int id = -1;
   //std::vector<pat::Electron> matchedElectron ;
   std::vector<const reco::Candidate*> matchedElectron ;
   for (std::vector<pat::Electron>::const_iterator it = electrons->begin(); it != electrons->end(); it++ ) {
       id++;
       int count =0 ;
       for (std::vector<int>::const_iterator it1= matchList.begin(); it1 != matchList.end(); it1++) { 
          if (id == (*it1) && count < 1 ) matchedElectron.push_back( &*it );
          count++; 
       }  
   }
   if ( matchedElectron.size() > 2 ) cout<<" WRONG MATCHING "<<endl;

   return matchedElectron ;
}

