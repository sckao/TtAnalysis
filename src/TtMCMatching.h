#ifndef TtMCMatching_H
#define TtMCMatching_H
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
// $Id: TtMCMatching.h,v 1.13 2009/10/15 06:47:16 sckao Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//include "FWCore/ParameterSet/interface/InputTag.h"
#include <FWCore/Utilities/interface/InputTag.h>


#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtAnalysisNtuple.h"
#include "TtAnalysisHisto.h"
#include "TtFormat.h"
#include "TtTools.h"

#include "TFile.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//

class TtMCMatching {
   public:
    /// Constructor
    explicit TtMCMatching();
    /// Destructor
    ~TtMCMatching();

    /// Perform the real analysis
    std::vector<const reco::Candidate*> GenTtCollection(edm::Handle<std::vector<reco::GenParticle> > genParticles ) ;

    std::vector<jmatch> matchJets(std::vector<const reco::Candidate*> genCollects, 
                    std::vector<ttCandidate>& selectedJets, HTOP7* histo7 = NULL, HTOP8* histo8 = NULL);

    std::vector<jmatch> matchWJets(edm::Handle<std::vector<reco::GenParticle> > genParticles,
                                   std::vector<const pat::Jet*> selectedWJets, HTOP8* histo8, bool fillhisto );

    std::vector<jmatch> matchbJets(edm::Handle<std::vector<reco::GenParticle> > genParticles,
                                   std::vector<const pat::Jet*> selectedbJets, HTOP7* histo7, bool fillhisto );

    std::vector<ttCandidate> matchMuon( std::vector<const reco::Candidate*> genCollects,
                                        std::vector<ttCandidate>& isoMuons, HTOP3* histo3 = NULL );

    std::vector<ttCandidate> matchElectron(std::vector<const reco::Candidate*> genCollects,
                                        std::vector<ttCandidate>& isoEle, HTOP4* histo4 = NULL );

    int matchLeptonicW(edm::Handle<std::vector<reco::GenParticle> > genParticles, std::vector<iReco> wSolution );

    bool matchingGeneral( LorentzVector p4_1, iTt ttObject, double& dR0, double& ptRes0 );

    bool matchingGeneral( LorentzVector aP4, LorentzVector bP4, double& dR0, double& ptRes0 );

    //std::vector<reco::Candidate> ttPartons( edm::Handle<std::vector<reco::GenParticle> > genParticles, int targetId ) ;
    std::vector<const reco::Candidate*> ttDecay( edm::Handle<std::vector<reco::GenParticle> > genParticles, int targetId ) ;
    std::vector<const reco::Candidate*> genMuonFromB( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::Candidate* genB ) ;

    void CheckGenParticle(  edm::Handle<std::vector<reco::GenParticle> > genParticles );

    double WBRCorrection_Ttbar( edm::Handle<std::vector<reco::GenParticle> > genParticles );

    //template<typename allT>
    //std::vector<const allT* > matchMuon( std::vector<const reco::Candidate*>  genCollects,
    //                                     std::vector<const allT*> isoMuons ,HTOP3* histo3 =NULL) ;
   private:

    TtTools*       tools;

};
/*
template<typename allT>
std::vector<const allT*> TtMCMatching::matchMuon( std::vector<const reco::Candidate*> genCollects,
                          std::vector<const allT*> isoMuons ,HTOP3* histo3 ){

   // Accumulate the leptonic dauaghters from W
   //std::vector<const reco::Candidate*> mcMuon = ttDecay(genParticles, 13) ;
   std::vector<const reco::Candidate*> mcMuon; 
   for (size_t i=0; i< genCollects.size(); i++) {
       if ( abs( genCollects[i]->pdgId() ) ==  13 ) mcMuon.push_back( genCollects[i] ); 
   }

   // find the matched muon 
   std::vector<int> matchList;
   for (size_t i=0; i< mcMuon.size(); i++) {

      // matching selected isoMuons with parton
      double dR1    = 9.;
      double ptRes1 = 9.;
      int usedMuon = -1;
      for (size_t j=0; j < isoMuons.size(); j++ ) {
          bool matched = matchingGeneral( mcMuon[i]->p4(), isoMuons[j]->p4(), dR1, ptRes1 );
          if ( matched && isoMuons[j]->charge() == mcMuon[i]->charge() ) usedMuon = static_cast<int>(j);
      }
      if (usedMuon != -1) {
         matchList.push_back( usedMuon );
         double ptRes2 = (isoMuons[usedMuon]->pt() /mcMuon[i]->pt() ) - 1. ;
         if ( histo3 != NULL ) histo3->Fill3b(          mcMuon[i]->eta(),          mcMuon[i]->phi(),
                                               isoMuons[usedMuon]->eta(), isoMuons[usedMuon]->phi(),
                                                mcMuon[i]->pt(), ptRes2 );
      }

   }
   // using one more loops to check whether double counting 
   std::vector<const allT*> matchedMuon ;

   std::vector<bool> used(isoMuons.size(), false  );

   for (std::vector<int>::const_iterator it= matchList.begin(); it != matchList.end(); it++) {
       if ( used[*it] ) continue;
       matchedMuon.push_back( isoMuons[*it] );
       used[*it] = true;
   }

   return matchedMuon ;
}
*/

#endif
