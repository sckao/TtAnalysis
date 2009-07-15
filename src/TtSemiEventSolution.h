#ifndef TtSemiEventSolution_H
#define TtSemiEventSolution_H
// -*- C++ -*-
//
// Package:    TtSemiEventSolution
// Class:      TtSemiEventSolution
// 
/**\class TtSemiEventSolution TtSemiEventSolution.cc PhysicsTools/TtAnalysis/src/TtSemiEventSolution.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtSemiEventSolution.h,v 1.11 2009/03/07 14:22:25 sckao Exp $
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
#include "FWCore/ParameterSet/interface/InputTag.h"


#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"

#include "DataFormats/MuonReco/interface/Muon.h" 
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/JetReco/interface/GenJet.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//#include "DataFormates/TrackReco/interface/Track.h"
//#include "DataFormates/TrackReco/interface/TrackFwd.h"

#include "TtAnalysisHisto.h"
#include "TtEvtSelector.h"
#include "TtMCMatching.h"
#include "TtMuon.h"
#include "TtElectron.h"
#include "TtMET.h"
#include "TtJet.h"
#include "TtEfficiency.h"
#include "TtFormat.h"
#include "TtTools.h"

#include "TFile.h"
#include "TVector3.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//
class TtEvtSelector;
class TtMCMatching;
class TtMuon;
class TtElectron;
class TtMET;
class TtJet;
class TtEfficiency;

typedef std::vector<int> Idx;

class TtSemiEventSolution {
   public:

    /// Constructor
    explicit TtSemiEventSolution( const edm::ParameterSet& );
    /// Destructor
    ~TtSemiEventSolution();

    /// Perform the real analysis
    void BuildSemiTt(const edm::Event & iEvent,   int topo, tHisto histos );
    void MCBuildSemiTt(const edm::Event & iEvent, int topo, tHisto histos );

    bool recoW( std::vector<const reco::Candidate*> wjets, std::vector<iReco>& wCandidate, HTOP10* histo10 = NULL, std::vector<bool>* btags = NULL );

    bool recoW( std::vector<const reco::Candidate*> lepton, LorentzVector metP4,
                std::vector<iReco>& wCandidate, HTOP10* hitso10 = NULL );
    bool recoW( std::vector<const reco::Candidate*> lepton, LorentzVector metP4,
                std::vector<iReco>& wCandidates, bool FoundWSolution, HTOP10* hitso10 = NULL );

    bool recoTop( std::vector<iReco> wCandidate, std::vector<const reco::Candidate*> bjets, std::vector<iReco>& tCandidate, std::vector<bool>* btags = NULL , HTOP11* histo11 = NULL);
    bool recoTop( std::vector<iReco> wCandidate, std::vector<const reco::Candidate*> bjets, std::vector<iReco>& tCandidate, bool btagging, HTOP11* histo11 = NULL );

    // with 2b tagging
    std::vector<iReco> recoSemiLeptonicTtEvent(int topo, std::vector<const reco::Candidate*> theWJets,
                  std::vector<const reco::Candidate*> thebJets, std::vector<const reco::Candidate*> theLep,
                  LorentzVector metP4, tHisto histos );
 
    // no b-tagging or 1 b-tagging
    std::vector<iReco> recoSemiLeptonicTtEvent(int topo, std::vector<const reco::Candidate*> theJets,
                  std::vector<const reco::Candidate*> theLep, LorentzVector metP4, tHisto histos, std::vector<bool>* btags = NULL   );
 
    void dmSortRecoObjects( std::vector<iReco>& objCand );

    void accuracySemiTt( std::vector<iReco> ttMC, std::vector<iReco> ttReco, HTOP9* histo9 ) ;

    void MCTruthCheck( std::vector<iReco> mcTt, std::vector<iReco> rcTt, int k,  std::vector<const reco::Candidate*> mcWJets,
                       std::vector<const reco::Candidate*> rcWJets, std::vector<const reco::Candidate*> mcBJets,
                       std::vector<const reco::Candidate*> rcBJets, HTOP6* histo6 ); 

    void McRecoCompare( int topo, int r, bool matchedpass, tHisto histos );

    void KeepBuildInfo( bool isData );

    void HadronicTopCombinatoric( std::vector<iReco>& tCandidates, std::vector<Idx>& tList, iReco t2, iReco w2, bool btagging );
    
    void Algo_dmMin( std::vector<iReco> lepTops, std::vector<iReco> hadTops, std::vector<iReco> hadWs, std::vector<Idx>& twIdx );
    void Algo_PtMin( std::vector<iReco> lepTops, std::vector<iReco> hadTops, std::vector<iReco> hadWs, std::vector<Idx>& twIdx );
    void Algo_Beta( std::vector<iReco> lepTops, std::vector<iReco> hadTops, std::vector<iReco> lepWs, std::vector<iReco> hadWs, std::vector<Idx>& twIdx );
    void Algo_Zero( std::vector<iReco> lepTops, std::vector<iReco> hadTops, std::vector<iReco> hadWs, std::vector<Idx>& twIdx, std::vector<bool>* bTags = NULL );
  
    void ResultRecord( int it, std::vector<Idx> twIdx, std::vector<iReco> lepTops, std::vector<bool>& usedLepT, 
                                                       std::vector<iReco> hadTops, std::vector<bool>& usedHadT,
                                          std::vector<iReco> lepWs, std::vector<iReco> hadWs , int type, tHisto histos );
   
   private:
      // ----------member data ---------------------------

    TtEvtSelector* evtSelected;
    TtMCMatching*  MCMatching;
    TtMuon*        ttMuon;
    TtElectron*    ttEle;
    TtMET*         ttMET;
    TtJet*         ttJet;
    TtEfficiency*  ttEff;
    TtTools*       tools;

    std::vector<iReco> semiTt;
    std::vector<const reco::Candidate*> isoLep;
    std::vector<const reco::Candidate*> selectedJets;
    std::vector<const reco::Candidate*> selectedbJets;
    std::vector<const reco::Candidate*> selectedWJets;

    std::vector<iReco> semiMCTt; 
    std::vector<const reco::Candidate*> mcLep;
    std::vector<const reco::Candidate*> mcWJets ;
    std::vector<const reco::Candidate*> mcbJets ;
   
    std::vector<TtResult> AllTt ;

    int pass;
    bool pure4Jet;

    // Switch for debug output
    bool debug;
    bool btag;

    std::string recoMuon;
    std::string algo;
    edm::InputTag muonSrc;
    edm::InputTag electronSrc;
    edm::InputTag metSrc;
    edm::InputTag tcmetSrc;
    edm::InputTag jetSrc;
    edm::InputTag jptSrc;
    edm::InputTag genJetSrc;
    edm::InputTag jetObj;
    edm::InputTag genSrc;
    edm::InputTag caloSrc;

};

#endif
