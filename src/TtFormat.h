#ifndef TtFormat_H
#define TtFormat_H
// -*- C++ -*-
//
// Package:    TtFormat
// Class:      TtFormat
// 
/**\class TtFormat TtFormat.cc PhysicsTools/TtAnalysis/src/TtFormat.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtFormat.h,v 1.2 2009/01/15 14:58:02 sckao Exp $
//
//


// system include files
#include <memory>

// user include files
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/EventHypothesis.h"
#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtAnalysisHisto.h"

#include "TFile.h"
#include "TVector3.h"
#include <vector>
#include <map>
#include <string>
#include <utility>

//
// class decleration
//
typedef math::XYZTLorentzVector LorentzVector;
typedef std::pair<int, LorentzVector> iParton;
// hfPos[0]:eta, hfPos[1]:phi, hfPos[2]:pt
typedef std::vector<double> hfPos ;
typedef std::vector<int> Idx;
typedef std::pair<int, double> IDPair;

enum ttChannel{ hadronic, semiMuon, dilep, semiElectron, semiTau, other };

//
struct iReco{
    LorentzVector p4;
    std::pair<int,int> from;
    std::pair<const reco::Candidate*, const reco::Candidate*> ptr;
    //std::pair<LorentzVector, LorentzVector> q4;
    std::vector<iParton> q4v ;
    double dm; // fill pt for letponic W
    double mt; // only filled for leptonic W
    double pt; 
};

// information of matched jets 
struct jmatch {
       int MomIdx ;
       double res_P;
       LorentzVector sumP4 ;
       std::vector<pat::Jet> assJets ;
       pat::Jet leadingJet ;
       const pat::Jet* trueJet ;
       reco::Particle mom ;
       bool hasMatched;
};

struct iTt {
    int pdgId;
    int momId;
    LorentzVector p4;
    GlobalVector  gv;
    double pt;
};

struct tHisto {

    HTOP1 *hJet;
    HTOP2 *hMET;
    HTOP3 *hMuon;
    HTOP4 *hEle;
    HTOP5 *hGam;
    HTOP6 *hMObj;
    HTOP7 *hBJet;
    HTOP8 *hWJet;
    HTOP9 *hTop;
    HTOP10 *hWs[4];   // 0: Lep-MC , 1: Lep-Reco , 2: Had-MC , 3: Had-Reco
    HTOP11 *hTops[4]; // 0: Lep-MC , 1: Lep-Reco , 2: Had-MC , 3: Had-Reco

};

struct TtResult {

    std::vector<iReco>                   Tt ;
    std::vector<const reco::Candidate*>  Ls ;
    std::vector<const pat::Jet*>         Js ;
    std::vector<const pat::Jet*>         bJ ;
    std::vector<const pat::Jet*>         WJ ;
    bool isData;
};

#endif
