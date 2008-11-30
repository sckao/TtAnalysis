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
// $Id$
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

//
struct iReco{
    LorentzVector p4;
    std::pair<int,int> from;
    std::pair<const reco::Candidate*, const reco::Candidate*> ptr;
    //std::pair<LorentzVector, LorentzVector> q4;
    std::vector<iParton> q4v ;
    double dm;
    double mt; // only filled for leptonic W
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

#endif
