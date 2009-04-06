#ifndef TtTools_H
#define TtTools_H
// -*- C++ -*-
//
// Package:    TtTools
// Class:      TtTools
// 
/**\class TtTools TtTools.cc PhysicsTools/TtAnalysis/src/TtTools.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Fri May 16 2008
// $Id: TtTools.h,v 1.3 2009/03/07 14:22:25 sckao Exp $
//
//


// system include files
#include <memory>

// For propagation
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Track Calo Mapping!
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"

//#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
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
#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/MuonReco/interface/Muon.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/Candidate/interface/Particle.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TtFormat.h"

#include "TFile.h"
#include "TVector3.h"
#include <vector>
#include <map>
#include <string>
#include <utility>


//
// class decleration
//

class TtTools {
   public:
    /// Constructor
    explicit TtTools();
    /// Destructor
    ~TtTools();

    /// Perform the real analysis
    FreeTrajectoryState getFTS(GlobalPoint GP, GlobalVector GV, int charge,
                               const AlgebraicSymMatrix66& cov, const MagneticField* field);

    double getEta(double vx, double vy, double vz );
    double getdPhi(  LorentzVector v1, LorentzVector v2 );
    double getdR(  LorentzVector v1, LorentzVector v2 );
    double getdRy(  LorentzVector v1, LorentzVector v2 );
    double getY( LorentzVector v1 );
    double getRelPt( LorentzVector a, LorentzVector b );
    double getBeta( LorentzVector a );
    double getInvMass( std::vector<LorentzVector> vlist );
    double getInvMass( LorentzVector lv );
    double getInvMass( LorentzVector lv1, LorentzVector lv2 );

    std::vector<const pat::Jet*> ReturnJetForm( std::vector<const reco::Candidate*> jCand, edm::Handle<std::vector<pat::Jet> > patJet );
    const pat::Jet* ReturnJetForm( const reco::Candidate* jCand, edm::Handle<std::vector<pat::Jet> > patJet, bool& goodmatching );

   private:


};

#endif
