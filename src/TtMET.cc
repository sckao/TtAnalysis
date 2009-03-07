// -*- C++ -*-
//
// Package:    TtMET
// Class:      TtMET
// 
/**\class TtMET TtMET.cc PhysicsTools/TtAnalysis/src/TtMET.cc

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
#include "TtMET.h"

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
TtMET::TtMET(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  /*
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  //recoJet           = iConfig.getUntrackedParameter<string> ("recoJets");
   */

  muonSrc         = iConfig.getParameter<edm::InputTag> ("muonSource");
  caloSrc         = iConfig.getParameter<edm::InputTag> ("caloSource");
  genSrc          = iConfig.getParameter<edm::InputTag> ("genParticles"); 

  ttMuon      = new TtMuon();
  tools           = new TtTools();
}


TtMET::~TtMET()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtMET Analysis] Destructor called" << endl;
   delete ttMuon;
   delete tools;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------

void TtMET::MetTreeFeeder(Handle<std::vector<pat::MET> > patMet, NJet* jtree, int eventId ) {

  for (std::vector<pat::MET>::const_iterator m1 = patMet->begin(); m1 != patMet->end(); m1++)
  {
     float emMET = (*m1).emEtInEB()  + (*m1).emEtInEE()  + (*m1).emEtInHF() ;
     float hdMET = (*m1).hadEtInHB() + (*m1).hadEtInHE() + (*m1).hadEtInHF() + (*m1).hadEtInHO() ;

     jtree->FillBpatNu( eventId, m1->eta(), m1->phi(), emMET, hdMET, m1->p(), m1->pt() );
  }
}

void TtMET::metAnalysis(Handle<std::vector<pat::MET> > patMet, const edm::Event& iEvent, HTOP2* histo2) {

 Handle<std::vector<pat::Muon> > muons;
 iEvent.getByLabel(muonSrc, muons);

 Handle<CaloTowerCollection>  caloTowers;
 iEvent.getByLabel(caloSrc, caloTowers);

 Handle<std::vector<reco::GenParticle> > genParticles;
 iEvent.getByLabel(genSrc, genParticles);


 for (std::vector<pat::MET>::const_iterator m1 = patMet->begin(); m1 != patMet->end(); m1++)
 {
     float emMET = (*m1).emEtInEB()  + (*m1).emEtInEE()  + (*m1).emEtInHF() ;
     float hdMET = (*m1).hadEtInHB() + (*m1).hadEtInHE() + (*m1).hadEtInHF() + (*m1).hadEtInHO() ;
     float emfCalo = 1. - (emMET/(emMET + hdMET)) ;

     float emf = (*m1).emEtFraction() ;
     //float hdf = (*m1).etFractionHadronic() ;
     histo2->Fill2a( (*m1).et(), emf, emfCalo, (*m1).sumEt() );

     // calculate the MET from CaloTowers
     LorentzVector calo = CaloMET(caloTowers);
     // calculate the muon correction
     std::vector<double> muPtCorr = ttMuon->MuonEtCorrection(muons);

     // the pure caloMET info
     double vsc =  sqrt( calo.Px()*calo.Px() + calo.Py()*calo.Py() );
     double phic = atan2( calo.Py(), calo.Px() );

     // caloMET + over corrected muon
     double calx1 = calo.Px() + muPtCorr[0];
     double caly1 = calo.Py() + muPtCorr[1];
     double vsc1  =  sqrt( calx1*calx1 + caly1*caly1 );
     double phic1 = atan2(caly1, calx1);

     // find neutrino from generator
     LorentzVector vP4 = METfromNeutrino(genParticles) ;

     double vPT = sqrt(vP4.Px()*vP4.Px() + vP4.Py()*vP4.Py());
     double vPhi= atan2(vP4.Py(), vP4.Px());
     double MET_Res[3] = {100.0, 100.0, 100.0};
     double Phi_Res[3] = {100.0, 100.0, 100.0};

     if ( vPT != 0 ) {
        MET_Res[0] = ((*m1).et() - vPT) / vPT ;
	MET_Res[1] = ( vsc - vPT) / vPT ;
	MET_Res[2] = ( vsc1 - vPT) / vPT ;

	Phi_Res[0] = (*m1).phi() - vPhi;
	Phi_Res[1] =  phic  - vPhi;
        Phi_Res[2] =  phic1 - vPhi;
        histo2->Fill2b( MET_Res[0],MET_Res[1],MET_Res[2],Phi_Res[0],Phi_Res[1],Phi_Res[2] );
     }
 }

}

LorentzVector TtMET::CaloMET( Handle<CaloTowerCollection> calotowers ) {

   double caloXY[2] = {0.0};
   for (CaloTowerCollection::const_iterator t1 = calotowers->begin(); t1 != calotowers->end(); t1++) {
       caloXY[0] -= t1->et()*cos(t1->phi());
       caloXY[1] -= t1->et()*sin(t1->phi());
   }
   double vsumEt = sqrt( caloXY[0]*caloXY[0] + caloXY[1]*caloXY[1] );

   // output MET info
   LorentzVector theMET( caloXY[0], caloXY[1], 0, vsumEt  );

   return theMET;
}

LorentzVector TtMET::CaloMET( const edm::Event& iEvent ) {

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

   double caloXY[2] = {0.0};
   for (CaloTowerCollection::const_iterator t1 = caloTowers->begin(); t1 != caloTowers->end(); t1++) {
       caloXY[0] -= t1->et()*cos(t1->phi());
       caloXY[1] -= t1->et()*sin(t1->phi());
   }
   double vsumEt = sqrt( caloXY[0]*caloXY[0] + caloXY[1]*caloXY[1] );

   // output MET info
   LorentzVector theMET( caloXY[0], caloXY[1], 0, vsumEt  );

   return theMET;
}

LorentzVector  TtMET::METfromObjects( std::vector<const reco::Candidate*> theLep, std::vector<const pat::Jet*> theJets ) {

  double Psum[3] ={ 0.0 };
  for (size_t i=0; i < theLep.size(); i++ ) {
      Psum[0] -= theLep[i]->p4().Px() ;
      Psum[1] -= theLep[i]->p4().Py() ;
  }
  for (size_t i=0; i < theJets.size(); i++ ) {
      Psum[0] -= theJets[i]->p4().Px() ;
      Psum[1] -= theJets[i]->p4().Py() ;
  }

  Psum[2] = sqrt( Psum[0]*Psum[0] + Psum[1]*Psum[1] );

  LorentzVector theMET( Psum[0], Psum[1], 0., Psum[2] );
  
  return theMET;

}

// return 4 momentum of neutrino
LorentzVector TtMET::METfromNeutrino( const edm::Event& iEvent ) {

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   double vp[4] = {0.0};
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       if ( abs((*it).pdgId()) != 12 && abs((*it).pdgId()) != 14 && abs((*it).pdgId()) != 16 ) continue;
       if ( it->status() != 3 ) continue;

       bool fromW = false;
       for (size_t v=0; v< (*it).numberOfMothers(); ++v) {
           const reco::Candidate *mom = (*it).mother(v) ;
           if( abs(mom->pdgId()) != 24 && abs(mom->pdgId()) != 6 ) continue;
           fromW = true;
       }
       if ( !fromW ) continue;
 
       vp[0] +=  it->px() ;
       vp[1] +=  it->py() ;
       vp[2] +=  it->pz() ;
       vp[3] +=  it->p() ;
   }

   LorentzVector vp4 = LorentzVector(vp[0],vp[1],vp[2],vp[3]);
   return vp4;
}

LorentzVector TtMET::METfromNeutrino( Handle<std::vector<reco::GenParticle> > genParticles ) {

   double vp[4] = {0.0};
   for (std::vector<reco::GenParticle>::const_iterator it = genParticles->begin(); it != genParticles->end(); it++ ){
       if ( abs((*it).pdgId()) != 12 && abs((*it).pdgId()) != 14 && abs((*it).pdgId()) != 16 ) continue;
       if ( it->status() != 3 ) continue;

       bool fromW = false;
       for (size_t v=0; v< (*it).numberOfMothers(); ++v) {
           const reco::Candidate *mom = (*it).mother(v) ;
           if( abs(mom->pdgId()) != 24 && abs(mom->pdgId()) != 6 ) continue;
           fromW = true;
       }
       if ( !fromW ) continue;
 
       vp[0] +=  it->px() ;
       vp[1] +=  it->py() ;
       vp[2] +=  it->pz() ;
       vp[3] +=  it->p() ;
   }

   LorentzVector vp4 = LorentzVector(vp[0],vp[1],vp[2],vp[3]);
   return vp4;
}

void TtMET::METandNeutrino( std::vector<const reco::Candidate*> theLep, std::vector<const pat::Jet*> theJets,
                            Handle<std::vector<pat::MET> > met, Handle<std::vector<reco::GenParticle> > genParticles, 
                            HOBJ2* histo ) {

     LorentzVector noMET(0.,0.,0.,0.);
     LorentzVector patMET = (met->size() > 0) ?  (*met)[0].p4() : noMET ;
     LorentzVector evtMET = METfromObjects( theLep, theJets );
     LorentzVector neuMET = METfromNeutrino( genParticles );

     double patResol = 9. ;
     double evtResol = 9. ;
     double theResol = 9. ;
     double dPhi_neu_pat = 9. ;
     double dPhi_neu_evt = 9. ;
     double dPhi_pat_evt = 9. ;

     if ( neuMET.Pt() != 0. ) {
        patResol = ( patMET.Pt()/neuMET.Pt() ) - 1. ;
	evtResol = ( evtMET.Pt()/neuMET.Pt() ) - 1. ;
        if ( patMET.Pt() != 0. && evtMET.Pt() != 0. ) {
	   dPhi_neu_pat = tools->getdPhi( patMET, neuMET );
	   dPhi_neu_evt = tools->getdPhi( evtMET, neuMET );
        }
        if ( patResol >= 1.98 ) patResol = 1.98 ;
        if ( evtResol >= 1.98 ) evtResol = 1.98 ;
        if ( dPhi_neu_pat >= 3.142 ) dPhi_neu_pat = 3.142 ;
        if ( dPhi_neu_evt >= 3.142 ) dPhi_neu_evt = 3.142 ;
     }
     if ( patMET.Pt() != 0. ) {
        theResol = ( evtMET.Pt()/patMET.Pt() ) - 1. ;
        dPhi_pat_evt = tools->getdPhi( patMET, evtMET );
     }

     histo->Fill_2a( patResol,  evtResol, theResol, dPhi_neu_pat, dPhi_neu_evt, dPhi_pat_evt );

}



void TtMET::MetAndMuon(Handle<std::vector<pat::MET> > met, std::vector<const reco::Candidate*> isoMu, HTOP2* histo2, int njets ) {

     if ( met->size() >  0 && isoMu.size() > 0) {

        double df = tools->getdPhi( (*met)[0].p4(), isoMu[0]->p4() );

        if ( njets >= 4 ) histo2->Fill2c0( (*met)[0].et(), df );
        if ( njets == 4 ) histo2->Fill2c1( (*met)[0].et(), df );
        if ( njets == 5 ) histo2->Fill2c2( (*met)[0].et(), df );
        if ( njets == 6 ) histo2->Fill2c3( (*met)[0].et(), df );
        if ( njets >= 7 ) histo2->Fill2c4( (*met)[0].et(), df );

        histo2->Fill2d( (*met)[0].et(), isoMu[0]->pt() );
 
     }

}

void TtMET::MetAndJets(Handle<std::vector<pat::MET> > met, std::vector<const pat::Jet*> theJets, HTOP2* histo2 ) {

   if ( met->size() >  0 && theJets.size() > 1) {

      LorentzVector J12 = theJets[0]->p4() + theJets[1]->p4();
      double df1  = tools->getdPhi( (*met)[0].p4(), theJets[0]->p4() );
      double df12 = tools->getdPhi( (*met)[0].p4(), J12 );

      histo2->Fill2e( (*met)[0].et(),  df1, df12 );

   }
}

