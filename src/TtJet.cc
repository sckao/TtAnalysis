// -*- C++ -*-
//
// Package:    TtJet
// Class:      TtJet
// 
/**\class TtJet TtJet.cc PhysicsTools/TtAnalysis/src/TtJet.cc

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
#include "TtJet.h"

#include "FWCore/Framework/interface/MakerMacros.h"

// constants, enums and typedefs

// static data member definitions


// constructors and destructor
using namespace edm;
using namespace std;
TtJet::TtJet(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  /*
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  rootFileName      = iConfig.getUntrackedParameter<string> ("rootFileName");
  leptonFlavour     = iConfig.getParameter<std::string>   ("leptonFlavour");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  */

  muonSrc         = iConfig.getParameter<edm::InputTag> ("muonSource");
  caloSrc         = iConfig.getParameter<edm::InputTag> ("caloSource");
  genSrc          = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  jptSrc          = iConfig.getParameter<edm::InputTag> ("jptSource");

  bCut            = iConfig.getUntrackedParameter<double> ("bTagCut");
  bTagAlgo        = iConfig.getUntrackedParameter<string> ("bTagAlgo");

  // get the associator parameters
  //tkParas =  iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  //theParameters.loadParameters( tkParas );

  tools        = new TtTools();
  ttMuon       = new TtMuon();
  JetMatching  = new TtMCMatching();

}


TtJet::~TtJet()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtJet Analysis] Destructor called" << endl;
   delete tools;
   delete ttMuon;
   delete JetMatching;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------
static bool EtDecreasing(const pat::Jet s1, const pat::Jet s2) { return ( s1.et() > s2.et() ); }
//static bool EtDecreasing1(const pat::Jet* s1, const pat::Jet* s2) { return ( s1->et() > s2->et() ); }
static bool EtDecreasing2(const reco::Candidate* s1, const reco::Candidate* s2) { return ( s1->et() > s2->et() ); }


void TtJet::JetTreeFeeder(Handle<std::vector<pat::Jet> > patJet, NJet* jtree, int eventId ) {

   for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++)
   {
       float emE0 = (*j1).emEnergyInEB()  + (*j1).emEnergyInEE()  + (*j1).emEnergyInHF() ;
       float hdE0 = (*j1).hadEnergyInHB() + (*j1).hadEnergyInHE() + (*j1).hadEnergyInHF()+ 
                    (*j1).hadEnergyInHO();

       jtree->FillBpatJ( eventId, j1->eta(), j1->phi(), emE0, hdE0, j1->p(), j1->pt() );

   }
}

void TtJet::jetAnalysis(Handle<std::vector<pat::Jet> > patJet, HTOP1* histo1){

   /// 1) overall jet information
   int nJets = patJet->size();
   for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++)
   {
       // Uncorrect Information
       // uncorrected pt from reco jet 
       // raw energy deposite in ECal and HCal
       float  emE0 = (*j1).emEnergyInEB() + (*j1).emEnergyInEE() + (*j1).emEnergyInHF() ;
       float  hdE0 = (*j1).hadEnergyInHB()+(*j1).hadEnergyInHE() + (*j1).hadEnergyInHF() + (*j1).hadEnergyInHO();
       double totalE  = emE0 +hdE0 ;
       double totalEt = totalE*sin(j1->theta() ) ;
       float  emF0 = emE0 / totalE ;

       // pat information -> suppose to be corrected 
       double corrPt = j1->pt();
       float  emF  = (*j1).emEnergyFraction() ;
       double EovH = EoverH(*j1);

       int nCon = (*j1).nConstituents() ;
       double rCon = NofJetConstituents( *j1 ) ;
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;

       if ( nCon == 0 || totalEt < 5. ) continue;

       histo1->Fill1a(  corrPt, (*j1).et(), EovH, nCon,(*j1).n60()/nCon, rCon, (*j1).towersArea(), assTk.size(), emF, emF0 ) ;

       if ( (*j1).pt() > 20.0 )  {
          histo1->Fill1h( EovH, (*j1).towersArea(), nCon, assTk.size() ,rCon, (*j1).n60()/nCon );
       }

       if ( (*j1).towersArea()/(*j1).pt() > 0.006 ) histo1->Fill1c( nCon, assTk.size(), EovH ); 

       if ( (*j1).towersArea()/(*j1).pt() < 0.004 ) histo1->Fill1d( nCon, assTk.size(), EovH );


   }
   histo1->Fill1f( nJets );

}


// leptonic w and jets
void TtJet::JetAndLepW( std::vector<const reco::Candidate*> jet_temp,  LorentzVector p1, HTOP1* histo1){

   int njet = 0;

   //double wEta    = tools->getEta( p1.Px(), p1.Py(), p1.Pz() );
   double wEta    = tools->getY( p1 );
   histo1->Fill1p( njet, wEta );

   if ( jet_temp.size() > 2) {


      // look at fisrt 3 jets
      double dPhi[3]= {-1.};
      double dR[3]  = {-1.};
      double j3[5] = { 0. };
      for (int i= 0; i<3; i++) {
          LorentzVector p2 = jet_temp[i]->p4() ;

          dPhi[i] = tools->getdPhi( p1, p2 );
          double jEta = tools->getY( p2 );
          double dEta = jEta - wEta ;

          dR[i] = sqrt( (dEta*dEta) + (dPhi[i]*dPhi[i]) );
          //if ( i > 1 ) continue ;
          j3[0] += p2.E() ;
          j3[1] += p2.Px() ;
          j3[2] += p2.Py() ;
          j3[3] += p2.Pz() ;
          j3[4] += jet_temp[i]->et() ;
      }

      LorentzVector J12( j3[1], j3[2], j3[3], j3[0] );

      double m3 = sqrt( (j3[0]*j3[0]) - (j3[1]*j3[1]) - (j3[2]*j3[2]) -  (j3[3]*j3[3]) );
      double dRwj = tools->getdR( p1, J12 );
      double dFwj = tools->getdPhi( p1, J12 );

      histo1->Fill1n( njet, wEta, dPhi[0], dPhi[1], dPhi[2], dR[0], dR[1], dR[2], dRwj, dFwj, m3 );
  
   }

}


void TtJet::genJetInfo(Handle<std::vector<reco::GenJet> > genJet, 
                       Handle<std::vector<reco::GenParticle> > genParticles, 
                       HTOP1* histo1, HTOP7* histo7,  HTOP8* histo8){

   std::vector<iTt> ttobjs = JetMatching->TtObjects(genParticles);
   // matching for bjets & wjets
   std::vector<int> bjets; 
   std::vector<int> wjets; 
   for ( size_t i=0; i<ttobjs.size(); i++ ) {
       if ( ttobjs[i].pdgId == 5 ) {
          double dR0 = 99.;
          double ptRes0 =99.;
          int bj = -1;
          for (size_t j=0; j < genJet->size(); j++) {
              bool matched = JetMatching->matchingGeneral( (*genJet)[j].p4(), ttobjs[i], dR0, ptRes0);
              if ( matched ) bj = static_cast<int>(j) ;
          }
          if ( bj != -1  ) bjets.push_back( bj );
       }
       if ( ttobjs[i].pdgId == -5 ) {
          double dR0 = 99.;
          double ptRes0 =99.;
          int bj = -1;
          for (size_t j=0; j < genJet->size(); j++) {
              bool matched = JetMatching->matchingGeneral( (*genJet)[j].p4(), ttobjs[i], dR0, ptRes0);
              if ( matched ) bj = static_cast<int>(j) ;
          }
          if ( bj != -1  ) bjets.push_back( bj );
       }
       if ( abs(ttobjs[i].pdgId) < 5 ) {
          double dR0 = 99.;
          double ptRes0 =99.;
          int wj = -1;
          for (size_t j=0; j < genJet->size(); j++) {
              bool matched = JetMatching->matchingGeneral( (*genJet)[j].p4(), ttobjs[i], dR0, ptRes0);
              if ( matched ) wj = static_cast<int>(j) ;
          }
          if ( wj != -1  ) wjets.push_back( wj );
       }
   }
   for (std::vector<int>::iterator i= bjets.begin(); i!= bjets.end(); i++ ) {
       double EovH = -1. ;
       if ( (*genJet)[*i].hadEnergy() != 0 ) EovH = (*genJet)[*i].emEnergy()/(*genJet)[*i].hadEnergy() ;
       histo7->Fill7f( (*genJet)[*i].pt(), (*genJet)[*i].eta(), EovH ) ; 
   }
   for (std::vector<int>::iterator i= wjets.begin(); i!= wjets.end(); i++ ) {
       double EovH = -1. ;
       if ( (*genJet)[*i].hadEnergy() != 0 ) EovH = (*genJet)[*i].emEnergy()/(*genJet)[*i].hadEnergy() ;
       histo8->Fill8f( (*genJet)[*i].pt(), (*genJet)[*i].eta(), EovH ) ; 
   }
   int nj[4] = {0};
   for (size_t j=0; j < genJet->size(); j++) {
       nj[0]++;
       if ( (*genJet)[j].nConstituents() ==0 ) continue;
       if ( (*genJet)[j].pt()  > 20. ) nj[1]++;
       if ( (*genJet)[j].pt()  > 40. ) nj[2]++;
       if (fabs( (*genJet)[j].eta() ) < 3.5 ) nj[3]++;

       double EovH = -1. ;
       if ( (*genJet)[j].hadEnergy() != 0 ) EovH = (*genJet)[j].emEnergy() / (*genJet)[j].hadEnergy() ;

       histo1->Fill1j( (*genJet)[j].eta(), (*genJet)[j].pt(), EovH );

   }
   histo1->Fill1k( nj[0], nj[1], nj[2], nj[3] ) ;
}


void TtJet::matchedWJetsAnalysis( std::vector<jmatch> mcwjets , std::vector<const reco::Candidate*> isoMuons,
                                  Handle<std::vector<pat::Jet> > patJets, HTOP8* histo8 ) {

  // look at matched jet properties
 
  std::vector<LorentzVector> wpjj;
  std::vector<LorentzVector> wnjj;
  double dR_wjj = -1. ; 
  std::vector<LorentzVector> wpqq;
  std::vector<LorentzVector> wnqq;
  double dR_wqq = -1. ; 
  for (std::vector<jmatch>::const_iterator j1 = mcwjets.begin(); j1 != mcwjets.end(); j1++) {

      if ( abs(j1->MomIdx)  == 5 ) continue;
      //pat::Jet truth =  *(j1->trueJet) ;
      const reco::Candidate* truth =  j1->trueJet ;
      reco::Particle jmom = j1->mom;

      int flv =  j1->MomIdx ;
      if ( flv == 2 || flv == 4 || flv == -1 || flv == -3 ) {
         if ( wpjj.size()==0 ) { 
            wpjj.push_back( truth->p4() );
            wpqq.push_back( jmom.p4() );
         } else if ( wpjj.size()==1 ) {
            wpjj.push_back( truth->p4() );
            wpqq.push_back( jmom.p4() );
         } 
      }
      if ( flv == -2 || flv == -4 || flv == 1 || flv == 3 ) { 
         if ( wnjj.size()==0 ) {
            wnjj.push_back( truth->p4() );
            wnqq.push_back( jmom.p4() );
         } else if ( wnjj.size()==1 ) {
            wnjj.push_back( truth->p4() );
            wnqq.push_back( jmom.p4() );
         }
      }

      if ( wpjj.size()==2 ) {
         dR_wjj = tools->getdR( wpjj[0], wpjj[1] );
         dR_wqq = tools->getdR( wpqq[0], wpqq[1] );
         histo8->Fill8j( dR_wjj , dR_wqq );
         wpjj.clear();
         wpqq.clear();
      }
      if ( wnjj.size()==2  ) {
         dR_wjj = tools->getdR( wnjj[0], wnjj[1] );
         dR_wqq = tools->getdR( wnjj[0], wnjj[1] );
         histo8->Fill8j( dR_wjj , dR_wqq );
         wnjj.clear();
         wnqq.clear();
      }
      double res_Pt  = (truth->pt() - jmom.pt()) / jmom.pt() ;

      bool goodTransfer = false;
      const pat::Jet* pat_truth = tools->ReturnJetForm( truth, patJets,  goodTransfer );
      //double jProb   = pat_truth->bDiscriminator("jetProbabilityBJetTags") ;
      //double tkCount = pat_truth->bDiscriminator("trackCountingHighEffBJetTags") ;
      double softMuTag = pat_truth->bDiscriminator( bTagAlgo ) ;

      double dR_WjMu = 999.;
      for (std::vector<const reco::Candidate*>::const_iterator m1= isoMuons.begin(); m1 != isoMuons.end(); m1++) {
          double dR = tools->getdR( pat_truth->p4() , (*m1)->p4() ) ;
          if ( dR < dR_WjMu )   dR_WjMu = dR; 
      }
      if (dR_WjMu != 999. ) { 
         histo8->Fill8c( res_Pt, softMuTag, dR_WjMu );
      }
  }

}

void TtJet::matchedbJetsAnalysis( std::vector<jmatch> mcjets, std::vector<const reco::Candidate*> isoMuons, 
                                  Handle<std::vector<pat::Jet> > patJets, HTOP7* histo7 ) {

  for (std::vector<jmatch>::const_iterator j1 = mcjets.begin(); j1 != mcjets.end(); j1++) {

      if ( abs(j1->MomIdx)  != 5 ) continue;

      //pat::Jet truth =  *(j1->trueJet) ;
      const reco::Candidate* reco_truth =  j1->trueJet ;
      bool goodTransfer = false;
      const pat::Jet* truth = tools->ReturnJetForm( reco_truth, patJets,  goodTransfer );
      reco::Particle jmom = j1->mom;
      
      int nCon = truth->nConstituents();

      double towerArea = 0.;
      double EovH = -1.;
      double emEF = -1.;
      int nTracks = -1;
      int n60 = -1;
      int n90 = -1;
      if ( truth->isCaloJet() ) {
         EovH = EoverH( *truth  );
         emEF = truth->emEnergyFraction();
         towerArea = truth->towersArea();
         edm::RefVector<reco::TrackCollection>  assTk = truth->associatedTracks() ;
         nTracks= assTk.size();
         n60 = truth->n60();
         n90 = truth->n90();
      }

      double Res_Pt = (truth->pt() - jmom.pt()) / jmom.pt() ; 
      double softMuTag = truth->bDiscriminator( bTagAlgo ) ;

      histo7->Fill7a( truth->pt(), truth->eta(), EovH, nCon, n60, n90, towerArea, nTracks, emEF, Res_Pt, softMuTag );

      // dR(closest isoMu, matched bjet )
      double dR_bmu = 9.;
      double RelPt_bmu = -1 ;
      for (std::vector<const reco::Candidate*>::const_iterator m1= isoMuons.begin(); m1 != isoMuons.end(); m1++) {
          double dR = tools->getdR( truth->p4(), (*m1)->p4() );
          if (dR < dR_bmu ) {
             dR_bmu = dR ;
             RelPt_bmu = tools->getRelPt( truth->p4(), (*m1)->p4() );
          }               
      } 
      double dR_bwj = 999.;
      for (std::vector<jmatch>::const_iterator j2 = mcjets.begin(); j2 != mcjets.end(); j2++ ) {
          if ( abs(j1->MomIdx) == 5 ) continue;
          double dh = truth->eta() - (j2->trueJet)->eta();
          double df = truth->phi() - (j2->trueJet)->phi();
          double dR = sqrt( dh*dh + df*df );
          if (dR < dR_bwj ) {
             dR_bwj = dR ;
          }               
      }
      if ( dR_bmu != 999. ) {  histo7->Fill7b( dR_bmu, RelPt_bmu );  }
      if ( dR_bwj != 999. ) {  histo7->Fill7c( dR_bwj );  }
  }

}

void TtJet::bTagAnalysis(Handle<std::vector<pat::Jet> > Jets, Handle<std::vector<pat::Muon> > patMu, HTOP7* histo7  ) {

  for (std::vector<pat::Jet>::const_iterator j1 = Jets->begin(); j1 != Jets->end(); j1++) {

      double bDis = j1->bDiscriminator( bTagAlgo ) ;

      if ( fabs(j1->eta()) > 2.7 ) continue; 
      if ( j1->et() < 20. ) continue; 

      histo7->Fill7g( bDis );

      if ( bDis < bCut ) continue;

      double dR = 999.;
      double d0 = 999.;
      double RelPt_bMu = 999.;
      int x = -1 ;
      int u = -1 ;
      for (std::vector<pat::Muon>::const_iterator m1 = patMu->begin(); m1!= patMu->end(); m1++) {
           x++;
           //bool isolated = ttMuon->IsoMuonID( *m1, 0.9 ) ;
           if ( !(m1->isTrackerMuon()) ) continue;
           double dR1 = tools->getdR( j1->p4() , m1->p4() ) ;
           if ( dR1 < dR ) {
              dR = dR1 ;
              if ( dR < 0.5 ) {
                 u = x ;
                 reco::TrackRef muTrack = m1->innerTrack() ;
                 d0 = muTrack->d0() ;
                 RelPt_bMu = tools->getRelPt( j1->p4(), m1->p4() );
              }
           }
      }
      if ( u != -1 ) histo7->Fill7d( dR, RelPt_bMu , j1->phi(), d0 );

      /*
      cout<<" bDis(softMu)   "<< j1->bDiscriminator("softMuonBJetTags")<<endl;
      cout<<" bDis(jetProb)  "<< j1->bDiscriminator("jetProbabilityBJetTags")       <<endl;
      cout<<" bDis(trkCount3)"<< j1->bDiscriminator("trackCountingHighPurBJetTag")  <<endl;
      cout<<" bDis(trkEff)   "<< j1->bDiscriminator("trackCountingHighEffBJetTags") <<endl;
      cout<<" bDis(softEle)  "<< j1->bDiscriminator("softElectronBJetTags")         <<endl;
      cout<<" bDis(2ndVtx)   "<< j1->bDiscriminator("impleSecondaryVertexBJetTags") <<endl;
      cout<<" bDis(cmbVtx)   "<< j1->bDiscriminator("combinedSVBJetTags")           <<endl;
      */
  }

}

void TtJet::selectedWJetsAnalysis(Handle<std::vector<pat::Jet> > patJet, std::vector<const reco::Candidate*> isoMuons, HTOP8* histo8  ) {

   std::vector<const reco::Candidate*> wjets = WJetSelection( patJet, isoMuons, 30. );

   // jets selection
   double jPt[2]  = {0.0, 0.0};
   double jPhi[2] = {-9.0, -9.0};
   double jEta[2] = {-9.0, -9.0};

   int idx=0;
   for (std::vector<const reco::Candidate*>::const_iterator j1 = wjets.begin(); j1 != wjets.end(); j1++)
   {
       if ( idx < 2){
          jPt[idx]  = (*j1)->pt() ;
          jPhi[idx] = (*j1)->phi() ;
          jEta[idx] = (*j1)->eta() ;
       }
       idx++;
   }
   histo8->Fill8a( jPt[0], jPt[1], jEta[0], jEta[1], wjets.size() );

   // rebuild W from selected jets!
   int jj1 =0;
   for (std::vector<const reco::Candidate*>::const_iterator j1 = wjets.begin(); j1 != wjets.end(); j1++)
   {
       // only semi-lep
       //if ( topo != 1 || !pass ) continue;
       jj1++ ;
       int jj2 =0;

       for (std::vector<const reco::Candidate*>::const_iterator j2 = wjets.begin(); j2 != wjets.end(); j2++)
       {
           jj2++;
           if (j2 <= j1) continue;

           //cout<<" flavor1: "<< (*j1)->partonFlavour()<<"  vtx:"<<(*j1)->vertex() <<endl;
           //cout<<" flavor2: "<< (*j2)->partonFlavour()<<"  vtx:"<<(*j2)->vertex() <<endl;
           //LorentzVector qm1 = (*j1)->energyFractionHadronic() * (*j1)->p4();
           //LorentzVector qm2 = (*j2)->energyFractionHadronic() * (*j2)->p4();
           LorentzVector qm1 = (*j1)->p4();
           LorentzVector qm2 = (*j2)->p4();
           LorentzVector pW = findW( qm1, qm2);
           double momW = sqrt( pW.Px()*pW.Px() + pW.Py()*pW.Py() + pW.Pz()*pW.Pz() );
           double massW = sqrt ( pW.E()*pW.E() - pW.Px()*pW.Px() - pW.Py()*pW.Py() - pW.Pz()*pW.Pz() );
           //cout<<" combination of "<<jj1<<" & "<<jj2<<endl;
           //cout<<" find W = "<< pW.Px() <<","<< pW.Py() <<","<< pW.Pz() <<","<< pW.P() <<endl;
           //cout<<" mass W = "<< massW <<endl;
           //if (!exclude) {
           histo8->Fill8b( massW, momW );
       }
   }

}


double TtJet::NofJetConstituents( pat::Jet theJet ) {

    std::vector< const reco::Candidate* >  jc =  theJet.getJetConstituentsQuick();

    //for (size_t i=0; i < jc.size(); i++ ) {
    double sumE =0 ;
    for (size_t i=0; i < jc.size(); i++ ) {
        //const reco::Candidate* jj = dynamic_cast<const reco::Candidate*>( &*(jc[i]) );
        //cout<<" jet constituents "<< jc[i]->p4() ;
        //cout<<" h:"<<jc[i]->eta()<<" f:"<<jc[i]->phi()<<endl;
        sumE = sumE + jc[i]->energy() ;
    }

    double sumP =0;
    edm::RefVector<reco::TrackCollection>  assTk = theJet.associatedTracks() ;
    for(size_t i=0; i < assTk.size(); i++ ) {
        sumP = sumP + (*assTk[i]).p() ;
        //cout <<" trk q:" <<(*assTk[i]).charge()<<" p4:  "<<(*assTk[i]).momentum() ;
        //cout <<" p: "<< (*assTk[i]).p() <<endl ;
    }

    if ( theJet.nConstituents() > 0) {
       double nC = static_cast<double>( theJet.nConstituents() );
       double nT = static_cast<double>( assTk.size() ) ;
       //cout <<" ratio = "<< nC/nT <<endl;
       return nT/nC ;
    } else {
       return -1 ;
    }
}

double TtJet::EoverH( pat::Jet theJet ) {

       float emE = theJet.emEnergyFraction() ;
       float hdE = theJet.energyFractionHadronic() ;
       double EovH = -1.0 ;
       if ( hdE != 0.0 ) {
          EovH = emE/hdE ;
       }
       return EovH;
}

LorentzVector TtJet::findW(LorentzVector qm1, LorentzVector qm2 ) {

     //exclude = false;
     double m1 = qm1.P();
     double m2 = qm2.P();
     double xW = qm1.Px() + qm2.Px() ;
     double yW = qm1.Py() + qm2.Py() ;
     double zW = qm1.Pz() + qm2.Pz() ;
     double EW = m1 + m2 ;
     LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
     //cout<<" pW: "<<xW<<","<<yW<<","<<zW<<","<<EW<<endl;
     //cout<<" mW: "<<mW.Px()<<","<<mW.Py()<<","<<mW.Pz()<<","<<mW.E()<<endl;
     return mW;
}

// using track associator
/*
void TtJet::JetMatchedMuon( Handle<std::vector<pat::Jet> > patJet , Handle<std::vector<pat::Muon> > patMuon, 
                            const Event& iEvent, const EventSetup& iSetup, const edm::ParameterSet& iConfig,  HTOP3* histo3){


  // get the associator parameters
  ParameterSet tkParas =  iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  TrackAssociatorParameters theParameters;
  theParameters.loadParameters( tkParas );

  TrackDetectorAssociator trackAssociator;
  trackAssociator.useDefaultPropagator();

  for (std::vector<pat::Muon>::const_iterator m1 = patMuon->begin(); m1 != patMuon->end(); m1++)
  {
       float vx = (m1->vertex()).x();
       float vy = (m1->vertex()).y();
       float vz = (m1->vertex()).z();
       float mx   = (m1->momentum()).x();
       float my   = (m1->momentum()).y();
       float mz   = (m1->momentum()).z();
       GlobalPoint gp = GlobalPoint(vx, vy, vz) ;
       GlobalVector gm = GlobalVector(mx, my, mz );
       TrackDetMatchInfo Info = trackAssociator.associate(iEvent, iSetup, gm, gp, m1->charge(), theParameters);

       DetId det_Id =  Info.findMaxDeposition( TrackDetMatchInfo::TowerTotal ) ;
       CaloTowerDetId caloId(det_Id); 
       //cout<<" DetId:"<<  caloId <<" Energy:"<<  Info.nXnEnergy(det_Id,TrackDetMatchInfo::TowerTotal,1) << endl;

       // flag muon isolation
       bool isoMuon = false;
       //reco::MuonIsolation R3 = m1->getIsolationR03();
       //reco::MuonIsolation R5 = m1->getIsolationR05();
       if ( ( m1->trackIso() == 1. )&&( m1->ecalIso() <= 3 )&& ( m1->hcalIso() <= 1 ) ) {
          isoMuon = true;
       }

       for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++)
       {
           float dh = j1->eta() - m1->eta();
           float df = j1->phi() - m1->phi();
           float dR = sqrt( (dh*dh) + (df*df) );
           if ( dR > 0.5 ) continue;

           CaloTowerPtr jcaloPtr;
	   for(int i=0; i< j1->nConstituents(); i++  ){
              jcaloPtr = j1->getCaloConstituent(static_cast<unsigned int>(i));
	      //cout<<"    JetTower"<<i<<" Id:"<< (*jcaloRef).id()  <<" E:"<< (*jcaloRef).energy() <<endl;
              if ( caloId != (*jcaloPtr).id() ) continue;
              if (isoMuon) {
                 histo3->Fill3j( (*jcaloPtr).energy(), m1->p(), j1->energy() );
              } else {
                 histo3->Fill3k( (*jcaloPtr).energy(), m1->p(), j1->energy() );
              }
           }
       }

  }

}
*/

// Using SteppingHelixPropogator 
void TtJet::JetMatchedMuon( Handle<std::vector<pat::Jet> > patJet , Handle<std::vector<pat::Muon> > patMuon, 
                            const edm::Event& iEvent, const EventSetup& eventSetup, HTOP3* histo3, bool Done ){

   // Magnetic field
   ESHandle<MagneticField> field;
   eventSetup.get<IdealMagneticFieldRecord>().get(field);
   ESHandle<Propagator> shPropAny;
   eventSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", shPropAny);
   ESHandle<CaloGeometry> caloGeomHandle;
   eventSetup.get<CaloGeometryRecord>().get(caloGeomHandle);
   //eventSetup.get<IdealGeometryRecord>().get(caloGeomHandle);
   const CaloGeometry& caloGeom = *caloGeomHandle;

   // method to use SteppingHelixPropagator       
   for (std::vector<pat::Muon>::const_iterator m1 = patMuon->begin(); m1 != patMuon->end(); m1++)
   {
       float vx = (m1->vertex()).x();
       float vy = (m1->vertex()).y();
       float vz = (m1->vertex()).z();
       float mx   = (m1->momentum()).x();
       float my   = (m1->momentum()).y();
       float mz   = (m1->momentum()).z();

       GlobalPoint gp = GlobalPoint(vx, vy, vz) ;
       GlobalVector gm = GlobalVector(mx, my, mz );
       AlgebraicSymMatrix66 covT;
       FreeTrajectoryState ftsStart = tools->getFTS(gp, gm, m1->charge(), covT, &*field);

       // flag muon isolation
       bool isoMuon = false;
       //reco::MuonIsolation R3 = m1->getIsolationR03();
       //reco::MuonIsolation R5 = m1->getIsolationR05();
       if ( ( m1->trackIso() == 1. )&&( m1->ecalIso() <= 3 )&& ( m1->hcalIso() <= 1 ) ) {
          isoMuon = true;
       }

       // looping the jets 
       for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++)
       {
           float dh = j1->eta() - m1->eta();
           float df = j1->phi() - m1->phi();
           float dR = sqrt( (dh*dh) + (df*df) );
           if ( dR > 0.5 ) continue;

           float dmut = 9998.;
           int muonhit = -1;
           //looping all its constituents => calotowers
           CaloTowerPtr jcaloRef;
           for(int i=0; i< j1->nConstituents(); i++  ){
              jcaloRef = j1->getCaloConstituent(i);

              //cout<<"JETConstitents"<<i<<" Id:"<< (*jcaloRef).id()  <<" E:"<< (*jcaloRef).energy();
              //cout<<"  h:"<<(*jcaloRef).eta()<<" f:"<<(*jcaloRef).phi();
              //cout<<" emE:"<<(*jcaloRef).emEnergy()<<" hadE:"<<(*jcaloRef).hadEnergy()<<endl;
              //cout<<"   # of subConstitents "<< (*jcaloRef).constituentsSize() << endl;
  
              float dl1 = 9999.;
              //looping all cells  => ECal or HCal cells
              for(size_t k=0; k< (*jcaloRef).constituentsSize(); k++  ){ 
                 GlobalPoint caloGP = caloGeom.getPosition(  (*jcaloRef).constituent(k) );
                 //double R = sqrt( (caloGP.x()*caloGP.x()) + (caloGP.y()*caloGP.y()) );
                 //cout<<"         GP:"<<caloGP<<"  R="<<R<<" h:"<<tools->getEta(caloGP.x(), caloGP.y(), caloGP.z())<<endl;
                // setup the final destination 
		 FreeTrajectoryState ftsDest1;
                 GlobalPoint pDest1 = caloGP;
		 const SteppingHelixPropagator* shPropAnyCPtr1 = dynamic_cast<const SteppingHelixPropagator*>(&*shPropAny);
		 ftsDest1 = shPropAnyCPtr1->propagate(ftsStart, pDest1);
		 GlobalPoint fGP = ftsDest1.position();
		 float dx = fGP.x() - pDest1.x() ;
		 float dy = fGP.y() - pDest1.y() ;
		 float dz = fGP.z() - pDest1.z() ;
		 double dl = sqrt( (dx*dx) + (dy*dy) + (dz*dz) ); 
		 //cout<<"         final gp"<<gp<<"=> "<<ftsDest1.position();
                 //cout<<" h:"<<tools->getEta(fGP.x(), fGP.y(), fGP.z()) << endl;
                 //cout<<"         muon E:"<< m1->p()<<endl;
		 //cout<<"         *** dL = "<< dl <<endl;
                 if (dl < dl1 ) dl1 = dl;         
              }
              if (dl1 >= dmut ) continue;
              dmut = dl1;
              muonhit = i ; 
           }
           if ( dmut > 3.0 ) continue;
           if (muonhit == -1 ) continue; 

           jcaloRef = j1->getCaloConstituent( muonhit );
           if (isoMuon  ) {
              histo3->Fill3j( (*jcaloRef).energy(), m1->p(), j1->energy() );
           } else {
              histo3->Fill3k( (*jcaloRef).energy(), m1->p(), j1->energy() );
           }

       }
   }
   Done = true ;
}

void TtJet::JetTrigger( Handle<std::vector<pat::Jet> > jets, Handle <edm::TriggerResults> triggers ) {

   for (std::vector<pat::Jet>::const_iterator it = jets->begin(); it!= jets->end(); it++) {
       std::vector<pat::TriggerPrimitive> trigInfo = it->triggerMatches() ;
       cout<<" === Jet trigger ==== "<<endl;
       for(size_t i=0; i< trigInfo.size(); i++) {
          cout<<"    ObjId:"<< trigInfo[i].triggerObjectId() ;
          cout<<"  ObjType:"<< trigInfo[i].triggerObjectType() <<endl;
          cout<<" filter name:"<<trigInfo[i].filterName() <<endl;
       }
   }
   /*
   edm::TriggerNames trigNames( *triggers );
   for (size_t i=0; i< triggers->size(); i++ ) {

       string triggered = triggers->accept(i) ? "Yes" : "No" ;

       cout<<" path("<<i<<") accepted ? "<< triggered ;
       cout<<" trigName: "<< trigNames.triggerName(i)<<endl;
   }*/

}

// JPT Correction
std::vector<const reco::Candidate*> TtJet::JetSelection( Handle<std::vector<reco::CaloJet> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold, HOBJ1* histo, std::vector<bool>* bTags ) {

   std::vector<const reco::Candidate*> jet_temp;

   for (std::vector<reco::CaloJet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( fabs(j1->eta()) > 2.7 ) continue; 
       if ( j1->et() < EtThreshold ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ 
                     j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;
       /*
       if ( bTags != NULL ) {
          double bDis = j1->bDiscriminator( "softMuonBJetTags" ) ;
          if ( bDis >= 0.9 ) (*bTags).push_back( true );
          if ( bDis <  0.9 ) (*bTags).push_back( false );
       }
       */

       jet_temp.push_back( &*j1 );
   }

   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtDecreasing2 );

   if ( histo != NULL ) { 
      int jetSize = static_cast<int>( jets->size() );
      if (jetSize > 20) jetSize = 20 ;
      // Cannot tag JPT, => set number of BJets = 0
      histo->Fill_1g( jetSize, jet_temp.size(), 0 );
   }  

   return jet_temp;

}
// gen jets
std::vector<const reco::Candidate* > TtJet::JetSelection( Handle<std::vector<reco::GenJet> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold, std::vector<bool>* bTags ) {

   std::vector<const reco::Candidate* > jet_temp ;
   for (std::vector<reco::GenJet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( fabs(j1->eta()) > 2.7 ) continue; 
       if ( j1->et() < EtThreshold ) continue; 

       double calE = j1->emEnergy() ;
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       if ( bTags != NULL ) {
          double bDis = abs( j1->pdgId() );
          if ( bDis == 5 ) (*bTags).push_back( true );
          if ( bDis != 5 ) (*bTags).push_back( false );
       }

       jet_temp.push_back( &*j1 );
   }

   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtDecreasing2 );
   return jet_temp;

}

// pat jets
std::vector<const reco::Candidate* > TtJet::JetSelection( Handle<std::vector<pat::Jet> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold, HOBJ1* histo, std::vector<bool>* bTags, string bTagAlgo1 ) {

   int nBJets = 0;
   std::vector<const reco::Candidate* > jet_temp ;
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( fabs(j1->eta()) > 2.7 ) continue; 
       if ( j1->et() < EtThreshold ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ 
                     j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       double bDis = j1->bDiscriminator( bTagAlgo1 ) ;
       if ( bDis >= bCut ) nBJets++ ;
       if ( bTags != NULL ) {
          if ( bDis >= bCut ) (*bTags).push_back( true );
          if ( bDis <  bCut ) (*bTags).push_back( false );
       }
       //if ( histo != NULL ) histo->Fill_1( bDis );
       jet_temp.push_back( &*j1 );
   }

   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtDecreasing2 );

   if ( histo != NULL ) { 
      int jetSize = static_cast<int>( jets->size() );
      if (jetSize > 20) jetSize = 20 ;
      histo->Fill_1g( jetSize, jet_temp.size(), nBJets );
   }  
   return jet_temp;

}

std::vector<const reco::Candidate* > TtJet::SoftJetSelection( Handle<std::vector<reco::CaloJet> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold,  HOBJ1* histo ) {

   std::vector<const reco::Candidate* > jet_temp ;
   for (std::vector<reco::CaloJet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( j1->et() >= EtThreshold ) continue; 
       if ( fabs(j1->eta()) > 2.7 ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ 
                     j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       jet_temp.push_back( &*j1 );
   }

   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtDecreasing2 );
   
   double leadingJetEt = ( jet_temp.size() > 0 ) ? jet_temp[0]->et() : 0.  ;

   histo->Fill_1i( jet_temp.size(), leadingJetEt );
   return jet_temp;

}

std::vector<const reco::Candidate* > TtJet::SoftJetSelection( Handle<std::vector<pat::Jet> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold, std::vector<bool>* bTags, string bTagAlgo, HOBJ1* histo ) {

   std::vector<const reco::Candidate* > jet_temp ;
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {

       if ( j1->et() >= EtThreshold ) continue; 
       if ( fabs(j1->eta()) > 2.7 ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ 
                     j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       if ( calEt < 4. ) continue;

       double EovH = EoverH(*j1);
       if ( EovH == -1. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           LorentzVector muP4 = IsoMuons[i]->p4() ;
           double dR_mu = tools->getdR( muP4, j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       if ( bTags != NULL ) {
          double bDis = j1->bDiscriminator( bTagAlgo ) ;
          if ( bDis >= bCut ) (*bTags).push_back( true );
          if ( bDis <  bCut ) (*bTags).push_back( false );
       }

       jet_temp.push_back( &*j1 );
   }

   if ( jet_temp.size() > 1 ) sort( jet_temp.begin(), jet_temp.end(), EtDecreasing2 );
   if ( histo != NULL ) {
      double leadingJetEt = ( jet_temp.size() > 0 ) ? jet_temp[0]->et() : 0.  ;
      histo->Fill_1i( jet_temp.size(), leadingJetEt );
   }
   return jet_temp;
}

/*
std::vector<const reco::Candidate* > TtJet::SoftJetSelection( Handle<std::vector<pat::Jet> > jets, std::vector<const reco::Candidate*> IsoMuons , double EtThreshold,  HOBJ1* histo, std::vector<bool>* bTags ) {

   std::vector<const reco::Candidate* > jet_temp = SoftJetSelection( jets, IsoMuons, EtThreshold, bTags ); 

   double leadingJetEt = ( jet_temp.size() > 0 ) ? jet_temp[0]->et() : 0.  ;
   histo->Fill_1i( jet_temp.size(), leadingJetEt );

   return jet_temp;
}
*/

std::vector<const reco::Candidate*> TtJet::WJetSelection( Handle<std::vector<pat::Jet> >  jets, std::vector<const reco::Candidate*> IsoMuons, double JetEtCut ) {

   std::vector<const reco::Candidate* > jCollection;
   jCollection.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {
       if ( fabs(j1->eta()) > 2.7 || j1->et() < JetEtCut ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       double EovH = EoverH(*j1);

       if ( EovH == -1. || calEt < 4. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           double dR_mu = tools->getdR( IsoMuons[i]->p4() , j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       double bDis = j1->bDiscriminator( bTagAlgo ) ;
       if ( bDis >= bCut ) continue;

       jCollection.push_back( &*j1 );
   }

   // sort the seeds by # of own segments
   sort(jCollection.begin(), jCollection.end(), EtDecreasing2 ) ;

   return jCollection;

}

std::vector<const reco::Candidate*> TtJet::bJetSelection( Handle<std::vector<pat::Jet> >  jets, std::vector<const reco::Candidate*> IsoMuons, double JetEtCut, string bTagAlgo ) {

   std::vector<const reco::Candidate* > jCollection;
   jCollection.clear();
   for (std::vector<pat::Jet>::const_iterator j1 = jets->begin(); j1 != jets->end(); j1++)
   {
       if ( fabs(j1->eta()) > 2.7 || j1->et() < JetEtCut ) continue; 

       double calE = j1->emEnergyInEB()  + j1->emEnergyInEE()  + j1->emEnergyInHF() +
                     j1->hadEnergyInHB() + j1->hadEnergyInHE() + j1->hadEnergyInHF()+ j1->hadEnergyInHO();
       double calEt = calE * sin( j1->theta() );
       double EovH = EoverH(*j1);

       if ( EovH == -1. || calEt < 4. ) continue;

       bool fakeJet = false;
       for ( size_t i =0; i < IsoMuons.size(); i++ ) {
           double dR_mu = tools->getdR( IsoMuons[i]->p4() , j1->p4() );  
           if ( dR_mu < 0.1 ) fakeJet = true ;
       }
       if ( fakeJet ) continue;

       edm::RefVector<reco::TrackCollection>  assTk = j1->associatedTracks() ;
       if (assTk.size() == 0) continue;

       double bDis = j1->bDiscriminator( bTagAlgo ) ;
       if ( bDis < bCut ) continue;

       jCollection.push_back( &*j1 );
   }

   // sort the seeds by # of own segments
   sort(jCollection.begin(), jCollection.end(), EtDecreasing2 ) ;

   return jCollection;
}

// isolated muon and jets
//void TtJet::MuonAndJet( std::vector<reco::Candidate> jet_temp,  LorentzVector p1, HTOP1* histo1){
void TtJet::MuonAndJet( std::vector<const reco::Candidate*> jet_temp,  const reco::Candidate* isoMu, HOBJ1* histo1){

   int njet = static_cast<int>( jet_temp.size() );
   //double muEta    = tools->getEta( p1.Px(), p1.Py(), p1.Pz() );
   double muEta     = isoMu->eta() ;
   LorentzVector p1 = isoMu->p4() ;

   // njets vs eta of muon 
   histo1->Fill_1d( njet, muEta );

   double dR = 9. ;
   double relPt = -1;
   for (int i= 0; i < njet ; i++) {
       LorentzVector p2 = jet_temp[i]->p4() ;

       double dR1    = tools->getdR( p1, p2 );
       double relPt1 = tools->getRelPt(p1, p2);
       if ( dR1 < dR ) {
          dR = dR1; 
          relPt = relPt1 ;
       }

       double dPhi[3] = { -1. };
       if ( jet_temp.size() > 2 ) {
          dPhi[i] = tools->getdPhi(p1, p2);
          histo1->Fill_1e(njet, dPhi[0], dPhi[1], dPhi[2]);
       }
   }
   // min dR(muon, jet)
   histo1->Fill_1b( njet, dR, relPt );

}

void TtJet::JetEtSpectrum( std::vector<const reco::Candidate*> jet_temp, HOBJ1* histo1){

   histo1->Fill_1f( jet_temp.size() );
   if ( jet_temp.size() > 3) {
      // get m3
      std::vector<LorentzVector> vlist;
      vlist.push_back( jet_temp[0]->p4() );
      vlist.push_back( jet_temp[1]->p4() );
      vlist.push_back( jet_temp[2]->p4() );
      double m3 = tools->getInvMass( vlist );

      histo1->Fill_1c(jet_temp[0]->et(), jet_temp[1]->et(), jet_temp[2]->et(), jet_temp[3]->et(), m3 ) ;
   }
}

void TtJet::JetEtSpectrum(Handle<std::vector<reco::GenJet> > genJet, HOBJ1* histo1){

   std::vector<reco::GenJet> jet_temp ;
   for (std::vector<reco::GenJet>::const_iterator j1 = genJet->begin(); j1 != genJet->end(); j1++)
   {
       jet_temp.push_back( *j1 );
   }

   histo1->Fill_1f( jet_temp.size() );
   if ( jet_temp.size() > 3) {
      sort(jet_temp.begin(), jet_temp.end(), EtDecreasing );
      std::vector<LorentzVector> vlist;
      vlist.push_back( jet_temp[0].p4() );
      vlist.push_back( jet_temp[1].p4() );
      vlist.push_back( jet_temp[2].p4() );
      double m3 = tools->getInvMass( vlist );

      histo1->Fill_1c(jet_temp[0].et(), jet_temp[1].et(), jet_temp[2].et(), jet_temp[3].et(), m3 ) ;
   }

}

void TtJet::JetdRAnalysis( std::vector<const reco::Candidate*>  patJet, HOBJ1* histo1 ){

    int jetSize = static_cast<int>( patJet.size() ) ;
    for(int i=0; i< jetSize; i++) {

       if ( i+1 >= jetSize  ) continue;
       for (int j=i+1; j< jetSize; j++) {

           double dR = tools->getdR( patJet[i]->p4(), patJet[j]->p4() );
           histo1->Fill_1a( dR );
       }
    }

}
