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
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  */

  muonSrc         = iConfig.getParameter<edm::InputTag> ("muonSource");
  caloSrc         = iConfig.getParameter<edm::InputTag> ("caloSource");
  genSrc          = iConfig.getParameter<edm::InputTag> ("genParticles"); 

  // get the associator parameters
  tkParas =  iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  theParameters.loadParameters( tkParas );

  fromTtMuon      = new TtMuon();
  JetMatching     = new TtMCMatching();
  evtSelected     = new TtEvtSelector();
}


TtJet::~TtJet()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //if (debug) cout << "[TtJet Analysis] Destructor called" << endl;
   delete fromTtMuon;
   delete JetMatching;
   delete evtSelected;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;

// ------------ method called to for each event  ------------

void TtJet::JetTreeFeeder(Handle<std::vector<pat::Jet> > patJet, NJet* jtree, int eventId ) {

   for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++)
   {
       float emE0 = (*j1).emEnergyInEB()  + (*j1).emEnergyInEE()  + (*j1).emEnergyInHF() ;
       float hdE0 = (*j1).hadEnergyInHB() + (*j1).hadEnergyInHE() + (*j1).hadEnergyInHF()+ 
                    (*j1).hadEnergyInHO();

       jtree->FillBpatJ( eventId, j1->eta(), j1->phi(), emE0, hdE0, j1->p(), j1->pt() );

   }
}
void TtJet::jetAnalysis(Handle<std::vector<pat::Jet> > patJet, const edm::Event& iEvent, HTOP1* histo1){

   /// 1) overall jet information
   int nJets = patJet->size();
   int nj[3] = {0};
   for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++)
   {
       // Uncorrect Information
       // uncorrected pt from reco jet 
       reco::Jet recj = j1->recJet();
       double uncorrPt = recj.pt() ;
       // raw energy deposite in ECal and HCal
       float emE0 = (*j1).emEnergyInEB() + (*j1).emEnergyInEE() + (*j1).emEnergyInHF() ;
       float hdE0 = (*j1).hadEnergyInHB() +(*j1).hadEnergyInHE() + (*j1).hadEnergyInHF() + (*j1).hadEnergyInHO();
       double totalE = emE0 +hdE0 ;
       float emF0 = emE0 /totalE ;

       // pat information -> suppose to be corrected 
       double corrPt = j1->pt();
       float emF  = (*j1).emEnergyFraction() ;
       double EovH = EoverH(*j1);

       int nCon = (*j1).nConstituents();
       double rCon = NofJetConstituents( *j1 ) ;
       edm::RefVector<reco::TrackCollection>  assTk = (*j1).associatedTracks() ;

       if (nCon == 0) continue;

       histo1->Fill1a( j1->eta(), corrPt, (*j1).et(), EovH, nCon,(*j1).n60()/nCon, (*j1).n90()/nCon, uncorrPt, rCon, (*j1).towersArea(), assTk.size(), emF, emF0 );

       if ( (*j1).pt() > 20.0  &&  assTk.size() > 0) {
          histo1->Fill1h( EovH, (*j1).eta(), (*j1).towersArea(), nCon, assTk.size() ,rCon, (*j1).n60()/nCon );
       }

       if ( (*j1).pt() > 20. && fabs((*j1).eta()) < 2.0) nj[0]++;
       if ( (*j1).pt() > 30. && fabs((*j1).eta()) < 2.0) nj[1]++;
       if ( (*j1).pt() > 40. && fabs((*j1).eta()) < 2.0) nj[2]++;

       /*
       cout <<"   NoCorr= "<<j1->correctionFactor( pat::Jet::NoCorrection ) <<endl;
       cout <<"  defCorr= "<<j1->correctionFactor( pat::Jet::DefaultCorrection ) <<endl;
       cout <<"  udsCorr= "<<j1->correctionFactor( pat::Jet::udsCorrection ) <<endl;
       cout <<"    bCorr= "<<j1->correctionFactor( pat::Jet::bCorrection ) <<endl;
       cout <<"    cCorr= "<<j1->correctionFactor( pat::Jet::cCorrection ) <<endl;
       cout <<"    gCorr= "<<j1->correctionFactor( pat::Jet::gCorrection ) <<endl;
       cout <<"   NrCorr= "<<j1->correctionFactor( pat::Jet::NrOfCorrections ) <<endl;
       cout <<" "<<endl;
       */

       // crapy jet scope
       if ((*j1).towersArea()/(*j1).pt() > 0.005 ) {
	  histo1->Fill1i( (*j1).eta(), EovH, nCon, assTk.size(), (*j1).n60()/nCon, (*j1).n90()/nCon, rCon );
       }
       // underlying events setup check .. 
       if ( fabs(j1->eta()) < 1 ) {
          histo1->Fill1c( corrPt, uncorrPt, j1->energy(), totalE);
       }
       if ( fabs(j1->eta()) > 2.8  ) {
          histo1->Fill1c2( corrPt, uncorrPt, j1->energy(), totalE);
       }
       if ( fabs(j1->eta()) > 2.0 &&  fabs(j1->eta()) < 2.8   ) {
          histo1->Fill1c3( corrPt, uncorrPt, j1->energy(), totalE);
       }

   }
   histo1->Fill1f( nJets, nj[0], nj[1], nj[2] );

}

void TtJet::matchedWJetsAnalysis( std::vector<jmatch> mcwjets , std::vector<const reco::Candidate*> isoMuons ,
                                  HTOP8* histo8 ) {

  // look at matched jet properties
  for (std::vector<jmatch>::const_iterator j1 = mcwjets.begin(); j1 != mcwjets.end(); j1++) {

      pat::Jet truth =  (*j1).trueJet ;
      reco::Particle jmom = j1->mom;

      int    nCon = truth.nConstituents();
      //if (nCon == 0) continue; 
      double EovH = EoverH( truth  );
      double rCon = NofJetConstituents( truth ) ;
      double res_Pt  = (truth.pt() - jmom.pt()) / jmom.pt() ;
      double res_Eta = (truth.eta() - jmom.eta()) / jmom.eta();

      edm::RefVector<reco::TrackCollection>  assTk = truth.associatedTracks() ;
      histo8->Fill8c( truth.pt(), truth.eta(), EovH, nCon, truth.n60()/nCon, truth.n90()/nCon, rCon , truth.towersArea(), assTk.size(), truth.emEnergyFraction(), res_Pt, res_Eta );

     // dR(closest isoMu, matched bjet )
     double dR_WjMu = 999.;
     for (std::vector<const reco::Candidate*>::const_iterator m1= isoMuons.begin(); m1 != isoMuons.end(); m1++) {
         double dh = truth.eta() - (*m1)->eta();
	 double df = truth.phi() - (*m1)->phi();
	 double dR = sqrt( (dh*dh) + (df*df) );
	 if (dR < dR_WjMu ) {  dR_WjMu = dR; }
     }
     histo8->Fill8d( dR_WjMu );

  }

  // access the 4 momentum of parton 
  /* 
  LorentzVector qm[4] ;
  std::vector<pat::Jet> Wjj1;
  std::vector<pat::Jet> Wjj2;
  std::vector<reco::Particle> Wqq1;
  std::vector<reco::Particle> Wqq2;
  for (std::vector<jmatch>::const_iterator j1 = mcwjets.begin(); j1 != mcwjets.end(); j1++) {

      if ( (*j1).MomIdx == 1 ) {
         qm[0] = (*j1).sumP4 ;
         Wjj1.push_back(j1->trueJet);
         Wqq1.push_back(j1->mom);
       }
      if ( (*j1).MomIdx == 2 ) {
         qm[1] = (*j1).sumP4 ;
         Wjj1.push_back(j1->trueJet);
         Wqq1.push_back(j1->mom);
      }
      if ( (*j1).MomIdx == 3 ) {
         qm[2] = (*j1).sumP4 ;
         Wjj2.push_back(j1->trueJet);
         Wqq2.push_back(j1->mom);
      }
      if ( (*j1).MomIdx == 4 ) { 
         qm[3] = (*j1).sumP4 ;
         Wjj2.push_back(j1->trueJet);
         Wqq2.push_back(j1->mom);
      }
  }

  // look at the W mass from matching jets 
  if (qm[0].E() != 0 && qm[1].E() != 0 ) {
     LorentzVector pW = findW( qm[0], qm[1]);
     double momW = sqrt( pW.Px()*pW.Px() + pW.Py()*pW.Py() + pW.Pz()*pW.Pz() );
     double massW = sqrt ( pW.E()*pW.E() - pW.Px()*pW.Px() - pW.Py()*pW.Py() - pW.Pz()*pW.Pz() );
     double dh = Wjj1[0].eta() - Wjj1[1].eta() ;
     double df = Wjj1[0].phi() - Wjj1[1].phi() ;
     double dhq = Wqq1[0].eta() - Wqq1[1].eta() ;
     double dfq = Wqq1[0].phi() - Wqq1[1].phi() ;
     double dRWjj = sqrt( pow(dh,2) + pow(df,2) );
     double dRWqq = sqrt( pow(dhq,2) + pow(dfq,2) );
     double Res_dR = (dRWjj - dRWqq) / dRWqq ;

     histo6->Fill6b( massW, momW, dRWjj, Res_dR );
  }
  if (qm[2].E() != 0 && qm[3].E() != 0 ) {
     LorentzVector pW = findW( qm[2], qm[3]);
     double momW = sqrt( pW.Px()*pW.Px() + pW.Py()*pW.Py() + pW.Pz()*pW.Pz() );
     double massW = sqrt ( pW.E()*pW.E() - pW.Px()*pW.Px() - pW.Py()*pW.Py() - pW.Pz()*pW.Pz() );
     double dh = Wjj2[0].eta() - Wjj2[1].eta() ;
     double df = Wjj2[0].phi() - Wjj2[1].phi() ;
     double dhq = Wqq2[0].eta() - Wqq2[1].eta() ;
     double dfq = Wqq2[0].phi() - Wqq2[1].phi() ;
     double dRWjj = sqrt( pow(dh,2) + pow(df,2) );
     double dRWqq = sqrt( pow(dhq,2) + pow(dfq,2) );
     double Res_dR = (dRWjj - dRWqq) / dRWqq ;

     histo6->Fill6b( massW, momW, dRWjj, Res_dR );
  }
  */
}

void TtJet::matchedbJetsAnalysis( std::vector<jmatch> mcbjets, std::vector<jmatch> mcwjets,
                                  std::vector<const reco::Candidate*> isoMuons, HTOP7* histo7 ) {


  for (std::vector<jmatch>::const_iterator j1 = mcbjets.begin(); j1 != mcbjets.end(); j1++) {

      pat::Jet truth =  j1->trueJet ;
      reco::Particle jmom = j1->mom;
      
      double Res_Pt = (truth.pt() - jmom.pt()) / jmom.pt() ; 
      double fk1 = truth.bDiscriminator("jetProbabilityBJetTags") ;
      double fk2 = truth.bDiscriminator("trackCountingHighEffBJetTags") ;
      histo7->Fill7a( Res_Pt, fk1 , fk2 );
      //histo7->Fill7a( Res_Pt, truth.bDiscriminator("trackCountingJetTags") );

      // dR(closest isoMu, matched bjet )
      double dR_bmu = 999.;
      for (std::vector<const reco::Candidate*>::const_iterator m1= isoMuons.begin(); m1 != isoMuons.end(); m1++) {
          double dh = truth.eta() - (*m1)->eta();
          double df = truth.phi() - (*m1)->phi();
          double dR = sqrt( dh*dh + df*df );
          if (dR < dR_bmu ) {
             dR_bmu = dR ;
          }               
      } 
      double dR_bwj = 999.;
      for (std::vector<jmatch>::const_iterator j2 = mcwjets.begin(); j2 != mcwjets.end(); j2++ ) {
          double dh = truth.eta() - (j2->trueJet).eta();
          double df = truth.phi() - (j2->trueJet).phi();
          double dR = sqrt( dh*dh + df*df );
          if (dR < dR_bwj ) {
             dR_bwj = dR ;
          }               
      }
      if ( dR_bmu != 999. ) {  histo7->Fill7b( dR_bmu );  }
      if ( dR_bwj != 999. ) {  histo7->Fill7c( dR_bwj );  }

  }


}

void TtJet::bTagAnalysis(Handle<std::vector<pat::Jet> > Jets, HTOP7* histo7  ) {

  for (std::vector<pat::Jet>::const_iterator j1 = Jets->begin(); j1 != Jets->end(); j1++) {
      double fk1 = j1->bDiscriminator("jetProbabilityBJetTags") ;
      double fk2 = j1->bDiscriminator("trackCountingHighEffBJetTags") ;
      //cout<<" bDis(jetProb)  "<< fk1 <<endl;
      //cout<<" bDis(trkCount2)"<< fk2 <<endl;
      histo7->Fill7d( fk1, fk2 );
      /*
      cout<<" bDis(trkCount3)"<<j1->bDiscriminator("trackCountingHighPurBJetTag")<<endl;
      cout<<" bDis(softMu)   "<<j1->bDiscriminator("softMuonBJetTags")<<endl;
      cout<<" bDis(softEle)  "<<j1->bDiscriminator("softElectronBJetTags")<<endl;
      cout<<" bDis(2ndVtx)   "<<j1->bDiscriminator("impleSecondaryVertexBJetTags")<<endl;
      cout<<" bDis(cmbVtx)   "<<j1->bDiscriminator("combinedSVBJetTags")<<endl;
      */
  }

}

void TtJet::selectedWJetsAnalysis(Handle<std::vector<pat::Jet> > patJet, HTOP8* histo8  ) {

   std::vector<pat::Jet> wjets = evtSelected->WJetSelection( patJet );

   // jets selection
   double jPt[2]  = {0.0, 0.0};
   double jPhi[2] = {-9.0, -9.0};
   double jEta[2] = {-9.0, -9.0};

   int idx=0;
   for (std::vector<pat::Jet>::const_iterator j1 = wjets.begin(); j1 != wjets.end(); j1++)
   {
       if ( idx < 2){
          jPt[idx] = j1->pt() ;
          jPhi[idx] = j1->phi() ;
          jEta[idx] = j1->eta() ;
       }
       idx++;
   }
   histo8->Fill8a( jPt[0], jPt[1], jEta[0], jEta[1], wjets.size() );

   // rebuild W from selected jets!
   int jj1 =0;
   for (std::vector<pat::Jet>::const_iterator j1 = wjets.begin(); j1 != wjets.end(); j1++)
   {
       // only semi-lep
       //if ( topo != 1 || !pass ) continue;

       jj1++ ;
       int jj2 =0;

       double EovH1 = EoverH(*j1) ;
       if ((*j1).nConstituents() < 5 || EovH1 > 20 || EovH1 < 0.01) continue;
       if ((*j1).towersArea()/(*j1).pt() > 0.005 ) continue;

       for (std::vector<pat::Jet>::const_iterator j2 = wjets.begin(); j2 != wjets.end(); j2++)
       {
           jj2++;
           if (j2 <= j1) continue;

           double EovH2 = EoverH(*j2) ;

           if ((*j2).nConstituents() < 5 || EovH2 > 20 || EovH2 < 0.01) continue;
           if ((*j2).towersArea()/(*j2).pt() > 0.005 ) continue;
           //cout<<" flavor1: "<< (*j1).partonFlavour()<<"  vtx:"<<(*j1).vertex() <<endl;
           //cout<<" flavor2: "<< (*j2).partonFlavour()<<"  vtx:"<<(*j2).vertex() <<endl;
           //LorentzVector qm1 = (*j1).energyFractionHadronic() * (*j1).p4();
           //LorentzVector qm2 = (*j2).energyFractionHadronic() * (*j2).p4();
           LorentzVector qm1 = (*j1).p4();
           LorentzVector qm2 = (*j2).p4();
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

    if ( assTk.size() > 0) {
       double nC = static_cast<double>( theJet.nConstituents() );
       double nT = static_cast<double>( assTk.size() ) ;
       //cout <<" ratio = "<< nC/nT <<endl;
       return nC/nT ;
    } else {
       return -1 ;
    }
}

double TtJet::EoverH( pat::Jet theJet ) {

       float emE = theJet.emEnergyFraction() ;
       float hdE = theJet.energyFractionHadronic() ;
       double EovH = -2.0 ;
       if (hdE == 0.0 || (emE/hdE) >= 100.0 ) {
          EovH = 100.0;
       }
       else {
          EovH = emE/hdE ;
       }
       return EovH;
}

LorentzVector TtJet::findW(LorentzVector qm1, LorentzVector qm2 ) {

     //exclude = false;
     double m1 = qm1.P();
     double m2 = qm2.P();
     /*
     double m1m2 = qm1.Px()*qm2.Px() + qm1.Py()*qm2.Py() + qm1.Pz()*qm2.Pz() ;
     double massW = sqrt( (2.0*m1*m2) - (2.0*m1m2) ) ;
     double angleqq = acos( m1m2/(m1*m2) ) ;
     cout<<" q1("<<m1<<") q2("<<m2<<") dP= "<< fabs(m1-m2) <<",open angle of qq = "<< angleqq << endl;
     cout<<" mass of W = "<< massW <<endl;
     if ( angleqq > 2.2 && (m1+m2)> 100. ) {
        cout<<" Candidate 1 !!! "<<endl;
     }
     else if ( angleqq < 1.7 && (m1-m2)< 20. && (m1+m2) > 100.) {
        cout<<" Candidate 2 !!! "<<endl;
     }
     else if ( (m1+m2) > 100. && (m1+m2)< 1200. ) {
        cout<<" Candidate 3 !!! "<<endl;
     }
     else {
        exclude = true;
        if (  massW < 100. && massW > 60. ) {
          cout<<" NO GOOD @$#%@ " <<endl;
        } else {
          cout<<" *** exclude **** "<<endl;
        }
     }
     */
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
void TtJet::JetMatchedMuon( Handle<std::vector<pat::Jet> > patJet , Handle<std::vector<pat::Muon> > patMuon, 
                            const Event& iEvent, const EventSetup& iSetup, HTOP3* histo3){

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
                 histo3->Fill3d( (*jcaloPtr).energy(), m1->p(), j1->energy() );
              } else {
                 histo3->Fill3e( (*jcaloPtr).energy(), m1->p(), j1->energy() );
              }
           }
       }

  }

}

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
       FreeTrajectoryState ftsStart = getFTS(gp, gm, m1->charge(), covT, &*field);

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
                 //cout<<"         GP:"<<caloGP<<"  R="<<R<<" h:"<<getEta(caloGP.x(), caloGP.y(), caloGP.z())<<endl;
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
                 //cout<<" h:"<<getEta(fGP.x(), fGP.y(), fGP.z()) << endl;
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
              histo3->Fill3d( (*jcaloRef).energy(), m1->p(), j1->energy() );
           } else {
              histo3->Fill3e( (*jcaloRef).energy(), m1->p(), j1->energy() );
           }

       }
   }
   Done = true ;
}

FreeTrajectoryState TtJet::getFTS(GlobalPoint GP, GlobalVector GV, int charge, 
                                  const AlgebraicSymMatrix66& cov, const MagneticField* field){


  GlobalTrajectoryParameters tPars(GP, GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

double TtJet::getEta(double vx, double vy, double vz ) {

      double va = sqrt( vx*vx + vy*vy + vz*vz );

      double theta = acos( vz/va );
      double eta = (-1.0)*log( tan(theta/2.0) )  ;
      return eta;
}

