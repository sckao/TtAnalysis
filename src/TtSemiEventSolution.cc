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
//
//


// system include files
#include <memory>

// user include files
#include "TtSemiEventSolution.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//
//static bool ptDecrease(const iReco t1, const iReco t2) { return ( t1.p4.Pt() > t2.p4.Pt() ); }
static bool probDecrease(const iReco t1, const iReco t2) { return ( t1.prob > t2.prob ); }

// constructors and destructor
using namespace edm;
using namespace std;
TtSemiEventSolution::TtSemiEventSolution(const edm::ParameterSet& iConfig ): evt_Id(0), ini_Id(0), ent_Id(0), ent_sz(0) {
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool> ("debug");
  btag              = iConfig.getUntrackedParameter<bool> ("btag");
  nbtagged          = iConfig.getUntrackedParameter<int> ("nbtagged");
  JEScale           = iConfig.getUntrackedParameter<double> ("JEScale");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  tcmetSrc          = iConfig.getParameter<edm::InputTag> ("tcMetSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  jptSrc            = iConfig.getParameter<edm::InputTag> ("jptSource");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  algo              = iConfig.getUntrackedParameter<string> ("recoAlgo");
  
  evtSelected = new TtEvtSelector( iConfig );
  MCMatching  = new TtMCMatching();
  ttMuon      = new TtMuon();
  ttEle       = new TtElectron();
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );
  ttEff       = new TtEfficiency();
  tools       = new TtTools();
 
}


TtSemiEventSolution::~TtSemiEventSolution()
{
 
   delete evtSelected;
   delete MCMatching;
   delete ttMuon;
   delete ttEle;
   delete ttMET;
   delete ttJet;
   delete ttEff; 
   delete tools;
}

//
// member functions
//

// ------------ method called to for each event  ------------
void TtSemiEventSolution::BuildSemiTt( const edm::Event& iEvent, int topo, int evtId, tHisto* histos ) {

   // ********************************************
   // *     Build the tt events from DATA!!!     *
   // ********************************************

   // 0. Event Selection Setup
   double JetEtCut = 30. ;

   // 1. Reset object containers
   evt_Id = evtId ;
   isoLep.clear();
   selectedWJets.clear();
   selectedbJets.clear();
   selectedJets.clear();
   solvedMetP4.clear();
   semiTt.clear();
   KProbability.clear();
   
   // 2. Event Selection
   std::vector<double> bDisList;
   std::vector<bool> bTags;
   pass = evtSelected->eventSelection( topo, JetEtCut, isoLep, selectedJets, solvedMetP4, iEvent, "patMet", "patJet", &bTags, &bDisList );
   cout<<" event selected w/ jet "<<selectedJets.size() <<endl;
   // 3. Reconstruct event
   if ( btag ) {
      int Nbtags = 0 ;
      for (size_t b=0; b< bTags.size(); b++) {
          if ( bTags[b] ) Nbtags ++ ;
      }
      cout <<"  bTager Sizer = "<< bTags.size() <<" w/ "<< Nbtags  <<endl;

      if ( pass > 3 && Nbtags >= 2 && nbtagged == 2 ) semiTt = recoSemiLeptonicTtEvent( selectedJets, isoLep, solvedMetP4, &bTags, &KProbability, histos );
      if ( pass > 3 && Nbtags == 1 && nbtagged == 1 ) semiTt = recoSemiLeptonicTtEvent( selectedJets, isoLep, solvedMetP4, &bTags, &KProbability, histos );
      if ( pass > 3 && nbtagged == 0 ) semiTt = recoSemiLeptonicTtEvent(selectedJets, isoLep, solvedMetP4, &bTags, &KProbability, histos );

   } else {
      if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent( selectedJets, isoLep, solvedMetP4, NULL, &KProbability, histos );
   }

   if ( semiTt.size() >= 2 )     KeepBuildInfo( true );

   cout<<" ------ reco finished --------------"<<endl;

}

// ntuple version
void TtSemiEventSolution::BuildSemiTt( const edm::Event& iEvent, int topo, int evtId, tNtuple* ntuples ) {

   // ********************************************
   // *     Build the tt events from DATA!!!     *
   // ********************************************

   // 0. Event Selection Setup
   double JetEtCut = 30. ;

   // 1. Reset object containers
   evt_Id = evtId ;
   isoLep.clear();
   selectedWJets.clear();
   selectedbJets.clear();
   selectedJets.clear();
   solvedMetP4.clear();
   semiTt.clear();
   KProbability.clear();
   
   // 2. Event Selection
   std::vector<double> bDisList;
   std::vector<bool> bTags;
   pass = evtSelected->eventSelection( topo, JetEtCut, isoLep, selectedJets, solvedMetP4, iEvent, "patMet", "patJet", &bTags, &bDisList );
   cout<<" event selected w/ jet "<<selectedJets.size() <<endl;
   // 3. Reconstruct event
   if ( btag ) {
      int Nbtags = 0 ;
      for (size_t b=0; b< bTags.size(); b++) {
          if ( bTags[b] ) Nbtags ++ ;
      }
      cout <<"  bTager Sizer = "<< bTags.size() <<" w/ "<< Nbtags  <<endl;

      if ( pass > 3 && Nbtags >= 2 && nbtagged == 2 ) semiTt = recoSemiLeptonicTtEvent( selectedJets, isoLep, solvedMetP4, &bTags, &KProbability, ntuples );
      if ( pass > 3 && Nbtags == 1 && nbtagged == 1 ) semiTt = recoSemiLeptonicTtEvent( selectedJets, isoLep, solvedMetP4, &bTags, &KProbability, ntuples );
      if ( pass > 3 && nbtagged == 0 ) semiTt = recoSemiLeptonicTtEvent(selectedJets, isoLep, solvedMetP4, &bTags, &KProbability, ntuples );

   } else {
      if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent( selectedJets, isoLep, solvedMetP4, NULL, &KProbability, ntuples );
   }

   // store information in ntuples
   if ( ntuples != NULL && pass == 4) {
      evt_cout++;
      for (size_t i=0; i< selectedJets.size(); i++) {
          double b_dis = ( bDisList.size() == selectedJets.size() ) ? bDisList[i] : -1 ;
          (*ntuples).jetTree->FillB( evtId, i, b_dis, selectedJets[i]->px(), selectedJets[i]->py(), selectedJets[i]->pz(), 
                                                      selectedJets[i]->energy(), selectedJets[i]->pt() );
      }

      for (size_t i=0; i< isoLep.size(); i++)  {
          (*ntuples).muTree->FillB( evtId, i, -1, isoLep[i]->px(), isoLep[i]->py(), isoLep[i]->pz(), isoLep[i]->energy(), isoLep[i]->pt() );
      }

      for (size_t i=0; i< solvedMetP4.size(); i++) {
          (*ntuples).neuTree->FillB( evtId, i, -1, solvedMetP4[i].px(), solvedMetP4[i].py(), solvedMetP4[i].pz(), 
                                                   solvedMetP4[i].energy(), solvedMetP4[i].pt() );
      }
   }
   cout<<" ------ reco finished --------------"<<endl;

}


// topoIdx : 4:tau+jets  3:electron+jest  2:di-lep,  1:muon+jets,  0:hadron,  -1:Non-Tt event
void TtSemiEventSolution::MCBuildSemiTt( const edm::Event& iEvent, int topoIdx, int evtId, tHisto* histos ) {

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   // check the topology
   //int topo = evtSelected->MCEvtSelection(genParticles);
   //bool build = ( topo ==  topoIdx ) ? true : false ;

   // 1. Reset object containers
   mcLep.clear();
   mcWJets.clear();
   mcbJets.clear();
   semiMCTt.clear();

   if ( pass < 4 ) return;

   if ( topoIdx == 1 )  mcLep = MCMatching->matchMuon(genParticles, isoLep, (*histos).hMuon );  
   if ( topoIdx == 3 )  mcLep = MCMatching->matchElectron(genParticles, isoLep, (*histos).hEle );  
   std::vector<jmatch> mc_jets = MCMatching->matchJets(genParticles, selectedJets, (*histos).hBJet, (*histos).hWJet );  

   // 2. classify wjets and bjets ..... parton Btagging
   int wQ = 0 ;
   std::vector<int> wjid;
   int bjh = -1;
   int bjl = -1;
   std::vector<iReco> hadWs;
   std::vector<const reco::Candidate*> bn;
   std::vector<const reco::Candidate*> bp;
   std::vector<bool> bTags(2, true) ;
   bool findhadW = false ;
   for (size_t i =0; i< mc_jets.size(); i++ ){

       if ( abs( mc_jets[i].MomIdx ) < 5  ) {
          mcWJets.push_back( mc_jets[i].trueJet ) ;
          wjid.push_back( mc_jets[i].Idx );
          wQ +=  mc_jets[i].MomIdx ;
          // reco hadronic W
          if ( mcWJets.size() == 2 ) findhadW = recoW( mcWJets, hadWs, (*histos).hWs[2] );
       }
       if ( mc_jets[i].MomIdx ==  5 ) bn.push_back( mc_jets[i].trueJet );
       if ( mc_jets[i].MomIdx == -5 ) bp.push_back( mc_jets[i].trueJet );
       if ( mcLep.size() == 0 ) continue;
       if ( mcLep[0]->charge() == -1 && mc_jets[i].MomIdx == -5 ) bjl = mc_jets[i].Idx ; 
       if ( mcLep[0]->charge() == -1 && mc_jets[i].MomIdx ==  5 ) bjh = mc_jets[i].Idx ; 
       if ( mcLep[0]->charge() ==  1 && mc_jets[i].MomIdx ==  5 ) bjl = mc_jets[i].Idx ; 
       if ( mcLep[0]->charge() ==  1 && mc_jets[i].MomIdx == -5 ) bjh = mc_jets[i].Idx ; 
   }
   // 3. reco leptonic W

   std::vector<iReco> wSols;
   std::vector<iReco> lepWs;
   bool findlepW = recoW( mcLep, solvedMetP4[0], wSols, (*histos).hWs[0] );
   // find the right solution

   if ( wSols.size() > 1) {
      int wl = MCMatching->matchLeptonicW( genParticles, wSols );
      if ( wl > -1 ) lepWs.push_back( wSols[wl] );
   }
   if ( wSols.size() == 1 )  lepWs.push_back( wSols[0] );

   // 4. reco leptonic Top

   std::vector<iReco> lepTs;
   bool findlepT = false;
   if ( findlepW ) {
      std::vector<const reco::Candidate*> b_temp;
      if ( mcLep[0]->charge() == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      if ( mcLep[0]->charge() ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      findlepT = recoTop( lepWs , b_temp , lepTs, &bTags, (*histos).hTops[0], true );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findlepT )   semiMCTt.push_back( lepTs[0] ) ;
   }
   // 5. reco hadronic Top

   std::vector<iReco> hadTs;
   bool findhadT = false ; 
   if ( findhadW ) {
      std::vector<const reco::Candidate*> b_temp;
      if ( wQ ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      if ( wQ == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      findhadT = recoTop( hadWs , b_temp , hadTs, &bTags, (*histos).hTops[2], true );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findhadT ) semiMCTt.push_back( hadTs[0] );
   }

   // 6. look at top mass distribution and results
   if ( semiMCTt.size() == 2 ) {

      Idx tws(4,0);
      std::vector<Idx> twIdx( 1 , tws ) ;
      std::vector<bool> usedHadT( hadTs.size() ,false);
      std::vector<bool> usedLepT( lepTs.size() ,false); 
      ResultRecord( 0, 0, twIdx, lepTs, usedLepT, hadTs, usedHadT, lepWs, hadWs, NULL, histos ) ;

      KeepBuildInfo( false );
   }

   bool matchedpass =  ( pass >= 4  ) ? true : false ;
   McRecoCompare( 1, 0, matchedpass , histos );
   cout<<" ====  MCMatching : size of lepT:"<< lepTs.size() <<" , hadT:"<<hadTs.size()  <<endl;

}

void TtSemiEventSolution::MCBuildSemiTt( const edm::Event& iEvent, int topoIdx, int evtId, tNtuple* ntuples ) {

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   // check the topology
   //int topo = evtSelected->MCEvtSelection(genParticles);
   //bool build = ( topo ==  topoIdx ) ? true : false ;

   // 1. Reset object containers
   mcLep.clear();
   mcWJets.clear();
   mcbJets.clear();
   semiMCTt.clear();

   if ( pass < 4 ) return;

   if ( topoIdx == 1 )  mcLep = MCMatching->matchMuon(genParticles, isoLep );  
   if ( topoIdx == 3 )  mcLep = MCMatching->matchElectron(genParticles, isoLep );  
   std::vector<jmatch> mc_jets = MCMatching->matchJets(genParticles, selectedJets );  

   // 2. classify wjets and bjets ..... parton Btagging
   int wQ = 0 ;
   std::vector<int> wjid;
   int bjh = -1;
   int bjl = -1;
   std::vector<iReco> hadWs;
   std::vector<const reco::Candidate*> bn;
   std::vector<const reco::Candidate*> bp;
   std::vector<bool> bTags(2, true) ;
   bool findhadW = false ;
   for (size_t i =0; i< mc_jets.size(); i++ ){

       if ( abs( mc_jets[i].MomIdx ) < 5  ) {
          mcWJets.push_back( mc_jets[i].trueJet ) ;
          wjid.push_back( mc_jets[i].Idx );
          wQ +=  mc_jets[i].MomIdx ;
          // reco hadronic W
          if ( mcWJets.size() == 2 ) findhadW = recoW( mcWJets, hadWs );
       }
       if ( mc_jets[i].MomIdx ==  5 ) bn.push_back( mc_jets[i].trueJet );
       if ( mc_jets[i].MomIdx == -5 ) bp.push_back( mc_jets[i].trueJet );
       if ( mcLep.size() == 0 ) continue;
       if ( mcLep[0]->charge() == -1 && mc_jets[i].MomIdx == -5 ) bjl = mc_jets[i].Idx ; 
       if ( mcLep[0]->charge() == -1 && mc_jets[i].MomIdx ==  5 ) bjh = mc_jets[i].Idx ; 
       if ( mcLep[0]->charge() ==  1 && mc_jets[i].MomIdx ==  5 ) bjl = mc_jets[i].Idx ; 
       if ( mcLep[0]->charge() ==  1 && mc_jets[i].MomIdx == -5 ) bjh = mc_jets[i].Idx ; 
   }

   // 3. reco leptonic W
   std::vector<iReco> wSols;
   std::vector<iReco> lepWs;
   bool findlepW = recoW( mcLep, solvedMetP4[0], wSols );
   // find the right solution

   int wl  = -1 ;
   if ( wSols.size() > 1) {
      wl = MCMatching->matchLeptonicW( genParticles, wSols );
      if ( wl > -1 )  lepWs.push_back( wSols[wl] );
   }
   if ( wSols.size() == 1 ) {
      lepWs.push_back( wSols[0] );
      wl = 0 ;
   }

   // 4. reco leptonic Top

   std::vector<iReco> lepTs;
   bool findlepT = false;
   if ( findlepW ) {
      std::vector<const reco::Candidate*> b_temp;
      if ( mcLep[0]->charge() == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      if ( mcLep[0]->charge() ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      findlepT = recoTop( lepWs , b_temp , lepTs, &bTags, NULL, true );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findlepT )   semiMCTt.push_back( lepTs[0] ) ;
   }
   // 5. reco hadronic Top

   std::vector<iReco> hadTs;
   bool findhadT = false ; 
   if ( findhadW ) {
      std::vector<const reco::Candidate*> b_temp;
      if ( wQ ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      if ( wQ == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      findhadT = recoTop( hadWs , b_temp , hadTs, &bTags, NULL, true );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findhadT ) semiMCTt.push_back( hadTs[0] );
   }

   // 6. look at top mass distribution and results
   if ( semiMCTt.size() == 2 &&  ntuples != NULL && pass == 4 ) {
      double lepW_mt = tools->getMt( lepWs[0].q4v[0].second, lepWs[0].q4v[1].second) ;
      (*ntuples).mcmTree->FillB( evtId, wjid[0], wjid[1], bjh, bjl, -1*wSols[wl].from.second, 
                                 ini_Id,  ent_sz, ent_Id, 
                                 pass,    1.,     1.,     1.,
                                 hadTs[0].p4.M(), lepTs[0].p4.M(),  hadWs[0].p4.M(), lepW_mt  );

   }
   cout<<" ====  MCMatching : size of lepT:"<< lepTs.size() <<" , hadT:"<<hadTs.size()  <<endl;

}


void TtSemiEventSolution::McRecoCompare( int topo, int r,  bool matchedpass, tHisto* histos ) {

   int sz = static_cast<int>( AllTt.size() );
   int m = sz-1 ;
   if ( sz > 1  && AllTt[m].isData == false  ) {

      int sz1 = static_cast<int>( AllTt[r].Tt.size() )/2 ;
      for (int i1=0; i1< sz1; i1++ ) {

          int i2 = i1*2  ;
          if ( sz > 1 && matchedpass ) MCTruthCheck( i2, AllTt[m].Tt, AllTt[r].Tt, AllTt[m].WJ, AllTt[r].Js, AllTt[m].bJ, AllTt[r].Js, (*histos).hMObj );

          accuracySemiTt( AllTt[m].Tt , AllTt[r].Tt , (*histos).hTop );
          if ( topo == 1 ) ttEff->IsoLeptonEfficiency( AllTt[r].Ls, AllTt[m].Ls, (*histos).hTop );
          if ( topo == 3 ) ttEff->IsoLeptonEfficiency( AllTt[r].Ls, AllTt[m].Ls, (*histos).hTop );
      }
   }
   // Reset "AllTt" containers
   AllTt.clear();

}

void TtSemiEventSolution::KeepBuildInfo( bool isData ) {

     if ( isData && semiTt.size() >= 2 ) {
        TtResult thisBuild;
        thisBuild.Tt = semiTt;
	thisBuild.Ls = isoLep;
	thisBuild.Js = selectedJets;
	thisBuild.bJ = selectedbJets;
	thisBuild.WJ = selectedWJets;
	thisBuild.metP4 = solvedMetP4;
        thisBuild.isData = true;
        AllTt.push_back( thisBuild );
     }
     if ( !isData && semiMCTt.size() == 2 ) {
        TtResult thisBuild;
        thisBuild.Tt = semiMCTt;
	thisBuild.Ls = mcLep;
	thisBuild.Js = selectedJets;
	thisBuild.bJ = mcbJets;
	thisBuild.WJ = mcWJets;
	thisBuild.metP4 = solvedMetP4;
        thisBuild.isData = false;
        AllTt.push_back( thisBuild );
     }

}
    
void TtSemiEventSolution::MCTruthCheck( int k, std::vector<iReco> mcTt, std::vector<iReco> rcTt, 
                     std::vector<const reco::Candidate*> mcWJets, std::vector<const reco::Candidate*> rcWJets,
                     std::vector<const reco::Candidate*> mcBJets, std::vector<const reco::Candidate*> rcBJets, HTOP6* histo6 ) {

   //1. mc matched reco W 
   std::vector<iReco> mcWs;
   bool mc = recoW( mcWJets, mcWs );
   //2. reco W
   std::vector<iReco> rcWs;
   bool rc = recoW( rcWJets, rcWs );

   if ( !mc || !rc ) cout<<" RECO ERROR !!!! "<<endl;

   // 3. from leptonic top
   double dRbj0 = tools->getdR( mcBJets[mcTt[0].from.second]->p4(), rcBJets[rcTt[k].from.second]->p4() );
   double dRw0  = tools->getdRy( mcWs[mcTt[0].from.first].p4, rcWs[rcTt[k].from.first].p4 );
   double dRt0  = tools->getdRy( mcTt[0].p4, rcTt[k].p4 );
   // 4. from hadronic top
   double dRbj1 = tools->getdR( mcBJets[mcTt[1].from.second]->p4(), rcBJets[rcTt[k+1].from.second]->p4() );
   double dRw1  = tools->getdRy( mcWs[mcTt[1].from.first].p4, rcWs[rcTt[k+1].from.first].p4 );
   double dRt1  = tools->getdRy( mcTt[1].p4, rcTt[k+1].p4 );
   
   histo6->Fill6a( dRt0, dRbj0, dRw0, dRt1, dRbj1, dRw1 );

   // 5. bjet and recoW from leptonic and hadronic top
   double dRbw_mch = tools->getdR( mcBJets[mcTt[1].from.second]->p4(), mcWs[mcTt[1].from.first].p4 );
   double dRbw_rch = tools->getdR( rcBJets[rcTt[k+1].from.second]->p4(), rcWs[rcTt[k+1].from.first].p4 );
   double dRbw_mcl = tools->getdR( mcBJets[mcTt[0].from.second]->p4(), mcWs[mcTt[0].from.first].p4 );
   double dRbw_rcl = tools->getdR( rcBJets[rcTt[k].from.second]->p4(), rcWs[rcTt[k].from.first].p4 );

   histo6->Fill6b( dRbw_mcl, dRbw_rcl, dRbw_mch, dRbw_rch );

   // 6. Record the matched reconstruction
   if (dRt1 < 0.01 && dRt0 < 0.01) {
       double mclepTMass = mcTt[0].p4.M() ; 
       double mchadTMass = mcTt[1].p4.M() ; 
       if ( pass  > 3 ) histo6->Fill6c(  mclepTMass , mchadTMass );
       if ( pass == 4 ) histo6->Fill6c0( mclepTMass , mchadTMass );
       if ( pass == 5 ) histo6->Fill6c1( mclepTMass , mchadTMass );
       if ( pass == 6 ) histo6->Fill6c2( mclepTMass , mchadTMass );
       if ( pass >= 7 ) histo6->Fill6c3( mclepTMass , mchadTMass );
   }

}

void TtSemiEventSolution::accuracySemiTt( std::vector<iReco> ttMC, std::vector<iReco> ttReco, HTOP9* histo9 ) {

     int nb = 0;
     int nw = 0;
     int nl = 0;
     if (ttMC.size() == 0 && ttReco.size() == 0 ) {
        histo9->Fill9j( 0, 0, 0) ;
        return ;
     }
     if (ttMC.size() == 0 && ttReco.size() > 0 ) {
        histo9->Fill9j( 0, 0, 0) ;
        return ;
     }
     if (ttMC.size() > 0 && ttReco.size() == 0 ) {
        histo9->Fill9j( 0, 0, 0) ;
        return ;
     }

     if ( ttMC[0].q4v[2] == ttReco[0].q4v[2] ) nb++ ; // b of leptonic top
     if ( ttMC[1].q4v[2] == ttReco[1].q4v[2] ) nb++ ; // b of hadronic top

     // w of hadronic top
     if ( ttMC[1].q4v[0] == ttReco[1].q4v[0] && ttMC[1].q4v[1] == ttReco[1].q4v[1] ) nw = 2 ;
     if ( ttMC[1].q4v[0] == ttReco[1].q4v[1] && ttMC[1].q4v[1] == ttReco[1].q4v[0] ) nw = 2 ;

     if ( ttMC[1].q4v[0] == ttReco[1].q4v[0] && ttMC[1].q4v[1] != ttReco[1].q4v[1] ) nw = 1 ;
     if ( ttMC[1].q4v[0] != ttReco[1].q4v[0] && ttMC[1].q4v[1] == ttReco[1].q4v[1] ) nw = 1 ;
     if ( ttMC[1].q4v[0] == ttReco[1].q4v[1] && ttMC[1].q4v[1] != ttReco[1].q4v[0] ) nw = 1 ;
     if ( ttMC[1].q4v[0] != ttReco[1].q4v[1] && ttMC[1].q4v[1] == ttReco[1].q4v[0] ) nw = 1 ;

     // leptons of hadronic top
     if ( ttMC[0].q4v[0] == ttReco[0].q4v[0] ) nl = 1;
   
     histo9->Fill9j( nb, nw, nl ) ;

}


// hadronic channel 
bool TtSemiEventSolution::recoW(  std::vector<const reco::Candidate*> wjets, std::vector<iReco>& wCandidate, HTOP10* histo, std::vector<bool>* btags  ){

    bool findW = false ; 
    if ( wjets.size() < 2 ) return findW ;

    double norm_p = 0;
    for (size_t i=0; i< wjets.size(); i++ ) {
       for (size_t j =i+1; j < wjets.size(); j++ ) {

           if ( btags != NULL ) {
              if ( (*btags)[i] || (*btags)[j] ) continue ;
           }

           LorentzVector m4W = wjets[i]->p4() + wjets[j]->p4() ;
	   if ( m4W.M2() < 0. ) continue;

           double chi = (m4W.M() - 80.4) / (80.4*0.1) ;
           double wProb = ( m4W.M() <= 160.8 ) ? exp( -0.5*chi*chi ) : 0. ;
           norm_p += wProb ;
          
           iReco candW;
           candW.p4   = m4W;
           candW.from = make_pair(i,j);
           candW.ptr  = make_pair( wjets[i], wjets[j] );
           candW.prob = exp( -0.5*chi*chi ) ;

	   iParton qi = make_pair( wjets[i]->pdgId() , wjets[i]->p4() );
	   iParton qj = make_pair( wjets[j]->pdgId() , wjets[j]->p4() );
	   candW.q4v.push_back( qi );
	   candW.q4v.push_back( qj );

           wCandidate.push_back( candW );

           findW = true;
       }
    }

    if ( findW ) {
       for (size_t i=0; i < wCandidate.size(); i++ ) {
           wCandidate[i].prob = (norm_p > 0.) ? wCandidate[i].prob/norm_p : 0. ;
           double wMass =  wCandidate[i].p4.M() ;
           if( wMass >= 160 ) wMass =  159.9 ;
           if ( histo != NULL ) histo->Fill10a( wMass );
        }
    }

    return findW;
}

// leptonic channel 
bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton,  LorentzVector metp4, 
                                 std::vector<iReco>& wCandidate, HTOP10* histo ){
  

    bool findW = false ; 
    if ( lepton.size() != 1 ) return findW ;

    // use w mass constrain to solve 4 momentum of W
    double xW  = lepton[0]->p4().Px() + metp4.Px() ;
    double yW  = lepton[0]->p4().Py() + metp4.Py() ;
    double nuPt2 = ( metp4.Px()*metp4.Px() ) + ( metp4.Py()*metp4.Py() ); 
 
    double zl = lepton[0]->p4().Pz() ;
    double El = lepton[0]->p4().E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    // EtW and PtW are only for mT of W calculation 
    double EtW  = lepton[0]->pt() + metp4.E() ;
    double PtW2 = (xW*xW) + (yW*yW);
    double MtW2 = (EtW*EtW) - PtW2 ;
    if ( MtW2 < 0. ) return findW ;

    if ( (B*B) < (4.*A*C) ) return findW ;

    // 2 solutions for  z momentum of neutrino
    double nz1 = (-1.*B + sqrt(B*B - 4.*A*C )) / (2.*A) ;
    double nz2 = (-1.*B - sqrt(B*B - 4.*A*C )) / (2.*A) ;
   
    // pick a better solution!
    for (int i=1; i<3; i++ ) {
       double nz = 0.0;  
       if (i==1) nz = nz1;
       if (i==2) nz = nz2;
       double ENu2 = ( metp4.Px()*metp4.Px() ) + ( metp4.Py()*metp4.Py() ) + (nz*nz);
       double zW = lepton[0]->p4().Pz() + nz ;
       double EW = lepton[0]->p4().E()  + sqrt(ENu2) ;
       LorentzVector np4 = LorentzVector( metp4.Px(), metp4.Py(), nz, sqrt(ENu2) );
       double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

       if (massTest > 0. ) {
          iReco candW; 
          LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
          candW.p4   = mW;
          candW.from = make_pair( -1, -1*i );      

          // probability must be 0.5 ...because 2 pz solution
          candW.prob = 0.5;

          iParton l1 = make_pair( lepton[0]->pdgId() , lepton[0]->p4() );
          int neuId =  ( lepton[0]->pdgId() > 0 ) ? lepton[0]->pdgId()+1 : lepton[0]->pdgId()-1 ;
          iParton l2 = make_pair( neuId , np4 );

          candW.q4v.push_back(l1);
          candW.q4v.push_back(l2);

          wCandidate.push_back( candW );
          findW = true;
       }
    }
    sort(wCandidate.begin(), wCandidate.end(), probDecrease );

    if ( findW ) {
       for (size_t i=0; i < wCandidate.size(); i++ ) {
           double wMt = ( wCandidate[i].p4.Mt() >= 160) ? 159.9 : wCandidate[i].p4.Mt() ;
           if ( histo != NULL ) histo->Fill10a( wMt );
       }
    }

    return findW;
}


// reco leptonic W with muon, MET agian ; assume pz of neutrino is zeor => MET is p4 of neutrino
bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton, LorentzVector metp4,
                                 std::vector<iReco>& wCandidate, bool FoundWSolution, HTOP10* histo ){
 
   if ( FoundWSolution ) return FoundWSolution;

   bool findW = false;
   if (  lepton.size() != 1 ) return FoundWSolution ;

   // Test the MT of lepton + MET
   double xW = lepton[0]->p4().Px() + metp4.Px() ;
   double yW = lepton[0]->p4().Py() + metp4.Py() ;
   double zW = lepton[0]->p4().Pz();

   double ptW = sqrt( (xW*xW) + (yW*yW) );
   double EtW = lepton[0]->pt() + metp4.E() ;
   double MtW2 = (EtW*EtW) - (ptW*ptW) ;
   if ( MtW2 < 0. ) return findW ;

   double pW2 =  (xW*xW) + (yW*yW) + (zW*zW) ;
   double EW2 =  lepton[0]->energy() + metp4.E() ;
   double masstest = EW2 - pW2 ;
   if ( masstest < 0 ) return findW ;

   // assume the neutrino has no Pz, only work for those exceed the Jacobian peak
   if ( MtW2 > 6464.16 ) {

      double EW = lepton[0]->p4().E() + metp4.E() ;
      double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

      if ( massTest > 0. ) {
	 iReco candW; 
         LorentzVector p4W = LorentzVector(xW,yW,zW,EW);
	 candW.p4   = p4W;
	 candW.from = make_pair(-1,-1);      

         candW.prob = 1. ;

         iParton l1 = make_pair( lepton[0]->pdgId() , lepton[0]->p4() );
         int neuId =  ( lepton[0]->pdgId() > 0 ) ? lepton[0]->pdgId()+1 : lepton[0]->pdgId()-1 ;
         iParton l2 = make_pair( neuId ,  metp4 );

	 candW.q4v.push_back(l1);
	 candW.q4v.push_back(l2);

	 wCandidate.push_back( candW );
	 findW = true;    
      } 
   }
   if ( findW ) {
      for (size_t i=0; i < wCandidate.size(); i++ ) {
          double wMt = ( wCandidate[i].p4.Mt() >= 160) ? 159.9 : wCandidate[i].p4.Mt() ;
          if ( histo != NULL ) histo->Fill10a( wMt );
      }
   }

   return findW;

}

bool TtSemiEventSolution::recoTop( std::vector<iReco> wCandidate, std::vector<const reco::Candidate*> bjets, std::vector<iReco>& tCandidates, std::vector<bool>* btags, HTOP11* histo, bool isMCMatching ) {

  bool findt = false ;
  tCandidates.clear();

  // check # of bjets
  int nBJet = 0;
  if ( btags != NULL ) {
     for (size_t i=0; i< btags->size(); i++ ) { 
         if ( (*btags)[i] ) nBJet++;
     }
  }

  double norm_p = 0.;
  std::vector<iReco> bWList;  
  for(size_t i=0; i < wCandidate.size(); i++ ){
     for (size_t j=0; j< bjets.size(); j++ ){

         // if btagging info is avaible and at least 2 jets are tagged , only apply the jets with btag
         if ( btags != NULL && nBJet  >= 2 &&  !(*btags)[j] && !isMCMatching ) continue;

         // avoid double use of bjets and w jets
         int bj = static_cast<int>(j);
         bool doubleUsed = (bj == wCandidate[i].from.first || bj == wCandidate[i].from.second) ? true : false ;
         if ( doubleUsed && !isMCMatching ) continue;

         LorentzVector m4T =  wCandidate[i].p4 + bjets[j]->p4() ;

	 if ( m4T.M2() <= 0. ) continue;
         iReco bWpair;
         bWpair.p4   = m4T ;
         bWpair.from = make_pair( i , j );
         double tProb = ( m4T.M() < 350. ) ? wCandidate[i].prob : 0. ;
         bWpair.prob = tProb ;
         norm_p += tProb ;

         iParton qb = make_pair( 5 , bjets[j]->p4() );
         bWpair.q4v.push_back( wCandidate[i].q4v[0] );
         bWpair.q4v.push_back( wCandidate[i].q4v[1] );
         bWpair.q4v.push_back( qb );
         HadronicTopCombinatoric( bWList, wCandidate, bWpair, wCandidate[i] );

	 findt = true;    
     }
  }

  sort(bWList.begin(), bWList.end(), probDecrease );
  tCandidates = bWList ;

  if ( findt ) { 
      for ( size_t i=0; i < tCandidates.size(); i++) {    
          tCandidates[i].prob = ( norm_p > 0. ) ? tCandidates[i].prob/norm_p :  0. ;
          double theMass = ( tCandidates[i].p4.M()  >= 480. ) ? 479.9 : tCandidates[i].p4.M() ;
	  double thePt   = ( tCandidates[i].p4.Pt() >= 480. ) ? 479.9 : tCandidates[i].p4.Pt() ;
          if ( histo != NULL ) histo->Fill11a( theMass, thePt );
          // only meaningful for no b-tagging    
	  if ( i == 0 && pass == 4 && histo != NULL ) histo->Fill11b( theMass, thePt );
      }
  }

  return findt;
}

// New general method for histogram, including b-tagging and no b-tagging scenario
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent( std::vector<const reco::Candidate*> theJets,
                      std::vector<const reco::Candidate*> theLep, std::vector<LorentzVector>& metp4,  
                      std::vector<bool>* btags, std::vector<iProb>* kProb, tHisto* histos ) {

   std::vector<iReco> TtCollection;
   std::vector<Idx> twIdx ;
   TtCollection.clear();
   
   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = recoW( theLep, metp4[0], lepWs, (*histos).hWs[1] );
   if ( leptonic ) metp4.clear()    ;
   for ( size_t i =0; i< lepWs.size(); i++ ) {
        metp4.push_back( lepWs[i].q4v[1].second ) ;
   }
   
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = recoW( theJets, hadWs, (*histos).hWs[3], btags );

   if ( leptonic && hadronic ) {
      std::vector<iReco> lepTs;
      std::vector<iReco> hadTs;
      bool findlepT  = recoTop( lepWs ,theJets, lepTs, btags, (*histos).hTops[1] );
      bool findhadT  = recoTop( hadWs ,theJets, hadTs, btags, (*histos).hTops[3] );

      if ( !findlepT || !findhadT ) return TtCollection ;
      if (algo == "ptMin" ) Algo_PtMin( lepTs, hadTs, hadWs, twIdx );
      if (algo == "beta"  ) Algo_Beta( lepTs, hadTs, lepWs, hadWs, twIdx );
      if (algo == "zero"  ) Algo_Zero( lepTs, hadTs, hadWs, twIdx, btags );
      if (algo == "kConstrain" ) Algo_KConstrain( lepTs, hadTs, hadWs, twIdx, kProb, btags );

      if ( theJets.size() >= 4 ) cout<<" lepT:"<< lepTs.size() <<" hadT:"<< hadTs.size() <<endl;

      std::vector<bool> usedHadT( hadTs.size() ,false);
      std::vector<bool> usedLepT( lepTs.size() ,false);
      double maxProb = -1;
      for ( size_t i=0; i< twIdx.size(); i++ ) {
          TtCollection.push_back( lepTs[ twIdx[i][0] ] );
	  TtCollection.push_back( hadTs[ twIdx[i][1] ] );
          ResultRecord( i, 1, twIdx, lepTs, usedLepT, hadTs, usedHadT, lepWs, hadWs, kProb, histos ) ;
          if ( algo == "kConstrain" ) maxProb = ( (*kProb)[i].Tt > maxProb ) ? (*kProb)[i].Tt : maxProb ;
      }
      if ( kProb != NULL && (*histos).hTop != NULL ) (*histos).hTop->Fill9a( maxProb );

      int nBJet = 0;
      if ( btags != NULL ) {
         for (size_t i=0; i< btags->size(); i++ ) { 
             if ( (*btags)[i] ) nBJet++;
         }
      }
      if ( (*histos).hTop != NULL ) (*histos).hTop->Fill9e( theJets.size(), nBJet );

   }

   return TtCollection ;

}

// New general method , for ntuple analysis 
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent( std::vector<const reco::Candidate*> theJets,
                      std::vector<const reco::Candidate*> theLep, std::vector<LorentzVector>& metp4,  
                      std::vector<bool>* btags, std::vector<iProb>* kProb, tNtuple* ntuples ) {

   std::vector<iReco> TtCollection;
   std::vector<Idx> twIdx ;
   TtCollection.clear();
   
   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = recoW( theLep, metp4[0], lepWs );
   if ( leptonic ) metp4.clear()    ;
   for ( size_t i =0; i< lepWs.size(); i++ ) {
        metp4.push_back( lepWs[i].q4v[1].second ) ;
   }
   
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = recoW( theJets, hadWs, NULL, btags );

   if ( leptonic && hadronic ) {
      std::vector<iReco> lepTs;
      std::vector<iReco> hadTs;
      bool findlepT  = recoTop( lepWs ,theJets, lepTs, btags );
      bool findhadT  = recoTop( hadWs ,theJets, hadTs, btags );

      if ( !findlepT || !findhadT ) return TtCollection ;
      if (algo == "ptMin" )      Algo_PtMin( lepTs, hadTs, hadWs, twIdx );
      if (algo == "beta"  )      Algo_Beta( lepTs, hadTs, lepWs, hadWs, twIdx );
      if (algo == "zero"  )      Algo_Zero( lepTs, hadTs, hadWs, twIdx, btags );
      if (algo == "kConstrain" ) Algo_KConstrain( lepTs, hadTs, hadWs, twIdx, kProb, btags );

      if ( theJets.size() >= 4 ) cout<<" lepT:"<< lepTs.size() <<" hadT:"<< hadTs.size() <<endl;

      for ( size_t i=0; i< twIdx.size(); i++ ) {
          TtCollection.push_back( lepTs[ twIdx[i][0] ] );
	  TtCollection.push_back( hadTs[ twIdx[i][1] ] );
          ResultRecord( i, twIdx, lepTs, hadTs, lepWs, hadWs, kProb, ntuples ) ;
      }

      int nBJet = 0;
      if ( btags != NULL ) {
         for (size_t i=0; i< btags->size(); i++ ) { 
             if ( (*btags)[i] ) nBJet++;
         }
      }

   }

   return TtCollection ;

}


void TtSemiEventSolution::Algo_PtMin( std::vector<iReco> lepTs,  std::vector<iReco> hadTs, std::vector<iReco> hadWs, std::vector<Idx>& twIdx ) {

      int lt = -1;
      int ht = -1;
      int lw = -1;
      int hw = -1;
      if ( lepTs.size() > 0  && hadTs.size() > 0 ) {
         double PtSum0 = 9999. ;
         for (size_t i=0; i< lepTs.size(); i++ ) {
             for (size_t j=0; j< hadTs.size(); j++) {

                 // first => W , second => b
                 if ( lepTs[i].from.second == hadTs[j].from.second ) continue;
                 int iw = hadTs[j].from.first; 
                 if ( lepTs[i].from.second == hadWs[iw].from.first  ) continue;
                 if ( lepTs[i].from.second == hadWs[iw].from.second ) continue;

                 double PxSum = lepTs[i].p4.Px() + hadTs[j].p4.Px();
                 double PySum = lepTs[i].p4.Py() + hadTs[j].p4.Py();
                 double PtSum = sqrt( (PxSum*PxSum) + (PySum*PySum) );
                 if ( PtSum > PtSum0 ) continue;
                 PtSum0 = PtSum ;
                 lt = i ;
                 ht = j ;
                 lw = lepTs[i].from.first ;
                 hw = hadTs[j].from.first ;
            }
         }
      }

      if ( lt != -1 && ht != -1 ) {
         Idx tws;
         tws.push_back(lt);
         tws.push_back(ht);
         tws.push_back(lw);
         tws.push_back(hw);

         twIdx.push_back(tws);
      }
}

void TtSemiEventSolution::Algo_Beta( std::vector<iReco> lepTs,  std::vector<iReco> hadTs, std::vector<iReco> lepWs, std::vector<iReco> hadWs, std::vector<Idx>& twIdx  ) {

      if ( lepTs.size() > 0  && hadTs.size() > 0 ) {

         for (size_t i=0; i< lepTs.size(); i++ ) {
             for (size_t j=0; j< hadTs.size(); j++) {

                 // first => W , second => b
                 if ( lepTs[i].from.second == hadTs[j].from.second ) continue;
                 int jw = hadTs[j].from.first; 
                 if ( lepTs[i].from.second == hadWs[jw].from.first  ) continue;
                 if ( lepTs[i].from.second == hadWs[jw].from.second ) continue;

                 int iw = lepTs[i].from.first; 
		 double lRelPt = tools->getRelPt( lepWs[iw].p4 , lepTs[i].p4 );
		 double hRelPt = tools->getRelPt( hadWs[jw].p4 , hadTs[j].p4 );

                 if ( lRelPt > 100. || hRelPt > 100. ) continue;

                 Idx tws;
                 tws.push_back( i );
                 tws.push_back( j );
                 tws.push_back( lepTs[i].from.first );
                 tws.push_back( hadTs[j].from.first );
                 twIdx.push_back( tws );
            }
         }
      }

}

void TtSemiEventSolution::Algo_Zero( std::vector<iReco> lepTs,  std::vector<iReco> hadTs, std::vector<iReco> hadWs, std::vector<Idx>& twIdx, std::vector<bool>* bTags ) {

    int nBJet = 0;
    if ( bTags != NULL ) {
       for (size_t i=0; i< bTags->size(); i++ ) { 
          if ( (*bTags)[i] ) nBJet++;
       }
    }

    int preComb = 0;
    for (size_t i=0; i< lepTs.size(); i++ ) {
        for (size_t j=0; j< hadTs.size(); j++) {

               // first => W , second => b
               if ( lepTs[i].from.second == hadTs[j].from.second ) continue;
               int jw = hadTs[j].from.first; 
	       if ( lepTs[i].from.second == hadWs[jw].from.first  ) continue;
	       if ( lepTs[i].from.second == hadWs[jw].from.second ) continue;

               preComb++;

               // if b jet exist, we must use it
               bool usedB = true ;
               int nbUsed = 0;
               if ( nBJet > 0 ) {
                  if ( (*bTags)[ hadWs[jw].from.first ] ) continue;
                  if ( (*bTags)[ hadWs[jw].from.second] ) continue;
                  if ( (*bTags)[ lepTs[i].from.second ] ) nbUsed++;
                  if ( (*bTags)[ hadTs[j].from.second ] ) nbUsed++;
                  if ( nbUsed == 0 && nBJet < 3 ) usedB = false;
                  if ( nbUsed != 2 && nBJet > 2 ) usedB = false;
               }
               if ( !usedB ) continue;
               // this excludes those no-btagged events if turn on b-tagger
               if ( bTags != NULL && nBJet == 0 ) continue;

               // record the result, 0:lepT, 1:hadT, 2:lepW, 3:hadW
	       Idx tws;
	       tws.push_back( i );
	       tws.push_back( j );
	       tws.push_back( lepTs[i].from.first );
	       tws.push_back( hadTs[j].from.first );

	       twIdx.push_back( tws );
        }
    }
    //cout<<" No of Combinatorics == "<< twIdx.size() <<" / "<< preComb <<endl;

}

void TtSemiEventSolution::Algo_KConstrain( std::vector<iReco> lepTs,  std::vector<iReco> hadTs, std::vector<iReco> hadWs, std::vector<Idx>& twIdx, std::vector<iProb>* kProb, std::vector<bool>* bTags ) {

    int nBJet = 0;
    if ( bTags != NULL ) {
       for (size_t i=0; i< bTags->size(); i++ ) { 
          if ( (*bTags)[i] ) nBJet++;
       }
    }

    double norm_W  = 0;
    double norm_dM = 0;
    double norm_Tt = 0;
    std::vector<double> prob_W ;
    std::vector<double> prob_dM;
    std::vector<double> prob_Tt;
    std::vector<iProb> probV ;
    for (size_t i=0; i< lepTs.size(); i++ ) {
        for (size_t j=0; j< hadTs.size(); j++) {

            // first => W , second => b
	    if ( lepTs[i].from.second == hadTs[j].from.second ) continue;
	    int jw = hadTs[j].from.first; 
	    if ( lepTs[i].from.second == hadWs[jw].from.first  ) continue;
	    if ( lepTs[i].from.second == hadWs[jw].from.second ) continue;

	    // if b jet exist, we must use it
	    bool usedB = true ;
	    int nbUsed = 0;
	    if ( nBJet > 0 ) {
               if ( (*bTags)[ hadWs[jw].from.first ] ) continue;
	       if ( (*bTags)[ hadWs[jw].from.second] ) continue;
	       if ( (*bTags)[ lepTs[i].from.second ] ) nbUsed++;
	       if ( (*bTags)[ hadTs[j].from.second ] ) nbUsed++;
	       if ( nbUsed == 0 && nBJet < 3 ) usedB = false;
	       if ( nbUsed != 2 && nBJet > 2 ) usedB = false;
            }
            if ( !usedB ) continue;

	    // this excludes those no-btagged events if turn on b-tagger
	    if ( bTags != NULL && nBJet == 0 ) continue;

	    // record the result, 0:lepT, 1:hadT, 2:lepW, 3:hadW
	    Idx tws;
	    tws.push_back( i );
	    tws.push_back( j );
	    tws.push_back( lepTs[i].from.first );
	    tws.push_back( hadTs[j].from.first );
	    twIdx.push_back( tws );
               
            // define the chi2 for ttbar from dm constrain
	    // input are come from mass-calib of algo_zero
            
	    double sgm_had = 22. ;
	    double sgm_lep = 26. ;
	    double sgm_tt  = sqrt( (sgm_had*sgm_had) + (sgm_lep*sgm_lep) ) ;
	    double dm     = hadTs[j].p4.M() - lepTs[i].p4.M() ;
	    double chi_tt = ( dm - 2.65 ) / sgm_tt ;
	    double prb_tt = exp( -0.5*chi_tt*chi_tt )  ;

	    norm_dM +=  prb_tt ;
	    norm_W  += ( lepTs[i].prob * hadTs[j].prob );
	    norm_Tt += lepTs[i].prob * hadTs[j].prob * prb_tt ;

            iProb probL ;
            probL.dM = prb_tt ;
            probL.W  = lepTs[i].prob * hadTs[j].prob ;
            probL.Tt = lepTs[i].prob * hadTs[j].prob * prb_tt ;
            probV.push_back( probL ) ;
        }
    }

    // normalized the probability of all permutation
    for (size_t i =0; i< probV.size(); i++ ) {
        probV[i].Tt = ( norm_Tt > 0.) ? probV[i].Tt / norm_Tt : 0. ;
        probV[i].W  = ( norm_W  > 0.) ? probV[i].W  / norm_W  : 0. ;
        probV[i].dM = ( norm_dM > 0.) ? probV[i].dM / norm_dM : 0. ;
        if (kProb != NULL ) kProb->push_back( probV[i] );
    }
    //cout<<"*** No of Combinatorics == "<< twIdx.size() <<" / "<< preComb<<"  Int P="<<allP <<endl;

}

// record in histograms
void TtSemiEventSolution::ResultRecord( int it, int type, std::vector<Idx> twIdx, 
                                        std::vector<iReco> lepTs, std::vector<bool>& usedLepT, 
                                        std::vector<iReco> hadTs, std::vector<bool>& usedHadT, 
                                        std::vector<iReco> lepWs, std::vector<iReco> hadWs, 
                                        std::vector<iProb>* kProb, tHisto* histos ) {

      Idx index = twIdx[it];
      int lt = index[0];
      int ht = index[1];
      int lw = index[2];
      int hw = index[3];
      double weight =  1 ;
      if ( algo == "kConstrain" && type == 1 )  weight = (*kProb)[it].Tt ;

      if ( lt != -1 && ht != -1 ) {

         double lepW_mt = tools->getMt( lepWs[lw].q4v[0].second, lepWs[lw].q4v[1].second) ;
	 double lepW_mth = ( lepW_mt >= 160 ) ? 159.9 : lepW_mt ;
         double hadW_m  = ( hadWs[hw].p4.M()  >= 160 ) ? 159.9 : hadWs[hw].p4.M() ;

         double lRelPt = tools->getRelPt( lepWs[ lw ].p4 , lepTs[lt].p4 );
         double hRelPt = tools->getRelPt( hadWs[ hw].p4 , hadTs[ht].p4 );

         double ltbeta = tools->getBeta( lepTs[lt].p4 );
         double htbeta = tools->getBeta( hadTs[ht].p4 );
         double lWdRab = tools->getdR( lepWs[ lw ].q4v[0].second, lepWs[ lw ].q4v[1].second );
         double hWdRab = tools->getdR( hadWs[ hw].q4v[0].second, hadWs[ hw].q4v[1].second );
	 double lepT_pt = ( lepTs[lt].p4.Pt() >= 480) ? 479.9 : lepTs[lt].p4.Pt() ;
	 double hadT_pt = ( hadTs[ht].p4.Pt() >= 480) ? 479.9 : hadTs[ht].p4.Pt() ;
       
         double lepTMass = lepTs[lt].p4.M() ; 
         double hadTMass = hadTs[ht].p4.M() ; 
         if (lepTMass >= 480. ) lepTMass = 479.99 ;
         if (hadTMass >= 480. ) hadTMass = 479.99 ;

         double PxTt = lepTs[lt].p4.Px() + hadTs[ht].p4.Px() ;
	 double PyTt = lepTs[lt].p4.Py() + hadTs[ht].p4.Py() ;
	 double PtTt = sqrt( (PxTt*PxTt) + (PyTt*PyTt) ) ;
	 if ( PtTt > 480. ) PtTt = 479.9 ;

         double Wtt = tools->getInvMass( lepTs[lt].p4, hadTs[ht].p4 );
         double dM = lepTMass - hadTMass ;

         // For Reco
         if ( type ==  1 ) {
            if ( algo == "kConstrian" || !usedLepT[lt] ) {
               (*histos).hWs[1]->Fill10b( lepW_mth, lRelPt, ltbeta, lWdRab );
	       (*histos).hTops[1]->Fill11c( lepTMass , lepT_pt , lepW_mth, weight );
	       (*histos).hTops[1]->Fill11d( lepTMass , hadTMass, PtTt, Wtt, dM );
	       if ( pass == 4 ) (*histos).hTops[1]->Fill11c0( lepTMass , lepT_pt, weight );
	       if ( pass == 5 ) (*histos).hTops[1]->Fill11c1( lepTMass , lepT_pt, weight );
	       if ( pass == 6 ) (*histos).hTops[1]->Fill11c2( lepTMass , lepT_pt, weight );
	       if ( pass >= 7 ) (*histos).hTops[1]->Fill11c3( lepTMass , lepT_pt, weight );
	       if ( pass == 4 ) (*histos).hTops[1]->Fill11d0( lepTMass , hadTMass, PtTt, Wtt, dM );
	       if ( pass == 5 ) (*histos).hTops[1]->Fill11d1( lepTMass , hadTMass, PtTt, Wtt, dM );
	       if ( pass == 6 ) (*histos).hTops[1]->Fill11d2( lepTMass , hadTMass, PtTt, Wtt, dM );
	       if ( pass >= 7 ) (*histos).hTops[1]->Fill11d3( lepTMass , hadTMass, PtTt, Wtt, dM );
               usedLepT[lt] = true;
            }
            if ( algo == "kConstrain" || !usedHadT[ht] ) {
               (*histos).hWs[3]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
	       (*histos).hTops[3]->Fill11c( hadTMass , hadT_pt , hadW_m, weight );
	       if ( pass == 4 ) (*histos).hTops[3]->Fill11c0( hadTMass , hadT_pt, weight );
	       if ( pass == 5 ) (*histos).hTops[3]->Fill11c1( hadTMass , hadT_pt, weight );
	       if ( pass == 6 ) (*histos).hTops[3]->Fill11c2( hadTMass , hadT_pt, weight );
	       if ( pass >= 7 ) (*histos).hTops[3]->Fill11c3( hadTMass , hadT_pt, weight );
	       usedHadT[ht] = true;
            }
         }
         // For MC matching 
         if ( type ==  0 ) {
            (*histos).hWs[0]->Fill10b( lepW_mth, lRelPt, ltbeta, lWdRab );
            (*histos).hTops[0]->Fill11c( lepTMass , lepT_pt , lepW_mth, weight );
            (*histos).hTops[0]->Fill11d( lepTMass , hadTMass, PtTt, Wtt, dM );

            if ( pass == 4 ) (*histos).hTops[0]->Fill11c0( lepTMass , lepT_pt, weight );
            if ( pass == 5 ) (*histos).hTops[0]->Fill11c1( lepTMass , lepT_pt, weight );
            if ( pass == 6 ) (*histos).hTops[0]->Fill11c2( lepTMass , lepT_pt, weight );
            if ( pass >= 7 ) (*histos).hTops[0]->Fill11c3( lepTMass , lepT_pt, weight );
            if ( pass == 4 ) (*histos).hTops[0]->Fill11d0( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 5 ) (*histos).hTops[0]->Fill11d1( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 6 ) (*histos).hTops[0]->Fill11d2( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass >= 7 ) (*histos).hTops[0]->Fill11d3( lepTMass , hadTMass, PtTt, Wtt, dM );

            (*histos).hWs[2]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            (*histos).hTops[2]->Fill11c( hadTMass , hadT_pt , hadW_m, weight );
            if ( pass == 4 ) (*histos).hTops[2]->Fill11c0( hadTMass , hadT_pt, weight );
            if ( pass == 5 ) (*histos).hTops[2]->Fill11c1( hadTMass , hadT_pt, weight );
            if ( pass == 6 ) (*histos).hTops[2]->Fill11c2( hadTMass , hadT_pt, weight );
            if ( pass >= 7 ) (*histos).hTops[2]->Fill11c3( hadTMass , hadT_pt, weight );
            
         }
      }
}

// record ntuples
void TtSemiEventSolution::ResultRecord( int it, std::vector<Idx> twIdx, 
                                        std::vector<iReco> lepTs, std::vector<iReco> hadTs, 
                                        std::vector<iReco> lepWs, std::vector<iReco> hadWs, 
                                        std::vector<iProb>* kProb, tNtuple* ntuples ) {

      Idx index = twIdx[it];
      int lt = index[0];
      int ht = index[1];
      int lw = index[2];
      int hw = index[3];
      double weight[3] = { 1 } ;
      if ( algo == "kConstrain" )  {
          weight[0] = (*kProb)[it].W ;
          weight[1] = (*kProb)[it].dM ;
          weight[2] = (*kProb)[it].Tt ;
      }
      ent_sz = static_cast<int>( twIdx.size() ) ;

      if ( lt != -1 && ht != -1 ) {
 
         double lepW_mt = tools->getMt( lepWs[lw].q4v[0].second, lepWs[lw].q4v[1].second) ;
         if ( ntuples != NULL && pass == 4 ) {

            ent_Id++;
            if ( it == 0 ) ini_Id = ent_Id ;
            (*ntuples).solTree->FillB( evt_Id, hadWs[hw].from.first,  hadWs[hw].from.second, 
                                               hadTs[ht].from.second, lepTs[lt].from.second, 
                                            -1*lepWs[lw].from.second, ini_Id,  ent_sz,    ent_Id, 
                                               pass,    weight[0], weight[1], weight[2],
                                               hadTs[ht].p4.M(), lepTs[lt].p4.M(), 
                                               hadWs[hw].p4.M(), lepW_mt  );
         }

      }
}

// updated version
void TtSemiEventSolution::HadronicTopCombinatoric( std::vector<iReco>& tCandidates, std::vector<iReco>& wCandidates, iReco t2, iReco w2 ) {

     // if it's from leptonic W, just accept it
     if ( w2.from.first == -1 ) {
        tCandidates.push_back(t2) ;
        return ;
     }

     if ( tCandidates.size() == 0 ) {
        tCandidates.push_back( t2 );
        return;
     }

     // flag the jets used
     Idx this_hadIdx(3, -1) ;
     this_hadIdx[0] = t2.from.second;  // bjet index
     this_hadIdx[1] = w2.from.first ;  // wjet index
     this_hadIdx[2] = w2.from.second ; // wjet index

     // count the same permutation => (1 2, 3) = (1 3, 2) = (2 3, 1)
     bool doubleUsed = false ;
     bool keepAll = ( algo == "kConstrain" ) ? true : false ;
     for(size_t k=0; k < tCandidates.size(); k++ ) {

        int repeat = 0;
        for(int i=0; i<3; i++) {
           int lastW = tCandidates[k].from.first ;
           if ( wCandidates[ lastW ].from.first  == this_hadIdx[i] ) repeat++ ;
           if ( wCandidates[ lastW ].from.second == this_hadIdx[i] ) repeat++ ;
           if ( tCandidates[k].from.second       == this_hadIdx[i] ) repeat++ ;
        }
        if ( repeat == 3 ) { 
            doubleUsed = true;
            if ( t2.prob > tCandidates[k].prob && !keepAll ) {
               tCandidates.erase( tCandidates.begin()+k );
               tCandidates.insert( tCandidates.begin()+k, t2 );
            }
            break;
        }
     }

     if ( !doubleUsed ||  keepAll ) {
        tCandidates.push_back( t2 );
        return;
     } 

}

