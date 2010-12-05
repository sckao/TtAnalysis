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
//#include "FWCore/Framework/interface/MakerMacros.h"

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
TtSemiEventSolution::TtSemiEventSolution(const edm::ParameterSet& iConfig ): ini_Id(0), ent_Id(0), ent_sz(0) {
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool> ("debug");
  jetSetup          = iConfig.getParameter<std::vector<double> >("jetSetup");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles");
  trigTag           = iConfig.getUntrackedParameter<string> ("trigTag");
  //algo              = iConfig.getUntrackedParameter<string> ("recoAlgo");
  
  evtSelected = new TtEvtSelector( iConfig );
  MCMatching  = new TtMCMatching();
  ttMuon      = new TtMuon( iConfig );
  ttEle       = new TtElectron( iConfig );
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );
  tools       = new TtTools();

  nGoodJet = -1 ;

  for (int i = 0; i< 9; i++) { counter[i] = 0 ; }

}


TtSemiEventSolution::~TtSemiEventSolution()
{

   cout <<" hltMu= "<<counter[0] <<" Vtx = "<<counter[1] <<" 1IsoM= "<<counter[2] <<" LM= "<<counter[3] <<" Le= "<<counter[4] <<" 1J= "<<counter[5] <<" 2J= "<<counter[6] <<" 3J= "<<counter[7] <<" 4J= "<<counter[8] <<endl;

   delete evtSelected;
   delete MCMatching;
   delete ttMuon;
   delete ttEle;
   delete ttMET;
   delete ttJet;
   delete tools;
}

//
// member functions
//

// ------------ method called to for each event  ------------

// ntuple version
void TtSemiEventSolution::RecordSolutions( const edm::Event& iEvent, int topo, int evtId, int njets, SolNtp2* solTree ) {

   // ********************************************
   // *     Build the events from DATA!!!     *
   // ********************************************

   // 1. Reset object containers
   isoLep.clear();
   selectedJets.clear();
   solvedMetP4.clear();
   pvInfo.clear();
   vetoInfo.clear();

   // temperary selection for hitfit people
   /*
   bool hitfit = false ;
   if ( iEvent.id().run() == 139786 && iEvent.id().event()== 30766631 )  hitfit = true; 
   if ( iEvent.id().run() == 140124 && iEvent.id().event()== 1749068 )  hitfit = true; 
   if ( iEvent.id().run() == 141960 && iEvent.id().event()== 36380903 )  hitfit = true; 
   if ( iEvent.id().run() == 142137 && iEvent.id().event()== 115049658 )  hitfit = true; 
   if ( iEvent.id().run() == 142524 && iEvent.id().event()== 119786558 )  hitfit = true; 
   if ( iEvent.id().run() == 142528 && iEvent.id().event()== 528865367 )  hitfit = true; 
   evtId =  iEvent.id().event() ;
   */

   // 2. Event Selection
   //bool passtrig = evtSelected->TriggerSelection( iEvent, "HLT_Mu15v1" ); 
   bool passtrig = evtSelected->TriggerSelection( iEvent, trigTag ); 
   //bool passtrig = true ; 

   bool goodVtx = evtSelected->VertexSelection( iEvent, pvInfo ) ;

   nGoodJet = evtSelected->eventSelection( topo, isoLep, selectedJets, solvedMetP4, vetoInfo, iEvent, "patMet" );

   //cout<<" event selected w/ jet "<<selectedJets.size()<<"  w/ lepton "<< isoLep.size() <<endl;
   int nLmu = 0;
   int nLel = 0;
   for ( size_t j=0 ; j< vetoInfo.size(); j++ ) {
       if ( vetoInfo[j].pdgId == 11 ) nLel++;
       if ( vetoInfo[j].pdgId == 13 ) nLmu++;
   }

   //int failIdx = 0 ;
   // sync exerc
   
   if ( passtrig )                                                     counter[0]++ ;
   if ( passtrig && goodVtx  )                                         counter[1]++ ;
   if ( passtrig && goodVtx && nGoodJet > -1 )                         counter[2]++ ;
   if ( passtrig && goodVtx && nGoodJet > -1 && nLmu == 0 )            counter[3]++ ;
   if ( passtrig && goodVtx && nGoodJet > -1 && vetoInfo.size() == 0 ) counter[4]++ ;
   if ( passtrig && goodVtx && nGoodJet >= 1 && vetoInfo.size() == 0 ) counter[5]++ ;
   if ( passtrig && goodVtx && nGoodJet >= 2 && vetoInfo.size() == 0 ) counter[6]++ ;
   if ( passtrig && goodVtx && nGoodJet >= 3 && vetoInfo.size() == 0 ) counter[7]++ ;
   if ( passtrig && goodVtx && nGoodJet >= 4 && vetoInfo.size() == 0 ) counter[8]++ ;
   
   /*
   if ( passtrig )                                                     failIdx++ ;
   if ( passtrig && goodVtx  )                                         failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet > -1 )                         failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet > -1 && nLmu == 0 )            failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet > -1 && vetoInfo.size() == 0 ) failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet >= 1 && vetoInfo.size() == 0 ) failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet >= 2 && vetoInfo.size() == 0 ) failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet >= 3 && vetoInfo.size() == 0 ) failIdx++ ;
   if ( passtrig && goodVtx && nGoodJet >= 4 && vetoInfo.size() == 0 ) failIdx++ ;
   */
    
   //double isoLepPt = -1 ;
   //if ( isoLep.size() > 0 ) isoLepPt = isoLep[0].p4.Pt() ;
   //cout<< iEvent.id().run() <<" : "<<iEvent.id().event()<<" : "<< iEvent.id().luminosityBlock() <<" : "<< isoLepPt <<
   //" : "<< failIdx+1 << endl;

   // selection for 2 jet sample
   bool getEvent = ( nGoodJet >= njets ) ? true : false ;
   if ( njets == 2 ) getEvent = ( nGoodJet == njets ) ? true : false ;

   if ( !goodVtx ) getEvent = false ;
   if ( !passtrig ) getEvent = false ;
   if ( vetoInfo.size() >  0 ) getEvent = false ;

   if ( getEvent ) {

      //cout<< iEvent.id().run() <<" : "<<iEvent.id().event()<<" : "<< iEvent.id().luminosityBlock() <<" : "<< isoLep[0].p4.Pt() <<endl;
      // 3. Reconstruct leptonic W 
     //cout<<"    -> "<<njets<<" jets event " <<endl ;
     std::vector<iReco> lepWs;
     bool leptonic = recoW( isoLep[0], solvedMetP4[0], lepWs );
     if ( leptonic ) solvedMetP4.clear()    ;
     for ( size_t i =0; i< lepWs.size(); i++ ) {
         solvedMetP4.push_back( lepWs[i].q4v[1].second ) ;
     }
     // save the events without considering neutrino pz solution => events with only 1 neutrino pz( = 0 )  
     solTree->FillB1( evtId, selectedJets, solvedMetP4, isoLep, pvInfo, vetoInfo );

   }
   //cout<<" ------ reco finished --------------"<<endl;

}

// topoIdx : 4:tau+jets  3:electron+jest  2:di-lep,  1:muon+jets,  0:hadron,  -1:Non-Tt event
// ntuple version
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

   if ( nGoodJet < 4 ) return;

   std::vector<const reco::Candidate*> genCollects = MCMatching->GenTtCollection( genParticles ) ;

   if ( topoIdx == 1 )  mcLep = MCMatching->matchMuon(genCollects, isoLep );  
   if ( topoIdx == 3 )  mcLep = MCMatching->matchElectron(genCollects, isoLep );  
   std::vector<jmatch> mc_jets = MCMatching->matchJets(genCollects, selectedJets );  

   // 2. classify wjets and bjets ..... parton Btagging
   int wQ = 0 ;
   std::vector<int> wjid;
   int bjh = -1;
   int bjl = -1;
   std::vector<iReco> hadWs;
   std::vector<LorentzVector> bn;
   std::vector<LorentzVector> bp;
   std::vector<bool> mcbTags(2, true) ;
   bool findhadW = false ;
   for (size_t i =0; i< mc_jets.size(); i++ ){

       if ( abs( mc_jets[i].MomIdx ) < 5  ) {
          mcWJets.push_back( selectedJets[ mc_jets[i].Idx ] ) ;
          wjid.push_back( mc_jets[i].Idx );
          wQ +=  mc_jets[i].MomIdx ;
          // reco hadronic W
          if ( mcWJets.size() == 2 ) findhadW = recoW( mcWJets, hadWs );
       }
       if ( mc_jets[i].MomIdx ==  5 ) bn.push_back( mc_jets[i].p4 );
       if ( mc_jets[i].MomIdx == -5 ) bp.push_back( mc_jets[i].p4 );
       if ( mcLep.size() == 0 ) continue;
       if ( mcLep[0].charge == -1 && mc_jets[i].MomIdx == -5 ) bjl = mc_jets[i].Idx ; 
       if ( mcLep[0].charge == -1 && mc_jets[i].MomIdx ==  5 ) bjh = mc_jets[i].Idx ; 
       if ( mcLep[0].charge ==  1 && mc_jets[i].MomIdx ==  5 ) bjl = mc_jets[i].Idx ; 
       if ( mcLep[0].charge ==  1 && mc_jets[i].MomIdx == -5 ) bjh = mc_jets[i].Idx ; 
   }

   // 3. reco leptonic W
   std::vector<iReco> wSols;
   std::vector<iReco> lepWs;
   bool findlepW = recoW( mcLep[0], solvedMetP4[0], wSols );
   // find the right solution

   int wl  = -1 ;
   if ( wSols.size() >= 1) {
      wl = MCMatching->matchLeptonicW( genParticles, wSols );
      if ( wl > -1 )  lepWs.push_back( wSols[wl] );
   }

   // 4. reco leptonic Top
   std::vector<iReco> lepTs;
   bool findlepT = false;
   if ( findlepW && wl > -1 ) {
      if ( mcLep[0].charge == -1 && bp.size() > 0 ) findlepT = true ; 
      if ( mcLep[0].charge ==  1 && bn.size() > 0 ) findlepT = true ; 
   }

   // 5. reco hadronic Top
   std::vector<iReco> hadTs;
   bool findhadT = false ; 
   if ( findhadW ) {
      if ( wQ ==  1 && bn.size() > 0 ) findhadT = true ; 
      if ( wQ == -1 && bp.size() > 0 ) findhadT = true ; 
   }

   // 6. look at top mass distribution and results
   if ( findlepT && findhadT &&  ntuples != NULL ) {
       //cout<<" matchID = ( "<< wjid[0]<<","<<wjid[1]<<","<<bjh<<","<<bjl<<","<<wSols[wl].from.second<<" )"<<endl;
       int matchId[6] = { wjid[0], wjid[1], bjh, bjl, -1*wSols[wl].from.second, 0 } ;
       /*
       cout<<" j0 = "<< selectedJets[ wjid[0] ].p4.Pt() <<" / "<< selectedJets[ wjid[0] ].p4.Pz() <<endl;
       cout<<" j1 = "<< selectedJets[ wjid[1] ].p4.Pt() <<" / "<< selectedJets[ wjid[1] ].p4.Pz() <<endl;
       cout<<" j2 = "<< selectedJets[ bjh ].p4.Pt() <<" / "<< selectedJets[ bjh ].p4.Pz() <<endl;
       cout<<" j3 = "<< selectedJets[ bjl ].p4.Pt() <<" / "<< selectedJets[ bjl ].p4.Pz() <<endl;
       cout<<" nu = "<< lepWs[0].q4v[1].second.Pt() <<" / "<< lepWs[0].q4v[1].second.Pz() <<endl;
       cout<<" mu = "<< mcLep[0].p4.Pt() <<" / "<< mcLep[0].p4.Pz() <<endl;
       */
       (*ntuples).mcmTree->FillB( evtId, matchId, selectedJets, lepWs[0].q4v[1].second,  mcLep[0].p4 );
       (*ntuples).genTree->FillB1( evtId, genCollects );

   }
   //cout<<" ====  MCMatching : size of lepT:"<< lepTs.size() <<" , hadT:"<<hadTs.size()  <<endl;

}

// hadronic channel 
bool TtSemiEventSolution::recoW(  std::vector<ttCandidate>& wjets, std::vector<iReco>& wCandidate, HTOP10* histo ){

    bool findW = false ; 
    if ( wjets.size() < 2 ) return findW ;

    double norm_p = 0;
    for (size_t i=0; i< wjets.size(); i++ ) {
       for (size_t j =i+1; j < wjets.size(); j++ ) {

           if ( wjets[i].cuts[0] > 5 || wjets[j].cuts[0] > 5 ) continue;

           LorentzVector m4W = wjets[i].p4 + wjets[j].p4 ;
	   if ( m4W.M2() < 0. ) continue;

           double chi = (m4W.M() - 80.4) / (80.4*0.1) ;
           double wProb = ( m4W.M() <= 160.8 ) ? exp( -0.5*chi*chi ) : 0. ;
           norm_p += wProb ;
          
           iReco candW;
           candW.p4   = m4W;
           candW.from = make_pair(i,j);
           candW.prob = exp( -0.5*chi*chi ) ;

	   iParton qi = make_pair( wjets[i].pdgId , wjets[i].p4 );
	   iParton qj = make_pair( wjets[j].pdgId , wjets[j].p4 );
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
bool TtSemiEventSolution::recoW( ttCandidate& lepton,  LorentzVector metp4, std::vector<iReco>& wCandidate, HTOP10* histo ){
  
    bool findW = false ; 

    // use w mass constrain to solve 4 momentum of W
    double xW  = lepton.p4.Px() + metp4.Px() ;
    double yW  = lepton.p4.Py() + metp4.Py() ;
    double nuPt2 = ( metp4.Px()*metp4.Px() ) + ( metp4.Py()*metp4.Py() ); 
 
    double zl = lepton.p4.Pz() ;
    double El = lepton.p4.E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    // EtW and PtW are only for mT of W calculation 
    double EtW  = lepton.p4.Pt() + metp4.E() ;
    double PtW2 = (xW*xW) + (yW*yW);
    double MtW2 = (EtW*EtW) - PtW2 ;
    if ( MtW2 < 0. ) return findW ;

    if ( (B*B) < (4.*A*C) ) return findW ;

    // 2 solutions for  z momentum of neutrino
    double nz1 = (-1.*B + sqrt(B*B - 4.*A*C )) / (2.*A) ;
    double nz2 = (-1.*B - sqrt(B*B - 4.*A*C )) / (2.*A) ;
   
    // get solutions!
    for (int i=1; i<3; i++ ) {
       double nz = ( i == 1) ? nz1 : nz2 ;  
       double ENu2 = ( metp4.Px()*metp4.Px() ) + ( metp4.Py()*metp4.Py() ) + (nz*nz);
       double zW = lepton.p4.Pz() + nz ;
       double EW = lepton.p4.E()  + sqrt(ENu2) ;
       LorentzVector np4 = LorentzVector( metp4.Px(), metp4.Py(), nz, sqrt(ENu2) );
       double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

       if (massTest > 0. ) {
          iReco candW; 
          LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
          candW.p4   = mW;
          candW.from = make_pair( -1, -1*i );      

          // probability must be 0.5 ...because 2 pz solution
          candW.prob = 0.5;

          iParton l1 = make_pair( lepton.pdgId , lepton.p4 );
          int neuId =  ( lepton.pdgId > 0 ) ? lepton.pdgId+1 : lepton.pdgId-1 ;
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
bool TtSemiEventSolution::recoW( ttCandidate& lepton, LorentzVector metp4,
                                 std::vector<iReco>& wCandidate, bool FoundWSolution, HTOP10* histo ){
 
   if ( FoundWSolution ) return FoundWSolution;

   bool findW = false;

   // Test the MT of lepton + MET
   double xW = lepton.p4.Px() + metp4.Px() ;
   double yW = lepton.p4.Py() + metp4.Py() ;
   double zW = lepton.p4.Pz();

   double ptW = sqrt( (xW*xW) + (yW*yW) );
   double EtW = lepton.p4.Pt() + metp4.E() ;
   double MtW2 = (EtW*EtW) - (ptW*ptW) ;
   if ( MtW2 < 0. ) return findW ;

   double pW2 =  (xW*xW) + (yW*yW) + (zW*zW) ;
   double EW2 =  lepton.p4.E() + metp4.E() ;
   double masstest = EW2 - pW2 ;
   if ( masstest < 0 ) return findW ;

   // assume the neutrino has no Pz, only work for those exceed the Jacobian peak
   if ( MtW2 > 6464.16 ) {

      double EW = lepton.p4.E() + metp4.E() ;
      double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

      if ( massTest > 0. ) {
	 iReco candW; 
         LorentzVector p4W = LorentzVector(xW,yW,zW,EW);
	 candW.p4   = p4W;
	 candW.from = make_pair(-1,-1);      

         candW.prob = 1. ;

         iParton l1 = make_pair( lepton.pdgId , lepton.p4 );
         int neuId =  ( lepton.pdgId > 0 ) ? lepton.pdgId+1 : lepton.pdgId-1 ;
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

