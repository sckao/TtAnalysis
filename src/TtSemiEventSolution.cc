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
static bool dmIncrease(const iReco t1, const iReco t2) { return ( t1.dm < t2.dm ); }
static bool ptDecrease(const iReco t1, const iReco t2) { return ( t1.pt > t2.pt ); }

// constructors and destructor
using namespace edm;
using namespace std;
TtSemiEventSolution::TtSemiEventSolution(const edm::ParameterSet& iConfig )
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool> ("debug");
  btag              = iConfig.getUntrackedParameter<bool> ("btag");
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
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource");
  
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
void TtSemiEventSolution::BuildSemiTt( const edm::Event& iEvent, int topo, tHisto histos ) {

   // ********************************************
   // *     Build the tt events from DATA!!!     *
   // ********************************************

   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<reco::MET> > tcmet;
   iEvent.getByLabel(tcmetSrc, tcmet);

   // 0. Event Selection Setup
   double JetEtCut = 30. ;

   // 1. Reset object containers
   isoLep.clear();
   selectedWJets.clear();
   selectedbJets.clear();
   selectedJets.clear();
   semiTt.clear();
   LorentzVector metp4;

   //2. Jet/MET selection and Reconstruct semi-leptonic Tt events 
   if ( btag ) {
      //pass = evtSelected->eventSelection( topo, JetEtCut, isoLep, selectedWJets, selectedbJets, metp4, iEvent,  "patMet" );
      //if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent(-1, selectedWJets, selectedbJets, isoLep, metp4, histos);

      std::vector<bool> bTags;
      pass = evtSelected->eventSelection( topo, JetEtCut, isoLep, selectedJets, metp4, iEvent, "patMet", "patJet", &bTags );
      cout<<" selJet Size = "<<selectedJets.size() <<"  bTager Sizer = "<< bTags.size() <<endl;
      if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent(-1, selectedJets, isoLep, metp4, histos, &bTags );

   } else {
      pass = evtSelected->eventSelection( topo, JetEtCut, isoLep, selectedJets, metp4, iEvent, "patMet", "patJet" );
      if ( algo =="ptMin" ) metp4 = (*mets)[0].p4() ;
      if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent(-1, selectedJets, isoLep, metp4, histos);
   }

   if ( semiTt.size() >= 2 )    KeepBuildInfo( true );

}


void TtSemiEventSolution::MCBuildSemiTt( const edm::Event& iEvent, int topoIdx, tHisto histos ) {

   // ********************************************
   // *   Build the tt events from MC !!!        *
   // ********************************************
 
   // retrieve the reco-objects
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<reco::MET> > tcmet;
   iEvent.getByLabel(tcmetSrc, tcmet);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::CaloJet> > jpts;
   iEvent.getByLabel(jptSrc, jpts);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   // 0. Event Pre-selection and Topology determination
   double JetEtCut = 30. ;
   int topo = evtSelected->MCEvtSelection(genParticles);
   bool build = ( topo ==  topoIdx ) ? true : false ;

   //  calculate the selection efficiency
   bool passSelect = ( pass > 3 ) ? true : false ;
   ttEff->EventEfficiency( topo, passSelect, histos.hTop );

   // 1. Reset object containers
   mcLep.clear();
   mcWJets.clear();
   mcbJets.clear();
   semiMCTt.clear();

   if ( !passSelect ) return ;
   // 2. gen-reco matching selection for leptons 
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> isoEle   = ttEle->IsoEleSelection( electrons );

   if ( topoIdx == 1 )  mcLep = MCMatching->matchMuon(genParticles, isoMuons, histos.hMuon, build);  
   if ( topoIdx == 3 )  mcLep = MCMatching->matchElectron(genParticles, isoEle, histos.hEle, build);  

   // 3. gen-reco matching selection for jets ;  For No B-tagging case only
   //    W and b jets selections 
   std::vector<const reco::Candidate*> theSelectedJets = ttJet->JetSelection( jets, isoMuons, JetEtCut );
   //std::vector<const reco::Candidate*> theSelectedJets = ttJet->JetSelection( jpts, isoMuons, JetEtCut );

   std::vector<const reco::Candidate*> additionalJets = ttJet->SoftJetSelection( jets, isoMuons, JetEtCut );
   if ( additionalJets.size() > 0  && additionalJets[0]->et() > 20. ) theSelectedJets.push_back( additionalJets[0] );

   std::vector<jmatch> mc_jets = MCMatching->matchJets(genParticles, theSelectedJets, histos.hBJet, histos.hWJet, build);  

   // 4. MET Selection
   //LorentzVector metP4 = ttMET->METfromObjects( mcLep, theSelectedJets );
   //LorentzVector metP4 = ttMET->METfromNeutrino( iEvent );
   //LorentzVector metP4 = (*tcmet)[0].p4();
   LorentzVector metP4 = (*mets)[0].p4();
   if ( algo =="ptMin" ) metP4 = (*mets)[0].p4() ;

   cout<<" sel Jet:"<< theSelectedJets.size() <<"   mc Jet:"<< mc_jets.size() <<endl;

   // 5. classify wjets and bjets
   int wQ = 0 ;
   std::vector<iReco> hadWs;
   std::vector<const reco::Candidate*> bn;
   std::vector<const reco::Candidate*> bp;
   bool findhadW = false ;
   for (size_t i =0; i< mc_jets.size(); i++ ){
       if ( abs(mc_jets[i].MomIdx) < 5 ) {
          mcWJets.push_back( mc_jets[i].trueJet ) ;
          wQ +=  mc_jets[i].MomIdx ;
          if ( mcWJets.size() == 2 ) findhadW = recoW( mcWJets, hadWs, histos.hWs[2] );
       }
       if ( mc_jets[i].MomIdx ==  5 ) bn.push_back( mc_jets[i].trueJet );
       if ( mc_jets[i].MomIdx == -5 ) bp.push_back( mc_jets[i].trueJet );
   }

   // 6. reco leptonic W
   std::vector<iReco> wSols;
   std::vector<iReco> lepWs;
   bool findlepW = recoW( mcLep, metP4, wSols, histos.hWs[0] );
   // find the right solution
   if ( wSols.size() > 1) {
      int wl = MCMatching->matchLeptonicW( genParticles, wSols );
      if ( wl > -1 ) lepWs.push_back( wSols[wl] );
   }
   if ( wSols.size() == 1 )  lepWs.push_back( wSols[0] );
   cout<<" wSol:"<< wSols.size() <<" => "<< lepWs.size() << endl;

   // 7. reco leptonic Top
   std::vector<iReco> lepTs;
   bool findlepT = false;
   if ( findlepW ) {
      std::vector<const reco::Candidate*> b_temp;
      if ( mcLep[0]->charge() == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      if ( mcLep[0]->charge() ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      findlepT = recoTop( lepWs , b_temp , lepTs, true, histos.hTops[0] );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findlepT )   semiMCTt.push_back( lepTs[0] ) ;
   }

   // 8. reco hadronic Top
   std::vector<iReco> hadTs;
   bool findhadT = false ; 
   if ( findhadW ) {
      std::vector<const reco::Candidate*> b_temp;
      if ( wQ ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      if ( wQ == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      findhadT = recoTop( hadWs , b_temp , hadTs, true, histos.hTops[2] );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findhadT ) semiMCTt.push_back( hadTs[0] );
   }

   // 9. look at top mass distribution and results
   if ( semiMCTt.size() == 2 ) {

      Idx tws(4,0);
      std::vector<Idx> twIdx( 1 , tws ) ;

      std::vector<bool> usedHadT( hadTs.size() ,false);
      std::vector<bool> usedLepT( lepTs.size() ,false);
      ResultRecord( 0, twIdx, lepTs, usedLepT, hadTs, usedHadT, lepWs, hadWs, 0, histos ) ;

      KeepBuildInfo( false );
   }

   bool matchedpass =  ( pass >= 4  ) ? true : false ;
   McRecoCompare( 1, 0, matchedpass , histos );

}

void TtSemiEventSolution::McRecoCompare( int topo, int r,  bool matchedpass, tHisto histos ) {

   int sz = static_cast<int>( AllTt.size() );
   int m = sz-1 ;
   if ( sz > 1  && AllTt[m].isData == false  ) {

      int sz1 = static_cast<int>( AllTt[r].Tt.size() )/2 ;

      for (int i1=0; i1< sz1; i1++ ) {

          int i2 = i1*2  ;

          if ( sz > 1 && matchedpass ) {
             //if (btag) {
             //   MCTruthCheck( AllTt[m].Tt , AllTt[r].Tt , i2, AllTt[m].WJ , AllTt[r].WJ , AllTt[m].bJ , AllTt[r].bJ, histos.hMObj );
             //} else {
                MCTruthCheck( AllTt[m].Tt , AllTt[r].Tt , i2, AllTt[m].WJ , AllTt[r].Js , AllTt[m].bJ , AllTt[r].Js, histos.hMObj );
             //}
          }

          accuracySemiTt( AllTt[m].Tt , AllTt[r].Tt , histos.hTop );

          /*
          if ( btag ) {
             ttEff->JetEfficiency( AllTt[r].bJ, AllTt[m].bJ, histos.hTop );
             ttEff->JetEfficiency( AllTt[r].WJ, AllTt[m].WJ, histos.hTop );
          }*/

          if ( topo == 1 ) ttEff->IsoLeptonEfficiency( AllTt[r].Ls, AllTt[m].Ls, histos.hTop );
          if ( topo == 3 ) ttEff->IsoLeptonEfficiency( AllTt[r].Ls, AllTt[m].Ls, histos.hTop );
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
        thisBuild.isData = false;
        AllTt.push_back( thisBuild );
     }

}
    
// hadronic channel 
bool TtSemiEventSolution::recoW(  std::vector<const reco::Candidate*> wjets, std::vector<iReco>& wCandidate, HTOP10* histo, std::vector<bool>* btags  ){

    bool findW = false ; 
    if ( wjets.size() < 2 ) return findW ;

    for (size_t i=0; i< wjets.size(); i++ ) {
       for (size_t j =i+1; j < wjets.size(); j++ ) {

           if ( btags != NULL ) {
              if ( (*btags)[i] || (*btags)[j] ) continue ;
           }

           double xW = wjets[i]->p4().Px() + wjets[j]->p4().Px() ;
           double yW = wjets[i]->p4().Py() + wjets[j]->p4().Py() ;
           double zW = wjets[i]->p4().Pz() + wjets[j]->p4().Pz() ;
           double EW = wjets[i]->p4().E()  + wjets[j]->p4().E() ;
	   double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;
	   if (massTest < 0. ) continue;
	   //if (massTest > 19600. ) continue;

           iReco candW;
           LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
           candW.p4   = mW;
           candW.from = make_pair(i,j);
           candW.ptr  = make_pair( wjets[i], wjets[j] );
           candW.dm   = fabs( sqrt(massTest) - 80.4 ) ;
           candW.pt   = sqrt( (xW*xW) + (yW*yW) );

	   iParton qi = make_pair( wjets[i]->pdgId() , wjets[i]->p4() );
	   iParton qj = make_pair( wjets[j]->pdgId() , wjets[j]->p4() );
	   candW.q4v.push_back( qi );
	   candW.q4v.push_back( qj );

           wCandidate.push_back( candW );

           findW = true;
       }
    }
    dmSortRecoObjects( wCandidate );

    if ( findW && histo != NULL ) {
       for (size_t i=0; i < wCandidate.size(); i++ ) {
           double wMass = tools->getInvMass( wCandidate[i].p4 );
           if( wMass >= 160 ) wMass =  159.9 ;
           histo->Fill10a( wMass );
       }
    }

    return findW;
}

// leptonic channel 
bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton, 
                          LorentzVector  metP4, std::vector<iReco>& wCandidate, HTOP10* histo ){
  

    bool findW = false ; 
    if ( lepton.size() != 1 ) return findW ;

    // use w mass constrain to solve 4 momentum of W
    double xW  = lepton[0]->p4().Px() + metP4.Px() ;
    double yW  = lepton[0]->p4().Py() + metP4.Py() ;
    double nuPt2 = ( metP4.Px()*metP4.Px() ) + ( metP4.Py()*metP4.Py() ); 
 
    double zl = lepton[0]->p4().Pz() ;
    double El = lepton[0]->p4().E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    // EtW and PtW are only for mT of W calculation 
    double EtW  = lepton[0]->pt() + metP4.E() ;
    double PtW2 = (xW*xW) + (yW*yW);
    double mtW2 = (EtW*EtW) - PtW2 ;
    if ( mtW2 < 0. ) mtW2 = 0. ;
    double mtW = sqrt( mtW2 );

    if ( (B*B) < (4.*A*C) ) return findW ;

    // 2 solutions for  z momentum of neutrino
    double nz1 = (-1.*B + sqrt(B*B - 4.*A*C )) / (2.*A) ;
    double nz2 = (-1.*B - sqrt(B*B - 4.*A*C )) / (2.*A) ;
   
    // pick a better solution!
    for (int i=1; i<3; i++ ) {
       double nz = 0.0;  
       if (i==1) nz = nz1;
       if (i==2) nz = nz2;
       double ENu2 = ( metP4.Px()*metP4.Px() ) + ( metP4.Py()*metP4.Py() ) + (nz*nz);
       double zW = lepton[0]->p4().Pz() + nz ;
       double EW = lepton[0]->p4().E()  + sqrt(ENu2) ;
       LorentzVector np4 = LorentzVector( metP4.Px(), metP4.Py(), nz, sqrt(ENu2) );
       double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

       if (massTest > 0. ) {
          iReco candW; 
          const reco::Candidate* metFakePtr = 0 ;   
          LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
          candW.p4   = mW;
          candW.from = make_pair(-1,-1);      
          candW.ptr  = make_pair( lepton[0] , metFakePtr );      
          candW.dm   = 0. ;
          candW.pt   = sqrt( (xW*xW) + (yW*yW) );
          candW.mt   = mtW ;

          iParton l1 = make_pair( lepton[0]->pdgId() , lepton[0]->p4() );
          iParton l2 = make_pair( 14 , np4 );
          candW.q4v.push_back(l1);
          candW.q4v.push_back(l2);

          wCandidate.push_back( candW );
          findW = true;
       }
    }
    //dmSortRecoObjects( wCandidate );
    sort(wCandidate.begin(), wCandidate.end(), ptDecrease );

    if ( findW && histo != NULL ) {
       for (size_t i=0; i < wCandidate.size(); i++ ) {
           double wMt = ( wCandidate[i].mt >= 160) ? 159.9 : wCandidate[i].mt ;
           histo->Fill10a( wMt );
       }
    }

    return findW;
}


// reco leptonic W with muon, MET agian ; assume pz of neutrino is zeor => MET is p4 of neutrino
bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton, LorentzVector metP4,
                                 std::vector<iReco>& wCandidate, bool FoundWSolution, HTOP10* histo ){
 
   if ( FoundWSolution ) return FoundWSolution;

   bool findW = false;
   if (  lepton.size() != 1 ) return FoundWSolution ;

   // Test the MT of lepton + MET
   double xW = lepton[0]->p4().Px() + metP4.Px() ;
   double yW = lepton[0]->p4().Py() + metP4.Py() ;
   double zW = lepton[0]->p4().Pz();

   double ptW = sqrt( (xW*xW) + (yW*yW) );
   double EtW = lepton[0]->pt() + metP4.E() ;
   double MtW2 = (EtW*EtW) - (ptW*ptW) ;
   if ( MtW2 < 0. ) return findW ;

   double pW2 =  (xW*xW) + (yW*yW) + (zW*zW) ;
   double EW2 =  lepton[0]->energy() + metP4.E() ;
   double masstest = EW2 - pW2 ;
   if ( masstest < 0 ) return findW ;

   // assume the neutrino has no Pz, only work for those exceed the Jacobian peak
   if ( MtW2 > 6464.16 ) {

      double EW = lepton[0]->p4().E() + metP4.E() ;
      double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

      if ( massTest > 0. ) {
         LorentzVector p4W = LorentzVector(xW,yW,zW,EW);
         const reco::Candidate* metFakePtr = 0 ;   
	 iReco candW; 
	 candW.p4   = p4W;
	 candW.from = make_pair(-1,-1);      
	 candW.ptr  = make_pair( lepton[0] , metFakePtr );      
	 candW.dm   = fabs( sqrt(massTest) - 80.4 ) ;
	 candW.pt   = sqrt( (xW*xW) + (yW*yW) );
	 candW.mt   = sqrt( MtW2 ) ;

	 iParton l1 = make_pair( lepton[0]->pdgId() , lepton[0]->p4() );
	 iParton l2 = make_pair( 14 , metP4 );
	 candW.q4v.push_back(l1);
	 candW.q4v.push_back(l2);

	 wCandidate.push_back( candW );
	 findW = true;    
      } 
   }

   if ( findW && histo != NULL ) {
      for (size_t i=0; i < wCandidate.size(); i++ ) {
          double wMt = ( wCandidate[i].mt >= 160) ? 159.9 : wCandidate[i].mt ;
          histo->Fill10a( wMt );
      }
   }

   return findW;

}

bool TtSemiEventSolution::recoTop( std::vector<iReco> wCandidate, std::vector<const reco::Candidate*> bjets, std::vector<iReco>& tCandidates, bool btagging, HTOP11* histo ) {

  bool findt = false ;
  tCandidates.clear();

  std::vector<iReco> bWpairList;  
  std::vector<Idx> tList;
  for(size_t i=0; i < wCandidate.size(); i++ ){
     for (size_t j=0; j< bjets.size(); j++ ){

         // avoid double use of bjets and w jets
         int bj = static_cast<int>(j);
         bool doubleUsed = (bj == wCandidate[i].from.first || bj == wCandidate[i].from.second) ? true : false ;
         if ( !btagging &&  doubleUsed ) continue;

         double xT = wCandidate[i].p4.Px() + bjets[j]->p4().Px() ;   
         double yT = wCandidate[i].p4.Py() + bjets[j]->p4().Py() ;
         double zT = wCandidate[i].p4.Pz() + bjets[j]->p4().Pz() ;
         double ET = wCandidate[i].p4.E()  + bjets[j]->p4().E()  ;
         double massTest =  (ET*ET) - (xT*xT) - (yT*yT) - (zT*zT) ;
	 LorentzVector mT = LorentzVector(xT,yT,zT,ET);

	 if ( massTest <= 0. ) continue;
         iReco bWpair;
         bWpair.p4   = mT ;
         bWpair.from = make_pair( i , j );
         bWpair.dm   = sqrt(massTest) - 171.2; 
         bWpair.pt   = sqrt( (xT*xT)+ (yT*yT) );

         iParton qb = make_pair( 5 , bjets[j]->p4() );
         bWpair.q4v.push_back( wCandidate[i].q4v[0] );
         bWpair.q4v.push_back( wCandidate[i].q4v[1] );
         bWpair.q4v.push_back( qb );
         HadronicTopCombinatoric( bWpairList, tList, bWpair, wCandidate[i], btagging );

	 findt = true;    
     }
  }

  cout<<"  sol Top:"<< bWpairList.size() <<endl;
  sort(bWpairList.begin(), bWpairList.end(), ptDecrease );
  tCandidates = bWpairList ;

  if ( findt && histo != NULL ) { 
      for ( size_t i=0; i < tCandidates.size(); i++) {    
          double theMass = ( tools->getInvMass( tCandidates[i].p4 ) >= 480. ) ? 479.9 : tools->getInvMass( tCandidates[i].p4 ) ;
	  double thePt   = ( tCandidates[i].pt  >= 480. ) ? 479.9 : tCandidates[i].pt ;
          histo->Fill11a( theMass, thePt );
          // only meaningful for no b-tagging    
	  if ( i == 0 && pass == 4 ) histo->Fill11b( theMass, thePt );
      }
  }

  return findt;
}

bool TtSemiEventSolution::recoTop( std::vector<iReco> wCandidate, std::vector<const reco::Candidate*> bjets, std::vector<iReco>& tCandidates, std::vector<bool>* btags, HTOP11* histo ) {

  bool findt = false ;
  tCandidates.clear();

  // check # of bjets
  int nBJet = 0;
  if ( btags != NULL ) {
     for (size_t i=0; i< btags->size(); i++ ) { 
         if ( (*btags)[i] ) nBJet++;
     }
  }

  std::vector<iReco> bWpairList;  
  std::vector<Idx> tList;
  for(size_t i=0; i < wCandidate.size(); i++ ){
     for (size_t j=0; j< bjets.size(); j++ ){

         // if btagging info is avaible and at least 2 jets are tagged , only apply the jets with btag
         if ( btags != NULL && nBJet  > 2 &&  !(*btags)[j]  ) continue;

         // avoid double use of bjets and w jets
         int bj = static_cast<int>(j);
         bool doubleUsed = (bj == wCandidate[i].from.first || bj == wCandidate[i].from.second) ? true : false ;
         if ( doubleUsed ) continue;

         double xT = wCandidate[i].p4.Px() + bjets[j]->p4().Px() ;   
         double yT = wCandidate[i].p4.Py() + bjets[j]->p4().Py() ;
         double zT = wCandidate[i].p4.Pz() + bjets[j]->p4().Pz() ;
         double ET = wCandidate[i].p4.E()  + bjets[j]->p4().E()  ;
         double massTest =  (ET*ET) - (xT*xT) - (yT*yT) - (zT*zT) ;
	 LorentzVector mT = LorentzVector(xT,yT,zT,ET);

	 if ( massTest <= 0. ) continue;
         iReco bWpair;
         bWpair.p4   = mT ;
         bWpair.from = make_pair( i , j );
         bWpair.dm   = sqrt(massTest) - 171.2; 
         bWpair.pt   = sqrt( (xT*xT)+ (yT*yT) );

         iParton qb = make_pair( 5 , bjets[j]->p4() );
         bWpair.q4v.push_back( wCandidate[i].q4v[0] );
         bWpair.q4v.push_back( wCandidate[i].q4v[1] );
         bWpair.q4v.push_back( qb );
         HadronicTopCombinatoric( bWpairList, tList, bWpair, wCandidate[i], false );

	 findt = true;    
     }
  }

  cout<<"  sol Top:"<< bWpairList.size() <<endl;
  sort(bWpairList.begin(), bWpairList.end(), ptDecrease );
  tCandidates = bWpairList ;

  if ( findt && histo != NULL ) { 
      for ( size_t i=0; i < tCandidates.size(); i++) {    
          double theMass = ( tools->getInvMass( tCandidates[i].p4 ) >= 480. ) ? 479.9 : tools->getInvMass( tCandidates[i].p4 ) ;
	  double thePt   = ( tCandidates[i].pt  >= 480. ) ? 479.9 : tCandidates[i].pt ;
          histo->Fill11a( theMass, thePt );
          // only meaningful for no b-tagging    
	  if ( i == 0 && pass == 4 ) histo->Fill11b( theMass, thePt );
      }
  }

  return findt;
}

// with b tagging
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<const reco::Candidate*> theWJets,
                            std::vector<const reco::Candidate*> thebJets, std::vector<const reco::Candidate*> theLep, 
                            LorentzVector metP4, tHisto histos ) {

   std::vector<iReco> TtCollection;
   std::vector<Idx> twIdx ;
   TtCollection.clear();
   if (topo != 1 && topo != -1 ) return TtCollection;

   int type = (topo == 1) ? 0 : 1; // MC : Reco
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = false;
   if ( type == 1 ) hadronic = recoW( theWJets, hadWs, histos.hWs[3] );
   if ( type == 0 ) hadronic = recoW( theWJets, hadWs, histos.hWs[2] );
   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = false;
   if ( type == 1 ) leptonic = recoW( theLep, metP4, lepWs, histos.hWs[1] );
   if ( type == 0 ) leptonic = recoW( theLep, metP4, lepWs, histos.hWs[0] );

   if ( hadronic && leptonic ) {
      std::vector<iReco> lepTs;
      bool findlepT = false;
      if( type == 0 ) findlepT = recoTop( lepWs ,thebJets, lepTs, true, histos.hTops[0] );
      if( type == 1 ) findlepT = recoTop( lepWs ,thebJets, lepTs, true, histos.hTops[1] );
      std::vector<iReco> hadTs;
      bool findhadT = false;
      if( type == 0 ) findhadT = recoTop( hadWs ,thebJets, hadTs, true, histos.hTops[2] );
      if( type == 1 ) findhadT = recoTop( hadWs ,thebJets, hadTs, true, histos.hTops[3] );

      if (algo == "dmMin" ) Algo_dmMin( lepTs, hadTs, hadWs, twIdx );
      if (algo == "ptMin" ) Algo_PtMin( lepTs, hadTs, hadWs, twIdx );
      if (algo == "beta" ) Algo_Beta( lepTs, hadTs, lepWs, hadWs, twIdx );
      if (algo == "zero" ) Algo_Zero( lepTs, hadTs, hadWs, twIdx );

      std::vector<bool> usedHadT( hadTs.size() ,false);
      std::vector<bool> usedLepT( lepTs.size() ,false);
      for ( size_t i=0; i< twIdx.size(); i++ ) {
          TtCollection.push_back( lepTs[ twIdx[i][0] ] );
	  TtCollection.push_back( hadTs[ twIdx[i][1] ] );
          ResultRecord( i, twIdx, lepTs, usedLepT, hadTs, usedHadT, lepWs, hadWs, type, histos ) ;
      }
   }
   return TtCollection ;

}

// with no b-tagging and 1 b-tagging ... maybe 2 b-tagging
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<const reco::Candidate*> theJets,
                   std::vector<const reco::Candidate*> theLep, LorentzVector metP4, tHisto histos, std::vector<bool>* btags ) {

   std::vector<iReco> TtCollection;
   std::vector<Idx> twIdx ;
   TtCollection.clear();
   if (topo != 1 && topo != -1) return TtCollection;
   

   // set MC-Reco or Real-Reco   0:MC , 1:Reco
   int type = ( topo == 1 ) ? 0:1 ;

   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = false;
   if ( type == 1 ) leptonic = recoW( theLep, metP4, lepWs, histos.hWs[1] );
   if ( type == 0 ) leptonic = recoW( theLep, metP4, lepWs, histos.hWs[0] );
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = false;
   if ( type == 1 ) hadronic = recoW( theJets, hadWs, histos.hWs[3], btags );
   if ( type == 0 ) hadronic = recoW( theJets, hadWs, histos.hWs[2], btags );

   if ( theJets.size() == 4 ) cout<<" ------ reco 4j info -------- "<<endl;
   if ( theJets.size() == 4 ) cout<<" lepW:"<< lepWs.size() <<" hadW:"<< hadWs.size() <<endl;

   if ( leptonic && hadronic ) {
      std::vector<iReco> lepTs;
      bool findlepT = false;
      if( type == 0 ) findlepT = recoTop( lepWs ,theJets, lepTs, btags, histos.hTops[0] ); // MC
      if( type == 1 ) findlepT = recoTop( lepWs ,theJets, lepTs, btags, histos.hTops[1] );
      std::vector<iReco> hadTs;
      bool findhadT = false;
      if( type == 0 ) findhadT = recoTop( hadWs ,theJets, hadTs, btags, histos.hTops[2] ); // MC
      if( type == 1 ) findhadT = recoTop( hadWs ,theJets, hadTs, btags, histos.hTops[3] );

      if (algo == "dmMin" ) Algo_dmMin( lepTs, hadTs, hadWs, twIdx );
      if (algo == "ptMin" ) Algo_PtMin( lepTs, hadTs, hadWs, twIdx );
      if (algo == "beta" ) Algo_Beta( lepTs, hadTs, lepWs, hadWs, twIdx );
      if (algo == "zero" ) Algo_Zero( lepTs, hadTs, hadWs, twIdx, btags );

      if ( theJets.size() == 4 ) cout<<" lepT:"<< lepTs.size() <<" hadT:"<< hadTs.size() <<endl;

      std::vector<bool> usedHadT( hadTs.size() ,false);
      std::vector<bool> usedLepT( lepTs.size() ,false);
      for ( size_t i=0; i< twIdx.size(); i++ ) {
          TtCollection.push_back( lepTs[ twIdx[i][0] ] );
	  TtCollection.push_back( hadTs[ twIdx[i][1] ] );
          ResultRecord( i, twIdx, lepTs, usedLepT, hadTs, usedHadT, lepWs, hadWs, type, histos ) ;
      }

      int nBJet = 0;
      if ( btags != NULL ) {
         for (size_t i=0; i< btags->size(); i++ ) { 
             if ( (*btags)[i] ) nBJet++;
         }
      }
      histos.hTop->Fill9e( theJets.size(), nBJet );

   }

   if ( theJets.size() == 4 ) cout<<" TtColl:"<<TtCollection.size() <<endl;
   return TtCollection ;

}

void TtSemiEventSolution::dmSortRecoObjects( std::vector<iReco>& objCandidates ) {

   sort(objCandidates.begin(), objCandidates.end(), dmIncrease );

}


void TtSemiEventSolution::accuracySemiTt( std::vector<iReco> ttMC, std::vector<iReco> ttReco, HTOP9* histo9 ) {

     int nb = 0;
     int nw = 0;
     int nl = 0;
     if (ttMC.size() == 0 && ttReco.size() == 0 ) {
        histo9->Fill9j(-1, 0, 0, 0) ;
        return ;
     }
     if (ttMC.size() == 0 && ttReco.size() > 0 ) {
        histo9->Fill9j(-2, 0, 0, 0) ;
        return ;
     }
     if (ttMC.size() > 0 && ttReco.size() == 0 ) {
        histo9->Fill9j(-3, 0, 0, 0) ;
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
   
     histo9->Fill9j( (ttReco[0].dm + ttReco[1].dm) , nb, nw, nl ) ;

}

void TtSemiEventSolution::MCTruthCheck( std::vector<iReco> mcTt, std::vector<iReco> rcTt, int k, 
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
       double mclepTMass = tools->getInvMass( mcTt[0].p4 ); 
       double mchadTMass = tools->getInvMass( mcTt[1].p4 ); 
       if ( pass  > 3 ) histo6->Fill6c(  mclepTMass , mchadTMass );
       if ( pass == 4 ) histo6->Fill6c0( mclepTMass , mchadTMass );
       if ( pass == 5 ) histo6->Fill6c1( mclepTMass , mchadTMass );
       if ( pass == 6 ) histo6->Fill6c2( mclepTMass , mchadTMass );
       if ( pass >= 7 ) histo6->Fill6c3( mclepTMass , mchadTMass );
   }
}

void TtSemiEventSolution::Algo_dmMin( std::vector<iReco> lepTs,  std::vector<iReco> hadTs, std::vector<iReco> hadWs, std::vector<Idx>& twIdx ) {

      int lt = -1;
      int ht = -1;
      int lw = -1;
      int hw = -1;
      if ( lepTs.size() > 0  && hadTs.size() > 0 ) {
         double dmtt = 999. ;
         for (size_t i=0; i< lepTs.size(); i++ ) {
             for (size_t j=0; j< hadTs.size(); j++) {

                 // first => W , second => b
                 if ( lepTs[i].from.second == hadTs[j].from.second ) continue;
                 int iw = hadTs[j].from.first; 
                 if ( lepTs[i].from.second == hadWs[iw].from.first  ) continue;
                 if ( lepTs[i].from.second == hadWs[iw].from.second ) continue;

                 double dmtt1 = fabs( lepTs[i].dm - hadTs[j].dm );
                 if ( dmtt1 > dmtt ) continue;
                 dmtt = dmtt1 ;
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
    cout<<" Total N of BJETs = "<< nBJet <<endl;

    int preComb = 0;
    if ( lepTs.size() > 0  && hadTs.size() > 0 ) {
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
               //cout<<" N of BJets Used => "<< nbUsed <<endl;
               if ( !usedB ) continue;
               if ( bTags != NULL && nBJet == 0 ) continue;
               //cout<<" ** keep this combination ** "<<endl;

	       Idx tws;
	       tws.push_back( i );
	       tws.push_back( j );
	       tws.push_back( lepTs[i].from.first );
	       tws.push_back( hadTs[j].from.first );
	       twIdx.push_back( tws );
           }
       }
    }
    //cout<<" No of Combinatorics == "<< twIdx.size() <<" / "<< preComb <<endl;
    cout<<" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "<<endl;

}

void TtSemiEventSolution::ResultRecord( int it, std::vector<Idx> twIdx, std::vector<iReco> lepTs, std::vector<bool>& usedLepT, 
                                        std::vector<iReco> hadTs,  std::vector<bool>& usedHadT, 
                                        std::vector<iReco> lepWs, std::vector<iReco> hadWs, int type, tHisto histos ) {

      Idx index = twIdx[it];
      int lt = index[0];
      int ht = index[1];
      int lw = index[2];
      int hw = index[3];

      if ( lt != -1 && ht != -1 ) {
	 double lepW_mt = (lepWs[ lw ].mt >= 160. ) ? 159.9 : lepWs[ lw ].mt ;
         double hadW_mass =  tools->getInvMass( hadWs[hw].p4) ;
         double hadW_m  = ( hadW_mass >= 160 ) ? 159.9 : hadW_mass ;

         double lRelPt = tools->getRelPt( lepWs[ lw ].p4 , lepTs[lt].p4 );
         double hRelPt = tools->getRelPt( hadWs[ hw].p4 , hadTs[ht].p4 );

         double ltbeta = tools->getBeta( lepTs[lt].p4 );
         double htbeta = tools->getBeta( hadTs[ht].p4 );
         double lWdRab = tools->getdR( lepWs[ lw ].q4v[0].second, lepWs[ lw ].q4v[1].second );
         double hWdRab = tools->getdR( hadWs[ hw].q4v[0].second, hadWs[ hw].q4v[1].second );
	 double lepT_pt = ( lepTs[lt].pt >= 480) ? 479.9 : lepTs[lt].pt ;
	 double hadT_pt = ( hadTs[ht].pt >= 480) ? 479.9 : hadTs[ht].pt ;
       
         double lepTMass = tools->getInvMass( lepTs[lt].p4 ); 
         double hadTMass = tools->getInvMass( hadTs[ht].p4 ); 
         if (lepTMass >= 480. ) lepTMass = 479.99 ;
         if (hadTMass >= 480. ) hadTMass = 479.99 ;

         double PxTt = lepTs[lt].p4.Px() + hadTs[ht].p4.Px() ;
	 double PyTt = lepTs[lt].p4.Py() + hadTs[ht].p4.Py() ;
	 double PtTt = sqrt( (PxTt*PxTt) + (PyTt*PyTt) ) ;
	 if ( PtTt > 480. ) PtTt = 479.9 ;

         double Wtt = tools->getInvMass( lepTs[lt].p4, hadTs[ht].p4 );
         double dM = lepTMass - hadTMass ;

         // For Reco
         if ( type ==  1 && !usedLepT[lt] ) {
            histos.hWs[1]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
            histos.hTops[1]->Fill11c( lepTMass , lepT_pt , lepW_mt);
            histos.hTops[1]->Fill11d( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 4 ) histos.hTops[1]->Fill11c0( lepTMass , lepT_pt );
            if ( pass == 5 ) histos.hTops[1]->Fill11c1( lepTMass , lepT_pt );
            if ( pass == 6 ) histos.hTops[1]->Fill11c2( lepTMass , lepT_pt );
            if ( pass >= 7 ) histos.hTops[1]->Fill11c3( lepTMass , lepT_pt );
            if ( pass == 4 ) histos.hTops[1]->Fill11d0( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 5 ) histos.hTops[1]->Fill11d1( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 6 ) histos.hTops[1]->Fill11d2( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass >= 7 ) histos.hTops[1]->Fill11d3( lepTMass , hadTMass, PtTt, Wtt, dM );
            usedLepT[lt] = true;
         }
         if ( type ==  1 && !usedHadT[ht] ) {
            histos.hWs[3]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            histos.hTops[3]->Fill11c( hadTMass , hadT_pt , hadW_m );
            if ( pass == 4 ) histos.hTops[3]->Fill11c0( hadTMass , hadT_pt );
            if ( pass == 5 ) histos.hTops[3]->Fill11c1( hadTMass , hadT_pt );
            if ( pass == 6 ) histos.hTops[3]->Fill11c2( hadTMass , hadT_pt );
            if ( pass >= 7 ) histos.hTops[3]->Fill11c3( hadTMass , hadT_pt );
            usedHadT[ht] = true;
         }
         // For MC matching 
         if ( type ==  0 ) {
            histos.hWs[0]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
            histos.hTops[0]->Fill11c( lepTMass , lepT_pt , lepW_mt );
            histos.hTops[0]->Fill11d( lepTMass , hadTMass, PtTt, Wtt, dM );

            if ( pass == 4 ) histos.hTops[0]->Fill11c0( lepTMass , lepT_pt );
            if ( pass == 5 ) histos.hTops[0]->Fill11c1( lepTMass , lepT_pt );
            if ( pass == 6 ) histos.hTops[0]->Fill11c2( lepTMass , lepT_pt );
            if ( pass >= 7 ) histos.hTops[0]->Fill11c3( lepTMass , lepT_pt );
            if ( pass == 4 ) histos.hTops[0]->Fill11d0( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 5 ) histos.hTops[0]->Fill11d1( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass == 6 ) histos.hTops[0]->Fill11d2( lepTMass , hadTMass, PtTt, Wtt, dM );
            if ( pass >= 7 ) histos.hTops[0]->Fill11d3( lepTMass , hadTMass, PtTt, Wtt, dM );

            histos.hWs[2]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            histos.hTops[2]->Fill11c( hadTMass , hadT_pt , hadW_m );
            if ( pass == 4 ) histos.hTops[2]->Fill11c0( hadTMass , hadT_pt );
            if ( pass == 5 ) histos.hTops[2]->Fill11c1( hadTMass , hadT_pt );
            if ( pass == 6 ) histos.hTops[2]->Fill11c2( hadTMass , hadT_pt );
            if ( pass >= 7 ) histos.hTops[2]->Fill11c3( hadTMass , hadT_pt );
         }

      }

}

void TtSemiEventSolution::HadronicTopCombinatoric( std::vector<iReco>& tCandidates, std::vector<Idx>& tList, iReco t2, iReco w2, bool btagging ) {

     if ( w2.from.first == -1 && w2.from.second == -1 ) {
        tCandidates.push_back(t2) ;
        return ;
     }

     if ( btagging ) {
        tCandidates.push_back(t2) ;
        return ;
     }

     Idx this_hadIdx(3, -1) ;
     this_hadIdx[0] = t2.from.second;
     this_hadIdx[1] = w2.from.first ;
     this_hadIdx[2] = w2.from.second ;

     if ( tList.size() == 0 ) {
        tList.push_back( this_hadIdx );
        tCandidates.push_back( t2 );
        return;
     }

     bool doubleUsed = false ;
     for(size_t k=0; k < tList.size(); k++ ) {

        std::vector<int> last_hadIdx = tList[k];
        int repeat = 0;
        for(int i=0; i<3; i++) {
           for(int j=0; j<3; j++) {
              if ( last_hadIdx[i] == -1 ) continue;
              if ( last_hadIdx[i] == this_hadIdx[j] ) repeat++ ;
           }
        }
        if ( repeat == 3 ) doubleUsed = true;
        if ( doubleUsed == true )  break;
     }

     if ( !doubleUsed ) {
        tList.push_back( this_hadIdx );
        tCandidates.push_back( t2 );
        return;
     }

}

