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
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
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

void TtSemiEventSolution::BuildSemiTt( const edm::Event& iEvent, const edm::EventSetup& iSetup, int topo, tHisto histos ) {

   // ********************************************
   // *     Build the tt events from DATA!!!     *
   // ********************************************

   // retrieve the reco-objects
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   // 0. Event Pre-selection
   double JetEtCut = 30. ;
   int pass = evtSelected->eventSelection(muons, electrons, jets, JetEtCut );

   // 1. Reset object containers
   isoLep.clear();
   selectedWJets.clear();
   selectedbJets.clear();
   selectedJets.clear();
   semiTt.clear();
 
   //2. Isolation Leptons Selections
   if ( topo == 1 )  isoLep = ttMuon->IsoMuonSelection( muons );
   if ( topo == 3 )  isoLep = ttEle->IsoEleSelection( electrons );

   //3. Jet selection and Reconstruct semi-leptonic Tt events 
   if ( btag ) {
      selectedWJets = ttJet->WJetSelection( jets, isoLep );
      selectedbJets = ttJet->bJetSelection( jets, isoLep );
      if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent(-1, selectedWJets, selectedbJets, isoLep, mets, histos);
   } else {
      selectedJets  = ttJet->JetSelection( jets, isoLep, JetEtCut );
      if ( pass > 3 ) semiTt = recoSemiLeptonicTtEvent(-1, selectedJets, isoLep, mets, histos);
   }

  
   // ********************************************
   // *     Check Results and Efficiency !!!     *
   // ********************************************
 
   if ( semiTt.size() == 2 ) {
      double lepTMass = tools->getInvMass( semiTt[0].p4 ); 
      double hadTMass = tools->getInvMass( semiTt[1].p4 ); 
      if (lepTMass >= 400. ) lepTMass = 399.99 ;
      if (hadTMass >= 400. ) hadTMass = 399.99 ;

      double PxTt = semiTt[0].p4.Px() + semiTt[1].p4.Px() ;
      double PyTt = semiTt[0].p4.Py() + semiTt[1].p4.Py() ;
      double PtTt = sqrt( (PxTt*PxTt) + (PyTt*PyTt) ) ;
      if ( PtTt > 400. ) PtTt = 399.9 ;

      double Wtt = tools->getInvMass( semiTt[0].p4, semiTt[1].p4 );
      double Wt2 = lepTMass + hadTMass ;

      if ( pass  > 3 ) histos.hTop->Fill9a0( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass == 4 ) histos.hTop->Fill9a1( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass == 5 ) histos.hTop->Fill9a2( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass == 6 ) histos.hTop->Fill9a3( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass >= 7 ) histos.hTop->Fill9a4( lepTMass, hadTMass, PtTt, Wtt, Wt2 );

   }

}


void TtSemiEventSolution::MCBuildSemiTt( const edm::Event& iEvent, const edm::EventSetup& iSetup, int topoIdx, tHisto histos ) {

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

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   // 0. Event Pre-selection and Topology determination
   double JetEtCut = 30. ;
   int pass = evtSelected->eventSelection(muons, electrons, jets, JetEtCut );
   int topo = evtSelected->MCEvtSelection(genParticles);
   bool build = ( topo ==  topoIdx ) ? true : false ;

   ///  calculate the selection efficiency
   bool passSelect = ( pass > 3 ) ? true : false ;
   ttEff->EventEfficiency( topo, passSelect, histos.hTop );

   // 1. Reset object containers
   mcLep.clear();
   mcWJets.clear();
   mcbJets.clear();
   semiMCTt.clear();

   // 2. gen-reco matching selection for leptons 
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection( muons );
   std::vector<const reco::Candidate*> isoEle   = ttEle->IsoEleSelection( electrons );
   if ( topoIdx == 1 )  mcLep = MCMatching->matchMuon(genParticles, isoMuons, histos.hMuon, build);  
   if ( topoIdx == 3 )  mcLep = MCMatching->matchElectron(genParticles, isoEle, histos.hEle, build);  

   // 3. gen-reco matching selection for jets ;  For No B-tagging case only
   //    W and b jets selections 
   std::vector<const pat::Jet*> selectedJets  = ttJet->JetSelection( jets, isoLep, JetEtCut );
   std::vector<jmatch> mc_jets  = MCMatching->matchJets(genParticles, selectedJets, histos.hBJet, histos.hWJet, build);  

   // 4. classify wjets and bjets
   int wQ = 0 ;
   std::vector<iReco> hadWs;
   std::vector<const pat::Jet*> bn;
   std::vector<const pat::Jet*> bp;
   bool hadronic = false ;
   for (size_t i =0; i< mc_jets.size(); i++ ){
       if ( abs(mc_jets[i].MomIdx) < 5 ) {
          mcWJets.push_back( mc_jets[i].trueJet ) ;
          wQ +=  mc_jets[i].MomIdx ;
          if ( mcWJets.size() == 2 ) hadronic = recoW( mcWJets, hadWs, histos.hWs[2] );
       }
       if ( mc_jets[i].MomIdx ==  5 ) bn.push_back( mc_jets[i].trueJet );
       if ( mc_jets[i].MomIdx == -5 ) bp.push_back( mc_jets[i].trueJet );
   }

   // 5. reco leptonic W and Top
   std::vector<iReco> wSols;
   std::vector<iReco> lepWs;
   bool leptonic = recoW( mcLep, mets, wSols, histos.hWs[0] );
   // find the right solution
   if ( wSols.size() > 1) {
      int wl = MCMatching->matchLeptonicW( genParticles, wSols );
      if ( wl > -1 ) lepWs.push_back( wSols[wl] );
   }
   if ( wSols.size() == 1 )  lepWs.push_back( wSols[0] );

   std::vector<iReco> lepTs;
   bool findlepT = false;
   if ( leptonic ) {
      std::vector<const pat::Jet*> b_temp;
      if ( mcLep[0]->charge() == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      if ( mcLep[0]->charge() ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      findlepT = recoTop( lepWs , b_temp , lepTs, true, histos.hTops[0] );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findlepT )   semiMCTt.push_back( lepTs[0] ) ;
   }

   // 6. reco hadronic W and Top
   std::vector<iReco> hadTs;
   bool findhadT = false ; 
   if (hadronic) {
      std::vector<const pat::Jet*> b_temp;
      if ( wQ ==  1 && bn.size() > 0 ) b_temp.push_back( bn[0] ) ; 
      if ( wQ == -1 && bp.size() > 0 ) b_temp.push_back( bp[0] ) ; 
      findhadT = recoTop( hadWs , b_temp , hadTs, true, histos.hTops[2] );
      if ( b_temp.size() > 0 ) mcbJets.push_back( b_temp[0] );
      if ( findhadT ) semiMCTt.push_back( hadTs[0] );
   }

   // 6. look at top mass distribution and results
   if ( semiMCTt.size() == 2 ) {

      double lepW_mt = (lepWs[0].mt >= 160) ? 159.9 : lepWs[0].mt ;
      double hadW_mass =  tools->getInvMass( hadWs[0].p4) ;
      double hadW_m  = ( hadW_mass >= 160 ) ? 159.9 : hadW_mass ;
      double lRelPt = tools->getRelPt( lepWs[0].p4 , lepTs[0].p4 );
      double hRelPt = tools->getRelPt( hadWs[0].p4 , hadTs[0].p4 );
      double ltbeta = tools->getBeta( lepTs[0].p4 );
      double htbeta = tools->getBeta( hadTs[0].p4 );
      double lWdRab = tools->getdR( lepWs[0].q4v[0].second, lepWs[0].q4v[1].second );
      double hWdRab = tools->getdR( hadWs[0].q4v[0].second, hadWs[0].q4v[1].second );
      double lepT_pt = ( lepTs[0].pt >= 400) ? 399.9 : lepTs[0].pt ;
      double hadT_pt = ( hadTs[0].pt >= 400) ? 399.9 : hadTs[0].pt ;
 
      histos.hWs[0]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
      histos.hWs[2]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );     
      histos.hTops[0]->Fill11c( lepTs[0].dm + 171.2 ,  lepT_pt );
      histos.hTops[2]->Fill11c( hadTs[0].dm + 171.2 ,  hadT_pt );

      double mcLepTMass = tools->getInvMass( semiMCTt[0].p4 ); 
      double mcHadTMass = tools->getInvMass( semiMCTt[1].p4 ); 

      double PxTt = semiMCTt[0].p4.Px() + semiMCTt[1].p4.Px() ;
      double PyTt = semiMCTt[0].p4.Py() + semiMCTt[1].p4.Py() ;
      double PtTt = sqrt( (PxTt*PxTt) + (PyTt*PyTt) ) ;

      double Wtt = tools->getInvMass( semiMCTt[0].p4, semiMCTt[1].p4 );
      double Wt2 = mcLepTMass + mcHadTMass ;

      histos.hTop->Fill9(mcLepTMass, mcHadTMass, PtTt, Wtt, Wt2 );
      if ( pass > 3 ) histos.hTop->Fill9g( mcLepTMass, mcHadTMass );    
   }
 
}

void TtSemiEventSolution::McRecoCompare( int topo, int r,  tHisto histos ) {

   int sz = static_cast<int>( AllTt.size() );
   int m = sz-1 ;

   if ( sz > 0 && AllTt[m].isData == false  ) {

      accuracySemiTt( AllTt[m].Tt , AllTt[r].Tt , histos.hTop );

      if ( sz > 1 ) {
         if (btag) {
            MCTruthCheck( AllTt[m].Tt , AllTt[r].Tt , AllTt[m].WJ , AllTt[r].WJ , AllTt[m].bJ , AllTt[r].bJ, histos.hMObj );
         } else {
            MCTruthCheck( AllTt[m].Tt , AllTt[r].Tt , AllTt[m].WJ , AllTt[r].Js , AllTt[m].bJ , AllTt[r].Js, histos.hMObj );
         }
      }

      if ( btag ) {
         ttEff->JetEfficiency( AllTt[r].bJ, AllTt[m].bJ, histos.hTop );
         ttEff->JetEfficiency( AllTt[r].WJ, AllTt[m].WJ, histos.hTop );
      }
      if ( topo == 1 ) ttEff->IsoLeptonEfficiency( AllTt[r].Ls, AllTt[m].Ls, histos.hTop );
      if ( topo == 3 ) ttEff->IsoLeptonEfficiency( AllTt[r].Ls, AllTt[m].Ls, histos.hTop );
   }
   // Reset "AllTt" containers
   AllTt.clear();

}

void TtSemiEventSolution::KeepBuildInfo( bool isData ) {

     if ( isData && semiTt.size() == 2 ) {
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
	thisBuild.bJ = mcbJets;
	thisBuild.WJ = mcWJets;
        thisBuild.isData = false;
        AllTt.push_back( thisBuild );
     }

}
    
// hadronic channel 
bool TtSemiEventSolution::recoW( std::vector<const pat::Jet*> wjets, std::vector<iReco>& wCandidate ) {

    bool findW = false ; 
    if ( wjets.size() < 2 ) return findW ;

    for (size_t i=0; i< wjets.size(); i++ ) {
       for (size_t j =i+1; j < wjets.size(); j++ ) {
           double xW = wjets[i]->p4().Px() + wjets[j]->p4().Px() ;
           double yW = wjets[i]->p4().Py() + wjets[j]->p4().Py() ;
           double zW = wjets[i]->p4().Pz() + wjets[j]->p4().Pz() ;
           double EW = wjets[i]->p4().E()  + wjets[j]->p4().E() ;
	   double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;
	   if (massTest < 0 ) continue;

           iReco candW;
           LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
           candW.p4   = mW;
           candW.from = make_pair(i,j);
           candW.ptr  = make_pair( wjets[i], wjets[j] );
           candW.dm   = fabs( sqrt(massTest) - 80.4 ) ;
           candW.pt   = sqrt( (xW*xW) + (yW*yW) );

	   iParton qi = make_pair( wjets[i]->partonFlavour() , wjets[i]->p4() );
	   iParton qj = make_pair( wjets[i]->partonFlavour() , wjets[j]->p4() );
	   candW.q4v.push_back( qi );
	   candW.q4v.push_back( qj );

           wCandidate.push_back( candW );

           findW = true;
       }
    }
    dmSortRecoObjects( wCandidate );
    return findW;
}

bool TtSemiEventSolution::recoW(  std::vector<const pat::Jet*> wjets, std::vector<iReco>& wCandidate, HTOP10* histo  ){
    
     bool hadronic = recoW( wjets, wCandidate );

     if ( hadronic ) {
        for (size_t i=0; i < wCandidate.size(); i++ ) {
            double wMass = tools->getInvMass( wCandidate[i].p4 );
            if( wMass >= 160 ) wMass =  159.9 ;
            histo->Fill10a( wMass );
        }
     }
     return hadronic ;

}


// leptonic channel 
bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton, 
                         Handle<std::vector<pat::MET> >  patMet, std::vector<iReco>& wCandidate ){
    
    std::vector<pat::MET> met ;
    for (std::vector<pat::MET>::const_iterator m1 = patMet->begin(); m1 != patMet->end(); m1++){
        met.push_back( *m1 );     
    }

    bool findW = false ; 
    if ( lepton.size() != 1 || met.size() != 1) return findW ;

    // use w mass constrain to solve 4 momentum of W
    double xW  = lepton[0]->p4().Px() + met[0].p4().Px() ;
    double yW  = lepton[0]->p4().Py() + met[0].p4().Py() ;
    double nuPt2 = (met[0].p4().Px()*met[0].p4().Px()) + (met[0].p4().Py()*met[0].p4().Py()); 
 
    double zl = lepton[0]->p4().Pz() ;
    double El = lepton[0]->p4().E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    // EtW and PtW are only for mT of W calculation 
    double EtW  = lepton[0]->et() + met[0].et() ;
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
       double ENu2 = (met[0].p4().Px()*met[0].p4().Px()) + (met[0].p4().Py()*met[0].p4().Py()) + (nz*nz);
       double zW = lepton[0]->p4().Pz() + nz ;
       double EW = lepton[0]->p4().E()  + sqrt(ENu2) ;
       LorentzVector np4 = LorentzVector( met[0].p4().Px(), met[0].p4().Py(), nz, sqrt(ENu2) );
       double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

       if (massTest > 0. ) {
          iReco candW; 
          LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
          candW.p4   = mW;
          candW.from = make_pair(-1,-1);      
          candW.ptr  = make_pair( lepton[0] , &met[0] );      
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

    return findW;
}


// reco leptonic W with muon, MET agian
bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton, Handle<std::vector<pat::MET> > met,
                                 std::vector<iReco>& wCandidate, bool FoundWSolution ){
 
   if ( FoundWSolution ) return FoundWSolution;
   bool findW = false;

   if (  lepton.size() != 1 || met->size() != 1 ) return findW ;

   // Test the MT of lepton + MET
   double xW = lepton[0]->p4().Px() + (*met)[0].p4().Px() ;
   double yW = lepton[0]->p4().Py() + (*met)[0].p4().Py() ;
   double zW = lepton[0]->p4().Pz();

   double ptW = sqrt( (xW*xW) + (yW*yW) );
   double EtW = lepton[0]->et() + (*met)[0].et() ;
   double MtW2 = (EtW*EtW) - (ptW*ptW) ;
   if ( MtW2 < 0. ) return findW ;

   // assume the neutrino has no Pz, only work for those exceed the Jacobian peak
   if ( MtW2 > 6464.16 ) {

      double EW = lepton[0]->p4().E() + (*met)[0].et() ;
      double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

      if ( massTest > 0. ) {
         LorentzVector p4W = LorentzVector(xW,yW,zW,EW);
	 iReco candW; 
	 candW.p4   = p4W;
	 candW.from = make_pair(-1,-1);      
	 candW.ptr  = make_pair( lepton[0] , &((*met)[0]) );      
	 candW.dm   = fabs( sqrt(massTest) - 80.4 ) ;
	 candW.pt   = sqrt( (xW*xW) + (yW*yW) );
	 candW.mt   = sqrt( MtW2 ) ;

	 iParton l1 = make_pair( lepton[0]->pdgId() , lepton[0]->p4() );
	 iParton l2 = make_pair( 14 , (*met)[0].p4() );
	 candW.q4v.push_back(l1);
	 candW.q4v.push_back(l2);

	 wCandidate.push_back( candW );
	 findW = true;    
      } else {
         findW = false;
      }
   }
   return findW;

}

bool TtSemiEventSolution::recoW( std::vector<const reco::Candidate*> lepton, 
                         Handle<std::vector<pat::MET> >  patMet, std::vector<iReco>& wCandidate, HTOP10* histo  ){
    
     bool foundWSolution = recoW( lepton, patMet, wCandidate );
     bool leptonic = recoW( lepton, patMet, wCandidate, foundWSolution );

     if ( leptonic ) {
        for (size_t i=0; i < wCandidate.size(); i++ ) {
            double wMt = ( wCandidate[i].mt >= 160) ? 159.9 : wCandidate[i].mt ;
            histo->Fill10a( wMt );
        }
     }
     return leptonic ;

}

bool TtSemiEventSolution::recoTop( std::vector<iReco> wCandidate, std::vector<const pat::Jet*> bjets, std::vector<iReco>& tCandidates, bool btagging ) {

  bool findt = false ;
  tCandidates.clear();

  std::vector<iReco> bWpairList;  
  for(size_t i=0; i < wCandidate.size(); i++ ){
     for (size_t j=0; j< bjets.size(); j++ ){

         int jj = static_cast<int>(j);
         bool doubleUsed = (jj == wCandidate[i].from.first || jj == wCandidate[i].from.second) ? true : false ;
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
         //bWpair.q4   = make_pair( wCandidate[i].p4 , bjets[j]->p4() );
         bWpair.dm   = sqrt(massTest) - 171.2; 
         bWpair.pt   = sqrt( (xT*xT)+ (yT*yT) );

         iParton qb = make_pair( 5 , bjets[j]->p4() );
         bWpair.q4v.push_back( wCandidate[i].q4v[0] );
         bWpair.q4v.push_back( wCandidate[i].q4v[1] );
         bWpair.q4v.push_back( qb );

	 bWpairList.push_back( bWpair );
	 findt = true;    
     }
  }
  //dmSortRecoObjects( bWpairList );
  sort(bWpairList.begin(), bWpairList.end(), ptDecrease );
  tCandidates = bWpairList ;

  return findt;
}
 

bool TtSemiEventSolution::recoTop( std::vector<iReco> wCandidate, std::vector<const pat::Jet*> bjets, std::vector<iReco>& tCandidates, bool btagging , HTOP11* histo ) {

  bool findTop = recoTop( wCandidate, bjets, tCandidates, btagging );

  if ( findTop ) { 
      for ( size_t i=0; i < tCandidates.size(); i++) {    
          double theMass = ( tools->getInvMass( tCandidates[i].p4 ) >= 400. ) ? 399.9 : tools->getInvMass( tCandidates[i].p4 ) ;
	  double thePt   = ( tCandidates[i].pt  >= 400. ) ? 399.9 : tCandidates[i].pt ;
          histo->Fill11a( theMass, thePt );    
	  if ( i == 0 ) histo->Fill11b( theMass, thePt );
      }
  }
  return findTop;
}


// with 2b tagging
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theWJets,
                            std::vector<const pat::Jet*> thebJets, std::vector<const reco::Candidate*> theLep, 
                            Handle<std::vector<pat::MET> >  mets, tHisto histos ) {

   std::vector<iReco> TtCollection;
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
   if ( type == 1 ) leptonic = recoW( theLep, mets, lepWs, histos.hWs[1] );
   if ( type == 0 ) leptonic = recoW( theLep, mets, lepWs, histos.hWs[0] );


   if ( hadronic && leptonic ) {
      std::vector<iReco> lepTs;
      bool findlepT = false;
      if( type == 0 ) findlepT = recoTop( lepWs ,thebJets, lepTs, true, histos.hTops[0] );
      if( type == 1 ) findlepT = recoTop( lepWs ,thebJets, lepTs, true, histos.hTops[1] );
      std::vector<iReco> hadTs;
      bool findhadT = false;
      if( type == 0 ) findhadT = recoTop( hadWs ,thebJets, hadTs, true, histos.hTops[2] );
      if( type == 1 ) findhadT = recoTop( hadWs ,thebJets, hadTs, true, histos.hTops[3] );

      std::pair<int, int> ttpair(-1, -1);
      std::pair<int, int> chosenW(-1, -1);
      if (findlepT && findhadT ) {
         double dmtt = 999. ;
         for (size_t i=0; i< lepTs.size(); i++ ) {
             for (size_t j=0; j< hadTs.size(); j++) {
                 if ( (lepTs[i].dm + hadTs[j].dm) > dmtt ) continue;
                 if ( lepTs[i].from.second == hadTs[j].from.second ) continue;

                 double dmtt1 = fabs( lepTs[i].dm - hadTs[j].dm );
                 if ( dmtt1 > dmtt ) continue;
                 dmtt = dmtt1 ;

                 ttpair = make_pair(i,j);
                 chosenW = make_pair(lepTs[i].from.first , hadTs[j].from.first );
            }
         }
      }
      if ( ttpair.first != -1 && ttpair.second != -1 ) {
         TtCollection.push_back( lepTs[ ttpair.first ] );
         TtCollection.push_back( hadTs[ ttpair.second ] );
	 double lepW_mt = (lepWs[ chosenW.first ].mt >= 160. ) ? 159.9 : lepWs[ chosenW.first ].mt ;
         double hadW_mass =  tools->getInvMass( hadWs[chosenW.second].p4) ;
         double hadW_m  = ( hadW_mass >= 160 ) ? 159.9 : hadW_mass ;
         double lRelPt = tools->getRelPt( lepWs[ chosenW.first ].p4 , lepTs[ttpair.first].p4 );
         double hRelPt = tools->getRelPt( hadWs[ chosenW.second].p4 , hadTs[ttpair.second].p4 );
         double ltbeta = tools->getBeta( lepTs[ttpair.first].p4 );
         double htbeta = tools->getBeta( hadTs[ttpair.second].p4 );
         double lWdRab = tools->getdR( lepWs[ chosenW.first ].q4v[0].second, lepWs[ chosenW.first ].q4v[1].second );
         double hWdRab = tools->getdR( hadWs[ chosenW.second].q4v[0].second, hadWs[ chosenW.second].q4v[1].second );
	 double lepT_pt = ( lepTs[ttpair.first].pt >= 400) ? 399.9 : lepTs[ttpair.first].pt ;
	 double hadT_pt = ( hadTs[ttpair.second].pt >= 400) ? 399.9 : hadTs[ttpair.second].pt ;

         if ( type ==  1 ) {
            histos.hWs[1]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
            histos.hWs[3]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            histos.hTops[1]->Fill11c( lepTs[ ttpair.first ].dm + 171.2 , lepT_pt );
            histos.hTops[3]->Fill11c( hadTs[ ttpair.second].dm + 171.2 , hadT_pt );
         }
         if ( type ==  0 ) {
            histos.hWs[0]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
            histos.hWs[2]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            histos.hTops[0]->Fill11c( lepTs[ ttpair.first ].dm + 171.2 , lepT_pt );
            histos.hTops[2]->Fill11c( hadTs[ ttpair.second].dm + 171.2 , hadT_pt );
         }
      }
   }
   return TtCollection ;

}

// with no b-tagging
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theJets,
                            std::vector<const reco::Candidate*> theLep, 
                            Handle<std::vector<pat::MET> >  mets, tHisto histos ) {

   std::vector<iReco> TtCollection;
   TtCollection.clear();
   if (topo != 1 && topo != -1) return TtCollection;

   // set MC-Reco or Real-Reco   0:MC , 1:Reco
   int type = ( topo == 1 ) ? 0:1 ;

   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = false;
   if ( type == 1 ) leptonic = recoW( theLep, mets, lepWs, histos.hWs[1] );
   if ( type == 0 ) leptonic = recoW( theLep, mets, lepWs, histos.hWs[0] );
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = false;
   if ( type == 1 ) hadronic = recoW( theJets, hadWs, histos.hWs[3] );
   if ( type == 0 ) hadronic = recoW( theJets, hadWs, histos.hWs[2] );

   if ( leptonic && hadronic ) {
      std::vector<iReco> lepTs;
      bool findlepT = false;
      if( type == 0 ) findlepT = recoTop( lepWs ,theJets, lepTs, true, histos.hTops[0] );
      if( type == 1 ) findlepT = recoTop( lepWs ,theJets, lepTs, true, histos.hTops[1] );
      std::vector<iReco> hadTs;
      bool findhadT = false;
      if( type == 0 ) findhadT = recoTop( hadWs ,theJets, hadTs, true, histos.hTops[2] );
      if( type == 1 ) findhadT = recoTop( hadWs ,theJets, hadTs, true, histos.hTops[3] );

      std::pair<int, int> chosenW(-1, -1);
      std::pair<int, int> ttpair(-1, -1);
      if (findlepT && findhadT ) {
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

                 ttpair = make_pair(i,j);
                 chosenW = make_pair(lepTs[i].from.first , hadTs[j].from.first );
            }
         }
      }
      if ( ttpair.first != -1 && ttpair.second != -1 ) {
         TtCollection.push_back( lepTs[ ttpair.first ] );
	 TtCollection.push_back( hadTs[ ttpair.second ] );
         double lepW_mt = (lepWs[ chosenW.first ].mt >= 160) ? 159.9 : lepWs[ chosenW.first ].mt ;
         double hadW_mass =  tools->getInvMass( hadWs[chosenW.second].p4) ;
         double hadW_m  = ( hadW_mass >= 160 ) ? 159.9 : hadW_mass ;
         double lRelPt = tools->getRelPt( lepWs[ chosenW.first ].p4 , lepTs[ttpair.first].p4 );
         double hRelPt = tools->getRelPt( hadWs[ chosenW.second ].p4 , hadTs[ttpair.second].p4 );
         double ltbeta = tools->getBeta( lepTs[ttpair.first].p4 );
         double htbeta = tools->getBeta( hadTs[ttpair.second].p4 );
         double lWdRab = tools->getdR( lepWs[ chosenW.first ].q4v[0].second, lepWs[ chosenW.first ].q4v[1].second );
         double hWdRab = tools->getdR( hadWs[ chosenW.second].q4v[0].second, hadWs[ chosenW.second].q4v[1].second );
	 double lepT_pt = ( lepTs[ttpair.first].pt >= 400) ? 399.9 : lepTs[ttpair.first].pt ;
	 double hadT_pt = ( hadTs[ttpair.second].pt >= 400) ? 399.9 : hadTs[ttpair.second].pt ;

         if ( type ==  1 ) {
            histos.hWs[1]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
            histos.hWs[3]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            histos.hTops[1]->Fill11c( lepTs[ ttpair.first ].dm + 171.2 , lepT_pt );
            histos.hTops[3]->Fill11c( hadTs[ ttpair.second].dm + 171.2 , hadT_pt );
         }
         if ( type ==  0 ) {
            histos.hWs[0]->Fill10b( lepW_mt, lRelPt, ltbeta, lWdRab );
            histos.hWs[2]->Fill10b( hadW_m , hRelPt, htbeta, hWdRab );
            histos.hTops[0]->Fill11c( lepTs[ ttpair.first ].dm + 171.2 , lepT_pt );
            histos.hTops[2]->Fill11c( hadTs[ ttpair.second].dm + 171.2 , hadT_pt );
         }
      }
   }

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


void TtSemiEventSolution::MCTruthCheck( std::vector<iReco> mcTt, std::vector<iReco> rcTt, std::vector<const pat::Jet*> mcWJets,
                                        std::vector<const pat::Jet*> rcWJets, std::vector<const pat::Jet*> mcBJets,
                                        std::vector<const pat::Jet*> rcBJets, HTOP6* histo6 ) {

   //1. mc matched reco W 
   std::vector<iReco> mcWs;
   bool mc = recoW( mcWJets, mcWs );
   //2. reco W
   std::vector<iReco> rcWs;
   bool rc = recoW( rcWJets, rcWs );

   if ( !mc || !rc ) cout<<" RECO ERROR !!!! "<<endl;

   // 3. from leptonic top
   double dRbj0 = tools->getdR( mcBJets[mcTt[0].from.second]->p4(), rcBJets[rcTt[0].from.second]->p4() ); 
   double dRt0  = tools->getdRy( mcTt[0].p4, rcTt[0].p4 );
   double dRw0  = tools->getdRy( mcWs[mcTt[0].from.first].p4, rcWs[rcTt[0].from.first].p4 );
   // 4. from hadronic top
   double dRbj1 = tools->getdR( mcBJets[mcTt[1].from.second]->p4(), rcBJets[rcTt[1].from.second]->p4() ); 
   double dRt1  = tools->getdRy( mcTt[1].p4, rcTt[1].p4 );
   double dRw1  = tools->getdRy( mcWs[mcTt[1].from.first].p4, rcWs[rcTt[1].from.first].p4 );

   histo6->Fill6a( dRt0, dRbj0, dRw0, dRt1, dRbj1, dRw1 );

   // 5. bjet and recoW from leptonic and hadronic top
   double dRbw_mch = tools->getdR( mcBJets[mcTt[1].from.second]->p4(), mcWs[mcTt[1].from.first].p4 );  
   double dRbw_rch = tools->getdR( rcBJets[rcTt[1].from.second]->p4(), rcWs[rcTt[1].from.first].p4 );  
   double dRbw_mcl = tools->getdR( mcBJets[mcTt[0].from.second]->p4(), mcWs[mcTt[0].from.first].p4 );  
   double dRbw_rcl = tools->getdR( rcBJets[rcTt[0].from.second]->p4(), rcWs[rcTt[0].from.first].p4 );  

   histo6->Fill6b( dRbw_mcl, dRbw_rcl, dRbw_mch, dRbw_rch );

}

