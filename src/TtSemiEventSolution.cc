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
static bool ptDecrease(const iReco t1, const iReco t2) { return ( t1.pt < t2.pt ); }
// constructors and destructor
using namespace edm;
using namespace std;
TtSemiEventSolution::TtSemiEventSolution(const edm::ParameterSet& iConfig)
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
  


  //evtSelected = new TtEvtSelector( iConfig );
  evtSelected = new TtEvtSelector( iConfig );
  MCMatching  = new TtMCMatching();
  ttMuon      = new TtMuon();
  ttEle       = new TtElectron();
  ttMET       = new TtMET( iConfig );
  ttJet       = new TtJet( iConfig );
  ttEff       = new TtEfficiency();

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

}

//
// member functions
//

// ------------ method called to for each event  ------------
void TtSemiEventSolution::BuildSemiTt(const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                HTOP1* histo1, HTOP3* histo3, HTOP4* histo4, HTOP7* histo7, HTOP8* histo8, HTOP9* histo9 ) {

   // retrieve the reco-objects
   
   Handle<std::vector<pat::Muon> > muons;
   iEvent.getByLabel(muonSrc, muons);

   Handle<std::vector<reco::Muon> > recomuons;
   if ( recoMuon == "paramMuons" ) {
      iEvent.getByLabel(recoMuon,"ParamGlobalMuons",recomuons);
   }
   if ( recoMuon == "muons" ) {
      iEvent.getByLabel(recoMuon,"",recomuons);
   }

   Handle<std::vector<pat::Electron> > electrons;
   iEvent.getByLabel(electronSrc, electrons);

   Handle<std::vector<pat::MET> > mets;
   iEvent.getByLabel(metSrc, mets);

   Handle<std::vector<pat::Jet> > jets;
   iEvent.getByLabel(jetSrc, jets);

   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel(genSrc, genParticles);

   Handle<std::vector<reco::GenJet> > genJets;
   iEvent.getByLabel(genJetSrc, genJets);

   Handle<CaloTowerCollection>  caloTowers;
   iEvent.getByLabel(caloSrc, caloTowers);

   // 0. select the semi-lep events and objects
   int pass = evtSelected->eventSelection(muons, electrons, jets, 20);
   int topo = evtSelected->MCEvtSelection(genParticles);
   /// calculate the selection efficiency
   bool passSelect = ( pass >= 4 ) ? true : false ;
   ttEff->EventEfficiency( topo, passSelect, histo9 );

   //1. Isolation Leptons Selections
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection(muons, histo3);
   std::vector<const reco::Candidate*> isoEle = ttEle->IsoEleSelection(electrons, histo4);
   // look at muon channel only
   isoEle.clear() ;
   //2. W and b jets selections
   std::vector<const pat::Jet*> selectedWJets = ttJet->WJetSelection( jets );
   std::vector<const pat::Jet*> selectedbJets = ttJet->bJetSelection( jets );
   /// For No B-tagging case only
   std::vector<const pat::Jet*> selectedJets  = ttJet->JetSelection( jets, isoMuons, 20., histo1 );
  
   if( isoMuons.size() == 1 ) histo4->Fill4i( isoEle.size(), selectedJets.size() );

   //3. reconstruct semi-leptonic Tt events 
   std::vector<iReco> semiTt;
   if ( pass > 3 ) {
      if (btag) {     
         semiTt = recoSemiLeptonicTtEvent(-1, selectedWJets, selectedbJets, isoMuons, isoEle, mets, histo9);
      } else {
         semiTt = recoSemiLeptonicTtEvent(-1, selectedJets, isoMuons, isoEle, mets, histo9);
      }
   }

   //bool build = ( topo == 1 || topo ==3 || topo ==4 ) ? true : false ;
   // Muon Channel Only
   bool build = ( topo == 3 ) ? true : false ;

   //4. gen-reco matching selection for leptons and jets 
   std::vector<const reco::Candidate*> mcMuons     = MCMatching->matchMuon(genParticles, muons, isoMuons, histo3, build);  
   ttMuon->matchedMuonAnalysis( mcMuons, histo3 );

   std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron(genParticles, electrons, isoEle, histo4, build);  
   ttEle->matchedElectronAnalysis( mcElectrons, histo4 );

   if ( topo == 3 ) mcElectrons.clear() ;  
   if ( topo == 4 ) mcMuons.clear() ;  

   std::vector<jmatch> mcwjets  = MCMatching->matchWJets(genParticles, jets, selectedWJets, histo8, build);
   std::vector<jmatch> mcbjets  = MCMatching->matchbJets(genParticles, jets, selectedbJets, histo7, build);
  
   //5. reconstruct semi-leptonic Tt events
   std::vector<const pat::Jet*> theWJets ; 
   for (size_t i =0; i< mcwjets.size(); i++ ){
       if ( mcwjets[i].hasMatched ) theWJets.push_back( mcwjets[i].trueJet ) ;
   }  
   std::vector<const pat::Jet*> thebJets ; 
   for (size_t i =0; i< mcbjets.size(); i++ ){
       if ( mcbjets[i].hasMatched ) thebJets.push_back( mcbjets[i].trueJet ) ;
   }

   //int topo1 = ( topo ==3 || topo ==4 ) ? 1 : 99 ;
   int topo1 = ( topo == 3 ) ? 1 : 99 ;
   std::vector<iReco> semiMCTt = recoSemiLeptonicTtEvent( topo1, theWJets, thebJets, mcMuons, mcElectrons, mets, histo9);


   //6. look at top mass distribution 
   if ( semiMCTt.size() == 2 ) {
      double mcLepTMass = getInvMass( semiMCTt[0].p4 ); 
      double mcHadTMass = getInvMass( semiMCTt[1].p4 ); 

      double PxTt = semiMCTt[0].p4.Px() + semiMCTt[1].p4.Px() ;
      double PyTt = semiMCTt[0].p4.Py() + semiMCTt[1].p4.Py() ;
      double PtTt = sqrt( (PxTt*PxTt) + (PyTt*PyTt) ) ;

      double Wtt = getInvMass( semiMCTt[0].p4, semiMCTt[1].p4 );
      double Wt2 = mcLepTMass + mcHadTMass ;

      histo9->Fill9(mcLepTMass, mcHadTMass, PtTt, Wtt, Wt2 );
      if ( semiTt.size() == 2 ) {
         if (btag) {
            MCTruthCheck(semiMCTt, semiTt, theWJets, selectedWJets, thebJets, selectedbJets, histo9 );
         } else {
            MCTruthCheck(semiMCTt, semiTt, theWJets, selectedJets, thebJets, selectedJets, histo9 );
         }
      }
   }
   if ( semiTt.size() == 2 ) {
      double lepTMass = getInvMass( semiTt[0].p4 ); 
      double hadTMass = getInvMass( semiTt[1].p4 ); 
      if (lepTMass >= 400. ) lepTMass = 399.99 ;
      if (hadTMass >= 400. ) hadTMass = 399.99 ;

      double PxTt = semiTt[0].p4.Px() + semiTt[1].p4.Px() ;
      double PyTt = semiTt[0].p4.Py() + semiTt[1].p4.Py() ;
      double PtTt = sqrt( (PxTt*PxTt) + (PyTt*PyTt) ) ;

      double Wtt = getInvMass( semiTt[0].p4, semiTt[1].p4 );
      double Wt2 = lepTMass + hadTMass ;

      if ( pass  > 3 ) histo9->Fill9a0( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass == 4 ) histo9->Fill9a1( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass == 5 ) histo9->Fill9a2( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass == 6 ) histo9->Fill9a3( lepTMass, hadTMass, PtTt, Wtt, Wt2 );
      if ( pass >= 7 ) histo9->Fill9a4( lepTMass, hadTMass, PtTt, Wtt, Wt2 );

   }

   accuracySemiTt( semiMCTt, semiTt, histo9 );

   ttEff->JetEfficiency( selectedbJets, thebJets, histo9 );
   ttEff->JetEfficiency( selectedWJets, theWJets, histo9 );
   ttEff->IsoLeptonEfficiency( isoMuons, mcMuons, histo9 );
   ttEff->IsoLeptonEfficiency( isoEle, mcElectrons, histo9 );

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
       if (i==1) nz =nz1;
       if (i==2) nz =nz2;
       double ENu2 = (met[0].p4().Px()*met[0].p4().Px()) + (met[0].p4().Py()*met[0].p4().Py()) + (nz*nz);
       double zW = lepton[0]->p4().Pz() + nz ;
       double EW = lepton[0]->p4().E() + sqrt(ENu2) ;
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
 
// with 2b tagging
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theWJets,
                            std::vector<const pat::Jet*> thebJets, std::vector<const reco::Candidate*> theMuons, 
                            std::vector<const reco::Candidate*> theElectrons, 
                            Handle<std::vector<pat::MET> >  mets, HTOP9* histo9 ) {

   std::vector<iReco> TtCollection;
   TtCollection.clear();
   if (topo != 1 && topo != -1) return TtCollection;

   // get the leptons 
   bool isMuon     = ( theMuons.size() ==1 && theElectrons.size() ==0 ) ? true : false ;
   bool isElectron = ( theMuons.size() ==0 && theElectrons.size() ==1 ) ? true : false ;
   std::vector<const reco::Candidate*> mclepton = ( isMuon && !isElectron  ) ? theMuons : theElectrons ;
   
   int type = (topo == 1 ) ? 0 : 1; // MC : Reco
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = recoW( theWJets, hadWs );
   if ( type ==  0 && hadronic ) histo9->Fill9d( getInvMass( hadWs[0].p4 ) ); // MC reco
   if ( type ==  1 && hadronic ) histo9->Fill9e( getInvMass( hadWs[0].p4 ) ); // real reco

   // reco leptonic W
   std::vector<iReco> lepWs;
   bool foundWSolution = recoW( mclepton, mets, lepWs );
   bool leptonic = recoW( mclepton, mets, lepWs, foundWSolution );
   if (type ==  0 && leptonic ) histo9->Fill9b( lepWs[0].mt ); // MC reco
   if (type ==  1 && leptonic ) histo9->Fill9c( lepWs[0].mt ); // real reco

   if ( hadronic && leptonic ) {
      std::vector<iReco> lepTs;
      bool findlepT = recoTop( lepWs ,thebJets, lepTs, true );
      std::vector<iReco> hadTs;
      bool findhadT = recoTop( hadWs ,thebJets, hadTs, true );

      if ( type ==  1 && findlepT ) histo9->Fill9q1( getInvMass( lepTs[0].p4 ) ); // real reco
      if ( type ==  0 && findlepT ) histo9->Fill9q2( getInvMass( lepTs[0].p4 ) ); // MC reco
      if ( type ==  1 && findhadT ) histo9->Fill9q3( getInvMass( hadTs[0].p4 ) ); // real reco
      if ( type ==  0 && findhadT ) histo9->Fill9q4( getInvMass( hadTs[0].p4 ) ); // MC reco

      std::pair<int, int> ttpair(-1, -1);
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
            }
         }
      }
      if ( ttpair.first != -1 && ttpair.second != -1 ) {
         TtCollection.push_back( lepTs[ ttpair.first ] );
         TtCollection.push_back( hadTs[ ttpair.second ] );
      }
   }
   return TtCollection ;

}

// with no b-tagging
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<const pat::Jet*> theJets,
                            std::vector<const reco::Candidate*> theMuons, 
                            std::vector<const reco::Candidate*> theElectrons, 
                            Handle<std::vector<pat::MET> >  mets, HTOP9* histo9 ) {

   std::vector<iReco> TtCollection;
   TtCollection.clear();
   if (topo != 1 && topo != -1) return TtCollection;

   // get the leptons 
   bool isMuon     = ( theMuons.size() ==1 && theElectrons.size() ==0 ) ? true : false ;
   bool isElectron = ( theMuons.size() ==0 && theElectrons.size() ==1 ) ? true : false ;
   std::vector<const reco::Candidate*> mclepton = ( isMuon && !isElectron  ) ? theMuons : theElectrons ;
   
   // set MC-Reco or Real-Reco
   int type = (topo == 1) ? 0:1 ;

   // reco leptonic W
   std::vector<iReco> lepWs;
   bool foundWSolution = recoW( mclepton, mets, lepWs );
   bool leptonic = recoW( mclepton, mets, lepWs, foundWSolution );
   if ( type ==  0 && leptonic ) histo9->Fill9b( lepWs[0].mt ); // MC reco
   if ( type ==  1 && leptonic ) histo9->Fill9c( lepWs[0].mt ); // real reco

   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = recoW( theJets, hadWs );
   if ( type ==  0 && hadronic ) histo9->Fill9d( getInvMass( hadWs[0].p4 ) ); // MC reco
   if ( type ==  1 && hadronic ) histo9->Fill9e( getInvMass( hadWs[0].p4 ) ); // real reco

   if ( leptonic && hadronic ) {
      std::vector<iReco> lepTs;
      bool findlepT = recoTop( lepWs ,theJets, lepTs, false );
      std::vector<iReco> hadTs;
      bool findhadT = recoTop( hadWs ,theJets, hadTs, false );

      if ( type ==  1 && findlepT ) histo9->Fill9q1( getInvMass( lepTs[0].p4 ) ); // real reco
      if ( type ==  0 && findlepT ) histo9->Fill9q2( getInvMass( lepTs[0].p4 ) ); // MC reco
      if ( type ==  1 && findhadT ) histo9->Fill9q3( getInvMass( hadTs[0].p4 ) ); // real reco
      if ( type ==  0 && findhadT ) histo9->Fill9q4( getInvMass( hadTs[0].p4 ) ); // MC reco

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
            }
         }
      }
      if ( ttpair.first != -1 && ttpair.second != -1 ) {
         TtCollection.push_back( lepTs[ ttpair.first ] );
         TtCollection.push_back( hadTs[ ttpair.second ] );
      }
   }

   return TtCollection ;

}


double TtSemiEventSolution::getInvMass( LorentzVector lv ) {

     double mom2  = (lv.Px()*lv.Px()) + (lv.Py()*lv.Py()) + (lv.Pz()*lv.Pz()) ;
     double mass2 = lv.E()*lv.E() - mom2;
     double mass = 0. ;
     if (mass2 < 0. ) mass2 = 0;
     mass = sqrt(mass2) ;
     return mass;

}

double TtSemiEventSolution::getInvMass( LorentzVector lv1, LorentzVector lv2 ) {

    double p1 = lv1.Px() + lv2.Px() ;
    double p2 = lv1.Py() + lv2.Py() ;
    double p3 = lv1.Pz() + lv2.Pz() ;
    double p4 = lv1.E() + lv2.E() ;
    double mass2 = p4*p4 - (p1*p1) - (p2*p2) - (p3*p3) ;
    if (mass2 < 0. ) mass2 = 0;

    double  mass = sqrt(mass2) ;
    return mass;
}

void TtSemiEventSolution::dmSortRecoObjects( std::vector<iReco>& objCandidates ) {

   sort(objCandidates.begin(), objCandidates.end(), dmIncrease );

}


void TtSemiEventSolution::accuracySemiTt( std::vector<iReco> ttMC, std::vector<iReco> ttReco, 
                                          HTOP9* histo9 ) {
//                                          std::vector<iReco> wMC, std::vector<iReco> wReco, HTOP9* histo9 ) {

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

/*
void TtSemiEventSolution::MCTruthCheckB( std::vector<iReco> mcTt, std::vector<iReco> rcTt, std::vector<const pat::Jet*> mcBJets,
                                        std::vector<const pat::Jet*> rcBJets, HTOP9* histo9 ) {

   // 1. bjet from leptonic top
   double dRbj0 = ttJet->getdR( mcBJets[mcTt[0].from.second]->p4(), rcBJets[rcTt[0].from.second]->p4() ); 
   double dRt0  = ttJet->getdR( mcTt[0].p4, rcTt[0].p4 );
   double dm0   = mcTt[0].dm - rcTt[0].dm ;
   // 2. bjet from hadronic top
   double dRbj1 = ttJet->getdR( mcBJets[mcTt[1].from.second]->p4(), rcBJets[rcTt[1].from.second]->p4() ); 
   double dRt1  = ttJet->getdR( mcTt[1].p4, rcTt[1].p4 );
   double dm1   = mcTt[1].dm - rcTt[1].dm ;
   histo9->Fill9n( dm0, dRt0, dRbj0, dm1, dRt1, dRbj1 );

}
*/

void TtSemiEventSolution::MCTruthCheck( std::vector<iReco> mcTt, std::vector<iReco> rcTt, std::vector<const pat::Jet*> mcWJets,
                                         std::vector<const pat::Jet*> rcWJets, std::vector<const pat::Jet*> mcBJets,
                                         std::vector<const pat::Jet*> rcBJets, HTOP9* histo9 ) {

   //1. mc matched reco W 
   std::vector<iReco> mcWs;
   bool mc = recoW( mcWJets, mcWs );
   //2. reco W
   std::vector<iReco> rcWs;
   bool rc = recoW( rcWJets, rcWs );

   if ( !mc || !rc ) cout<<" RECO ERROR !!!! "<<endl;

   double dRw = ttJet->getdR( mcWs[mcTt[1].from.first].p4, rcWs[rcTt[1].from.first].p4 );
   double dRt = ttJet->getdR( mcTt[1].p4, rcTt[1].p4 );
   double dm  = mcTt[1].dm - rcTt[1].dm ;
   int j1 = rcWs[rcTt[1].from.first].from.first ;
   int j2 = rcWs[rcTt[1].from.first].from.second ;
   double dRWjj_rc = ttJet->getdR( rcWJets[j1]->p4(), rcWJets[j2]->p4() );
   histo9->Fill9o( dRw, dRt, dm, dRWjj_rc );

   // 3. bjet from leptonic top
   double dRbj0 = ttJet->getdR( mcBJets[mcTt[0].from.second]->p4(), rcBJets[rcTt[0].from.second]->p4() ); 
   double dRt0  = ttJet->getdR( mcTt[0].p4, rcTt[0].p4 );
   double dm0   = mcTt[0].dm - rcTt[0].dm ;
   // 4. bjet from hadronic top
   double dRbj1 = ttJet->getdR( mcBJets[mcTt[1].from.second]->p4(), rcBJets[rcTt[1].from.second]->p4() ); 
   double dRt1  = ttJet->getdR( mcTt[1].p4, rcTt[1].p4 );
   double dm1   = mcTt[1].dm - rcTt[1].dm ;
   histo9->Fill9n( dm0, dRt0, dRbj0, dm1, dRt1, dRbj1 );

   // 5. bjet and recoW from leptonic and hadronic top
   double dRbw_mch = ttJet->getdR( mcBJets[mcTt[1].from.second]->p4(), mcWs[mcTt[1].from.first].p4 );  
   double dRbw_rch = ttJet->getdR( rcBJets[rcTt[1].from.second]->p4(), rcWs[rcTt[1].from.first].p4 );  
   double dRbw_mcl = ttJet->getdR( mcBJets[mcTt[0].from.second]->p4(), mcWs[mcTt[0].from.first].p4 );  
   double dRbw_rcl = ttJet->getdR( rcBJets[rcTt[0].from.second]->p4(), rcWs[rcTt[0].from.first].p4 );  
   double dRbb_mc  = ttJet->getdR( mcBJets[mcTt[0].from.second]->p4(), mcBJets[mcTt[1].from.second]->p4() );
   double dRbb_rc  = ttJet->getdR( rcBJets[rcTt[0].from.second]->p4(), rcBJets[rcTt[1].from.second]->p4() );
   histo9->Fill9p( dRbw_mcl, dRbw_rcl, dRbw_mch, dRbw_rch, dRbb_mc, dRbb_rc );

}


