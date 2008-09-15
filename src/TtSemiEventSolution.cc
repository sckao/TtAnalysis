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

// constructors and destructor
using namespace edm;
using namespace std;
TtSemiEventSolution::TtSemiEventSolution(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  debug             = iConfig.getUntrackedParameter<bool>   ("debug");
  muonSrc           = iConfig.getParameter<edm::InputTag> ("muonSource");
  electronSrc       = iConfig.getParameter<edm::InputTag> ("electronSource");
  metSrc            = iConfig.getParameter<edm::InputTag> ("metSource");
  jetSrc            = iConfig.getParameter<edm::InputTag> ("jetSource");
  genJetSrc         = iConfig.getParameter<edm::InputTag> ("genJetSource");
  genSrc            = iConfig.getParameter<edm::InputTag> ("genParticles"); 
  recoMuon          = iConfig.getUntrackedParameter<string> ("recoMuons");
  caloSrc           = iConfig.getParameter<edm::InputTag> ("caloSource"); 


  evtSelected = new TtEvtSelector();
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
                                      HTOP3* histo3, HTOP4* histo4, HTOP7* histo7, HTOP8* histo8, HTOP9* histo9 ) {

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
   bool pass = evtSelected->eventSelection(muons, electrons, jets);
   int  topo = evtSelected->MCEvtSelection(genParticles);
   /// calculate the selection efficiency
   ttEff->EventEfficiency( topo, pass, histo9 );

   //1. Isolation Leptons Selections
   std::vector<const reco::Candidate*> isoMuons = ttMuon->IsoMuonSelection(muons, histo3);
   std::vector<const reco::Candidate*> isoEle = ttEle->IsoEleSelection(electrons, histo4);
   //2. W and b jets selections
   std::vector<pat::Jet> selectedWJets = evtSelected->WJetSelection(jets);
   std::vector<pat::Jet> selectedbJets = evtSelected->bJetSelection(jets);
   
   //3. reconstruct semi-leptonic Tt events 
   std::vector<iReco> semiTt;
   if ( pass ) {
      semiTt = recoSemiLeptonicTtEvent(-1, selectedWJets, selectedbJets, isoMuons, isoEle, mets, histo9);
      ttEff->EventShape(topo, isoMuons.size(), isoEle.size(), selectedbJets.size(), selectedWJets.size(), histo9);
   }

   bool build = ( topo == 1 ) ? true : false ;

   //4. gen-reco matching selection for leptons and jets 
   std::vector<const reco::Candidate*> mcMuons     = MCMatching->matchMuon(genParticles, muons, isoMuons, histo3, build);  
   std::vector<const reco::Candidate*> mcElectrons = MCMatching->matchElectron(genParticles, electrons, isoEle, histo4, build);  
   std::vector<jmatch> mcwjets  = MCMatching->matchWJets(genParticles, jets, selectedWJets, histo8, build);
   std::vector<jmatch> mcbjets  = MCMatching->matchbJets(genParticles, jets, selectedbJets, histo7, build);
 
  
   //5. reconstruct semi-leptonic Tt events
   std::vector<pat::Jet> theWJets ; 
   for (size_t i =0; i< mcwjets.size(); i++ ){
       if ( mcwjets[i].hasMatched ) theWJets.push_back( mcwjets[i].trueJet ) ;
   }  
   std::vector<pat::Jet> thebJets ; 
   for (size_t i =0; i< mcbjets.size(); i++ ){
       if ( mcbjets[i].hasMatched ) thebJets.push_back( mcbjets[i].trueJet ) ;
   }
   std::vector<iReco> semiMCTt = recoSemiLeptonicTtEvent( topo, theWJets, thebJets, mcMuons, mcElectrons, mets, histo9);

   //6. look at top mass distribution 
   for(size_t i=0; i< semiMCTt.size(); i++ ){
      double tMass = getInvMass( semiMCTt[i].p4 ); 
      histo9->Fill9(tMass);
   }
   for(size_t i=0; i< semiTt.size(); i++ ){
      double tMass = getInvMass( semiTt[i].p4 ); 
      histo9->Fill9a(tMass);
   }

   ttEff->JetEfficiency( selectedbJets, thebJets, histo9 );
   ttEff->JetEfficiency( selectedWJets, theWJets, histo9 );
   ttEff->IsoLeptonEfficiency( isoMuons, mcMuons, histo9 );
   ttEff->IsoLeptonEfficiency( isoEle, mcElectrons, histo9 );

}

// hadronic channel 
bool TtSemiEventSolution::recoW( int type, std::vector<pat::Jet> wjets, std::vector<iReco>& wCandidate, HTOP9* histo9 ) {

    bool findW = false ; 
    if ( wjets.size() < 2 ) return findW ;

    for (size_t i=0; i< wjets.size(); i++ ) {
       for (size_t j =i+1; j < wjets.size(); j++ ) {
           double xW = wjets[i].p4().Px() + wjets[j].p4().Px() ;
           double yW = wjets[i].p4().Py() + wjets[j].p4().Py() ;
           double zW = wjets[i].p4().Pz() + wjets[j].p4().Pz() ;
           double EW = wjets[i].p4().E()  + wjets[j].p4().E() ;
	   double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;
	   if (massTest < 0 ) continue;
           iReco candW;
           LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
           std::pair<int,int> jetPair(i,j);
           candW.p4   = mW;
           candW.from = jetPair ;
           candW.dm   = fabs( sqrt(massTest) - 80.4 ) ;
           wCandidate.push_back( candW );
           if (type ==  0 )   histo9->Fill9d( sqrt(massTest) ); // MC reco
           if (type ==  1 )   histo9->Fill9e( sqrt(massTest) ); // real reco
           findW = true;
       }
    }
    dmSortRecoObjects( wCandidate );
    return findW;
}

// semi-leptonic channel 
bool TtSemiEventSolution::recoW( int type, std::vector<const reco::Candidate*> lepton, 
                         Handle<std::vector<pat::MET> >  patMet, std::vector<iReco>& wCandidate, HTOP9* histo9 ){
    
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
       double massTest =  (EW*EW) - (xW*xW) - (yW*yW) - (zW*zW) ;

       if (massTest > 0. ) {
          iReco candW; 
          LorentzVector mW = LorentzVector(xW,yW,zW,EW) ;
          std::pair<int,int> lvPair(0,0);
          candW.p4   = mW;
          candW.from = lvPair;      
          candW.dm   = fabs( sqrt(massTest) - 80.4 );
          candW.mt   = mtW ;
          wCandidate.push_back( candW );
          findW = true;
       }
    }
    
    dmSortRecoObjects( wCandidate );
    if (type == 0 )histo9->Fill9b( wCandidate[0].mt );  // mc reco
    if (type == 1 )histo9->Fill9c( wCandidate[0].mt );  // real reco

    return findW;
}

bool TtSemiEventSolution::recoTop( std::vector<iReco> wCandidate, std::vector<pat::Jet> bjets, std::vector<iReco>& tCandidates ) {

  bool findt = false ;
  tCandidates.clear();
  //if ( wCandidate.size() < 2 || bjets.size() < 2 ) return findt;

  //std::vector<iReco> bW1;
  //std::vector<iReco> bW2;
  std::vector<iReco> bWpairList;  
  int widx = -1;
  for(std::vector<iReco>::const_iterator i1 = wCandidate.begin(); i1 != wCandidate.end(); i1++ ){
     widx++;
     int bj = -1;
     for (std::vector<pat::Jet>::const_iterator i2 = bjets.begin(); i2 != bjets.end(); i2++ ){
         bj++;    
         double xT = i1->p4.Px() + (i2->p4()).Px() ;   
         double yT = i1->p4.Py() + (i2->p4()).Py() ;
         double zT = i1->p4.Pz() + (i2->p4()).Pz() ;
         double ET = i1->p4.E()  + (i2->p4()).E()  ;
	 LorentzVector mT = LorentzVector(xT,yT,zT,ET);
         double massTest =  (ET*ET) - (xT*xT) - (yT*yT) - (zT*zT) ;

	 if ( massTest <= 0. ) continue;
         iReco bWpair;
         bWpair.p4   = mT ;
         bWpair.from = std::pair<int,int>( widx , bj );
         bWpair.dm   = fabs(sqrt(massTest) - 174.); 
	 bWpairList.push_back( bWpair );
	 findt = true;    
     }
     //if (widx == 1) bW1 = bWpairList;
     //if (widx == 2) bW2 = bWpairList;
  }
  dmSortRecoObjects( bWpairList );
  tCandidates = bWpairList ;
  /*
  int sz = ( bW1.size() > bW2.size() ) ? bW1.size() : bW2.size() ;
  for(int i=0; i < sz; i++) {

     if (bW1[i].from.second != bW2[i].from.second ) {
        tCandidates.push_back(bW1[i]);
        tCandidates.push_back(bW2[i]);
        break;
     } else {
        if ( bW1[i].dm < bW2[i].dm ) {
           tCandidates.push_back(bW1[i]);
           tCandidates.push_back(bW2[i+1]);
           break;
        }
        if ( bW1[i].dm > bW2[i].dm ) {
           tCandidates.push_back(bW1[i+1]);
           tCandidates.push_back(bW2[i]);
           break;
        }
     }
     if ( tCandidates.size() > 1 ) break;

  }
  */
 
  return findt;
}
 
std::vector<iReco> TtSemiEventSolution::recoSemiLeptonicTtEvent(int topo, std::vector<pat::Jet> theWJets,
                            std::vector<pat::Jet> thebJets, std::vector<const reco::Candidate*> theMuons, 
                            std::vector<const reco::Candidate*> theElectrons, 
                            Handle<std::vector<pat::MET> >  mets, HTOP9* histo9 ) {

   std::vector<iReco> TtCollection;
   if (topo != 1 && topo != -1) return TtCollection;

   // get the leptons 
   bool isMuon     = ( theMuons.size() ==1 && theElectrons.size() ==0 ) ? true : false ;
   bool isElectron = ( theMuons.size() ==0 && theElectrons.size() ==1 ) ? true : false ;
   std::vector<const reco::Candidate*> mclepton = ( isMuon && !isElectron  ) ? theMuons : theElectrons ;
   
   int type = 1;
   if (topo == 1 ) type = 0;
   // reco hadronic W 
   std::vector<iReco> hadWs;
   bool hadronic = recoW( type ,theWJets, hadWs, histo9 );
   // reco leptonic W
   std::vector<iReco> lepWs;
   bool leptonic = recoW( type ,mclepton, mets, lepWs, histo9 );

   if ( hadronic && leptonic ) {
      std::vector<iReco> lepTs;
      bool findlepT = recoTop( lepWs ,thebJets, lepTs );
      std::vector<iReco> hadTs;
      bool findhadT = recoTop( hadWs ,thebJets, hadTs );

      std::pair<int, int> ttpair(-1, -1);
      if (findlepT && findhadT ) {
         double dmtt = 999. ;
         for (size_t i=0; i< lepTs.size(); i++ ) {
             for (size_t j=0; j< hadTs.size(); j++) {
                 if ( (lepTs[i].dm + hadTs[j].dm) > dmtt ) continue;
                 if ( lepTs[i].from.second == hadTs[j].from.second ) continue;
                 dmtt = lepTs[i].dm + hadTs[j].dm ;
                 ttpair = make_pair(i,j);
            }
         }
      }
      if ( ttpair.first != -1 && ttpair.second != -1 ) {
         TtCollection.push_back( lepTs[ ttpair.first ] );
         TtCollection.push_back( hadTs[ ttpair.second ] );
      }
   }
   /*
   bool findtops = false ;
   if ( hadronic && leptonic ) {
      std::vector<iReco> wCandidates;
      std::vector<iReco> Tt_temp;
      for (size_t i=0; i< lepWs.size(); i++) {
         double dm_had = 999.;
         for (size_t j=0; j< hadWs.size(); j++) {
            wCandidates.push_back(lepWs[i]);
            wCandidates.push_back(hadWs[j]);
            findtops = recoTop( wCandidates ,thebJets, Tt_temp );
	    wCandidates.clear();
            if ( !findtops ) continue;
            // pick the best hadronic top
            if (Tt_temp[1].dm < dm_had) TtCollection = Tt_temp ;
         } 
      }
   }
   */
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


void TtSemiEventSolution::dmSortRecoObjects( std::vector<iReco>& objCandidates ) {

   sort(objCandidates.begin(), objCandidates.end(), dmIncrease );

}
