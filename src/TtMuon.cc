// -*- C++ -*-
//
// Package:    TtMuon
// Class:      TtMuon
// 
/**\class TtMuon TtMuon.cc PhysicsTools/TtAnalysis/src/TtMuon.cc

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
#include "TtMuon.h"

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
TtMuon::TtMuon(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muSetup  = iConfig.getParameter<std::vector<double> >("muSetup");
  //muonSrc = iConfig.getParameter<edm::InputTag> ("muonSource");
  tools        = new TtTools();

}


TtMuon::~TtMuon()
{
   delete tools ;
}

//
// member functions
//
//typedef std::pair<double, pat::Jet> ptjet ;
static bool PtDecreasing(const reco::Candidate* s1, const reco::Candidate* s2) { return ( s1->pt() > s2->pt() ); }
//static bool PtDecreasing1(const pat::Muon* s1, const pat::Muon* s2) { return ( s1->pt() > s2->pt() ); }
static bool PtDecreasing2(ttCandidate s1, ttCandidate s2) { return ( s1.p4.pt() > s2.p4.pt() ); }

// ------------ method called to for each event  ------------
void TtMuon::PromptMuonScope( Handle<std::vector<pat::Muon> > patMu, HOBJ3* histo ){

   double thePt = 0 ;
   double theEta = 999 ;
   double theIso = 999 ;
   for (std::vector<pat::Muon>::const_iterator m1 = patMu->begin(); m1!= patMu->end(); m1++) {
       if ( m1->pt() > thePt ) {
          thePt = m1->pt() ;
          theEta = m1->eta() ;
          double emIso = m1->isolationR03().emEt ;
          double hdIso = m1->isolationR03().hadEt ;
	  double tkIso = m1->isolationR03().sumPt ;
	  double RelIso = ( tkIso + emIso + hdIso )/ m1->pt() ;
          theIso = RelIso ;
        }
   }
   histo->Fill_3d( thePt, theEta, theIso, 0., 0., 0. );

}

void TtMuon::PatMuonScope( Handle<std::vector<pat::Muon> > patMu, const edm::Event& iEvent, HOBJ3* histo ){

     // get the beam spot information
     //edm::Handle<reco::BeamSpot> theBeamSpot;
     //iEvent.getByLabel("offlineBeamSpot", theBeamSpot );

     //math::XYZPoint pos_beamspot( 0, 0, 0 );
     //if ( theBeamSpot->isValid() )  pos_beamspot = pos_beamspot( beamSpot.x0(), beamSpot.y0(), beamSpot.z0() );
     //math::XYZPoint pos_beamspot( theBeamSpot->x0(), theBeamSpot->y0(), theBeamSpot->z0() );

     int N_glbMu = 0;
     for (std::vector<pat::Muon>::const_iterator m1 = patMu->begin(); m1!= patMu->end(); m1++) {

         //1. Muon Isolation
         /*
         const reco::IsoDeposit* ecalIso  = (*m1)->ecalIsoDeposit();
	 const reco::IsoDeposit* hcalIso  = (*m1)->hcalIsoDeposit();
	 const reco::IsoDeposit* trackIso = (*m1)->trackIsoDeposit();
	 std::pair<double, int> emR3 = ecalIso->depositAndCountWithin(0.3);
	 std::pair<double, int> hdR3 = hcalIso->depositAndCountWithin(0.3);
	 std::pair<double, int> tkR3 = trackIso->depositAndCountWithin(0.3);
         double RelIso =  (emR3.first + hdR3.first + tkR3.first)/(*m1)->pt()  ;
         */
         double emIso = m1->isolationR03().emEt ;
	 double hdIso = m1->isolationR03().hadEt ;
	 double tkIso = m1->isolationR03().sumPt ;
	 double RelIso = ( tkIso + emIso + hdIso )/ m1->pt() ;

         //histo->Fill_3d( (*m1)->pt(), (*m1)->eta(), RelIso, emR3.first, hdR3.first, tkR3.first );
         histo->Fill_3d( m1->pt(), m1->eta(), RelIso, emIso, hdIso, tkIso );

         // track quality
         int nSeg = 0 ;
         double stdX2 = -1 ;
         if ( m1->isStandAloneMuon() ) {
            reco::TrackRef outTrack = m1->outerTrack();
            nSeg  = outTrack->numberOfValidHits();
            stdX2 = outTrack->normalizedChi2();
            histo->Fill_3e( nSeg, stdX2 );
         }

         int nHit = 0 ;
         double trkX2 = -1 ;
         if ( m1->isTrackerMuon() ) {
            reco::TrackRef inTrack = m1->innerTrack();
            nHit  = inTrack->numberOfValidHits();
            trkX2 = inTrack->normalizedChi2();
            // beam spot study , according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideFindingBeamSpot
            /*
            double d0    = -1.* m1->innerTrack()->dxy( pos_beamspot );
            double beamWdthX = theBeamSpot->BeamWidthX() ;
            double beamWdthY = theBeamSpot->BeamWidthY() ;
            double d0Err = sqrt( inTrack->d0Error()*inTrack->d0Error() + 0.5*beamWdthX*beamWdthX + 0.5* beamWdthY*beamWdthY );
            */
            double d0    = m1->innerTrack()->d0();
            double d0Err = inTrack->d0Error()  ;
            double ipSig = d0 / d0Err ;
            histo->Fill_3f( nHit, trkX2, ipSig );
         }

         int aHit = 0 ;
         double glbX2 = -1 ;
         if ( m1->isGlobalMuon() ) {
            reco::TrackRef inTrack = m1->innerTrack();
            reco::TrackRef glbTrack = m1->globalTrack();
            aHit  = glbTrack->numberOfValidHits();
            glbX2 = glbTrack->normalizedChi2();

            double Pt_Res = (glbTrack->pt()/inTrack->pt()) - 1  ;
            histo->Fill_3g( aHit, glbX2, glbTrack->pt(), m1->eta(), RelIso, inTrack->pt(), Pt_Res );
            N_glbMu++;
         }

     }
     int muSize = static_cast<int>( patMu->size() );
     histo->Fill_3c( muSize, N_glbMu );

}

std::vector<const reco::Candidate*> TtMuon::IsoMuonSelection( Handle<std::vector<pat::Muon> > patMu) {

 std::vector<const reco::Candidate*> isoMuons;
 isoMuons.clear();
 bool isolated = false ;
 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     if ( !(it->isGlobalMuon()) ) continue ;
     isolated = IsoMuonID( *it, muSetup[2] );    

     if ( !(it->isTrackerMuon()) ) continue ;
     reco::TrackRef inTrack = it->innerTrack();
     int    nHit  = inTrack->numberOfValidHits();
     double d0    = inTrack->d0();
     double d0Err = inTrack->d0Error()  ;
     double ipSig = d0 / d0Err ;
     

     bool pass = true ;
     if ( it->pt()         < muSetup[0]  ) pass = false ;  
     if ( fabs(it->eta())  > muSetup[1]  ) pass = false ;  
     if ( nHit             < muSetup[3]  ) pass = false ;  
     if ( it->normChi2()   > muSetup[4]  ) pass = false ;
     if ( ipSig            > muSetup[5]  ) pass = false ;

     if ( isolated && pass ) isoMuons.push_back( &*it );
 }
 if ( isoMuons.size() > 1 ) sort( isoMuons.begin(), isoMuons.end(), PtDecreasing );
 return isoMuons ;

}

std::vector<ttCandidate> TtMuon::IsoMuonSelection1( Handle<std::vector<pat::Muon> > patMu,  Handle<std::vector<pat::Jet> > patJet, Handle<reco::BeamSpot> bSpot_, double PVz, std::vector<ttCandidate>& vetoInfo ) {

 std::vector<ttCandidate> isoMuons;
 isoMuons.clear();
 for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {

     if ( !(it->isGlobalMuon()) ) continue ;
     double relIso = IsoMuonID( *it );    

     bool pass = true ;
     
     reco::TrackRef glbTrack = it->globalTrack();
     double nChi2 = glbTrack->normalizedChi2() ;

     reco::MuonEnergy muE = it->calEnergy() ;
     double muCalE = muE.em + muE.had ;

     // beam Spot information
     reco::BeamSpot bSpot = *bSpot_ ;
     //double beamWdthX = bSpot.BeamWidthX() ;
     //double beamWdthY = bSpot.BeamWidthY() ;

     // dB = d0(Bsp) if set process.patMuons.usePV = False 
     double d0 = it->dB() ;
     int  nHit = 0 ; 
     int  nSiLayer = 0 ;
     if ( it->isTrackerMuon() ) {
        reco::TrackRef inTrack = it->innerTrack();
	nHit  = inTrack->numberOfValidHits();
	// d0    = inTrack->d0();
	// double d0Err = inTrack->d0Error()  ;
	d0    = -1.*inTrack->dxy( bSpot.position() ) ;
	//double d0Err = sqrt( inTrack->d0Error()*inTrack->d0Error() + 0.5*beamWdthX*beamWdthX + 0.5* beamWdthY*beamWdthY );
	//double ipSig = d0 / d0Err ;
        nSiLayer = inTrack->hitPattern().pixelLayersWithMeasurement() ;
     } else {
       pass = false ;
     }

     double dZ = fabs( it->vertex().z() - PVz ) ;

     if ( !it->muonID("GlobalMuonPromptTight") ) pass = false ;
     if ( it->pt()         <= muSetup[0]  ) pass = false ;  
     if ( fabs(it->eta())  >= muSetup[1]  ) pass = false ;  
     if ( relIso           >= muSetup[2]  ) pass = false ;  
     if ( nHit             <= muSetup[3]  ) pass = false ;  
     if ( nChi2            >= muSetup[4]  ) pass = false ;
     if ( fabs( d0 )       >= muSetup[5]  ) pass = false ;
     if ( nSiLayer         <           1  ) pass = false ;
     if ( it->numberOfMatches() <      2  ) pass = false ;
     if ( dZ               >=          1  ) pass = false ;
     if ( glbTrack->hitPattern().numberOfValidMuonHits() < 1 ) pass = false ;

     //double mindR = 999.;
     /*    
     for (std::vector<pat::Jet>::const_iterator j1 = patJet->begin(); j1 != patJet->end(); j1++) {

         if ( !pass ) break ;
         
         //if ( !j1->isCaloJet() ) continue ;
         //if ( j1->pt() < 30. || fabs(j1->eta()) > 2.4 ||  j1->emEnergyFraction() < 0.01 || j1->jetID().n90Hits <= 1 ) continue; 
         //if ( j1->jetID().fHPD >= 0.98 ) continue;
         
         if ( !j1->isPFJet() ) continue;
         if ( j1->pt() < 25. || fabs(j1->eta()) > 2.4 ) continue ;
	 int    NConstituents = j1->numberOfDaughters() ;
	 double CEF = j1->chargedEmEnergyFraction() ;
	 double NHF = j1->neutralHadronEnergyFraction() ;
	 double NEF = j1->neutralEmEnergyFraction() ;
	 double CHF = j1->chargedHadronEnergyFraction() ;
	 int    NCH = j1->chargedMultiplicity() ;
         if ( NConstituents <= 1 )                        continue ;
         if ( CEF >= 0.99 || NHF >= 0.99 || NEF >= 0.99 ) continue ;
         if ( fabs( j1->eta()) < 2.4 && CHF <= 0 )        continue ;
         if ( fabs( j1->eta()) < 2.4 && NCH <= 0 )        continue ;

         double dR = deltaR( it->p4(), j1->p4() ) ;
         if ( dR <= 0.3 )  pass = false ;
         //if ( dR < mindR )  mindR = dR ;
     }
     */
     if ( pass ) {
        ttCandidate ttMuon ;
        ttMuon.p4 = it->p4() ;
        ttMuon.eta = it->eta() ;
        ttMuon.iso =  relIso ;
        ttMuon.cuts[0] = d0 ;
        ttMuon.cuts[1] = nChi2 ;
        ttMuon.cuts[2] = muCalE  ;
        ttMuon.charge = it->charge() ;
        ttMuon.nHits = nHit ;
        ttMuon.pdgId = it->pdgId() ;
        isoMuons.push_back( ttMuon );
     } else {
        if ( it->pt() < 10.  ) continue ;  
	if ( fabs(it->eta())  > 2.5  ) continue ;  
	if ( relIso           > 0.2  ) continue ;  
	ttCandidate vetoMuon ;
	vetoMuon.p4 = it->p4() ;
	vetoMuon.eta = it->eta() ;
	vetoMuon.iso =  relIso ;
	vetoMuon.cuts[0] = it->dB() ;
	vetoMuon.cuts[1] = nChi2 ;
	vetoMuon.cuts[2] = muCalE  ;
	vetoMuon.charge = it->charge() ;
	vetoMuon.nHits = it->numberOfValidHits() ;
	vetoMuon.pdgId = 13 ;
	vetoInfo.push_back( vetoMuon );
     }

 }
 if ( isoMuons.size() > 1 ) sort( isoMuons.begin(), isoMuons.end(), PtDecreasing2 );
 return isoMuons ;

}

bool TtMuon::IsoMuonID( pat::Muon muon, double isoCut ) {

    bool isolation = false;
    //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
    /*
    const reco::IsoDeposit* ecalIso  = muon.ecalIsoDeposit();
    const reco::IsoDeposit* hcalIso  = muon.hcalIsoDeposit();
    const reco::IsoDeposit* trackIso = muon.trackIsoDeposit();
    std::pair<double, int> emR3 = ecalIso->depositAndCountWithin(0.3);
    std::pair<double, int> hdR3 = hcalIso->depositAndCountWithin(0.3);
    std::pair<double, int> tkR3 = trackIso->depositAndCountWithin(0.3);
    //double RelIso =  (emR3.first + hdR3.first + tkR3.first)/muon.pt()  ;
    */
    double emIso = muon.isolationR03().emEt ;
    double hdIso = muon.isolationR03().hadEt ;
    double tkIso = muon.isolationR03().sumPt ;
    if ( isoCut > 0.3 ) {
       emIso = muon.isolationR05().emEt ;
       hdIso = muon.isolationR05().hadEt ;
       tkIso = muon.isolationR05().sumPt ;
    }
    //double RelIso = muon.pt() / ( muon.pt() + emR3.first + hdR3.first + tkR3.first ) ;
    //if ( RelIso >= isoCut ) isolation = true;
    // new Iso definition
    double RelIso = ( tkIso + emIso + hdIso )/ muon.pt() ;
    // Isolation Cut
    if ( RelIso < isoCut ) isolation = true;

    return isolation;
}

double TtMuon::IsoMuonID( pat::Muon muon ) {

    //const reco::IsoDeposit* AllIso  = it->isoDeposit( pat::TrackerIso );
    /*
    const reco::IsoDeposit* ecalIso  = muon.ecalIsoDeposit();
    const reco::IsoDeposit* hcalIso  = muon.hcalIsoDeposit();
    const reco::IsoDeposit* trackIso = muon.trackIsoDeposit();
    std::pair<double, int> emR3 = ecalIso->depositAndCountWithin(0.3);
    std::pair<double, int> hdR3 = hcalIso->depositAndCountWithin(0.3);
    std::pair<double, int> tkR3 = trackIso->depositAndCountWithin(0.3);
    double RelIso =  (emR3.first + hdR3.first + tkR3.first)/muon.pt()  ;

    double emIso = muon.ecalIso() ;
    double hdIso = muon.hcalIso() ;
    double tkIso = muon.trackIso() ;
    */
    double emIso = muon.isolationR03().emEt ;
    double hdIso = muon.isolationR03().hadEt ;
    double tkIso = muon.isolationR03().sumPt ;

    double RelIso = ( tkIso + emIso + hdIso )/ muon.pt() ;

    return RelIso;
}

std::vector<double> TtMuon::MuonEtCorrection( Handle<std::vector<pat::Muon> > mu ) {

     std::vector<double> ptcorr;
     ptcorr.clear();
     double pxy[2] ={0.0};
     for (std::vector<pat::Muon>::const_iterator u1 = mu->begin(); u1 != mu->end(); u1++)
     {
         if ( !(*u1).isGlobalMuon() ) continue;
	 const reco::IsoDeposit* caloE = u1->ecalIsoDeposit(); 
	 const reco::IsoDeposit* caloH = u1->hcalIsoDeposit(); 
         // energy deposite in ECal and HCal
	 double Uem = caloE->candEnergy() ;
	 double Uhd = caloH->candEnergy() ;
         double ex = (Uem+Uhd) * sin((*u1).theta()) * cos((*u1).phi()) ;
         double ey = (Uem+Uhd) * sin((*u1).theta()) * sin((*u1).phi()) ;
         double et = sqrt( ex*ex + ey*ey );

         // over-correction
         bool badCorr =  ( et > (*u1).pt() ? true:false ) ;

         if ( ex > 5. || ey > 5. || badCorr ) { 
            pxy[0] += ((*u1).px() - ex) ;
            pxy[1] += ((*u1).py() - ey) ;
         }
     }
     ptcorr.push_back(pxy[0]);
     ptcorr.push_back(pxy[1]);
     return ptcorr ;
}

double TtMuon::PatMuonRelIso( const pat::Muon& patMu, double ConeSize,  double& ecalDepo, double& hcalDepo, double& trkDepo ){

       /*
       reco::MuonIsolation isoR3 = (*m1)->isolationR03();
       float emIso  = isoR3.emVetoEt ;
       float hadIso = isoR3.hadVetoEt + isoR3.hoVetoEt;
       float trkIso = isoR3.trackerVetoPt ;
       double RelIso = ( emIso + hadIso + trkIso ) / (*m1)->pt() ;
       */
       const reco::IsoDeposit* ecalIso  = patMu.ecalIsoDeposit();
       const reco::IsoDeposit* hcalIso  = patMu.hcalIsoDeposit();
       const reco::IsoDeposit* trackIso = patMu.trackIsoDeposit();
       std::pair<double, int> emR3 = ecalIso->depositAndCountWithin( ConeSize );
       std::pair<double, int> hdR3 = hcalIso->depositAndCountWithin( ConeSize );
       std::pair<double, int> tkR3 = trackIso->depositAndCountWithin( ConeSize );
       double RelIso =  (emR3.first + hdR3.first + tkR3.first) / patMu.pt()  ;
       ecalDepo = emR3.first ;
       hcalDepo = hdR3.first ;
       trkDepo = tkR3.first ;
       return RelIso ;
}

/*
void TtMuon::MuonTrigger( Handle<std::vector<pat::Muon> >patMu, Handle <edm::TriggerResults> triggers ) {

   for (std::vector<pat::Muon>::const_iterator it = patMu->begin(); it!= patMu->end(); it++) {
       std::vector<pat::TriggerPrimitive> trigInfo = it->triggerMatches() ;
       cout<<" === Muon trigger ==== "<<endl;
       for(size_t i=0; i< trigInfo.size(); i++) {
          cout<<"    ObjId:"<< trigInfo[i].triggerObjectId() ;
          cout<<"  ObjType:"<< trigInfo[i].triggerObjectType() <<endl;
          cout<<" filter name:"<<trigInfo[i].filterName() <<endl; 
       }
       if ( trigInfo.size() < 1) continue;
       edm::TriggerNames trigNames( *triggers );
       for (size_t i=0; i< triggers->size(); i++ ) {
           string triggered = triggers->accept(i) ? "Yes" : "No" ;
           cout<<" path("<<i<<") accepted ? "<< triggered ;
           cout<<" trigName: "<< trigNames.triggerName(i)<<endl;
       }
   }

}

*/

