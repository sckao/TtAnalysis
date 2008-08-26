// -*- C++ -*-
//
// Package:    TtEfficiency
// Class:      TtEfficiency
// 
/**\class TtEfficiency TtEfficiency.cc PhysicsTools/TtAnalysis/src/TtEfficiency.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Thur Aug 16 2008
//
//


// system include files
#include <memory>

// user include files
#include "TtEfficiency.h"
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
//TtEfficiency::TtEfficiency(const edm::ParameterSet& iConfig)
TtEfficiency::TtEfficiency()
{
   //now do what ever initialization is needed

}


TtEfficiency::~TtEfficiency()
{
   // do anything here that needs to be done at desctruction time
}

//
// member functions
//

// ------------ method called to for each event  ------------
void TtEfficiency::EventEfficiency( int topo , bool pass , HTOP9* histo9 ) {

    // hadronic channel
    if      (topo == 0 &&  pass) { histo9->Fill9f( 0.  ); }
    else if (topo == 0 && !pass) { histo9->Fill9f( 0.5 ); } 
    // semi-leptonic channel
    else if (topo == 1 &&  pass) { histo9->Fill9f( 1.  ); }
    else if (topo == 1 && !pass) { histo9->Fill9f( 1.5 ); }
    // di-leptonic channel
    else if (topo == 2 &&  pass) { histo9->Fill9f( 2.  ); }
    else if (topo == 2 && !pass) { histo9->Fill9f( 2.5 ); }
    else                         { histo9->Fill9f( -1. ); }

}

void TtEfficiency::EventShape( int topo, size_t isoMu, size_t isoE , size_t nBJ, size_t nWJ, HTOP9* histo9 ){

    int  Channel = 0;
    int nBJets =  static_cast<int>(nBJ) ;
    int nWJets =  static_cast<int>(nWJ) ;
    if( isoE  > 1 && isoMu == 0 )  Channel = 1;
    if( isoE == 0 && isoMu  > 1 )  Channel = 2;
    if( isoE  > 0 && isoMu  > 0 )  Channel = 3;
    if( isoE == 0 && isoMu == 1 )  Channel = 5; 
    if( isoE == 1 && isoMu == 0 )  Channel = 7; 
    if( isoE == 0 && isoMu == 0 )  Channel = 9;

    int id = (topo*10) + Channel ;

    if( topo > 2) id = -1 ;

    histo9->Fill9g( id, nBJets, nWJets );
}

void TtEfficiency::IsoLeptonEfficiency(std::vector<const reco::Candidate*> isolep, std::vector<const reco::Candidate*> mclep, HTOP9* histo9 ){ 

 int mc = 0; 
 int rc = 0; 
 for(size_t i=0; i< mclep.size(); i++ ){
    LorentzVector mcP = mclep[i]->p4();
    int pid1 = mclep[i]->pdgId();
    if (abs(pid1) == 11) {
       mc = 1;
       rc = 2;
    }
    if (abs(pid1) == 13) {
       mc = 3;
       rc = 4;
    }
    histo9->Fill9i( mc ); 

    for (size_t j=0; j< isolep.size(); j++  ) {
        LorentzVector isoP = isolep[j]->p4();
        int pid2 = isolep[j]->pdgId();
        if (pid1 == pid2 && mcP == isoP ) {
           histo9->Fill9i( rc ); 
        }
    }
 } 



}


void TtEfficiency::JetEfficiency(std::vector<pat::Jet> recojets, std::vector<pat::Jet> mcjets, HTOP9* histo9 ){ 

 int mc = 0; 
 int rc = 0; 
 for(size_t i=0; i< mcjets.size(); i++ ){
    LorentzVector mcP = mcjets[i].p4();
    int pid1 = mcjets[i].partonFlavour();
    if ( abs(pid1) == 5 ) {
       mc = 5;
       rc = 6;
    }
    if ( abs(pid1)  < 5 ) {
       mc = 7;
       rc = 8;
    }
    //cout <<" flavor1:"<<pid1<<"  p4="<<mcP<<endl;
    histo9->Fill9i( mc ); 

    for (size_t j=0; j< recojets.size(); j++  ) {
        LorentzVector recoP = recojets[j].p4();
        int pid2 = recojets[j].partonFlavour();
        if (pid1 == pid2 && mcP == recoP ) { 
           histo9->Fill9i( rc ); 
        }
        //cout<<" "<<endl;
    }
 } 

}
//define this as a plug-in
