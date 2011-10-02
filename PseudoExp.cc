#include "PseudoExp.h"

PseudoExp::PseudoExp(){

  fitInput = new MassAnaInput();
  fitInput->GetParameters( "MuonCuts", &muonCuts );
  fitInput->GetParameters( "JetCuts",  &jetCuts );

}

PseudoExp::~PseudoExp(){

  delete fitInput ;

}

bool ItIncreasing( int s1, int s2) { return ( s1 < s2 ); }
bool EtDescending( TLorentzVector v1, TLorentzVector v2 ) { return ( v1.Et() > v2.Et() ) ;  }

vector<int> PseudoExp::GetEnsemble( string fileName, TString treeName, double pMean, int RandomSeed  ) {

  // get files and trees
  TFile* file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName, treeName, file );

  // event solution tree
  //int evtId1 ;
  //tr1->SetBranchAddress( "evtId" , &evtId1 );
  Int_t tsz = tr1->GetEntries();
  cout<<" All Evt Size = "<< tsz << endl;

  // set up the random number function
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( RandomSeed );  

  int NEvts = tRan->Poisson( pMean );
  cout<<" N of Evts needed = "<< NEvts <<endl ;
  
  vector<int> ensembles ;
  for (int k=0; k < NEvts; k++) {
      int selEntry = tRan->Integer( tsz );

      // double check the selected entries, in case any repeat selection
      bool repeat = false ;
      for (size_t j = 0; j< ensembles.size(); j++ ) {
          repeat = ( ensembles[j] ==  selEntry ) ? true : false ;
          if ( repeat ) break ;
      }
      if ( repeat ) {
         k-- ;
         continue;
      }

      ensembles.push_back( selEntry );
      cout<<k<<"  == selEntry : " << selEntry <<endl;

  }
  sort( ensembles.begin(), ensembles.end(), ItIncreasing ) ;
  return ensembles ;

}

vector<int> PseudoExp::EventShuffle( int theSize, int RandomSeed ) {

  const int tsz = theSize;

  // set up the random number function
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( RandomSeed );  

  vector<int> evtpool;
  for(int i=0; i<tsz; i++) {
     evtpool.push_back( i ) ;
  }
  for(int i=tsz-1; i>1; i--) {
     int k = tRan->Integer( i );
     swap( evtpool[k], evtpool[i]);//swaps the randomly picks character with n
  }

  delete tRan ;
  return evtpool ;

}

void PseudoExp::JetEtReSort( vector<TLorentzVector>& vs ) {

  // re-sort the Jets
  vector<TLorentzVector> newVs ;
  for (size_t i=0; i< vs.size()- 2; i++ ) {
      newVs.push_back( vs[i] );
  }

  sort( newVs.begin(), newVs.end(), EtDescending ) ;
  newVs.push_back( vs[ vs.size() -2 ] );
  newVs.push_back( vs[ vs.size() -1 ] );

  vs.clear();
  for (size_t i=0; i< newVs.size(); i++ ) {
      vs.push_back( newVs[i] ) ;
  }

}


void PseudoExp::PhaseSmearing( vector<TLorentzVector>& vs, int RandomSeed, double jes, bool ReMET ) {

  //cout<<" Start Smearing ~~~ "<<endl;
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( RandomSeed );  
  double rf = tRan->Uniform( -3.1415, 3.1415 ) ;
  double bz = tRan->Uniform( -0.2, 0.2 ) ;
  //double bx = tRan->Gaus( -0.1, 0.1 ) ;
  //double by = tRan->Gaus( -0.1, 0.1 ) ;
  TVector3 bv( 0, 0, bz ) ;
  //cout<<" bv  ("<<bx<<","<<by<<","<<bz<<")"<<"  rotate : "<< rf <<endl;

  //const size_t sz = vs.size() ;
  //double sc[sz] ;
  //tRan->RndmArray( sz, sc ) ;

  double npx = 0.;
  double npy = 0.;
  for (size_t i=0; i< vs.size()- 2; i++ ) {
      //cout<<" "<< i <<"*  vs :"<<vs[i].Pt()<<" / "<< vs[i].Pz()<<endl;
      //mA = mA + vs[i] ;
      //if ( i < 3 ) m3A = m3A + vs[i] ;
      int nS = 0 ;
      //double minPt = ( vs[i].Pt() >= jetCuts[0]  ) ?  jetCuts[0] : 20. ;

      TLorentzVector vs1( vs[i].Px(), vs[i].Py(), vs[i].Pz(), vs[i].E() );
      if ( vs1.Pt() == 0 && vs1.Pz() == 0 ) continue ;

      while( nS == 0 || vs1.Pt() < jetCuts[0] ) {
          vs1 = TLorentzVector( vs[i].Px(), vs[i].Py(), vs[i].Pz(), vs[i].E() ) ;
          tRan->SetSeed( i*17 + nS + RandomSeed );  // just change the random number seed
          double sc = tRan->Rndm();
          double scl = ( jes == 0 ) ? (sc-0.5)/10. : tRan->Gaus( 0, jes ) ;
          if ( nS == 0 ) { 
             vs1.RotateZ( rf );
             vs1 = vs1*( 1. +  scl );
             vs1.Boost( -bv );
          } else  {
             vs1.RotateZ( rf );
          } 
          nS++ ;
      }
      vs[i] = vs1 ;
      npx = npx - vs[i].Px() ;
      npy = npy - vs[i].Py() ;
      //cout<<" -> vs :"<<vs[i].Pt()<<" / "<< vs[i].Pz() <<endl;
      //mB = mB + vs[i] ;
      //if ( i < 3 ) m3B = m3B + vs[i] ;
  }

  // smearing the muon and MET
  //double sc0 = tRan->Rndm();
  //double scl0 = (sc0-0.5)/10. ;
  double scl0 = tRan->Gaus( 0, 0.01 ) ;
  int im = vs.size() - 2 ;
  int in = vs.size() - 1 ;
  double dphi = vs[ im ].DeltaPhi( vs[ in ] ) ;
  double Mt2  = 2.* vs[im].Pt() * vs[in].Pt() *( 1- cos(dphi) );
  
  for (int i= im ; i < in+1 ; i++ ) {
      //cout<<" "<< i <<"*  vs :"<<vs[i].Pt()<<" / "<< vs[i].Pz()<<" M = "<< vs[i].M() <<endl;

      double theM = vs[i].M() ;
      TLorentzVector vs1( vs[i].Px(), vs[i].Py(), vs[i].Pz(), vs[i].E() );
      if ( vs1.Pt() == 0 && vs1.Pz() == 0 ) continue ;

      vs1.RotateZ( rf );
      vs1.Boost( -bv );
      if ( i == im ) { 
         vs1 = vs1*( 1. +  scl0 );
         if ( vs1.Pt() < muonCuts[0] || fabs(vs1.Eta()) > muonCuts[1] ) vs1 = vs[i] ;
      }
      if ( i == in ) {
         double sc = tRan->Rndm();
         scl0 = ( jes == 0 ) ? (sc-0.5)/10. : tRan->Gaus( 0,  jes ) ;
         dphi = vs[ im ].DeltaPhi( vs1 ) ;
	 double new_met =  Mt2 / ( 2.*vs[ im ].Pt()*( 1- cos(dphi) ) );
	 double msc = new_met*(1.+scl0) / vs1.Pt() ;
	 //cout <<" new MET = "<< new_met <<" scale = "<< msc <<endl ;
	 double theE = msc*sqrt( vs1.Pt()*vs1.Pt() + theM*theM ) ;
	 vs1 = TLorentzVector( vs1.Px()* msc, vs1.Py()*msc, 0, theE ) ;
         if ( vs1.Pt() < jetCuts[2] ) vs1 = vs1 * ( jetCuts[2]/vs1.Pt()  )  ;
      }
      vs[i] = vs1 ;
      npx = npx - vs[i].Px() ;
      npy = npy - vs[i].Py() ;
      //cout<<" -> vs :"<<vs[i].Pt()<<" / "<< vs[i].Pz() <<" M = "<< vs[i].M() <<endl;
  }

  // Re-do the MET, if desire ...
  if ( ReMET ) { 
     // smearing the phi 
     TLorentzVector newMET( npx, npy, 0., sqrt( (npx*npx) + (npy*npy) ) ) ;
     double s_phi = 0.69*pow( newMET.Pt(), -0.64 ) ;
     double sf = tRan->Gaus( 0., s_phi ) ;
     newMET.RotateZ( sf ) ;
     // smearing the MET 
     double s_met = 0.7*sqrt( (5.6*5.6)/(newMET.Pt()*newMET.Pt() ) + ( 1.25*1.25/newMET.Pt() ) + (0.033*0.033) );
     double sE = tRan->Gaus( 0., s_met ) ;
     newMET = newMET*( 1. + sE ) ;
     vs.push_back( newMET );
  }
 
  delete tRan ;

}

