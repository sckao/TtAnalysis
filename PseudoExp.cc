#include "PseudoExp.h"

PseudoExp::PseudoExp( double massL, double massH ){

  fitInput = new MassAnaInput( "had", massL, massH );

}

PseudoExp::~PseudoExp(){

  delete fitInput ;

}

bool ItIncreasing( int s1, int s2) { return ( s1 < s2 ); }

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
  TRandom* tRan = new TRandom();
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
  TRandom* tRan = new TRandom();
  tRan->SetSeed( RandomSeed );  

  vector<int> evtpool;
  for(int i=0; i<tsz; i++) {
     evtpool.push_back( i ) ;
  }
  for(int i=tsz-1; i>1; i--) {
     int k = tRan->Integer( i );
     swap( evtpool[k], evtpool[i]);//swaps the randomly picks character with n
  }
  return evtpool ;

}

void PseudoExp::PhaseSmearing( vector<TLorentzVector>& vs, int RandomSeed  ) {

  TRandom* tRan = new TRandom();
  tRan->SetSeed( RandomSeed );  
  double rf = tRan->Uniform( -3.1415, 3.1415 ) ;
  double bz = tRan->Uniform( -0.5, 0.5 ) ;
  double bx = tRan->Gaus( 0., 0.15 ) ;
  double by = tRan->Gaus( 0., 0.15 ) ;
  const size_t sz = vs.size() ;
  double sc[sz] ;
  tRan->RndmArray( sz, sc ) ;
  if ( bx >=  0.5 ) bz =  0.5 ;
  if ( by >=  0.5 ) by =  0.5 ;
  if ( bx <= -0.5 ) bz = -0.5 ;
  if ( by <= -0.5 ) by = -0.5 ;
  TVector3 bv( bx, by, bz ) ;
  
  //cout<<" bv  ("<<bx<<","<<by<<","<<bz<<")"<<"  rotate : "<< rf <<endl;
  TLorentzVector mA ;
  TLorentzVector mB ;
  /*
  TLorentzVector m3A ;
  TLorentzVector m3B ;
  */
  double npx = 0.;
  double npy = 0.;
  for (int i=0; i< vs.size(); i++ ) {
      //cout<<" *  vs :"<<vs[i].Pt()<<" / "<< vs[i].Pz()<<" scale : "<< 1 + (sc[i] - 0.5 )/5 <<endl;
      mA = mA + vs[i] ;
      //if ( i < 3 ) m3A = m3A + vs[i] ;
      vs[i] = vs[i]*( 1. +  (sc[i]-0.5)/5. );
      vs[i].RotateZ( rf );
      vs[i].Boost( -bv );
      npx = npx - vs[i].Px() ;
      npy = npy - vs[i].Py() ;
      //cout<<" -> vs :"<<vs[i].Pt()<<" / "<< vs[i].Pz() <<endl;
      mB = mB + vs[i] ;
      //if ( i < 3 ) m3B = m3B + vs[i] ;
  }
  // smearing the phi 
  TLorentzVector newMET( npx, npy, 0., sqrt( (npx*npx) + (npy*npy) ) ) ;
  double s_phi = 0.69*pow( newMET.Pt(), -0.64 ) ;
  double sf = tRan->Gaus( 0., s_phi ) ;
  newMET.RotateZ( sf ) ;
  // smearing the MET 
  double s_met = 0.7*sqrt( (5.6*5.6)/(newMET.Pt()*newMET.Pt() ) +  ( 1.25*1.25/newMET.Pt() ) +  ( 0.033*0.033 ) );
  double sE = tRan->Gaus( 0., s_met ) ;
  newMET = newMET*( 1. + sE ) ;
  vs.push_back( newMET );
  //cout<<" mA:"<< mA.M() <<" pt: "<< mA.Pt() <<" pz: "<< mA.Pz() <<endl;
  //cout<<" mB:"<< mB.M() <<" pt: "<< mB.Pt() <<" pz: "<< mB.Pz() <<" sc = "<< mB.M()/ mA.M() <<endl;
  //cout<<" newMET  px = "<< newMET.Px() <<" py =  "<< newMET.Py()<<endl ;
  /*
  cout<<" m3A:"<< m3A.M() <<" pt: "<< m3A.Pt() <<" pz: "<< m3A.Pz() <<endl;
  cout<<" m3B:"<< m3B.M() <<" pt: "<< m3B.Pt() <<" pz: "<< m3B.Pz()<<" sc = "<< m3B.M()/ m3A.M()  <<endl;
  */
}

//void PseudoExp::PhaseSmearing( TLorentzVector& vs, double scale  ) {
//}
