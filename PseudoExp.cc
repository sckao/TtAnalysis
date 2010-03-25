#include "PseudoExp.h"

PseudoExp::PseudoExp( double massL, double massH ){

  fitInput = new MassAnaInput( "had", massL, massH );

}

PseudoExp::~PseudoExp(){

  delete fitInput ;

}

bool IdIncreasing( pair<int,int> s1, pair<int,int> s2) { return ( s1.first < s2.first ); }
bool ItIncreasing( int s1, int s2) { return ( s1 < s2 ); }

vector< pair<int,int> > PseudoExp::GetEnsemble( string fileName, double pMean, int RandomSeed  ) {

  // get files and trees
  TFile* file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName, "solTt", file );

  // event solution tree
  int evtId1, entId1, iniId1, entSz1;
  tr1->SetBranchAddress( "evtId" , &evtId1 );
  tr1->SetBranchAddress( "entId" , &entId1 );
  tr1->SetBranchAddress( "iniId" , &iniId1 );
  tr1->SetBranchAddress( "entSz" , &entSz1 );

  // accumulate all the events
  vector< pair<int,int> > evts ;
  int evtId0 = 0;
  for ( int j= 0 ; j< tr1->GetEntries() ; j++ ) {

      tr1->GetEntry(j);
      if ( evtId0 != evtId1 ) {
         pair<int,int> evtIdx = make_pair( evtId1, iniId1 );
         evts.push_back( evtIdx ) ;
         //cout<<" EvtId list = "<< evtId1 <<endl;
         evtId0 = evtId1;
      }
  }
  cout<<" All Evt Size = "<< evts.size() << endl;

  TRandom* tRan = new TRandom();
  tRan->SetSeed( RandomSeed );  

  int NEvts = tRan->Poisson( pMean );
  cout<<" N of Evts needed = "<< NEvts <<endl ;
  
  vector< pair<int,int> > ensembles ;
  for (int k=0; k < NEvts; k++) {
      int selEntry = tRan->Integer( evts.size() );
      int selId = evts[ selEntry ].first ;
      int selIt = evts[ selEntry ].second ;
      // double check the selected entries, in case any repeat selection
      bool repeat = false ;
      for (size_t j = 0; j< ensembles.size(); j++ ) {
          repeat = ( ensembles[j].first ==  selId ) ? true : false ;
          if ( repeat ) break ;
      }
      if ( repeat ) {
         k-- ;
         continue;
      }

      ensembles.push_back( evts[ selEntry ] );
      cout<<k<<"  == selEntry : " << selEntry <<" sel ID: "<< selId <<" iniId: "<< selIt <<endl;

  }
  sort( ensembles.begin(), ensembles.end(), IdIncreasing ) ;

  return ensembles ;

}

vector<int> PseudoExp::GetEnsemble( string fileName, TString treeName, double pMean, int RandomSeed  ) {

  // get files and trees
  TFile* file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName, treeName, file );

  // event solution tree
  int evtId1 ;
  tr1->SetBranchAddress( "evtId" , &evtId1 );
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
