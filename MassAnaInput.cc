#include "MassAnaInput.h"

MassAnaInput::MassAnaInput() {

  probName = "probM";
  weighting = true;
 
  string decayType ;
  GetParameters( "DecayType", &decayType );
  GetParameters( "MassLBound", &mL );
  GetParameters( "MassLBound", &mH );

  ch_name = decayType;
  hname = decayType + "TM";
  //n_btag = NBTag;
  GetParameters( "bThreshold", &bTh);
  GetParameters( "n_btag", &n_btag);
  //N_tt = 8720;  // for OctX, MET_Workshop
  neu_str = 0;


}

MassAnaInput::~MassAnaInput(){


}

vector<TTree*> forest2J ;
vector<TTree*> forest4J ;
vector<TTree*> forestXJ ;
vector<TTree*> forestData ;

void MassAnaInput::Initialize( TString* hfolder ) {

   if ( n_btag == -1 ) *hfolder = "WFitter_"+ch_name+"/";
   if ( n_btag == 0 )  *hfolder = "AllTags0_"+ch_name+"/";
   if ( n_btag == 1 )  *hfolder = "AllTags1_"+ch_name+"/";
   if ( n_btag == 2 )  *hfolder = "AllTags2_"+ch_name+"_10/";

}

void MassAnaInput::LinkForests( TString treeName ){

  cout<<" Linking all trees !!!"<<endl;
  vector<string> fNames2J ;
  GetParameters( "2JSamples", &fNames2J );
  for ( size_t i =0 ; i< fNames2J.size(); i++ ) {
      TTree* tr = GetTree( fNames2J[i], treeName ) ;
      forest2J.push_back( tr );
  }

  vector<string> fNames4J ;
  GetParameters( "4JSamples", &fNames4J );
  for ( size_t i =0 ; i< fNames4J.size(); i++ ) {
      TTree* tr = GetTree( fNames4J[i], treeName ) ;
      forest4J.push_back( tr );
  }

  vector<string> fNamesData ;
  GetParameters( "TheData", &fNamesData );
  for ( size_t i =0 ; i< fNamesData.size(); i++ ) {
      TTree* tr = GetTree( fNamesData[i], treeName ) ;
      forestData.push_back( tr );
  }

  vector<string> fNamesXJ ;
  GetParameters( "FakeData", &fNamesXJ );
  for ( size_t i =0 ; i< fNamesXJ.size(); i++ ) {
      TTree* tr = GetTree( fNamesXJ[i], treeName ) ;
      forestXJ.push_back( tr );
  }
}

TTree* MassAnaInput::TreeMap( string fileName ){
 
    vector<string> f0Names ;
    GetParameters( "TheData", &f0Names );
    vector<string> f1Names ;
    GetParameters( "FakeData", &f1Names );
    vector<string> fNames2J ;
    GetParameters( "2JSamples", &fNames2J );
    vector<string> fNames4J ;
    GetParameters( "4JSamples", &fNames4J );

    TTree* theTr = 0;
    for ( int i=0; i< f0Names.size(); i++ ) {
        if ( f0Names[i] == fileName ) theTr = forestData[i] ;
    }
    for ( int i=0; i< f1Names.size(); i++ ) {
        if ( f1Names[i] == fileName ) theTr = forestXJ[i] ;
    }
    for ( int i=0; i< fNames2J.size(); i++ ) {
        if ( fNames2J[i] == fileName ) theTr = forest2J[i] ;
    }
    for ( int i=0; i< fNames4J.size(); i++ ) {
        if ( fNames4J[i] == fileName ) theTr = forest4J[i] ;
    }
    return theTr ;
}

int MassAnaInput::TreeSize( string fileName ){
 
    vector<string> f0Names ;
    GetParameters( "TheData", &f0Names );
    vector<string> f1Names ;
    GetParameters( "FakeData", &f1Names );
    vector<string> fNames2J ;
    GetParameters( "2JSamples", &fNames2J );
    vector<string> fNames4J ;
    GetParameters( "4JSamples", &fNames4J );

    int treeSize = 0;
    for ( int i=0; i< f0Names.size(); i++ ) {
        if ( f0Names[i] == fileName ) treeSize = forestData[i]->GetEntries() ;
    }
    for ( int i=0; i< f1Names.size(); i++ ) {
        if ( f1Names[i] == fileName ) treeSize = forestXJ[i]->GetEntries() ;
    }
    for ( int i=0; i< fNames2J.size(); i++ ) {
        if ( fNames2J[i] == fileName ) treeSize = forest2J[i]->GetEntries() ;
    }
    for ( int i=0; i< fNames4J.size(); i++ ) {
        if ( fNames4J[i] == fileName ) treeSize = forest4J[i]->GetEntries() ;
    }

    return treeSize ;
}

vector<TTree*> MassAnaInput::GetForest( string DataSet, TString treeName ) {

    cout<<"  =>>> getting a forest of "<< treeName <<endl ;
    vector<string> fileList;
    GetParameters( DataSet , &fileList );

    vector<TTree*> forest ;
    for ( size_t i =0 ; i< fileList.size(); i++ ) {
        TTree* tr = GetTree( fileList[i], treeName ) ;
        forest.push_back( tr );
    }
    return forest ;
}

TTree* MassAnaInput::GetTree( string chName, TString treeName, TFile* file  ) {
  
  TTree* tr = 0;

  string filePath ;
  GetParameters( "RootFiles", &filePath );

  TString theFileName ;
  TChain* theChain = new TChain( treeName ) ;

  if ( chName[ chName.size()-1 ] == '+'  ) {
     string ChainName = chName.substr( 0, chName.size()-1 ) + "Chain"  ;
     vector<string> chainlist;
     GetParameters( ChainName, &chainlist );
     cout<<" * fileName+ = "<< ChainName <<endl;
     for ( size_t j=0; j< chainlist.size(); j++) {
         theFileName = filePath + chainlist[j]+".root" ;
         //cout<<" ** fileName = "<< theFileName <<endl;
         theChain->Add( theFileName );
     }
     tr = theChain ;
  } else {
    theFileName = filePath + chName+".root" ;
    cout<<" * fileName = "<< theFileName <<endl;
    if ( file == NULL ) file = TFile::Open( theFileName );
    tr = (TTree*) file->Get( treeName );
    //theChain->Add( theFileName );
    //tr = theChain ;
  }

  return tr ;
}

TTree* MassAnaInput::GetTree( string chName, TString treeName, string fSuffix ) {
  
  TTree* tr = 0;

  TString theFileName ;
  TChain* theChain = new TChain( treeName ) ;

  string ChainName = chName + fSuffix + "Chain"  ;
  vector<string> chainlist;
  GetParameters( ChainName, &chainlist );
  cout<<" * fileName+ = "<< ChainName <<endl;
  for ( size_t j=0; j< chainlist.size(); j++) {
      theFileName = chainlist[j]+".root" ;
      theChain->Add( theFileName );
  }
  tr = theChain ;

  return tr ;
}

// general method
void MassAnaInput::get_h1Obj(TString fname, TString TName, TString BName, TH1D* h1, double theScale, bool weight) {

  TFile* file = TFile::Open( fname );
  TTree*  tr1 = (TTree*) file->Get( TName );

  // retrieve the variables
  double brV;
  tr1->SetBranchAddress( BName ,&brV );
  double prb;
  tr1->SetBranchAddress( probName ,&prb );

  std::vector<bool> blist;
  std::vector<double> bProb;
  ListWithBTag( fname, n_btag, &blist, &bProb );

  // get entries for branch
  Int_t tsz = tr1->GetEntries();
  //cout<<" size of tree = "<< tsz << endl;

  for (int k=0; k<tsz; k++) {
      tr1->GetEntry(k);
      if (  weight ) h1->Fill( brV, bProb[k] );
      if ( !weight ) h1->Fill( brV );
  }

  h1->Scale( theScale ) ;
  file->Close();
}

// method for TChain
void MassAnaInput::get_h1Obj(TChain* tr1, TString BName, TH1D* h1, double theScale, bool weight) {

  // retrieve the variables
  double brV;
  tr1->SetBranchAddress( BName ,&brV );
  double prb;
  tr1->SetBranchAddress( probName ,&prb );

  //std::vector<bool> blist;
  //std::vector<double> bProb;
  //ListWithBTag( fname, n_btag, &blist, &bProb );

  // get entries for branch
  Int_t tsz = tr1->GetEntries();
  //cout<<" size of tree = "<< tsz << endl;

  for (int k=0; k<tsz; k++) {
      tr1->GetEntry(k);
      //if (  weight ) h1->Fill( brV, bProb[k] );
      if ( !weight ) h1->Fill( brV );
  }

  h1->Scale( theScale ) ;

}

void MassAnaInput::getMcMatching( TString mName, TString brName, TH1D* h1, double theScale ) {

  //GetFileName( mName, 1 );
  vector<string> msets;
  GetParameters( "TMassAssumption", &msets );
  vector<string> sglist;
  GetParameters( "Signals", &sglist );

  TString theFileName ;
  for (size_t i=0; i< msets.size(); i++) {
      TString massump = msets[i].substr(0,3) ;
      if ( massump == mName ) theFileName = sglist[i] + ".root" ;
  }
  
  //TFile* file = TFile::Open( theSG );
  TFile* file = TFile::Open( theFileName );
  TTree*  tr0 = (TTree*) file->Get( "mcmTt" );
  TTree*  tr1 = (TTree*) file->Get( "solTt" );

  int evtId0, entId0, iniId0, entSz0, jid0[5] ;
  tr0->SetBranchAddress( "evtId" ,&evtId0 );
  tr0->SetBranchAddress( "entId" ,&entId0 );
  tr0->SetBranchAddress( "entSz" ,&entSz0 );
  tr0->SetBranchAddress( "iniId" ,&iniId0 );
  tr0->SetBranchAddress( "JId"   ,&jid0   );

  int evtId1, entId1, iniId1, entSz1, jid1[5] ;
  tr1->SetBranchAddress( "evtId" ,&evtId1 );
  tr1->SetBranchAddress( "entId" ,&entId1 );
  tr1->SetBranchAddress( "entSz" ,&entSz1 );
  tr1->SetBranchAddress( "iniId" ,&iniId1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );

  double brV, prob ;
  tr1->SetBranchAddress( probName ,&prob );
  tr1->SetBranchAddress( brName ,&brV );

  int startPoint = 0;
  int endPoint   = 0;
  for (int j=0; j< tr0->GetEntries() ; j++) {
      tr0->GetEntry(j);
      startPoint = iniId0 -1 ;
      endPoint   = entId0 -1 ;
      for (int k= startPoint ; k< endPoint ; k++) {
          tr1->GetEntry(k);
          if ( evtId0  != evtId1 )  continue;
          if ( jid0[4] != jid1[4] ) continue;
          if ( jid0[3] != jid1[3] ) continue;
          if ( jid0[2] != jid1[2] ) continue;
          bool isW1 = ( jid0[0] == jid1[0] && jid0[1] == jid1[1]  ) ? true : false ;
          bool isW2 = ( jid0[0] == jid1[1] && jid0[1] == jid1[0]  ) ? true : false ;
          if ( !isW1 && !isW2 ) continue ;
          h1->Fill( brV, prob );
      }
  }
  h1->Scale( theScale ) ;
  //NormalizeComponents( luminosity, N_tt, 1., h1 );
  NormalizeComponents( "tt", h1 );
  file->Close();

}

// Get all hadronic
void MassAnaInput::getHadPermutation( TString fName, TString brName, TH1D* h1, double theScale, std::vector<bool>* blist ) {

  TFile* file = TFile::Open( fName );
  TTree*  tr1 = (TTree*) file->Get( "solTt" );

  int entId1, iniId1, entSz1, jid1[5] ;
  tr1->SetBranchAddress( "entId" ,&entId1 );
  tr1->SetBranchAddress( "entSz" ,&entSz1 );
  tr1->SetBranchAddress( "iniId" ,&iniId1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );

  double brV;
  tr1->SetBranchAddress( brName ,&brV );

  std::vector<int> j0;
  std::vector<int> j1;
  std::vector<int> j2;
  for (int j=0; j< tr1->GetEntries() ; j++) {
      tr1->GetEntry(j);

      if ( blist != NULL ) {
         if ( !(*blist)[j] ) continue;
      } 

      if ( j0.size() == 0 ) {
         j0.push_back( jid1[0] ) ;
         j1.push_back( jid1[1] ) ;
         j2.push_back( jid1[2] ) ;
         h1->Fill( brV );
      }
      int sameM = 0;
      for ( size_t k = 0; k < j0.size() ; k++ ) {
          for (int i= 0; i< 3; i++) {
              if ( jid1[i] == j0[k] ) sameM++ ;
              if ( jid1[i] == j1[k] ) sameM++ ;
              if ( jid1[i] == j2[k] ) sameM++ ;
          }
          if ( sameM == 3 ) break;
          if ( sameM  < 3 ) sameM = 0;
      }
      if ( sameM < 3 ) {
         j0.push_back( jid1[0] ) ;
         j1.push_back( jid1[1] ) ;
         j2.push_back( jid1[2] ) ;
         h1->Fill( brV ); 
      }
      if ( (entId1 + 1) == (iniId1 + entSz1) ) {
         j0.clear() ;
         j1.clear() ;
         j2.clear() ;
      }
  }

  h1->Scale( theScale ) ;
  //NormalizeComponents( luminosity, N_tt, 1., h1 );
  NormalizeComponents( "tt", h1 );
  file->Close();
}

//  check all combination with BTagging 
void MassAnaInput::ListWithBTag( TString fName, int nbtags, std::vector<bool>* blist, std::vector<double>* bProb ) {

  TFile* file = TFile::Open( fName );
  TTree*  tr1 = (TTree*) file->Get( "solTt"  );
  TTree*  tr2 = (TTree*) file->Get( "selJet" );

  int entId1, iniId1, entSz1, jid1[5], nj ;
  tr1->SetBranchAddress( "entId" ,&entId1 );
  tr1->SetBranchAddress( "iniId" ,&iniId1 );
  tr1->SetBranchAddress( "entSz" ,&entSz1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );
  tr1->SetBranchAddress( "nJ"    ,&nj     );
  double prb;
  tr1->SetBranchAddress( probName ,&prb );

  double bDis;
  tr2->SetBranchAddress( "qCut"  ,&bDis   );

  std::vector<int> Bj ;
  int evtCount =  0 ;
  double bNorm = 0;
  for (int j=0; j< tr1->GetEntries() ; j++) {
      tr1->GetEntry(j);

      // count the B jets
      if ( entId1 == iniId1 ) {
         evtCount += nj ;
         for (int i = evtCount-nj ; i < evtCount; i++ ) {
             tr2->GetEntry(i);
             if ( bDis > bTh ) Bj.push_back( i - evtCount + nj );
         }
      }
      // check # of b jets used
      int usedB = 0;
      for ( size_t k = 0; k < Bj.size() ; k++ ) {
          if ( jid1[2] == Bj[k] ) usedB++;
          if ( jid1[3] == Bj[k] ) usedB++;
      }
      // record # of b jets used and calculate bjet used probability
      //if ( bProb != NULL && nbtags > 0 ) { 
      if ( bProb != NULL ) { 
         double totalPrb = usedB*prb ; 
         bProb->push_back( totalPrb );
         bNorm += totalPrb ;
      }
      if ( bProb != NULL && nbtags == -1 )  bProb->push_back( prb );

      // flag the btagging scenario
      if ( usedB == nbtags )  blist->push_back(true) ;
      blist->push_back( false );

      if ( (entId1 + 1) == (iniId1 + entSz1) ) {
          Bj.clear();
          // normalize the btagging use probability
          for (int i=iniId1-1; i<iniId1+entSz1-1; i++ ) {
              if ( bProb == NULL || nbtags == -1 ) break;
              if ( bNorm != 0 ) (*bProb)[i] = (*bProb)[i] / bNorm ;
              if ( bNorm == 0 ) (*bProb)[i] = 0 ;
              //(*bProb)[i] = 1. ;
          }
          bNorm = 0;
      }
      
  }
  file->Close();

}

void MassAnaInput::getLepPermutation( TString fName, TString brName, TH1D* h1, double theScale, std::vector<bool>* blist ) {

  TFile* file = TFile::Open( fName );
  TTree*  tr1 = (TTree*) file->Get( "solTt" );

  int entId1, iniId1, entSz1, jid1[5] ;
  tr1->SetBranchAddress( "entId" ,&entId1 );
  tr1->SetBranchAddress( "entSz" ,&entSz1 );
  tr1->SetBranchAddress( "iniId" ,&iniId1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );

  double brV;
  tr1->SetBranchAddress( brName ,&brV );

  std::vector<int> j3;
  std::vector<int> j4;
  for (int j=0; j< tr1->GetEntries() ; j++) {
      tr1->GetEntry(j);

      if ( blist != NULL ) {
         if ( !(*blist)[j] ) continue;
      }
 
      if ( j3.size() == 0 ) {
         j3.push_back( jid1[3] ) ;
         j4.push_back( jid1[4] ) ;
         h1->Fill( brV );
      }
      int sameM = 0;
      for ( size_t k = 0; k < j3.size() ; k++ ) {
          for (int i= 3; i< 5; i++) {
              if ( jid1[i] == j3[k] ) sameM++ ;
              if ( jid1[i] == j4[k] ) sameM++ ;
          }
      }
      if ( sameM < 2 ) {
         j3.push_back( jid1[3] ) ;
         j4.push_back( jid1[4] ) ;
         h1->Fill( brV ); 
      }
      if ( (entId1 + 1) == (iniId1 + entSz1) ) {
         j3.clear() ;
         j4.clear() ;
      }
  }

  h1->Scale( theScale ) ;
  //NormalizeComponents( luminosity, N_tt, 1., h1 );
  NormalizeComponents( "tt", h1 );
  file->Close();

}


void MassAnaInput::getMostProb( TString mName, TString brName, TH1D* h1,  double theScale ) {

  //GetFileName( mName, 1 );

  vector<string> msets;
  GetParameters( "TMassAssumption", &msets );
  vector<string> sglist;
  GetParameters( "Signals", &sglist );
  TString theFileName ;
  for (size_t i=0; i< msets.size(); i++) {
      TString massump = msets[i].substr(0,3) ;
      if ( massump == mName ) theFileName = sglist[i] + ".root" ;
  }

  TFile* file = TFile::Open( theFileName );
  TTree*  tr1 = (TTree*) file->Get( "solTt" );

  int evtId, entId, iniId, entSz, jid[5] ;
  tr1->SetBranchAddress( "evtId" ,&evtId );
  tr1->SetBranchAddress( "entId" ,&entId );
  tr1->SetBranchAddress( "entSz" ,&entSz );
  tr1->SetBranchAddress( "iniId" ,&iniId );
  tr1->SetBranchAddress( "JId"   ,&jid   );

  double brV, prob ;
  tr1->SetBranchAddress( probName ,&prob );
  tr1->SetBranchAddress( brName ,&brV );

  double mostProb = 0 ;
  double mostV = 0 ;
  for (int j=0; j< tr1->GetEntries() ; j++) {
      tr1->GetEntry(j);

      if ( prob > mostProb ) {
          mostProb = prob ;
          mostV    = brV ;
      }
      if ( (entId + 1) == (iniId + entSz) ) {
         h1->Fill( mostV, mostProb );
         mostProb = 0;
         mostV = 0 ;
      }
  }
  h1->Scale( theScale ) ;
  NormalizeComponents( "tt", h1 );
  file->Close();

}

// get signal component for fitter parameterization
void MassAnaInput::getTt( TH1D* hData, double tmass, bool matched ) {

  //int nbin = ( mH - mL )/ rbin ;

  vector<double> mlist;
  GetParameters( "TMassAssumption", &mlist );
  vector<string> flist;
  GetParameters( "Signals", &flist );

  TString fileName ;
  for (size_t i=0; i < mlist.size(); i++) {
      if ( mlist[i] == tmass ) {
         fileName = flist[i]+".root" ;
      }
  }

  TH1D* h_sg = (TH1D*) hData->Clone("h_sg") ;
  //TH1D* h_sg = new TH1D("h_sg","", nbin, mL, mH );
  get_h1Obj( fileName, "mcmTt", hname, h_sg );

  if ( !matched ) {
     get_h1Obj( fileName, "solTt", hname, hData );
     hData->Add(h_sg, -1 );
     hData->SetFillColor(kPink-2);
  }  
  if ( matched ) { 
     hData->Add( h_sg, 1 ) ;
     hData->SetFillColor(2);
  }

  NormalizeComponents( "tt" , hData );

  delete h_sg;
}

// No tt wrong permutation 
void MassAnaInput::getBackground( TH1D* allbg, THStack* bgstk, vector<TH1D*>& hlist, int groupId ){

     vector<string> channel;
     GetParameters( "channel", &channel );
     vector<int> colors;
     GetParameters("ColorCode", &colors );
     vector<string> bglist;
     GetParameters( "Backgrounds", &bglist );
     vector<int> bgGrp;
     GetParameters( "BgGroups", &bgGrp ); 

     TH1D* hini = (TH1D*) allbg->Clone("hini") ;
     TString fileName ;
     for (size_t i=0; i < bglist.size(); i++) {
         // exclude un-wanted component 
         if ( bgGrp[i] != groupId ) continue;

         TH1D* htmp = (TH1D*) hini->Clone("htmp") ;
         // chain the files if necessary
         if ( bglist[i][ bglist[i].size() -1 ] == '+'  ) {
            string ChainName = bglist[i].substr( 0, bglist[i].size()-1 ) + "Chain"  ;
            vector<string> chainlist;
            GetParameters( ChainName, &chainlist );
            TChain* fchain = new TChain("solTt") ;
            for ( size_t j=0; j< chainlist.size(); j++) {
                fileName = chainlist[j]+".root" ;
                fchain->Add( fileName );
            }
            get_h1Obj( fchain, hname, htmp );

         } else {
            fileName = bglist[i]+".root" ;
            get_h1Obj( fileName, "solTt", hname, htmp );
         }
         NormalizeComponents(  channel[i+1], htmp );

         allbg->Add( htmp, 1 );
         htmp->SetFillColor( colors[i+1] );
         hlist.push_back( htmp );
         fileName.Clear() ;
     }
     for (size_t i = hlist.size(); i > 0; i--) {
         size_t j = i-1;
         bgstk->Add( hlist[j] );
     }
}

// get the fake data 
void MassAnaInput::getFakeData( int rbin, TH1D* ttadd, THStack* ttstk, vector<TH1D*>& hlist ){

  vector<string> fNameList;
  GetParameters("FakeData", &fNameList );
  vector<string> channel;
  GetParameters("channel", &channel );
  vector<int> colors;
  GetParameters("ColorCode", &colors );

  int nbin = ( mH - mL )/ rbin ;
  for (size_t i=0; i< fNameList.size(); i++) {
      TH1D* htmp = new TH1D( channel[i].c_str() ,"", nbin, mL, mH );
      hlist.push_back( htmp );
      TString f_n( fNameList[i] );
      f_n = f_n + ".root" ;
      get_h1Obj( f_n, "solTt", hname, hlist[i], 1 );
      NormalizeComponents(  channel[i], hlist[i] );
      ttadd->Add( hlist[i], 1);
      hlist[i]->SetFillColor( colors[i] );
  }
  for (size_t i = hlist.size(); i > 0; i--) {
      size_t j = i-1;
      ttstk->Add( hlist[j] );
  }

}

vector<TLorentzVector> MassAnaInput::GetNeutrinos( TTree* tr, int evtIdx, int& synId, bool nextEvt ){

  double px,py,pz,E ;
  int evtId,objId ;
  tr->SetBranchAddress("px"    ,&px);
  tr->SetBranchAddress("py"    ,&py);
  tr->SetBranchAddress("pz"    ,&pz);
  tr->SetBranchAddress("E"     ,&E);
  tr->SetBranchAddress("evtId" ,&evtId);
  tr->SetBranchAddress("objId" ,&objId);
  vector<TLorentzVector> neuV;
  // need to be fixed , solution Id = 1 or 2; objId = 0 or 1
 
  if (synId ==0 ) neu_str = 0 ;
  if ( nextEvt ) {
     for (int i=neu_str; i<= tr->GetEntries() ; i++) {
         tr->GetEntry(i);

         //cout<<"  evtIdx = "<< evtIdx <<" evtId = "<< evtId <<" neu_str = "<< neu_str <<" objId = "<< objId ;
         //cout<<" synId = "<< synId << endl;
         if ( evtId < evtIdx ) continue;

         if ( evtId == evtIdx )  {
            //cout<<"  -> evtId = "<< evtId <<" objId ="<< objId <<" pz = "<< pz<<" px = "<< px << endl;
            TLorentzVector neuP4( px, py, pz, E );
            neuV.push_back( neuP4 );
            if ( synId != evtIdx ) { 
               synId = evtIdx ;
               neu_str = i ;
            }
         }

         if ( evtId > evtIdx ) break;
      }
      neu_str = neu_str + neuV.size()  ;
  }
  //cout<<"   final str = "<< neu_str <<endl;
  return neuV;

}

vector<TLorentzVector> MassAnaInput::GetJets( TTree* tr, int evtIdx, vector<double>& bCut, int& synId, bool nextEvt ){

  double px,py,pz,E, bDis ;
  int evtId,objId ;
  tr->SetBranchAddress("px"    ,&px);
  tr->SetBranchAddress("py"    ,&py);
  tr->SetBranchAddress("pz"    ,&pz);
  tr->SetBranchAddress("E"     ,&E);
  tr->SetBranchAddress("evtId" ,&evtId);
  tr->SetBranchAddress("objId" ,&objId);
  tr->SetBranchAddress("qCut"  ,&bDis );

  vector<TLorentzVector> jetV ;
  bCut.clear(); 

  if ( synId ==0 ) jet_str = 0 ;
  if ( nextEvt ) { 
     for (int i=jet_str; i<= tr->GetEntries() ; i++) {
         tr->GetEntry(i);

         //cout<<" evtIdx = "<< evtIdx <<" evtId = "<< evtId <<" jet_str = "<< jet_str <<" objId = "<< objId ;
         //cout<<" synId = "<< synId << endl;
         if ( evtId < evtIdx ) continue;

         if ( evtId == evtIdx )  {
            //cout<<"  -> evtId = "<< evtId <<" objId ="<< objId <<" pz = "<< pz<<" px = "<< px <<endl;
            TLorentzVector jetP4( px, py, pz, E );
            jetV.push_back( jetP4 );
            bCut.push_back( bDis );
            if ( synId != evtIdx ) { 
               synId = evtIdx ;
               jet_str = i ;
            }
         }

         if ( evtId > evtIdx ) break;
     }
     jet_str = jet_str + jetV.size()  ;
  }
  //cout<<"   final str = "<< jet_str <<endl;
  return jetV;
}

vector<TLorentzVector> MassAnaInput::GetMuons( TTree* tr, int evtIdx, int& synId, bool nextEvt ){

  double px,py,pz,E ;
  int evtId,objId ;
  tr->SetBranchAddress("px"    ,&px);
  tr->SetBranchAddress("py"    ,&py);
  tr->SetBranchAddress("pz"    ,&pz);
  tr->SetBranchAddress("E"     ,&E);
  tr->SetBranchAddress("evtId" ,&evtId);
  tr->SetBranchAddress("objId" ,&objId);

  vector<TLorentzVector> muV ;
 
  if ( synId ==0 ) mu_str = 0 ;
  if ( nextEvt ) { 
     for (int i=mu_str; i<= tr->GetEntries() ; i++) {
         tr->GetEntry(i);

         //cout<<"  evtIdx = "<< evtIdx <<" evtId = "<< evtId <<" mu_str = "<< mu_str <<" objId = "<< objId ;
         //cout<<" synId = "<< synId << endl;
         if ( evtId < evtIdx ) continue;

         if ( evtId == evtIdx )  {
            //cout<<"  -> evtId = "<< evtId <<" objId ="<< objId <<" pz = "<< pz<<" px = "<< px <<endl;
            TLorentzVector muP4( px, py, pz, E );
            muV.push_back( muP4 );
            if ( synId != evtIdx ) { 
               synId = evtIdx ;
               mu_str = i ;
            }
         }

         if ( evtId > evtIdx ) break;
     }
     mu_str = mu_str + muV.size()  ;
  }
  //cout<<"   final str = "<< mu_str <<endl;
  return muV;
}

void MassAnaInput::NormalizeComponents( double lumi, double nEvents, int theChannel, TH1D* tmp ){

  //**old for 100 /pb => ttbar signal = 8992 events ; wjets = 30393 ; QCD = 2841112
  //**new for 100 /pb => ttbar signal = 9067 events ; wjets = 34000 ; 
  //                     st_tw = 2730, st_t = 6360 
  // channel #   => ttbar signal =  1 ; wjets = 2 ; Signle Top t_ch = 3 ; Single Top tW_ch = 4 ; QCD = 5

  int idx = theChannel - 1;
  // 10TeV cross-section
  double xsec[5] = {   414,  40000,   130,    29, 121675 };  // unit:pb, for MuEnrichedQCD already applied filter efficiency
  double Eff[5]  = { 0.219, 0.0085, 0.195, 0.183, 0.2335 };  // HLT and Data Skim efficiency
  //
  // 7 TeV cross-section !!! Eff needed to be fixed !!!
  //double xsec[5] = {   187,  24000,    42,    11, 109853 };  // unit:pb, for MuEnrichedQCD already applied filter efficiency
  //double Eff[5]  = { 0.207, 0.0085, 0.195, 0.183, 0.2335 };  // HLT and Data Skim efficiency
  double nBase = xsec[idx]*Eff[idx];

  double Scal = (nBase*lumi) / nEvents ;
  tmp->Scale(Scal);

}

void MassAnaInput::NormalizeComponents(  string theChannel, TH1D* tmp ){

  double lumi ;
  GetParameters("Lumi", &lumi );

  vector<double> nEvents ;
  GetParameters( "nEvents" , &nEvents );
  vector<double> xsec;
  GetParameters("xsec", &xsec);
  vector<double> Eff;
  GetParameters("EffHLT", &Eff)  ;
  vector<string>  channel;
  GetParameters( "channel" , &channel );

  int idx = -1;
  for (size_t i =0; i< channel.size(); i++) {
      if ( channel[i] == theChannel )  idx = i ;
  }

  if ( idx >= 0 ) {
     double nBase = xsec[idx]*Eff[idx];
     double Scal = (nBase*lumi) / nEvents[idx] ;
     cout<<"-- Scal of "<< channel[idx]<< " = " << Scal <<endl;
     tmp->Scale(Scal);
  
  } else {
     cout <<" No matched componenet !! " <<endl;
  }

}

double MassAnaInput::NormalizeComponents(  string theChannel ){

  double lumi ;
  double Scal = 1 ;
  GetParameters("Lumi", &lumi );

  vector<double> nEvents ;
  GetParameters( "nEvents" , &nEvents );
  vector<double> xsec;
  GetParameters("xsec", &xsec);
  vector<double> Eff;
  GetParameters("EffHLT", &Eff)  ;
  vector<string>  channel;
  GetParameters( "channel" , &channel );

  int idx = -1;
  for (size_t i =0; i< channel.size(); i++) {
      if ( channel[i] == theChannel )  idx = i ;
  }

  if ( idx >= 0 ) {
     double nBase = xsec[idx]*Eff[idx];
     Scal = (nBase*lumi) / nEvents[idx] ;
     cout<<" Scal of "<< channel[idx]<< " = " << Scal <<endl;
  
  } else {
     cout <<" No matched componenet !! " <<endl;
  }

  return Scal;
}

void MassAnaInput::GetPermutes( int njets, vector<jlist>& jlistV ) {

    vector<int> pools ;
    for (int i=0; i < njets; i++) {
        pools.push_back(i) ;
    }

    jlistV.clear() ;

    do {
       for (int i=0; i< njets-1 ; i++) {
           for (int j=i+1; j < njets; j++ ) {
               if ( pools[0] == i && pools[1] == j ) {
                  jlist theEntry ;
                  theEntry.w1 = pools[0] ;
                  theEntry.w2 = pools[1] ;
                  theEntry.bh = pools[2] ;
                  if ( pools.size() >= 4 ) theEntry.bl = pools[3] ;
                  if ( pools.size() < 4  ) theEntry.bl = 0 ;
                  jlistV.push_back( theEntry ) ;
               }
           }
       }
    } while ( next_permutation( pools.begin() ,pools.end() ) ) ;

}


// Methods to read DataCard.txt
void MassAnaInput::GetParameters(string paraName, int* thePara ){

     fstream paraFile("DataCard.txt");
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              getValue = line.substr( vpos );
              *thePara = atoi( getValue.c_str() );
              //cout<< paraName <<" = "<< *thePara << endl;
           }
     }
     paraFile.close();
}

void MassAnaInput::GetParameters(string paraName, double* thePara ){

     fstream paraFile("DataCard.txt");
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              getValue = line.substr( vpos );
              *thePara = atof( getValue.c_str() );
              //cout<< paraName <<" = "<< *thePara << endl;
           }
     }
     paraFile.close();
}

void MassAnaInput::GetParameters(string paraName, string* thePara ){

     fstream paraFile("DataCard.txt");
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     size_t  pos ;
     size_t  vpos ;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 2;
           if ( pos < line.npos ) {
              //cout<<" pos = "<< pos <<endl;
              getName  = line.substr( pos, paraName.size() );
              //*thePara = line.substr( vpos );
              //cout<< paraName <<" = "<< *thePara << endl;
              string strTmp = line.substr( vpos );
              for (string::iterator it = strTmp.begin(); it< strTmp.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') thePara->push_back( *it );
              }
           }
     }
     paraFile.close();
}

void MassAnaInput::GetParameters(string paraName, vector<double>* thePara ){

     fstream paraFile("DataCard.txt");
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<double>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     vvec.push_back( atof( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

} 

void MassAnaInput::GetParameters(string paraName, vector<string>* thePara ){

     fstream paraFile("DataCard.txt");

     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<string>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     vvec.push_back( vtemp ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

}
 
void MassAnaInput::GetParameters(string paraName, vector<int>* thePara ){

     fstream paraFile("DataCard.txt");
     if ( !paraFile.is_open() )  cout<<" file opened error => check file path and the folder "<<endl;
     string  line;
     string  getName;
     string  getValue;
     size_t  pos ;
     size_t  vpos ;
     vector<int>  vvec;

     while ( getline(paraFile, line) ) {
           if ( line[0] == '#' ) continue ;

           pos = line.find( paraName );
           vpos = pos + paraName.size() + 1;
           if ( pos < line.npos ) {
              getName  = line.substr( pos, paraName.size() );
              string arrVal = line.substr( vpos );
	      int vidx = 0;
	      string vtemp ;
	      //cout<< paraName <<" = ( " ;
              for (string::iterator it = arrVal.begin(); it< arrVal.end(); it++) {
                  if ( (*it) != ',' && (*it) != ' ' && (*it) != '(' && (*it) != ')' && (*it) != '=') vtemp.push_back( *it );
                  if ( (*it) == ',' || (*it) == ')' ) { 
                     vvec.push_back( atoi( vtemp.c_str() ) ) ;
		     //cout<< vtemp << *it;
		     vidx++ ;
		     vtemp.clear() ;
                  }
              }
              *thePara = vvec ;
           }
     }
     paraFile.close();

}
 

