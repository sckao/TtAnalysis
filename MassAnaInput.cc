#include "MassAnaInput.h"

MassAnaInput::MassAnaInput( TString channel, int NBTag, double massL, double massH ) {

  hname = channel+"TM";
  probName = "probTt";

  mL = massL;
  mH = massH;
  ch_name = channel;
  n_btag = NBTag;
  weighting = true;
 
  bTh = 5;
  luminosity = 100;    //  pb^-1
  N_tt = 12000;
  N_wj = 120000;
  N_stt = 45018;
  N_stw = 25935;
  N_qcd = 2841112;

}

MassAnaInput::~MassAnaInput(){


}

void MassAnaInput::Initialize( TString* hfolder ) {

   if ( n_btag == -1 ) *hfolder = "NoTag_"+ch_name+"/";
   if ( n_btag == 0 )  *hfolder = "AllTags0_"+ch_name+"/";
   if ( n_btag == 1 )  *hfolder = "AllTags1_"+ch_name+"/";
   if ( n_btag == 2 )  *hfolder = "AllTags2_"+ch_name+"_10/";

}

void MassAnaInput::GetFileName( TString mName, int type, TString* fNameList ) {

  if ( n_btag == -1 ) subtag = "NB";
  if ( n_btag ==  0 ) subtag = "B0";
  if ( n_btag ==  1 ) subtag = "B1";
  if ( n_btag ==  2 ) subtag = "B2";
  TString algotag = "kcon"; 
  // File names for fake data 
  if ( type == 0 ) {
     theSG  = "pseud"+mName +"_"+algotag+"_"+subtag+".root"  ;
     theBG1 = "pseud_WJets_"+algotag+"_"+subtag+".root"  ;
     theBG2 = "pseud_STT_"  +algotag+"_"+subtag+".root"  ;
     theBG3 = "pseud_STTW_" +algotag+"_"+subtag+".root"  ;
     theBG4 = "pseud_QCD_"  +algotag+"_"+subtag+".root"  ;
  }
  // File names for templates 
  if ( type == 1 ) {
     theSG  = algotag + mName+"_"+subtag+".root"  ;
     theBG1 = algotag + "_WJets_" +subtag+".root"  ;
     theBG2 = algotag + "_STT_"   +subtag+".root"  ;
     theBG3 = algotag + "_STTW_"  +subtag+".root"  ;
     theBG4 = algotag + "_QCD_"   +subtag+".root"  ;
  }
  if ( fNameList != NULL ) {
     fNameList[0] = theSG ;
     fNameList[1] = theBG1 ;
     fNameList[2] = theBG2 ;
     fNameList[3] = theBG3 ;
     fNameList[4] = theBG4 ;
  }
}

// general method
void MassAnaInput::get_h1Obj(TString fname, TString TName, TString BName, TH1D* h1, double theScale ) {

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
  cout<<" size of tree = "<< tsz << endl;

  for (int k=0; k<tsz; k++) {
      tr1->GetEntry(k);
      h1->Fill( brV, bProb[k] );
      //h1->Fill( brV, prb );
  }

  h1->Scale( theScale ) ;
  file->Close();
 
}

void MassAnaInput::getMcMatching( TString mName, TString brName, TH1D* h1, double theScale ) {

  GetFileName( mName, 1 );
  TFile* file = TFile::Open( theSG );
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
  NormalizeComponents( luminosity, N_tt, 1., h1 );
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
  NormalizeComponents( luminosity, N_tt, 1., h1 );
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
      if ( bProb != NULL && nbtags > 0 ) { 
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
  NormalizeComponents( luminosity, N_tt, 1., h1 );
  file->Close();

}


void MassAnaInput::getMostProb( TString mName, TString brName, TH1D* h1,  double theScale ) {

  GetFileName( mName, 1 );
  TFile* file = TFile::Open( theSG );
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
  NormalizeComponents( luminosity, N_tt, 1., h1 );
  file->Close();

}


// combined the fake data
/*
void MassAnaInput::getFakeData( TString mName, int rbin, TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, TH1D* dth4, TH1D* dth5 ){

  // get the file names of fake data
  GetFileName( mName, 0 );
  // get all fake data 
  if ( dth0 != NULL && dth1 != NULL ) {
     getFakeData(dth0, theSG,  true,  rbin,  1. );                  // tt-signal
     getFakeData(dth1, theSG,  false, rbin,  1. );                  // tt-wrong combinatorics
  }
  if (dth1 != NULL && dth0 == NULL ) get_h1Obj( theSG, "solTt", hname, dth1 ); // for kinematic constrain  
  if (dth2 != NULL ) get_h1Obj( theBG1, "solTt", hname, dth2 );            // w+jets
  if (dth3 != NULL ) get_h1Obj( theBG2, "solTt", hname, dth3,  0.195 );    // single top t-channel
  if (dth4 != NULL ) get_h1Obj( theBG3, "solTt", hname, dth4,  0.183 );    // single top tW-channel
  if (dth5 != NULL ) get_h1Obj( theBG4, "solTt", hname, dth5,  2.01  );    // QCD

  // mix all fake samples
  if (dth0 != NULL ) ttadd->Add(dth0, 1);
  if (dth1 != NULL ) ttadd->Add(dth1, 1);
  if (dth2 != NULL ) ttadd->Add(dth2, 1);
  if (dth3 != NULL ) ttadd->Add(dth3, 1);
  if (dth4 != NULL ) ttadd->Add(dth4, 1);
  if (dth5 != NULL ) ttadd->Add(dth5, 1);

  // give different color for samples
  if (dth5 != NULL ) dth5->SetFillColor(5);  // QCD
  if (dth4 != NULL ) dth4->SetFillColor(6);  // single top tw channel
  if (dth3 != NULL ) dth3->SetFillColor(4);  // single top t channel
  if (dth2 != NULL ) dth2->SetFillColor(2);  // w+jets
  if (dth1 != NULL ) dth1->SetFillColor(7);
  if (dth0 != NULL ) dth0->SetFillColor(3);

  // stack them
  if (dth5 != NULL ) ttstk->Add( dth5 );
  if (dth4 != NULL ) ttstk->Add( dth4 );
  if (dth3 != NULL ) ttstk->Add( dth3 );
  if (dth2 != NULL ) ttstk->Add( dth2 );
  if (dth1 != NULL ) ttstk->Add( dth1 );
  if (dth0 != NULL ) ttstk->Add( dth0 );

}

// for tt signal and backgrounds
void MassAnaInput::getFakeData( TH1D* fkdat, TString thefileName, bool isSignal, int rbin, double theScale ){

  // all combinations
  int nbin = ( mH - mL ) / rbin ;
  TH1D* bg1 = new TH1D("bg1","", nbin, mL, mH );
  get_h1Obj( thefileName, "solTt", hname, bg1 );

  // matched combinations
  TH1D* sg1 = new TH1D("sg1","", nbin, mL, mH );
  get_h1Obj( thefileName, "mcmTt", hname, sg1 );

  if ( !isSignal ) fkdat->Add( bg1, sg1, theScale, -1.*theScale );
  if (  isSignal ) fkdat->Add( sg1, theScale );

  delete sg1;
  delete bg1;
}
*/

void MassAnaInput::getSignal(TH1D* h_Sg, int type, TString mName ) {

  GetFileName( mName, 1 );
  if ( type == 0 ) get_h1Obj(theSG, "solTt", hname, h_Sg );
  if ( type == 1 ) get_h1Obj(theSG, "mcmTt", hname, h_Sg );
  NormalizeComponents( luminosity, N_tt, 1., h_Sg );

}

void MassAnaInput::getBackground( TH1D* h_Bg, int type, int nbin, TString mName ){

  GetFileName( mName, 1 );

  TString thefileName ;
  if (type == 1 ) thefileName = theSG;
  if (type == 2 ) thefileName = theBG1 ;
  if (type == 3 ) thefileName = theBG2 ;
  if (type == 4 ) thefileName = theBG3 ;
  if (type == 5 ) thefileName = theBG4 ;

  get_h1Obj( thefileName, "solTt", hname, h_Bg );

  if( type == 1 ) {
    TH1D* sgl = new TH1D("sgl","", nbin, mL, mH );
    get_h1Obj( theSG, "mcmTt", hname, sgl );
    h_Bg->Add(sgl, -1 );
    delete sgl;
  }

  // scale to 100 /pb
  if ( type == 1 ) NormalizeComponents( luminosity, N_tt,   1, h_Bg );  // ttbar wrong permutation
  if ( type == 2 ) NormalizeComponents( luminosity, N_wj,   2, h_Bg );  // wjets
  if ( type == 3 ) NormalizeComponents( luminosity, N_stt,  3, h_Bg );  // single top t
  if ( type == 4 ) NormalizeComponents( luminosity, N_stw,  4, h_Bg );  // single top tW
  if ( type == 5 ) NormalizeComponents( luminosity, N_qcd,  5, h_Bg );  // QCD

}

void MassAnaInput::NormalizeComponents( double lumi, double nEvents, int theChannel, TH1D* tmp ){

  //**old for 100 /pb => ttbar signal = 8992 events ; wjets = 30393 ; QCD = 2841112
  //**new for 100 /pb => ttbar signal = 9067 events ; wjets = 34000 ; 
  //                     st_tw = 2730, st_t = 6360 
  // channel #   => ttbar signal =  1 ; wjets = 2 ; Signle Top t_ch = 3 ; Single Top tW_ch = 4 ; QCD = 5

  int idx = theChannel - 1;
  double xsec[5] = {   414,  40000,   130,    29, 121675 };  // unit:pb, for MuEnrichedQCD already applied filter efficiency
  double Eff[5]  = { 0.219, 0.0085, 0.195, 0.183, 0.2335 };  // HLT and Data Skim efficiency
  double nBase = xsec[idx]*Eff[idx];

  double Scal = (nBase*lumi) / nEvents ;
  tmp->Scale(Scal);

}


