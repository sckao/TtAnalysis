#include "BgEstimation.h"

BgEstimation::BgEstimation() {

  fitInput  = new MassAnaInput();
  objInfo   = new ObjectInfo();

  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "PlotType", &plotType );

  string phaseSmear ;
  fitInput->GetParameters( "PhaseSmear", &phaseSmear );
  smearing = ( phaseSmear == "ON" ) ? true : false ;

}

BgEstimation::~BgEstimation(){

  delete fitInput ;
  delete objInfo ;

}

// use 1 inclusive file 
// index : 0->all , 1->nonQCD ,  2->QCD     => using EvtSelector, MC truth selection
// index :          3->nonQCD ,  4->QCD     => using QCDSelector, control region selection
vector<double> BgEstimation::RatioXY( int njx, int njy, vector<string>& fNames, int index, bool doPlot, bool inclX, bool inclY, bool normMC ){

      ostringstream nameStr ;
      nameStr << "R" ;
      nameStr << index ;
      nameStr << "_" ;
      nameStr << "J" ;
      nameStr << njx ;
      nameStr << "J" ;
      nameStr << njy ;

      if ( fNames[0].substr(0,4) == "data" ) nameStr << "_data" ;
      cout<<" ===  Measuring Ratio "<< njx <<"/"<<njy<<" === "<< nameStr.str() <<endl;
      if ( index == 3 || index == 1 )  cout<<" **** non-QCD Region **** "<<endl;
      if ( index == 4 || index == 2 )  cout<<" ****   QCD Region   **** "<<endl;

      if ( index == 3 ) objInfo->ResetBGCuts(0, 50 );
      if ( index == 3 ) objInfo->ResetBGCuts(1, 30 );
      if ( index == 4 ) objInfo->ResetBGCuts(0, 35 );
      if ( index == 4 ) objInfo->ResetBGCuts(1, 20 );

      if ( normMC ) objInfo->Reset(1,true) ;  // Data-MC normalization
      int scaleMode =  ( normMC ) ? 2 : 0  ;  // normalize QCD and Vjets

      double theScale = 1 ;

      objInfo->Reset(0, inclX );
      objInfo->Reset(0, njx) ;
      hLepM2* nc_xj = new hLepM2( nameStr.str().substr(3,2) ) ;
      vector<double> n_tmp4 ;
      double n4Jtmp[2] = {0,0} ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          if ( index == 1 && i == 6 ) continue ; 
          if ( index == 2 && i != 6 ) continue ; 
          if ( fNames[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNames[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNames[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNames[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNames[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( fNames[i].substr(0,2) == "ww" ) theScale = fitInput->NormalizeComponents( "ww" )  ;
          if ( fNames[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" ) ;
          if ( index <= 2 && njx  > 2 && i != 0 )  objInfo->EvtSelector( fNames[i], nc_xj, smearing, theScale );
          if ( index <= 2 && njx <= 2           )  objInfo->EvtSelector( fNames[i], nc_xj, smearing, theScale );
          if ( index == 3  )  objInfo->QCDSelector( fNames[i], nc_xj, smearing, theScale, false, scaleMode );
          if ( index == 4  )  objInfo->QCDSelector( fNames[i], nc_xj, smearing, theScale, true,  scaleMode );
          n4Jtmp[0] = n4Jtmp[1] ;
          nc_xj->CounterVec( n_tmp4 );
          n4Jtmp[1] = n_tmp4[0]  ;
          cout<<njx<<" jets "<<fNames[i].substr(0,2) <<" = "<< n4Jtmp[1] - n4Jtmp[0] <<"  -> "<< n_tmp4[0]<< endl; 
          n_tmp4.clear();
      }
      vector<double> nData_xj ;
      nc_xj->CounterVec( nData_xj );

      if ( doPlot ) {
         TString chtag = nameStr.str().substr(0,5) ; 
         vector<TH1D*> hs_xj ;
	 nc_xj->Fill1DVec( hs_xj );
	 LepM2Plotter( hs_xj, chtag );
      }

      objInfo->Reset(0, inclY );
      objInfo->Reset(0, njy) ;
      hLepM2* nc_yj = new hLepM2( nameStr.str().substr(5,2) ) ;
      vector<double> n_tmp ;
      double n2Jtmp[2] = {0,0} ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          if ( index == 1 && i == 6 ) continue ; 
          if ( index == 2 && i != 6 ) continue ; 
          if ( fNames[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNames[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNames[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNames[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNames[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( fNames[i].substr(0,2) == "ww" ) theScale = fitInput->NormalizeComponents( "ww" )  ;
          if ( fNames[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" ) ;
          if ( index <= 2  )  objInfo->EvtSelector( fNames[i], nc_yj, smearing, theScale );
          if ( index == 3  )  objInfo->QCDSelector( fNames[i], nc_yj, smearing, theScale, false, scaleMode );
          if ( index == 4  )  objInfo->QCDSelector( fNames[i], nc_yj, smearing, theScale, true,  scaleMode );
          n2Jtmp[0] = n2Jtmp[1] ;
          nc_yj->CounterVec( n_tmp );
          n2Jtmp[1] = n_tmp[0]  ;
          cout<<njy<<" jets "<<fNames[i].substr(0,2) <<" = "<< n2Jtmp[1] - n2Jtmp[0] <<"  -> "<< n_tmp[0]<< endl; 
          n_tmp.clear();
      }
      vector<double> nData_yj ;
      nc_yj->CounterVec( nData_yj );

      if ( doPlot ) {
         TString chtag = nameStr.str().substr(0,3) + nameStr.str().substr(5,2);
         vector<TH1D*> hs_yj ;
	 nc_yj->Fill1DVec( hs_yj );
	 LepM2Plotter( hs_yj, chtag );
      }

      if ( doPlot ) {
         TString chtag = nameStr.str()  ;
         RatioPlotter( nData_xj, nData_yj, chtag );
      }
  
      vector<double> Rxy ;
      for ( int i=0; i< 5; i++) {
          Rxy.push_back( nData_xj[i]/nData_yj[i] ); 
      }
      for ( int i=0; i< 5; i++) {
          double RxyErr = MassFitFunction::ErrAovB( nData_xj[i], nData_yj[i], nData_xj[i+5], nData_yj[i+5] );
          Rxy.push_back( RxyErr ); 
          cout<<" Ratio"<<i<<": "<< nData_xj[i] <<" / "<<nData_yj[i]<<" = "<< nData_xj[i]/nData_yj[i] ;
          cout<<" +/- "<<RxyErr <<endl;
      }


      delete nc_xj ;
      delete nc_yj ;
      return Rxy ;
}

vector<double> BgEstimation::VQFraction( vector<string>& fNames, vector<double>& scaleV ){

      // no normalization, 
      objInfo->Reset(1, false) ;
      // set to 2 jet events only
      objInfo->Reset(0, false );
      objInfo->Reset(0, 2) ;

      hLepM2* c_vj  = new hLepM2("v2j") ;
      hLepM2* c_qcd = new hLepM2("q2j") ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          if ( fNames[i].substr(0,3) == "qcd" ) objInfo->EvtSelector( fNames[i], c_qcd, smearing, scaleV[i] );
          if ( fNames[i].substr(0,3) != "qcd" ) objInfo->EvtSelector( fNames[i], c_vj, smearing, scaleV[i] );
      }
      vector<double> n_vj ;
      c_vj->CounterVec( n_vj );
      vector<double> n_qcd ;
      c_qcd->CounterVec( n_qcd );

      double f_vj  = n_vj[0]  / (n_vj[0] + n_qcd[0] ) ;
      double f_qcd = n_qcd[0] / (n_vj[0] + n_qcd[0] ) ;

      vector<double> fracV;
      fracV.push_back( f_vj );
      fracV.push_back( f_qcd );
      cout<<"  VJ% = "<< f_vj <<"  QCD% = "<< f_qcd <<endl;

      delete c_vj ;
      delete c_qcd ;
      return fracV ;
}

vector<double> BgEstimation::VQFraction( vector<string>& fNames ){

      // normalize the MC for signal and background counter
      string mcNormalization ;
      fitInput->GetParameters( "MCNormalization", &mcNormalization );
      bool mcNorm = ( mcNormalization == "YES" ) ? true : false ;
      objInfo->Reset(1, mcNorm ) ;

      // set to 2 jet events only
      //objInfo->Reset(0, false );
      //objInfo->Reset(0, 2) ;

      double theScale = 1 ;
      hLepM2* c_vj  = new hLepM2("v2j") ;
      hLepM2* c_qcd = new hLepM2("q2j") ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          if ( fNames[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNames[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNames[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNames[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNames[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( fNames[i].substr(0,2) == "ww" ) theScale = fitInput->NormalizeComponents( "ww" )  ;
          if ( fNames[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" ) ;

          if ( fNames[i].substr(0,3) == "qcd" ) objInfo->EvtSelector( fNames[i], c_qcd, smearing, theScale );
          if ( fNames[i].substr(0,3) != "qcd" ) objInfo->EvtSelector( fNames[i], c_vj, smearing, theScale );
      }
      vector<double> n_vj ;
      c_vj->CounterVec( n_vj );
      vector<double> n_qcd ;
      c_qcd->CounterVec( n_qcd );

      double f_vj  = n_vj[0]  / (n_vj[0] + n_qcd[0] ) ;
      double f_qcd = n_qcd[0] / (n_vj[0] + n_qcd[0] ) ;

      vector<double> fracV;
      fracV.push_back( f_vj );
      fracV.push_back( f_qcd );
      cout<<"  VJ% = "<< f_vj <<"  QCD% = "<< f_qcd <<endl;

      delete c_vj ;
      delete c_qcd ;
      return fracV ;
}

vector<double> BgEstimation::AlphaSFromR21( vector<string>& fNames, vector<double>& scaleV2, vector<double>& scaleV1 ){

      // no normalization, 
      objInfo->Reset(1, false) ;
      // set to exclusive jet only
      objInfo->Reset(0, false );

      // get the alpha_s for vjets process
      objInfo->ResetBGCuts(0, 50 );
      objInfo->ResetBGCuts(1, 30 );
      hLepM2* c_v2j  = new hLepM2("v2j") ;
      hLepM2* c_v1j  = new hLepM2("v1j") ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          objInfo->Reset(0, 2) ;
          objInfo->QCDSelector( fNames[i], c_v2j, smearing, scaleV2[i], false, 0 );
          objInfo->Reset(0, 1) ;
          objInfo->QCDSelector( fNames[i], c_v1j, smearing, scaleV1[i], false, 0 );
      }
      vector<double> n_v2j ;
      c_v2j->CounterVec( n_v2j );
      vector<double> n_v1j ;
      c_v1j->CounterVec( n_v1j );
      
      double alphaV = n_v2j[0] / n_v1j[0] ;
      double berend_Vs = pow( alphaV, 2 ) + pow( alphaV, 3 ) + pow( alphaV, 4 ) ;

      // get the alpha_s for qcd process
      objInfo->ResetBGCuts(0, 35 );
      objInfo->ResetBGCuts(1, 20 );
      hLepM2* c_q2j = new hLepM2("q2j") ;
      hLepM2* c_q1j = new hLepM2("q1j") ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          objInfo->Reset(0, 2) ;
          objInfo->QCDSelector( fNames[i], c_q2j, smearing, scaleV2[i], true, 0 );
          objInfo->Reset(0, 1) ;
          objInfo->QCDSelector( fNames[i], c_q1j, smearing, scaleV1[i], true, 0 );
      }
      vector<double> n_q2j ;
      c_q2j->CounterVec( n_q2j );
      vector<double> n_q1j ;
      c_q1j->CounterVec( n_q1j );

      double alphaQ = n_q2j[0] / n_q1j[0] ;
      double berend_Qs = pow( alphaQ, 2 ) + pow( alphaQ, 3 ) + pow( alphaQ, 4 ) ;
      cout<<" a_V = "<< alphaV << " a_Q = "<< alphaQ << endl; 
      cout<<" berend_V = "<< berend_Vs << " berend_Q = "<< berend_Qs << endl; 

      vector<double>  berends ;
      berends.push_back( berend_Vs );
      berends.push_back( berend_Qs );
      berends.push_back( alphaV );
      berends.push_back( alphaQ );

      delete c_v2j ;
      delete c_v1j ;
      delete c_q2j ;
      delete c_q1j ;
      return berends ;
}

// calculate the ratio by using the alpha_s
vector<double> BgEstimation::RatioX2( vector<string>& MCNames, vector<string>& DataNames ){

   // get the Data-MC normalization
   vector<double> normV ;
   fitInput->GetParameters( "vjNorm", &normV );
   vector<double> normQ ;
   fitInput->GetParameters( "qcdNorm", &normQ );
   // get the nominal scaling for all MC samples
   vector<double> scaleV2 ;
   vector<double> scaleV1 ;
   vector<double> scaleData ;
   vector<string> channels ;
   fitInput->GetParameters( "channel", &channels, "DataCard.txt" );
   for ( size_t i=0; i< channels.size(); i++ ) {
       double theScale = fitInput->NormalizeComponents( channels[i] )  ;
       if ( i != 6 ) scaleV2.push_back( theScale*normV[1] ) ;
       if ( i == 6 ) scaleV2.push_back( theScale*normQ[1] ) ;
       if ( i != 6 ) scaleV1.push_back( theScale*normV[0] ) ;
       if ( i == 6 ) scaleV1.push_back( theScale*normQ[0] ) ;
       scaleData.push_back(1.) ;
   }

   vector<double> frac_vq = VQFraction( MCNames, scaleV2 );
   bool isRealData = ( DataNames[0].substr(0,4) == "data" ) ? true : false ;
   vector<double> berends = ( isRealData ) ? AlphaSFromR21( DataNames, scaleData, scaleData ) : AlphaSFromR21( DataNames, scaleV2, scaleV1 );

   double RX2 = frac_vq[0]*berends[0] + frac_vq[1]*berends[1] ;
   cout<<" RX2 = "<< RX2 << endl;
   vector<double> ratios ;
   ratios.push_back( RX2 );
   return ratios ;
}

vector<double> BgEstimation::RatioX2( vector<string>& MCNames, vector<string>& DataNames, vector<double>& scaleV2, vector<double>& scaleV1 ){

   // get the nominal scaling for all MC samples
   vector<double> scaleData ;
   for ( size_t i=0; i< DataNames.size(); i++ ) {
       scaleData.push_back(1.) ;
   }

   vector<double> frac_vq = VQFraction( MCNames, scaleV2 );
   bool isRealData = ( DataNames[0].substr(0,4) == "data" ) ? true : false ;
   vector<double> berends = ( isRealData ) ? AlphaSFromR21( DataNames, scaleData, scaleData ) : AlphaSFromR21( DataNames, scaleV2, scaleV1 );

   double RX2 = frac_vq[0]*berends[0] + frac_vq[1]*berends[1] ;
   vector<double> ratios ;
   ratios.push_back( RX2 );
   return ratios ;
}


// use 2 file lists which already have differnt jet multiplicity
vector<double> BgEstimation::RatioXY( vector<string>& fNameX, vector<string>& fNameY, int index, bool includeTt, bool doPlot ){

      cout<<" =====  Measuring Ratio from 2 Sets of files ===== "<<endl;
      double theScale = 1 ;
      //objInfo->Reset(0,true) ;
      string chtag = "R42" ; 

      objInfo->Reset(0, 4) ;
      hLepM2* nc_xj = new hLepM2( "Nxj", 20. ) ;
      int fNameXSize = static_cast<int>( fNameX.size() ) ;
      for ( int i = 0; i < fNameXSize; i++) {
          if ( i > 3 ) continue ; 
          if ( !includeTt && i == 0 && index <  0 ) continue; 
          if ( index > -1 && i != index ) continue ; 
          if ( fNameX[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNameX[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNameX[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNameX[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNameX[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( fNameX[i].substr(0,2) == "ww" ) theScale = fitInput->NormalizeComponents( "ww" )  ;
          if ( fNameX[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" ) ;
          objInfo->EvtSelector( fNameX[i], nc_xj, smearing, theScale );

      }
      vector<double> nData_xj ;
      nc_xj->CounterVec( nData_xj );

      if ( doPlot ) {
         vector<TH1D*> hs_xj ;
	 nc_xj->Fill1DVec( hs_xj );
	 LepM2Plotter( hs_xj, chtag );
      }

      objInfo->Reset(0, 2) ;
      hLepM2* nc_yj = new hLepM2( "Nyj", 20. ) ;
      int fNameYSize = static_cast<int>( fNameY.size() ) ;
      for ( int i = 0; i < fNameYSize ; i++) {
          if ( i > 3 ) continue ; 
          if ( !includeTt && i == 0 ) continue; 
          if ( index > -1 && i != index ) continue ; 
          if ( fNameY[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNameY[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNameY[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNameY[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNameY[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( fNameY[i].substr(0,2) == "ww" ) theScale = fitInput->NormalizeComponents( "ww" )  ;
          if ( fNameY[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" )  ;
          objInfo->EvtSelector( fNameY[i], nc_yj, smearing, theScale );
      }
      vector<double> nData_yj ;
      nc_yj->CounterVec( nData_yj );

      if ( doPlot ) {
         vector<TH1D*> hs_yj ;
	 nc_yj->Fill1DVec( hs_yj );
	 LepM2Plotter( hs_yj, chtag );
      }

      chtag = "" ;
      if ( doPlot ) RatioPlotter( nData_xj, nData_yj, chtag );

      vector<double> R42 ;
      for ( int i=0; i< 5; i++) {
          R42.push_back( nData_xj[i]/nData_yj[i] ); 
      }
      for ( int i=0; i< 5; i++) {
          double R42Err = MassFitFunction::ErrAovB( nData_xj[i], nData_yj[i], nData_xj[i+5], nData_yj[i+5] );
          R42.push_back( R42Err ); 
          cout<<" Ratio"<<i<<": "<< nData_xj[i] <<" / "<<nData_yj[i]<<" = "<< nData_xj[i]/nData_yj[i] ;
          cout<<" +/- "<<R42Err <<endl;
      }

      delete nc_xj ;
      delete nc_yj ;
      return R42 ;
}

// ratio used for systematic studies
vector<double> BgEstimation::RatioXY( vector<string>& fNames, vector<double>& scaleX, vector<double>& scaleY, int njx, int njy, bool inclX, bool inclY ){

      cout<<" ===  Measuring Ratio "<< njx <<"/"<<njy<<" === "<<endl;

      objInfo->Reset(1, false ); // no event weighting
      objInfo->Reset(0, inclX );
      objInfo->Reset(0, njx) ;
      hLepM2* nc_xj = new hLepM2( "NXJ" ) ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          if ( i != 0 ) {
             objInfo->EvtSelector( fNames[i], nc_xj, smearing, scaleX[i] );
           }
      }
      vector<double> nData_xj ;
      nc_xj->CounterVec( nData_xj );

      objInfo->Reset(0, inclY );
      objInfo->Reset(0, njy) ;
      hLepM2* nc_yj = new hLepM2( "NYJ" ) ;
      for ( size_t i = 0; i < fNames.size(); i++) {
          if ( njy > 2 && i == 0 ) continue;
          objInfo->EvtSelector( fNames[i], nc_yj, smearing, scaleY[i] );
      }
      vector<double> nData_yj ;
      nc_yj->CounterVec( nData_yj );

      
      vector<double> Rxy ;
      for ( int i=0; i< 5; i++) {
          Rxy.push_back( nData_xj[i]/nData_yj[i] ); 
      }
      for ( int i=0; i< 5; i++) {
          double RxyErr = MassFitFunction::ErrAovB( nData_xj[i], nData_yj[i], nData_xj[i+5], nData_yj[i+5] );
          Rxy.push_back( RxyErr ); 
          cout<<" Ratio"<<i<<": "<< nData_xj[i] <<" / "<<nData_yj[i]<<" = "<< nData_xj[i]/nData_yj[i] ;
          cout<<" +/- "<<RxyErr <<endl;
      }

      delete nc_xj ;
      delete nc_yj ;
      return Rxy ;
}

void BgEstimation::RatioForSystematic( string cfgFile ) {

  // open a log file
  string fpName = hfolder+"ratio_sys.log" ;
  FILE *ffile = fopen( fpName.c_str() ,"a");

  // get the nominal scaling for all MC samples
  vector<double> scaleV ;
  vector<string> channels ;
  fitInput->GetParameters( "channel", &channels, "DataCard.txt" );
  for ( size_t i=0; i< channels.size(); i++ ) {
      double theScale = fitInput->NormalizeComponents( channels[i] )  ;
      scaleV.push_back( theScale ) ;
  }

  // get the type of systematic
  string sysType ;
  fitInput->GetParameters( "SysType", &sysType, cfgFile );
  vector<double> cutV ;
  fitInput->GetParameters( sysType, &cutV, cfgFile );
  int ntest = cutV.size() ;
  cout<<" N of test = "<< ntest  <<" for "<< sysType << endl;

  // get the systematic samples
  vector<string> sysch ;
  fitInput->GetParameters( "channel", &sysch, cfgFile );

  for ( int k=0; k< ntest; k++ ) {

      // get the files for systematic combinations
      int sampleID = ( sysch.size() == 0 ) ? 0 : k+1 ;
      ostringstream fnameStr ;
      fnameStr << "sample" ;
      fnameStr << sampleID ;
      vector<string> SysFiles ;
      fitInput->GetParameters( fnameStr.str(), &SysFiles, cfgFile );
      cout<<" Testing "<< fnameStr.str() <<" -> "<<cutV[k]<<endl;

      // get the corresponding Data-MC normalization
      ostringstream normStr ;
      normStr << "Norm_" ;
      normStr << k+1 ;
      vector<double> normV ;
      fitInput->GetParameters( normStr.str(), &normV, cfgFile );

      // replace the scaling for the systematic samples
      double sysScale = ( sysch.size() == 0 ) ? 1. : fitInput->NormalizeComponents( sysch[k], cfgFile ) ;
      for ( size_t j=0; j< scaleV.size(); j++ ) { 
          cout<<"-> scale["<<j<<"]  = "<<scaleV[j] <<endl;
      }

      // current situ : sc_norm4 == sc_norm2
      vector<double> sc_norm4 ;
      vector<double> sc_norm2 ;
      for( size_t j=0; j< scaleV.size() ; j++ ) {
         double theScale = ( sysType == "Replace" && cutV[k] == j ) ?  sysScale : scaleV[j] ;
         if ( j== 1 || j==2 || j== 3 || j==4 || j ==5 ) {  
            //sc_norm4.push_back( (theScale*normV[1]*normV[1]*normV[1])/(normV[0]*normV[0]) );
            sc_norm4.push_back( theScale*normV[1] );
            sc_norm2.push_back( theScale*normV[1] );
         } else if ( j== 6 ) {
            //sc_norm4.push_back( (theScale*normV[3]*normV[3]*normV[3])/(normV[2]*normV[2]) );
            sc_norm4.push_back( theScale*normV[3] );
            sc_norm2.push_back( theScale*normV[3] );
         } else {
            sc_norm4.push_back( theScale  );
            sc_norm2.push_back( theScale  );
         }
      }
   
      for ( size_t j=0; j< scaleV.size(); j++ ) { 
          cout<<"=>  scale4["<<j<<"]  = "<<sc_norm4[j] << "=>  scale2["<<j<<"]  = "<<sc_norm2[j] <<endl;
      }

      if (sysType == "JetPt" )   objInfo->Reset( 0, cutV[k] );
      if (sysType == "MuonPt" )  objInfo->Reset( 2, cutV[k] );
      if (sysType == "MuonIso" ) objInfo->Reset( 4, cutV[k] );
      if (sysType == "SetJES"  ) objInfo->Reset( 1, static_cast<int>(cutV[k]) );
  
      //                                       inclX   inclY 
      vector<double> rxy = RatioXY( SysFiles, sc_norm4, sc_norm2, 4, 2, true,  false ) ;
      fprintf(ffile," %d -> %.4f , %.4f ,  %.4f ,  %.4f ,  %.4f ,  %.4f ,  %.4f ,  %.4f ,  %.4f ,  %.4f \n", 
                       k,  rxy[0], rxy[1], rxy[2], rxy[3], rxy[4], rxy[5], rxy[6], rxy[7], rxy[8], rxy[9] );

  }
  fclose(ffile);
}


// estimate background in different W pt regions
vector<double> BgEstimation::BgEstimate( vector<double>& R42, vector<double>& nData_2j ){

     //vector<string> fNames ;
     //fitInput->GetParameters( "FakeData", &fNames );
     //vector<double> R42 = RatioXY( 4, 2, fNames, -1, false, false ) ;

     // 2. Calculate the 4Jets background and the uncertainty 
     double bg4J = 0 ; 
     double sbg4J[2] = {0,0} ;
     double bg4JinPt[4] = { 0. } ;

     cout<<" BG Estimated : ( ";
     for(int i=1; i<5; i++) {

        bg4JinPt[i-1] = R42[i]*nData_2j[i] ;

        bg4J += R42[i]*nData_2j[i] ;
        cout<<R42[i]*nData_2j[i] <<" , " ;
        double pErr4J = MassFitFunction::ErrAxB( R42[i], nData_2j[i], R42[i+5], -1, true );
        double nErr4J = MassFitFunction::ErrAxB( R42[i], nData_2j[i], R42[i+5], -1, false );

        sbg4J[0] += pErr4J*pErr4J ;
        sbg4J[1] += nErr4J*nErr4J ;
     }
     cout<<" ) == "<< bg4J <<endl;
     vector<double> bg4Jv ;
     bg4Jv.push_back( bg4J ); 
     bg4Jv.push_back( sqrt( sbg4J[0] ) ); // upward 
     bg4Jv.push_back( sqrt( sbg4J[1] ) ); // downward
     for ( int i=0; i<4; i++){
         bg4Jv.push_back( bg4JinPt[i]  );
     }

     return bg4Jv ;
}

// estimate background using overall ratio4/2
vector<double> BgEstimation::BgEstimate( vector<double>& R42, double nData_2j ){

     //  Calculate the 4Jets background and the uncertainty 
     double bg4J = R42[0]*nData_2j ;

     double pErr4J = MassFitFunction::ErrAxB( R42[0], nData_2j, R42[5], -1, true );
     double nErr4J = MassFitFunction::ErrAxB( R42[0], nData_2j, R42[5], -1, false );

     cout<<" BG Estimated = "<< bg4J <<endl;

     vector<double> bg4Jv ;
     bg4Jv.push_back( bg4J ); 
     bg4Jv.push_back( pErr4J ); // upward 
     bg4Jv.push_back( nErr4J ); // downward
     return bg4Jv ;
}

void BgEstimation::LepM2Plotter( vector<TH1D*>& h1Ds, TString Tag  ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  TString theSubFolder = "hLepM2Info/"+Tag ;
  gSystem->mkdir( "hLepM2Info" );
  gSystem->cd( "hLepM2Info/" );
  gSystem->mkdir( Tag );
  gSystem->cd( "../../" );

  TString hNames[7] = { "W_Pt",  "W_Mt",   "SumJetPt",
                         "Mt_SumJetPt1", "Mt_SumJetPt2", "Mt_SumJetPt3", "Mt_SumJetPt4" };


  gStyle->SetOptStat("io");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetStatFontSize( 0.08);
  for (int i=0; i< 3 ; i++) {

      gStyle->SetStatFontSize( 0.08);
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd();

      h1Ds[i]->Draw();
      c1->Update();

      TString plotname = hfolder + theSubFolder +"/"+ hNames[i]+ Tag + "."+ plotType ;
      c1->Print( plotname );
      delete c1;
  }

  for (int i=0; i< 4 ; i++) {

      int j = i*2 + 3 ;
      gStyle->SetStatFontSize( 0.08);
      gStyle->SetStatFontSize( 0.08);
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->Divide(1,2);

      c1->cd(1);
      h1Ds[j]->Draw();
      c1->Update();
      c1->cd(2);
      h1Ds[j+1]->Draw();
      c1->Update();

      TString plotname = hfolder + theSubFolder +"/"+ hNames[i+3]+ "."+ plotType ;
      c1->Print( plotname );
      delete c1;
  }

}

void BgEstimation::RatioPlotter( vector<double>&  nW1, vector<double>&  nW2, TString Tag, double scale  ){

  TCanvas* c1 = new TCanvas("c1","", 750, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->cd();

  TF1 *func = new TF1("func", MassFitFunction::fitPoly, 0, 200, 2);
  func->FixParameter(0, nW1[0]/nW2[0] ) ;
  func->FixParameter(1,    0  ) ;
  func->SetLineColor(2);
  func->SetLineWidth(2);

  double ratio[4] = { nW1[1]/nW2[1], nW1[2]/nW2[2], nW1[3]/nW2[3], nW1[4]/nW2[4]  };
  double ptcut[4] = {            10,            30,            50,            70  };
  double rErr[4] = {0.};
  double ptRange[4] = {0.};
  for (int i=0; i<4; i++) {
      rErr[i] = MassFitFunction::ErrAovB( nW1[i+1], nW2[i+1], nW1[i+6], nW2[i+6] );
      ptRange[i] = 10. ;
      //if (i==4) ptRange[i] = 60. ;
  }

  TGraphErrors* hRatio = new TGraphErrors( 4, ptcut, ratio, ptRange, rErr );
  hRatio->SetMaximum( 0.6 );
  hRatio->SetMinimum( 0.05 );
  hRatio->GetXaxis()->SetLimits(0,80);
  hRatio->SetMarkerColor(4);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(2.0);
  hRatio->SetLineWidth(2);
  hRatio->SetTitle("Ratio vs Pt of Mt");

  hRatio->Draw("AP");
  func->Draw("SAME");

  c1->Update();

  TString plotname = hfolder + "hLepM2Info/"+ Tag + "."+plotType ;
  c1->Print( plotname );

  delete hRatio ;
  delete func ;
  delete c1 ;
}


vector<double> BgEstimation::MtFitter( string& DataName, vector<string>& fNames, int phaseIdx ){

     int k = 1 + phaseIdx*2 ;

     // Get Data
     hLepM2* hs_data = new hLepM2( "hs_data" ) ;
     objInfo->EvtSelector( DataName, hs_data, false, 1 );
     vector<TH1D*> hData_v ;
     hs_data->Fill1DVec( hData_v );
     double N_Data = hData_v[k]->Integral() ;
   
     int nbin = 30 ;
     if ( hData_v[k]->Integral() < 500. && phaseIdx < 3) {
        nbin = 15 ;
        hData_v[k]->Rebin(2);
     }

     vector<TH1D*> hTemplates ;


     // template 1 : W + Z + Tt
     double scaleWj = fitInput->NormalizeComponents( "wj" )  ;
     double scaleZj = fitInput->NormalizeComponents( "zj" )  ;
     double scaleTt = fitInput->NormalizeComponents( "tt" )  ;
     double scaleQCD = fitInput->NormalizeComponents( "qcd" )  ;

     //objInfo->Reset(0, true );
     //objInfo->Reset(0, 0 ) ;
     hLepM2* h_wj = new hLepM2( "h_wj", nbin ) ;
     objInfo->EvtSelector( fNames[1], h_wj, false, scaleWj );
     vector<TH1D*> hW_xj ;
     h_wj->Fill1DVec( hW_xj );
   
     hLepM2* h_zj = new hLepM2( "h_zj", nbin ) ;
     objInfo->EvtSelector( fNames[3], h_zj, false, scaleZj );
     vector<TH1D*> hZ_xj ;
     h_zj->Fill1DVec( hZ_xj );

     hLepM2* h_tt = new hLepM2( "h_tt", nbin ) ;
     objInfo->EvtSelector( fNames[0], h_tt, false, scaleTt );
     vector<TH1D*> hT_xj ;
     h_tt->Fill1DVec( hT_xj );

     hW_xj[k]->Add( hT_xj[k], 1 );
     hW_xj[k]->Add( hZ_xj[k], 1 );
     double N_VJ = hW_xj[k]->Integral() ;

     hW_xj[k]->Scale(  150000./hW_xj[k]->Integral()  );

    // template 2 : QCD
     hLepM2* h_qcd = new hLepM2( "h_qcd", nbin ) ;
     objInfo->EvtSelector( fNames[2], h_qcd, false, scaleQCD );
     vector<TH1D*> hQ_xj ;
     h_qcd->Fill1DVec( hQ_xj );

     double N_QCD = hQ_xj[k]->Integral() ;

     hQ_xj[k]->Scale(  150000./hQ_xj[k]->Integral()  );

     hTemplates.push_back( hQ_xj[k] );
     hTemplates.push_back( hW_xj[k] );

     // Fit 
     int fbin =  ( nbin != 30 || phaseIdx > 2 ) ? 11 : 22 ;
     vector<double> fRatio = JacobianFitter( hData_v[k], hTemplates, phaseIdx, 1, fbin );
     vector<double> fracN ;
     for ( size_t i=0; i< fRatio.size(); i++){
         fracN.push_back( N_Data*fRatio[i] ) ;
     }
     cout <<" QCD     scale = "<< N_Data*fRatio[0]/N_QCD <<endl;
     cout <<" Non-QCD scale = "<< N_Data*fRatio[1]/N_VJ  <<endl;

     delete h_wj ;
     delete h_zj ;
     delete h_tt ;
     delete h_qcd ;
     delete hs_data ;
     return fracN ;
}

vector<double> BgEstimation::JacobianFitter( TH1D* hData, vector<TH1D*>& hTemplates, int phaseIdx, int fbin1, int fbin2 ){

  TString theFolder = hfolder ;
  TString theSubFolder = "JbFitter/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );   

   gStyle->SetOptFit(111);
   gStyle->SetOptStat("ieo");
   TCanvas* c1 = new TCanvas("c1","", 900, 700);
   c1->SetGrid();
   c1->SetFillColor(10);
   c1->SetFillColor(10);
   c1->cd();

   int nTemp = hTemplates.size() ;
   TObjArray* temps = new TObjArray( nTemp );
   for ( int i=0; i< nTemp; i++ ) {
      temps->Add( hTemplates[i] ) ;
   }

   TFractionFitter* jbfit = new TFractionFitter( hData, temps );  
   jbfit->SetRangeX( fbin1, fbin2 ) ;
   Int_t status = jbfit->Fit();
   cout<<" fit status = "<< status <<endl;

   vector<double> fRatio ;
   if ( status == 0 || status == 4 ) {         
      double ratio[10];
      double err[10];
      for ( int i=0; i< nTemp; i++ ) {
          jbfit->GetResult( i, ratio[i], err[i] );
          fRatio.push_back( ratio[i] );
      }
     
      double theScale = hData->Integral() / hTemplates[0]->Integral() ;

      THStack* hStk = new THStack("hStk", "LepM2T_Fit" );
      hTemplates[0]->SetFillColor(kYellow);
      hTemplates[1]->SetFillColor(kGreen);

      hTemplates[0]->Scale( ratio[0]*theScale );
      hTemplates[1]->Scale( ratio[1]*theScale );
      hStk->Add( hTemplates[1] );
      hStk->Add( hTemplates[0] );
      hData->SetMarkerSize(1);
      hData->SetMarkerStyle(21);
      hData->Draw("PE");
      c1->Update();
      hStk->Draw("same");
      c1->Update();

      hData->Draw("PE SAME");
      c1->Update();
   }
   
   char xchar[2];
   sprintf( xchar, "%d", phaseIdx ) ;
   string xtag = xchar ;

   TString plotname = hfolder + theSubFolder + "LepWFit_" + xtag + "." +plotType ;
   c1->Print( plotname );
   

   delete jbfit ;
   //delete temps ;
   delete c1 ;
   return fRatio ;
}

vector<double> BgEstimation::MeasureScale( string& dataFile, vector<string>& fNames  ){

   cout<<" ---  MC Scaling Measuring --- "<<endl;
       
   double scale1 = fitInput->NormalizeComponents( "wj" );
   double scale2 = fitInput->NormalizeComponents( "qcd" );
   double scale3 = fitInput->NormalizeComponents( "zj" );

   bgCounter* hwj = new bgCounter("wj_") ;
   objInfo->EvtSelector( fNames[1], hwj, false, scale1 );
   vector<double> h_wj ;
   hwj->CounterVec( h_wj );

   bgCounter* hzj = new bgCounter("zj_") ;
   objInfo->EvtSelector( fNames[3], hzj, false, scale3 );
   vector<double> h_zj ;
   hzj->CounterVec( h_zj );

   bgCounter* hqcd = new bgCounter("qcd_") ;
   objInfo->EvtSelector( fNames[2], hqcd, false, scale2 );
   vector<double> h_qcd ;
   hqcd->CounterVec( h_qcd );

   vector<double> mcScale ;
   for ( int i=0; i< 5; i++) {
       double nqcd0 = h_qcd[i] ;
       double nV0   = h_wj[i] + h_zj[i] ;
       vector<double> fitR = MtFitter( dataFile, fNames, i );
       cout<<" Phase "<<i<<" QCDNorm = "<< fitR[0]/ nqcd0 <<"  VNorm = "<< fitR[1]/nV0 <<endl ;
       mcScale.push_back( fitR[0]/ nqcd0 );
       mcScale.push_back( fitR[1]/ nV0 );
   }

   delete hwj;
   delete hzj;
   delete hqcd;
   return mcScale ;
}

double BgEstimation::MeasureScale2D( string& DataName, vector<string>& fakeData, double MtCut, double METCut, bool doQCD, bool isVjNorm ){

  objInfo->ResetBGCuts(0, MtCut);
  objInfo->ResetBGCuts(1, METCut);
  objInfo->Reset(1, isVjNorm );  

  int nbin_ = 2 ;

  hObjs* hdt = new hObjs("data", nbin_ ) ;
  objInfo->QCDSelector( DataName, hdt, false, 1, doQCD );
  vector<TH2D*> h_dt ;
  hdt->Fill2DVec( h_dt );

  if ( h_dt[0]->Integral() < 500. ) {
     nbin_ = 15 ;
     h_dt[0]->RebinY(2);
  }

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  hObjs* htt = new hObjs("tt", nbin_ ) ;
  objInfo->QCDSelector( fakeData[0], htt, false, scale0, doQCD );
  vector<TH2D*> h_tt ;
  htt->Fill2DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[1], hwj, false, scale1, doQCD );
  vector<TH2D*> h_wj ;
  hwj->Fill2DVec( h_wj );

  hObjs* hzj = new hObjs("zj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[2], hzj, false, scale2, doQCD );
  vector<TH2D*> h_zj ;
  hzj->Fill2DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin_ ) ;
  objInfo->QCDSelector( fakeData[3], htq, false, scale3, doQCD );
  vector<TH2D*> h_tq ;
  htq->Fill2DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin_ ) ;
  objInfo->QCDSelector( fakeData[4], htw, false, scale4, doQCD );
  vector<TH2D*> h_tw ;
  htw->Fill2DVec( h_tw );

  hObjs* hww = new hObjs("ww", nbin_ ) ;
  objInfo->QCDSelector( fakeData[5], hww, false, scale5, doQCD );
  vector<TH2D*> h_ww ;
  hww->Fill2DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd", nbin_ ) ;
  objInfo->QCDSelector( fakeData[6], hqcd, false, scale6, doQCD );
  vector<TH2D*> h_qcd ;
  hqcd->Fill2DVec( h_qcd );

  TH2D* hMC   = new TH2D("hMC", " MET(X) vs lep Mt(Y) ",   30, 0, 150,  nbin_, 0, 150 );
  TH2D* hMC_A = new TH2D("hMC_A", " MET(X) vs lep Mt(Y) ", 30, 0, 150,  nbin_, 0, 150 );
  TH2D* hData = new TH2D("hData", " MET(X) vs lep Mt(Y) ", 30, 0, 150,  nbin_, 0, 150 );

  hMC_A->Add( h_qcd[0] );
  hMC_A->Add( h_wj[0] );
  hMC_A->Add( h_zj[0] );
  hMC_A->Add( h_tt[0] );
  hMC_A->Add( h_tq[0] );
  hMC_A->Add( h_tw[0] );
  hMC_A->Add( h_ww[0] );

  if ( doQCD ) { 
     hMC->Add( h_qcd[0] );
     hData->Add( h_dt[0],  1. );
     hData->Add( h_tt[0], -1. );
     hData->Add( h_wj[0], -1. );
     hData->Add( h_zj[0], -1. );
     hData->Add( h_tq[0], -1. );
     hData->Add( h_tw[0], -1. );
     hData->Add( h_ww[0], -1. );
  }
  if ( !doQCD ) { 
     hMC->Add( h_wj[0] );
     hMC->Add( h_zj[0] );
     hMC->Add( h_tt[0] );
     hMC->Add( h_tq[0] );
     hMC->Add( h_tw[0] );
     hMC->Add( h_ww[0] );
     hData->Add( h_dt[0],  1. );
     hData->Add( h_qcd[0], -1. );
  }
  cout<<" N of MC QCD = "<< h_qcd[0]->Integral() <<endl;
  cout<<" N of MC WJ  = "<< h_wj[0]->Integral() <<endl;
  cout<<" N of MC zJ  = "<< h_zj[0]->Integral() <<endl;
  cout<<" N of MC tt  = "<< h_tt[0]->Integral() <<endl;
  cout<<" N of MC tq  = "<< h_tq[0]->Integral() <<endl;
  cout<<" N of MC tw  = "<< h_tw[0]->Integral() <<endl;
  cout<<" N of MC tw  = "<< h_ww[0]->Integral() <<endl;
  cout<<" N of Data = "<< h_dt[0]->Integral() <<" -> "<< hData->Integral() <<endl ;

  cout<<" Event # normalization = "<< hData->Integral()/hMC->Integral() <<endl ;
 
  int wbin = 150/nbin_ ;
  int bx1 = ( doQCD ) ?        1 :  METCut/5 ;
  int bx2 = ( doQCD ) ? METCut/5 :        20 ;
  int by1 = ( doQCD ) ?        1 :  MtCut/wbin ;
  int by2 = ( doQCD ) ?  MtCut/5 :    100/wbin ;
  double bestNorm =  Chi2Normalization( hData, hMC, bx1, bx2, by1, by2, 1, scale2 );
  cout<<" The best QCD Normalization  = "<< bestNorm <<endl;

  TString theFolder = hfolder ;
  TString theSubFolder = "hQCDMu18_Ex/" ;
  if ( !doQCD ) theSubFolder = "hSGMu18/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  //gStyle->SetOptStat("nieuo");
  //gStyle->SetPalette(1);
  //gStyle->SetStatX(0.85);

  gStyle->SetOptStat("");
  //gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPalette(1);
  gStyle->SetTitleFontSize(0.08);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  double xR = 120. ;
  double yR = 120. ;

  c1->cd(1);
  gStyle->SetNumberContours(10);
  //h_dt[0]->SetName("Data");
  h_dt[0]->SetTitle("Data");
  h_dt[0]->SetLabelSize( 0.08, "X");  
  h_dt[0]->SetLabelSize( 0.08, "Y");  
  h_dt[0]->SetLabelSize( 0.08, "Z");  
  h_dt[0]->SetAxisRange(0,xR,"X") ;
  h_dt[0]->SetAxisRange(0,yR,"Y") ;
  h_dt[0]->Draw("COLZ");
  c1->Update();

  c1->cd(2);
  gStyle->SetNumberContours(4);
  //h_qcd[0]->SetName("QCD");
  h_qcd[0]->SetTitle("QCD MC");
  h_qcd[0]->SetLabelSize( 0.08, "X");
  h_qcd[0]->SetLabelSize( 0.08, "Y");  
  h_qcd[0]->SetLabelSize( 0.08, "Z");  
  h_qcd[0]->SetAxisRange(0,xR,"X") ;
  h_qcd[0]->SetAxisRange(0,yR,"Y") ;
  h_qcd[0]->Draw("COLZ");
  c1->Update();

  c1->cd(3);
  gStyle->SetNumberContours(10);
  //h_wj[0]->SetName("WJets");
  h_wj[0]->SetTitle("WJets MC");
  h_wj[0]->SetLabelSize( 0.08, "X");
  h_wj[0]->SetLabelSize( 0.08, "Y");  
  h_wj[0]->SetLabelSize( 0.08, "Z");  
  h_wj[0]->SetAxisRange(0,xR,"X") ;
  h_wj[0]->SetAxisRange(0,yR,"Y") ;
  h_wj[0]->Draw("COLZ");
  c1->Update();

  c1->cd(4);
  gStyle->SetNumberContours(10);
  //hMC_A->SetName("MC");
  hMC_A->SetTitle("Mix MC");
  hMC_A->SetLabelSize( 0.08, "X");
  hMC_A->SetLabelSize( 0.08, "Y");  
  hMC_A->SetLabelSize( 0.08, "Z");  
  hMC_A->SetAxisRange(0,xR,"X") ;
  hMC_A->SetAxisRange(0,yR,"Y") ;
  hMC_A->Draw("COLZ");
  c1->Update();

  gStyle->SetStatX(0.95);

  ostringstream pNameStr ;
  pNameStr << "_" ;
  pNameStr << MtCut ;
  pNameStr << "-" ;
  pNameStr << METCut ;
  TString pNameTStr = pNameStr.str() ;
  TString plotname1 = theFolder + theSubFolder +  "ScaleCheck" + pNameTStr + "." + plotType ;
  c1->Print( plotname1 );

  delete c1;
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;
  delete hMC_A ;  
  delete hMC ;  
  delete hData ;
  return bestNorm ;

}

double BgEstimation::Chi2Normalization( TH2D* hData, TH2D* hMC, int Bx1, int Bx2, int By1, int By2, double s1, double s2, bool doFit ){

   // 
   double Itgr1 = hData->Integral() ;
   double Itgr2 = hMC->Integral() ;

   double d_norm = Itgr1/Itgr2 ;
   double norm = d_norm - (0.34*d_norm) ;
   double step = (d_norm*0.34)/25. ;
   double bestNorm = d_norm ;
 
   cout<<" delta Norm = "<< d_norm*0.34 <<"  Norm start from = "<< norm <<endl ;

   double minX2 = 99999. ;

   double nm[50];
   double chi2[50];
   for( int k=0; k<50; k++ ) {
      norm = norm + step ;
      // chi2 sum
      double X2 = 0 ;
      for (int i= Bx1; i<= Bx2; i++){
          for (int j= By1; j<= By2; j++){
              double N1 = hData->GetBinContent(i,j);
              double N2 = hMC->GetBinContent(i,j);
              if ( N1 < 0 ) N1 = 0 ;
              if ( N2 <= 0. ) continue;
              double x2_ij = (N1 - (N2*norm) )*(N1 - (N2*norm) ) / (s1*N1 + s2*N2*norm*norm ) ;
              X2 += x2_ij ;
          }
      }
      nm[k] = norm ;
      chi2[k] = X2 ;
      //cout<<" Norm =  "<<norm <<"   X2 = "<< X2 <<"   step : "<< step <<endl; ;

      if ( k==0 ) minX2 = X2 ;
      if ( X2 < minX2 && k > 0 ) { 
         minX2 = X2 ;
         bestNorm = norm ;
      }
   }

   if ( doFit ) {
      gStyle->SetOptFit(111);
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd();

      TString theFolder = hfolder ;
      TString theSubFolder = "hQCDBG/" ;
      gSystem->mkdir( theFolder );
      gSystem->cd( theFolder );
      gSystem->mkdir( theSubFolder );
      gSystem->cd( "../" );

      TGraph* hNormX2 = new TGraph(50, nm, chi2 );

      double fitMin = ( nm[0] < nm[49] ) ? nm[0] : nm[49] ;
      double fitMax = ( fitMin == nm[0] ) ? nm[49] : nm[0] ;
      TF1 *func1 = new TF1("func1",MassFitFunction::fitParabola, fitMin, fitMax, 3 );
      func1->SetParLimits(0, fitMin, fitMax );

      hNormX2->SetMarkerSize(1);
      hNormX2->SetMarkerColor(4);
      hNormX2->SetMarkerStyle(21);
      hNormX2->Draw("AP");
      c1->Update();

      hNormX2->Fit( func1, "R", "sames", fitMin, fitMax );
      double fitNorm = func1->GetParameter(0);
      double fChi2 = func1->GetChisquare();
      int fNDF  = func1->GetNDF() ;
      double nChi2 = fChi2 / static_cast<double>( fNDF ) ;
      cout<<" The Fitted Normalization = "<< fitNorm <<" the best norm = "<< bestNorm <<endl;
      if ( nChi2 < 10. ) bestNorm = fitNorm ;
      c1->Update();

      TString plotname1 = theFolder + theSubFolder +  "X2Fit2D."+plotType ;
      c1->Print( plotname1 );

      delete hNormX2;
      delete func1 ;
      delete c1;
   }
   cout<<"  Final Normalization = "<< bestNorm << endl;

   return bestNorm ;
}

double BgEstimation::MeasureScale1D( string& DataName, vector<string>& fakeData, int hID, double MtCut, double METCut, bool doQCD, bool isVjNorm, int njets ){

  objInfo->Reset(0, njets );
  objInfo->ResetBGCuts(0, MtCut);
  objInfo->ResetBGCuts(1, METCut);
  objInfo->Reset(0, false );  // set to Exclusive
  objInfo->Reset(1, isVjNorm );  
  int scaleMode =  ( isVjNorm ) ? 1 : 0 ;

  int nbin_ = 2 ;

  hObjs* hdt = new hObjs("data", nbin_ ) ;
  objInfo->QCDSelector( DataName, hdt, false, 1., doQCD );
  vector<TH1D*> h_dt ;
  hdt->Fill1DVec( h_dt );

  if ( h_dt[13]->Integral() < 500. ) {
     nbin_ = 15 ;
     h_dt[13]->Rebin(2);
     h_dt[14]->Rebin(2);
     h_dt[15]->Rebin(2);
     cout<<" bin width changes -> "<< h_dt[13]->GetBinWidth(1) << endl;
  }

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  hObjs* htt = new hObjs("tt", nbin_ ) ;
  objInfo->QCDSelector( fakeData[0], htt, false, scale0, doQCD, scaleMode );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[1], hwj, false, scale1, doQCD, scaleMode );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hzj = new hObjs("zj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[2], hzj, false, scale2, doQCD, scaleMode );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin_ ) ;
  objInfo->QCDSelector( fakeData[3], htq, false, scale3, doQCD, scaleMode );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin_ ) ;
  objInfo->QCDSelector( fakeData[4], htw, false, scale4, doQCD, scaleMode );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  hObjs* hww = new hObjs("ww", nbin_ ) ;
  objInfo->QCDSelector( fakeData[5], hww, false, scale5, doQCD, scaleMode );
  vector<TH1D*> h_ww ;
  hww->Fill1DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd", nbin_ ) ;
  objInfo->QCDSelector( fakeData[6], hqcd, false, scale6, doQCD, scaleMode );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  int nBin = h_dt[hID]->GetNbinsX() ;
  cout<<" N bin set = "<< nBin <<endl;
  TH1D* hMC   = new TH1D("hMC",   "MC   lep Mt ", nBin, 0, 150 );
  TH1D* hData = new TH1D("hData", "Data lep Mt ", nBin, 0, 150 );
  if ( doQCD ) { 
     hMC->Add( h_qcd[hID] );
     hData->Add( h_dt[hID],  1. );
     hData->Add( h_tt[hID], -1. );
     hData->Add( h_wj[hID], -1. );
     hData->Add( h_zj[hID], -1. );
     hData->Add( h_tq[hID], -1. );
     hData->Add( h_tw[hID], -1. );
     hData->Add( h_ww[hID], -1. );
  }
  if ( !doQCD ) { 
     hMC->Add( h_wj[hID] );
     hMC->Add( h_zj[hID] );
     hMC->Add( h_tt[hID] );
     hMC->Add( h_tq[hID] );
     hMC->Add( h_tw[hID] );
     hMC->Add( h_ww[hID] );
     hData->Add( h_dt[hID],  1. );
     hData->Add( h_qcd[hID], -1. );
  }

  cout<<" N of MC QCD = "<< h_qcd[hID]->Integral() <<endl;
  cout<<" N of MC WJ  = "<< h_wj[hID]->Integral() <<endl;
  cout<<" N of MC zJ  = "<< h_zj[hID]->Integral() <<endl;
  cout<<" N of MC tt  = "<< h_tt[hID]->Integral() <<endl;
  cout<<" N of MC tq  = "<< h_tq[hID]->Integral() <<endl;
  cout<<" N of MC tw  = "<< h_tw[hID]->Integral() <<endl;
  cout<<" N of MC ww  = "<< h_ww[hID]->Integral() <<endl;
  cout<<" N of Data = "<< h_dt[hID]->Integral() <<" -> "<< hData->Integral() <<endl ;

  cout<<" Event # normalization = "<< hData->Integral()/hMC->Integral() <<endl ;
  int wbin = 150/nBin ;
  int bin1 = ( doQCD ) ?         1    : (MtCut/wbin) + 1 ;
  int bin2 = ( doQCD ) ? (MtCut/wbin) :         100/wbin ;
  cout<<" bin1 = "<< bin1 <<"   bin2 = "<<bin2<<endl;

  ostringstream pNameStr ;
  pNameStr << MtCut ;
  pNameStr << hID ;
  pNameStr << "_" ;
  pNameStr << njets ;
  TString pNameTStr = pNameStr.str() ;
  double bestNorm =  Chi2Normalization( hData, hMC, bin1, bin2, 1, 1, pNameTStr, true );

  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hww;
  delete hzj;
  delete hqcd;
  delete hMC ;  
  delete hData ;

  return bestNorm ;
}

// MeasureScale for systematic samples
void BgEstimation::MeasureScale1D( string& DataName, string cfgFile, int hID, double MtCut, double METCut, bool isQCD, bool isVjNorm, int njets ) {

  // open a log file
  string fpName = hfolder+"norm_VJ.log" ;
  FILE *ffile = fopen( fpName.c_str() ,"a");

  // get the nominal scaling for all MC samples, *** Be Careful the Order ****
  vector<double> scaleV ;
  vector<string> channels ;
  fitInput->GetParameters( "channel", &channels, "DataCard.txt" );
  for ( size_t i=0; i< channels.size(); i++ ) {
      double theScale = fitInput->NormalizeComponents( channels[i] )  ;
      scaleV.push_back( theScale ) ;
  }

  // get the type of systematic
  string sysType ;
  fitInput->GetParameters( "SysType", &sysType, cfgFile );
  vector<double> cutV ;
  fitInput->GetParameters( sysType, &cutV, cfgFile );
  int nTest = cutV.size() ;

  vector<string> sysch ;
  fitInput->GetParameters( "channel", &sysch, cfgFile );
   
  for ( int k=0; k< nTest; k++ ) {

      // get the files for systematic combinations -> only for Vjets and ttbar
      int sampleID = ( sysch.size() == 0 ) ? 0 : k+1 ;
      ostringstream fnameStr ;
      fnameStr << "sample" ;
      fnameStr << sampleID ;
      vector<string> SysFiles ;
      fitInput->GetParameters( fnameStr.str(), &SysFiles, cfgFile );
      cout<<" Testing "<< sysType <<" -> "<< fnameStr.str()<<"  => "<< cutV[k] <<endl;

      // replace the scaling for the systematic samples -> only for Vjets and ttbar
      double sysScale = ( sysch.size() == 0 ) ? 1. : fitInput->NormalizeComponents( sysch[k], cfgFile ) ;

      // get the corresponding Data-MC normalization
      ostringstream normStr ;
      normStr << "Norm_" ;
      normStr << k+1 ;
      vector<double> normV ;
      fitInput->GetParameters( normStr.str(), &normV, cfgFile );
      double normv1 = ( isVjNorm ) ? normV[ njets-1 ] : 1. ;

      // Reset all the scaling*normalization values
      vector<double> sc_norm ;
      for( size_t j=0; j< scaleV.size() ; j++ ) {
         double theScale = ( sysType == "Replace" && cutV[k] == j ) ?  sysScale : scaleV[j] ;
         if ( j== 1 || j==2 || j==3 || j== 4 || j ==5 ) {  
            sc_norm.push_back( theScale*normv1 );
         } else {
            sc_norm.push_back( theScale  );
         }
      }

      if (sysType == "MuonPt"  ) objInfo->Reset( 2, cutV[k] );
      if (sysType == "MuonIso" ) objInfo->Reset( 4, cutV[k] );
      if (sysType == "JetPt"  ) objInfo->Reset( 0, cutV[k] );
      if (sysType == "SetJES"  ) objInfo->Reset( 1, static_cast<int>(cutV[k]) );

      objInfo->Reset(0, njets );
      objInfo->ResetBGCuts(0, MtCut);
      objInfo->ResetBGCuts(1, METCut);
      objInfo->Reset(0, false );  // set to Exclusive
      objInfo->Reset(1, false );  // no event weighting , it's done by scaling*normalization
     
      int nbin_ = 2 ;

      hObjs* hdt = new hObjs("data", nbin_ ) ;
      objInfo->QCDSelector( DataName, hdt, false, 1., isQCD );
      vector<TH1D*> h_dt ;
      hdt->Fill1DVec( h_dt );

      if ( h_dt[13]->Integral() < 500. ) {
         nbin_ = 15 ;
         h_dt[13]->Rebin(2);
         h_dt[14]->Rebin(2);
         h_dt[15]->Rebin(2);
         cout<<" bin width changes -> "<< h_dt[13]->GetBinWidth(1) << endl;
      }

      hObjs* htt = new hObjs("tt", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[0], htt, false, sc_norm[0], isQCD );
      vector<TH1D*> h_tt ;
      htt->Fill1DVec( h_tt );

      hObjs* hwj = new hObjs("wj", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[1], hwj, false, sc_norm[1], isQCD );
      vector<TH1D*> h_wj ;
      hwj->Fill1DVec( h_wj );

      hObjs* hzj = new hObjs("zj", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[2], hzj, false, sc_norm[2], isQCD );
      vector<TH1D*> h_zj ;
      hzj->Fill1DVec( h_zj );

      hObjs* htq = new hObjs("tq", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[3], htq, false, sc_norm[3], isQCD );
      vector<TH1D*> h_tq ;
      htq->Fill1DVec( h_tq );

      hObjs* htw = new hObjs("tw", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[4], htw, false, sc_norm[4], isQCD );
      vector<TH1D*> h_tw ;
      htw->Fill1DVec( h_tw );

      hObjs* hww = new hObjs("ww", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[5], hww, false, sc_norm[5], isQCD );
      vector<TH1D*> h_ww ;
      hww->Fill1DVec( h_ww );

      hObjs* hqcd = new hObjs("qcd", nbin_ ) ;
      objInfo->QCDSelector( SysFiles[6], hqcd, false, sc_norm[6], isQCD );
      vector<TH1D*> h_qcd ;
      hqcd->Fill1DVec( h_qcd );

      int nBin = h_dt[hID]->GetNbinsX() ;
      cout<<" N bin set = "<< nBin <<endl;
      TH1D* hMC   = new TH1D("hMC",   "MC   lep Mt ", nBin, 0, 150 );
      TH1D* hData = new TH1D("hData", "Data lep Mt ", nBin, 0, 150 );
      if ( isQCD ) { 
         hMC->Add( h_qcd[hID] );
	 hData->Add( h_dt[hID],  1. );
	 hData->Add( h_tt[hID], -1. );
	 hData->Add( h_wj[hID], -1. );
	 hData->Add( h_zj[hID], -1. );
	 hData->Add( h_tq[hID], -1. );
	 hData->Add( h_tw[hID], -1. );
	 hData->Add( h_ww[hID], -1. );
      }
      if ( !isQCD ) { 
         hMC->Add( h_wj[hID] );
	 hMC->Add( h_zj[hID] );
	 hMC->Add( h_tt[hID] );
	 hMC->Add( h_tq[hID] );
	 hMC->Add( h_tw[hID] );
	 hMC->Add( h_ww[hID] );
	 hData->Add( h_dt[hID],  1. );
	 hData->Add( h_qcd[hID], -1. );
      }

      cout<<" N of MC QCD = "<< h_qcd[hID]->Integral() <<endl;
      cout<<" N of MC WJ  = "<< h_wj[hID]->Integral() <<endl;
      cout<<" N of MC zJ  = "<< h_zj[hID]->Integral() <<endl;
      cout<<" N of MC tt  = "<< h_tt[hID]->Integral() <<endl;
      cout<<" N of MC tq  = "<< h_tq[hID]->Integral() <<endl;
      cout<<" N of MC tw  = "<< h_tw[hID]->Integral() <<endl;
      cout<<" N of MC ww  = "<< h_ww[hID]->Integral() <<endl;
      cout<<" N of Data = "<< h_dt[hID]->Integral() <<" -> "<< hData->Integral() <<endl ;

      double evtNorm = hData->Integral()/hMC->Integral() ;
      cout<<" Event # normalization = "<< hData->Integral()/hMC->Integral() <<endl ;
      int wbin = 150/nBin ;
      int bin1 = ( isQCD ) ?         1    : (MtCut/wbin) + 1 ;
      int bin2 = ( isQCD ) ? (MtCut/wbin) :         100/wbin ;
      cout<<" bin1 = "<< bin1 <<"   bin2 = "<<bin2<<endl;

      ostringstream pNameStr ;
      pNameStr << MtCut ;
      pNameStr << hID ;
      pNameStr << "_" ;
      pNameStr << njets ;
      TString pNameTStr = pNameStr.str() ;
      double bestNorm = Chi2Normalization( hData, hMC, bin1, bin2, 1, 1, pNameTStr, false );
      fprintf(ffile," %.3f , %.3f ,  %.3f , %d \n", cutV[k], evtNorm, bestNorm, njets );

      delete htt;
      delete htq;
      delete htw;
      delete hwj;
      delete hww;
      delete hzj;
      delete hqcd;
      delete hMC ;  
      delete hData ;
      delete hdt;
  }
  fclose(ffile);

}

//MC -> MC normalization
double BgEstimation::MeasureScale1D( vector<string>& fakeData, int hID, double MtCut, double METCut, bool doQCD, bool isVjNorm ){

  objInfo->ResetBGCuts(0, MtCut);
  objInfo->ResetBGCuts(1, METCut);

  int nbin_ = 2 ;

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  // get the damn pseduo-data
  int idx = 2 ;
  double ddx = fabs( idx*1.0 ) ;
  // change the W sample for systematic study
  vector<string> SysSample ;
  fitInput->GetParameters( "SysSamples", &SysSample );
  double uScaleW = (3.1*28049*0.2177) / 348887 ;
  double dScaleW = (3.1*28049*0.2177) / 1234550 ;
  double scaleW  = 1.;
  if ( SysSample[0].substr(3,1) == "+" ) {
      scaleW = fitInput->NormalizeComponents( "wj" ) ;
  } else {
      scaleW  = ( SysSample[0].substr(7,2) == "up" ) ? uScaleW : dScaleW ;
  }
  // create the pseudo data
  hObjs* hdt = new hObjs("dt", nbin_ ) ;

  objInfo->QCDSelector( fakeData[0], hdt, false, ddx*scale0, doQCD, 0, idx );
  objInfo->QCDSelector( fakeData[1], hdt, false, ddx*scale1, doQCD, 0, idx );
  objInfo->QCDSelector( fakeData[2], hdt, false, ddx*scale2, doQCD, 0, idx );
  objInfo->QCDSelector( fakeData[3], hdt, false, ddx*scale3, doQCD, 0, idx );
  objInfo->QCDSelector( fakeData[4], hdt, false, ddx*scale4, doQCD, 0, idx );
  objInfo->QCDSelector( fakeData[5], hdt, false, ddx*scale5, doQCD, 0, idx );
  objInfo->QCDSelector( fakeData[6], hdt, false, ddx*scale6, doQCD, 0, idx );

  vector<TH1D*> h_dt ;
  hdt->Fill1DVec( h_dt );

  if ( h_dt[13]->Integral() < 500. ) {
     nbin_ = 15 ;
     h_dt[13]->Rebin(2);
     h_dt[14]->Rebin(2);
     h_dt[15]->Rebin(2);
     cout<<" bin width changes -> "<< h_dt[13]->GetBinWidth(1) << endl;
  }

  objInfo->Reset(1, isVjNorm );  
  int scaleMode =  ( isVjNorm ) ? 1 : 0 ;

  // Get the MC distribution
  hObjs* htt = new hObjs("tt", nbin_ ) ;
  objInfo->QCDSelector( fakeData[0], htt, false, scale0, doQCD, scaleMode );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin_ ) ;
  //objInfo->QCDSelector( fakeData[1], hwj, false, scale1, doQCD, scaleMode );
  objInfo->QCDSelector( SysSample[0], hwj, false, scaleW, doQCD, scaleMode );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hzj = new hObjs("zj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[2], hzj, false, scale2, doQCD, scaleMode );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin_ ) ;
  objInfo->QCDSelector( fakeData[3], htq, false, scale3, doQCD, scaleMode );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin_ ) ;
  objInfo->QCDSelector( fakeData[4], htw, false, scale4, doQCD, scaleMode );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  hObjs* hww = new hObjs("ww", nbin_ ) ;
  objInfo->QCDSelector( fakeData[5], hww, false, scale5, doQCD, scaleMode );
  vector<TH1D*> h_ww ;
  hww->Fill1DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd", nbin_ ) ;
  objInfo->QCDSelector( fakeData[6], hqcd, false, scale6, doQCD, scaleMode );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  int nBin = h_dt[hID]->GetNbinsX() ;
  cout<<" N bin set = "<< nBin <<endl;
  TH1D* hMC   = new TH1D("hMC",   "MC   lep Mt ", nBin, 0, 150 );
  TH1D* hData = new TH1D("hData", "Data lep Mt ", nBin, 0, 150 );
  if ( doQCD ) { 
     hMC->Add( h_qcd[hID] );
     hData->Add( h_dt[hID],  1. );
     hData->Add( h_tt[hID], -1. );
     hData->Add( h_wj[hID], -1. );
     hData->Add( h_zj[hID], -1. );
     hData->Add( h_tq[hID], -1. );
     hData->Add( h_tw[hID], -1. );
     hData->Add( h_ww[hID], -1. );
  }
  if ( !doQCD ) { 
     hMC->Add( h_wj[hID] );
     hMC->Add( h_zj[hID] );
     hMC->Add( h_tt[hID] );
     hMC->Add( h_tq[hID] );
     hMC->Add( h_tw[hID] );
     hMC->Add( h_ww[hID] );
     hData->Add( h_dt[hID],  1. );
     hData->Add( h_qcd[hID], -1. );
  }

  cout<<" N of MC QCD = "<< h_qcd[hID]->Integral() <<endl;
  cout<<" N of MC WJ  = "<< h_wj[hID]->Integral() <<endl;
  cout<<" N of MC zJ  = "<< h_zj[hID]->Integral() <<endl;
  cout<<" N of MC tt  = "<< h_tt[hID]->Integral() <<endl;
  cout<<" N of MC tq  = "<< h_tq[hID]->Integral() <<endl;
  cout<<" N of MC tw  = "<< h_tw[hID]->Integral() <<endl;
  cout<<" N of MC ww  = "<< h_ww[hID]->Integral() <<endl;
  cout<<" N of Data = "<< h_dt[hID]->Integral() <<" -> "<< hData->Integral() <<endl ;

  cout<<" Event # normalization = "<< hData->Integral()/hMC->Integral() <<endl ;
  int wbin = 150/nBin ;
  int bin1 = ( doQCD ) ?         1    : (MtCut/wbin) + 1 ;
  int bin2 = ( doQCD ) ? (MtCut/wbin) :         100/wbin ;
  cout<<" bin1 = "<< bin1 <<"   bin2 = "<<bin2<<endl;

  ostringstream pNameStr ;
  pNameStr << MtCut ;
  pNameStr << hID ;
  TString pNameTStr = pNameStr.str() ;
  double bestNorm =  Chi2Normalization( hData, hMC, bin1, bin2, 1, scale2, pNameTStr, true );

  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hww;
  delete hzj;
  delete hqcd;
  delete hMC ;  
  delete hData ;

  return bestNorm ;
}

double BgEstimation::Chi2Normalization( TH1D* hData, TH1D* hMC, int Bx1, int Bx2, double s1, double s2, TString plotname, bool doFit ){

   //
   double Itgr1 = hData->Integral() ;
   double Itgr2 = hMC->Integral() ;
  
   double d_norm = Itgr1/Itgr2 ;
   double norm = d_norm - (0.34*d_norm) ;
   double step = (d_norm*0.34)/25 ;

   double minX2 = 99999. ;
   double bestNorm = d_norm ;

   //cout<<" delta Norm = "<< d_norm*0.34 <<"  Norm start from = "<< norm <<endl ;

   double nm[50];
   double chi2[50];
   for( int k=0; k<50; k++ ) {
      norm = norm + step ;
      // chi2 sum
      double X2 = 0 ;
      for (int i= Bx1; i<= Bx2; i++){
          double N1 = hData->GetBinContent(i);
	  double N2 = hMC->GetBinContent(i);
          if ( N1 < 0 ) N1 = 0 ;
	  if ( N2 <= 0. ) continue;
	  double x2_i = (N1 - (N2*norm) )*(N1 - (N2*norm) ) / (s1*N1 + s2*N2*norm*norm ) ;
          //cout<<"   N1:"<<N1<<" N2: "<<N2<<" x2i:"<<x2_i<<" X2:"<<X2<<endl;
	  X2 += x2_i ;
      }
      nm[k]   = norm ;
      chi2[k] = X2 ;
      //cout<<" Norm =  "<<nm[k] <<"   X2 = "<< chi2[k] <<"   step : "<< step <<endl; ;

      if ( k==0 ) minX2 = X2 ;
      if ( X2 < minX2 && k > 0 ) { 
         minX2 = X2 ;
         bestNorm = norm ;
      }
   }

   if ( doFit ) { 
      gStyle->SetOptFit(111);
      gStyle->SetOptTitle(0);
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd();

      TString theFolder = hfolder ;
      TString theSubFolder = "hFitting/" ;
      gSystem->mkdir( theFolder );
      gSystem->cd( theFolder );
      gSystem->mkdir( theSubFolder );
      gSystem->cd( "../" );

      TGraph* hNormX2 = new TGraph(50, nm, chi2 );

      double fitMin = ( nm[0] < nm[49] ) ? nm[0] : nm[49] ;
      double fitMax = ( fitMin == nm[0] ) ? nm[49] : nm[0] ;
      TF1 *func1 = new TF1("func1",MassFitFunction::fitParabola, fitMin, fitMax, 3 );
      //double nL = bestNorm - (bestNorm*0.1) ;
      //double nH = bestNorm + (bestNorm*0.1) ;
      //func1->SetParLimits(0, nL, nH );
      func1->FixParameter(0, bestNorm );
      func1->SetParLimits(1, 0.00001, 10. );
      func1->SetParLimits(2, minX2*0.9 , minX2*1.1 );

      hNormX2->SetMarkerSize(1);
      hNormX2->SetMarkerColor(4);
      hNormX2->SetMarkerStyle(21);
      hNormX2->Draw("AP");
      c1->Update();

      hNormX2->Fit( func1, "R", "sames", fitMin, fitMax );
      double fitNorm = func1->GetParameter(0);
      double fChi2 = func1->GetChisquare();
      int fNDF  = func1->GetNDF() ;
      double nChi2 = fChi2 / static_cast<double>( fNDF ) ;
      cout<<" The Fitted Normalization = "<< fitNorm <<"    the best norm = "<< bestNorm <<endl;
      if ( nChi2 < 5 ) bestNorm = fitNorm ;


      c1->Update();

      TString plotname1 = theFolder + theSubFolder +  "X2Fit" + plotname + "." + plotType ;
      c1->Print( plotname1 );
      delete c1;
      delete hNormX2 ;
      delete func1 ;
   }
   cout<<"  Final Normalization = "<< bestNorm << endl;
   return bestNorm ;

}


void BgEstimation::RatioScan( vector<string>& DataFiles, vector<string>& fNames, bool normMC ){

       // V+Jets
       vector<double> v10 = RatioXY( 1, 0, fNames,  1, false, false, false, normMC ) ;
       vector<double> v21 = RatioXY( 2, 1, fNames,  1, false, false, false, normMC ) ;
       vector<double> v32 = RatioXY( 3, 2, fNames,  1, false, false, false, normMC ) ;
       vector<double> v43 = RatioXY( 4, 3, fNames,  1, false, false, false, normMC ) ;
       vector<double> v54 = RatioXY( 5, 4, fNames,  1, false, false, false, normMC ) ;

       double vR[5] = { v10[0], v21[0], v32[0], v43[0], v54[0] }; 
       double evR[5] = { v10[5], v21[5], v32[5], v43[5], v54[5] }; 
       vector<double> v21c = RatioXY( 2, 1, fNames,  3, false, false, false, normMC ) ;
       vector<double> v21d = RatioXY( 2, 1, DataFiles,  3, false, false, false, normMC ) ;

       // QCD
       /*
       vector<double> q10 = RatioXY( 1, 0, fNames,  2, false, false, false, normMC ) ;
       vector<double> q21 = RatioXY( 2, 1, fNames,  2, false, false, false, normMC ) ;
       vector<double> q32 = RatioXY( 3, 2, fNames,  2, false, false, false, normMC ) ;
       vector<double> q43 = RatioXY( 4, 3, fNames,  2, false, false, false, normMC ) ;
       vector<double> q54 = RatioXY( 5, 4, fNames,  2, false, false, false, normMC ) ;

       double qR[5] = { q10[0], q21[0], q32[0], q43[0], q54[0] }; 
       double eqR[5] = { q10[5], q21[5], q32[5], q43[5], q54[5] }; 

       vector<double> q21c = RatioXY( 2, 1, fNames,  4, false, false, false, normMC ) ;
       vector<double> q21d = RatioXY( 2, 1, DataFiles, 4, false, false, false, normMC ) ;
       */

       //double x[5] = { 1, 2, 3, 4, 5} ;
       //double ex[5] = { 0. } ;

       gStyle->SetOptStat("");
       gStyle->SetOptTitle(0);
       gStyle->SetPadLeftMargin(0.15);
       gStyle->SetPadBottomMargin(0.15);
       gStyle->SetErrorX(0) ;

       TCanvas* c1 = new TCanvas("c1","", 800, 700);
       c1->SetGrid();
       c1->SetFillColor(10);
       c1->SetFillColor(10);
       c1->cd();

       TH1D* vjR = new TH1D("vjR"," VJets ratio", 5, 0.5, 5.5 );
       //TH1D* qjR = new TH1D("qjR"," QCD ratio  ", 5, 0.5, 5.5 );
       for ( int i=0; i<5; i++){
           vjR->Fill( i+1. , vR[i] );
           vjR->SetBinError( i+1, evR[i] );
           //qjR->Fill( i+1. , qR[i] );
           //qjR->SetBinError( i+1, eqR[i] );
       }
       
       vjR->SetAxisRange(0.01, 0.26, "Y");
       vjR->SetXTitle(" Ratio between different jet multiplicity");
       vjR->SetYTitle(" #frac{Number of N jets}{Number of N-1 jets} " );
       vjR->SetTitleOffset(2.0, "Y") ;
       vjR->SetTitleOffset(1.5, "X") ;

       vjR->SetLabelSize(0.05, "X");
       vjR->GetXaxis()->SetBinLabel(1, "R(1/0)");
       vjR->GetXaxis()->SetBinLabel(2, "R(2/1)");
       vjR->GetXaxis()->SetBinLabel(3, "R(3/2)");
       vjR->GetXaxis()->SetBinLabel(4, "R(4/3)");
       vjR->GetXaxis()->SetBinLabel(5, "R(5/4)");

       vjR->SetMarkerColor(4);
       vjR->SetMarkerSize(2);
       vjR->SetMarkerStyle(20);

       vjR->Draw("PE1");
       c1->Update();
       
       TF1 *fvR = new TF1("fvR", MassFitFunction::fitPoly, 0, 6, 1);
       fvR->FixParameter(0, v21c[0] ) ;
       fvR->SetLineColor(4);
       fvR->SetLineWidth(3);

       TF1 *dvR = new TF1("dvR", MassFitFunction::fitPoly, 0, 6, 1);
       dvR->FixParameter(0, v21d[0] ) ;
       dvR->SetLineColor(4);
       dvR->SetLineWidth(3);
       dvR->SetLineStyle(2);

       fvR->Draw("same");
       dvR->Draw("same");
        
       /*
       qjR->SetAxisRange(0.01, 0.26, "Y");
       qjR->SetXTitle(" Ratio between different jet multiplicity");
       qjR->SetYTitle(" #frac{Number of N jets}{Number of N-1 jets} " );
       qjR->SetTitleOffset(1.5, "Y") ;
       qjR->SetTitleOffset(1.5, "X") ;

       qjR->SetLabelSize(0.05, "X");
       qjR->GetXaxis()->SetBinLabel(1, "R(1/0)");
       qjR->GetXaxis()->SetBinLabel(2, "R(2/1)");
       qjR->GetXaxis()->SetBinLabel(3, "R(3/2)");
       qjR->GetXaxis()->SetBinLabel(4, "R(4/3)");
       qjR->GetXaxis()->SetBinLabel(5, "R(5/4)");

       qjR->SetMarkerColor(6);
       qjR->SetMarkerSize(2);
       qjR->SetMarkerStyle(21);

       TF1 *fqR = new TF1("fqR", MassFitFunction::fitPoly, 0, 6, 1);
       fqR->FixParameter(0, q21c[0] ) ;
       fqR->SetLineColor(6);
       fqR->SetLineWidth(3);

       TF1 *dqR = new TF1("dqR", MassFitFunction::fitPoly, 0, 6, 1);
       dqR->FixParameter(0, q21d[0] ) ;
       dqR->SetLineColor(2);
       dqR->SetLineWidth(3);
       dqR->SetLineStyle(2);
       
       qjR->Draw("PE1");       
       c1->Update();

       fqR->Draw("same");
       dqR->Draw("same");
       */

       TLegend *leg = new TLegend(.18, .16, .65, .50 );
       leg->AddEntry(vjR,  "V+Jets , single top",  "P");
       leg->AddEntry(fvR,  "R2/1 from control region - MC", "L");
       leg->AddEntry(dvR,  "R2/1 from contral region - Data", "L");
       /*
       leg->AddEntry(qjR,  "QCD",     "P");
       leg->AddEntry(fqR,  "R2/1 from contral region - MC", "L");
       leg->AddEntry(dqR,  "R2/1 from contral region - Data", "L");
        */
       leg->Draw("same");

       c1->Update();

       TString plotname = hfolder +"Ratios."+ plotType ;
       c1->Print(plotname);
       delete c1;
}
