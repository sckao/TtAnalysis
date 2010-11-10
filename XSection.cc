#include "XSection.h"
#include "WFormat.h"

XSection::XSection(){

  fitInput  = new MassAnaInput();
  wmfitter  = new HadWMassFitter();
  pseudoExp = new PseudoExp();
  bgEst     = new BgEstimation();
  objInfo   = new ObjectInfo();


  fitInput->GetParameters( "4JSamples", &File4J );
  fitInput->GetParameters( "2JSamples", &File2J );
  fitInput->GetParameters( "FakeData",  &File0J );
  fitInput->GetParameters( "TheData",   &dataFile );

  string phaseSmear ;
  fitInput->GetParameters( "PhaseSmear", &phaseSmear );
  smearing = ( phaseSmear == "ON" ) ? true : false ;

  fitInput->GetParameters( "Path", &hfolder );
  theFolder = hfolder ;
  gSystem->mkdir( theFolder );

  fitInput->GetParameters( "InputMean4J", &inputMean4J );
  fitInput->GetParameters( "InputMean2J", &inputMean2J );

  fitInput->GetParameters( "PlotType", &plotType );

  // normalize the MC for signal and background counter
  string mcNormalization ;
  fitInput->GetParameters( "MCNormalization", &mcNormalization );
  mcNorm = ( mcNormalization == "YES" ) ? true : false ;
  wmfitter->SetMCNormalization( mcNorm ) ;
  objInfo->Reset(1, mcNorm ) ;   
}

XSection::~XSection(){

  delete fitInput ;
  delete wmfitter ;
  delete pseudoExp ;
  delete bgEst ;
  delete objInfo ;

}

// Exercise for background estimation

vector<double> XSection::RunEnsembles( int nJets, int tsz, double pMean, int nRun, int RandomSeed, string fileName, TH1D* hGen ){

  cout<<"  <<< Getting Ensembles for "<< nRun <<" >>> "<<endl;

  // set up the random number function
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( RandomSeed );

  int  nRun_ = 0 ;
  vector<double> nCountV ;
  vector<int> ensembles ;
  int nRepeat = 0 ;
  while ( nRun_ < nRun) {

       // shuffle the events
       vector<int> evtline = pseudoExp->EventShuffle( tsz, RandomSeed );
       //if ( nRepeat > 0 ) cout<<"   smearing samples "<<endl;
       if ( nRepeat >= 0 ) smearing = true  ;
       bool nextRun = true ;
       int  RunStop = 0 ;
       int  NEvts   = 0 ;
       int evtSize = static_cast<int>( evtline.size() ) ;
       for (int k=0; k < evtSize; k++) {
 
           if ( nRun_ == nRun ) break; 
           if ( nextRun ) {
	      NEvts = tRan->Poisson( pMean );
	      nextRun = false ;
	      RunStop = k + NEvts ;
	      ensembles.clear() ;
              if ( RunStop >= evtSize ) break ;
	      cout<<"  Run["<<nRun_<<"] " ; 
              hGen->Fill( NEvts );
           }
           if ( k < RunStop ) {
              ensembles.push_back( evtline[k] );
           }
           if ( k == (RunStop - 1) || RunStop == 0 ) {
              ACounter* counter = new ACounter();
	      wmfitter->ReFitSolution( fileName, counter , nJets, 1, &ensembles, 0, smearing );
	      vector<double> nCount = counter->Output();
	      nCountV.push_back( nCount[0] );
	      nextRun = true ;
	      nRun_++ ;
	      cout<<" ("<<k<<") N of Evts needed For this Run = "<< NEvts <<" passed # = "<< nCount[0] <<endl ;
              delete counter;
           }
       }
       nRepeat++ ;
  }
  cout<<" Number of Repeat = "<< nRepeat <<endl ;

  delete tRan ;
  return nCountV ;
}

vector<double> XSection::CutEff( int nj ){

  vector<double> nEvts ;
  fitInput->GetParameters( "nEvents" , &nEvts );

  double scaleTt  = fitInput->NormalizeComponents( "tt" );
  double scaleWJ  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD = fitInput->NormalizeComponents( "qcd" );
  double scaleZJ  = fitInput->NormalizeComponents( "zj" );
  double scaleTq  = fitInput->NormalizeComponents( "tq" );
  double scaleTw  = fitInput->NormalizeComponents( "tw" );

  // reset the cuts -- use the default value from DataCard
  // Without cuts
  wmfitter->ResetCuts( 0, 999, 0, 999, 0, 999 );

  ACounter* count4jSg = new ACounter();
  wmfitter->ReFitSolution( File0J[0], count4jSg, nj, scaleTt, NULL, 0, smearing );
  vector<double> n4JSg= count4jSg->Output();

  ACounter* count4jBg = new ACounter();
  wmfitter->ReFitSolution( File0J[1], count4jBg, nj, scaleWJ,  NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[2], count4jBg, nj, scaleQCD, NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[3], count4jBg, nj, scaleZJ, NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[4], count4jBg, nj, scaleTq, NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[5], count4jBg, nj, scaleTw, NULL, 0, smearing );
  vector<double> n4JBg = count4jBg->Output();

  // With cuts
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );

  ACounter* count4jSg_c = new ACounter();
  wmfitter->ReFitSolution( File0J[0], count4jSg_c, nj, scaleTt, NULL, 0, smearing );
  vector<double> n4JSg_cuts = count4jSg_c->Output();

  ACounter* count4jBg_c = new ACounter();
  wmfitter->ReFitSolution( File0J[1], count4jBg_c, nj, scaleWJ,  NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[2], count4jBg_c, nj, scaleQCD, NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[3], count4jBg_c, nj, scaleZJ, NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[4], count4jBg_c, nj, scaleTq, NULL, 0, smearing );
  wmfitter->ReFitSolution( File0J[5], count4jBg_c, nj, scaleTw, NULL, 0, smearing );
  vector<double> n4JBg_cuts = count4jBg_c->Output();

  // signal efficiency with cuts
  double  Eff = n4JSg_cuts[0] / (nEvts[0]*scaleTt) ;
  double sEff = MassFitFunction::ErrAovB( n4JSg_cuts[0], nEvts[0], n4JSg_cuts[1], -1 ) / scaleTt ;

  // 4J background cut efficiency 
  double Eff_Bg  =  n4JBg_cuts[0] / n4JBg[0] ;
  double sEff_Bg = MassFitFunction::ErrAovB( n4JBg_cuts[0], n4JBg[0], n4JBg_cuts[1], n4JBg[1] );

  cout<<" ===== Estimate Cut Efficiency ===== "<<endl;
  cout<<" MC Signal  "<< nj <<"j = "<< n4JSg[0] <<" Eff = "<< Eff <<" +/- "<< sEff <<  endl;
  cout<<" MC Background "<< nj <<"j = "<< n4JBg[0] <<" Eff = "<< Eff_Bg  <<" +/- "<< sEff_Bg <<endl;
  cout<<" ------------------------------------- "<<endl;
  delete count4jSg ;
  delete count4jBg ;
  delete count4jSg_c ;
  delete count4jBg_c ;

  vector<double> Effv ;
  Effv.push_back(  Eff );
  Effv.push_back(  Eff_Bg );
  Effv.push_back( sEff );
  Effv.push_back( sEff_Bg );
  return Effv ;

}

// calculate cross-section by estimating background in different W pt region
vector<double> XSection::CrossSection( double nData_4J, vector<double>& nData_2J, vector<double>& R42, vector<double>& EffCut ){

  double lumi ;
  fitInput->GetParameters( "Lumi" , &lumi );
  vector<double> EffHLT ;
  fitInput->GetParameters( "EffHLT" , &EffHLT );

  vector<double> nBg_4J = bgEst->BgEstimate( R42, nData_2J ) ;

  double xsecTt = ( nData_4J - (nBg_4J[0]*EffCut[1]) ) / ( lumi*EffHLT[0]*EffCut[0] );

  // calculate the uncertainty of cross-section 
  vector<double> s_n4J = MassFitFunction::StatErr( nData_4J );
  double sn_Bg = MassFitFunction::ErrAxB( nBg_4J[0], EffCut[1],  nBg_4J[2], EffCut[3] ) ;
  double sp_Bg = MassFitFunction::ErrAxB( nBg_4J[0], EffCut[1],  nBg_4J[1], EffCut[3] ) ;

  double sn_Tt = sqrt( (s_n4J[0]*s_n4J[0]) + ( sn_Bg*sn_Bg ) );
  double sp_Tt = sqrt( (s_n4J[1]*s_n4J[1]) + ( sp_Bg*sp_Bg ) );
  double n_Tt = nData_4J - (nBg_4J[0]*EffCut[1]) ;

  double sn_xsec0 = MassFitFunction::ErrAovB( n_Tt, EffCut[0], sn_Tt, EffCut[2] );
  double sp_xsec0 = MassFitFunction::ErrAovB( n_Tt, EffCut[0], sp_Tt, EffCut[2] );
  double sn_xsec = sn_xsec0/(lumi*EffHLT[0]) ;
  double sp_xsec = sp_xsec0/(lumi*EffHLT[0]) ;

  cout<<" ========== Cross-Section Report ============ "<<endl;
  cout<<" N of All 4J = "<< nData_4J <<endl;
  cout<<" N of All 2J = "<<  nData_2J[0] ;
  cout<<" -> "<< nData_2J[1]<<","<<nData_2J[2]<<","<<nData_2J[3]<<","<<nData_2J[4]<<","<<nData_2J[5]<<endl;
  cout<<" N of Tt = "<< n_Tt <<" - "<< sn_Tt<<" + "<< sp_Tt <<endl;
  cout<<" N of Bg "<< nBg_4J[0] <<" Eff= "<<EffCut[1] <<" - "<<sn_Bg<<" + "<< sp_Bg << endl;
  cout<<" Tt Xsec = "<< xsecTt <<" - "<< sn_xsec <<" + "<< sp_xsec << endl;
  cout<<" ============================================ "<<endl;

  vector<double> xVs ;
  xVs.push_back( xsecTt );                
  xVs.push_back( sn_xsec );
  xVs.push_back( sp_xsec );
  xVs.push_back( nBg_4J[0]*EffCut[1] );
  return xVs ;

}

// calculate cross-section by estimating background in one phase space
vector<double> XSection::CrossSection( double nData_4J, double nData_2J, vector<double>& R42, vector<double>& EffCut ){

  double lumi ;
  fitInput->GetParameters( "Lumi" , &lumi );
  vector<double> EffHLT ;
  fitInput->GetParameters( "EffHLT" , &EffHLT );

  vector<double> nBg_4J = bgEst->BgEstimate( R42, nData_2J ) ;

  double xsecTt = ( nData_4J - (nBg_4J[0]*EffCut[1]) ) / ( lumi*EffHLT[0]*EffCut[0] );

  // calculate the uncertainty of cross-section 
  vector<double> s_n4J = MassFitFunction::StatErr( nData_4J );
  double sn_Bg = MassFitFunction::ErrAxB( nBg_4J[0], EffCut[1],  nBg_4J[2], EffCut[3] ) ;
  double sp_Bg = MassFitFunction::ErrAxB( nBg_4J[0], EffCut[1],  nBg_4J[1], EffCut[3] ) ;

  double sn_Tt = sqrt( (s_n4J[0]*s_n4J[0]) + ( sn_Bg*sn_Bg ) );
  double sp_Tt = sqrt( (s_n4J[1]*s_n4J[1]) + ( sp_Bg*sp_Bg ) );
  double n_Tt = nData_4J - (nBg_4J[0]*EffCut[1]) ;

  double sn_xsec0 = MassFitFunction::ErrAovB( n_Tt, EffCut[0], sn_Tt, EffCut[2] );
  double sp_xsec0 = MassFitFunction::ErrAovB( n_Tt, EffCut[0], sp_Tt, EffCut[2] );
  double sn_xsec = sn_xsec0/(lumi*EffHLT[0]) ;
  double sp_xsec = sp_xsec0/(lumi*EffHLT[0]) ;

  cout<<" ========== Cross-Section Report ============ "<<endl;
  cout<<" N of All 4J = "<< nData_4J <<" N of All 2J = "<<  nData_2J <<endl;
  cout<<" N of Tt = "<< n_Tt <<" - "<< sn_Tt<<" + "<< sp_Tt <<endl;
  cout<<" N of Bg "<< nBg_4J[0]*EffCut[1] <<" - "<<sn_Bg<<" + "<< sp_Bg <<endl;
  cout<<" Tt Xsec = "<< xsecTt <<" - "<< sn_xsec <<" + "<< sp_xsec << endl;
  cout<<" ============================================ "<<endl;

  vector<double> xVs ;
  xVs.push_back( xsecTt );
  xVs.push_back( sn_xsec );
  xVs.push_back( sp_xsec );
  xVs.push_back( nBg_4J[0]*EffCut[1] );
  return xVs ;

}

void XSection::MethodTest1( int nj ){

  string useDefault_ ;
  fitInput->GetParameters( "UseDefault", &useDefault_ ); 
  bool useDefault = ( useDefault_ == "True" ) ? true : false ;

  vector<double> R42 ;
  fitInput->GetParameters( "Ratio42", &R42 );
  vector<double> EffCut;
  if ( useDefault ) {
     fitInput->GetParameters( "EffCut",  &EffCut );
     cout<<" R42 size = "<< R42.size() <<"  Eff size = "<< EffCut.size() <<endl;
  } else {
     EffCut = CutEff( nj );
  }

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ0  = fitInput->NormalizeComponents( "zj" );
  double scaleTq0  = fitInput->NormalizeComponents( "tq" );
  double scaleTW0  = fitInput->NormalizeComponents( "tw" );

  int idx = -2 ;
  double scaleTt   = scaleTt0 * 2. ;
  double scaleWJ   = scaleWJ0 * 2. ;
  double scaleQCD  = scaleQCD0* 2. ;
  double scaleZJ   = scaleZJ0 * 2. ;
  double scaleTq   = scaleTq0 * 2. ;
  double scaleTW   = scaleTW0 * 2. ;

  // run the execise - use the cuts in DataCard
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );

  ACounter* count4j_sg = new ACounter();
  wmfitter->ReFitSolution( File0J[0], count4j_sg, nj, scaleTt, NULL, idx, smearing );
  vector<double> n4J_sg= count4j_sg->Output();

  ACounter* count4j_bg = new ACounter();
  wmfitter->ReFitSolution( File0J[1], count4j_bg, nj, scaleWJ,  NULL, idx, smearing );
  wmfitter->ReFitSolution( File0J[2], count4j_bg, nj, scaleQCD, NULL, idx, smearing );
  wmfitter->ReFitSolution( File0J[3], count4j_bg, nj, scaleZJ,  NULL, idx, smearing );
  wmfitter->ReFitSolution( File0J[4], count4j_bg, nj, scaleTq,  NULL, idx, smearing );
  wmfitter->ReFitSolution( File0J[5], count4j_bg, nj, scaleTW,  NULL, idx, smearing );
  vector<double> n4J_bg = count4j_bg->Output();
  double totalN4J = n4J_sg[0] + n4J_bg[0] ;

  // method for different W pt space ...
  objInfo->Reset(0, 2) ;      // set 2 jets
  objInfo->Reset(0, false) ;  // set exclusive 
  bgCounter* count2j_bg = new bgCounter("eff") ;
  objInfo->EvtSelector( File0J[0], count2j_bg, smearing, scaleTt,  NULL, idx );
  objInfo->EvtSelector( File0J[1], count2j_bg, smearing, scaleWJ,  NULL, idx );
  objInfo->EvtSelector( File0J[2], count2j_bg, smearing, scaleQCD, NULL, idx );
  objInfo->EvtSelector( File0J[3], count2j_bg, smearing, scaleZJ,  NULL, idx );
  objInfo->EvtSelector( File0J[4], count2j_bg, smearing, scaleTq,  NULL, idx );
  objInfo->EvtSelector( File0J[5], count2j_bg, smearing, scaleTW,  NULL, idx );
  vector<double> n2J_bg ;
  count2j_bg->CounterVec( n2J_bg );

  cout<<" ===== Method Test ===== "<<endl;
  cout<<"  Generated 4J Tt = "<< n4J_sg[0] <<" Bg = "<<n4J_bg[0] <<endl;

  vector<double> xtt = CrossSection(  totalN4J, n2J_bg[0], R42, EffCut ) ;
 
  delete count4j_sg ;
  delete count4j_bg ;
  delete count2j_bg ;

}

// ** turn off the MCNormalization
void XSection::MethodTest2( string cfgFile, int nj ){

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ0  = fitInput->NormalizeComponents( "zj" );
  double scaleTq0  = fitInput->NormalizeComponents( "tq" );
  double scaleTW0  = fitInput->NormalizeComponents( "tw" );

  double lumi ;
  fitInput->GetParameters( "Lumi" , &lumi );
  vector<double> EffHLT ;
  fitInput->GetParameters( "EffHLT" , &EffHLT );
  vector<int> testPara;
  fitInput->GetParameters( "TestPara", &testPara, cfgFile ); 

  // scan muon pt and different cuts
  int sysType         = testPara[0] ;
  const static int sz = testPara[1] ;
  int centralValue    = testPara[2] ;
  int cut_topo        = testPara[3] ;

  string sfx = "" ;
  if ( sysType == 0 ) sfx = "_" ;
  if ( sysType == 1 ) sfx = "_i" ;
  if ( sysType == 2 ) sfx = "_j" ;
  if ( sysType == 3 ) sfx = "_q" ;
  TString gNames[4] = { "muon Pt", "muon Isolation", "JES %", " q^{2} scaling of W " } ;
  double  varScal[4] = { 1, 0.001, 1, 1 };

  double xV[sz], xErrn[sz], xErrp[sz], vars[sz], varsErr[sz] ;
  double nAll[sz], nBg[sz], Effs[sz]  ;
  ACounter*  count4j_sg = new ACounter();
  ACounter*  count4j_bg = new ACounter();
  bgCounter* count2j_bg = new bgCounter("eff") ;

  for (int k=0; k<3; k++ ) {
 
      int cut_idx = ( cut_topo == 9 ) ? k : cut_topo ;
      if ( k != cut_idx ) continue ;
      ostringstream cutIdStr ;
      cutIdStr << "_" ;
      cutIdStr << cut_idx ;

      for (int i=0; i<sz; i++ ) {
 
          double muonPt  = ( sysType == 0 ) ? (16.0+(i*2))   : 20. ;
	  double muonIso = ( sysType == 1 ) ? (0.05+(i*0.025)) : 0.1 ;
	  int var = 0 ;
	  if ( sysType == 0 ) var = 16 + (i*2) ;
	  if ( sysType == 1 ) var = 50 + (i*25) ;
	  if ( sysType == 2 ) var = 95 + (i*5) ;
	  if ( sysType == 3 ) var = 1 + i ;

	  ostringstream normStr ;
	  normStr << "Norm" ;
	  normStr << sfx ;
	  normStr << var ;

	  ostringstream ratioStr ;
	  ratioStr << "Ratio42" ;
	  ratioStr << sfx ;
	  ratioStr << var ;

	  ostringstream cutStr ;
	  cutStr << "TopoCut" ;
	  cutStr << cut_idx ;

	  ostringstream effStr ;
	  effStr << "Eff" ;
	  effStr << cut_idx ;
	  effStr << sfx ;
	  effStr << var ;

	  vector<double> normV ;
	  fitInput->GetParameters( normStr.str(), &normV, cfgFile ); 
	  vector<double> R42 ;
	  fitInput->GetParameters( ratioStr.str(), &R42, cfgFile );
	  vector<double> topoCuts;
	  fitInput->GetParameters( cutStr.str(),  &topoCuts, cfgFile );
	  vector<double> EffCuts;
	  fitInput->GetParameters( effStr.str(),  &EffCuts, cfgFile );

	  cout<<" *****  Systematic test "<< var <<" *********** "<<endl ;
	  cout<<" Norm(V,Q) = ( "<< normV[1] <<" , "<<normV[3]<<" )"<< endl;
	  cout<<" Ratio42   = ( "<< R42[0]<<" , "<< R42[5] <<" )"<< endl;
	  cout<<" Cut Eff   = ( "<< EffCuts[0]<<", "<< EffCuts[1]<<", "<< EffCuts[2]<<", "<< EffCuts[3] <<" )"<< endl;
	  cout<<" Cut Vaule = ( "<< topoCuts[0]<<", "<< topoCuts[1]<<", "<< topoCuts[2]<<", "<< topoCuts[3]<<", ";
	  cout<< topoCuts[4]<<", "<< topoCuts[5] <<" )"<< endl;

	  int idx = -2 ;
	  double scaleTt   = scaleTt0 * 2. ;
	  double scaleWJ   = scaleWJ0 * 2. * normV[1] ;
	  double scaleZJ   = scaleZJ0 * 2. * normV[1] ;
	  double scaleQCD  = scaleQCD0* 2. * normV[3] ;
	  double scaleTq   = scaleTq0 * 2. ;
	  double scaleTW   = scaleTW0 * 2. ;

	  // run the execise - use the cuts in DataCard
	  wmfitter->ResetCuts( topoCuts[0], topoCuts[1], topoCuts[2], topoCuts[3], topoCuts[4], topoCuts[5] );
	  wmfitter->SetMuonCuts( muonPt, -1, muonIso );

	  wmfitter->ReFitSolution( File0J[0], count4j_sg, nj, scaleTt, NULL, idx, smearing );
	  vector<double> n4J_sg= count4j_sg->Output();

	  wmfitter->ReFitSolution( File0J[1], count4j_bg, nj, scaleWJ,  NULL, idx, smearing );
	  wmfitter->ReFitSolution( File0J[2], count4j_bg, nj, scaleQCD, NULL, idx, smearing );
	  wmfitter->ReFitSolution( File0J[3], count4j_bg, nj, scaleZJ,  NULL, idx, smearing );
	  wmfitter->ReFitSolution( File0J[4], count4j_bg, nj, scaleTq,  NULL, idx, smearing );
	  wmfitter->ReFitSolution( File0J[5], count4j_bg, nj, scaleTW,  NULL, idx, smearing );
	  vector<double> n4J_bg = count4j_bg->Output();
	  double totalN4J = n4J_sg[0] + n4J_bg[0] ;

	  // method for different W pt space ...
	  objInfo->Reset(0, 2) ;      // set 2 jets
	  objInfo->Reset(0, false) ;  // set exclusive 
	  objInfo->Reset(2, muonPt);  // set muon pt
	  objInfo->Reset(4, muonIso);  // set muon isolation
	  objInfo->EvtSelector( File0J[0], count2j_bg, smearing, scaleTt,  NULL, idx );
	  objInfo->EvtSelector( File0J[1], count2j_bg, smearing, scaleWJ,  NULL, idx );
	  objInfo->EvtSelector( File0J[2], count2j_bg, smearing, scaleQCD, NULL, idx );
	  objInfo->EvtSelector( File0J[3], count2j_bg, smearing, scaleZJ,  NULL, idx );
	  objInfo->EvtSelector( File0J[4], count2j_bg, smearing, scaleTq,  NULL, idx );
	  objInfo->EvtSelector( File0J[5], count2j_bg, smearing, scaleTW,  NULL, idx );
	  vector<double> n2J_bg ;
	  count2j_bg->CounterVec( n2J_bg );

	  cout<<" ===== Method Test ===== "<<endl;
	  cout<<"  Generated 4J Tt = "<< n4J_sg[0] <<" Bg = "<<n4J_bg[0] <<endl;

	  vector<double> xtt = CrossSection(  totalN4J, n2J_bg[0], R42, EffCuts ) ;

	  vars[i]    = var*varScal[sysType] ; 
	  varsErr[i] = 0.0 ; 
	  xV[i]     = xtt[0] ;
	  xErrn[i]  = xtt[1] ;
	  xErrp[i]  = xtt[2] ;
	  nBg[i]    = xtt[3]  ;
	  nAll[i]   = totalN4J ;
	  Effs[i]   = EffCuts[0] ;

	  count4j_sg->Reset();
	  count4j_bg->Reset();
	  count2j_bg->Reset();
      }
      // calculate the percent error
      double dSg[sz] ;
      for (int i=0; i<sz; i++ ) {
          double d_nAll = nAll[i] - nAll[centralValue] ;
	  double d_nBg  = nBg[i]  - nBg[centralValue]  ;
	  dSg[i] = sqrt( (d_nAll*d_nAll) + (d_nBg*d_nBg) ) / ( lumi*EffHLT[0]*Effs[i]);
          if ( sysType > 1 ) dSg[i] = fabs( xV[i] - xV[ centralValue ] );
	  double percentErr = (xV[i] - xV[ centralValue ]) / xV[ centralValue ]  ;
	  cout<<" dX "<< xV[i] - xV[ centralValue ] << " %Err "<< percentErr ;
          cout<<"   dXsec ="<< dSg[i] <<"  Err% = "<< dSg[i]/xV[centralValue]  << endl;
	  if ( dSg[i]/xV[i] > 1. )  dSg[i] = xV[i] ;
      }

      TString theFolder = hfolder ;
      gSystem->cd( theFolder );
      gSystem->mkdir( "Systematic"  );
      gSystem->cd("../");

      gStyle->SetOptTitle(0);
      gStyle->SetOptStat("");
      gStyle->SetLabelSize( 0.05, "X");
      gStyle->SetLabelSize( 0.05, "Y");
      gStyle->SetPadLeftMargin(0.20);
      gStyle->SetPadBottomMargin(0.15);

      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd();

      TGraphAsymmErrors* gXs = new TGraphAsymmErrors( sz, vars, xV,  varsErr, varsErr, xErrn, xErrp );
      gXs->SetMarkerColor(4);
      gXs->SetMarkerStyle(21);
      gXs->SetMarkerSize(2.0);
      gXs->GetYaxis()->SetTitle("#sigma  #pm #Delta#sigma_{stat} (pb^{-1})" );
      gXs->GetYaxis()->SetTitleOffset(2.0);
      gXs->GetXaxis()->SetTitle( gNames[sysType] );
      gXs->GetXaxis()->SetTitleOffset(1.7);
      gXs->Draw("AP");
      c1->Update();

      TString plotname1 = hfolder + "Systematic/Xsection"+ sfx + "Sys"+ cutIdStr.str() + "." +plotType ;
      c1->Print( plotname1 );

      TCanvas* c2 = new TCanvas("c2","", 800, 600);
      c2->SetGrid();
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      c2->cd();

      TGraphErrors* gdX = new TGraphErrors( sz, vars, xV, varsErr, dSg );
      gdX->SetMarkerColor(4);
      gdX->SetMarkerStyle(21);
      gdX->SetMarkerSize(2.0);
      gdX->GetYaxis()->SetTitle("#sigma  #pm #Delta#sigma_{stat} (pb^{-1})" );
      gdX->GetYaxis()->SetTitleOffset(2.0);
      gdX->GetXaxis()->SetTitleOffset(1.7);
      gdX->GetXaxis()->SetTitle( gNames[sysType] );
      gdX->Draw("AP");
      c2->Update();

      TString plotname2 = hfolder + "Systematic/dXsec" + sfx + "Sys"+ cutIdStr.str() + "." +plotType ;
      c2->Print( plotname2 );

      delete c1;
      delete c2;
      delete gXs;
      delete gdX;
  }
  delete count4j_sg ;
  delete count4j_bg ;
  delete count2j_bg ;

	
}

void XSection::MethodTest3( int nj ){

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );

  // scan muon pt and different cuts
  int idx = -2 ;
  for (int i=0; i<3; i++ ) {

     ACounter* count4j_sg = new ACounter();
     wmfitter->ReFitSolution( File0J[0], count4j_sg, nj, scaleTt0, NULL, idx, smearing );
     ACounter* count4j_bg = new ACounter();
     wmfitter->ReFitSolution( File0J[1], count4j_sg, nj, scaleWJ0, NULL, idx, smearing );
     vector<double> n4J_sg= count4j_sg->Output();

     delete count4j_sg ;
     delete count4j_bg ;
     cout<<" n4J = "<< n4J_sg[0] <<endl;
  }
  for (int i=0; i<3; i++ ) {
     bgCounter* count2j_bg = new bgCounter("eff") ;
     objInfo->Reset(0, 4) ;      // set 4 jets
     objInfo->Reset(0, true) ;  // set exclusive 
     objInfo->EvtSelector( File0J[0], count2j_bg, smearing, scaleTt0,  NULL, idx );
     delete count2j_bg ;
  }
  for (int i=0; i<3; i++ ) {

     ACounter* count4j_sg = new ACounter();
     wmfitter->ReFitSolution( File0J[0], count4j_sg, nj, scaleTt0, NULL, idx, smearing );
     ACounter* count4j_bg = new ACounter();
     wmfitter->ReFitSolution( File0J[1], count4j_sg, nj, scaleWJ0, NULL, idx, smearing );
     vector<double> n4J_sg= count4j_sg->Output();

     delete count4j_sg ;
     delete count4j_bg ;
     cout<<" n4J = "<< n4J_sg[0] <<endl;
  }

}

void XSection::BgClosureTest( int nX, int nY, bool inclX, bool inclY ){

  // get the Ratio42
  vector<double> RXY = bgEst->RatioXY( nX, nY, File0J, 0, false, inclX, inclY, true  ) ;

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ0  = fitInput->NormalizeComponents( "zj" );
  double scaleTq0  = fitInput->NormalizeComponents( "tq" );
  double scaleTw0  = fitInput->NormalizeComponents( "tw" );

  int idx = -2 ;
  double scaleTt   = scaleTt0* 2.;
  double scaleWJ   = scaleWJ0* 2. ;
  double scaleQCD  = scaleQCD0* 2.  ;
  double scaleZJ   = scaleZJ0* 2. ;
  double scaleTq   = scaleTq0* 2. ;
  double scaleTw   = scaleTw0* 2. ;
  
  // actural XJ background
  objInfo->Reset(0, nX ) ;
  objInfo->Reset(0, inclX ) ;
  bgCounter* countxj_bg = new bgCounter("bgxj") ;
  //objInfo->EvtSelector( File0J[0], countxj_bg, smearing, scaleTt, NULL, idx );
  objInfo->EvtSelector( File0J[1], countxj_bg, smearing, scaleWJ,  NULL, idx );
  objInfo->EvtSelector( File0J[2], countxj_bg, smearing, scaleQCD, NULL, idx );
  objInfo->EvtSelector( File0J[3], countxj_bg, smearing, scaleZJ,  NULL, idx );
  objInfo->EvtSelector( File0J[4], countxj_bg, smearing, scaleTq,  NULL, idx );
  objInfo->EvtSelector( File0J[5], countxj_bg, smearing, scaleTw,  NULL, idx );
  vector<double> nXJ_bg ;
  countxj_bg->CounterVec( nXJ_bg );

  // Input YJ background
  objInfo->Reset(0, nY ) ;
  objInfo->Reset(0, inclY ) ;
  bgCounter* countyj_bg = new bgCounter("bgyj") ;
  objInfo->EvtSelector( File0J[0], countyj_bg, smearing, scaleTt,  NULL, idx );
  objInfo->EvtSelector( File0J[1], countyj_bg, smearing, scaleWJ,  NULL, idx );
  objInfo->EvtSelector( File0J[2], countyj_bg, smearing, scaleQCD, NULL, idx );
  objInfo->EvtSelector( File0J[3], countyj_bg, smearing, scaleZJ,  NULL, idx );
  objInfo->EvtSelector( File0J[4], countyj_bg, smearing, scaleTq,  NULL, idx );
  objInfo->EvtSelector( File0J[5], countyj_bg, smearing, scaleTw,  NULL, idx );
  vector<double> nYJ_bg ;
  countyj_bg->CounterVec( nYJ_bg );

  vector<double> expBg4J  = bgEst->BgEstimate( RXY, nYJ_bg ) ;
  vector<double> expBg4J0 = bgEst->BgEstimate( RXY, nYJ_bg[0] ) ;

  cout<<" ===== Background Closure Test ===== "<< endl;
  cout<<"  Generated "<<nX<<" Background = "<< nXJ_bg[0] ;
  cout<<" ( "<< nXJ_bg[1] <<" , "<< nXJ_bg[2] <<" , "<< nXJ_bg[3] <<" , "<< nXJ_bg[4] <<" )"<< endl;

  cout<<"  Measured "<<nX<<" Background0 = "<< expBg4J0[0] <<" + "<< expBg4J0[1] <<" - "<< expBg4J0[2] << endl; 

  cout<<"  Measured "<<nX<<" BG in different Pt = "<< expBg4J[0] <<" + "<< expBg4J[1] <<" - "<< expBg4J[2] << endl;
  cout<<" ( "<< expBg4J[3] <<" , "<< expBg4J[4] <<" , "<< expBg4J[5] <<" , "<< expBg4J[6] <<" )"<< endl;
  cout<<" =================================== "<< endl;

  delete countxj_bg ;
  delete countyj_bg ;

}

void XSection::EnsembleTest1( int nRun, int randomSeed ){

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  // 4J Samples
  int tsz0 = fitInput->TreeSize( File4J[0] );
  int tsz1 = fitInput->TreeSize( File4J[1] );
  int tsz2 = fitInput->TreeSize( File4J[2] );
  int tsz3 = fitInput->TreeSize( File4J[3] );
  int tsz4 = fitInput->TreeSize( File4J[4] );
  int tsz5 = fitInput->TreeSize( File4J[5] );

  TH1D* gTt4J = new TH1D("gTt4J", " Tt 4J         ",  25, 0, 50 );
  TH1D* gBg4J = new TH1D("gBg4J", " Background 4J ",  15, 0, 30 );

  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  vector<double> n4JTt = RunEnsembles( 4, tsz0, inputMean4J[0],  nRun, randomSeed, File4J[0], gTt4J );
  vector<double> nWj   = RunEnsembles( 4, tsz1, inputMean4J[1],  nRun, randomSeed, File4J[1], gBg4J );
  vector<double> nQCD  = RunEnsembles( 4, tsz2, inputMean4J[2],  nRun, randomSeed, File4J[2], gBg4J );
  vector<double> nZj   = RunEnsembles( 4, tsz3, inputMean4J[3],  nRun, randomSeed, File4J[3], gBg4J );
  vector<double> nTq   = RunEnsembles( 4, tsz4, inputMean4J[4],  nRun, randomSeed, File4J[4], gBg4J );
  vector<double> nTw   = RunEnsembles( 4, tsz5, inputMean4J[5],  nRun, randomSeed, File4J[5], gBg4J );
  cout<<"  size of 4J tt:"<< tsz0 <<"  WJ:"<< tsz1 <<"  QCD:"<< tsz2 <<" ZJ: "<< tsz3<<endl;

  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  // 2J Samples
  int t2sz0 = fitInput->TreeSize( File2J[0] );
  int t2sz1 = fitInput->TreeSize( File2J[1] );
  int t2sz2 = fitInput->TreeSize( File2J[2] );
  int t2sz3 = fitInput->TreeSize( File2J[3] );

  TH1D* g2J = new TH1D("gW2J", " W+2J  ", 125, 50, 300 );

  //vector< vector<double> > vW2j   = RunEnsembles1( 2, t2sz1, inputMean[4], nRun, randomSeed, File2J[1], gW2J );
  //vector< vector<double> > vQCD2j = RunEnsembles1( 2, t2sz2, inputMean[5], nRun, randomSeed, File2J[2], gQ2J );
  //vector< vector<double> > vZ2j   = RunEnsembles1( 2, t2sz3, inputMean[6], nRun, randomSeed, File2J[3], gZ2J );
  vector<double> nT2j   = RunEnsembles1( 2, t2sz0, inputMean2J[0], nRun, randomSeed, File2J[0], g2J );
  vector<double> nW2j   = RunEnsembles1( 2, t2sz1, inputMean2J[1], nRun, randomSeed, File2J[1], g2J );
  vector<double> nQCD2j = RunEnsembles1( 2, t2sz2, inputMean2J[2], nRun, randomSeed, File2J[2], g2J );
  vector<double> nZ2j   = RunEnsembles1( 2, t2sz3, inputMean2J[3], nRun, randomSeed, File2J[3], g2J );
  
  cout<<"  size of 2J tt:"<< t2sz0 <<"  WJ:"<< t2sz1 <<"  QCD:"<< t2sz2 <<endl;

  // get the Cut Efficiency
  vector<double> EffCut = CutEff();
  // get the Ratio42
  //vector<double> R42 = bgEst->RatioXY( File4J, File2J, 0, true, false) ;
  vector<double> R42 = bgEst->RatioXY( 4, 2, File0J, 0, false, true, false, true  ) ;

  TH1D* hXtt = new TH1D("hXtt", " tt cross-section @ pb ", 24, 0, 360 );
  TH1D* hsXn = new TH1D("hsXn", " -Error of tt x-section @ pb ", 20, 0, 100 );
  TH1D* hsXp = new TH1D("hsXp", " +Error of tt x-section @ pb ", 20, 0, 100 );

  TH1D* hW2J = new TH1D("hW2J", " W+2J after cuts ", 50,  0, 300 );
  TH1D* hQ2J = new TH1D("hQ2J", " QCD2J after cuts ", 50,  0, 300 );
  TH1D* hW4J = new TH1D("hW4J", " W+4J after cuts ",  15, 0, 30 );
  TH1D* hQ4J = new TH1D("hQ4J", " QCD4J after cuts ",  15, 0, 30 );

  TH1D* hT4J = new TH1D("hT4J", " Tt+4J after cuts ", 25, 0, 50 );
  TH1D* eT4J = new TH1D("eT4J", " Expected Tt+4J ", 25, 0, 50 );

  TH1D* hB4J = new TH1D("hB4J", " BG+4J after cuts ", 50, 0, 50 );
  TH1D* eB4J = new TH1D("eB4J", " Expected BG+4J ", 50, 0, 50 );

  double n4JBg = 0; 
  double n4J   = 0;
  double n2J   = 0;
  for ( int i=0; i< nRun ; i++) {
      n4JBg = nWj[i]   + nQCD[i] + nZj[i] + nTq[i] + nTw[i] ;
      n4J   = n4JTt[i] + n4JBg     ;
      n2J   = nW2j[i]  + nQCD2j[i] + nZ2j[i] + nT2j[i] ;

      cout<<" ********** Run "<< i <<" ******************"<<endl;
      cout<<"  Ntt = "<<n4JTt[i] <<" Nbg = "<< n4JBg <<"   2JBg = "<< n2J  << endl;
      cout<<"  Ensembles 4J tt:"<< n4JTt[i] <<"  WJ:"<< nWj[i] <<"  QCD:"<< nQCD[i] <<" ZJ:"<< nZj[i] <<endl;
      
      vector<double> xtt = CrossSection(  n4J, n2J, R42, EffCut ) ;
      hXtt->Fill( xtt[0] );
      hsXn->Fill( xtt[1] );
      hsXp->Fill( xtt[2] );

      eB4J->Fill( xtt[3] );
      eT4J->Fill( n4J - xtt[3] );
      hW2J->Fill(   nW2j[0] );
      hQ2J->Fill( nQCD2j[0] );
      hW4J->Fill( nWj[i] );
      hQ4J->Fill( nQCD[i] );
      hT4J->Fill( n4JTt[i] );
      hB4J->Fill( n4JBg  );
  }

  gStyle->SetOptStat("neiruom");
  gStyle->SetOptFit(111);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->cd();
  hXtt->Draw();

  TF1* func0 = new TF1("func0", MassFitFunction::fitGS , 0, 300, 3);
  func0->SetParLimits(1, 50, 200);
  func0->SetParLimits(2, 20, 80);

  hXtt->Fit( func0, "R", "", 50., 300. );
  c1->Update();
  c1->Print(  theFolder+"ttXsection.gif"  );

  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->Divide(1,2);
  c2->cd(1);
  hsXn->Draw();
  c2->Update();
  c2->cd(2);
  hsXp->Draw();
  c2->Update();
  c2->Print(  theFolder+"ttXsectionErr.gif"  );

  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->Divide(2,2);
  c3->cd(1);
  hW2J->Draw();
  c3->Update();
  c3->cd(2);
  hQ2J->Draw();
  c3->Update();
  c3->cd(3);
  hW4J->Draw();
  c3->Update();
  c3->cd(4);
  hQ4J->Draw();
  c3->Update();
  c3->Print(  theFolder+"BGPoisson.gif"  );

  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  c4->Divide(2,2);
  c4->cd(1);
  hT4J->Draw();
  c4->Update();
  c4->cd(2);
  hB4J->Draw();
  c4->Update();
  c4->cd(3);
  eT4J->Draw();
  c4->Update();
  c4->cd(4);
  eB4J->Draw();
  c4->Update();
  c4->Print(  theFolder+"4JPoisson.gif"  );

  TCanvas* c5 = new TCanvas("c5","", 800, 600);
  c5->SetGrid();
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->Divide(2,2);
  c5->cd(1);
  gTt4J->Draw();
  c5->Update();

  c5->cd(2);
  gBg4J->Draw();
  c5->Update();

  c5->cd(3);
  g2J->Draw();
  c5->Update();

  c5->Print(  theFolder+"NGenInfo.gif"  );

  cout<<" !!! Ensemble Test Done !!! "<<endl;
}


vector<double> XSection::RunEnsembles1( int nJets, int tsz, double pMean, int nRun, int RandomSeed, string fileName, TH1D* hGen ){

  cout<<"  <<< Getting Ensembles for "<< nRun <<" >>> "<<endl;
  objInfo->Reset(0, nJets) ;

  // set up the random number function
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( RandomSeed );

  int  nRun_ = 0 ;
  vector<double> nCountV ;
  vector<int> ensembles ;
  int nRepeat = 0 ;
  while ( nRun_ < nRun) {

       // shuffle the events
       vector<int> evtline = pseudoExp->EventShuffle( tsz, RandomSeed );
       //if ( nRepeat > 0 ) cout<<"   smearing samples "<<endl;
       if ( nRepeat >= 0 ) smearing = true  ;
       bool nextRun = true ;
       int  RunStop = 0 ;
       int  NEvts   = 0 ;
       int  evtSize = static_cast<int>( evtline.size() );
       for (int k=0; k < evtSize; k++) {
 
           if ( nRun_ == nRun ) break; 
           if ( nextRun ) {
              //int PoiSeed = tRan->Integer( 1000+k );
	      //tRan->SetSeed( PoiSeed );
	      NEvts = tRan->Poisson( pMean );
	      nextRun = false ;
	      RunStop = k + NEvts ;
	      ensembles.clear() ;
              if ( RunStop >= evtSize ) break ;
	      cout<<"  Run["<<nRun_<<"] " ; 
              hGen->Fill( NEvts );
           }
           if ( k < RunStop ) {
              ensembles.push_back( evtline[k] );
           }
           if ( k == (RunStop - 1) || RunStop == 0 ) {

              bgCounter* counter = new bgCounter();
              objInfo->EvtSelector( fileName, counter, smearing, 1, &ensembles );
              vector<double> nCount;
              counter->CounterVec( nCount );
	      nCountV.push_back( nCount[0] );

	      nextRun = true ;
	      nRun_++ ;
	      cout<<" ("<<k<<") N of Evts needed For this Run = "<< NEvts <<" passed # = "<<nCount[0]  <<endl;
              //cout<<"-> "<< nCount[1]<<","<<nCount[2]<<","<<nCount[3]<<","<<nCount[4]<<","<<nCount[5] <<endl ;
              delete counter;
           }
       }
       nRepeat++ ;
  }
  cout<<" Number of Repeat = "<< nRepeat <<endl ;

  delete tRan ;
  return nCountV ;
}

void XSection::RealDataAnalysis( int nX, int nY ){

  string useDefault_ ;
  fitInput->GetParameters( "UseDefault", &useDefault_ ); 
  bool useDefault = ( useDefault_ == "True" ) ? true : false ;

  vector<double> RXY ;
  fitInput->GetParameters( "Ratio42", &RXY );

  vector<double> EffCut;
  if ( useDefault ) {
     fitInput->GetParameters( "EffCut",  &EffCut );
  } else {
     //RXY = bgEst->RatioXY( nX, nY, File0J, 0, false, inclX, inclY, true  ) ;
     EffCut = CutEff( nX );
  }

  // Data YJ background, Y = 2, always be exclusive
  objInfo->Reset(0, nY ) ;
  objInfo->Reset(0, false ) ;
  bgCounter* countyj_data = new bgCounter("datayj") ;
  objInfo->EvtSelector( dataFile[0], countyj_data, false, 1. );
  vector<double> nYJ_data ;
  countyj_data->CounterVec( nYJ_data );

  // Data XJ background , X = 4
  ACounter* countxj_data = new ACounter();
  wmfitter->ReFitSolution( dataFile[0], countxj_data, nX, 1, NULL );
  vector<double> nXJ_data = countxj_data->Output();

  vector<double> expBg4J  = bgEst->BgEstimate( RXY, nYJ_data ) ;
  vector<double> expBg4J0 = bgEst->BgEstimate( RXY, nYJ_data[0] ) ;

  cout<<" ====== Cross-section Calculation ====== "<< endl;

  cout<<"  Measured "<<nX<<"J Background0 = "<< expBg4J0[0] <<" + "<< expBg4J0[1] <<" - "<< expBg4J0[2] << endl; 

  cout<<"  Measured "<<nX<<"J BG in different Pt = "<< expBg4J[0] <<" + "<< expBg4J[1] <<" - "<< expBg4J[2] ;
  cout<<" | ( "<< expBg4J[3] <<", "<< expBg4J[4] <<", "<< expBg4J[5] <<", "<< expBg4J[6] <<", "<< expBg4J[7] <<" )"<< endl;
  cout<<" ----------------------------------------"<<endl;
  cout<<"  Data "<<nX<<"J  = "<< nXJ_data[0] <<"  +/- "<< sqrt( nXJ_data[1] )<<endl;
  cout<<"  Data "<<nY<<"J  = "<< nYJ_data[0] ;
  cout<<" ( "<< nYJ_data[1] <<" , "<< nYJ_data[2] <<" , "<< nYJ_data[3] <<" , "<< nYJ_data[4] <<" ) "<<endl;
  cout<<" ---------------------------------------- "<< endl;

  CrossSection( nXJ_data[0], nYJ_data[0], RXY, EffCut );
  cout<<" ======================================= "<< endl;

  delete countxj_data ;
  delete countyj_data ;

}

void XSection::DataForSystematic( string cfgFile, int nX  ){


  TString theFolder = hfolder ;
  gSystem->cd( theFolder );
  gSystem->mkdir( "Systematic"  );
  gSystem->cd("../");

  double lumi ;
  fitInput->GetParameters( "Lumi" , &lumi );
  vector<double> EffHLT ;
  fitInput->GetParameters( "EffHLT" , &EffHLT );
  vector<int> testPara;
  fitInput->GetParameters( "TestPara", &testPara, cfgFile ); 

  // scan muon pt and different cuts
  int sysType         = testPara[0] ;
  const static int sz = testPara[1] ;
  int centralValue    = testPara[2] ;
  int cut_topo        = testPara[3] ;

  string sfx = "" ;
  if ( sysType == 0 ) sfx = "_" ;
  if ( sysType == 1 ) sfx = "_i" ;
  if ( sysType == 2 ) sfx = "_j" ;
  if ( sysType == 3 ) sfx = "_q" ;
  TString gNames[4] = { "muon Pt", "muon Isolation", "JES %", " q^{2} scaling of W " } ;
  double  varScal[4] = { 1, 0.001, 1, 1 };

  double xV[sz], xErrn[sz], xErrp[sz], vars[sz], varsErr[sz] ;
  double nAll[sz], nBg[sz], Effs[sz] ;
  for (int k=0; k<3; k++ ) {
 
      int cut_idx = ( cut_topo == 9 ) ? k : cut_topo ;
      if ( k != cut_idx ) continue ;
      ostringstream cutIdStr ;
      cutIdStr << "_" ;
      cutIdStr << cut_idx ;

      for (int i=0; i<sz; i++ ) {
 
          double muonPt  = ( sysType == 0 ) ? (16.0+(i*2))   : 20. ;
	  double muonIso = ( sysType == 1 ) ? (0.05+(i*0.025)) : 0.1 ;
	  int var = 0 ;
	  if ( sysType == 0 ) var = 16 + (i*2) ;
	  if ( sysType == 1 ) var = 50 + (i*25) ;
	  if ( sysType == 2 ) var = 95 + (i*5) ;
	  if ( sysType == 3 ) var =  4 + i ;

	  ostringstream ratioStr ;
	  ratioStr << "Ratio42" ;
	  ratioStr << sfx ;
	  ratioStr << var ;

	  ostringstream effStr ;
	  effStr << "Eff" ;
	  effStr << cut_idx ;
	  effStr << sfx ;
	  effStr << var ;

	  ostringstream cutStr ;
	  cutStr << "TopoCut" ;
	  cutStr << cut_idx ;

	  vector<double> R42 ;
	  fitInput->GetParameters( ratioStr.str(), &R42, cfgFile );
	  vector<double> topoCuts;
	  fitInput->GetParameters( cutStr.str(),  &topoCuts, cfgFile );
	  vector<double> EffCuts;
	  fitInput->GetParameters( effStr.str(),  &EffCuts, cfgFile );

	  cout<<" **** Systematic test "<< var <<" *** "<<"muon pt = "<< muonPt <<"  Iso = "<< muonIso <<endl;
	  cout<<" Ratio42   = ( "<< R42[0]<<" , "<< R42[5] <<" )"<< endl;
	  cout<<" Cut Eff   = ( "<< EffCuts[0]<<", "<< EffCuts[1]<<", "<< EffCuts[2]<<", "<< EffCuts[3] <<" )"<< endl;
	  cout<<" Cut Vaule = ( "<< topoCuts[0]<<", "<< topoCuts[1]<<", "<< topoCuts[2]<<", "<< topoCuts[3]<<", ";
	  cout<< topoCuts[4]<<", "<< topoCuts[5] <<" )"<< endl;

	  // run the execise - use the cuts in DataCard
	  wmfitter->ResetCuts( topoCuts[0], topoCuts[1], topoCuts[2], topoCuts[3], topoCuts[4], topoCuts[5] );
	  wmfitter->SetMuonCuts( muonPt, -1, muonIso );

	  ACounter* countxj_data = new ACounter();
	  wmfitter->ReFitSolution( dataFile[0], countxj_data, nX, 1, NULL );
	  vector<double> nXJ_data = countxj_data->Output();

	  // Data YJ background, Y = 2, always be exclusive
	  objInfo->Reset(0, 2 ) ;
	  objInfo->Reset(0, false ) ;
	  objInfo->Reset(2, muonPt);  // set muon pt
	  objInfo->Reset(4, muonIso);  // set muon Iso
	  bgCounter* countyj_data = new bgCounter("datayj") ;
	  objInfo->EvtSelector( dataFile[0], countyj_data, false, 1. );
	  vector<double> nYJ_data ;
	  countyj_data->CounterVec( nYJ_data );
	  vector<double> xtt = CrossSection( nXJ_data[0], nYJ_data[0], R42, EffCuts );
	  cout<<" ======================================= "<< endl;

	  xV[i]     = xtt[0] ;
	  xErrn[i]  = xtt[1] ;
	  xErrp[i]  = xtt[2] ;
	  nAll[i]   = nXJ_data[0] ;
	  nBg[i]    = xtt[3] ;
	  Effs[i]   = EffCuts[0] ;
	  vars[i]   = var * varScal[sysType] ; 
	  varsErr[i] = 0.0 ; 

          delete countxj_data ;
          delete countyj_data ;
      }
      double dSg[sz] ;
      for ( int i=0; i< sz; i++) {
          double d_nAll = nAll[i] - nAll[centralValue] ;
	  double d_nBg  = nBg[i]  - nBg[centralValue]   ;
	  dSg[i] = sqrt( (d_nAll*d_nAll) + (d_nBg*d_nBg) ) / ( lumi*EffHLT[0]*Effs[i]);
          if ( sysType > 1 ) dSg[i] = fabs( xV[i] - xV[ centralValue ] ) ;
	  double percentErr = ( xV[i] - xV[ centralValue ] ) / xV[ centralValue ]  ;
	  cout<<" dX "<< xV[i] - xV[ centralValue ] << " %Err "<< percentErr ;
          cout<<"   dXsec ="<< dSg[i] <<"  Err% = "<< dSg[i]/xV[centralValue]  << endl;
          if ( dSg[i]/xV[i] > 1. )  dSg[i] = xV[i] ;
      }

      gStyle->SetOptTitle(0);
      gStyle->SetOptStat("");
      gStyle->SetLabelSize( 0.05, "X");
      gStyle->SetLabelSize( 0.05, "Y");
      gStyle->SetPadLeftMargin(0.20);
      gStyle->SetPadBottomMargin(0.15);
      //gStyle->SetErrorX(0) ;

      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      //c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->cd();

      TGraphAsymmErrors* gXs = new TGraphAsymmErrors( sz, vars, xV,  varsErr, varsErr, xErrn, xErrp );
      gXs->GetYaxis()->SetTitle("#sigma  #pm #Delta#sigma_{stat} (pb^{-1})" );
      gXs->GetYaxis()->SetTitleOffset(2.0);
      gXs->GetXaxis()->SetTitle( gNames[sysType] );
      gXs->GetXaxis()->SetTitleOffset(1.7);
      gXs->SetMarkerColor(4);
      gXs->SetMarkerStyle(21);
      gXs->SetMarkerSize(2.0);
      gXs->Draw("AP");
      c1->Update();

      TString plotname1 = hfolder + "Systematic/Xsec" + sfx + "DataSys" + cutIdStr.str() + "." + plotType ;
      c1->Print( plotname1 );

      TCanvas* c2 = new TCanvas("c2","", 800, 600);
      //c2->SetGrid();
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      c2->cd();

      TGraphErrors* gdX = new TGraphErrors( sz, vars, xV, varsErr, dSg );
      gdX->GetYaxis()->SetTitle("#sigma #pm #Delta#sigma_{syst}  (pb^{-1}) ");
      gdX->GetYaxis()->SetTitleOffset(2.0);
      gdX->GetXaxis()->SetTitleOffset(1.7);
      gdX->GetXaxis()->SetTitle( gNames[sysType] );
      gdX->SetMarkerColor(4);
      gdX->SetMarkerStyle(21);
      gdX->SetMarkerSize(2.0);
      gdX->Draw("AP");
      c2->Update();

      TString plotname2 = hfolder + "Systematic/dXsec" + sfx + "DataSys"+ cutIdStr.str() + "." + plotType ;
      c2->Print( plotname2 );

      delete c1;
      delete c2;
      delete gXs;
      delete gdX;
  }
}



vector<double> XSection::ErrAxB( double A, double s_A, double B, double s_B ){

	vector<double> sA = MassFitFunction::StatErr( A ) ;
	double sAp = ( s_A != -1 ) ? s_A : sA[1]; 
    double sAn = ( s_A != -1 ) ? s_A : -1*sA[0]; 
    vector<double> sB = MassFitFunction::StatErr( B ) ;
    double sBp = ( s_B != -1 ) ? s_B : sB[1]; 
    double sBn = ( s_B != -1 ) ? s_B : -1*sB[0]; 

    double s_fp = sqrt( B*B*sAp*sAp + A*A*sBp*sBp ) ;
    double s_fn = sqrt( B*B*sAn*sAn + A*A*sBn*sBn ) ;
 
    vector<double> sf ;
    sf.push_back( s_fn );
    sf.push_back( s_fp );
    return sf ;

}

