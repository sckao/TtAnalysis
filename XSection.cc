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

  fitInput->GetParameters( "InputMean", &inputMean );

}

XSection::~XSection(){

  delete fitInput ;
  delete wmfitter ;
  delete pseudoExp ;
  delete bgEst ;
  delete objInfo ;

}

// Exercise for background estimation
vector<double> XSection::Ratio42( int statIdx ){

  cout<<" ---  Measuring Ratio(4J/2J) --- "<<endl;
  vector<double> nEvts ;
  fitInput->GetParameters( "nEvents" , &nEvts );

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleTt   = ( statIdx != 0 ) ? scaleTt0*static_cast<double>( statIdx )    : scaleTt0  ;
  double scaleQCD  = ( statIdx != 0 ) ? scaleQCD0*static_cast<double>( statIdx )   : scaleQCD0 ;
  double scaleWJ   = ( statIdx != 0 ) ? scaleWJ0*static_cast<double>( statIdx )   : scaleWJ0  ;
  cout<<" renew scaleTt = "<< scaleTt <<"  scaleWJ = "<< scaleWJ <<"  scaleQCD = "<< scaleQCD <<endl; 

  // reset cuts -- remove m2m3 and leptonic m2t constrains
  wmfitter->ResetCuts( 0, 300, 0, 350, 0, 999 );

  // ratio(4J/2J) 
  ACounter* count4j = new ACounter();
  wmfitter->ReFitSolution( File4J[1], count4j, 4, scaleWJ,  NULL, statIdx, smearing );
  vector<double> nW4JEvts = count4j->Output();
  wmfitter->ReFitSolution( File4J[2], count4j, 4, scaleQCD, NULL, statIdx, smearing );
  vector<double> n4JEvts = count4j->Output();
  double rW4J = nW4JEvts[0]/n4JEvts[0] ;
  double rQ4J = 1. - rW4J  ;
  cout<<" NoCut MC BG4J = "<< n4JEvts[0] <<"  W4J = "<< nW4JEvts[0] <<" QCD4J "<< n4JEvts[0] - nW4JEvts[0] <<endl;
  cout<<"       Ratio_W4J = "<< rW4J  <<" Ratio_QCD4J = "<< rQ4J <<endl;

  ACounter* count2j = new ACounter();
  wmfitter->ReFitSolution( File2J[1], count2j, 2, scaleWJ, NULL,  statIdx, smearing );
  wmfitter->ReFitSolution( File2J[2], count2j, 2, scaleQCD, NULL, statIdx, smearing );
  vector<double> n2JEvts = count2j->Output();
  cout<<" MC BG2J = "<< n2JEvts[0] << endl;

  /// understand the signal contamination
  ACounter* count2jtt = new ACounter();
  wmfitter->ReFitSolution( File2J[0], count2jtt, 2, scaleTt, NULL, statIdx, smearing );
  vector<double> n2JttEvts = count2jtt->Output();
  cout<<" MC Tt2J  = "<< n2JttEvts[0] << endl;
  
  // ratio(4/2) without cuts
  vector<double> ratio42;
  ratio42.push_back(  n4JEvts[0] / n2JEvts[0]  ) ;
  double sR = sqrt( n4JEvts[1] + ratio42[0]*ratio42[0]*n2JEvts[1]  ) / n2JEvts[0] ;
  ratio42.push_back( sR ) ;
  cout<<" [0,1] MC R4/2 = "<< ratio42[0] <<" +/- "<< ratio42[1]   <<endl;

  // reset the cuts -- use the default value from DataCard
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );

  // signal efficiency with cuts
  ACounter* count4jSg = new ACounter();
  wmfitter->ReFitSolution( File4J[0], count4jSg, 4, scaleTt, NULL, statIdx, smearing );
  vector<double> n4JSg= count4jSg->Output();
  double  Eff = n4JSg[0] / (scaleTt0* nEvts[0]) ;
  double sEff = sqrt(  n4JSg[1] +  (Eff*Eff*scaleTt0*nEvts[0]) ) / (scaleTt0*nEvts[0]) ;
  ratio42.push_back( Eff );
  ratio42.push_back( sEff );
  cout<<" [2,3] Signal MC 4J = "<< n4JSg[0] <<" Eff = "<< Eff <<" +/- "<< sEff <<  endl;

  // 4J background cut efficiency 
  ACounter* count4jc = new ACounter();
  wmfitter->ReFitSolution( File4J[1], count4jc, 4, scaleWJ,  NULL, statIdx, smearing );
  vector<double> nW4JcEvts = count4jc->Output();
  wmfitter->ReFitSolution( File4J[2], count4jc, 4, scaleQCD, NULL, statIdx, smearing );
  vector<double> n4JcEvts = count4jc->Output();
  double cutEffW4J =  nW4JcEvts[0]/nW4JEvts[0] ;
  double cutEffQ4J = (n4JcEvts[0] -nW4JcEvts[0])/ (n4JEvts[0] - nW4JEvts[0]) ;
 
  cout<<" Cuts BG 4J = "<< n4JcEvts[0] <<"   W4J= "<< nW4JcEvts[0] <<" QCD4J= "<< n4JcEvts[0] - nW4JcEvts[0] <<endl;
  cout<<"     Cut Eff W4J= "<< cutEffW4J  <<" QCD4J= "<< cutEffQ4J  <<endl;

  // background efficiency => Ratio (4J_cuts/4J_all)
  double cutEff = n4JcEvts[0]/n4JEvts[0] ;
  double sEffBg_p = MassFitFunction::ErrAovB( n4JcEvts[0], n4JEvts[0], -1, -1, true );
  double sEffBg_n = MassFitFunction::ErrAovB( n4JcEvts[0], n4JEvts[0], -1, -1, false );

  //double cutEff1 = (rW4J*cutEffW4J) + (rQ4J*cutEffQ4J) ;
  //cout<<" cutEff1 = "<< cutEff1 <<endl ;
  ratio42.push_back( cutEff );
  ratio42.push_back( sEffBg_n );
  ratio42.push_back( sEffBg_p );

  cout<<" [4,5,6] CutEff(4J_cuts/4J_all) ="<< n4JcEvts[0]/n4JEvts[0] ;
  cout<<" - "<< sEffBg_n <<" + "<< sEffBg_p <<endl;

  delete count4j;
  delete count2j;
  delete count4jSg;
  delete count4jc;

  cout<<" Ratio Measured !! "<<endl;
  
  return ratio42 ;  
}

void XSection::MethodTest(){

  // 0:Ratio42 , 1: sigma_Ratio42 , 2 : Signal Eff , 3: sigma of Signal Eff
  vector<double> R42 = Ratio42( 2 ) ; 

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD  = scaleQCD0* 2.  ;
  double scaleWJ   = scaleWJ0* 2. ;
  double scaleTt   = scaleTt0* 2.;
  cout<<" renew scaleWJ = "<< scaleWJ <<"  scaleQCD = "<< scaleQCD <<endl; 

  // run the execise - use the cuts in DataCard
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  ACounter* count4jdata = new ACounter();
  wmfitter->ReFitSolution( File4J[1], count4jdata, 4, scaleWJ, NULL, -2, smearing );
  vector<double> nW4JData= count4jdata->Output();
  wmfitter->ReFitSolution( File4J[2], count4jdata, 4, scaleQCD, NULL, -2, smearing );
  vector<double> n4JData= count4jdata->Output();

  ACounter* count4jSgData = new ACounter();
  wmfitter->ReFitSolution( File4J[0], count4jSgData, 4, scaleTt, NULL, -2, smearing );
  vector<double> n4JSgData= count4jSgData->Output();

  wmfitter->ResetCuts( 0, 300, 0, 350, 0, 999 );
  ACounter* count2jdata = new ACounter();
  wmfitter->ReFitSolution( File2J[1], count2jdata, 2, scaleWJ,  NULL, -2, smearing );
  wmfitter->ReFitSolution( File2J[2], count2jdata, 2, scaleQCD, NULL, -2, smearing );
  vector<double> n2JData = count2jdata->Output();
 
  cout<<" ========== M2M3 =============== "<<endl;
  cout<<" Data Background 4J M2M3 = "<< n4JData[0] <<" WJ= "<< nW4JData[0] <<" QCD= "<< n4JData[0] - nW4JData[0]<< endl;
  cout<<" Signal 4J M2M3 = "<< n4JSgData[0] << endl;
  cout<<" Data Background 2J = "<< n2JData[0] << endl;
  GetXSection( R42, n2JData[0], n4JData[0]+n4JSgData[0] ) ;

  cout<<" =============================== "<<endl;

  delete count4jdata;
  delete count2jdata;
  delete count4jSgData;

  cout<<" Count !! "<<endl;
}

vector<double> XSection::GetXSection( vector<double>& R42, double n2J, double n4J ){

  vector<double> EffHLT ;
  fitInput->GetParameters( "EffHLT" , &EffHLT );
  double lumi ;
  fitInput->GetParameters( "Lumi" , &lumi );

  cout<<" ========== Cross-Section Calculation =============== "<<endl;
  cout<<" MC R4/2      = "<< R42[0] <<" +/- "<< R42[1]   <<endl;
  cout<<" Data 2J M2M3 = "<< n2J <<" +/- "<< sqrt( n2J ) <<endl;

  double  Exp4JBG_noCut = R42[0]*n2J ;
  vector<double> sExp4JBG_noCut = ErrAxB( R42[0], R42[1], n2J, -1 );
  cout<<" Expected 4J Background w/o cuts = "<< Exp4JBG_noCut <<" - "<<  sExp4JBG_noCut[0] <<" + "<< sExp4JBG_noCut[1] << endl;

  double  Exp4JBG = ( R42[4] < 1. ) ? Exp4JBG_noCut*R42[4] : Exp4JBG_noCut ; 
  vector<double> sExp4JBG ;

  if ( R42[4] < 1. ) {
     vector<double> sExp4JBG_n = ErrAxB( Exp4JBG_noCut, sExp4JBG_noCut[0], R42[4], R42[5] );
     vector<double> sExp4JBG_p = ErrAxB( Exp4JBG_noCut, sExp4JBG_noCut[1], R42[4], R42[6] );
     sExp4JBG.push_back( sExp4JBG_n[0] );
     sExp4JBG.push_back( sExp4JBG_p[0] );
  } else {
     sExp4JBG = sExp4JBG_noCut ;
  }
  cout<<" Expected 4J Background w/  cuts = "<< Exp4JBG <<" - "<<  sExp4JBG[0] <<" + "<< sExp4JBG[1] << endl;

  double xsecTt = ( n4J - Exp4JBG ) / ( R42[2]*lumi*EffHLT[0] );

  vector<double> s_n4J = MassFitFunction::StatErr( n4J );
  double s2_Ttn = (s_n4J[0]*s_n4J[0]) + (sExp4JBG[0]*sExp4JBG[0]);
  double s2_Ttp = (s_n4J[1]*s_n4J[1]) + (sExp4JBG[1]*sExp4JBG[1]);
  double nExp  = ( n4J - Exp4JBG ) / R42[2] ;
  double s_xsec_n = sqrt( s2_Ttn + nExp*nExp*R42[3]*R42[3] ) / ( R42[2]*lumi*EffHLT[0] );
  double s_xsec_p = sqrt( s2_Ttp + nExp*nExp*R42[3]*R42[3] ) / ( R42[2]*lumi*EffHLT[0] );

  cout<<" Expected 4J Signal M2M3     = "<< n4J - Exp4JBG <<" - "<< sqrt( s2_Ttn )<<" + "<< sqrt( s2_Ttp ) << endl;

  cout<<" *** Final XSection of Tt = "<< xsecTt <<" - "<< s_xsec_n <<" + "<< s_xsec_p <<endl ;
  cout<<" =============================== "<<endl;

  vector<double> xVs ;
  xVs.push_back( xsecTt );
  xVs.push_back( s_xsec_n );
  xVs.push_back( s_xsec_p );
  xVs.push_back( Exp4JBG );
  return xVs ;

}

void XSection::EnsembleTest( int nRun, int randomSeed ){

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  // 4J Samples
  //vector<string> File4J ;
  //fitInput->GetParameters( "4JSamples", &File4J );

  int tsz0 = fitInput->TreeSize( File4J[0] );
  int tsz1 = fitInput->TreeSize( File4J[1] );
  int tsz2 = fitInput->TreeSize( File4J[2] );

  TH1D* gT4J = new TH1D("gT4J", " Tt+4J ", 25, 0, 50 );
  TH1D* gW4J = new TH1D("gW4J", " W+4J  ",  15, 0, 30 );
  TH1D* gQ4J = new TH1D("gQ4J", " QCD4J ",  15, 0, 30 );

  vector<double> n4JTt = RunEnsembles( 4, tsz0, inputMean[0],  nRun, randomSeed, File4J[0], gT4J );
  vector<double> nWj   = RunEnsembles( 4, tsz1, inputMean[1],  nRun, randomSeed, File4J[1], gW4J );
  vector<double> nQCD  = RunEnsembles( 4, tsz2, inputMean[2], nRun, randomSeed, File4J[2], gQ4J );
  cout<<"  size of 4J tt:"<< tsz0 <<"  WJ:"<< tsz1 <<"  QCD:"<< tsz2 <<endl;

  // 2J Samples
  //vector<string> File2J ;
  //fitInput->GetParameters( "2JSamples", &File2J );

  int t2sz0 = fitInput->TreeSize( File2J[0] );
  int t2sz1 = fitInput->TreeSize( File2J[1] );
  int t2sz2 = fitInput->TreeSize( File2J[2] );

  TH1D* gW2J = new TH1D("gW2J", " W+2J  ", 80, 50, 210 );
  TH1D* gQ2J = new TH1D("gQ2J", " QCD2J ", 80, 50, 210 );

  wmfitter->ResetCuts( 0, 300, 0, 350, 0, 999 );
  
  vector<double> nW2j   = RunEnsembles( 2, t2sz1, inputMean[3], nRun, randomSeed, File2J[1], gW2J );
  vector<double> nQCD2j = RunEnsembles( 2, t2sz2, inputMean[4], nRun, randomSeed, File2J[2], gQ2J );
  
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  cout<<"  size of 2J tt:"<< t2sz0 <<"  WJ:"<< t2sz1 <<"  QCD:"<< t2sz2 <<endl;

  // get the ration42 and measure the cross-secton
  vector<double> R42 = Ratio42( 0 ) ;
  double n4JBg = 0; 
  double n2J   = 0;
  double n4J   = 0;

  TH1D* hXtt = new TH1D("hXtt", " tt cross-section @ 5/pb ", 24, 0, 360 );
  TH1D* hsXn = new TH1D("hsXn", " -Error of tt x-section @ 5/pb ", 20, 0, 100 );
  TH1D* hsXp = new TH1D("hsXp", " +Error of tt x-section @ 5/pb ", 20, 0, 100 );

  TH1D* hW2J = new TH1D("hW2J", " W+2J after cuts ", 40,  0, 200 );
  TH1D* hQ2J = new TH1D("hQ2J", " QCD2J after cuts ", 40,  0, 200 );
  TH1D* hW4J = new TH1D("hW4J", " W+4J after cuts ",  15, 0, 30 );
  TH1D* hQ4J = new TH1D("hQ4J", " QCD4J after cuts ",  15, 0, 30 );

  TH1D* hT4J = new TH1D("hT4J", " Tt+4J after cuts ", 25, 0, 50 );
  TH1D* eT4J = new TH1D("eT4J", " Expected Tt+4J ", 25, 0, 50 );

  TH1D* hB4J = new TH1D("hB4J", " BG+4J after cuts ", 50, 0, 50 );
  TH1D* eB4J = new TH1D("eB4J", " Expected BG+4J ", 50, 0, 50 );


  for ( size_t i=0; i< nRun ; i++) {
      n4JBg = nWj[i]   + nQCD[i]   ;
      n2J   = nW2j[i]  + nQCD2j[i] ;
      //n2J   = n2JEvts[0]  ;
      n4J   = n4JTt[i] + n4JBg     ;
      cout<<" ********** Run "<< i <<" ******************"<<endl;
      cout<<"  Ntt = "<<n4JTt[i] <<" Nbg = "<< n4JBg <<"   2JBg = "<< n2J << endl;
      cout<<"   ensembles 4J tt:"<< n4JTt[i] <<"  WJ:"<< nWj[i] <<"  QCD:"<< nQCD[i] <<endl;
      cout<<"   ensembles 2J wj:"<< nW2j[i] <<"  QCD:"<< nQCD2j[i] <<endl;
      //cout<<"   ensembles 2J wj:"<< nW2J[0] <<"  QCD:"<< nQCD2J <<endl;
      vector<double> xtt = GetXSection( R42, n2J, n4J ) ;
      hXtt->Fill( xtt[0] );
      hsXn->Fill( xtt[1] );
      hsXp->Fill( xtt[2] );

      eB4J->Fill( xtt[3] );
      eT4J->Fill( n4J - xtt[3] );
      hW2J->Fill( nW2j[i] );
      hQ2J->Fill( nQCD2j[i] );
      //hW2J->Fill( nW2J[0] );
      //hQ2J->Fill( nQCD2J );
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
  gT4J->Draw();
  c5->Update();

  c5->cd(2);
  gW4J->SetLineColor(2);
  gW4J->Draw();
  c5->Update();
  gStyle->SetStatY(0.60);
  gQ4J->SetLineColor(4);
  gQ4J->DrawCopy("sames");
  c5->Update();

  c5->cd(3);
  gStyle->SetStatY(0.95);
  gW2J->Draw();
  c5->Update();
  
  c5->cd(4);
  gQ2J->Draw();
  c5->Update();

  c5->Print(  theFolder+"NGenInfo.gif"  );

  cout<<" !!! Ensemble Test Done !!! "<<endl;
}


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
       for (int k=0; k < evtline.size(); k++) {
 
           if ( nRun_ == nRun ) break; 
           if ( nextRun ) {
              //int PoiSeed = tRan->Integer( 1000+k );
	      //tRan->SetSeed( PoiSeed );
	      NEvts = tRan->Poisson( pMean );
	      nextRun = false ;
	      RunStop = k + NEvts ;
	      ensembles.clear() ;
              if ( RunStop >= evtline.size() ) break ;
	      cout<<"  Run["<<nRun_<<"] " ; 
              hGen->Fill( NEvts );
           }
           if ( k < RunStop ) {
              ensembles.push_back( evtline[k] );
           }
           if ( k == (RunStop - 1) ) {
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

vector<double> XSection::CutEff(){

  vector<double> nEvts ;
  fitInput->GetParameters( "nEvents" , &nEvts );

  double scaleTt  = fitInput->NormalizeComponents( "tt" );
  double scaleWJ  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD = fitInput->NormalizeComponents( "qcd" );
  double scaleZJ  = fitInput->NormalizeComponents( "zj" );

  // reset the cuts -- use the default value from DataCard
  // Without cuts
  wmfitter->ResetCuts( 0, 999, 0, 999, 0, 999 );

  ACounter* count4jSg = new ACounter();
  wmfitter->ReFitSolution( File4J[0], count4jSg, 4, scaleTt, NULL, 0, smearing );
  vector<double> n4JSg= count4jSg->Output();

  /*
  ACounter* count4jBg = new ACounter();
  wmfitter->ReFitSolution( File4J[1], count4jBg, 4, scaleWJ,  NULL, 0, smearing );
  wmfitter->ReFitSolution( File4J[2], count4jBg, 4, scaleQCD, NULL, 0, smearing );
  wmfitter->ReFitSolution( File4J[3], count4jBg, 4, scaleZJ, NULL, 0, smearing );
  vector<double> n4JBg = count4jBg->Output();
  */

  objInfo->Reset(0, 4) ;
  bgCounter* count4jBg = new bgCounter("eff") ;
  objInfo->EvtSelector( File4J[1], count4jBg, smearing, scaleWJ );
  objInfo->EvtSelector( File4J[2], count4jBg, smearing, scaleQCD );
  objInfo->EvtSelector( File4J[3], count4jBg, smearing, scaleZJ );
  vector<double> n4JBg ;
  count4jBg->CounterVec( n4JBg );

  // With cuts
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );

  ACounter* count4jSg_c = new ACounter();
  wmfitter->ReFitSolution( File4J[0], count4jSg_c, 4, scaleTt, NULL, 0, smearing );
  vector<double> n4JSg_cuts = count4jSg_c->Output();

  ACounter* count4jBg_c = new ACounter();
  wmfitter->ReFitSolution( File4J[1], count4jBg_c, 4, scaleWJ,  NULL, 0, smearing );
  wmfitter->ReFitSolution( File4J[2], count4jBg_c, 4, scaleQCD, NULL, 0, smearing );
  wmfitter->ReFitSolution( File4J[3], count4jBg_c, 4, scaleZJ, NULL, 0, smearing );
  vector<double> n4JBg_cuts = count4jBg_c->Output();

  // signal efficiency with cuts
  double  Eff = n4JSg_cuts[0] / (nEvts[0]*scaleTt) ;
  double sEff = MassFitFunction::ErrAovB( n4JSg_cuts[0], nEvts[0], n4JSg_cuts[1], -1 ) / scaleTt ;

  // 4J background cut efficiency 
  double Eff_Bg  =  n4JBg_cuts[0] / n4JBg[0] ;
  double sEff_Bg = MassFitFunction::ErrAovB( n4JBg_cuts[0], n4JBg[0], n4JBg_cuts[1], n4JBg[6] );

  cout<<" ===== Estimate Cut Efficiency ===== "<<endl;
  cout<<" MC Signal 4J = "<< n4JSg[0] <<" Eff = "<< Eff <<" +/- "<< sEff <<  endl;
  cout<<" MC Background 4J = "<< n4JBg[0] <<" Eff = "<< Eff_Bg  <<" +/- "<< sEff_Bg <<endl;
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

vector<double> XSection::CrossSection( double nData_4J, vector<double>& nData_2J, vector<double>& R42, vector<double>& EffCut ){

  double lumi ;
  fitInput->GetParameters( "Lumi" , &lumi );
  vector<double> EffHLT ;
  fitInput->GetParameters( "EffHLT" , &EffHLT );

  vector<double> nBg_4J = bgEst->BgEstimate( R42, nData_2J ) ;

  double xsecTt = ( nData_4J - (nBg_4J[0]*EffCut[1])*1.2 ) / ( lumi*EffHLT[0]*EffCut[0] );

  // calculate the uncertainty of cross-section 
  vector<double> s_n4J = MassFitFunction::StatErr( nData_4J );
  double sn_Bg = MassFitFunction::ErrAxB( nBg_4J[0], EffCut[1],  nBg_4J[2], EffCut[3] ) * 1.2;
  double sp_Bg = MassFitFunction::ErrAxB( nBg_4J[0], EffCut[1],  nBg_4J[1], EffCut[3] ) * 1.2;

  double sn_Tt = sqrt( (s_n4J[0]*s_n4J[0]) + ( sn_Bg*sn_Bg ) );
  double sp_Tt = sqrt( (s_n4J[1]*s_n4J[1]) + ( sp_Bg*sp_Bg ) );
  double n_Tt = nData_4J - (nBg_4J[0]*EffCut[1])*1.2 ;

  double sn_xsec0 = MassFitFunction::ErrAovB( n_Tt, EffCut[0], sn_Tt, EffCut[2] );
  double sp_xsec0 = MassFitFunction::ErrAovB( n_Tt, EffCut[0], sp_Tt, EffCut[2] );
  double sn_xsec = sn_xsec0/(lumi*EffHLT[0]) ;
  double sp_xsec = sp_xsec0/(lumi*EffHLT[0]) ;

  cout<<" ========== Cross-Section Report ============ "<<endl;
  cout<<" N of All 4J = "<< nData_4J <<endl;
  cout<<" N of All 2J = "<<  nData_2J[0] ;
  cout<<" -> "<< nData_2J[1]<<","<<nData_2J[2]<<","<<nData_2J[3]<<","<<nData_2J[4]<<","<<nData_2J[5]<<endl;
  cout<<" N of Tt = "<< n_Tt <<" - "<< sn_Tt<<" + "<< sp_Tt <<endl;
  cout<<" N of Bg "<< nBg_4J[0]*1.2 <<" Eff= "<<EffCut[1] <<" - "<<sn_Bg<<" + "<< sp_Bg << endl;
  cout<<" Tt Xsec = "<< xsecTt <<" - "<< sn_xsec <<" + "<< sp_xsec << endl;
  cout<<" ============================================ "<<endl;

  vector<double> xVs ;
  xVs.push_back( xsecTt );
  xVs.push_back( sn_xsec );
  xVs.push_back( sp_xsec );
  xVs.push_back( nBg_4J[0]*EffCut[1] );
  return xVs ;

}

void XSection::MethodTest1(){

  // get the Cut Efficiency
  vector<double> EffCut = CutEff();
  // get the Ratio42
  vector<double> R42 = bgEst->RatioXY( File4J, File2J, -1, false, false ) ;

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ0  = fitInput->NormalizeComponents( "zj" );

  int idx = -3 ;
  double scaleTt   = scaleTt0* 3.;
  double scaleWJ   = scaleWJ0* 3. ;
  double scaleQCD  = scaleQCD0* 3.  ;
  double scaleZJ   = scaleWJ0* 3. ;
  

  // run the execise - use the cuts in DataCard
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );

  ACounter* count4j_sg = new ACounter();
  wmfitter->ReFitSolution( File4J[0], count4j_sg, 4, scaleTt, NULL, idx, smearing );
  vector<double> n4J_sg= count4j_sg->Output();

  ACounter* count4j_bg = new ACounter();
  wmfitter->ReFitSolution( File4J[1], count4j_bg, 4, scaleWJ,  NULL, idx, smearing );
  wmfitter->ReFitSolution( File4J[2], count4j_bg, 4, scaleQCD, NULL, idx, smearing );
  wmfitter->ReFitSolution( File4J[3], count4j_bg, 4, scaleZJ,  NULL, idx, smearing );
  vector<double> n4J_bg = count4j_bg->Output();
  double totalN4J = n4J_sg[0] + n4J_bg[0] ;

  objInfo->Reset(0, 2) ;
  bgCounter* count2j_bg = new bgCounter("eff") ;
  objInfo->EvtSelector( File2J[1], count2j_bg, smearing, scaleWJ, NULL, idx );
  objInfo->EvtSelector( File2J[2], count2j_bg, smearing, scaleQCD, NULL, idx );
  objInfo->EvtSelector( File2J[3], count2j_bg, smearing, scaleZJ, NULL, idx );
  vector<double> n2J_bg ;
  count2j_bg->CounterVec( n2J_bg );

  cout<<" ===== Method Test ===== "<<endl;
  cout<<"  Generated 4J Tt = "<< n4J_sg[0] <<" Bg = "<<n4J_bg[0] <<endl;

  vector<double> xtt = CrossSection(  totalN4J, n2J_bg, R42, EffCut ) ;
 
  delete count4j_sg ;
  delete count4j_bg ;
  delete count2j_bg ;

}

void XSection::BgClosureTest( int nX, int nY, bool inclX, bool inclY ){

  // get the Ratio41
  vector<double> RXY = bgEst->RatioXY( nX, nY, File0J, -1, false, inclX, inclY  ) ;

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ0  = fitInput->NormalizeComponents( "zj" );

  int idx = -2 ;
  double scaleTt   = scaleTt0* 2.;
  double scaleWJ   = scaleWJ0* 2. ;
  double scaleQCD  = scaleQCD0* 2.  ;
  double scaleZJ   = scaleWJ0* 2. ;
  
  // actural XJ background
  objInfo->Reset(0, nX ) ;
  objInfo->Reset(0, inclX ) ;
  bgCounter* count4j_bg = new bgCounter("bgxj") ;
  //objInfo->EvtSelector( File0J[0], count4j_bg, smearing, scaleTt, NULL, idx );
  objInfo->EvtSelector( File0J[1], count4j_bg, smearing, scaleWJ, NULL, idx );
  objInfo->EvtSelector( File0J[2], count4j_bg, smearing, scaleQCD, NULL, idx );
  objInfo->EvtSelector( File0J[3], count4j_bg, smearing, scaleZJ, NULL, idx );
  vector<double> n4J_bg ;
  count4j_bg->CounterVec( n4J_bg );

  // Input YJ background
  objInfo->Reset(0, nY ) ;
  objInfo->Reset(0, inclY ) ;
  bgCounter* count1j_bg = new bgCounter("bgyj") ;
  objInfo->EvtSelector( File0J[0], count1j_bg, smearing, scaleTt, NULL, idx );
  objInfo->EvtSelector( File0J[1], count1j_bg, smearing, scaleWJ, NULL, idx );
  objInfo->EvtSelector( File0J[2], count1j_bg, smearing, scaleQCD, NULL, idx );
  objInfo->EvtSelector( File0J[3], count1j_bg, smearing, scaleZJ, NULL, idx );
  vector<double> n1J_bg ;
  count1j_bg->CounterVec( n1J_bg );

  vector<double> expBg4J  = bgEst->BgEstimate( RXY, n1J_bg ) ;
  vector<double> expBg4J0 = bgEst->BgEstimate( RXY, n1J_bg[0] ) ;

  cout<<" ===== Background Closure Test ===== "<< endl;
  cout<<"  Generated "<<nX<<" Background = "<< n4J_bg[0] ;
  cout<<" ( "<< n4J_bg[1] <<" , "<< n4J_bg[2] <<" , "<< n4J_bg[3] <<" , "<< n4J_bg[4] <<" )"<< endl;

  cout<<"  Measured "<<nX<<" Background0 = "<< expBg4J0[0] <<" + "<< expBg4J0[1] <<" - "<< expBg4J0[2] << endl; 

  cout<<"  Measured "<<nX<<" BG in different Pt = "<< expBg4J[0] <<" + "<< expBg4J[1] <<" - "<< expBg4J[2] << endl;
  cout<<" ( "<< expBg4J[3] <<" , "<< expBg4J[4] <<" , "<< expBg4J[5] <<" , "<< expBg4J[6] <<" )"<< endl;
  cout<<" =================================== "<< endl;

  delete count4j_bg ;
  delete count1j_bg ;

}

void XSection::EnsembleTest1( int nRun, int randomSeed ){

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  // 4J Samples
  int tsz0 = fitInput->TreeSize( File4J[0] );
  int tsz1 = fitInput->TreeSize( File4J[1] );
  int tsz2 = fitInput->TreeSize( File4J[2] );
  int tsz3 = fitInput->TreeSize( File4J[3] );

  TH1D* gT4J = new TH1D("gT4J", " Tt+4J ", 25, 0, 50 );
  TH1D* gW4J = new TH1D("gW4J", " W+4J  ",  15, 0, 30 );
  TH1D* gQ4J = new TH1D("gQ4J", " QCD4J ",  15, 0, 30 );
  TH1D* gZ4J = new TH1D("gZ4J", " Z+4J ",  15, 0, 30 );

  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  vector<double> n4JTt = RunEnsembles( 4, tsz0, inputMean[0],  nRun, randomSeed, File4J[0], gT4J );
  vector<double> nWj   = RunEnsembles( 4, tsz1, inputMean[1],  nRun, randomSeed, File4J[1], gW4J );
  vector<double> nQCD  = RunEnsembles( 4, tsz2, inputMean[2],  nRun, randomSeed, File4J[2], gQ4J );
  vector<double> nZj   = RunEnsembles( 4, tsz3, inputMean[3],  nRun, randomSeed, File4J[3], gZ4J );
  cout<<"  size of 4J tt:"<< tsz0 <<"  WJ:"<< tsz1 <<"  QCD:"<< tsz2 <<" ZJ: "<< tsz3<<endl;

  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  // 2J Samples
  int t2sz0 = fitInput->TreeSize( File2J[0] );
  int t2sz1 = fitInput->TreeSize( File2J[1] );
  int t2sz2 = fitInput->TreeSize( File2J[2] );
  int t2sz3 = fitInput->TreeSize( File2J[3] );

  TH1D* gW2J = new TH1D("gW2J", " W+2J  ", 80, 50, 210 );
  TH1D* gQ2J = new TH1D("gQ2J", " QCD2J ", 80, 50, 210 );
  TH1D* gZ2J = new TH1D("gZ2J", " Z+2J ", 80, 50, 210 );

  vector< vector<double> > vW2j   = RunEnsembles1( 2, t2sz1, inputMean[4], nRun, randomSeed, File2J[1], gW2J );
  vector< vector<double> > vQCD2j = RunEnsembles1( 2, t2sz2, inputMean[5], nRun, randomSeed, File2J[2], gQ2J );
  vector< vector<double> > vZ2j   = RunEnsembles1( 2, t2sz3, inputMean[6], nRun, randomSeed, File2J[3], gZ2J );
  
  cout<<"  size of 2J tt:"<< t2sz0 <<"  WJ:"<< t2sz1 <<"  QCD:"<< t2sz2 <<endl;

  double n4JBg = 0; 
  double n4J   = 0;
  vector<double> n2J ;

  // get the Cut Efficiency
  vector<double> EffCut = CutEff();
  // get the Ratio42
  vector<double> R42 = bgEst->RatioXY( File4J, File2J, -1, false, false ) ;

  TH1D* hXtt = new TH1D("hXtt", " tt cross-section @ 5/pb ", 24, 0, 360 );
  TH1D* hsXn = new TH1D("hsXn", " -Error of tt x-section @ 5/pb ", 20, 0, 100 );
  TH1D* hsXp = new TH1D("hsXp", " +Error of tt x-section @ 5/pb ", 20, 0, 100 );

  TH1D* hW2J = new TH1D("hW2J", " W+2J after cuts ", 40,  0, 200 );
  TH1D* hQ2J = new TH1D("hQ2J", " QCD2J after cuts ", 40,  0, 200 );
  TH1D* hW4J = new TH1D("hW4J", " W+4J after cuts ",  15, 0, 30 );
  TH1D* hQ4J = new TH1D("hQ4J", " QCD4J after cuts ",  15, 0, 30 );

  TH1D* hT4J = new TH1D("hT4J", " Tt+4J after cuts ", 25, 0, 50 );
  TH1D* eT4J = new TH1D("eT4J", " Expected Tt+4J ", 25, 0, 50 );

  TH1D* hB4J = new TH1D("hB4J", " BG+4J after cuts ", 50, 0, 50 );
  TH1D* eB4J = new TH1D("eB4J", " Expected BG+4J ", 50, 0, 50 );


  for ( size_t i=0; i< nRun ; i++) {
      n4JBg = nWj[i]   + nQCD[i] + nZj[i]  ;
      n4J   = n4JTt[i] + n4JBg     ;

      vector<double> nW2j   = vW2j[i] ;
      vector<double> nQCD2j = vQCD2j[i] ;
      vector<double> nZ2j   = vZ2j[i] ;
      n2J.clear() ;
      for ( size_t k=0; k < 6; k++) {
          n2J.push_back(  nW2j[k]  + nQCD2j[k] + nZ2j[k] );
      }

      cout<<" ********** Run "<< i <<" ******************"<<endl;
      cout<<"  Ntt = "<<n4JTt[i] <<" Nbg = "<< n4JBg <<"   2JBg = "<< n2J[0]  << endl;
      cout<<"  Ensembles 4J tt:"<< n4JTt[i] <<"  WJ:"<< nWj[i] <<"  QCD:"<< nQCD[i] <<" ZJ:"<< nZj[i] <<endl;
      
      vector<double> xtt = CrossSection(  n4J, n2J, R42, EffCut ) ;
      hXtt->Fill( xtt[0] );
      hsXn->Fill( xtt[1] );
      hsXp->Fill( xtt[2] );

      eB4J->Fill( xtt[3] );
      eT4J->Fill( n4J - xtt[3] );
      hW2J->Fill(   nW2j[0] );
      hQ2J->Fill( nQCD2j[0] );
      //hW2J->Fill( nW2J[0] );
      //hQ2J->Fill( nQCD2J );
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
  gT4J->Draw();
  c5->Update();

  c5->cd(2);
  gW4J->SetLineColor(2);
  gW4J->Draw();
  c5->Update();
  gStyle->SetStatY(0.60);
  gQ4J->SetLineColor(4);
  gQ4J->DrawCopy("sames");
  c5->Update();

  c5->cd(3);
  gStyle->SetStatY(0.95);
  gW2J->Draw();
  c5->Update();
  
  c5->cd(4);
  gQ2J->Draw();
  c5->Update();

  c5->Print(  theFolder+"NGenInfo.gif"  );

  cout<<" !!! Ensemble Test Done !!! "<<endl;
}


vector< vector<double> > XSection::RunEnsembles1( int nJets, int tsz, double pMean, int nRun, int RandomSeed, string fileName, TH1D* hGen ){

  cout<<"  <<< Getting Ensembles for "<< nRun <<" >>> "<<endl;
  objInfo->Reset(0, nJets) ;

  // set up the random number function
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( RandomSeed );

  int  nRun_ = 0 ;
  vector< vector<double> > nCountV ;
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
       for (int k=0; k < evtline.size(); k++) {
 
           if ( nRun_ == nRun ) break; 
           if ( nextRun ) {
              //int PoiSeed = tRan->Integer( 1000+k );
	      //tRan->SetSeed( PoiSeed );
	      NEvts = tRan->Poisson( pMean );
	      nextRun = false ;
	      RunStop = k + NEvts ;
	      ensembles.clear() ;
              if ( RunStop >= evtline.size() ) break ;
	      cout<<"  Run["<<nRun_<<"] " ; 
              hGen->Fill( NEvts );
           }
           if ( k < RunStop ) {
              ensembles.push_back( evtline[k] );
           }
           if ( k == (RunStop - 1) ) {

              bgCounter* counter = new bgCounter();
              objInfo->EvtSelector( fileName, counter, smearing, 1, &ensembles );
              vector<double> nCount;
              counter->CounterVec( nCount );
	      nCountV.push_back( nCount );

	      nextRun = true ;
	      nRun_++ ;
	      cout<<" ("<<k<<") N of Evts needed For this Run = "<< NEvts <<" passed # = "<<nCount[0]  ;
              cout<<"-> "<< nCount[1]<<","<<nCount[2]<<","<<nCount[3]<<","<<nCount[4]<<","<<nCount[5] <<endl ;
              delete counter;
           }
       }
       nRepeat++ ;
  }
  cout<<" Number of Repeat = "<< nRepeat <<endl ;

  delete tRan ;
  return nCountV ;
}

void XSection::RealBackground( int nX, int nY, bool inclX, bool inclY ){

  
  // get the RatioXY
  vector<double> RXY = bgEst->RatioXY( nX, nY, File0J, -1, false, inclX, inclY  ) ;

  double scaleTt  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ  = fitInput->NormalizeComponents( "wj" );
  double scaleZJ  = fitInput->NormalizeComponents( "zj" );

  // M.C. XJ background
  objInfo->Reset(0, nX ) ;
  objInfo->Reset(0, inclX ) ;
  bgCounter* countxj_bg = new bgCounter("mcxj") ;
  //objInfo->EvtSelector( File0J[0], countxj_bg, smearing, scaleTt, NULL );
  objInfo->EvtSelector( File0J[1], countxj_bg, smearing, scaleWJ, NULL );
  objInfo->EvtSelector( File0J[2], countxj_bg, smearing, scaleQCD, NULL );
  objInfo->EvtSelector( File0J[3], countxj_bg, smearing, scaleZJ, NULL );
  vector<double> mcXJ_bg ;
  countxj_bg->CounterVec( mcXJ_bg );

  // Data YJ background
  objInfo->Reset(0, nY ) ;
  objInfo->Reset(0, inclY ) ;
  bgCounter* countyj_data = new bgCounter("datayj") ;
  objInfo->EvtSelector( dataFile[0], countyj_data, false, 1. );
  vector<double> nYJ_data ;
  countyj_data->CounterVec( nYJ_data );

  // Data XJ background
  objInfo->Reset(0, nX ) ;
  objInfo->Reset(0, inclX ) ;
  bgCounter* countxj_data = new bgCounter("dataxj") ;
  objInfo->EvtSelector( dataFile[0], countxj_data, false, 1. );
  vector<double> nXJ_data ;
  countxj_data->CounterVec( nXJ_data );

  vector<double> expBg4J  = bgEst->BgEstimate( RXY, nYJ_data ) ;
  vector<double> expBg4J0 = bgEst->BgEstimate( RXY, nYJ_data[0] ) ;

  cout<<" ===== Real Background Calculation ===== "<< endl;
  cout<<"  MC "<<nX<<"J Background = "<< mcXJ_bg[0] ;
  cout<<" ( "<< mcXJ_bg[1] <<" , "<< mcXJ_bg[2] <<" , "<< mcXJ_bg[3] <<" , "<< mcXJ_bg[4] <<" , "<< mcXJ_bg[5] <<" )"<< endl;

  cout<<"  Measured "<<nX<<"J Background0 = "<< expBg4J0[0] <<" + "<< expBg4J0[1] <<" - "<< expBg4J0[2] << endl; 

  cout<<"  Measured "<<nX<<"J BG in different Pt = "<< expBg4J[0] <<" + "<< expBg4J[1] <<" - "<< expBg4J[2] << endl;
  cout<<" ( "<< expBg4J[3] <<" , "<< expBg4J[4] <<" , "<< expBg4J[5] <<" , "<< expBg4J[6] <<" , "<< expBg4J[7] <<" )"<< endl;
  cout<<" ----------------------------------------"<<endl;
  cout<<"  Data "<<nX<<"J in different Pt = "<< nXJ_data[0] << endl;
  cout<<" ( "<< nXJ_data[1] <<" , "<< nXJ_data[2] <<" , "<< nXJ_data[3] <<" , "<< nXJ_data[4] <<" , "<< nXJ_data[5] <<" )"<< endl;
  cout<<" ======================================= "<< endl;

  cout<<" ====== Cross-section Calculation ====== "<< endl;
  vector<double> EffCut = CutEff();
  CrossSection( nXJ_data[0], nYJ_data, RXY, EffCut );
  cout<<" ======================================= "<< endl;

  delete countxj_bg ;
  delete countxj_data ;
  delete countyj_data ;

}

vector<double> XSection::ErrAxB( double A, double s_A, double B, double s_B ){

    vector<double> sA = MassFitFunction::StatErr( A ) ;
    double sAp = ( s_A != -1 ) ? s_A : sA[1]; 
    double sAn = ( s_A != -1 ) ? s_A : -1*sA[0]; 
    vector<double> sB = MassFitFunction::StatErr( B ) ;
    double sBp = ( s_B != -1 ) ? s_B : sB[1]; 
    double sBn = ( s_B != -1 ) ? s_B : -1*sB[0]; 

    double f = A * B ;
    double s_fp = sqrt( B*B*sAp*sAp + A*A*sBp*sBp ) ;
    double s_fn = sqrt( B*B*sAn*sAn + A*A*sBn*sBn ) ;
 
    vector<double> sf ;
    sf.push_back( s_fn );
    sf.push_back( s_fp );
    return sf ;

}

