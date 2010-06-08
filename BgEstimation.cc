#include "BgEstimation.h"
#include "WFormat.h"

BgEstimation::BgEstimation(){

  fitFunc   = new MassFitFunction();
  fitInput  = new MassAnaInput();
  wmfitter  = new HadWMassFitter();
  pseudoExp = new PseudoExp();

  string phaseSmear ;
  fitInput->GetParameters( "PhaseSmear", &phaseSmear );
  smearing = ( phaseSmear == "ON" ) ? true : false ;

  fitInput->GetParameters( "Path", &hfolder );
  theFolder = hfolder ;
  gSystem->mkdir( theFolder );

  fitInput->GetParameters( "InputMean", &inputMean );

}

BgEstimation::~BgEstimation(){

  delete fitFunc;
  delete fitInput ;
  delete wmfitter ;
  delete pseudoExp ;

}

// Exercise for background estimation
vector<double> BgEstimation::Ratio42( int statIdx ){

  cout<<" ---  Measuring Ratio(4J/2J) --- "<<endl;
  vector<string> File4J ;
  fitInput->GetParameters( "FakeData", &File4J );

  vector<string> File2J ;
  fitInput->GetParameters( "2JSamples", &File2J );

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
  bgCounter* count4j = new bgCounter();
  wmfitter->ReFitSolution( File4J[1], count4j, scaleWJ,  NULL, statIdx, smearing );
  vector<double> nW4JEvts = count4j->Output();
  wmfitter->ReFitSolution( File4J[2], count4j, scaleQCD, NULL, statIdx, smearing );
  vector<double> n4JEvts = count4j->Output();
  double rW4J = nW4JEvts[0]/n4JEvts[0] ;
  double rQ4J = 1. - rW4J  ;
  cout<<" NoCut MC BG4J = "<< n4JEvts[0] <<"  W4J = "<< nW4JEvts[0] <<" QCD4J "<< n4JEvts[0] - nW4JEvts[0] <<endl;
  cout<<"       Ratio_W4J = "<< rW4J  <<" Ratio_QCD4J = "<< rQ4J <<endl;

  bgCounter* count2j = new bgCounter();
  wmfitter->ReFitSolution( File2J[1], count2j, scaleWJ, NULL,  statIdx, smearing );
  wmfitter->ReFitSolution( File2J[2], count2j, scaleQCD, NULL, statIdx, smearing );
  vector<double> n2JEvts = count2j->Output();
  cout<<" MC BG2J = "<< n2JEvts[0] << endl;

  /// understand the signal contamination
  bgCounter* count2jtt = new bgCounter();
  wmfitter->ReFitSolution( File2J[0], count2jtt, scaleTt, NULL, statIdx, smearing );
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
  bgCounter* count4jSg = new bgCounter();
  wmfitter->ReFitSolution( File4J[0], count4jSg, scaleTt, NULL, statIdx, smearing );
  vector<double> n4JSg= count4jSg->Output();
  double  Eff = n4JSg[0] / (scaleTt0* nEvts[0]) ;
  double sEff = sqrt(  n4JSg[1] +  (Eff*Eff*scaleTt0*nEvts[0]) ) / (scaleTt0*nEvts[0]) ;
  ratio42.push_back( Eff );
  ratio42.push_back( sEff );
  cout<<" [2,3] Signal MC 4J = "<< n4JSg[0] <<" Eff = "<< Eff <<" +/- "<< sEff <<  endl;

  // 4J background cut efficiency 
  bgCounter* count4jc = new bgCounter();
  wmfitter->ReFitSolution( File4J[1], count4jc, scaleWJ,  NULL, statIdx, smearing );
  vector<double> nW4JcEvts = count4jc->Output();
  wmfitter->ReFitSolution( File4J[2], count4jc, scaleQCD, NULL, statIdx, smearing );
  vector<double> n4JcEvts = count4jc->Output();
  double cutEffW4J =  nW4JcEvts[0]/nW4JEvts[0] ;
  double cutEffQ4J = (n4JcEvts[0] -nW4JcEvts[0])/ (n4JEvts[0] - nW4JEvts[0]) ;
 
  cout<<" Cuts BG 4J = "<< n4JcEvts[0] <<"   W4J= "<< nW4JcEvts[0] <<" QCD4J= "<< n4JcEvts[0] - nW4JcEvts[0] <<endl;
  cout<<"     Cut Eff W4J= "<< cutEffW4J  <<" QCD4J= "<< cutEffQ4J  <<endl;

  // background efficiency => Ratio (4J_cuts/4J_all)
  double cutEff = n4JcEvts[0]/n4JEvts[0] ;
  vector<double> sEffBg = ErrAovB( n4JcEvts[0], -1, n4JEvts[0], -1 );

  //double cutEff1 = (rW4J*cutEffW4J) + (rQ4J*cutEffQ4J) ;
  //cout<<" cutEff1 = "<< cutEff1 <<endl ;
  ratio42.push_back( cutEff );
  ratio42.push_back( sEffBg[0] );
  ratio42.push_back( sEffBg[1] );

  cout<<" [4,5,6] CutEff(4J_cuts/4J_all) ="<< n4JcEvts[0]/n4JEvts[0] ;
  cout<<" - "<< sEffBg[0] <<" + "<< sEffBg[1] <<endl;

  delete count4j;
  delete count2j;
  delete count4jSg;
  delete count4jc;

  cout<<" Ratio Measured !! "<<endl;
  
  return ratio42 ;  
}

void BgEstimation::MethodTest(){

  // 0:Ratio42 , 1: sigma_Ratio42 , 2 : Signal Eff , 3: sigma of Signal Eff
  vector<double> R42 = Ratio42( 0 ) ; 

  vector<string> File4J ;
  fitInput->GetParameters( "FakeData", &File4J );

  vector<string> File2J ;
  fitInput->GetParameters( "2JSamples", &File2J );

  string phaseSmear ;
  fitInput->GetParameters( "PhaseSmear", &phaseSmear );
  bool smearing = ( phaseSmear == "ON" ) ? true : false ;

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD  = scaleQCD0* 3.  ;
  double scaleWJ   = scaleWJ0* 3. ;
  double scaleTt   = scaleTt0* 3.;
  cout<<" renew scaleWJ = "<< scaleWJ <<"  scaleQCD = "<< scaleQCD <<endl; 

  // run the execise - use the cuts in DataCard
  wmfitter->ResetCuts( -1, -1, -1, -1, -1, -1, true );
  bgCounter* count4jdata = new bgCounter();
  wmfitter->ReFitSolution( File4J[1], count4jdata, scaleWJ, NULL,  -3, smearing );
  vector<double> nW4JData= count4jdata->Output();
  wmfitter->ReFitSolution( File4J[2], count4jdata, scaleQCD, NULL, -3, smearing );
  vector<double> n4JData= count4jdata->Output();

  bgCounter* count4jSgData = new bgCounter();
  wmfitter->ReFitSolution( File4J[0], count4jSgData, scaleTt, NULL, -3, smearing );
  vector<double> n4JSgData= count4jSgData->Output();

  wmfitter->ResetCuts( 0, 300, 0, 350, 0, 999 );
  bgCounter* count2jdata = new bgCounter();
  wmfitter->ReFitSolution( File2J[1], count2jdata, scaleWJ,  NULL, -3, smearing );
  wmfitter->ReFitSolution( File2J[2], count2jdata, scaleQCD, NULL, -3, smearing );
  vector<double> n2JData = count2jdata->Output();
 
  cout<<" ========== M2M3 =============== "<<endl;
  cout<<" Data Background 4J M2M3 = "<< n4JData[0] <<" WJ= "<< nW4JData[0] <<" QCD= "<< n4JData[0] - nW4JData[0]<< endl;
  cout<<" Signal 4J M2M3 = "<< n4JSgData[0] << endl;
  cout<<" Data Background 2J = "<< n2JData[0] << endl;
  XSection( R42, n2JData[0], n4JData[0]+n4JSgData[0] ) ;

  cout<<" =============================== "<<endl;

  delete count4jdata;
  delete count2jdata;
  delete count4jSgData;

  cout<<" Count !! "<<endl;
}

vector<double> BgEstimation::XSection( vector<double>& R42, double n2J, double n4J ){

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

  vector<double> s_n4J = StatErr( n4J );
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

void BgEstimation::EnsembleTest( int randomSeed ){

  cout<<" @@ Ensemble Test Start @@ "<<endl;

  TString treeName = "muJets";

  // 4J Samples
  vector<string> flist;
  fitInput->GetParameters( "FakeData", &flist );

  bgCounter* sg = new bgCounter();
  vector<int> sglist = pseudoExp->GetEnsemble( flist[0], treeName, inputMean[0], randomSeed );
  wmfitter->ReFitSolution( flist[0], sg,  1, &sglist );
  vector<double> n4JTt = sg->Output();

  bgCounter* bg = new bgCounter();
  vector<int> bg1list = pseudoExp->GetEnsemble( flist[1], treeName, inputMean[1], randomSeed );
  wmfitter->ReFitSolution( flist[1], bg , 1, &bg1list );

  vector<int> bg2list = pseudoExp->GetEnsemble( flist[2], treeName, inputMean[2], randomSeed );
  wmfitter->ReFitSolution( flist[2], bg , 1, &bg2list );
  vector<double> n4JBg = bg->Output();

  // 2J Samples
  vector<string> File2J ;
  fitInput->GetParameters( "2JSamples", &File2J );

  bgCounter* bg2j = new bgCounter();

  vector<int> w2jlist = pseudoExp->GetEnsemble( File2J[1], treeName, inputMean[3], randomSeed );
  wmfitter->ReFitSolution( File2J[1], bg2j , 1, &w2jlist );

  vector<int> qcd2jlist = pseudoExp->GetEnsemble( File2J[2], treeName, inputMean[4], randomSeed );
  wmfitter->ReFitSolution( File2J[2], bg2j , 1, &qcd2jlist );

  vector<double> n2J = bg2j->Output();

  double n4J = n4JTt[0] + n4JBg[0] ;

  // get the ration42 and measure the cross-secton
  cout<<" Ntt = "<<n4JTt[0] <<" Nbg = "<< n4JBg[0] << endl;
  vector<double> R42 = Ratio42( 0 ) ;
  XSection( R42, n2J[0], n4J ) ;

  delete sg ;
  delete bg ;
  delete bg2j ;

  cout<<" Ensemble Test Done !!! "<<endl;
}

void BgEstimation::EnsembleTest( int nRun, int randomSeed ){

  cout<<" @@@@ Ensemble Test Start @@@@ "<<endl;

  // 4J Samples
  vector<TTree*> Tr4J = fitInput->GetForest( "FakeData", "muJets" );

  int tsz0 = Tr4J[0]->GetEntries() ;
  int tsz1 = Tr4J[1]->GetEntries() ;
  int tsz2 = Tr4J[2]->GetEntries() ;

  TH1D* gT4J = new TH1D("gT4J", " Tt+4J (mean ~ 25 )", 25, 0, 50 );
  TH1D* gW4J = new TH1D("gW4J", " W+4J (mean ~ 7 ) ",  15, 0, 30 );
  TH1D* gQ4J = new TH1D("gQ4J", " QCD4J (mean ~ 9 )",  15, 0, 30 );

  vector<double> n4JTt = RunEnsembles( tsz0, "", inputMean[0],  nRun, randomSeed, Tr4J[0], gT4J );
  vector<double> nWj   = RunEnsembles( tsz1, "", inputMean[1],  nRun, randomSeed, Tr4J[1], gW4J );
  vector<double> nQCD  = RunEnsembles( tsz2, "", inputMean[2], nRun, randomSeed, Tr4J[2], gQ4J );
  cout<<"  size of 4J tt:"<< tsz0 <<"  WJ:"<< tsz1 <<"  QCD:"<< tsz2 <<endl;

  // 2J Samples
  vector<TTree*> Tr2J = fitInput->GetForest( "2JSamples", "muJets" );

  int t2sz0 = Tr2J[0]->GetEntries() ;
  int t2sz1 = Tr2J[1]->GetEntries() ;
  int t2sz2 = Tr2J[2]->GetEntries() ;

  TH1D* gW2J = new TH1D("gW2J", " W+2J (mean = 513 ) ", 80, 300, 700 );
  TH1D* gQ2J = new TH1D("gQ2J", " QCD2J (mean = 503 )", 80, 300, 700 );

  wmfitter->ResetCuts( 0, 300, 0, 350, 0, 999 );
  
  vector<double> nW2j   = RunEnsembles( t2sz1, "", inputMean[3], nRun, randomSeed, Tr2J[1], gW2J );
  vector<double> nQCD2j = RunEnsembles( t2sz2, "", inputMean[4], nRun, randomSeed, Tr2J[2], gQ2J );
  
  /*
  double scaleWJ  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD = fitInput->NormalizeComponents( "qcd" );
  bgCounter* bg2j = new bgCounter();
  wmfitter->ReFitSolution( "", bg2j, scaleWJ*2, NULL,  -1, false, Tr2J[1] );
  vector<double> nW2J = bg2j->Output();
  wmfitter->ReFitSolution( "", bg2j, scaleQCD*2, NULL, -1, false, Tr2J[2] );
  vector<double> n2JEvts = bg2j->Output();
  double nQCD2J = n2JEvts[0] - nW2J[0] ;
  */

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
      vector<double> xtt = XSection( R42, n2J, n4J ) ;
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


vector<double> BgEstimation::RunEnsembles( int tsz, string fileName, double pMean, int nRun, int RandomSeed, TTree* theTree, TH1D* hGen ){

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
       bool smearing = ( nRepeat >= 0 ) ? true : false ;
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
	      wmfitter->ReFitSolution( fileName, counter , 1, &ensembles, 0, smearing, theTree );
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

vector<double> BgEstimation::StatErr( double m ){

  vector<double> pErr ;
  if ( m < 1. ) {
     pErr.push_back( -1*m ) ;
     pErr.push_back( m ) ;
  }
  else if ( m > 25. ) {
     pErr.push_back( -1*sqrt(m) ) ;
     pErr.push_back( sqrt(m) ) ;
  } 
  else {

     double step = 0.01 ;

     // -34%
     double k = m ;
     double lm = 0. ;
     double pp = 0. ;
     while (lm <= 0.34 || k < 0 ) {
          k = k - step ;
	  pp = TMath::Poisson( k, m );
	  lm = lm + (pp*step) ;
	  //cout<<" k = "<< k <<" , p = "<< pp <<" int_P = "<< lm <<endl;
     } 
     // +34%
     double j = m ;
     double hm = 0 ;
     double hp = 0 ;
     while ( hm <=0.34 || j < 0 ) {
           j = j + step ;
	   hp = TMath::Poisson( j, m );
	   hm = hm + (hp*step) ;
	   //cout<<" j = "<< j <<" , p = "<< hp <<" int_P = "<< hm <<endl;
     }
     pErr.push_back( k - m );
     pErr.push_back( j - m );
  }
  return pErr ;

}

vector<double> BgEstimation::ErrAovB( double A, double s_A, double B, double s_B ){

    vector<double> sA = StatErr( A ) ;
    double sAp = ( s_A != -1 ) ? s_A : sA[1]; 
    double sAn = ( s_A != -1 ) ? s_A : -1*sA[0]; 
    vector<double> sB = StatErr( B ) ;
    double sBp = ( s_B != -1 ) ? s_B : sB[1]; 
    double sBn = ( s_B != -1 ) ? s_B : -1*sB[0]; 

    double f = A / B ;
    double s_fp = sqrt( sAp*sAp + f*f*sBp*sBp ) / B ;
    double s_fn = sqrt( sAn*sAn + f*f*sBn*sBn ) / B ;
 
    vector<double> sf ;
    sf.push_back( s_fn );
    sf.push_back( s_fp );
    return sf ;

}

vector<double> BgEstimation::ErrAxB( double A, double s_A, double B, double s_B ){

    vector<double> sA = StatErr( A ) ;
    double sAp = ( s_A != -1 ) ? s_A : sA[1]; 
    double sAn = ( s_A != -1 ) ? s_A : -1*sA[0]; 
    vector<double> sB = StatErr( B ) ;
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

