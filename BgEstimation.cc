#include "BgEstimation.h"
#include "WFormat.h"

BgEstimation::BgEstimation( double massL, double massH ){

  fitFunc  = new MassFitFunction();
  fitInput = new MassAnaInput( "had", massL, massH );
  wmfitter  = new HadWMassFitter( massL, massH );
  pseudoExp = new PseudoExp( massL, massH );

}

BgEstimation::~BgEstimation(){

  delete fitFunc;
  delete fitInput ;
  delete wmfitter ;
  delete pseudoExp ;

}

// Exercise for background estimation
vector<double> BgEstimation::Ratio42( bool isHalf){

  cout<<" ---  Measuring Ratio(4J/2J) --- "<<endl;
  vector<string> File4J ;
  fitInput->GetParameters( "FakeData", &File4J );

  vector<string> File2J ;
  fitInput->GetParameters( "2JSamples", &File2J );

  vector<double> nEvts ;
  fitInput->GetParameters( "nEvents" , &nEvts );

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD  = ( isHalf ) ? scaleQCD0* 2.  : scaleQCD0 ;
  double scaleWJ   = ( isHalf ) ? scaleWJ0* 2.   : scaleWJ0  ;
  double scaleTt   = ( isHalf ) ? scaleTt0* 2.   : scaleWJ0  ;
  cout<<" renew scaleWJ = "<< scaleWJ <<"  scaleQCD = "<< scaleQCD <<endl; 

  int statIdx = ( isHalf ) ? -1 : 0 ;
  // ratio(4J/2J) 
  bgCounter* count4j = new bgCounter();
  wmfitter->ReFitSolution( File4J[1], count4j, scaleWJ,  NULL, statIdx );
  wmfitter->ReFitSolution( File4J[2], count4j, scaleQCD, NULL, statIdx );
  vector<double> n4JEvts = count4j->Output();
  cout<<" MC 4J M2M3 = "<< n4JEvts[0] <<endl;

  bgCounter* count2j = new bgCounter();
  wmfitter->ReFitSolution( File2J[1], count2j, scaleWJ, NULL,  statIdx );
  wmfitter->ReFitSolution( File2J[2], count2j, scaleQCD, NULL, statIdx );
  vector<double> n2JEvts = count2j->Output();
  cout<<" MC 2J M2M3= "<< n2JEvts[0] << endl;

  vector<double> ratio42;
  ratio42.push_back(  n4JEvts[0] / n2JEvts[0]  ) ;
  double sR = sqrt( n4JEvts[1] + ratio42[0]*ratio42[0]*n2JEvts[1]  ) / n2JEvts[0] ;
  ratio42.push_back( sR ) ;
  cout<<" MC R4/2 = "<< ratio42[0] <<" +/- "<< ratio42[1]   <<endl;

  // signal efficiency 
  bgCounter* count4jSg = new bgCounter();
  wmfitter->ReFitSolution( File4J[0], count4jSg, scaleTt, NULL, statIdx );
  vector<double> n4JSg= count4jSg->Output();
  double  Eff = n4JSg[0] / (scaleTt0* nEvts[0]) ;
  double sEff = sqrt(  n4JSg[1] +  (Eff*Eff*scaleTt0*nEvts[0]) ) / (scaleTt0*nEvts[0]) ;
  ratio42.push_back( Eff );
  ratio42.push_back( sEff );
  cout<<" Signal MC 4J M2M3= "<< n4JSg[0] <<" Eff = "<< Eff <<" +/- "<< sEff <<  endl;

  delete count4j;
  delete count2j;
  delete count4jSg;

  cout<<" Ratio Measured !! "<<endl;

  return ratio42 ;  
}

void BgEstimation::MethodTest(){

  // 0:Ratio42 , 1: sigma_Ratio42 , 2 : Signal Eff , 3: sigma of Signal Eff
  vector<double> R42 = Ratio42( true ) ; 

  vector<string> File4J ;
  fitInput->GetParameters( "FakeData", &File4J );

  vector<string> File2J ;
  fitInput->GetParameters( "2JSamples", &File2J );

  double scaleTt0  = fitInput->NormalizeComponents( "tt" );
  double scaleQCD0 = fitInput->NormalizeComponents( "qcd" );
  double scaleWJ0  = fitInput->NormalizeComponents( "wj" );
  double scaleQCD  = scaleQCD0* 2.  ;
  double scaleWJ   = scaleWJ0* 2. ;
  double scaleTt   = scaleTt0* 2.;
  cout<<" renew scaleWJ = "<< scaleWJ <<"  scaleQCD = "<< scaleQCD <<endl; 

  // run the execise
  bgCounter* count4jdata = new bgCounter();
  wmfitter->ReFitSolution( File4J[1], count4jdata, scaleWJ, NULL,  1 );
  wmfitter->ReFitSolution( File4J[2], count4jdata, scaleQCD, NULL, 1 );
  vector<double> n4JData= count4jdata->Output();

  bgCounter* count2jdata = new bgCounter();
  wmfitter->ReFitSolution( File2J[1], count2jdata, scaleWJ,  NULL, 1 );
  wmfitter->ReFitSolution( File2J[2], count2jdata, scaleQCD, NULL, 1 );
  vector<double> n2JData = count2jdata->Output();
 
  bgCounter* count4jSgData = new bgCounter();
  wmfitter->ReFitSolution( File4J[0], count4jSgData, scaleTt, NULL, 1 );
  vector<double> n4JSgData= count4jSgData->Output();

  cout<<" ========== M2M3 =============== "<<endl;
  cout<<" Data Background 4J M2M3 = "<< n4JData[0] <<" Signal 4J M2M3 = "<< n4JSgData[0] << endl;
 
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

  cout<<" ========== Ratio 4J/2J  =============== "<<endl;
  cout<<" MC R4/2      = "<< R42[0] <<" +/- "<< R42[1]   <<endl;
  cout<<" Data 2J M2M3 = "<< n2J <<" +/- "<< sqrt( n2J ) <<endl;
  double  Exp4JBG = R42[0]*n2J ;
  //double sExp4JBG = R42[0]*R42[0]*n2J + n2J*n2J*R42[1]*R42[1] ;
  vector<double> sExp4JBG = ErrAxB ( R42[0], R42[1], n2J, -1 );
  cout<<" Expected 4J Background M2M3 = "<< Exp4JBG <<" - "<<  sExp4JBG[0] <<" + "<< sExp4JBG[1] << endl;

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

  double inputMean_tt = 18.57*1.365 ;
  double inputMean_wj = 4.985*1.308 ;
  double inputMean_qcd = 8.331*1.077 ;

  bgCounter* sg = new bgCounter();
  vector<int> sglist = pseudoExp->GetEnsemble( flist[0], treeName, inputMean_tt, randomSeed );
  wmfitter->ReFitSolution( flist[0], sg,  1, &sglist );
  vector<double> n4JTt = sg->Output();

  bgCounter* bg = new bgCounter();
  vector<int> bg1list = pseudoExp->GetEnsemble( flist[1], treeName, inputMean_wj, randomSeed );
  wmfitter->ReFitSolution( flist[1], bg , 1, &bg1list );

  vector<int> bg2list = pseudoExp->GetEnsemble( flist[2], treeName, inputMean_qcd, randomSeed );
  wmfitter->ReFitSolution( flist[2], bg , 1, &bg2list );
  vector<double> n4JBg = bg->Output();

  // 2J Samples
  vector<string> File2J ;
  fitInput->GetParameters( "2JSamples", &File2J );

  bgCounter* bg2j = new bgCounter();

  double mean_w2j =  71.3*7.203 ;
  double mean_qcd2j = 109.9*4.58 ;

  vector<int> w2jlist = pseudoExp->GetEnsemble( File2J[1], treeName, mean_w2j, randomSeed );
  wmfitter->ReFitSolution( File2J[1], bg2j , 1, &w2jlist );

  vector<int> qcd2jlist = pseudoExp->GetEnsemble( File2J[2], treeName, mean_qcd2j, randomSeed );
  wmfitter->ReFitSolution( File2J[2], bg2j , 1, &qcd2jlist );

  vector<double> n2J = bg2j->Output();

  double n4J = n4JTt[0] + n4JBg[0] ;

  // get the ration42 and measure the cross-secton
  cout<<" Ntt = "<<n4JTt[0] <<" Nbg = "<< n4JBg[0] << endl;
  vector<double> R42 = Ratio42( true ) ;
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

  double inputMean_tt  = 18.57*1.365 ;
  double inputMean_wj  = 4.985*1.308 ;
  double inputMean_qcd = 8.331*1.077 ;

  int tsz0 = Tr4J[0]->GetEntries() ;
  int tsz1 = Tr4J[1]->GetEntries() ;
  int tsz2 = Tr4J[2]->GetEntries() ;
  cout<<"  size of 4J tt:"<< tsz0 <<"  WJ:"<< tsz1 <<"  QCD:"<< tsz2 <<endl;

  vector<double> n4JTt = RunEnsembles( tsz0, "", inputMean_tt,  nRun, randomSeed, Tr4J[0] );
  vector<double> nWj   = RunEnsembles( tsz1, "", inputMean_wj,  nRun, randomSeed, Tr4J[1] );
  vector<double> nQCD  = RunEnsembles( tsz2, "", inputMean_qcd, nRun, randomSeed, Tr4J[2] );

  // 2J Samples
  vector<TTree*> Tr2J = fitInput->GetForest( "2JSamples", "muJets" );

  bgCounter* bg2j = new bgCounter();

  double mean_w2j   = 71.3*7.203 ;
  double mean_qcd2j = 109.9*4.58 ;

  int t2sz0 = Tr2J[0]->GetEntries() ;
  int t2sz1 = Tr2J[1]->GetEntries() ;
  int t2sz2 = Tr2J[2]->GetEntries() ;
  cout<<"  size of 2J tt:"<< t2sz0 <<"  WJ:"<< t2sz1 <<"  QCD:"<< t2sz2 <<endl;

  vector<double> nW2j   = RunEnsembles( t2sz1, "", mean_w2j, nRun, randomSeed, Tr2J[1] );
  vector<double> nQCD2j = RunEnsembles( t2sz2, "", mean_qcd2j, nRun, randomSeed, Tr2J[2] );

  // get the ration42 and measure the cross-secton
  vector<double> R42 = Ratio42( true ) ;
  double n4JBg = 0; 
  double n2J   = 0;
  double n4J   = 0;

  TH1D* hXtt = new TH1D("hXtt", " tt cross-section @ 5/pb ", 24, 0, 360 );
  TH1D* hsXn = new TH1D("hsXn", " -Error of tt x-section @ 5/pb ", 20, 0, 100 );
  TH1D* hsXp = new TH1D("hsXp", " +Error of tt x-section @ 5/pb ", 20, 0, 100 );

  TH1D* hB4J = new TH1D("hB4J", " Expected BG+4J ", 10, 4, 24 );
  TH1D* hT4J = new TH1D("hT4J", " Tt+4J (mean = 19 )", 18, 0, 36 );
  TH1D* hW2J = new TH1D("hW2J", " W+2J (mean = 71 ) ", 20, 20, 120 );
  TH1D* hQ2J = new TH1D("hQ2J", " QCD2J (mean = 110 )", 30, 50, 200 );
  TH1D* hW4J = new TH1D("hW4J", " W+4J (mean = 5 ) ", 8, 0, 16 );
  TH1D* hQ4J = new TH1D("hQ4J", " QCD4J (mean = 8 )", 10, 0, 20 );

  for ( size_t i=0; i< nRun ; i++) {
      n4JBg = nWj[i]   + nQCD[i]   ;
      n2J   = nW2j[i]  + nQCD2j[i] ;
      n4J   = n4JTt[i] + n4JBg     ;
      cout<<" ********** Run "<< i <<" ******************"<<endl;
      cout<<"  Ntt = "<<n4JTt[i] <<" Nbg = "<< n4JBg << endl;
      vector<double> xtt = XSection( R42, n2J, n4J ) ;
      hXtt->Fill( xtt[0] );
      hsXn->Fill( xtt[1] );
      hsXp->Fill( xtt[2] );

      hB4J->Fill( xtt[3] );
      hT4J->Fill( n4JTt[i] );
      hW2J->Fill( nW2j[i] );
      hQ2J->Fill( nQCD2j[i] );
      hW4J->Fill( nWj[i] );
      hQ4J->Fill( nQCD[i] );
  }

  gStyle->SetOptStat("neirom");
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
  c1->Print(  "NoJESNoB_4JPF/ttXsection.gif"  );

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
  c2->Print(  "NoJESNoB_4JPF/ttXsectionErr.gif"  );

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
  c3->Print(  "NoJESNoB_4JPF/BGPoisson.gif"  );

  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  c4->Divide(1,2);
  c4->cd(1);
  hT4J->Draw();
  c4->Update();
  c4->cd(2);
  hB4J->Draw();
  c4->Update();
  c4->Print(  "NoJESNoB_4JPF/4JPoisson.gif"  );

  cout<<" !!! Ensemble Test Done !!! "<<endl;
}

// under construction 
vector<double> BgEstimation::RunEnsembles( int tsz, string fileName, double pMean, int nRun, int RandomSeed, TTree* theTree ){

  cout<<"  <<< Getting Ensembles >>> "<<endl;

  // event solution tree
  /*
  TFile*  file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName , "muJets", file );
  int tsz = tr1->GetEntries();
  cout<<" All Evt Size = "<< tsz << endl;
  delete tr1 ;
  delete file ;
  */

  // set up the random number function
  TRandom* tRan = new TRandom();
  tRan->SetSeed( RandomSeed );

  // shuffle the events
  //vector<int> evtline = pseudoExp->EventShuffle( tsz, RandomSeed );
  double nRepeat_ = ( pMean*nRun ) / static_cast<double>( tsz ) ;
  int nRepeat = ( nRepeat_ > 1. ) ? static_cast<int>(nRepeat_)+1 : 1 ;
  vector<int> evtline;
  for (int i =0; i < nRepeat; i++ ) {
      vector<int> evtmp = pseudoExp->EventShuffle( tsz, RandomSeed );
      for ( size_t j=0; j< evtmp.size(); j++) {
          evtline.push_back( evtmp[j] ) ;
      }
  }
  cout<<" nRepeat = "<< nRepeat <<" -> "<< tsz << " evtline size = "<< evtline.size() <<" / "<< pMean <<endl ;

  bool smearing = ( nRepeat > 1 ) ? true : false ;
  vector<double> nCountV ;
  vector<int> ensembles ;
  bool nextRun = true ;
  int  RunStop = 0 ;
  int  NEvts = 0 ;
  int  nRun_ = 0 ;
  for (int k=0; k < evtline.size(); k++) {
 
      if ( nRun_ == nRun ) break; 
      if ( nextRun ) {
         int PoiSeed = tRan->Integer( 1000+k );
         tRan->SetSeed( PoiSeed );
         NEvts = tRan->Poisson( pMean );
         //cout<<" ("<<k<<") N of Evts needed = "<< NEvts <<endl ;
         nextRun = false ;
         RunStop = k + NEvts ;
         ensembles.clear() ;
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
      }

  }

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

