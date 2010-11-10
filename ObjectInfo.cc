#include "ObjectInfo.h"

ObjectInfo::ObjectInfo() {

  fitInput  = new MassAnaInput();
  pseudoExp = new PseudoExp();

  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "JetCuts",  &jetCuts );
  fitInput->GetParameters( "MuonCuts", &muonCuts );
  fitInput->GetParameters( "n_Jets",   &n_Jets );
  fitInput->GetParameters( "PlotType", &plotType );
  fitInput->GetParameters( "Inclusive", &Inclusive );
  inclu = ( Inclusive == "YES" ) ? true : false ;
  if (  inclu ) cout<<" Inclusive !!" <<endl;
  if ( !inclu ) cout<<" Exclusive !!" <<endl;
  isWeight = false ;
  
  bgMtCut  = 35;
  bgMETCut = 20;

}

ObjectInfo::~ObjectInfo(){

  delete fitInput ;
  delete pseudoExp;

}

// inclusive = true :  >=njet  
void ObjectInfo::EvtSelector( string fileName, recoObj* histos, bool smearing, double scale, vector<int>* evtlist, int evtSplit ) {

  TTree* tr = fitInput->TreeMap( fileName );

  // this is only for MC Matching
  double jpx[15],jpy[15],jpz[15],jE[15], bDis[15];
  int    n90[10] , nHits[2];
  double mpx[2],mpy[2],mpz[2],mE[2], mIso[2], d0[2], X2[2] ;
  double npx[3],npy[3],npz[3],nE[3] ;
  int nNu, nMu, nJ ;
  tr->SetBranchAddress("jpx"    ,&jpx);
  tr->SetBranchAddress("jpy"    ,&jpy);
  tr->SetBranchAddress("jpz"    ,&jpz);
  tr->SetBranchAddress("jE"     ,&jE);
  tr->SetBranchAddress("nJ"     ,&nJ);
  tr->SetBranchAddress("n90"    ,&n90);
  tr->SetBranchAddress("mpx"    ,&mpx);
  tr->SetBranchAddress("mpy"    ,&mpy);
  tr->SetBranchAddress("mpz"    ,&mpz);
  tr->SetBranchAddress("mE"     ,&mE);
  tr->SetBranchAddress("mIso"   ,&mIso);
  tr->SetBranchAddress("nHits"  ,&nHits);
  tr->SetBranchAddress("X2"     ,&X2);
  tr->SetBranchAddress("d0"     ,&d0);
  tr->SetBranchAddress("nHits"  ,&nHits);
  tr->SetBranchAddress("npx"    ,&npx);
  tr->SetBranchAddress("npy"    ,&npy);
  tr->SetBranchAddress("npz"    ,&npz);
  tr->SetBranchAddress("nE"     ,&nE);
  tr->SetBranchAddress("bTh"    ,&bDis );
  tr->SetBranchAddress("nNu"    ,&nNu );
  tr->SetBranchAddress("nMu"    ,&nMu );

  int idx = 0;
  int loopEnd = ( evtlist == NULL ) ? tr->GetEntries() : evtlist->size() ;

  //int fk = 0;
  //int ffk = 0;
  //int fffk = 0;
  for (int k=0; k< loopEnd ; k++) {
      if ( evtlist != NULL && evtlist->size() == 0 ) break;
      idx = ( evtlist == NULL ) ? k : (*evtlist)[k] ;

      if ( evtSplit == -1 && (k%2) == 0 ) continue ;
      if ( evtSplit ==  1 && (k%2) == 1 ) continue ;
      if ( evtSplit >  1 && (k%evtSplit ) != 0 ) continue ;
      if ( evtSplit < -1 && (k%evtSplit ) != 1 ) continue ;

      tr->GetEntry( idx );
      if ( nJ < n_Jets ) continue ;
      
      TLorentzVector m0( mpx[0], mpy[0], mpz[0], mE[0] );
      if ( m0.Pt()  < muonCuts[0] ) continue ;
      if ( m0.Eta() > muonCuts[1] ) continue ;
      if ( mIso[0]  > muonCuts[2] ) continue ;

      vector<TLorentzVector> objlist;
      int NCountedJets = 0;
      for ( int i = 0; i< nJ ; i ++) {
          TLorentzVector jn( jpx[i], jpy[i], jpz[i], jE[i] );
          // sync cuts
          if ( jn.Pt()          < jetCuts[0] ) continue ;
          if ( fabs( jn.Eta() ) > jetCuts[1] ) continue ;
          objlist.push_back( jn );
          NCountedJets++ ;
      }
      //cout<<" nJets? = "<< NCountedJets <<" sz = "<< objlist.size()<<" nJ = "<< nJ <<endl;
      if ( NCountedJets != n_Jets && !inclu ) continue ;
      if ( NCountedJets > 6 || NCountedJets < n_Jets ) continue ;

      while ( objlist.size() < 4 ) {
            TLorentzVector j0( 0., 0., 0., 0. );
            objlist.push_back( j0 );
      }
      TLorentzVector n0( npx[0], npy[0], npz[0], nE[0] );
      objlist.push_back( m0 );
      objlist.push_back( n0 );
      //if ( npz[0] == 0 ) fk++ ;
      //if ( npz[0] != 0 ) ffk++ ;
      //if ( nNu == 1 ) fffk++ ;

      if ( smearing )  {
         pseudoExp->PhaseSmearing( objlist, 0 );
         pseudoExp->JetEtReSort( objlist );
      }

      int sz = objlist.size() ;
      //cout<<" obj size = "<< sz <<endl ;
      histos->gethad( objlist[0], objlist[1], objlist[2] );   
      histos->getlep( objlist[3], objlist[sz-2], objlist[sz-1] );  
      double muonQv[3] = { mIso[0], d0[0], X2[0] };
      int intQv[2] = { NCountedJets, nHits[0] };
      histos->getFloats( muonQv );
      histos->getIntegrals( intQv );
      //histos->getOther( NCountedJets, nHits[0], mIso[0], d0[0], X2[0] ) ;
      //
      //TLorentzVector vM2 = objlist[sz-2] + objlist[sz-1] ;
      //double vPt = vM2.Pt() ;
      //double MtdPhi = m0.DeltaPhi( n0 );
      //if ( fabs(MtdPhi) < 0.5 ) continue ;
      //double weight = ( isWeight ) ? EvtScaling(vPt, fileName ) : 1. ;
      double weight = ( isWeight ) ? EvtScaling( NCountedJets, fileName ) : 1. ;
      histos->Fillh( weight, scale ) ;
  }
  //cout<<" neu fk = "<<fk <<" ffk = "<< ffk <<"  -> "<<  fffk <<endl;

}

// mode = 1 : measuring MC normalization, mode = 2 : measuring the ratio
void ObjectInfo::QCDSelector( string fileName, recoObj* histos, bool smearing, double scale, bool doQCD, int mode, int evtSplit ) {

  TTree* tr = fitInput->TreeMap( fileName );

  // this is only for MC Matching
  double jpx[10],jpy[10],jpz[10],jE[10], bDis[10];
  int    n90[10] , nHits[2];
  double mpx[2],mpy[2],mpz[2],mE[2], mIso[2], d0[2], X2[2] ;
  double npx[2],npy[2],npz[2],nE[2] ;
  int nNu, nMu, nJ ;
  tr->SetBranchAddress("jpx"    ,&jpx);
  tr->SetBranchAddress("jpy"    ,&jpy);
  tr->SetBranchAddress("jpz"    ,&jpz);
  tr->SetBranchAddress("jE"     ,&jE);
  tr->SetBranchAddress("nJ"     ,&nJ);
  tr->SetBranchAddress("n90"    ,&n90);
  tr->SetBranchAddress("mpx"    ,&mpx);
  tr->SetBranchAddress("mpy"    ,&mpy);
  tr->SetBranchAddress("mpz"    ,&mpz);
  tr->SetBranchAddress("mE"     ,&mE);
  tr->SetBranchAddress("mIso"   ,&mIso);
  tr->SetBranchAddress("nHits"  ,&nHits);
  tr->SetBranchAddress("X2"     ,&X2);
  tr->SetBranchAddress("d0"     ,&d0);
  tr->SetBranchAddress("nHits"  ,&nHits);
  tr->SetBranchAddress("npx"    ,&npx);
  tr->SetBranchAddress("npy"    ,&npy);
  tr->SetBranchAddress("npz"    ,&npz);
  tr->SetBranchAddress("nE"     ,&nE);
  tr->SetBranchAddress("bTh"    ,&bDis );
  tr->SetBranchAddress("nNu"    ,&nNu );
  tr->SetBranchAddress("nMu"    ,&nMu );
  cout<<" Current Mt cut = "<< bgMtCut <<"  MET cut = "<< bgMETCut << endl; 

  for (int k=0; k< tr->GetEntries() ; k++) {

      if ( evtSplit == -1 && (k%2) == 0 ) continue ;
      if ( evtSplit ==  1 && (k%2) == 1 ) continue ;
      if ( evtSplit >  1 && (k%evtSplit ) != 0 ) continue ;
      if ( evtSplit < -1 && (k%evtSplit ) != 1 ) continue ;

      tr->GetEntry(k);

      if ( nJ < n_Jets ) continue ;
      
      TLorentzVector m0( mpx[0], mpy[0], mpz[0], mE[0] );
      if ( m0.Pt()  < muonCuts[0] ) continue ;
      if ( m0.Eta() > muonCuts[1] ) continue ;
      if ( mIso[0]  > muonCuts[2] ) continue ;

      vector<TLorentzVector> objlist;
      int NCountedJets = 0;
      double dRmin = 999 ;
      double ptRel = 999 ;
      for ( int i = 0; i< nJ ; i ++) {
          TLorentzVector jn( jpx[i], jpy[i], jpz[i], jE[i] );
          // sync cuts
          if ( jn.Pt()          < jetCuts[0] ) continue ;
          if ( fabs( jn.Eta() ) > jetCuts[1] ) continue ;
          objlist.push_back( jn );
          NCountedJets++ ;
          double dRmj = m0.DeltaR( jn );
          if ( dRmj < dRmin ) {
             dRmin = dRmj ;
             ptRel = PtRel( m0, jn );
          }
      }
      //cout<<" nJets? = "<< NCountedJets <<" sz = "<< objlist.size()<<" nJ = "<< nJ <<endl;
      if ( NCountedJets != n_Jets && !inclu ) continue ;
      if ( NCountedJets > 6 || NCountedJets < n_Jets ) continue ;

      while ( objlist.size() < 4 ) {
            TLorentzVector j0( 0., 0., 0., 0. );
            objlist.push_back( j0 );
      }
      TLorentzVector n0( npx[0], npy[0], npz[0], nE[0] );
      objlist.push_back( m0 );
      objlist.push_back( n0 );

      double dPhi = m0.DeltaPhi( n0 ) ;
      double theMt = sqrt( 2.*m0.Pt()*n0.Pt()*( 1. - cos(dPhi) ) );

      //TLorentzVector vM2 = n0 + m0 ;
      //double vPt = vM2.Pt() ;

      if ( n0.Pt() > bgMETCut && doQCD ) continue;
      if ( theMt   > bgMtCut  && doQCD ) continue ;
      if ( n0.Pt() < bgMETCut && !doQCD ) continue;
      if ( theMt   < bgMtCut  && !doQCD ) continue ;

      //if ( n0.Pt() < 25. && theMt < 35 && !doQCD ) continue;

      if ( smearing )  {
         pseudoExp->PhaseSmearing( objlist, 0 );
         pseudoExp->JetEtReSort( objlist );
      }

      int sz = objlist.size() ;
      //cout<<" obj size = "<< sz <<endl ;
      histos->gethad( objlist[0], objlist[1], objlist[2] );   
      histos->getlep( objlist[3], objlist[sz-2], objlist[sz-1] );  
      double muonQv[5] = { mIso[0], d0[0], X2[0], dRmin, ptRel };
      int intQv[2] = { NCountedJets, nHits[0] };
      histos->getFloats( muonQv );
      histos->getIntegrals( intQv );

      double weight = 1. ;
      // used for normalization study
      if ( isWeight && mode == 1 ) weight = QCDScaling( NCountedJets, fileName ) ;
      // used for Ratio study
      if ( isWeight && mode == 2 ) weight = EvtScaling( NCountedJets, fileName ) ;

      histos->Fillh( weight, fabs(scale) ) ;
  }
  
}

double ObjectInfo::MtFitter( TH1D* hs1, TF1& fjb1, TCanvas* c1, double minMt, double maxMt, bool doPlot ){

  TF1* fJb = new TF1("fJb", MassFitFunction::ConvJacob, 0, 160, 7);

  // Leptonic W transverse mass 
  fJb->SetParLimits(0,  0.,  10000.);
  fJb->SetParLimits(1,  0.,  10.);
  fJb->SetParLimits(2, -50.,  50.);
  fJb->SetParLimits(3,  0.0001,  1000.);
  fJb->SetParLimits(4,  45.4,    115.4);
  fJb->SetParLimits(5,  0.,      50.);
  fJb->SetParLimits(6,  0.0001,   100.);

  std::vector<double> PJ ;
  fitInput->GetParameters( "PJacob", &PJ );
  std::vector<int> PS ;
  fitInput->GetParameters( "PSwitch", &PS );
  
  if ( PS[0] == 1 ) fJb->SetParameter(0,  PJ[0] );
  if ( PS[1] == 1 ) fJb->SetParameter(1,  PJ[1] );
  if ( PS[2] == 1 ) fJb->SetParameter(2,  PJ[2] );
  if ( PS[3] == 1 ) fJb->SetParameter(3,  PJ[3] );
  if ( PS[4] == 1 ) fJb->SetParameter(4,  PJ[4] );
  if ( PS[5] == 1 ) fJb->SetParameter(5,  PJ[5] );
  if ( PS[6] == 1 ) fJb->SetParameter(6,  PJ[6] );
  
  if ( PS[0] == 0 ) fJb->FixParameter(0,  PJ[0] );
  if ( PS[1] == 0 ) fJb->FixParameter(1,  PJ[1] );
  if ( PS[2] == 0 ) fJb->FixParameter(2,  PJ[2] );
  if ( PS[3] == 0 ) fJb->FixParameter(3,  PJ[3] );
  if ( PS[4] == 0 ) fJb->FixParameter(4,  PJ[4] );
  if ( PS[5] == 0 ) fJb->FixParameter(5,  PJ[5] );
  if ( PS[6] == 0 ) fJb->FixParameter(6,  PJ[6] );

  hs1->Fit( fJb, "R","", minMt, maxMt );

  // fast integral
  double nWJ = 0 ; 
  double lowMt =   0;
  double upMt  = 140;
  double step = ( upMt - lowMt )/100. ;
  for (int i=0; i<100; i++){
      double theMt = lowMt + step*( i + 0.5 ) ; 
      nWJ += step*fJb->Eval( theMt );
  }
  nWJ = nWJ / hs1->GetBinWidth(1) ;
  cout<<" The W+Jets content = "<< nWJ <<" / Total = "<< hs1->Integral() << endl ;

  // plotting the result
  if ( doPlot ) {
     gStyle->SetOptFit(111);
     gStyle->SetOptStat("ieo");
     TCanvas* c1a = new TCanvas("c1a","", 800, 600);
     c1a->SetGrid();
     c1a->SetFillColor(10);
     c1a->SetFillColor(10);
     c1a->cd();

     hs1->SetLineWidth(2);
     hs1->Draw();
     fJb->SetLineColor(2);
     fJb->Draw("sames");
     
     c1a->Update();

     TString plotName1a = hfolder + "hObjects/FitJacobianTest."+plotType ;
     c1a->Print( plotName1a );
     delete c1a;
  }
  if ( c1 != NULL ) {
     c1->cd();
     fJb->SetLineColor(2);
     fJb->Draw("sames");
     c1->Update();
  }

  fjb1 = *fJb ;

  delete fJb ;
  return nWJ ;

}

void ObjectInfo::ObjHistoPlotter( string fileName, bool smearing,  bool doWMtfit ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hObjects" );
  gSystem->cd( "../" );

  double scale = 1. ;
  if ( fileName.substr(0,2) == "tt" )  scale = fitInput->NormalizeComponents( "tt" );
  if ( fileName.substr(0,2) == "wj" )  scale = fitInput->NormalizeComponents( "wj" );
  if ( fileName.substr(0,2) == "qc" )  scale = fitInput->NormalizeComponents( "qcd" );
  if ( fileName.substr(0,2) == "zj" )  scale = fitInput->NormalizeComponents( "zj" );
  if ( fileName.substr(0,2) == "tq" )  scale = fitInput->NormalizeComponents( "tq" );
  if ( fileName.substr(0,2) == "tw" )  scale = fitInput->NormalizeComponents( "tw" );

  hObjs* hs = new hObjs( fileName.substr(0,2) ) ;
  EvtSelector( fileName, hs, smearing, scale );

  vector<TH1D*> h1Ds ;
  hs->Fill1DVec( h1Ds );

  gStyle->SetOptStat("ieruom");
  gStyle->SetLabelSize( 0.06, "X");
  gStyle->SetLabelSize( 0.06, "Y");
  gStyle->SetLabelSize( 0.06, "Z");

  // Jet Et Spectrum
  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetLogy();
  c1->cd();

  //gStyle->SetTextSize(0.6);
  //gStyle->SetStatFontSize(0.6);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(7);
  h1Ds[3]->SetLineColor(7);
  h1Ds[3]->SetLineWidth(2);
  h1Ds[3]->Draw();
  c1->Update();

  gStyle->SetStatX(0.75);
  gStyle->SetStatTextColor(6);
  h1Ds[2]->SetLineColor(6);
  h1Ds[2]->SetLineWidth(2);
  h1Ds[2]->DrawCopy("sames");
  c1->Update();

  gStyle->SetStatY(0.70);
  gStyle->SetStatX(0.95);
  gStyle->SetStatTextColor(4);
  h1Ds[1]->SetLineColor(4);
  h1Ds[1]->SetLineWidth(2);
  h1Ds[1]->DrawCopy("sames");
  c1->Update();

  gStyle->SetStatX(0.75);
  gStyle->SetStatTextColor(2);
  h1Ds[0]->SetLineColor(2);
  h1Ds[0]->SetLineWidth(2);
  h1Ds[0]->DrawCopy("sames");
  c1->Update();

  TString plotname1 = hfolder + "hObjects/"+ fileName + "_JetEt."+plotType ;
  c1->Print( plotname1 );

  // Muon and MET Spectrum
  TCanvas* c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->cd();
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(2);
  h1Ds[4]->SetLineColor(2);
  h1Ds[4]->SetLineWidth(2);
  h1Ds[4]->Draw();
  c2->Update();

  gStyle->SetStatY(0.7);
  gStyle->SetStatTextColor(4);
  h1Ds[5]->SetLineColor(4);
  h1Ds[5]->SetLineWidth(2);
  h1Ds[5]->DrawCopy("sames");
  c2->Update();

  TString plotname2 = hfolder + "hObjects/"+ fileName + "_LepEt."+plotType ;
  c2->Print( plotname2 );


  // Muon Eta distribution
  TCanvas* c3 = new TCanvas("c3","", 800, 600);
  c3->SetGrid();
  c3->SetFillColor(10);
  c3->SetFillColor(10);
  c3->cd();

  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(2);
  h1Ds[6]->SetLineColor(1);
  h1Ds[6]->SetLineWidth(2);
  h1Ds[6]->Draw();
  c3->Update();

  TString plotname3 = hfolder + "hObjects/"+ fileName + "_MuonEta."+plotType ;
  c3->Print( plotname3 );

  // Leptonic W transverse mass 
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("ieo");
  TCanvas* c4 = new TCanvas("c4","", 800, 600);
  c4->SetGrid();
  c4->SetFillColor(10);
  c4->SetFillColor(10);
  c4->cd();

  h1Ds[13]->SetLineWidth(2);
  h1Ds[13]->Draw();

  if ( doWMtfit ) {
     TF1* fjb1 = new TF1("fjb1", MassFitFunction::ConvJacob, 0, 160, 7);
     MtFitter( h1Ds[13], *fjb1, c4, 5, 120, true ) ;
     fjb1->SetLineColor(2);
     fjb1->Draw("sames");
     delete fjb1;
  }
  c4->Update();

  TString plotname4 = hfolder + "hObjects/"+ fileName + "_LepWMt."+plotType ;
  c4->Print( plotname4 );

  TCanvas* c5 = new TCanvas("c5","", 800, 600);
  c5->SetGrid();
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->cd();

  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(2);
  h1Ds[12]->SetLineColor(1);
  h1Ds[12]->SetLineWidth(2);
  h1Ds[12]->Draw();
  c5->Update();

  TString plotname5 = hfolder + "hObjects/"+ fileName + "_lepW_Pt."+plotType ;
  c5->Print( plotname5 );

  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
}

void ObjectInfo::JacobTester(){

  std::vector<double> PJ ;
  fitInput->GetParameters( "PJacob", &PJ );

  TCanvas* c5 = new TCanvas("c5","", 800, 600);
  c5->SetGrid();
  c5->SetFillColor(10);
  c5->SetFillColor(10);
  c5->cd();

  TF1 *fJb1   = new TF1("fJb1", MassFitFunction::ConvJacob, 0, 140, 7);
  //TF1 *fJb1   = new TF1("fJb1", MassFitFunction::fitJacob, 0., 140., 5);
  fJb1->FixParameter(0, PJ[0] );
  fJb1->FixParameter(1, PJ[1] );
  fJb1->FixParameter(2, PJ[2] );
  fJb1->FixParameter(3, PJ[3] );
  fJb1->FixParameter(4, PJ[4] );
  fJb1->FixParameter(5, PJ[5] );
  fJb1->FixParameter(6, PJ[6] );
 
  fJb1->Draw();
  c5->Update();
  TString plotname5 = hfolder + "hObjects/Jacobian."+plotType ;
  c5->Print( plotname5 );

  delete c5;
  delete fJb1 ;
}

void ObjectInfo::DataPlotter( string DataName, vector<string>& fakeData, bool doScale ){

  TString theFolder = hfolder ;
  TString theSubFolder = ( doScale == false ) ?  "hData18/" : "hData18_Scale/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 30 ;

  hObjs* hdt = new hObjs("data", 30) ;
  EvtSelector( DataName, hdt );
  vector<TH1D*> h_dt ;
  hdt->Fill1DVec( h_dt );

  if ( h_dt[13]->Integral() < 500. ) {
     nbin = 15 ;
     h_dt[13]->Rebin(2);
     h_dt[14]->Rebin(2);
     h_dt[15]->Rebin(2);
  }

  //vector<string> flist;
  //fitInput->GetParameters( "FakeData", &flist );
  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  Reset( 1, doScale ) ;

  hObjs* htt = new hObjs("tt", nbin ) ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  EvtSelector( fakeData[2], hqcd, false, scale2 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  EvtSelector( fakeData[3], hzj, false, scale3 );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  EvtSelector( fakeData[4], htq, false, scale4 );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  EvtSelector( fakeData[5], htw, false, scale5 );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  gStyle->SetOptStat("ieruom");
  //gStyle->SetOptStat("");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");
  //gStyle->SetOptTitle(0);

  TString hNames[18] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",    "Mt_4" } ;


  //TF1* fjb1 = new TF1("fjb1", MassFitFunction::ConvJacob, 0, 140, 7);

  for (int i=0; i<18; i++) {
      if ( i>0 && i <4 ) continue;
      THStack* hStk = new THStack("hStk", hNames[i] );
      h_wj[i]->SetFillColor(kGreen);
      hStk->Add( h_wj[i] );
      h_zj[i]->SetFillColor(kAzure-2);
      hStk->Add( h_zj[i] );
      h_tq[i]->SetFillColor(kMagenta+2);
      hStk->Add( h_tq[i] );
      h_tw[i]->SetFillColor(kMagenta);
      hStk->Add( h_tw[i] );
      h_tt[i]->SetFillColor(kRed+1);
      hStk->Add( h_tt[i] );
      h_qcd[i]->SetFillColor(kYellow);
      hStk->Add( h_qcd[i] );

      // Jet1 Et Spectrum
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c1->SetLogy();
      if ( i==8 ) h_dt[i]->SetMinimum(10.);
      c1->cd();

      //if ( i > 4 ) h_dt[i]->SetMaximum( yMax[i-5] );
      h_dt[i]->SetMarkerSize(1);
      h_dt[i]->SetMarkerStyle(21);
      h_dt[i]->Draw("PE");
      
      if ( i== 7 ) {
         //double n_W = MtFitter( h_dt[i], *fjb1, c1, 50, 120  ) ;
         cout<<" orignal QCD MC = "<< h_qcd[i]->Integral() <<endl;
      }
      
      c1->Update();

      hStk->Draw("same");
      c1->Update();
      h_dt[i]->Draw("PE SAME");
      //if ( i == 7 ) fjb1->Draw("sames") ;
      c1->Update();

      TString plotname1 = hfolder + theSubFolder+ DataName + hNames[i]+ "."+plotType ;
      c1->Print( plotname1 );
      delete c1;
      delete hStk ;
  }
 
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;
  //delete fjb1;

}

void ObjectInfo::CombinedMCPlotter( vector<string>& fakeData, bool doWMtFit, bool doScale ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "CombinedMC" );
  gSystem->cd( "../" );

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  Reset( 1, doScale ) ;

  hObjs* htt = new hObjs() ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj") ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd") ;
  EvtSelector( fakeData[2], hqcd, false, scale2 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  hObjs* hzj = new hObjs("zj") ;
  EvtSelector( fakeData[3], hzj, false, scale3 );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq") ;
  EvtSelector( fakeData[4], htq, false, scale4 );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw") ;
  EvtSelector( fakeData[5], htw, false, scale5 );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  TString hNames[18] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",    "Mt_4" } ;

  TF1* fjb1 = new TF1("fjb1", MassFitFunction::ConvJacob, 0, 140, 7);
  //gStyle->SetOptStat("ieruom");
  //gStyle->SetOptFit(111);
  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize( 0.06, "X");
  gStyle->SetLabelSize( 0.06, "Y");
  gStyle->SetLabelSize( 0.06, "Z");

  for (int i=0; i<18; i++) {
      if ( i>0 && i <4 ) continue;
      //TH1D* All_MC  = new TH1D("All_MC", "Combined MC ", 15, 0, 150);
      TH1D* All_MC  = (TH1D*) h_wj[i]->Clone("combined MC"); 
      THStack* hStk = new THStack("hStk", hNames[i] );
      h_qcd[i]->SetFillColor(kYellow);
      hStk->Add( h_qcd[i] );
      h_zj[i]->SetFillColor(kAzure-2);
      hStk->Add( h_zj[i] );
      h_wj[i]->SetFillColor(kGreen);
      hStk->Add( h_wj[i] );
      h_tq[i]->SetFillColor(kMagenta+2);
      hStk->Add( h_tq[i] );
      h_tw[i]->SetFillColor(kMagenta);
      hStk->Add( h_tw[i] );
      h_tt[i]->SetFillColor(kRed+1);
      hStk->Add( h_tt[i] );

      //All_MC->Add( h_wj[i] );
      All_MC->Add( h_zj[i] );
      All_MC->Add( h_tt[i] );
      All_MC->Add( h_tw[i] );
      All_MC->Add( h_tq[i] );
      All_MC->Add( h_qcd[i] );

      if ( i == 7 ) {
         cout<<" => ";
         for( int k=1; k<8; k++ ) {
            cout<<" bin"<<k<<" : "<< h_wj[i]->GetBinContent(k);
         }
         cout<<" "<<endl;
      }

      // Jet1 Et Spectrum
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      //c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c1->SetLogy();
      c1->cd();

      hStk->Draw();
      All_MC->Draw("sames");
      if ( i== 7 && doWMtFit )  MtFitter( All_MC, *fjb1, c1, 60, 120  ) ;
      //if ( i == 7 && doWMtFit ) fjb1->Draw("sames") ;
      c1->Update();

      TString plotname1 = hfolder + "CombinedMC/MC_"+ hNames[i]+ "."+plotType ;
      c1->Print( plotname1 );
      delete c1;
      delete All_MC ;
      delete hStk ;
  }
 
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;
  delete fjb1;
}

void ObjectInfo::MCPlotter1( vector<string>& fakeData, double norm  ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hMCAll" );
  gSystem->cd( "../" );

  double scale0 = ( norm == 0 ) ? fitInput->NormalizeComponents( "tt" ) : norm / fitInput->TreeSize( fakeData[0] ); 
  double scale1 = ( norm == 0 ) ? fitInput->NormalizeComponents( "wj" ) : norm / fitInput->TreeSize( fakeData[1] );
  double scale2 = ( norm == 0 ) ? fitInput->NormalizeComponents( "qcd" ): norm / fitInput->TreeSize( fakeData[2] );
  double scale3 = ( norm == 0 ) ? fitInput->NormalizeComponents( "zj" ) : norm / fitInput->TreeSize( fakeData[3] );
  //double scale4 = ( norm == 0 ) ? fitInput->NormalizeComponents( "tq" ) : norm / fitInput->TreeSize( fakeData[4] );
  //double scale5 = ( norm == 0 ) ? fitInput->NormalizeComponents( "tw" ) : norm / fitInput->TreeSize( fakeData[5] );

  /*
  double scale0 = fitInput->NormalizeComponents( "tt" )  ; 
  double scale1 = fitInput->NormalizeComponents( "wj" )  ;
  double scale2 = fitInput->NormalizeComponents( "qcd" ) ;
  double scale3 = fitInput->NormalizeComponents( "zj" )  ;
  */

  hObjs* htt = new hObjs() ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs( "wj" ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hqcd = new hObjs( "qcd" ) ;
  EvtSelector( fakeData[2], hqcd, false, scale2 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  hObjs* hzj = new hObjs( "zj" ) ;
  EvtSelector( fakeData[3], hzj, false, scale3 );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  /* 
  hObjs* htq = new hObjs( "tq" ) ;
  EvtSelector( fakeData[4], htq, false, scale4 );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs( "tw" ) ;
  EvtSelector( fakeData[5], htw, false, scale5 );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );
  */

  TString hNames[18] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",    "Mt_4" } ;

  double yMax[18]    = { 1000., 1000., 1000., 1000., 1000.,
                                 500.,  100., 1500., 1000.,
                                 300.,  150.,  150,  500., 
                                 500.,  500.,  500., 800.,  800. };

  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize( 0.06, "X");
  gStyle->SetLabelSize( 0.06, "Y");
  gStyle->SetLabelSize( 0.06, "Z");

  for (int i=0; i<18; i++) {
      if ( i>0 && i <4 ) continue;

      h_wj[i]->SetLineWidth(3);
      h_zj[i]->SetLineWidth(3);
      h_qcd[i]->SetLineWidth(3);
      h_tt[i]->SetLineWidth(3);
      //h_tq[i]->SetLineWidth(2);
      //h_tw[i]->SetLineWidth(2);

      h_wj[i]->SetLineColor(kGreen);
      h_zj[i]->SetLineColor(kAzure-2);
      h_tt[i]->SetLineColor(kRed+1);
      h_qcd[i]->SetLineColor(kBlack);
      //h_tq[i]->SetLineColor(kMagenta+2);
      //h_tw[i]->SetLineColor(kMagenta);

      
      if ( norm > 0 ) {   
         double nwj = h_wj[i]->Integral();
	 h_wj[i]->Scale( norm /nwj );
	 double nzj = h_zj[i]->Integral();
	 h_zj[i]->Scale( norm /nzj );
	 double nqcd = h_qcd[i]->Integral();
	 h_qcd[i]->Scale( norm /nqcd );
	 double ntt = h_tt[i]->Integral();
	 h_tt[i]->Scale( norm /ntt );
	 cout<<" scale wj = "<< norm /nwj <<" zj = "<< norm/nzj << " scale qcd = "<< norm/nqcd <<" scale tt = "<< norm/ntt << endl; 
      }

      // making plots
      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      //c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      if ( i== 0 || i == 7 || i== 8 ) c1->SetLogy();
      c1->cd();

      yMax[i] = ( norm > 0. ) ? yMax[i] : h_wj[i]->Integral()/2 ;
      h_wj[i]->SetMaximum( yMax[i] );

      h_wj[i]->Draw();
      h_zj[i]->Draw("same");
      //h_tq[i]->Draw("same");
      //h_tw[i]->Draw("same");
      h_tt[i]->Draw("same");
      h_qcd[i]->Draw("same");

      TLegend *leg = new TLegend(.75, .65, .90, .90 );
      leg->AddEntry(h_tt[i],  "Ttbar",  "L");
      leg->AddEntry(h_wj[i],  "W+Jets", "L");
      leg->AddEntry(h_zj[i],  "Z+Jets", "L");
      leg->AddEntry(h_qcd[i], "QCD",    "L");
      leg->Draw("same");

      c1->Update();

      TString plotname1 = hfolder + "hMCAll/MC_"+ hNames[i]+ "."+plotType ;
      c1->Print( plotname1 );
      delete c1;
      delete leg;
  }
 
  delete htt;
  //delete htq;
  //delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;
 

}

vector<double> ObjectInfo::BinErr( double m ){

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
      }
     // +34%
     double j = m ;
     double hm = 0 ;
     double hp = 0 ;
     while ( hm <=0.34 || j < 0 ) {
           j = j + step ;
           hp = TMath::Poisson( j, m );
           hm = hm + (hp*step) ;
     }
     pErr.push_back( k - m );
     pErr.push_back( j - m );
  }
  return pErr ;

}

void ObjectInfo::Reset( int idx, bool newValue ){
  if ( idx != -1 ){
     if ( idx == 0 ) inclu = newValue ;
     if ( idx == 1 ) isWeight = newValue ;
  } else {
     if ( idx == 0 ) fitInput->GetParameters( "Inclusive", &Inclusive );
     if ( idx == 1 ) isWeight = false ;
  }
}

void ObjectInfo::Reset( int idx, int newValue ){
  if ( newValue != -1 ){
     if ( idx == 0 ) n_Jets =  newValue ;
  } else {
     if ( idx == 0 ) fitInput->GetParameters( "n_Jets",   &n_Jets );
  }
}
void ObjectInfo::Reset( int idx, double newValue ){
  if ( newValue != -1 ){
     if ( idx == 0 ) jetCuts[0] = newValue ; 
     if ( idx == 1 ) jetCuts[1] = newValue ; 
     if ( idx == 2 ) muonCuts[0] = newValue ; 
     if ( idx == 3 ) muonCuts[1] = newValue ; 
     if ( idx == 4 ) muonCuts[2] = newValue ; 
  } else {
     if ( idx == 0 || idx == 1 ) {
        fitInput->GetParameters( "JetCuts",  &jetCuts );
     }
     if ( idx == 2 || idx == 3 || idx == 4 ) {
        fitInput->GetParameters( "MuonCuts", &muonCuts );
     }
  }
}

void ObjectInfo::ResetBGCuts( int idx, double newValue ){

   if ( idx == 0 ) bgMtCut  = newValue ;
   if ( idx == 1 ) bgMETCut = newValue ;

   if ( idx == - 1 ) {
      bgMtCut  = 35 ;
      bgMETCut = 20 ;
   }

}

double ObjectInfo::EvtScaling( double w_pt, string fileName ){

    // Mu18, 0J
    //double qSet[5] = { 1.72, 1.85, 1.48, 1.26, 1.03 };
    //double vSet[5] = { 1.06, 1.00, 1.46, 1.39, 1.44 };
    // PF
    // double qSet[5] = { 1.45, 2.15, 1.20, 1.40, 5.91  }; // with V+Jets scaling
    //double qSet[5] = { 1.46, 2.16, 1.20, 1.43, 6.2  };
    //double vSet[5] = { 1.02, 1.03, 0.82, 1.16, 1.36 };
    // Calo
    //double qSet[5] = { 1.69, 1.78, 1.57, 0.97, 42.65  }; // with V+Jets scaling
    double qSet[5] = { 1.70, 1.79, 1.57, 0.96, 24.41  };
    double vSet[5] = { 1.02, 0.98, 1.40, 1.33, 1.10 };

    
    double theScale = 1 ;
    if ( fileName.substr(0,2) == "wj" || fileName.substr(0,2) == "zj" || fileName.substr(0,1) == "t" ) {
       if ( w_pt < 20 )               theScale =  vSet[1] ;
       if ( w_pt >= 20 && w_pt < 40 ) theScale =  vSet[2] ;
       if ( w_pt >= 40 && w_pt < 60 ) theScale =  vSet[3] ;
       if ( w_pt >= 60 )              theScale =  vSet[4] ;
    }
    if ( fileName.substr(0,2) == "qc" ) {
       if ( w_pt < 20 )               theScale =  qSet[1] ;
       if ( w_pt >= 20 && w_pt < 40 ) theScale =  qSet[2] ;
       if ( w_pt >= 40 && w_pt < 60 ) theScale =  qSet[3] ;
       if ( w_pt >= 60 )              theScale =  qSet[4] ;
    }

    return theScale ;
}

double ObjectInfo::EvtScaling( int NJets, string fileName ){

    //Calo Mu18,         1J    2J+
    //double qSet[2] =  { 1.43, 1.00 } ;
    //double vSet[2] =  { 1.02, 1.32 } ;
    //Calo Mu18,         1J    2J+   , Template fitting
    //double qSet[2] =  { 1.41, 1.17 } ;
    //double vSet[2] =  { 1.10, 1.46 } ;
    //PF   Mu18,         1J    2J+   , Template fitting
    //double qSet[2] =  { 1.41, 1.47 } ;
    //double vSet[2] =  { 1.08, 1.35 } ;

    vector<double> qSet;
    fitInput->GetParameters( "qcdNorm", &qSet );
    vector<double> vSet;
    fitInput->GetParameters( "vjNorm", &vSet );

    double theScale = 1 ;
    if ( fileName.substr(0,2) == "wj" || fileName.substr(0,2) == "zj" ) {
       if ( NJets == 1 ) theScale =  vSet[0] ;
       if ( NJets >= 2 ) theScale =  vSet[1] ;
       //if ( NJets == 3 ) theScale =  (vSet[1]*vSet[1])/vSet[0] ;
       //if ( NJets >= 4 ) theScale =  (vSet[1]*vSet[1]*vSet[1])/(vSet[0]*vSet[0]) ;
    }
    if ( fileName.substr(0,2) == "qc" ) {
       if ( NJets == 1 ) theScale =  qSet[0] ;
       if ( NJets >= 2 ) theScale =  qSet[1] ;
       //if ( NJets == 3 ) theScale =  (qSet[1]*qSet[1])/qSet[0] ;
       //if ( NJets >= 4 ) theScale =  (qSet[1]*qSet[1]*qSet[1])/(qSet[0]*qSet[0]) ;
    }

    return theScale ;
}

double ObjectInfo::QCDScaling( double w_pt ){

    // Mu18, 0J inclusive
    //double qSet[4] = { 1.88, 1.55, 1.11, 12.17 };
    // Mu18, 1J inclusive
    double qSet[4] = { 1.32, 1.33, 1.05, 13.1 };

    double qScale = 1 ;
    if ( w_pt < 20 )               qScale = qSet[0];
    if ( w_pt >= 20 && w_pt < 40 ) qScale = qSet[1];
    if ( w_pt >= 40 && w_pt < 60 ) qScale = qSet[2];
    if ( w_pt >= 60 )              qScale = qSet[3];

    return qScale ;
}

double ObjectInfo::QCDScaling( int NJets, string fileName ){
    // calo
    //double vSet[2] =  { 1.04, 1.38 } ;
    // pf
    // double vSet[2] =  { 1.07, 1.24 } ; Mu18, 1.93 /nb
    // double vSet[2] =  { 0.95, 1.16 } ; // Mu20, 2.88 /nb
    // double vSet[2] =  { 1.17, 1.39 } ; // Mu20, 2.88 /nb , JES + 5%
    vector<double> vSet;
    fitInput->GetParameters( "vjNorm", &vSet );

    double theScale = 1 ;
    if ( fileName.substr(0,2) == "wj" || fileName.substr(0,2) == "zj" ) {
       if ( NJets == 1 ) theScale =  vSet[0] ;
       if ( NJets >= 2 ) theScale =  vSet[1] ;
       //if ( NJets == 3 ) theScale =  (vSet[1]*vSet[1])/vSet[0] ;
       //if ( NJets >= 4 ) theScale =  (vSet[1]*vSet[1]*vSet[1])/(vSet[0]*vSet[0]) ;
    }
    return theScale ;
}

void ObjectInfo::QCDBGPlotter( string DataName, vector<string>& fakeData, bool doQCD, double MtCut, double METCut, double scale ){

  ResetBGCuts(0, MtCut );
  ResetBGCuts(1, METCut );

  TString theFolder = hfolder ;
  TString theSubFolder = "hQCDMu18/" ;
  if ( !doQCD ) theSubFolder = "hSGMu18/" ;
  if ( scale > 1 ) theSubFolder = "hSGMu18_Scale1/" ;
  if ( scale == -1 ) theSubFolder = "hSGMu18_Scale4/" ;

  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 30 ;

  hObjs* hdt = new hObjs("data", nbin ) ;
  QCDSelector( DataName, hdt, false, 1, doQCD );
  vector<TH1D*> h_dt ;
  hdt->Fill1DVec( h_dt );

  if ( h_dt[13]->Integral() < 500. ) {
     nbin = 15 ;
     h_dt[13]->Rebin(2);
     h_dt[14]->Rebin(2);
     h_dt[15]->Rebin(2);
  }
  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = scale * fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  hObjs* htt = new hObjs("tt", nbin ) ;
  QCDSelector( fakeData[0], htt, false, scale0, doQCD );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  QCDSelector( fakeData[1], hwj, false, scale1, doQCD );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  QCDSelector( fakeData[2], hqcd, false, scale2, doQCD );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  QCDSelector( fakeData[3], hzj, false, scale3, doQCD );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  QCDSelector( fakeData[4], htq, false, scale4, doQCD );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  QCDSelector( fakeData[5], htw, false, scale5, doQCD );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  //gStyle->SetOptStat("ieruom");
  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");

  TString hNames[20] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",    "Mt_4", "dRmin", "PtRel" } ;

  //for PF
  /*
  double   yMax[13] = {    900,    300,   4000,   4000,
                           900,    450,    600,   2500, 
                           550,    450,    150,     40,   30 } ;
   */
  //double yMax[13] = { 400,  70, 1000, 500, 
  //                    200,  90,  150, 500,
  //                    200, 180,  150,  60, 30  } ;

  for (int i=0; i<20; i++) {
      if ( i>0 && i <4 ) continue;
      THStack* hStk = new THStack("hStk", hNames[i] );
      h_wj[i]->SetFillColor(kGreen);
      hStk->Add( h_wj[i] );
      h_zj[i]->SetFillColor(kAzure-2);
      hStk->Add( h_zj[i] );
      h_tq[i]->SetFillColor(kMagenta+2);
      hStk->Add( h_tq[i] );
      h_tw[i]->SetFillColor(kMagenta);
      hStk->Add( h_tw[i] );
      h_tt[i]->SetFillColor(kRed+1);
      hStk->Add( h_tt[i] );
      h_qcd[i]->SetFillColor(kYellow);
      hStk->Add( h_qcd[i] );
 
      double wQCD[5] = {1,1,1,1,1} ;
      if ( i > 12 ) {
         double N_oMC = h_wj[i]->Integral() + h_zj[i]->Integral() + h_tt[i]->Integral() 
                        + h_tq[i]->Integral() + h_tw[i]->Integral() ;
         double N_data  = h_dt[i]->Integral() ;

         wQCD[i-13] = ( N_data - N_oMC ) / h_qcd[i]->Integral() ;
         cout<<" weighting "<< i-13 <<" = "<< wQCD[i-13]<<endl;
      }

      TCanvas* c1 = new TCanvas("c1","", 800, 600);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c1->SetLogy();
      if ( i==8 ) h_dt[i]->SetMinimum(1.);
      c1->cd();

      //if ( i > 4 ) h_dt[i]->SetMaximum( yMax[i-5] );
      
      h_dt[i]->SetMarkerSize(1);
      h_dt[i]->SetMarkerStyle(21);
      h_dt[i]->Draw("PE");
      c1->Update();

      hStk->Draw("same");
      c1->Update();
      h_dt[i]->Draw("PE SAME");
      c1->Update();

      TString plotname1 = hfolder + theSubFolder + DataName + hNames[i]+ "."+plotType ;
      c1->Print( plotname1 );
    
      delete c1;
      delete hStk ;
  }
 
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;

}

void ObjectInfo::QCDBG2DPlotter( string DataName, vector<string>& fakeData, bool doQCD, double MtCut, double METCut, double scale ){

  ResetBGCuts(0, MtCut );
  ResetBGCuts(1, METCut );

  TString theFolder = hfolder ;
  TString theSubFolder = "hQCDMu18/" ;
  if ( !doQCD ) theSubFolder = "hSGMu18/" ;
  if ( scale > 1 ) theSubFolder = "hSGMu18_Scale1/" ;
  if ( scale == -1 ) theSubFolder = "hSGMu18_Scale4/" ;

  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 15 ;

  hObjs* hdt = new hObjs("data", nbin ) ;
  QCDSelector( DataName, hdt, false, 1, doQCD );

  vector<TH2D*> h_dt ;
  hdt->Fill2DVec( h_dt );

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = scale * fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  hObjs* htt = new hObjs("tt", nbin ) ;
  QCDSelector( fakeData[0], htt, false, scale0, doQCD );
  vector<TH2D*> h_tt ;
  htt->Fill2DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  QCDSelector( fakeData[1], hwj, false, scale1, doQCD );
  vector<TH2D*> h_wj ;
  hwj->Fill2DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  QCDSelector( fakeData[2], hqcd, false, scale2, doQCD );
  vector<TH2D*> h_qcd ;
  hqcd->Fill2DVec( h_qcd );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  QCDSelector( fakeData[3], hzj, false, scale3, doQCD );
  vector<TH2D*> h_zj ;
  hzj->Fill2DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  QCDSelector( fakeData[4], htq, false, scale4, doQCD );
  vector<TH2D*> h_tq ;
  htq->Fill2DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  QCDSelector( fakeData[5], htw, false, scale5, doQCD );
  vector<TH2D*> h_tw ;
  htw->Fill2DVec( h_tw );

  gStyle->SetOptStat("ieruom");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");

  TH2D* hMetMt0 = new TH2D("hMetMt0", " MET(X) vs lep Mt(Y) ", 30, 0, 150,  nbin, 0, 150 );
  hMetMt0->Add( h_tt[0] );
  hMetMt0->Add( h_tq[0] );
  hMetMt0->Add( h_tw[0] );
  hMetMt0->Add( h_wj[0] );
  hMetMt0->Add( h_zj[0] );
  hMetMt0->Add( h_qcd[0] );

  TH2D* hWPtMt0 = new TH2D("hWptMt0", " W Pt(X) vs lep Mt(Y) ", 20, 0, 200,  nbin, 0, 150 );
  hWPtMt0->Add( h_tt[1] );
  hWPtMt0->Add( h_tq[1] );
  hWPtMt0->Add( h_tw[1] );
  hWPtMt0->Add( h_wj[1] );
  hWPtMt0->Add( h_zj[1] );
  hWPtMt0->Add( h_qcd[1] );

  TH2D* hWPtdF0 = new TH2D("hWptdF0", " W Pt(X) vs dPhi(mu, MET)(Y) ", 20, 0, 200,  32, 0, 3.2 );
  hWPtdF0->Add( h_tt[2] );
  hWPtdF0->Add( h_tq[2] );
  hWPtdF0->Add( h_tw[2] );
  hWPtdF0->Add( h_wj[2] );
  hWPtdF0->Add( h_zj[2] );
  hWPtdF0->Add( h_qcd[2] );

  TString hNames[3] = { "MetMt", "WPtMt", "WPtdF" } ;
  gStyle->SetStatX(0.85);
  gStyle->SetPalette(1);
  gStyle->SetOptStat("nieuo");
  for (int i=0; i<3; i++) {

      TCanvas* c1 = new TCanvas("c1","", 800, 700);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->Divide(1,2);

      c1->cd(1);
      gStyle->SetNumberContours(8);
      h_dt[i]->SetName("Data");
      h_dt[i]->Draw("COLZ");
      c1->Update();

      c1->cd(2);
      gStyle->SetNumberContours(8);
      if ( i ==0 )  hMetMt0->Draw("COLZ");
      if ( i ==1 )  hWPtMt0->Draw("COLZ");
      if ( i ==2 )  hWPtdF0->Draw("COLZ");
      c1->Update();

      TString plotname1 = hfolder + theSubFolder + hNames[i] + "."+plotType ;
      c1->Print( plotname1 );
    
      delete c1;
  }
 
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;
 
  delete hMetMt0;
  delete hWPtMt0;
  delete hWPtdF0;
 
}

void ObjectInfo::WScaleStudy( vector<string>& fakeData  ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "CombinedMC" );
  gSystem->cd( "../" );

  // for the damn systematic studies
  vector<string> SysSample ;
  fitInput->GetParameters( "OtherSamples", &SysSample );
  double uScaleW  = (3.1*28049*0.2177) /  348887;
  double dScaleW  = (3.1*28049*0.2177) / 1234550;
  double nScaleW = fitInput->NormalizeComponents( "wj" );

  Reset( 1, false ) ;

  hObjs* hwj = new hObjs("wj") ;
  EvtSelector( fakeData[1], hwj, false, nScaleW );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hwju = new hObjs("wj_up") ;
  EvtSelector( SysSample[0], hwju, false, uScaleW );
  vector<TH1D*> h_wju ;
  hwju->Fill1DVec( h_wju );

  hObjs* hwjd = new hObjs("wj_down") ;
  EvtSelector( SysSample[1], hwjd, false, dScaleW );
  vector<TH1D*> h_wjd ;
  hwjd->Fill1DVec( h_wjd );


  //gStyle->SetOptStat("ieruom");
  gStyle->SetOptStat("");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetErrorX(0) ;

  double Rn[5] = {0} ;
  double Ru[5] = {0} ;
  double Rd[5] = {0} ;
  double eRn[5] = {0} ;
  double eRu[5] = {0} ;
  double eRd[5] = {0} ;
  for(int i=0; i<5; i++){
     int j = i+1 ;
     Rn[i] =  h_wj[7]->GetBinContent(j+1) /  h_wj[7]->GetBinContent(j)  ;
     Ru[i] = h_wju[7]->GetBinContent(j+1) / h_wju[7]->GetBinContent(j)  ;
     Rd[i] = h_wjd[7]->GetBinContent(j+1) / h_wjd[7]->GetBinContent(j)  ;
     eRn[i] = MassFitFunction::ErrAovB(h_wj[7]->GetBinContent(j+1)/nScaleW,   h_wj[7]->GetBinContent(j)/nScaleW ) ;
     eRu[i] = MassFitFunction::ErrAovB(h_wju[7]->GetBinContent(j+1)/uScaleW, h_wju[7]->GetBinContent(j)/uScaleW ) ;
     eRd[i] = MassFitFunction::ErrAovB(h_wjd[7]->GetBinContent(j+1)/dScaleW, h_wjd[7]->GetBinContent(j)/dScaleW ) ;
  }

  /*
  for( int k=1; k<8; k++ ) {
      cout<<" bin"<<k<<" : "<< h_wj[7]->GetBinContent(k);
  }
  cout<<"  Int = "<< h_wj[7]->Integral() <<endl;
  for( int k=1; k<8; k++ ) {
      cout<<" bin"<<k<<" : "<< h_wju[7]->GetBinContent(k);
  }
  cout<<"  Int = "<< h_wju[7]->Integral() <<endl;
  for( int k=1; k<8; k++ ) {
      cout<<" bin"<<k<<" : "<< h_wjd[7]->GetBinContent(k);
  }
  cout<<"  Int = "<< h_wjd[7]->Integral() <<endl;
  */

  // Jet1 Et Spectrum
  TCanvas* c1 = new TCanvas("c1","", 850, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetLogy();
  c1->cd();


  h_wj[7]->SetMarkerColor(1);
  h_wj[7]->SetMarkerSize(1.5);
  h_wj[7]->SetMarkerStyle(20);

  h_wjd[7]->SetMarkerColor(4);
  h_wjd[7]->SetMarkerSize(1.5);
  h_wjd[7]->SetMarkerStyle(23);

  h_wju[7]->SetMarkerColor(2);
  h_wju[7]->SetMarkerSize(1.5);
  h_wju[7]->SetMarkerStyle(22);

  h_wj[7]->Draw("PE1");
  c1->Update();
  h_wju[7]->Draw("samePE1");
  h_wjd[7]->Draw("samePE1");

  TLegend *leg = new TLegend(.22, .16, .45, .35 );
  leg->AddEntry(h_wj[7],  "Nominal",    "P");
  leg->AddEntry(h_wju[7], "Scale-Up",   "P");
  leg->AddEntry(h_wjd[7], "Scale-Down", "P");
  leg->Draw("same");

  c1->Update();
  TString plotname1 = hfolder + "CombinedMC/MC_WJ."+plotType ;
  c1->Print( plotname1 );

  TH1D* hRn = new TH1D("hRn"," VJets ratio1 ", 5, 0.5, 5.5 );
  TH1D* hRu = new TH1D("hRu"," VJets ratio2 ", 5, 0.5, 5.5 );
  TH1D* hRd = new TH1D("hRd"," VJets ratio3 ", 5, 0.5, 5.5 );

  for ( int i=0; i<5; i++){
	  hRn->Fill( i+1. , Rn[i] );
	  hRu->Fill( i+1. , Ru[i] );
	  hRd->Fill( i+1. , Rd[i] );
	  hRn->SetBinError( i+1, eRn[i] );
	  hRu->SetBinError( i+1, eRu[i] );
	  hRd->SetBinError( i+1, eRd[i] );
  }

  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize( 0.04, "X");
  gStyle->SetLabelSize( 0.04, "Y");
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetErrorX(0) ;

  TCanvas* c2 = new TCanvas("c2","", 850, 700);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->cd();

  hRn->SetAxisRange(0.1, 0.26, "Y");
  hRn->SetXTitle(" Ratio between different jet multiplicity");
  hRn->SetYTitle(" #frac{Number of N jets}{Number of N-1 jets} " );
  hRn->SetTitleOffset(1.5, "X") ;
  hRn->SetTitleOffset(2.0, "Y") ;
  hRn->SetLabelSize(0.04, "Y");

  hRn->SetMarkerColor(1);
  hRn->SetMarkerSize(2);
  hRn->SetMarkerStyle(20);

  hRn->GetXaxis()->SetBinLabel(1, "R(1/0)");
  hRn->GetXaxis()->SetBinLabel(2, "R(2/1)");
  hRn->GetXaxis()->SetBinLabel(3, "R(3/2)");
  hRn->GetXaxis()->SetBinLabel(4, "R(4/3)");
  hRn->GetXaxis()->SetBinLabel(5, "R(5/4)");

  hRn->Draw("PE1");
  c2->Update();
  
  hRu->SetMarkerColor(2);
  hRu->SetMarkerSize(2);
  hRu->SetMarkerStyle(22);

  hRd->SetMarkerColor(4);
  hRd->SetMarkerSize(2);
  hRd->SetMarkerStyle(23);

  hRu->Draw("samePE1");
  hRd->Draw("samePE1");

  TLegend *leg1 = new TLegend(.75, .75, .99, .99 );
  leg1->AddEntry(hRn, "Nominal",    "P");
  leg1->AddEntry(hRu, "Scale-Up",   "P");
  leg1->AddEntry(hRd, "Scale-Down", "P");

  leg1->Draw("same");

  c2->Update();
  
  TString plotname2 = hfolder + "CombinedMC/RatioMC_WJ." + plotType ;
  c2->Print( plotname2 );


  delete hwj;
  delete hwju;
  delete hwjd;
  delete hRn;
  delete hRu;
  delete hRd;
  delete leg;
  delete leg1;
  delete c1;
  delete c2;
}


void ObjectInfo::Data2DPlotter( string DataName, vector<string>& fakeData, bool doScale ){

  TString theFolder = hfolder ;
  TString theSubFolder = ( doScale == false ) ?  "hData18/" : "hData18_Scale/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 30 ;

  hObjs* hdt = new hObjs("data", 30) ;
  EvtSelector( DataName, hdt );
  vector<TH2D*> h_dt ;
  hdt->Fill2DVec( h_dt );

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  Reset( 1, doScale ) ;

  hObjs* htt = new hObjs("tt", nbin ) ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH2D*> h_tt ;
  htt->Fill2DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH2D*> h_wj ;
  hwj->Fill2DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  EvtSelector( fakeData[2], hqcd, false, scale2 );
  vector<TH2D*> h_qcd ;
  hqcd->Fill2DVec( h_qcd );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  EvtSelector( fakeData[3], hzj, false, scale3 );
  vector<TH2D*> h_zj ;
  hzj->Fill2DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  EvtSelector( fakeData[4], htq, false, scale4 );
  vector<TH2D*> h_tq ;
  htq->Fill2DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  EvtSelector( fakeData[5], htw, false, scale5 );
  vector<TH2D*> h_tw ;
  htw->Fill2DVec( h_tw );

  //gStyle->SetOptStat("ieruom");
  //gStyle->SetStatX(0.85);
  gStyle->SetOptStat("");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(8);

  TString hNames[4] = { "MET_Mt", "WPt_Mt", "WPt_dPhi",  "MuMET" } ;

  for (int i=0; i<4; i++) {

      TH2D* h_all =  (TH2D*) h_tw[i]->Clone("combinedMC");
      h_all->Add( h_qcd[i] );
      h_all->Add( h_wj[i] );
      h_all->Add( h_zj[i] );
      h_all->Add( h_tt[i] );
      h_all->Add( h_tq[i] );
      h_all->Add( h_tw[i] );

      // Jet1 Et Spectrum
      TCanvas* c1 = new TCanvas("c1","", 800, 700);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->Divide(2,2);

      c1->cd(1);
      h_dt[i]->Draw("COLZ");
      c1->Update();

      c1->cd(2);
      h_wj[i]->Draw("COLZ");
      c1->Update();

      c1->cd(3);
      h_qcd[i]->Draw("COLZ");
      c1->Update();

      c1->cd(4);
      h_all->Draw("COLZ");
      c1->Update();

      TString plotname1 = hfolder + theSubFolder+ DataName + hNames[i]+ "."+plotType ;
      c1->Print( plotname1 );
      delete c1;
      delete h_all;
  }
 
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hzj;
  delete hqcd;

}

double  ObjectInfo::PtRel( TLorentzVector v1, TLorentzVector v2 ){

       double px = ( v1.Y()*v2.Z() ) - ( v1.Z()*v2.Y() );
       double py = ( v1.Z()*v2.X() ) - ( v1.X()*v2.Z() );
       double pz = ( v1.X()*v2.Y() ) - ( v1.Y()*v2.X() );

       double p = sqrt( (px*px) + (py*py) + (pz*pz) ) ;
       double p_v2 = sqrt( (v2.X()*v2.X())  + (v2.Y()*v2.Y()) + (v2.Z()*v2.Z()) );

       double ptRel =  p / p_v2  ;
       return ptRel ;

}

