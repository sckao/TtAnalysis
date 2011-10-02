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
  fitInput->GetParameters( "JESType", &JESType );
  inclu = ( Inclusive == "YES" ) ? true : false ;
  fitInput->GetParameters( "DataLike",   &dataLike );

  if (  inclu ) cout<<" Inclusive !!" <<endl;
  if ( !inclu ) cout<<" Exclusive !!" <<endl;
  isWeight = false ;
  
  bgMtCut  = 35;
  bgMETCut = 20;
  maxNJets = static_cast<int>( jetCuts[3] ) ;

}

ObjectInfo::~ObjectInfo(){

  delete fitInput ;
  delete pseudoExp;

}

// inclusive = true :  >=njet  
void ObjectInfo::EvtSelector( string fileName, recoObj* histos, bool smearing, double scale, vector<int>* evtlist, int evtSplit ) {

  TTree* tr = fitInput->TreeMap( fileName );
  bool isData = ( fileName.substr(0,4) == "data") ? true : false ;
  if ( dataLike == "True" ) isData = true ;

  // this is only for MC Matching
  double jpx[15],jpy[15],jpz[15],jE[15], jes[15], jer[15], bDis[15];
  int    n90[15] , nHits[2];
  double mpx[2],mpy[2],mpz[2],mE[2], mIso[2], d0[2], X2[2] ;
  double npx[3],npy[3],npz[3],nE[3] ;
  double metErrX, metErrY, WBRc ;
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
  if ( !isData  ) {
     tr->SetBranchAddress("jes"    ,&jes);
     tr->SetBranchAddress("jer"    ,&jer);
     tr->SetBranchAddress("metErrX", &metErrX);
     tr->SetBranchAddress("metErrY", &metErrY);
     tr->SetBranchAddress("WBRc",    &WBRc);
  } else {
     metErrX = 1. ;
     metErrY = 1. ;
     WBRc    = 1. ;
  }
  int idx = 0;
  int loopEnd = ( evtlist == NULL ) ? tr->GetEntries() : evtlist->size() ;

  //int fk = 0;
  //int ffk = 0;
  //int fffk = 0;
  for (int k=0; k< loopEnd ; k++) {
      double HTlep = 0 ;
      double HT = 0 ;
      if ( evtlist != NULL && evtlist->size() == 0 ) break;
      idx = ( evtlist == NULL ) ? k : (*evtlist)[k] ;

      if ( evtSplit == -1 && (k%2) == 0 ) continue ;
      if ( evtSplit ==  1 && (k%2) == 1 ) continue ;
      if ( evtSplit >  1 && (k%evtSplit ) != 0 ) continue ;
      if ( evtSplit < -1 && (k%evtSplit ) != 1 ) continue ;

      tr->GetEntry( idx );
      //if ( nJ < n_Jets ) continue ;
      
      TLorentzVector m0( mpx[0], mpy[0], mpz[0], mE[0] );
      if ( m0.Pt()        < muonCuts[0] ) continue ;
      if ( fabs(m0.Eta()) > muonCuts[1] ) continue ;
      if ( mIso[0]        > muonCuts[2] ) continue ;

      vector<TLorentzVector> objlist;
      int NCountedJets = 0;
      double _jes = 1. ;
      double metCorrVx = 0.;
      double metCorrVy = 0.;
      bool badMuon = false ;
      double mindRmj = 99 ;
      double jec_tot = 0 ;
      for ( int i = 0; i< nJ ; i ++) {
          // get jet p4, jes will be applied if the switch is turn one
          TLorentzVector jn( jpx[i], jpy[i], jpz[i], jE[i] );
          double jec_a   = (0.75*0.8*2.2) / jn.Pt()   ;
          double jec_b   = ( bDis[i] > 3 ) ? 0.03 : 0. ;
          if ( bDis[i] > 3 && jn.Pt() > 50. &&  jn.Pt() < 200. && fabs(jn.Eta()) < 2.  ) jec_b = 0.02 ;
          jec_tot = sqrt( (jes[i]*jes[i]) + (jec_a*jec_a) + (jec_b*jec_b) ) ;
          if ( !isData )                 _jes = jer[i] ;
          if ( JESType == 2 && !isData ) _jes = jer[i] + sqrt( (jes[i]*jes[i]) + (jec_a*jec_a) + (jec_b*jec_b) ) ;
          if ( JESType == 3 && !isData ) _jes = jer[i] - sqrt( (jes[i]*jes[i]) + (jec_a*jec_a) + (jec_b*jec_b) ) ;
          if ( JESType == 4 && !isData ) _jes = 1. ;
          if ( JESType == 5 && !isData ) _jes = 1 + 2*( jer[i] -1 ) ;
          if ( JESType == 8 && !isData ) _jes = 1 + jetCuts[4]*( jer[i] -1 ) ;

          metCorrVx += jn.Px() ;
          metCorrVy += jn.Py() ;
          jn = jn*_jes;
          metCorrVx -= jn.Px() ;
          metCorrVy -= jn.Py() ;

          // sync cuts
          if ( jn.Pt()          < jetCuts[0] ) continue ;
          if ( fabs( jn.Eta() ) > jetCuts[1] ) continue ;
          double dRmj = m0.DeltaR( jn );
          if ( dRmj < muonCuts[3] ) badMuon = true ;
          mindRmj = ( dRmj < mindRmj ) ? dRmj : mindRmj ;
          if ( badMuon ) break ;

          objlist.push_back( jn );
          NCountedJets++ ;
          HT += jn.Pt() ;
          HTlep += jn.Pt() ;
      }
      if ( badMuon ) continue ;

      if ( JESType == 6 && !isData ) metCorrVx += metErrX ;
      if ( JESType == 6 && !isData ) metCorrVy += metErrY ;
      if ( JESType == 7 && !isData ) metCorrVx -= metErrX ;
      if ( JESType == 7 && !isData ) metCorrVy -= metErrY ;

      //cout<<" nJets? = "<< NCountedJets <<" sz = "<< objlist.size()<<" nJ = "<< nJ <<endl;
      if ( NCountedJets != n_Jets && !inclu ) continue ;
      if ( NCountedJets > maxNJets || NCountedJets < n_Jets ) continue ;

      while ( objlist.size() < 4 ) {
            TLorentzVector j0( 0., 0., 0., 0. );
            objlist.push_back( j0 );
      }
      TLorentzVector n0( npx[0], npy[0], npz[0], nE[0] );
      double neuPx = npx[0] + metCorrVx ;
      double neuPy = npy[0] + metCorrVy ;
      if ( n0.Pt() < jetCuts[2] ) continue;

      if ( !isData ) n0 = TLorentzVector( neuPx , neuPy, 0., sqrt( (neuPx*neuPx) +  (neuPy*neuPy) ) );

      objlist.push_back( m0 );
      objlist.push_back( n0 );
      HTlep += m0.Pt() ;
      HT += m0.Pt() ;
      HT += n0.Pt() ;
      //if ( npz[0] == 0 ) fk++ ;
      //if ( npz[0] != 0 ) ffk++ ;
      //if ( nNu == 1 ) fffk++ ;
      /*
      double MtdPhi = m0.DeltaPhi( n0 );
      TLorentzVector MtP4 = m0 + n0 ;
      double dPhi_upBound = 2 - ( MtP4.Pt() / 50 ) ;
      if ( fabs(MtdPhi) > dPhi_upBound  ) continue ;
      */

      /*
      if ( smearing )  {
         //pseudoExp->PhaseSmearing( objlist, 0, jec_tot );
         pseudoExp->PhaseSmearing( objlist, 0 );
         pseudoExp->JetEtReSort( objlist );
      }
      */

      int sz = objlist.size() ;
      //cout<<" obj size = "<< sz <<endl ;
      histos->gethad( objlist[0], objlist[1], objlist[2] );   
      histos->getlep( objlist[3], objlist[sz-2], objlist[sz-1] );  
      double muonQv[7] = { mIso[0], d0[0], X2[0], mindRmj, 0, HT, HTlep };
      int intQv[2] = { NCountedJets, nHits[0] };
      histos->getFloats( muonQv );
      histos->getIntegrals( intQv );
      //histos->getOther( NCountedJets, nHits[0], mIso[0], d0[0], X2[0] ) ;
      //
      //TLorentzVector vM2 = objlist[sz-2] + objlist[sz-1] ;
      //double vPt = vM2.Pt() ;
      //double weight = ( isWeight ) ? EvtScaling(vPt, fileName ) : 1. ;
      double weight = ( isWeight ) ? EvtScaling( NCountedJets, fileName ) : 1. ;
      histos->Fillh( 1., scale*weight*WBRc ) ;
  }
  //cout<<" neu fk = "<<fk <<" ffk = "<< ffk <<"  -> "<<  fffk <<endl;

}

// mode = 1 : measuring MC normalization, mode = 2 : measuring the ratio
void ObjectInfo::QCDSelector( string fileName, recoObj* histos, bool smearing, double scale, bool doQCD, int mode, int evtSplit ) {

  TTree* tr = fitInput->TreeMap( fileName );
  bool isData = ( fileName.substr(0,4) == "data") ? true : false ;
  if ( dataLike == "True" ) isData = true ;

  // this is only for MC Matching
  double jpx[15],jpy[15],jpz[15],jE[15], jes[15], jer[15], bDis[15];
  int    n90[15] , nHits[2];
  double mpx[2],mpy[2],mpz[2],mE[2], mIso[2], d0[2], X2[2] ;
  double npx[3],npy[3],npz[3],nE[3] ;
  double metErrX, metErrY, WBRc ;
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
  if ( !isData ) { 
     tr->SetBranchAddress("jes"    ,&jes);
     tr->SetBranchAddress("jer"    ,&jer);
     tr->SetBranchAddress("metErrX", &metErrX);
     tr->SetBranchAddress("metErrY", &metErrY);
     tr->SetBranchAddress("WBRc",    &WBRc);
  } else {
     metErrX = 1. ;
     metErrY = 1. ;
     WBRc    = 1. ;
  }

  cout<<" Current Mt cut = "<< bgMtCut <<"  MET cut = "<< bgMETCut << endl; 
  double Ht = 0 ;
  for (int k=0; k< tr->GetEntries() ; k++) {
      Ht =0 ;
      if ( evtSplit == -1 && (k%2) == 0 ) continue ;
      if ( evtSplit ==  1 && (k%2) == 1 ) continue ;
      if ( evtSplit >  1 && (k%evtSplit ) != 0 ) continue ;
      if ( evtSplit < -1 && (k%evtSplit ) != 1 ) continue ;

      tr->GetEntry(k);
      //if ( nJ < n_Jets ) continue ;
      
      TLorentzVector m0( mpx[0], mpy[0], mpz[0], mE[0] );
      if ( m0.Pt()  < muonCuts[0] ) continue ;
      if ( fabs( m0.Eta()) > muonCuts[1] ) continue ;
      if ( mIso[0]  > muonCuts[2] ) continue ;

      vector<TLorentzVector> objlist;
      int NCountedJets = 0;
      double dRmin = 999 ;
      double ptRel = 999 ;
      double _jes = 1. ;
      double metCorrVx = 0.;
      double metCorrVy = 0.;
      bool   badMuon = false ;
      for ( int i = 0; i< nJ ; i ++) {
          TLorentzVector jn( jpx[i], jpy[i], jpz[i], jE[i] );
          double jec_a = (0.75*0.8*2.2) / jn.Pt()   ;
          double jec_b = ( bDis[i] > 3 ) ? 0.03 : 0. ;
          if ( bDis[i] > 3 && jn.Pt() > 50. &&  jn.Pt() < 200. && fabs(jn.Eta()) < 2.  ) jec_b = 0.02 ;
          if (  !isData                ) _jes = jer[i] ;
          if ( JESType == 2 && !isData ) _jes = jer[i] + sqrt( (jes[i]*jes[i]) + (jec_a*jec_a) + (jec_b*jec_b) ) ;
          if ( JESType == 3 && !isData ) _jes = jer[i] - sqrt( (jes[i]*jes[i]) + (jec_a*jec_a) + (jec_b*jec_b) ) ;
          if ( JESType == 4 && !isData ) _jes = 1 ;
          if ( JESType == 5 && !isData ) _jes = 1 + 2*( jer[i] -1 ) ;
          if ( JESType == 8 && !isData ) _jes = 1 + jetCuts[4]*( jer[i] -1 ) ;

          metCorrVx += jn.Px() ;
          metCorrVy += jn.Py() ;
          jn = jn*_jes;
          metCorrVx -= jn.Px() ;
          metCorrVy -= jn.Py() ;

          // sync cuts
          if ( jn.Pt()          < jetCuts[0] ) continue ;
          if ( fabs( jn.Eta() ) > jetCuts[1] ) continue ;
          double dRmj = m0.DeltaR( jn );
          if ( dRmj < muonCuts[3] ) badMuon = true ; 
          if ( badMuon ) break ;
          if ( dRmj < dRmin ) {
             dRmin = dRmj ;
             ptRel = PtRel( m0, jn );
          }
          Ht += jn.Pt() ;
          objlist.push_back( jn );
          NCountedJets++ ;
      }
      if ( badMuon ) continue ;

      if ( JESType == 6 ) metCorrVx +=  metErrX ;
      if ( JESType == 6 ) metCorrVy +=  metErrY ;
      if ( JESType == 7 ) metCorrVx -=  metErrX ;
      if ( JESType == 7 ) metCorrVy -=  metErrY ;

      //cout<<" nJets? = "<< NCountedJets <<" sz = "<< objlist.size()<<" nJ = "<< nJ <<endl;
      if ( NCountedJets != n_Jets && !inclu ) continue ;
      if ( NCountedJets > maxNJets || NCountedJets < n_Jets ) continue ;

      while ( objlist.size() < 4 ) {
            TLorentzVector j0( 0., 0., 0., 0. );
            objlist.push_back( j0 );
      }
      TLorentzVector n0( npx[0], npy[0], npz[0], nE[0] );
      double neuPx = npx[0]+metCorrVx ;
      double neuPy = npy[0]+metCorrVy ;
      if ( !isData ) n0 = TLorentzVector( neuPx , neuPy, 0., sqrt( (neuPx*neuPx) +  (neuPy*neuPy) ) );
      objlist.push_back( m0 );
      objlist.push_back( n0 );
      Ht += m0.Pt() ;
      Ht += n0.Pt() ;
      if ( n0.Pt() < jetCuts[2] ) continue;

      double dPhi = m0.DeltaPhi( n0 ) ;
      double theMt = sqrt( 2.*m0.Pt()*n0.Pt()*( 1. - cos(dPhi) ) );
      //if ( fabs(dPhi) > 2.5 && m0.Pt() < 25. && n0.Et() < 25. ) continue ;

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
      double muonQv[6] = { mIso[0], d0[0], X2[0], dRmin, ptRel, Ht };
      int intQv[2] = { NCountedJets, nHits[0] };
      histos->getFloats( muonQv );
      histos->getIntegrals( intQv );

      double weight = 1. ;
      // used for normalization study
      if ( isWeight && mode == 1 ) weight = QCDScaling( NCountedJets, fileName ) ;
      // used for Ratio study
      if ( isWeight && mode == 2 ) weight = EvtScaling( NCountedJets, fileName ) ;

      histos->Fillh( 1, fabs(scale)*weight*WBRc ) ;
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
  if ( fileName.substr(0,2) == "zj" )  scale = fitInput->NormalizeComponents( "zj" );
  if ( fileName.substr(0,2) == "tq" )  scale = fitInput->NormalizeComponents( "tq" );
  if ( fileName.substr(0,2) == "tw" )  scale = fitInput->NormalizeComponents( "tw" );
  if ( fileName.substr(0,2) == "ww" )  scale = fitInput->NormalizeComponents( "ww" );
  if ( fileName.substr(0,2) == "qc" )  scale = fitInput->NormalizeComponents( "qcd" );

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

  TCanvas* c5 = new TCanvas("c5","", 800, 700);
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

void ObjectInfo::JESPlotter( vector<string>& fakeData, double norm ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hJES" );
  gSystem->cd( "../" );

  // only use the scale up value for different jes_type 
  int jes_type = JESType ;

  double scale0 = ( norm == 0 ) ? fitInput->NormalizeComponents( "tt" ) : norm / fitInput->TreeSize( fakeData[0] );
  double scale1 = ( norm == 0 ) ? fitInput->NormalizeComponents( "wj" ) : norm / fitInput->TreeSize( fakeData[1] );
  double scale6 = ( norm == 0 ) ? fitInput->NormalizeComponents( "qcd" ): norm / fitInput->TreeSize( fakeData[6] );

  TString tName = "JES_" ;
  if ( jes_type == 5 ) tName = "JER_" ; 
  if ( jes_type == 6 ) tName = "MET_" ; 

  Reset(1, jes_type );   // scale up JES, Unclustered Energy, scale-down JER

  hObjs* httu = new hObjs("ttu" ) ;
  EvtSelector( fakeData[0], httu, false, scale0 );
  vector<TH1D*> h_ttu ;
  httu->Fill1DVec( h_ttu );

  hObjs* hwju = new hObjs( "wju" ) ;
  EvtSelector( fakeData[1], hwju, false, scale1 );
  vector<TH1D*> h_wju ;
  hwju->Fill1DVec( h_wju );

  hObjs* hqcdu = new hObjs( "qcdu" ) ;
  EvtSelector( fakeData[6], hqcdu, false, scale6 );
  vector<TH1D*> h_qcdu ;
  hqcdu->Fill1DVec( h_qcdu );

  Reset(1, 0 );   // no JES scale

  hObjs* htt = new hObjs("tt" ) ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs( "wj" ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hqcd = new hObjs( "qcd" ) ;
  EvtSelector( fakeData[6], hqcd, false, scale6 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  if ( jes_type == 2 || jes_type == 6 ) Reset(1, jes_type+1 );   // scale down JES
  if ( jes_type == 5                  ) Reset(1, jes_type-1 );   // scale down JER

  hObjs* httd = new hObjs("ttd" ) ;
  EvtSelector( fakeData[0], httd, false, scale0 );
  vector<TH1D*> h_ttd ;
  httd->Fill1DVec( h_ttd );

  hObjs* hwjd = new hObjs( "wjd" ) ;
  EvtSelector( fakeData[1], hwjd, false, scale1 );
  vector<TH1D*> h_wjd ;
  hwjd->Fill1DVec( h_wjd );

  hObjs* hqcdd = new hObjs( "qcdd" ) ;
  EvtSelector( fakeData[6], hqcdd, false, scale6 );
  vector<TH1D*> h_qcdd ;
  hqcdd->Fill1DVec( h_qcdd );

  gStyle->SetOptStat("");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");
  gStyle->SetOptTitle(0);

  TString hNames[21] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",  
                                   "Mt_4",     "dRmj",    "RelPt",  "Ht_lep" } ;


  for (int i=0; i<21; i++) {
      if ( i !=0 && i !=5 && i !=7  && i != 13 && i!= 15 && i != 21 ) continue;
      if ( i == 13 ) h_tt[i]->SetMinimum( 0.);
      if ( i == 13 ) h_wj[i]->SetMinimum( 0.);
      if ( i == 13 ) h_qcd[i]->SetMinimum( 0.);

      h_ttu[i]->SetLineWidth(3);
      h_ttu[i]->SetLineColor(kRed);
      h_tt[i]->SetLineColor(1);
      h_ttd[i]->SetLineWidth(3);
      h_ttd[i]->SetLineColor(kBlue);

      h_wju[i]->SetLineWidth(3);
      h_wju[i]->SetLineColor(kRed);
      h_wj[i]->SetLineColor(1);
      h_wjd[i]->SetLineWidth(3);
      h_wjd[i]->SetLineColor(kBlue);

      h_qcdu[i]->SetLineWidth(3);
      h_qcdu[i]->SetLineColor(kRed);
      h_qcd[i]->SetLineColor(1);
      h_qcdd[i]->SetLineWidth(3);
      h_qcdd[i]->SetLineColor(kBlue);


      // Tt jet multiplicity
      TCanvas* c0 = new TCanvas("c0","", 800, 700);
      c0->SetFillColor(10);
      c0->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c0->SetLogy();
      c0->cd();

      double tMax = 1.5* h_tt[i]->GetBinContent(h_tt[i]->GetMaximumBin()) ;
      h_tt[i]->SetMaximum( tMax );

      h_tt[i]->SetMarkerSize(1);
      h_tt[i]->SetMarkerStyle(21);
      h_tt[i]->Draw("P");
      c0->Update();
      h_ttu[i]->Draw("SAME");
      c0->Update();
      h_ttd[i]->Draw("SAME");
      c0->Update();

      TLegend *tleg = new TLegend(.62, .7, .9, .9 );
      tleg->AddEntry(h_ttu[i], "Scale-up",   "L");
      tleg->AddEntry(h_tt[i],  "Nominal",    "P");
      tleg->AddEntry(h_ttd[i], "Scale-down", "L");
      tleg->Draw("same");
      c0->Update();

      TString plotname0 = hfolder + "hJES/Tt_"+ tName + hNames[i]+ "."+plotType ;
      c0->Print( plotname0 );

      // WJ jet multiplicity
      TCanvas* c1 = new TCanvas("c1","", 800, 700);
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c1->SetLogy();
      c1->cd();

      double wMax = 1.5* h_wj[i]->GetBinContent(h_wj[i]->GetMaximumBin()) ;
      h_wj[i]->SetMaximum( wMax );

      h_wj[i]->SetMarkerSize(1);
      h_wj[i]->SetMarkerStyle(21);
      h_wj[i]->Draw("P");
      c1->Update();
      h_wju[i]->Draw("SAME");
      c1->Update();
      h_wjd[i]->Draw("SAME");
      c1->Update();

      TLegend *wleg = new TLegend(.62, .7, .9, .9 );
      wleg->AddEntry(h_wju[i], "Scale-up",   "L");
      wleg->AddEntry(h_wj[i],  "Nominal",    "P");
      wleg->AddEntry(h_wjd[i], "Scale-down", "L");
      wleg->Draw("same");
      c1->Update();

      TString plotname1 = hfolder + "hJES/WJ_"+ tName + hNames[i]+ "."+plotType ;
      c1->Print( plotname1 );

      // QCD jet multiplicity
      TCanvas* c2 = new TCanvas("c2","", 800, 700);
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c2->SetLogy();
      c2->cd();

      double qMax = 1.5* h_qcd[i]->GetBinContent(h_qcd[i]->GetMaximumBin()) ;
      h_qcd[i]->SetMaximum( qMax );

      h_qcd[i]->SetMarkerSize(1);
      h_qcd[i]->SetMarkerStyle(21);
      h_qcd[i]->Draw("P");
      c2->Update();
      h_qcdu[i]->Draw("SAME");
      c2->Update();
      h_qcdd[i]->Draw("SAME");
      c2->Update();

      TLegend *qleg = new TLegend(.62, .7, .9, .9 );
      qleg->AddEntry(h_qcdu[i], "Scale-up",   "L");
      qleg->AddEntry(h_qcd[i],  "Nominal",    "P");
      qleg->AddEntry(h_qcdd[i], "Scale-down", "L");
      qleg->Draw("same");
      c2->Update();

      TString plotname2 = hfolder + "hJES/QCD_" + tName + hNames[i]+ "."+plotType ;
      c2->Print( plotname2 );

      delete c0;
      delete c1;
      delete c2;
  }

}

void ObjectInfo::DataPlotter( string DataName, vector<string>& fakeData, bool doScale ){

  TString theFolder = hfolder ;
  TString theSubFolder = ( doScale == false ) ?  "hData/" : "hData_Scale/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 2 ; // scale of bin number

  hObjs* hdt = new hObjs("data", nbin) ;
  EvtSelector( DataName, hdt );
  vector<TH1D*> h_dt ;
  hdt->Fill1DVec( h_dt );

  /*
  if ( h_dt[13]->Integral() < 600. ) {
     nbin = 1 ;
     h_dt[13]->Rebin(2);
     h_dt[14]->Rebin(2);
     h_dt[15]->Rebin(2);
     h_dt[16]->Rebin(2);
  }
  */

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  Reset( 1, doScale ) ;

  hObjs* htt = new hObjs("tt", nbin ) ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  EvtSelector( fakeData[2], hzj, false, scale2 );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  EvtSelector( fakeData[3], htq, false, scale3 );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  EvtSelector( fakeData[4], htw, false, scale4 );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  hObjs* hww = new hObjs("ww", nbin ) ;
  EvtSelector( fakeData[5], hww, false, scale5 );
  vector<TH1D*> h_ww ;
  hww->Fill1DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  EvtSelector( fakeData[6], hqcd, false, scale6 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  //gStyle->SetOptStat("ieruom");
  gStyle->SetOptStat("");
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");
  gStyle->SetOptTitle(0);

  TString hNames[24] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",  
                                   "Mt_4",     "dRmj",    "RelPt",  "Ht_lep",
                                   "Ht_tot",   "Jet1Eta", "Jet2Eta"          } ;

  //TF1* fjb1 = new TF1("fjb1", MassFitFunction::ConvJacob, 0, 140, 7);
  for (int i=0; i<24; i++) {
      //if ( i>1 && i <4 ) continue;
      if ( i == 19 ) continue;

      THStack* hStk = new THStack("hStk", hNames[i] );
      h_wj[i]->SetFillColor(kGreen);
      hStk->Add( h_wj[i] );
      h_zj[i]->SetFillColor(kAzure-2);
      hStk->Add( h_zj[i] );
      h_tq[i]->SetFillColor(kMagenta+2);
      hStk->Add( h_tq[i] );
      h_tw[i]->SetFillColor(kMagenta);
      hStk->Add( h_tw[i] );
      h_ww[i]->SetFillColor(kWhite);
      hStk->Add( h_ww[i] );
      h_qcd[i]->SetFillColor(kYellow);
      hStk->Add( h_qcd[i] );
      h_tt[i]->SetFillColor(kRed+1);
      hStk->Add( h_tt[i] );

      
      // construction for KS test
      TH1D* All_MC  = (TH1D*) h_tt[i]->Clone("combined MC"); 
      All_MC->Add( h_wj[i] );
      All_MC->Add( h_zj[i] );
      All_MC->Add( h_tw[i] );
      All_MC->Add( h_tq[i] );
      All_MC->Add( h_ww[i] );
      All_MC->Add( h_qcd[i] );

      // Jet1 Et Spectrum
      TCanvas* c1 = new TCanvas("c1","", 800, 700);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      double hMax = 1.5* h_dt[i]->GetBinContent(h_dt[i]->GetMaximumBin()) ;
      h_dt[i]->SetMaximum( hMax );
      //if ( i < 4 || i==7 || i== 8 ) c1->SetLogy();
      if ( i==7 || i== 8 ) c1->SetLogy();
      if ( i < 4 || i==7 || i== 8 ) h_dt[i]->SetMaximum( 1.5*hMax );

      if ( i==8 ) h_dt[i]->SetMinimum(10.);
      if ( i==7 ) {
         //h_dt[i]->SetMaximum(1800.);
         //h_dt[i]->SetMinimum(0.);
	 double ndt = h_dt[i]->GetEntries();
	 double ntt = h_tt[i]->GetEntries();
	 double nwj = h_wj[i]->GetEntries();
	 double nzj = h_zj[i]->GetEntries();
	 double ntq = h_tq[i]->GetEntries();
	 double ntw = h_tw[i]->GetEntries();
	 double nww = h_ww[i]->GetEntries();
	 double nqc = h_qcd[i]->GetEntries();
         cout<<" Data= "<< ndt <<" tt= "<< ntt<<" wj= "<<nwj<<" zj= "<< nzj ;
         cout<<" tq= "<<ntq<<" tw= "<<ntw<<" ww= "<<nww<<" qcd= "<<nqc<<endl; 
         
         //for (int k=1; k< 2; k++) {
         //    int kf =  7 ;
         for (int k=1; k< 6; k++) {
             int kf = ( k > 4 ) ? 7 : k ;
             double n_dt = h_dt[i]->Integral(k,kf);
	     double n_tt = h_tt[i]->Integral(k,kf);
	     double n_wj = h_wj[i]->Integral(k,kf);
	     double n_zj = h_zj[i]->Integral(k,kf);
	     double n_tq = h_tq[i]->Integral(k,kf);
	     double n_tw = h_tw[i]->Integral(k,kf);
	     double n_ww = h_ww[i]->Integral(k,kf);
	     double n_qc = h_qcd[i]->Integral(k,kf);
             double n_mc = n_tt + n_wj + n_zj + n_tq + n_tw + n_ww + n_qc  ;
             //double n_bg = n_wj + n_zj + n_tq + n_tw + n_ww + n_qc ;
             //double n_vjbg = n_wj + n_zj + n_tq + n_tw + n_ww  ;
             cout<<" "<<k-1<<"J Data= "<< n_dt <<" MC = "<< n_mc<<" tt = "<< n_tt;
             //cout<<" MC_Bg = "<< n_bg <<" Vj_Bg = "<< n_vjbg <<" QCD_Bg = "<< n_qc <<endl;
             cout<<" wj= "<<n_wj<<" zj= "<< n_zj <<" tq= "<<n_tq<<" tw= "<<n_tw<<" ww= "<<n_ww<<" qcd= "<<n_qc<<endl; 
         }
         
      }
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
      h_dt[i]->Draw("APE SAME");
      //if ( i == 7 ) fjb1->Draw("sames") ;
      c1->Update();

      TLegend *leg = new TLegend(.6, .6, .9, .9 );
      leg->AddEntry(h_tt[i],  "Ttbar",  "F");
      leg->AddEntry(h_qcd[i], "QCD",    "F");
      leg->AddEntry(h_tw[i],  "Single Top(tW)", "F");
      leg->AddEntry(h_tq[i],  "Single Top(t)", "F");
      leg->AddEntry(h_ww[i],  "WW",     "F");
      leg->AddEntry(h_zj[i],  "Z+Jets", "F");
      leg->AddEntry(h_wj[i],  "W+Jets", "F");
      leg->Draw("same");
      c1->Update();

      // KS test block
      /*
      double ksP = All_MC->KolmogorovTest( h_dt[i]  ) ;
      double ks = All_MC->KolmogorovTest( h_dt[i], "M" ) ;
      //double nA = All_MC->Integral();
      //double nD = h_dt[i]->Integral();
      //cout<<" KSP = "<< ksP <<"  KS = "<<ks <<" ks*nA = "<< ks*sqrt(nA) <<" ks*nD = "<< ks*sqrt(nD)<<endl ;
      ostringstream ksStr ;
      ksStr << "KS = ";
      ksStr << setprecision(4) << ksP ;
      ksStr << " D = " ;
      ksStr << setprecision(4) << ks ;
      TString ksReport = ksStr.str() ;
      TPaveText *pvtxt = new TPaveText( 0.25, 0.85, .6, .9,"NDC" );
      //pvtxt->SetFillColor(0);
      pvtxt->SetTextSize(0.03) ;
      TText *tx1 = pvtxt->AddText( ksReport ) ;
      tx1->SetTextAlign(11) ;

      pvtxt->Draw("same");
      c1->Update();
      */

      TString plotname1 = hfolder + theSubFolder+ DataName + hNames[i]+ "."+plotType ;
      //if ( doScale ) plotname1 = hfolder + theSubFolder+ DataName + hNames[i]+ "sc."+plotType ;
      c1->Print( plotname1 );
      delete c1;
      delete hStk ;
  }
 
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hww;
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
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw");
  double scale5 = fitInput->NormalizeComponents( "ww");
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  Reset( 1, doScale ) ;

  hObjs* htt = new hObjs() ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj") ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hzj = new hObjs("zj") ;
  EvtSelector( fakeData[2], hzj, false, scale2);
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq") ;
  EvtSelector( fakeData[3], htq, false, scale3 );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw") ;
  EvtSelector( fakeData[4], htw, false, scale4 );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  hObjs* hww = new hObjs("ww") ;
  EvtSelector( fakeData[5], hww, false, scale5 );
  vector<TH1D*> h_ww ;
  hww->Fill1DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd") ;
  EvtSelector( fakeData[6], hqcd, false, scale6 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

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
      h_ww[i]->SetFillColor(kWhite);
      hStk->Add( h_ww[i] );
      h_tt[i]->SetFillColor(kRed+1);
      hStk->Add( h_tt[i] );

      //All_MC->Add( h_wj[i] );
      All_MC->Add( h_zj[i] );
      All_MC->Add( h_tt[i] );
      All_MC->Add( h_tw[i] );
      All_MC->Add( h_ww[i] );
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
      if ( i==7 ) hStk->SetMinimum(1.);
      c1->cd();

      hStk->Draw();
      All_MC->Draw("sames");
      if ( i== 7 && doWMtFit )  MtFitter( All_MC, *fjb1, c1, 60, 120  ) ;
      //if ( i == 7 && doWMtFit ) fjb1->Draw("sames") ;
      c1->Update();

      TLegend *leg = new TLegend(.6, .6, .88, .88 );
      leg->AddEntry(h_tt[i],  "Ttbar",  "F");
      leg->AddEntry(h_wj[i],  "W+Jets", "F");
      leg->AddEntry(h_zj[i],  "Z+Jets", "F");
      leg->AddEntry(h_tq[i],  "Single Top(t)", "F");
      leg->AddEntry(h_tw[i],  "Single Top(tW)", "F");
      leg->AddEntry(h_ww[i],  "WW",     "F");
      leg->AddEntry(h_qcd[i], "QCD",    "F");
      leg->Draw("same");
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
  delete hww;
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
  double scale2 = ( norm == 0 ) ? fitInput->NormalizeComponents( "zj" ) : norm / fitInput->TreeSize( fakeData[2] );
  double scale6 = ( norm == 0 ) ? fitInput->NormalizeComponents( "qcd" ): norm / fitInput->TreeSize( fakeData[6] );
  /*
  double scale3 = ( norm == 0 ) ? fitInput->NormalizeComponents( "tq" ) : norm / fitInput->TreeSize( fakeData[4] );
  double scale4 = ( norm == 0 ) ? fitInput->NormalizeComponents( "tw" ) : norm / fitInput->TreeSize( fakeData[5] );
  double scale5 = ( norm == 0 ) ? fitInput->NormalizeComponents( "ww" ) : norm / fitInput->TreeSize( fakeData[6] );
  */
  hObjs* htt = new hObjs() ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs( "wj" ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hzj = new hObjs( "zj" ) ;
  EvtSelector( fakeData[2], hzj, false, scale2 );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* hqcd = new hObjs( "qcd" ) ;
  EvtSelector( fakeData[6], hqcd, false, scale6 );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  /*
  hObjs* htq = new hObjs( "tq" ) ;
  EvtSelector( fakeData[3], htq, false, scale3 );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs( "tw" ) ;
  EvtSelector( fakeData[4], htw, false, scale4 );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  hObjs* hww = new hObjs( "ww" ) ;
  EvtSelector( fakeData[5], hww, false, scale5 );
  vector<TH1D*> h_ww ;
  hww->Fill1DVec( h_ww );
  */

  TString hNames[18] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",    "Mt_4" } ;

  double yMax[18]    = { 1000., 1000., 1000., 1000., 1000.,
                                 500.,  100., 3000., 1000.,
                                 300.,  150.,  150,  500., 
                                 600.,  600.,  600., 600.,  600. };

  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize( 0.06, "X");
  gStyle->SetLabelSize( 0.06, "Y");
  gStyle->SetLabelSize( 0.06, "Z");

  for (int i=0; i<18; i++) {
      if ( i>0 && i <4 ) continue;

      h_wj[i]->SetLineWidth(5);
      h_zj[i]->SetLineWidth(5);
      h_qcd[i]->SetLineWidth(3);
      h_tt[i]->SetLineWidth(3);

      h_wj[i]->SetLineColor(kGreen+2);
      //h_zj[i]->SetLineColor(kAzure-2);
      h_zj[i]->SetLineColor(kBlue);

      h_tt[i]->SetLineColor(kRed+1);

      h_qcd[i]->SetLineColor(kBlack);
      /*
      h_tq[i]->SetLineWidth(2);
      h_tw[i]->SetLineWidth(2);
      h_ww[i]->SetLineWidth(2);

      h_tq[i]->SetLineColor(kMagenta+2);
      h_tw[i]->SetLineColor(kMagenta);
      h_ww[i]->SetLineColor(kWhite);
      */
      
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
      if ( i== 0 || i == 4 || i == 7 || i== 8  ) c1->SetLogy();
      //if ( i== 1 || i == 2 || i == 3           ) c1->SetLogy();
      c1->cd();

      yMax[i] = ( norm > 0. ) ? yMax[i] : h_wj[i]->Integral()/2 ;
      h_wj[i]->SetMaximum( yMax[i] );

      h_wj[i]->Draw();
      h_zj[i]->Draw("same");
      //h_tq[i]->Draw("same");
      //h_tw[i]->Draw("same");
      //h_ww[i]->Draw("same");
      h_tt[i]->Draw("same");
      h_qcd[i]->Draw("same");

      
      TLegend *leg = new TLegend(.70, .68, .89, .89 );
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
  delete hwj;
  delete hzj;
  delete hqcd;
 
  //delete htq;
  //delete htw;
  //delete hww;

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

// set 0 : inclusive/exclusive of jet multiplicity  1: ON/OFF of event weighting 
void ObjectInfo::Reset( int idx, bool newValue ){
  if ( idx != -1 ){
     if ( idx == 0 ) inclu = newValue ;
     if ( idx == 1 ) isWeight = newValue ;
  } else {
     fitInput->GetParameters( "Inclusive", &Inclusive );
     isWeight = false ;
  }
}

// set 0: min number of jet  1: JES type
void ObjectInfo::Reset( int idx, int newValue ){
  if ( newValue != -1 ){
     if ( idx == 0 ) n_Jets =  newValue ;
     if ( idx == 1 ) JESType =  newValue ;
  } else {
     if ( idx == 0 ) fitInput->GetParameters( "n_Jets",   &n_Jets );
     if ( idx == 1 ) fitInput->GetParameters( "JESType",   &JESType );
  }
}

void ObjectInfo::Reset( int idx, double newValue ){
  if ( newValue != -1 ){
     if ( idx == 0 ) jetCuts[0] = newValue ;   // pt
     if ( idx == 1 ) jetCuts[1] = newValue ;   // eta
     if ( idx == 2 ) muonCuts[0] = newValue ;  // pt
     if ( idx == 3 ) muonCuts[1] = newValue ;  // eta
     if ( idx == 4 ) muonCuts[2] = newValue ;  // iso
  } else {
        fitInput->GetParameters( "JetCuts",  &jetCuts );
        fitInput->GetParameters( "MuonCuts", &muonCuts );
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

    // PF
    // double qSet[5] = { 1.45, 2.15, 1.20, 1.40, 5.91  }; // with V+Jets scaling
    //double qSet[5] = { 1.46, 2.16, 1.20, 1.43, 6.2  };
    //double vSet[5] = { 1.02, 1.03, 0.82, 1.16, 1.36 };
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

    vector<double> qSet;
    fitInput->GetParameters( "qcdNorm", &qSet );
    vector<double> vSet;
    fitInput->GetParameters( "vjNorm", &vSet );

    double theScale = 1 ;
    if ( fileName.substr(0,2) == "wj" || fileName.substr(0,2) == "zj" || fileName.substr(0,2) == "ww" ) {
       if ( NJets == 1 ) theScale =   vSet[0] ;
       if ( NJets >= 2 ) theScale =   vSet[1] ;
       //if ( NJets == 3 ) theScale =   1.01 ;
       //if ( NJets >= 4 ) theScale =  vSet[1]*1.0346 ;
    }
    if ( fileName.substr(0,2) == "tq" || fileName.substr(0,2) == "tw"  ) {
       if ( NJets == 1 ) theScale =   vSet[0] ;
       if ( NJets >= 2 ) theScale =   vSet[1] ;
       //if ( NJets >= 4 ) theScale =  vSet[1]*1.0346 ;
       //if ( NJets == 3 ) theScale =   1.01 ;
       //if ( NJets >= 4 ) theScale =  (vSet[1]*vSet[1]*vSet[1])/(vSet[0]*vSet[0]) ;
    }
    if ( fileName.substr(0,2) == "qc" ) {
       if ( NJets == 1 ) theScale =   qSet[0] ;
       if ( NJets >= 2 ) theScale =   qSet[1] ;
       //if ( NJets == 3 ) theScale =  (qSet[1]*qSet[1])/qSet[0] ;
       //if ( NJets >= 4 ) theScale =  qSet[1]*1.0346 ;
    }
    //if ( fileName.substr(0,2) == "tt" ) {
    //   if ( NJets == 3 ) theScale =  (160.4/157.5)  ;
    //}
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

    vector<double> vSet;
    fitInput->GetParameters( "vjNorm", &vSet );

    double theScale = 1 ;
    if ( fileName.substr(0,2) == "wj" || fileName.substr(0,2) == "zj" || fileName.substr(0,2) == "ww" ) {
       if ( NJets == 1 ) theScale =   vSet[0] ;
       if ( NJets >= 2 ) theScale =   vSet[1] ;
       //if ( NJets == 3 ) theScale =   vSet[1]/vSet[0] ;
       //if ( NJets >= 4 ) theScale =  (vSet[1]*vSet[1])/(vSet[0]*vSet[0]) ;
    }
    return theScale ;
}

// scaleMode = 0 : no scaling , scaleMode = 1 : only scale Vjets componenet , scaleMode = 2 : regular scaling 
void ObjectInfo::QCDBGPlotter( string DataName, vector<string>& fakeData, bool doQCD, double MtCut, double METCut, int scaleMode ){

  ResetBGCuts(0, MtCut );
  ResetBGCuts(1, METCut );
  if (scaleMode > 0 ) Reset(1, true ) ; // normalize MC with Data

  TString theFolder = hfolder ;
  TString theSubFolder = "hQCDMu/" ;
  if ( !doQCD ) theSubFolder = "hSGMu/" ;

  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 2 ;

  hObjs* hdt = new hObjs("data", nbin ) ;
  QCDSelector( DataName, hdt, false, 1, doQCD );
  vector<TH1D*> h_dt ;
  hdt->Fill1DVec( h_dt );

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  hObjs* htt = new hObjs("tt", nbin ) ;
  QCDSelector( fakeData[0], htt, false, scale0, doQCD, scaleMode );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  QCDSelector( fakeData[1], hwj, false, scale1, doQCD, scaleMode );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  QCDSelector( fakeData[2], hzj, false, scale2, doQCD, scaleMode );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  QCDSelector( fakeData[3], htq, false, scale3, doQCD, scaleMode );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  QCDSelector( fakeData[4], htw, false, scale4, doQCD, scaleMode );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

  hObjs* hww = new hObjs("ww", nbin ) ;
  QCDSelector( fakeData[5], hww, false, scale5, doQCD, scaleMode );
  vector<TH1D*> h_ww ;
  hww->Fill1DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  QCDSelector( fakeData[6], hqcd, false, scale6, doQCD, scaleMode );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  //gStyle->SetOptStat("ieruom");
  gStyle->SetOptStat("");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize( 0.05, "X");
  gStyle->SetLabelSize( 0.05, "Y");
  gStyle->SetLabelSize( 0.05, "Z");

  TString hNames[21] = { "Jet1Pt", "Jet2Pt",   "Jet3Pt",  "Jet4Pt", "MuPt",
                                   "MET",      "MuEta",   "NJets",  "MuIso",
                                   "MuNHits",  "MuD0",    "MuX2",   "lepM2Pt", 
                                   "Mt",       "Mt_1",    "Mt_2",   "Mt_3",    
                                   "Mt_4",    "dRmin",    "PtRel",  "Ht_lep" } ;

   /* double hMax[21] =    { 20000,  3500, 20000, 20000,  9500,
                                 5000,  3700, 45000, 25000, 
                                 7500,  4000,  8000, 10000,
                                 4000,  2800,   600,   140,
                                  140,  5000, 10000, 10000 } ;
  */
  for (int i=0; i<21; i++) {
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
 
      // construction for KS test
      TH1D* All_MC  = (TH1D*) h_tt[i]->Clone("combined MC"); 
      All_MC->Add( h_wj[i] );
      All_MC->Add( h_zj[i] );
      All_MC->Add( h_tw[i] );
      All_MC->Add( h_tq[i] );
      All_MC->Add( h_ww[i] );
      All_MC->Add( h_qcd[i] );

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
      double hMax = 1.5* h_dt[i]->GetBinContent(h_dt[i]->GetMaximumBin()) ;
      h_dt[i]->SetMaximum( hMax );
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

      TLegend *leg = new TLegend(.6, .6, .9, .9 );
      leg->AddEntry(h_tt[i],  "Tt",      "F");
      leg->AddEntry(h_wj[i],  "WJets",   "F");
      leg->AddEntry(h_zj[i],  "ZJets",   "F");
      leg->AddEntry(h_tq[i],  "Single Top(t)",  "F");
      leg->AddEntry(h_tw[i],  "Single Top(tW)", "F");
      leg->AddEntry(h_ww[i],  "WW",      "F");
      leg->AddEntry(h_qcd[i], "QCD",     "F");
      leg->Draw("same");
      c1->Update();

      // KS test block
      /*
      double ksP = All_MC->KolmogorovTest( h_dt[i]  ) ;
      //double P_val = h_dt[i]->Chi2Test( All_MC, "UW" ) ;
      //cout<<" KSP = "<< ksP <<"  P-value = "<< P_val << endl ;
      ostringstream ksStr ;
      ksStr << "KS Prob. = ";
      ksStr << setprecision(4) << ksP ;
      //ksStr << " P-Val = " ;
      //ksStr << P_val ;
      TString ksReport = ksStr.str() ;
      TPaveText *pvtxt = new TPaveText( 0.25, 0.85, .6, .9,"NDC" );
      //pvtxt->SetFillColor(0);
      pvtxt->SetTextSize(0.03) ;
      pvtxt->AddText( ksReport ) ;
      pvtxt->SetTextAlign(11) ;
      pvtxt->Draw("same");
      c1->Update();
      */
      TString plotname1 = hfolder + theSubFolder + DataName + hNames[i]+ "."+plotType ;
      if (scaleMode > 0 ) plotname1 = hfolder + theSubFolder + DataName + hNames[i]+ "_scaled."+plotType ;
      c1->Print( plotname1 );
    
      delete c1;
      delete hStk ;
      delete leg ;
      delete All_MC ;
  }
 
  delete hdt;
  delete htt;
  delete htq;
  delete htw;
  delete hwj;
  delete hww;
  delete hzj;
  delete hqcd;

}

void ObjectInfo::QCDBG2DPlotter( string DataName, vector<string>& fakeData, bool doQCD, double MtCut, double METCut, double scale ){

  ResetBGCuts(0, MtCut );
  ResetBGCuts(1, METCut );

  TString theFolder = hfolder ;
  TString theSubFolder = "hQCDMu/" ;
  if ( !doQCD ) theSubFolder = "hSGMu/" ;
  if ( scale > 1 ) theSubFolder = "hSGMu_Scale1/" ;
  if ( scale == -1 ) theSubFolder = "hSGMu_Scale4/" ;

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

void ObjectInfo::Data2DPlotter( string DataName, vector<string>& fakeData, bool doScale ){

  TString theFolder = hfolder ;
  TString theSubFolder = ( doScale == false ) ?  "hData/" : "hData_Scale/" ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( theSubFolder );
  gSystem->cd( "../" );

  int nbin = 2 ;

  hObjs* hdt = new hObjs("data", nbin ) ;
  EvtSelector( DataName, hdt );
  vector<TH2D*> h_dt ;
  hdt->Fill2DVec( h_dt );

  double scale0 = fitInput->NormalizeComponents( "tt" );
  double scale1 = fitInput->NormalizeComponents( "wj" );
  double scale2 = fitInput->NormalizeComponents( "zj" );
  double scale3 = fitInput->NormalizeComponents( "tq" );
  double scale4 = fitInput->NormalizeComponents( "tw" );
  double scale5 = fitInput->NormalizeComponents( "ww" );
  double scale6 = fitInput->NormalizeComponents( "qcd" );

  Reset( 1, doScale ) ;

  hObjs* htt = new hObjs("tt", nbin ) ;
  EvtSelector( fakeData[0], htt, false, scale0 );
  vector<TH2D*> h_tt ;
  htt->Fill2DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin ) ;
  EvtSelector( fakeData[1], hwj, false, scale1 );
  vector<TH2D*> h_wj ;
  hwj->Fill2DVec( h_wj );

  hObjs* hzj = new hObjs("zj", nbin ) ;
  EvtSelector( fakeData[2], hzj, false, scale2 );
  vector<TH2D*> h_zj ;
  hzj->Fill2DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin ) ;
  EvtSelector( fakeData[3], htq, false, scale3 );
  vector<TH2D*> h_tq ;
  htq->Fill2DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin ) ;
  EvtSelector( fakeData[4], htw, false, scale4 );
  vector<TH2D*> h_tw ;
  htw->Fill2DVec( h_tw );

  hObjs* hww = new hObjs("ww", nbin ) ;
  EvtSelector( fakeData[5], hww, false, scale5 );
  vector<TH2D*> h_ww ;
  hww->Fill2DVec( h_ww );

  hObjs* hqcd = new hObjs("qcd", nbin ) ;
  EvtSelector( fakeData[6], hqcd, false, scale6 );
  vector<TH2D*> h_qcd ;
  hqcd->Fill2DVec( h_qcd );

  //gStyle->SetOptStat("ieruom");
  //gStyle->SetStatX(0.85);
  gStyle->SetOptStat("");
  //gStyle->SetOptTitle(1);
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(8);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleX(0.14);
  gStyle->SetTitleW(0.36);

  TString hNames[5] = { "MET_Mt", "WPt_Mt", "WPt_dPhi",  "MuMET", "Mt_dPhi" } ;
  TString xtitle[5] = { "MET",    "Pt (#mu , #nu)", "Pt (#mu , #nu)", "Pt of #mu", "M_{T}"} ;
  TString ytitle[5] = { " M_{T} ", "M_{T} ", "#Delta#phi (#mu , #nu) ", "MET ", "#Delta#phi (#mu , #nu) " } ;
  //TPaveText* hTitle = new TPaveText( 0.7, 0.7, 0.85, 0.85, "NDC" );

  for (int i=0; i<5; i++) {

      TH2D* h_all =  (TH2D*) h_tw[i]->Clone("combinedMC");
      h_all->Add( h_qcd[i] );
      h_all->Add( h_wj[i] );
      h_all->Add( h_zj[i] );
      h_all->Add( h_tt[i] );
      h_all->Add( h_tq[i] );
      h_all->Add( h_ww[i] );

      // Jet1 Et Spectrum
      TCanvas* c1 = new TCanvas("c1","", 900, 700);
      c1->SetGrid();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->Divide(2,2);

      c1->cd(1);
      h_dt[i]->SetTitle( "Data" );
      h_dt[i]->SetXTitle( xtitle[i] );
      h_dt[i]->SetYTitle( ytitle[i] );
      h_dt[i]->SetTitleSize(0.07, "X") ;
      h_dt[i]->SetTitleSize(0.07, "Y") ;
      h_dt[i]->SetTitleOffset(1.0, "X") ;
      h_dt[i]->SetTitleOffset(1.0, "Y") ;
      h_dt[i]->SetLabelSize(0.06,"X");
      h_dt[i]->SetLabelSize(0.06,"Y");
      h_dt[i]->SetLabelSize(0.06,"Z");
      h_dt[i]->GetXaxis()->SetNdivisions(507) ;
      h_dt[i]->GetYaxis()->SetNdivisions(507) ;
      h_dt[i]->Draw("COLZ");
      c1->Update();
      /*
      TText* hTitle1 = new TText( 0.6, 0.75, "Data" );
      hTitle1->SetNDC(kTRUE);
      hTitle1->SetTextSize(0.1);
      hTitle1->Draw("same");
      */

      c1->cd(2);
      h_wj[i]->SetTitle( "W+Jets" );
      h_wj[i]->SetXTitle( xtitle[i] );
      h_wj[i]->SetYTitle( ytitle[i] );
      h_wj[i]->SetTitleSize(0.07, "X") ;
      h_wj[i]->SetTitleSize(0.07, "Y") ;
      h_wj[i]->SetTitleOffset(1.0, "X") ;
      h_wj[i]->SetTitleOffset(1.0, "Y") ;
      h_wj[i]->SetLabelSize(0.06,"X");
      h_wj[i]->SetLabelSize(0.06,"Y");
      h_wj[i]->SetLabelSize(0.06,"Z");
      h_wj[i]->GetXaxis()->SetNdivisions(507) ;
      h_wj[i]->GetYaxis()->SetNdivisions(507) ;
      h_wj[i]->Draw("COLZ");
      c1->Update();

      c1->cd(3);
      h_qcd[i]->SetTitle( "QCD" );
      h_qcd[i]->SetXTitle( xtitle[i] );
      h_qcd[i]->SetYTitle( ytitle[i] );
      h_qcd[i]->SetTitleSize(0.07, "X") ;
      h_qcd[i]->SetTitleSize(0.07, "Y") ;
      h_qcd[i]->SetTitleOffset(1.0, "X") ;
      h_qcd[i]->SetTitleOffset(1.0, "Y") ;
      h_qcd[i]->SetLabelSize(0.06,"X");
      h_qcd[i]->SetLabelSize(0.06,"Y");
      h_qcd[i]->SetLabelSize(0.06,"Z");
      h_qcd[i]->GetXaxis()->SetNdivisions(507) ;
      h_qcd[i]->GetYaxis()->SetNdivisions(507) ;
      h_qcd[i]->Draw("COLZ");
      c1->Update();

      /*
      c1->cd(4);
      h_all->SetTitle( "Mixed MC " );
      h_all->SetXTitle( xtitle[i] );
      h_all->SetYTitle( ytitle[i] );
      h_all->SetTitleSize(0.07, "X") ;
      h_all->SetTitleSize(0.07, "Y") ;
      h_all->SetTitleOffset(1.0, "X") ;
      h_all->SetTitleOffset(1.0, "Y") ;
      h_all->SetLabelSize(0.06,"X");
      h_all->SetLabelSize(0.06,"Y");
      h_all->SetLabelSize(0.06,"Z");
      h_all->GetXaxis()->SetNdivisions(507) ;
      h_all->GetYaxis()->SetNdivisions(507) ;
      h_all->Draw("COLZ");
      c1->Update();
      */

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

