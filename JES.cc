#include "JES.h"

JES::JES() {

  fitInput = new MassAnaInput();
  fitFunc  = new MassFitFunction();

  fitInput->Initialize( &hfolder  );
  fitInput->GetParameters( "MassLBound", &mL);
  fitInput->GetParameters( "MassHBound", &mH);

}

JES::~JES(){

  delete fitInput ;
  delete fitFunc;

}

double JES::WMSC( TLorentzVector v1, TLorentzVector v2, double jes ){

   TLorentzVector sjv1( v1.Px()*jes, v1.Py()*jes, v1.Pz()*jes, v1.E()*jes );
   TLorentzVector sjv2( v2.Px()*jes, v2.Py()*jes, v2.Pz()*jes, v2.E()*jes );

   TLorentzVector wp4 = sjv1 + sjv2;

   return wp4.M() ;
}

double JES::deltaR( TLorentzVector v1, TLorentzVector v2 ){

   double dA = v1.Eta() - v2.Eta();
   double dF = v1.Phi() - v2.Phi();
   double dR = sqrt( (dA*dA) + (dF*dF) ) ;

   return dR ;
}

void JES::WMassSpectrum( TString mName, TH1D* hWM0, TH1D* hWM, TH1D* hMCW, double jes ) {

  // get the file names of fake data
  //TString fNameList[5];
  //fitInput->GetFileName( mName, 1, fNameList );
  //TFile* file = TFile::Open( fNameList[0] );

  TString theFileName ;
  vector<string> msets;
  fitInput->GetParameters( "TMassAssumption", &msets );
  vector<string> sglist;
  fitInput->GetParameters( "Signals", &sglist );
  for (size_t i=0; i< msets.size(); i++) {
      TString massump = msets[i].substr(0,3) ;
      if ( massump == mName ) theFileName = sglist[i] + ".root" ;
  }
  TFile* file = TFile::Open( theFileName );

  TTree*  tr0 = (TTree*) file->Get( "mcmTt" );  
  TTree*  tr1 = (TTree*) file->Get( "solTt" );
  TTree*  tr2 = (TTree*) file->Get( "selJet" );

  // this is only for MC Matching
  int evtId0, entSz0, jid0[5] ;
  tr0->SetBranchAddress( "evtId" ,&evtId0 );
  tr0->SetBranchAddress( "JId"   ,&jid0   );
  tr0->SetBranchAddress( "entSz" ,&entSz0);
  
  int evtId1, entId1, iniId1, entSz1, jid1[5] ;
  tr1->SetBranchAddress( "evtId" ,&evtId1 );
  tr1->SetBranchAddress( "entId" ,&entId1 );
  tr1->SetBranchAddress( "entSz" ,&entSz1 );
  tr1->SetBranchAddress( "iniId" ,&iniId1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );

  double px,py,pz,E, bDis;
  int evtId2 ;
  tr2->SetBranchAddress("px"    ,&px);
  tr2->SetBranchAddress("py"    ,&py);
  tr2->SetBranchAddress("pz"    ,&pz);
  tr2->SetBranchAddress("E"     ,&E);
  tr2->SetBranchAddress("evtId" ,&evtId2);
  tr2->SetBranchAddress("qCut"  ,&bDis   );

  int evtId = 0 ;
  int nbtag = 0 ;
  std::vector<TLorentzVector> jpv ;
  std::vector<double> blist ;
  int it1 = 1;
  int it0 = 1;
  tr0->GetEntry(it0);
  tr1->GetEntry(it1);
  for (int k=1; k<= tr2->GetEntries() ; k++) {
      tr2->GetEntry(k);
      //cout<<" objs event id = "<< evtId2 <<endl;

      if ( evtId != evtId2 ) {
         //cout<<" evtId = "<< evtId<<" evtId0 = "<< evtId0<<"  it1: "<<it1 <<endl;
         // looking for MC Matched permutation
         
         tr0->GetEntry(it0);
         if ( evtId !=0 && evtId == evtId0 && hMCW != NULL ) {
            //cout<<" n of permu = "<< entSz0 <<" j0:"<<jid0[0]<<" j1:"<<jid0[1]<<endl;
            if ( jpv[ jid0[0] ].Pt() >= 30 && jpv[ jid0[1] ].Pt() >= 30 ) {
            //if ( entSz0 == 24 ) {
               double mcwmass = WMSC(  jpv[ jid0[0] ], jpv[ jid0[1] ], jes ) ;
               hMCW->Fill( mcwmass ); 
            }
            it0++ ;
         }

         // skip no solution events
         tr1->GetEntry(it1);
         bool noSolution = true;
         if ( evtId == evtId1 ) {
            it1 = it1 + entSz1 ;
            noSolution = false;
         }
 
         // W mass spectrum
         int nW = 0 ;
         for (size_t i1=0; i1 < jpv.size()-1; i1++ ) {
             if ( noSolution ) break;
             for (size_t i2=i1+1; i2<jpv.size(); i2++ ) {
                 double wmass = WMSC(  jpv[i1], jpv[i2], jes ) ;
                 hWM0->Fill( wmass );
                 if ( blist[i1] < 5 && blist[i2] < 5 && nbtag > 0 ) {
                    hWM->Fill( wmass );
                    nW++;
                 }
             }
         }
         //cout<<" nW = "<< nW <<" nbtag = "<< nbtag<< endl;
         evtId = evtId2 ;
         jpv.clear();
         blist.clear();
         nbtag = 0;

      }
       
      TLorentzVector jp1( px, py, pz, E );
      jpv.push_back( jp1 );
      blist.push_back( bDis );
      if ( bDis > 5 ) nbtag++;
  }

}

void JES::GetJESFromW( TString mName, Double_t* para, double jes ){

  
  int jesdigi = static_cast<int>(jes*100) ;
  cout<<" jes digig = "<< jesdigi <<"  -> "<<jes <<endl;
  char jname[6] = { 0 };
  sprintf( jname, "_%d_", jesdigi );
  

  TH1D* ttW0 = new TH1D("ttW0"," W Mass from Tt w/o Btag  ",   40, 0, 200 );
  TH1D* ttW1 = new TH1D("ttW1"," W Mass from Tt w/  Btag  ",   20, 0, 200 );
  TH1D* ttW2 = new TH1D("ttW2"," W Mass from Tt MC matched",   40, 0, 200 );
  TH1D* wjW0 = new TH1D("wjW0"," W Mass from Wj w/o Btag  ",   40, 0, 200 );
  TH1D* wjW1 = new TH1D("wjW1"," W Mass from Wj w/  Btag  ",   20, 0, 200 );

  WMassSpectrum(    mName, ttW0, ttW1, ttW2, jes );
  WMassSpectrum( "_WJets", wjW0, wjW1, NULL, jes );
  
  fitInput->NormalizeComponents( "tt", ttW0 );
  fitInput->NormalizeComponents( "wj", wjW0 );
  fitInput->NormalizeComponents( "tt", ttW1 );
  fitInput->NormalizeComponents( "wj", wjW1 );
  fitInput->NormalizeComponents( "tt", ttW2 );
  /*
  fitInput->NormalizeComponents( 100,  12000,  1, ttW0 );
  fitInput->NormalizeComponents( 100, 120000,  2, wjW0 );
  fitInput->NormalizeComponents( 100,  12000,  1, ttW1 );
  fitInput->NormalizeComponents( 100, 120000,  2, wjW1 );
  fitInput->NormalizeComponents( 100,  12000,  1, ttW2 );
  */
  // add them
  TH1D* allW0 = new TH1D("allW0","", 40, 0, 200 );
  allW0->Add( ttW0 );
  allW0->Add( wjW0 );

  TH1D* allW1 = new TH1D("allW1","", 20, 0, 200 );
  allW1->Add( ttW1 );
  allW1->Add( wjW1 );

  // stack them
  THStack* wStk0 = new THStack("wStk0", "Combined W ");
  wjW0->SetFillColor(2);
  ttW0->SetFillColor(3);
  wStk0->Add( wjW0 );
  wStk0->Add( ttW0 );

  THStack* wStk1 = new THStack("wStk1", "Combined W ");
  wjW1->SetFillColor(2);
  ttW1->SetFillColor(3);
  wStk1->Add( wjW1 );
  wStk1->Add( ttW1 );

  cout<<" wj0 # = "<< wjW0->Integral() << " wj1 # = "<< wjW1->Integral() << endl;

  gStyle->SetOptStat("neim");
  gStyle->SetOptFit(111);
  gStyle->SetStatFontSize(0.05);
  TCanvas* c1 = new TCanvas("c1","", 800, 700);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);

  c1->cd(1);

  wStk0->Draw(); 
  TLegend *leg = new TLegend(.70, .6, .95, .65);
  leg->AddEntry(ttW0, "Tt", "f");
  leg->AddEntry(wjW0, "wjets", "f");
  leg->Draw("same");
  c1->Update();

  TF1* funcW0 = new TF1("funcW0", MassFitFunction::fitLD , 40, 180, 3);
  funcW0->SetParLimits(0, 10., allW0->Integral() );
  funcW0->SetParLimits(1, 45, 125.);
  funcW0->SetParLimits(2,  10., 200.);
  funcW0->SetLineColor(1);
  allW0->Fit( funcW0, "R","sames", 40, 180 );
  c1->Update();

  if ( para != NULL ) para[0] = funcW0->GetParameter(1);
  if ( para != NULL ) para[1] = funcW0->GetParameter(2);

  c1->cd(2);
  wStk1->Draw(); 

  TF1* funcW1 = new TF1("funcW1", MassFitFunction::fitGS , 40, 140, 3);
  funcW1->SetParLimits(0, 10., allW1->Integral() );
  funcW1->SetParLimits(1, 45, 125.);
  funcW1->SetParLimits(2,  10., 200.);
  funcW1->SetLineColor(1);
  allW1->Fit( funcW1, "R","sames", 40, 140 );
  c1->Update();

  if ( para != NULL ) para[2] = funcW1->GetParameter(1);
  if ( para != NULL ) para[3] = funcW1->GetParameter(2);

  c1->cd(3);
  ttW2->Draw();

  TF1* funcW2 = new TF1("funcW2", MassFitFunction::fitGS , 40, 140, 3);
  funcW2->SetParLimits(0, 5., ttW2->Integral() );
  funcW2->SetParLimits(1, 50, 110.);
  funcW2->SetParLimits(2, 5., 100.);
  ttW2->Fit( funcW2, "R","sames", 50*jes,110*jes );
  c1->Update();

  if ( para != NULL ) para[4] = funcW2->GetParameter(1);
  if ( para != NULL ) para[5] = funcW2->GetParameter(2);

  c1->Print(hfolder+mName+jname+"_AllWMass.gif");
  //c1->Print(hfolder+mName+"_AllWMass.gif");
 
  delete c1;
  delete funcW0;
  delete funcW1;
  delete funcW2;

}

void JES::CalibJES( TString mName ){

     Double_t para[6];
     Double_t p0[7];
     Double_t p1[7];
     Double_t p2[7];
     Double_t p3[7];
     Double_t p4[7];
     Double_t p5[7];
     Double_t x[7];
     Double_t xErr[7] = { 0.0 };
     for ( int i =0; i < 7; i ++) {
         double theJES = 0.94 + i*0.02 ;
         GetJESFromW( mName, para, theJES  );
         p0[i] = para[0] ;
         p1[i] = para[1] ;
         p2[i] = para[2] ;
         p3[i] = para[3] ;
         p4[i] = para[4] ;
         p5[i] = para[5] ;
         x[i] = theJES ;
     }

     c2 = new TCanvas("c2","", 800, 600);
     c2->SetFillColor(10);
     c2->SetFillColor(10);
     gPad->SetGridx();
     gPad->SetGridy();
     c2->Divide(2,2);

     c2->cd(1);
     gStyle->SetOptStat(kTRUE);
     gStyle->SetOptFit(111);

     TGraphErrors* hWJES0 = new TGraphErrors( 7, x, p0, xErr, p1 );
     hWJES0->SetMarkerColor(4);
     hWJES0->SetMarkerStyle(20);
     hWJES0->SetLineWidth(2);
     hWJES0->SetLineColor(4);
     hWJES0->SetTitle("No B Tagging");
     hWJES0->Draw("AP");

     TF1* func0 = new TF1("func0", MassFitFunction::fitPoly , 0.9, 1.1, 2);
     hWJES0->Fit( func0, "RQ","sames", x[0]-0.2, x[6]+0.2 );

     c2->Update();

     c2->cd(2);
     TGraphErrors* hWJES1 = new TGraphErrors( 7, x, p2, xErr, p3 );
     hWJES1->SetMarkerColor(4);
     hWJES1->SetMarkerStyle(20);
     hWJES1->SetLineWidth(2);
     hWJES1->SetLineColor(4);
     hWJES1->SetTitle("B Tagged Events");
     hWJES1->Draw("AP");

     hWJES1->Fit( func0, "RQ","sames", x[0]-0.2, x[6]+0.2 );

     c2->Update();

     c2->cd(3);
     TGraphErrors* hWJES2 = new TGraphErrors( 7, x, p4, xErr, p5 );
     hWJES2->SetMarkerColor(4);
     hWJES2->SetMarkerStyle(20);
     hWJES2->SetLineWidth(2);
     hWJES2->SetLineColor(4);
     hWJES2->SetTitle("MC Matched");
     hWJES2->Draw("AP");

     hWJES2->Fit( func0, "RQ","sames", x[0]-0.2, x[6]+0.2 );

     c2->Update();
     c2->Print(hfolder+mName+"_JESTest.gif");

}

void JES::Matching( TString mName, double jes, TH1D* hMatchW ){

  // get the file names of fake data
  //TString fNameList[5];
  //fitInput->GetFileName( mName, 1, fNameList );
  //TFile* file = TFile::Open( fNameList[0] );

   TString theFileName ;
   vector<string> msets;
   fitInput->GetParameters( "TMassAssumption", &msets );
   vector<string> sglist;
   fitInput->GetParameters( "Signals", &sglist );
   for (size_t i=0; i< msets.size(); i++) {
       TString massump = msets[i].substr(0,3) ;
       if ( massump == mName ) theFileName = sglist[i] + ".root" ;
   }
  TFile* file = TFile::Open( theFileName );
  TTree*  tr0 = (TTree*) file->Get( "mcmTt" );  
  TTree*  tr1 = (TTree*) file->Get( "selJet" );
  TTree*  tr2 = (TTree*) file->Get( "gen" );

  // this is only for MC Matching
  int evtId0, entSz0, jid0[5] ;
  tr0->SetBranchAddress( "evtId" ,&evtId0 );
  tr0->SetBranchAddress( "JId"   ,&jid0   );
  tr0->SetBranchAddress( "entSz" ,&entSz0);
  
  double px1,py1,pz1,E1 ;
  int evtId1, objId1 ;
  tr1->SetBranchAddress("px"    ,&px1);
  tr1->SetBranchAddress("py"    ,&py1);
  tr1->SetBranchAddress("pz"    ,&pz1);
  tr1->SetBranchAddress("E"     ,&E1);
  tr1->SetBranchAddress("evtId" ,&evtId1);
  tr1->SetBranchAddress("objId" ,&objId1);

  // parton p4
  double px2,py2,pz2,E2 ;
  int evtId2, objId2 ;
  tr2->SetBranchAddress("px"    ,&px2);
  tr2->SetBranchAddress("py"    ,&py2);
  tr2->SetBranchAddress("pz"    ,&pz2);
  tr2->SetBranchAddress("E"     ,&E2);
  tr2->SetBranchAddress("evtId" ,&evtId2);
  tr2->SetBranchAddress("objId" ,&objId2);

  std::vector<TLorentzVector> jpv ;
  std::vector<TLorentzVector> qpv ;
  int j1=1;
  int j2=1;
  for (int k=1; k<= tr0->GetEntries() ; k++) {
      tr0->GetEntry(k);
      //cout<<" objs event id = "<< evtId2 <<endl;

      for (int j=j1; j<= tr1->GetEntries() ; j++) {
          tr1->GetEntry(j);
          if ( evtId0 != evtId1 ) continue;
          if ( objId1 == jid0[0]) {
             TLorentzVector jp0( px1, py1, pz1, E1 );
             jpv.push_back( jp0 );
          }
          if ( objId1 == jid0[1]) {
             TLorentzVector jp1( px1, py1, pz1, E1 );
             jpv.push_back( jp1 );
          }
          if ( jpv.size() == 2 ) {
             j1 = j;
             break ;
          }
      }
   
      for (int j=j2; j<= tr2->GetEntries() ; j++) {
          tr2->GetEntry(j);
          if ( evtId0 != evtId2 ) continue;
          if ( abs(objId2) < 5 ) {
             TLorentzVector qp0( px2, py2, pz2, E2 );
             qpv.push_back( qp0 );
          }
          if ( qpv.size() == 2 ) {
             j2 = j;
             break ;
          }
      }

      double w1 = WMSC( jpv[0], jpv[1], jes ) ;
      //double w2 = WMSC( qpv[0], qpv[1], jes ) ;
      double dR11 = deltaR( jpv[0],  qpv[0] ); 
      double dR22 = deltaR( jpv[1],  qpv[1] ); 
      double dR12 = deltaR( jpv[0],  qpv[1] ); 
      double dR21 = deltaR( jpv[1],  qpv[0] ); 
      double DR1 =  sqrt( dR11*dR11 + dR22*dR22 ) ;    
      double DR2 =  sqrt( dR12*dR12 + dR21*dR21 ) ;

      /*
      cout<<" evt:"<< evtId1 <<" wMass1 = "<< w1 <<endl;
      cout<<" jpv0:("<< jpv[0].Px()<<","<< jpv[0].Py()<<"," << jpv[0].Pz()<<")"<<endl;
      cout<<" jpv1:("<< jpv[1].Px()<<","<< jpv[1].Py()<<"," << jpv[1].Pz()<<")"<<endl;
      cout<<" evt:"<< evtId2 <<" wMass2 = "<< w2 <<endl;
      cout<<" qpv0:("<< qpv[0].Px()<<","<< qpv[0].Py()<<"," << qpv[0].Pz()<<")"<<endl;
      cout<<" qpv1:("<< qpv[1].Px()<<","<< qpv[1].Py()<<"," << qpv[1].Pz()<<")"<<endl;
      cout<<" dR11:"<<dR11<<" dR22:"<<dR22 <<"  DR => "<< DR1 <<endl;     
      cout<<" dR12:"<<dR12<<" dR21:"<<dR21 <<"  DR => "<< DR2 <<endl;     
      if ( DR1 < 0.7 || DR2 < 0.7 ) cout <<" *** good W matching ("<< w1 <<" , "<<w2 <<" )"<<endl ;
      cout<<" -------------------------------------------------------- "<<endl;
      */
      bool matched = ( DR1 < 0.7 || DR2 < 0.7 ) ? true : false ;
      if ( matched && hMatchW != NULL ) hMatchW->Fill( w1 );

      jpv.clear();
      qpv.clear();

  }

}

void JES::Smearing( TString mName, double sgm, TH1D* hSmearW ){

     // get the file names of fake data
     //TString fNameList[5];
     //fitInput->GetFileName( mName, 1, fNameList );

     TString theFileName ;
     vector<string> msets;
     fitInput->GetParameters( "TMassAssumption", &msets );
     vector<string> sglist;
     fitInput->GetParameters( "Signals", &sglist );
     for (size_t i=0; i< msets.size(); i++) {
         TString massump = msets[i].substr(0,3) ;
         if ( massump == mName ) theFileName = sglist[i] + ".root" ;
     }

     TFile* file = TFile::Open( theFileName );
     TTree* tr = (TTree*) file->Get( "gen" );

     double px,py,pz,E ;
     int evtId, objId ;
     tr->SetBranchAddress("px"    ,&px);
     tr->SetBranchAddress("py"    ,&py);
     tr->SetBranchAddress("pz"    ,&pz);
     tr->SetBranchAddress("E"     ,&E);
     tr->SetBranchAddress("evtId" ,&evtId);
     tr->SetBranchAddress("objId" ,&objId);

     TH1D* hGen = new TH1D("hGen"," Generated E ", 50 ,  0., 250. );
     TH1D* hSmr = new TH1D("hSmr"," Smearing E  ", 50 ,  0., 250. );
     TH1D* hRes = new TH1D("hRes"," E Resolution", 50 ,  -2., 2. );
     TH1D* hgW  = new TH1D("hgW"," generated W ", 40 ,  30., 150. );
     TH1D* hsW  = new TH1D("hsW"," Smeared W   ", 40 ,  30., 150. );
     TRandom *gRan = new TRandom();

     int evtId0 = 0;
     std::vector<TLorentzVector> wjj ;
     std::vector<TLorentzVector> wqq ;
     int wid = 0 ;
     bool nlep = 0;
     for (int i=1; i<= tr->GetEntries() ; i++) {
         tr->GetEntry(i);
         if ( evtId == 1 ) evtId0 = evtId ;
         // check W mass
         bool realW = ( wid == 3 || wid == 7 ) ? true : false ;
         if ( evtId != evtId0 && wjj.size() == 2 && wqq.size()==2 && realW && nlep == 1 ) {
            TLorentzVector swP4( wjj[0].Px()+wjj[1].Px(),wjj[0].Py()+wjj[1].Py(),wjj[0].Pz()+wjj[1].Pz(),wjj[0].E()+wjj[1].E() );
            TLorentzVector gwP4( wqq[0].Px()+wqq[1].Px(),wqq[0].Py()+wqq[1].Py(),wqq[0].Pz()+wqq[1].Pz(),wqq[0].E()+wqq[1].E() );
            hsW->Fill( swP4.M() );
            hgW->Fill( gwP4.M() );
            if ( hSmearW != NULL ) hSmearW->Fill( swP4.M() ); 
            //cout<<" record "<< evtId0 <<endl;
         }
         if ( evtId != evtId0 ) {
            wjj.clear();
            wqq.clear(); 
            wid = 0;
            evtId0 = evtId ;
            nlep = 0;
         }
         if ( abs(objId) ==  13 ) nlep++ ;
         if ( abs(objId) > 4 || objId == 0 ) continue;
         //cout<<"objId:"<< objId <<" evtId:"<<evtId<<" wjj size="<< wjj.size()<<" wqq size="<< wqq.size()<<" wid:"<<wid<< endl;

         TLorentzVector qP4( px, py, pz, E );
         double gval = gRan->Gaus(0,1) ;

         double dEt = sgm*sqrt( (5.6*5.6) +  ( 1.25*1.25*qP4.Et() ) +  ( 0.033*0.033*qP4.Et()*qP4.Et() ) ); 
         double dE = gval*dEt*( qP4.P()/ qP4.Pt() ) ;
         double Scale = ( E + dE ) / E ;
         TLorentzVector jP4( Scale*px, Scale*py, Scale*pz, Scale*E );
         double ResE = 1 -(jP4.E()/E)  ;

         if ( fabs(jP4.Eta()) > 2.7 ) continue ; 
         if ( jP4.Et() < 30 ) continue ; 
         wjj.push_back( jP4 );
         wqq.push_back( qP4 );
         wid += abs( objId ) ;

         hGen->Fill( E );
         hSmr->Fill( jP4.E() );
         hRes->Fill( ResE );
     }

     if ( hSmearW == NULL ) {

        TCanvas* c3 = new TCanvas("c3","", 900, 800);
	gStyle->SetStatFontSize(0.04);
	gStyle->SetStatW(0.25);
	gStyle->SetOptStat("neim");
	gStyle->SetOptFit(111);
	c3->SetGrid();
	c3->SetFillColor(10);
	c3->Divide(2,2);
	c3->cd(1);
	gPad->SetGridx();
	hGen->Draw(); 
	c3->Update();

	c3->cd(2);
	gPad->SetGridx();
	hSmr->SetLineColor(2);
	hSmr->Draw();
	c3->Update();

	c3->cd(3);
	hRes->Draw();
	hRes->Fit( "gaus", "R","sames", -1, 1);
	c3->Update();

	c3->cd(4);
	hgW->Draw();
	c3->Update();
	hsW->SetLineColor(2);
	gStyle->SetStatY(0.75);
	gStyle->SetStatTextColor(2);
	hsW->Draw("sames");
	c3->Update();

	c3->Print(hfolder+"gRandom1.gif");

	TCanvas* c4 = new TCanvas("c4","", 700, 800);
	c4->SetFillColor(10);
	c4->Divide(1,2);
	c4->cd(1);
	gPad->SetGridx();
	hgW->Draw(); 
	c4->Update();

	c4->cd(2);
	gPad->SetGridx();
	gStyle->SetStatY(0.95);
	hsW->Draw(); 
	c4->Update();
	c4->Print(hfolder+"wmass1.gif");
     }
     
}

void JES::SmearAndMatch( TString mName ){

     TH1D* hsm0 = new TH1D("hsm0"," W Mass w/ 0.7  ", 40, 0, 200 );
     TH1D* hsm1 = new TH1D("hsm1"," W Mass w/ 0.9  ", 40, 0, 200 );
     TH1D* hsm2 = new TH1D("hsm2"," W Mass w/ 1.1  ", 40, 0, 200 );
     TH1D* hsm3 = new TH1D("hsm3"," W Mass w/ 1.3  ", 40, 0, 200 );
     TH1D* hsm4 = new TH1D("hsm4"," W Mass w/ 1.5  ", 40, 0, 200 );
     TH1D* ttW2 = new TH1D("ttW2"," W Mass from Tt MC matched  ", 40, 0, 200 );

     Smearing( mName, 0.7, hsm0 );
     Smearing( mName, 0.9, hsm1 );
     Smearing( mName, 1.1, hsm2 );
     Smearing( mName, 1.3, hsm3 );
     Smearing( mName, 1.5, hsm4 );
     //WMassSpectrum( mName, ttW0, ttW1, ttW2, 1. );
     Matching( mName, 1., ttW2 );

     double s0 = hsm0->Integral();
     double s1 = hsm1->Integral();
     double s2 = hsm2->Integral();
     double s3 = hsm3->Integral();
     double s4 = hsm4->Integral();
     double n1 = ttW2->Integral();

     hsm0->Scale(2500/s0) ;
     hsm1->Scale(2500/s1) ;
     hsm2->Scale(2500/s2) ;
     hsm3->Scale(2500/s3) ;
     hsm4->Scale(2500/s4) ;
     ttW2->Scale(2500/n1) ;

     double mX[5] ={0.7, 0.9, 1.1, 1.3, 1.5};
     double eX[5] ={0};
     double mW[5] ={0};
     double eW[5] ={0};

     TCanvas* c5 = new TCanvas("c5","", 900, 800);
     gStyle->SetStatW(0.25);
     gStyle->SetOptStat("neim");
     gStyle->SetOptFit(111);
     c5->SetGrid();
     c5->SetFillColor(10);

     hsm0->Draw();
     c5->Update();

     TF1* fGs = new TF1("fGs", MassFitFunction::fitGS , 30, 150, 3);
     fGs->SetParLimits(0, 5., hsm0->Integral() );
     fGs->SetParameter(1, 80.4 );
     fGs->SetParLimits(2, 5., 100.);
     fGs->SetLineColor(5);
     hsm0->SetLineColor(5);
     hsm0->Fit( fGs, "R","sames", 30, 150 );
     mW[0] = fGs->GetParameter(1) ;
     eW[0] = fGs->GetParameter(2) ;
     c5->Update();

     fGs->SetParLimits(0, 5., hsm1->Integral() );
     fGs->SetParameter(1, 80.4 );
     fGs->SetParLimits(2, 5., 100.);
     fGs->SetLineColor(2);
     hsm1->SetLineColor(2);
     hsm1->Fit( fGs, "R","sames", 30, 150 );
     mW[1] = fGs->GetParameter(1) ;
     eW[1] = fGs->GetParameter(2) ;
     c5->Update();

     fGs->SetParLimits(0, 5., hsm2->Integral() );
     fGs->SetParameter(1, 80.4 );
     fGs->SetParLimits(2, 5., 100.);
     fGs->SetLineColor(4);
     hsm2->SetLineColor(4);
     hsm2->Fit( fGs, "R","sames", 30, 150 );
     mW[2] = fGs->GetParameter(1) ;
     eW[2] = fGs->GetParameter(2) ;
     c5->Update();

     fGs->SetParLimits(0, 5., hsm3->Integral() );
     fGs->SetParameter(1, 80.4 );
     fGs->SetParLimits(2, 5., 100.);
     fGs->SetLineColor(6);
     hsm3->SetLineColor(6);
     hsm3->Fit( fGs, "R","sames", 30, 150 );
     mW[3] = fGs->GetParameter(1) ;
     eW[3] = fGs->GetParameter(2) ;
     c5->Update();

     fGs->SetParLimits(0, 5., hsm4->Integral() );
     fGs->SetParameter(1, 80.4 );
     fGs->SetParLimits(2, 5., 100.);
     fGs->SetLineColor(8);
     hsm4->SetLineColor(8);
     hsm4->Fit( fGs, "R","sames", 30, 150 );
     mW[4] = fGs->GetParameter(1) ;
     eW[4] = fGs->GetParameter(2) ;
     c5->Update();

     fGs->SetParLimits(0, 5., ttW2->Integral() );
     fGs->SetParameter(1, 80.4 );
     fGs->SetParLimits(2, 5., 100.);
     fGs->SetLineColor(1);
     ttW2->SetLineColor(1);
     ttW2->Fit( fGs, "R","sames", 30, 150 );
     c5->Update();

     c5->Print(hfolder+"SmearEffb.gif");
      
     TCanvas* c6 = new TCanvas("c6","", 800, 600);
     c6->SetFillColor(10);
     gPad->SetGridx();
     gPad->SetGridy();
     gStyle->SetOptStat(kTRUE);
     gStyle->SetOptFit(111);
     c6->cd();

     TGraphErrors* hSmEff = new TGraphErrors( 5, mX, mW, eX, eW );
     hSmEff->SetMarkerColor(4);
     hSmEff->SetMarkerStyle(20);
     hSmEff->SetLineWidth(2);
     hSmEff->SetLineColor(4);
     hSmEff->SetTitle(" Smearing Effect");
     hSmEff->Draw("AP");

     TF1* func0 = new TF1("func0", MassFitFunction::fitPoly , 0.9, 1.1, 2);
     hSmEff->Fit( func0, "RQ","sames", mX[0]-0.1, mX[4]+0.1 );

     c6->Update();
     c6->Print(hfolder+"SmearEff1.gif");
}

