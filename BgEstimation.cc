
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
// index : 0->all , 1->nonQCD ,  2->QCD     => using EvtSelector
// index :          3->nonQCD ,  4->QCD     => using QCDSelector
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

      if ( index == 3 ) objInfo->ResetBGCuts(0, 50 );
      if ( index == 3 ) objInfo->ResetBGCuts(1, 30 );
      if ( index == 4 ) objInfo->ResetBGCuts(0, 35 );
      if ( index == 4 ) objInfo->ResetBGCuts(1, 20 );
      if ( normMC ) objInfo->Reset(1,true) ;
      int scaleMode =  ( normMC ) ? 2 : 0  ;

      double theScale = 1 ;

      objInfo->Reset(0, inclX );
      objInfo->Reset(0, njx) ;
      hLepM2* nc_xj = new hLepM2( nameStr.str().substr(3,2) ) ;
      for ( int i = 0; i < fNames.size(); i++) {
          if ( index == 1 && i == 2 ) continue ; 
          if ( index == 2 && i != 2 ) continue ; 
          if ( fNames[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNames[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNames[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNames[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" )  ;
          if ( fNames[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNames[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( index <= 2 )  objInfo->EvtSelector( fNames[i], nc_xj, smearing, theScale );
          if ( index == 3 )  objInfo->QCDSelector( fNames[i], nc_xj, smearing, theScale, false, scaleMode );
          if ( index == 4 )  objInfo->QCDSelector( fNames[i], nc_xj, smearing, theScale, true,  scaleMode );
          
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
      for ( int i = 0; i < fNames.size(); i++) {
          if ( index == 1 && i == 2 ) continue ; 
          if ( index == 2 && i != 2 ) continue ; 
          if ( fNames[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNames[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNames[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNames[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" )  ;
          if ( fNames[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNames[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
          if ( index <= 2 )  objInfo->EvtSelector( fNames[i], nc_yj, smearing, theScale );
          if ( index == 3 )  objInfo->QCDSelector( fNames[i], nc_yj, smearing, theScale, false, scaleMode );
          if ( index == 4 )  objInfo->QCDSelector( fNames[i], nc_yj, smearing, theScale, true,  scaleMode );
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

// use 2 file lists which already have differnt jet multiplicity
vector<double> BgEstimation::RatioXY( vector<string>& fNameX, vector<string>& fNameY, int index, bool includeTt, bool doPlot ){

      cout<<" =====  Measuring Ratio from 2 Sets of files ===== "<<endl;
      double theScale = 1 ;
      //objInfo->Reset(0,true) ;
      string chtag = "R42" ; 

      objInfo->Reset(0, 4) ;
      hLepM2* nc_xj = new hLepM2( "Nxj", 20. ) ;
      for ( int i = 0; i < fNameX.size(); i++) {
          if ( i > 3 ) continue ; 
          if ( !includeTt && i == 0 && index <  0 ) continue; 
          if ( index > -1 && i != index ) continue ; 
          if ( fNameX[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNameX[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNameX[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNameX[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" ) ;
          if ( fNameX[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNameX[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
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
      for ( int i = 0; i < fNameY.size(); i++) {
          if ( i > 3 ) continue ; 
          if ( !includeTt && i == 0 ) continue; 
          if ( index > -1 && i != index ) continue ; 
          if ( fNameY[i].substr(0,2) == "tt" ) theScale = fitInput->NormalizeComponents( "tt" )  ;
          if ( fNameY[i].substr(0,2) == "wj" ) theScale = fitInput->NormalizeComponents( "wj" )  ;
          if ( fNameY[i].substr(0,2) == "zj" ) theScale = fitInput->NormalizeComponents( "zj" )  ;
          if ( fNameY[i].substr(0,2) == "qc" ) theScale = fitInput->NormalizeComponents( "qcd" )  ;
          if ( fNameY[i].substr(0,2) == "tq" ) theScale = fitInput->NormalizeComponents( "tq" )  ;
          if ( fNameY[i].substr(0,2) == "tw" ) theScale = fitInput->NormalizeComponents( "tw" )  ;
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

vector<double> BgEstimation::BgEstimate( vector<double>& R42, double nData_2j ){

     //vector<string> fNames ;
     //fitInput->GetParameters( "FakeData", &fNames );
     //vector<double> R42 = RatioXY( 4, 2, fNames, -1, false, false ) ;

     // 2. Calculate the 4Jets background and the uncertainty 
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

/*
void BgEstimation::CombinedRatio( int nx, int ny, bool inclX, bool inclY ){

      char rchar[5];
      sprintf( rchar, "R%d%d", nx, ny) ;
      string rTag = rchar ;

     vector<string> theFiles ;
     fitInput->GetParameters( "FakeData", &theFiles );

     bool doPlots = false ;
     vector< vector<double> > Rxy_V ;
     //vector<double> Rxy_tt  = RatioXY( nx, ny, theFiles,  0, true, doPlots,  inclX, inclY ) ;
     vector<double> Rxy_qcd = RatioXY( nx, ny, theFiles,  2, false,  doPlots,  inclX, inclY ) ;
     vector<double> Rxy_wj  = RatioXY( nx, ny, theFiles,  1, false,  doPlots,  inclX, inclY ) ;
     vector<double> Rxy_zj  = RatioXY( nx, ny, theFiles,  3, false,  doPlots,  inclX, inclY ) ;
     vector<double> Rxy_all = RatioXY( nx, ny, theFiles, -1, withTt, doPlots,  inclX, inclY ) ;
     //Rxy_V.push_back( Rxy_tt );
     Rxy_V.push_back( Rxy_qcd );
     Rxy_V.push_back( Rxy_wj );
     Rxy_V.push_back( Rxy_zj );
     Rxy_V.push_back( Rxy_all );

     TF1 *func = new TF1("func", MassFitFunction::fitPoly, 0, 200, 2);

     TCanvas* c1 = new TCanvas("c1","", 800, 600);
     c1->SetGrid();
     c1->SetFillColor(10);
     c1->SetFillColor(10);
     c1->cd();

     //TString nametag[5] = { "tt", "wj", "zj", "qcd", "all" } ;
     //int mcolor[5] = { 2, 3, 4, 5, 1 };
     TString nametag[4] = { "qcd", "wj", "zj", "all" } ;
     int mcolor[4]      = {     5,    3,    4,     1 };

     double aPt[4] = { 3, 23, 43, 63 } ;
     //double aPt[4] = { 10, 30, 50, 70 } ;
     double ePt[4] = { 0. } ;
     double aRxy[4] ={ 0. } ;
     double eRxy[4] ={ 0. } ;
     for ( size_t i=0; i< Rxy_V.size() ; i++) {
         for ( int j=1; j< 5; j++ ) {
             aPt[j-1]  = aPt[j-1] + 3 ;
             aRxy[j-1] = Rxy_V[i][j] ;
             eRxy[j-1] = Rxy_V[i][j+5] ;
         }

         TGraphErrors* hRatio = new TGraphErrors( 4, aPt, aRxy, ePt, eRxy );
	 //hRatio->SetMaximum( 16 );
	 hRatio->GetXaxis()->SetLimits(0,100);
	 hRatio->SetMarkerColor( mcolor[i] );
	 hRatio->SetMarkerStyle(20);
	 hRatio->SetMarkerSize(1.5);
	 hRatio->SetLineWidth(2);
	 hRatio->SetTitle("Ratio vs M2 Pt");

	 //hRatio->Draw("AP");
         //func->Draw("SAME");
	 if (i == 0) {
            hRatio->Draw("AP");
            func->FixParameter(0, Rxy_V[3][0] ) ;
            func->FixParameter(1,    0  ) ;
            func->SetLineColor( 2 );
            func->Draw("SAME");
         }
	 if (i  > 0) hRatio->Draw("P");
     
	 c1->Update();
         //TString plotname = hfolder + "hLepM2Info/Ratio_"+ nametag[i] +"."+plotType ;
         //c1->Print( plotname );
         //delete hRatio ;
         //delete c1;
     }
     TString plotname = hfolder + "hLepM2Info/" + rTag + "_bgCombined."+plotType ;
     c1->Print( plotname );
     delete func ;
     delete c1;
}
*/

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

  int nbin_ = 30 ;

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
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  hObjs* htt = new hObjs("tt", nbin_ ) ;
  objInfo->QCDSelector( fakeData[0], htt, false, scale0, doQCD );
  vector<TH2D*> h_tt ;
  htt->Fill2DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[1], hwj, false, scale1, doQCD );
  vector<TH2D*> h_wj ;
  hwj->Fill2DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd", nbin_ ) ;
  objInfo->QCDSelector( fakeData[2], hqcd, false, scale2, doQCD );
  vector<TH2D*> h_qcd ;
  hqcd->Fill2DVec( h_qcd );

  hObjs* hzj = new hObjs("zj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[3], hzj, false, scale3, doQCD );
  vector<TH2D*> h_zj ;
  hzj->Fill2DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin_ ) ;
  objInfo->QCDSelector( fakeData[4], htq, false, scale4, doQCD );
  vector<TH2D*> h_tq ;
  htq->Fill2DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin_ ) ;
  objInfo->QCDSelector( fakeData[5], htw, false, scale5, doQCD );
  vector<TH2D*> h_tw ;
  htw->Fill2DVec( h_tw );

  TH2D* hMC   = new TH2D("hMC", " MET(X) vs lep Mt(Y) ",   30, 0, 150,  nbin_, 0, 150 );
  TH2D* hMC_A = new TH2D("hMC_A", " MET(X) vs lep Mt(Y) ", 30, 0, 150,  nbin_, 0, 150 );
  TH2D* hData = new TH2D("hData", " MET(X) vs lep Mt(Y) ", 30, 0, 150,  nbin_, 0, 150 );

  hMC_A->Add( h_qcd[0] );
  hMC_A->Add( h_wj[0] );
  hMC_A->Add( h_zj[0] );
  hMC_A->Add( h_tt[0] );
  hMC_A->Add( h_tq[0] );
  hMC_A->Add( h_tw[0] );

  if ( doQCD ) { 
     hMC->Add( h_qcd[0] );
     hData->Add( h_dt[0],  1. );
     hData->Add( h_tt[0], -1. );
     hData->Add( h_wj[0], -1. );
     hData->Add( h_zj[0], -1. );
     hData->Add( h_tq[0], -1. );
     hData->Add( h_tw[0], -1. );
  }
  if ( !doQCD ) { 
     hMC->Add( h_wj[0] );
     hMC->Add( h_zj[0] );
     hMC->Add( h_tt[0] );
     hMC->Add( h_tq[0] );
     hMC->Add( h_tw[0] );
     hData->Add( h_dt[0],  1. );
     hData->Add( h_qcd[0], -1. );
  }
  cout<<" N of MC QCD = "<< h_qcd[0]->Integral() <<endl;
  cout<<" N of MC WJ  = "<< h_wj[0]->Integral() <<endl;
  cout<<" N of MC zJ  = "<< h_zj[0]->Integral() <<endl;
  cout<<" N of MC tt  = "<< h_tt[0]->Integral() <<endl;
  cout<<" N of MC tq  = "<< h_tq[0]->Integral() <<endl;
  cout<<" N of MC tw  = "<< h_tw[0]->Integral() <<endl;
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
  gStyle->SetOptStat("nieuo");
  gStyle->SetPalette(1);
  gStyle->SetStatX(0.85);

  TCanvas* c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  double xR = 120. ;
  double yR = 120. ;

  c1->cd(1);
  gStyle->SetNumberContours(10);
  h_dt[0]->SetName("Data");
  h_dt[0]->SetAxisRange(0,xR,"X") ;
  h_dt[0]->SetAxisRange(0,yR,"Y") ;
  h_dt[0]->Draw("COLZ");
  c1->Update();

  c1->cd(2);
  gStyle->SetNumberContours(4);
  h_qcd[0]->SetName("QCD");
  h_qcd[0]->SetAxisRange(0,xR,"X") ;
  h_qcd[0]->SetAxisRange(0,yR,"Y") ;
  h_qcd[0]->Draw("COLZ");
  c1->Update();

  c1->cd(3);
  gStyle->SetNumberContours(10);
  h_wj[0]->SetName("WJets");
  h_wj[0]->SetAxisRange(0,xR,"X") ;
  h_wj[0]->SetAxisRange(0,yR,"Y") ;
  h_wj[0]->Draw("COLZ");
  c1->Update();

  c1->cd(4);
  gStyle->SetNumberContours(10);
  hMC_A->SetName("MC");
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
      if ( nChi2 < 2 ) bestNorm = fitNorm ;
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

double BgEstimation::MeasureScale1D( string& DataName, vector<string>& fakeData, int hID, double MtCut, double METCut, bool doQCD, bool isVjNorm ){

  objInfo->ResetBGCuts(0, MtCut);
  objInfo->ResetBGCuts(1, METCut);
  objInfo->Reset(1, isVjNorm );  
  int scaleMode =  ( isVjNorm ) ? 1 : 0 ;

  int nbin_ = 30 ;

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
  double scale2 = fitInput->NormalizeComponents( "qcd" );
  double scale3 = fitInput->NormalizeComponents( "zj" );
  double scale4 = fitInput->NormalizeComponents( "tq" );
  double scale5 = fitInput->NormalizeComponents( "tw" );

  hObjs* htt = new hObjs("tt", nbin_ ) ;
  objInfo->QCDSelector( fakeData[0], htt, false, scale0, doQCD, scaleMode );
  vector<TH1D*> h_tt ;
  htt->Fill1DVec( h_tt );

  hObjs* hwj = new hObjs("wj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[1], hwj, false, scale1, doQCD, scaleMode );
  vector<TH1D*> h_wj ;
  hwj->Fill1DVec( h_wj );

  hObjs* hqcd = new hObjs("qcd", nbin_ ) ;
  objInfo->QCDSelector( fakeData[2], hqcd, false, scale2, doQCD, scaleMode );
  vector<TH1D*> h_qcd ;
  hqcd->Fill1DVec( h_qcd );

  hObjs* hzj = new hObjs("zj", nbin_ ) ;
  objInfo->QCDSelector( fakeData[3], hzj, false, scale3, doQCD, scaleMode );
  vector<TH1D*> h_zj ;
  hzj->Fill1DVec( h_zj );

  hObjs* htq = new hObjs("tq", nbin_ ) ;
  objInfo->QCDSelector( fakeData[4], htq, false, scale4, doQCD, scaleMode );
  vector<TH1D*> h_tq ;
  htq->Fill1DVec( h_tq );

  hObjs* htw = new hObjs("tw", nbin_ ) ;
  objInfo->QCDSelector( fakeData[5], htw, false, scale5, doQCD, scaleMode );
  vector<TH1D*> h_tw ;
  htw->Fill1DVec( h_tw );

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
  }
  if ( !doQCD ) { 
     hMC->Add( h_wj[hID] );
     hMC->Add( h_zj[hID] );
     hMC->Add( h_tt[hID] );
     hMC->Add( h_tq[hID] );
     hMC->Add( h_tw[hID] );
     hData->Add( h_dt[hID],  1. );
     hData->Add( h_qcd[hID], -1. );
  }

  cout<<" N of MC QCD = "<< h_qcd[hID]->Integral() <<endl;
  cout<<" N of MC WJ  = "<< h_wj[hID]->Integral() <<endl;
  cout<<" N of MC zJ  = "<< h_zj[hID]->Integral() <<endl;
  cout<<" N of MC tt  = "<< h_tt[hID]->Integral() <<endl;
  cout<<" N of MC tq  = "<< h_tq[hID]->Integral() <<endl;
  cout<<" N of MC tw  = "<< h_tw[hID]->Integral() <<endl;
  cout<<" N of Data = "<< h_dt[hID]->Integral() <<" -> "<< hData->Integral() <<endl ;

  cout<<" Event # normalization = "<< hData->Integral()/hMC->Integral() <<endl ;
  int wbin = 150/nBin ;
  int bin1 = ( doQCD ) ?         1    : (MtCut/wbin) + 1 ;
  int bin2 = ( doQCD ) ? (MtCut/wbin) :         100/wbin ;
  cout<<" bin1 = "<< bin1 <<"   bin2 = "<<bin2<<endl;
  char pNameChar[10];

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

   cout<<" delta Norm = "<< d_norm*0.34 <<"  Norm start from = "<< norm <<endl ;

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
      double nL = bestNorm - (bestNorm*0.1) ;
      double nH = bestNorm + (bestNorm*0.1) ;
      func1->SetParLimits(0, nL, nH );
      func1->SetParLimits(1, 0.00001, 100. );
      func1->SetParLimits(2, minX2*0.8 , minX2*1.2 );

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
      if ( nChi2 < 2 ) bestNorm = fitNorm ;

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
