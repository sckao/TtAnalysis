#include "SystInfo.h"
#include "WFormat.h"

SystInfo::SystInfo(){

  fitInput = new MassAnaInput();
  objInfo  = new ObjectInfo();
  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "PlotType", &plotType );
  fitInput->GetParameters( "JESType", &JESType );

}

SystInfo::~SystInfo(){

  delete fitInput ;
  delete objInfo  ;
}

void SystInfo::SystPlotter(){

  string SystId ;
  fitInput->GetParameters( "SystPlotSet", &SystId );
  TString subFolder = SystId ;

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hSys" );
  gSystem->cd( "hSys") ;
  gSystem->mkdir( subFolder );
  gSystem->cd( "../../" );

  string  flist[3]  ;
  TString legStr[3] ;
  double scale0 = 1  ;
  double scale1 = 1  ;
  double scale2 = 1  ;

  if ( SystId == "wjME" ) {
     string  flistA[3]  = { "wjMEu_387_0j", "wj0+", "wjMEd_387_0j" } ;
     TString legStrA[3] = { "Matching Threshold Up", "Nominal", "Matching Threshold Down" } ;
     scale0 = fitInput->NormalizeComponents( "wjMEu", "VJetSystematic.txt" )  ;
     scale1 = fitInput->NormalizeComponents( "wj" )  ;
     scale2 = fitInput->NormalizeComponents( "wjMEd", "VJetSystematic.txt" )  ;
     for (int i=0; i<3; i++) {
         flist[i]  = flistA[i] ;
         legStr[i] = legStrA[i] ;
     }
  }

  if ( SystId == "wjQ2" ) {
     string flistA[3]   = { "wjQ2u_387_0j", "wj0+", "wjQ2d_387_0j" } ;
     TString legStrA[3] = { "Q2 Up", "Nominal", "Q2 Down" } ;
     scale0 = fitInput->NormalizeComponents( "wjQ2u", "VJetSystematic.txt" )  ;
     scale1 = fitInput->NormalizeComponents( "wj" )  ;
     scale2 = fitInput->NormalizeComponents( "wjQ2d", "VJetSystematic.txt" )  ;
     for (int i=0; i<3; i++) {
         flist[i]  = flistA[i] ;
         legStr[i] = legStrA[i] ;
     }
  }

  if ( SystId == "ttR" ) {
     string flistA[3]   = { "ttRb_387_0j", "tt_387_0j", "ttRs_387_0j" } ;
     TString legStrA[3] = { "Larger ISR/FSR", "Nominal", "Smaller ISR/FSR" } ;
     scale0 = fitInput->NormalizeComponents( "ttRb", "TtSystematic.txt" )  ;
     scale1 = fitInput->NormalizeComponents( "tt" )  ;
     scale2 = fitInput->NormalizeComponents( "ttRs", "TtSystematic.txt" )  ;
     for (int i=0; i<3; i++) {
         flist[i]  = flistA[i] ;
         legStr[i] = legStrA[i] ;
     }
  }

  if ( SystId == "ttQ2" ) {
     string flistA[3]   = { "ttQ2u_387_0j", "tt_387_0j", "ttQ2d_387_0j" } ;
     TString legStrA[3] = { "Q2 up", "Nominal", "Q2 down" };
     scale0 = fitInput->NormalizeComponents( "ttQ2u", "TtSystematic.txt" )  ;
     scale1 = fitInput->NormalizeComponents( "tt" )  ;
     scale2 = fitInput->NormalizeComponents( "ttQ2d", "TtSystematic.txt" )  ;
     for (int i=0; i<3; i++) {
         flist[i]  = flistA[i] ;
         legStr[i] = legStrA[i] ;
     }
  }

  if ( SystId == "ttME" ) {
     string flistA[3]   = { "ttMEu_387_0j", "tt_387_0j", "ttMEd_387_0j" } ;
     TString legStrA[3] = { "Matching Threshold Up", "Nominal", "Matching Threshold Down" };
     scale0 = fitInput->NormalizeComponents( "ttMEu", "TtSystematic.txt" )  ;
     scale1 = fitInput->NormalizeComponents( "tt" )  ;
     scale2 = fitInput->NormalizeComponents( "ttMEd", "TtSystematic.txt" )  ;
     for (int i=0; i<3; i++) {
         flist[i]  = flistA[i] ;
         legStr[i] = legStrA[i] ;
     }
  }
  //double scale0 = norm / fitInput->TreeSize( flist[0] )  ;
  //double scale1 = norm / fitInput->TreeSize( flist[1] )  ;
  //double scale2 = norm / fitInput->TreeSize( flist[2] )  ;

  // scale-up sample
  hObjs* hsysu = new hObjs( "sysu" ) ;
  objInfo->EvtSelector( flist[0], hsysu, false, scale0 );
  vector<TH1D*> h_sysu ;
  hsysu->Fill1DVec( h_sysu );

  // nominal sample
  hObjs* hsys = new hObjs( "sys" ) ;
  objInfo->EvtSelector( flist[1], hsys, false, scale1 );
  vector<TH1D*> h_sys ;
  hsys->Fill1DVec( h_sys );

  // scale down sample
  hObjs* hsysd = new hObjs( "sysd" ) ;
  objInfo->EvtSelector( flist[2], hsysd, false, scale2 );
  vector<TH1D*> h_sysd ;
  hsysd->Fill1DVec( h_sysd );


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
      if ( i == 13 ) h_sys[i]->SetMinimum( 0.);

      h_sysu[i]->SetLineWidth(3);
      h_sysu[i]->SetLineColor(kRed);
      h_sys[i]->SetLineColor(1);
      h_sysd[i]->SetLineWidth(3);
      h_sysd[i]->SetLineColor(kBlue);

      // Tt jet multiplicity
      TCanvas* c0 = new TCanvas("c0","", 800, 700);
      c0->SetFillColor(10);
      c0->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c0->SetLogy();
      c0->cd();

      double tMax = 1.7* h_sys[i]->GetBinContent(h_sys[i]->GetMaximumBin()) ;
      if ( SystId.substr(0,2) == "tt" && (i==0 || i==7) ) tMax = 4.0*h_sys[i]->GetBinContent(h_sys[i]->GetMaximumBin()) ;
      if ( SystId.substr(0,2) == "wjQ2" && (i==13 || i== 15) ) tMax = 2.0*h_sys[i]->GetBinContent(h_sys[i]->GetMaximumBin()) ;
      h_sys[i]->SetMaximum( tMax );

      h_sys[i]->SetMarkerSize(1);
      h_sys[i]->SetMarkerStyle(21);
      h_sys[i]->Draw("P");
      c0->Update();
      h_sysu[i]->Draw("SAME");
      c0->Update();
      h_sysd[i]->Draw("SAME");
      c0->Update();

      double legXmin =  ( SystId.substr(2,4) == "Q2"  ) ? 0.58 : 0.52 ;
      TLegend *tleg = new TLegend( legXmin, .72, .9, .9 );
      tleg->AddEntry(h_sysu[i], legStr[0],   "L");
      tleg->AddEntry(h_sys[i],  legStr[1],   "P");
      tleg->AddEntry(h_sysd[i], legStr[2],   "L");
      tleg->Draw("same");
      c0->Update();

      TString plotname0 = hfolder + "hSys/"+subFolder+"/Sys_"+ hNames[i]+ "."+plotType ;
      c0->Print( plotname0 );

      delete c0;
      delete tleg;
  }

}

// JESType = 2(JES), 5(JER), 6(MET)
void SystInfo::JESPlotter( vector<string>& fakeData, double norm ){

  TString theFolder = hfolder ;
  gSystem->mkdir( theFolder );
  gSystem->cd( theFolder );
  gSystem->mkdir( "hJES" );
  gSystem->cd( "../" );

  int SystId = 1 ;
  fitInput->GetParameters( "SystPlotSet", &SystId );
  int k = SystId ;
  // only use the scale up value for different jes_type 
  int jes_type = JESType ;

  double scale[7] = { 1. } ; 
  scale[0] = ( norm == 0 ) ? fitInput->NormalizeComponents( "tt" ) : norm / fitInput->TreeSize( fakeData[0] );
  scale[1] = ( norm == 0 ) ? fitInput->NormalizeComponents( "wj" ) : norm / fitInput->TreeSize( fakeData[1] );
  scale[2] = ( norm == 0 ) ? fitInput->NormalizeComponents( "zj" ) : norm / fitInput->TreeSize( fakeData[2] );
  scale[3] = ( norm == 0 ) ? fitInput->NormalizeComponents( "tq" ) : norm / fitInput->TreeSize( fakeData[3] );
  scale[4] = ( norm == 0 ) ? fitInput->NormalizeComponents( "tw" ) : norm / fitInput->TreeSize( fakeData[4] );
  scale[5] = ( norm == 0 ) ? fitInput->NormalizeComponents( "ww" ) : norm / fitInput->TreeSize( fakeData[5] );
  scale[6] = ( norm == 0 ) ? fitInput->NormalizeComponents( "qcd" ): norm / fitInput->TreeSize( fakeData[6] );

  TString tName = "JES_" ;
  if ( jes_type == 5 ) tName = "JER_" ;
  if ( jes_type == 6 ) tName = "MET_" ;

  objInfo->Reset(1, jes_type );   // scale up JES, Unclustered Energy, scale-down JER

  hObjs* hjesu = new hObjs("jesu" ) ;
  objInfo->EvtSelector( fakeData[k], hjesu, false, scale[k] );
  vector<TH1D*> h_jesu ;
  hjesu->Fill1DVec( h_jesu );

  objInfo->Reset(1, 0 );   // no JES scale

  hObjs* hjes0 = new hObjs("jes0" ) ;
  objInfo->EvtSelector( fakeData[k], hjes0, false, scale[k] );
  vector<TH1D*> h_jes0 ;
  hjes0->Fill1DVec( h_jes0 );

  if ( jes_type == 2 || jes_type == 6 ) objInfo->Reset(1, jes_type+1 );   // scale down JES
  if ( jes_type == 5                  ) objInfo->Reset(1, jes_type-1 );   // scale down JER

  hObjs* hjesd = new hObjs("jesd" ) ;
  objInfo->EvtSelector( fakeData[k], hjesd, false, scale[k] );
  vector<TH1D*> h_jesd ;
  hjesd->Fill1DVec( h_jesd );

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
      if ( i == 13 ) h_jes0[i]->SetMinimum( 0.);

      h_jesu[i]->SetLineWidth(3);
      h_jesu[i]->SetLineColor(kRed);
      h_jes0[i]->SetLineColor(1);
      h_jesd[i]->SetLineWidth(3);
      h_jesd[i]->SetLineColor(kBlue);

      // Tt jet multiplicity
      TCanvas* c0 = new TCanvas("c0","", 800, 700);
      c0->SetFillColor(10);
      c0->SetFillColor(10);
      if ( i==0 || i==7 || i== 8 ) c0->SetLogy();
      c0->cd();

      double yMaxScale = (k == 0 && (i==0 || i==7 ) ) ? 3 : 1.5 ;
      double tMax = yMaxScale* h_jes0[i]->GetBinContent(h_jes0[i]->GetMaximumBin()) ;
      h_jes0[i]->SetMaximum( tMax );
      h_jes0[i]->SetMarkerSize(1);
      h_jes0[i]->SetMarkerStyle(21);
      h_jes0[i]->Draw("P");
      c0->Update();
      h_jesu[i]->Draw("SAME");
      c0->Update();
      h_jesd[i]->Draw("SAME");
      c0->Update();

      TLegend *tleg = new TLegend(.62, .7, .9, .9 );
      tleg->AddEntry(h_jesu[i], "Scale-up",   "L");
      tleg->AddEntry(h_jes0[i],  "Nominal",    "P");
      tleg->AddEntry(h_jesd[i], "Scale-down", "L");
      tleg->Draw("same");
      c0->Update();

      TString chName = fakeData[k].substr(0,2);
      TString plotname0 = hfolder + "hJES/" + chName +"_"+ tName + hNames[i]+ "."+plotType ;
      c0->Print( plotname0 );

      delete c0;
  }

}

