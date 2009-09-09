#include "MassFitInput.h"

MassFitInput::MassFitInput( TString channel, int NBTag ){
 
  hname = channel+"_selTmass_pt0";
  sname = "mcTtMass0";
  ch_name = channel;
  n_btag = NBTag;
  luminosity = 100;
  N_tt = 12000;
  N_wj = 108980;
  N_stt = 45018;
  N_stw = 25935;
  N_qcd = 2841112;

}

MassFitInput::~MassFitInput(){
 
  delete ttmass;

}

void MassFitInput::Initialize( TString* hfolder ) {

   if ( n_btag == -1 ) *hfolder = "NoTag_"+ch_name+"/";
   if ( n_btag == 0 )  *hfolder = "AllTags0_"+ch_name+"/";
   if ( n_btag == 1 )  *hfolder = "AllTags1_"+ch_name+"/";
   if ( n_btag == 2 )  *hfolder = "AllTags2_"+ch_name+"_10/";

}

void MassFitInput::GetFileName( TString mName, int type ) {

  if ( n_btag == -1 ) subtag = "NB";
  if ( n_btag ==  0 ) subtag = "B0";
  if ( n_btag ==  1 ) subtag = "B1";
  if ( n_btag ==  2 ) subtag = "B2";

  // File names for fake data 
  if ( type == 0 ) {
     theSG  = "pseud"+mName+subtag  ;
     theBG1 = "pseud_WJets_"+subtag ;
     theBG2 = "pseud_STT_"+subtag    ;
     theBG3 = "pseud_STTW_"+subtag   ;
     theBG4 = "pseud_QCD_"+subtag    ;
  }
  // File names for templates 
  if ( type == 1 ) {
     theSG  = "TCHE5_"+mName+"_"+subtag ;
     theBG1 = "TCHE5_WJets_"+subtag     ;
     theBG2 = "TCHE5_STT_"+subtag       ;
     theBG3 = "TCHE5_STTW_"+subtag      ;
     theBG4 = "TCHE5_QCD_"+subtag       ;
  }

}

void MassFitInput::getFakeData( TString mName, int rbin, TH1D* ttadd, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, TH1D* dth4, TH1D* dth5 ){


  // get the file names of fake data
  GetFileName( mName, 0 );
  // get all fake data 
  getFakeData(dth0, theSG,  true, rbin, 1. );             // tt-signal
  getFakeData(dth1, theSG,  false, rbin, 1. );            // tt-wrong combinatorics
  if (dth2 != NULL ) getFakeData(dth2, theBG1, rbin );          // w+jets
  if (dth3 != NULL ) getFakeData(dth3, theBG2, rbin, 0.195 );   // single top t-channel
  if (dth4 != NULL ) getFakeData(dth4, theBG3, rbin, 0.183 );   // single top tW-channel
  if (dth5 != NULL ) getFakeData(dth5, theBG4, rbin, 2.01 );    // QCD

  // mix all fake samples
  ttadd->Add(dth0, 1);
  ttadd->Add(dth1, 1);
  if (dth2 != NULL ) ttadd->Add(dth2, 1);
  if (dth3 != NULL ) ttadd->Add(dth3, 1);
  if (dth4 != NULL ) ttadd->Add(dth4, 1);
  if (dth5 != NULL ) ttadd->Add(dth5, 1);

  // give different color for samples
  if (dth5 != NULL ) dth5->SetFillColor(5);  // QCD
  if (dth4 != NULL ) dth4->SetFillColor(6);  // single top tw channel
  if (dth3 != NULL ) dth3->SetFillColor(4);  // single top t channel
  if (dth2 != NULL ) dth2->SetFillColor(2);  // w+jets
  dth1->SetFillColor(7);
  dth0->SetFillColor(3);

  // stack them
  if (dth5 != NULL ) ttstk->Add( dth5 );
  if (dth4 != NULL ) ttstk->Add( dth4 );
  if (dth3 != NULL ) ttstk->Add( dth3 );
  if (dth2 != NULL ) ttstk->Add( dth2 );
  ttstk->Add( dth1 );
  ttstk->Add( dth0 );

}

// for backgrounds
void MassFitInput::getFakeData( TH1D* ttadd, TString thefileName, int rbin, double theScale ){

  int nbin = 480./rbin ;
  TFile *datafile0  = TFile::Open(thefileName+".root");

  TH2D* fdata  = (TH2D *) datafile0->Get("Tops/"+hname);
  TH1D* fdata_pj = fdata->ProjectionX("fdata_pj",1,480,"");
  fdata_pj->Rebin(rbin);

  ttadd->Add(fdata_pj, theScale);
  ttadd->SetBinContent(nbin, 0.);

  delete fdata_pj;
  delete fdata;

  datafile0->Close();

}
// for tt signal and backgrounds
void MassFitInput::getFakeData( TH1D* fkdat, TString thefileName, bool isSignal, int rbin, double theScale ){

  int nbin = 480./rbin ;
  TFile *datafile1 = TFile::Open(thefileName+".root");

  // all combinations
  TH2D* hdata  = (TH2D *) datafile1->Get("Tops/"+hname);
  TH1D* hdata_pj = hdata->ProjectionX("hdata_pj",1,480,"");
  hdata_pj->Rebin(rbin);

  // matched combinations
  TH2D* hMC   = (TH2D *) datafile1->Get("MObjs/"+sname);
  TH1D* hMC_px = hMC->ProjectionX("hMC_px",1,480,"");
  TH1D* hMC_py = hMC->ProjectionY("hMC_py",1,480,"");
  hMC_px->Rebin(rbin);
  hMC_py->Rebin(rbin);

  TH1D* sg1 = new TH1D("sg1","", nbin, 0., 480.);
  double theN = 0. ;
  for (int i=1; i< nbin; i++) {
      if ( ch_name == "Had" ) theN = hMC_py->GetBinContent(i);
      if ( ch_name == "Lep" ) theN = hMC_px->GetBinContent(i);
      sg1->SetBinContent(i,theN);
  }

  if ( !isSignal ) fkdat->Add(hdata_pj,sg1, theScale, -1.*theScale );
  if (  isSignal ) fkdat->Add(sg1, theScale );
  fkdat->SetBinContent(nbin, 0.);

  delete sg1;
  delete hMC_px;
  delete hMC_py;
  delete hMC;
  delete hdata_pj;
  delete hdata;

  datafile1->Close();

}


void MassFitInput::getSignal(TH1D* h_Sg, int rbin, TString mName ) {

  int nbin = 480./rbin ;

  GetFileName( mName, 1 );

  TFile* thefile  = TFile::Open( theSG+".root" );
  TH2D* ttMC  = (TH2D *) thefile->Get("MObjs/"+sname);

  TH1D* ttMC_px = ttMC->ProjectionX("ttMC_px",0,-1,"");
  TH1D* ttMC_py = ttMC->ProjectionY("ttMC_py",0,-1,"");
  ttMC_px->Rebin(rbin);
  ttMC_py->Rebin(rbin);
  //ttMC_px->Scale(0.428);
  //ttMC_py->Scale(0.428);

  double theN = 0. ;
  for (int i=1; i< nbin; i++) {
      if ( ch_name == "Had" ) theN = ttMC_py->GetBinContent(i);
      if ( ch_name == "Lep" ) theN = ttMC_px->GetBinContent(i);
      h_Sg->SetBinContent(i,theN);
  }

  NormalizeComponents( luminosity, N_tt, 1., h_Sg );

  delete ttMC_px;
  delete ttMC_py;
  delete ttMC;
  thefile->Close();
}

void MassFitInput::getBackground(TH1D* h_Bg, int type, int rbin, TString mName ){

  int nbin = 480./rbin ;

  GetFileName( mName, 1 );

  if (type == 1) file1  = TFile::Open( theSG+".root" );  // tt-wrong combinatorics
  if (type == 2) file2  = TFile::Open( theBG1+".root" ); // w+jets
  if (type == 3) file2  = TFile::Open( theBG2+".root" ); // single top t-channel
  if (type == 4) file2  = TFile::Open( theBG3+".root" ); // single top tW-channel
  if (type == 5) file2  = TFile::Open( theBG3+".root" ); // QCD

  if ( type == 1 ) ttmass = (TH2D *) file1->Get("Tops/"+hname);
  if ( type >= 2 ) ttmass = (TH2D *) file2->Get("Tops/"+hname);

  TH1D* ttmass_pj = ttmass->ProjectionX("ttmass_pj",1,480,"");
  ttmass_pj->Rebin(rbin);

  if( type == 1 ) {

    TH1D* sgl = new TH1D("sgl","", nbin, 0., 480.);
    TH2D* sgnl   = (TH2D *) file1->Get("MObjs/"+sname);
    TH1D* sgnl_px = sgnl->ProjectionX("sgnl_px",1,480,"");
    TH1D* sgnl_py = sgnl->ProjectionY("sgnl_py",1,480,"");
    sgnl_px->Rebin(rbin);
    sgnl_py->Rebin(rbin);
    double theN = 0. ;
    for (int i=1; i< nbin; i++) {
        if ( ch_name == "Had" ) theN = sgnl_py->GetBinContent(i);
        if ( ch_name == "Lep" ) theN = sgnl_px->GetBinContent(i);
        sgl->SetBinContent(i,theN);
    }

    ttmass_pj->Add(sgl, -1 );
    delete sgnl;
    delete sgnl_px;
    delete sgnl_py;
    delete sgl;
  }

  for (int i=1; i< nbin; i++) {
      double theN = ttmass_pj->GetBinContent(i);
      h_Bg->SetBinContent(i,theN);
  }

  // scale to 100 /pb
  if ( type == 1 ) NormalizeComponents( luminosity, N_tt,  1, h_Bg );  // ttbar
  if ( type == 2 ) NormalizeComponents( luminosity, N_wj, 2, h_Bg );  // wjets
  if ( type == 3 ) NormalizeComponents( luminosity, N_stt,  3, h_Bg );  // single top t
  if ( type == 4 ) NormalizeComponents( luminosity, N_stw,  4, h_Bg );  // single top tW
  if ( type == 5 ) NormalizeComponents( luminosity, N_qcd,5, h_Bg );  // QCD

  delete ttmass_pj;

  if (type == 1 ) file1->Close();
  if (type >= 2 ) file2->Close();

}
 
void MassFitInput::combineBG( TString mName, TH1D* allbg, int rbin ) {

  int nbin = 480./rbin ;
  TH1D* tt = new TH1D("tt","", nbin, 0., 480.);
  TH1D* wj = new TH1D("wj","", nbin, 0., 480.);
  TH1D* stt = new TH1D("stt","", nbin, 0., 480.);
  TH1D* stw = new TH1D("stw","", nbin, 0., 480.);

  getBackground( tt, 1, rbin, mName );
  getBackground( wj, 2, rbin, mName );
  getBackground( stt, 3, rbin, mName );
  getBackground( stw, 4, rbin, mName );
  allbg->Add(tt, 1.);
  allbg->Add(wj, 1.);
  allbg->Add(stt, 1.);
  allbg->Add(stw, 1.);

  delete tt;
  delete wj;
  delete stt;
  delete stw;
}
 
void MassFitInput::getData( TH1D* h_data, TString thefileName, int rbin ){

  int nbin = 480./rbin ;

  TFile *file0  = TFile::Open(thefileName+".root");
  TH2D* hdata  = (TH2D *) file0->Get("Tops/"+hname);

  TH1D* hdata_pj = hdata->ProjectionX("hdata_pj",1,480,"");
  hdata_pj->Rebin(rbin);

  double theN = 0. ;
  for (int i=1; i< nbin; i++) {
      theN  = hdata_pj->GetBinContent(i);
      h_data->SetBinContent(i,theN);
  }

  delete hdata_pj;
  delete hdata;
  file0->Close();

}

void MassFitInput::NormalizeComponents( double lumi, double nEvents, int theChannel, TH1D* tmp ){

  //**old for 100 /pb => ttbar signal = 8992 events ; wjets = 30393 ; QCD = 2841112
  //**new for 100 /pb => ttbar signal = 9067 events ; wjets = 34000 ; 
  //                     st_tw = 2730, st_t = 6360 
  // channel #   => ttbar signal =  1 ; wjets = 2 ; Signle Top t_ch = 3 ; Single Top tW_ch = 4 ; QCD = 5
  
  int idx = theChannel - 1;
  double xsec[5] = {   414,  40000,   130,    29, 121675 };  // unit:pb, for MuEnrichedQCD already applied filter efficiency
  double Eff[5]  = { 0.219, 0.0085, 0.195, 0.183, 0.2335 };  // HLT and Data Skim efficiency
  double nBase = xsec[idx]*Eff[idx];
 
  double Scal = (nBase*lumi) / nEvents ;
  tmp->Scale(Scal);

}

