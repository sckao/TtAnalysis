#include "HadWMassFitter.h"

HadWMassFitter::HadWMassFitter( double massL, double massH ): eventId(0){

  fitInput = new MassAnaInput( "had", massL, massH );
  //fitTools = new MassAna( "had", massL, massH );

  mL = massL;
  mH = massH;

  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "bThreshold", &bTh);
  fitInput->GetParameters( "n_btag", &n_btag);

}

HadWMassFitter::~HadWMassFitter(){

  delete fitInput ;
  //delete fitTools ;

}

TMinuit minuit(6);

void HadWMassFitter::FitW( TLorentzVector jetv1, TLorentzVector jetv2, Double_t* pars, Double_t* errs, bool isJES ) {

   // Initialize TMinuit via generic fitter interface with a maximum of 6 params
   FitVectors fvs  ;
   fvs.jv1 = jetv1 ;
   fvs.jv2 = jetv2 ;
   minuit.SetObjectFit( &fvs );

   pars[0] = 0.;
   pars[1] = 0.;
   pars[2] = 0.;
   pars[3] = 0.;
   pars[4] = 1.;
   pars[5] = 1.;
   
   //TLorentzVector jv12 = jetv1 + jetv2 ;
   //cout<<" jet1 p4: ( "<< jetv1.Px()<<" , "<< jetv1.Py()<<" , "<<jetv1.Pz()<<" , "<< jetv1.Pt()<<" ) " <<endl;
   //cout<<" jet2 p4: ( "<< jetv2.Px()<<" , "<< jetv2.Py()<<" , "<<jetv2.Pz()<<" , "<< jetv1.Pt()<<" ) " <<endl;
   //cout<<" Wmass 1 = "<< jv12.M() <<endl;
   minuit.SetFCN( WMFCN );

   // correction for PT
   sigma1 = 0.7*sqrt( (5.6*5.6) +  ( 1.25*1.25*jetv1.Pt() ) +  ( 0.033*0.033*jetv1.Pt()*jetv1.Pt() ) );
   sigma2 = 0.7*sqrt( (5.6*5.6) +  ( 1.25*1.25*jetv2.Pt() ) +  ( 0.033*0.033*jetv2.Pt()*jetv2.Pt() ) );
   minuit.DefineParameter(0, "RES1",   0,   0.5,   -1*sigma1,  sigma1 );
   minuit.DefineParameter(1, "RES2",   0,   0.5,   -1*sigma2,  sigma2 );

   zigma1 = sigma1*fabs( jetv1.Pz() / jetv1.Pt() );
   zigma2 = sigma2*fabs( jetv2.Pz() / jetv2.Pt() );
   // correction for PZ
   minuit.DefineParameter(2, "RESz1",  0,   0.5,   -1*zigma1,  zigma1 );
   minuit.DefineParameter(3, "RESz2",  0,   0.5,   -1*zigma2,  zigma2 );

   // Re-Scale Jet Energy
   // JES1 for |eta| < 1.3 , JES2 for |eta| >= 1.3
   minuit.DefineParameter(4, "JES1",   1.,   0.1,   0.90,  1.1 );
   minuit.DefineParameter(5, "JES2",   1.,   0.1,   0.90,  1.1 );

   // -1  quiet (also suppresse all warnings), 0: normal, 1: verbose
   minuit.SetPrintLevel(-1);
   //	1 for Chi2 , 0.5 for negative log likelihood
   minuit.SetErrorDef(0.5);
   
   if ( isJES ) {
      minuit.FixParameter(0);
      minuit.FixParameter(1);
      minuit.FixParameter(2);
      minuit.FixParameter(3);
   } else {
      minuit.FixParameter(4);
      minuit.FixParameter(5);
   }
   /*
   minuit.Migrad();
   minuit.Release(2);
   */
   minuit.Migrad(); 

   minuit.GetParameter(0, pars[0], errs[0]);
   minuit.GetParameter(1, pars[1], errs[1]);
   minuit.GetParameter(2, pars[2], errs[2]);
   minuit.GetParameter(3, pars[3], errs[3]);
   minuit.GetParameter(4, pars[4], errs[4]);
   minuit.GetParameter(5, pars[5], errs[5]);

}


void HadWMassFitter::WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag)
{

    FitVectors * fVec = (FitVectors*) minuit.GetObjectFit();
    TLorentzVector jv1 = fVec->jv1 ;
    TLorentzVector jv2 = fVec->jv2 ;

    double JES1 =  ( fabs( jv1.Eta() ) < 1.3 ) ? par[4] : par[5] ;
    double JES2 =  ( fabs( jv2.Eta() ) < 1.3 ) ? par[4] : par[5] ;

    double cos1 = jv1.Px()/jv1.Pt() ;
    double cos2 = jv2.Px()/jv2.Pt() ;
    double sin1 = jv1.Py()/jv1.Pt() ;
    double sin2 = jv2.Py()/jv2.Pt() ;

    double px1 = (jv1.Px()+par[0]*cos1)*JES1 ;
    double py1 = (jv1.Py()+par[0]*sin1)*JES1 ;
    double px2 = (jv2.Px()+par[1]*cos2)*JES2 ;
    double py2 = (jv2.Py()+par[1]*sin2)*JES2 ;

    double pz1 = (jv1.Pz()+par[2])*JES1 ;
    double pz2 = (jv2.Pz()+par[3])*JES2 ;

    double jE1 = sqrt( px1*px1 + py1*py1 + pz1*pz1 ) ;
    double jE2 = sqrt( px2*px2 + py2*py2 + pz2*pz2 ) ;

    TLorentzVector wp4( px1+px2, py1+py2, pz1+pz2, jE1+jE2 );
    Double_t dwmass = wp4.M() - 80.4 ;

    wChi2 =  dwmass*dwmass ;
}

double HadWMassFitter::ReFitWMass( TLorentzVector v1, TLorentzVector v2, Double_t *par, Double_t* ProbW ){

   double JES1 =  ( fabs( v1.Eta() ) < 1.3 ) ? par[4] : par[5] ;
   double JES2 =  ( fabs( v2.Eta() ) < 1.3 ) ? par[4] : par[5] ;

   double cos1 = v1.Px()/v1.Pt() ;
   double cos2 = v2.Px()/v2.Pt() ;
   double sin1 = v1.Py()/v1.Pt() ;
   double sin2 = v2.Py()/v2.Pt() ;

   double px1 = (v1.Px()+par[0]*cos1)*JES1 ;
   double py1 = (v1.Py()+par[0]*sin1)*JES1 ;
   double px2 = (v2.Px()+par[1]*cos2)*JES2 ;
   double py2 = (v2.Py()+par[1]*sin2)*JES2 ;

   double pz1 = (v1.Pz()+par[2])*JES1 ;
   double pz2 = (v2.Pz()+par[3])*JES2 ;

   double jE1 = sqrt( px1*px1 + py1*py1 + pz1*pz1 ) ;
   double jE2 = sqrt( px2*px2 + py2*py2 + pz2*pz2 ) ;

   TLorentzVector sjv1( px1, py1, pz1, jE1 );
   TLorentzVector sjv2( px2, py2, pz2, jE2 );

   TLorentzVector wp4( px1+px2, py1+py2, pz1+pz2, jE1+jE2 );

   /*
   double dP1 = jE1 - v1.E() ;
   double dP2 = jE2 - v2.E() ;
   double sW  = sqrt(dP1*dP1 + dP2*dP2 ) ;
   double chi2 = (wp4.M() - 80.4)*(wp4.M() - 80.4) / ( 2*dP1*dP1 + 2*dP2*dP2 );
   ProbW[0] = 1. / exp(chi2)*sW*sqrt(6.283) ;
   */
   double chi2 = fabs( wp4.M() - 80.4 );
   //ProbW[0] = ( chi2 < 1.1 ) ? 1. : 0. ;
   ProbW[0] = 1. ;

   return wp4.M() ;
}

double HadWMassFitter::ReFitHadTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par ) {

   double cos1 = v1.Px()/v1.Pt() ;
   double cos2 = v2.Px()/v2.Pt() ;
   double sin1 = v1.Py()/v1.Pt() ;
   double sin2 = v2.Py()/v2.Pt() ;

   double JES1 =  ( fabs( v1.Eta() ) < 1.3 ) ? par[4] : par[5] ;
   double JES2 =  ( fabs( v2.Eta() ) < 1.3 ) ? par[4] : par[5] ;
   double JES3 =  ( fabs( v3.Eta() ) < 1.3 ) ? par[4] : par[5] ;

   double px1 = (v1.Px()+par[0]*cos1)*JES1;
   double py1 = (v1.Py()+par[0]*sin1)*JES1;
   double px2 = (v2.Px()+par[1]*cos2)*JES2;
   double py2 = (v2.Py()+par[1]*sin2)*JES2;

   double pz1 = (v1.Pz()+par[2])*JES1;
   double pz2 = (v2.Pz()+par[3])*JES2;

   double px3 = v3.Px()*JES3;
   double py3 = v3.Py()*JES3;
   double pz3 = v3.Pz()*JES3;

   double jE1 = sqrt( px1*px1 + py1*py1 + pz1*pz1 ) ;
   double jE2 = sqrt( px2*px2 + py2*py2 + pz2*pz2 ) ;
   double jE3 = sqrt( px3*px3 + py3*py3 + pz3*pz3 ) ;

   TLorentzVector tp4( px1+px2+px3, py1+py2+py3, pz1+pz2+pz3, jE1+jE2+jE3 );

   //cout<<" paras = "<<par[0]<<" , "<<par[1]<<" , "<<par[2]<<" , "<<par[3]<<" , "<<par[4]<<" , "<<par[5]<<endl;
   return tp4.M() ;
}

// use v1: bjet, v2: neutrino, v3: muon
double HadWMassFitter::ReFitLepTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par ) {

   double JES1 =  ( fabs( v1.Eta() ) < 1.3 ) ? par[4] : par[5] ;

   double px1 = v1.Px()*JES1 ;
   double py1 = v1.Py()*JES1 ;
   double pz1 = v1.Pz()*JES1 ;
   double jE1 = v1.E()*JES1 ;

   TLorentzVector tp4( px1+v2.Px()+v3.Px(), py1+v2.Py()+v3.Py(), pz1+v2.Pz()+v3.Pz(), jE1+v2.E()+v3.E() );

   //cout<<" paras = "<<par[0]<<" , "<<par[1]<<" , "<<par[2]<<" , "<<par[3]<<" , "<<par[4]<<" , "<<par[5]<<endl;
   return tp4.M() ;
}

Double_t HadWMassFitter::jetAngle( TLorentzVector v1, TLorentzVector v2 ){
 
   double v1v2 = v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz() ;
   double v1l  = sqrt( v1.Px()*v1.Px() + v1.Py()*v1.Py() + v1.Pz()*v1.Pz() );
   double v2l  = sqrt( v2.Px()*v2.Px() + v2.Py()*v2.Py() + v2.Pz()*v2.Pz() );
   double cosA = v1v2 / ( v1l*v2l ) ;
   double angle = acos( cosA ) ;
   return angle ;
}

bool HadWMassFitter::HadPermutation( int i1, int i2, int i3, int evtId ) {

   bool duplicated = false ;

   if ( evtId == eventId ) {
      int samej = 0;
      for ( size_t k = 0; k < jd1.size() ; k++ ) {
          if ( i3 == jd3[k] ) samej++ ;
          if ( i1 == jd1[k] && i2 == jd2[k] ) samej++ ;
          if ( i1 == jd2[k] && i2 == jd1[k] ) samej++ ;

          if ( samej == 2 ) break;
      }
      if ( samej < 2 ) {
         jd1.push_back( i1 );
	 jd2.push_back( i2 );
	 jd3.push_back( i3 );
         samej = 0 ;
      }
      if ( samej ==2 ) {
          duplicated = true ;
          samej = 0 ;
      }
   } 

   if ( evtId != eventId  ) {
      eventId = evtId ;
      jd1.clear();
      jd2.clear();
      jd3.clear();
      jd1.push_back( i1 );
      jd2.push_back( i2 );
      jd3.push_back( i3 );
      duplicated = false ;
   }

   return duplicated ;
}

// jet re-selection, added more pt constrains
bool HadWMassFitter::JetPtFilter( vector<TLorentzVector> jpv ) {

     bool pass = true ;

     for( size_t i = 0; i < jpv.size(); i++ ) {
        if ( i >= 4 ) break ;
        if ( fabs( jpv[i].Eta() ) > 2.4 ) pass = false ;
        if ( jpv[i].Pt() < 30 ) pass = false;
     }
     return pass ;
}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void HadWMassFitter::ReFitSolution( string fileName, TH2D* hM2M3, TH2D* hM3M3, int type, bool isMCMatched ) {

  // get files and trees
  TFile*  file = NULL ;
  TTree* tr1;
  if ( isMCMatched ) {
     tr1 = fitInput->GetTree( fileName, "mcmTt",   file );
  } else {
     tr1 = fitInput->GetTree( fileName, "solTt",   file );
  }
  TTree* tr2 = fitInput->GetTree( fileName, "selJet",  file );
  TTree* tr3 = fitInput->GetTree( fileName, "solNeu",  file );
  TTree* tr4 = fitInput->GetTree( fileName, "selMu",   file );

  bool isJES = ( type == 2 ) ? true : false ;

  // event solution tree
  double wm, tm , prob;
  int evtId1, jid1[5] ;
  tr1->SetBranchAddress( "evtId" ,&evtId1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );
  tr1->SetBranchAddress( "hadWM" ,&wm   );
  tr1->SetBranchAddress( "hadTM" ,&tm   );
  tr1->SetBranchAddress( "probW" ,&prob   );

  int neu_str = 0 ;
  int jet_str = 0 ;
  int mu_str  = 0 ;
  int evtId0  = 0 ;
  vector<double> bCut ;
  vector<TLorentzVector> jvs;
  vector<TLorentzVector> nvs;
  vector<TLorentzVector> mvs;
  bool nextEvt[3] = {false, false, false} ;

  double norm = 0;
  std::vector<double> hTMlist ;
  std::vector<double> lTMlist ;
  std::vector<double> WMlist ;
  std::vector<double> wProb ;

  int bCounter = 0;
  int EvtCounter = 0;
  for ( int j= 0 ; j< tr1->GetEntries() ; j++ ) {
      tr1->GetEntry(j);

     // cout<<" -------------------------------------------- "<<endl;
     // cout<<" Event "<< evtId1 <<" $$ jid("<<jid1[0]<<","<<jid1[1]<<","<<jid1[2]<<","<<jid1[3]<<","<<jid1[4]<<")"<<endl;
     //
      // reset all containers at the beginning of the event
      if ( evtId0 != evtId1 ) {
         // record
         //cout <<" N of Top Solutions = "<< hTMlist.size() <<" Normal = "<<norm<< endl;
         EvtCounter++;
         for (size_t i=0; i< WMlist.size(); i++) {
             if ( norm == 0 ) continue;
             hM2M3->Fill( hTMlist[i], WMlist[i], wProb[i]/norm );
             hM3M3->Fill( hTMlist[i], lTMlist[i], wProb[i]/norm );
         }
         int bnu =0 ;
         for (size_t k=0; k< bCut.size(); k++) {
             if (bCut[k] > bTh) bnu++ ;
         }
         if ( bnu > 0 ) bCounter++;
         evtId0 = evtId1 ;
         nextEvt[0] = true ;
         nextEvt[1] = true ;
         nextEvt[2] = true ;
	 hTMlist.clear();
	 lTMlist.clear();
	 WMlist.clear();
	 wProb.clear();
	 jvs.clear();
	 nvs.clear();
	 mvs.clear();
         bCut.clear();
         norm = 0;

         // accumulate objects in the event
         jvs = fitInput->GetJets( tr2, evtId1, bCut, jet_str, nextEvt[0] );
         nvs = fitInput->GetNeutrinos( tr3, evtId1, neu_str, nextEvt[1] );
         mvs = fitInput->GetMuons( tr4, evtId1, mu_str, nextEvt[2] );

      }
      // exclude the 5th jet
      if ( jid1[0] == 4 || jid1[1] == 4 || jid1[2] == 4 || jid1[3] == 4 ) continue;
      // jet pt reselection
      if ( !JetPtFilter( jvs ) ) continue;

      //cout<<"  jvs size = "<< jvs.size() <<"  nvs size = "<< nvs.size() <<"  mvs size = "<< mvs.size() <<endl;
      // number of b tagged
      int NofB = 0;
      for (size_t k=0; k< bCut.size(); k++) {
          if (bCut[k] > bTh) NofB++ ;
      }
      //cout <<" N of B ="<< NofB <<endl;
      // btagger switch
      if ( NofB == 0  && n_btag >= 0 ) continue;

      // b Cut  -- need more studies
      double probb = 0 ;
      if ( n_btag ==2 && NofB >= 2 && bCut[ jid1[2] ] >= bTh && bCut[ jid1[3]] >= bTh ) probb = 1;
      if ( n_btag ==0 && NofB >= 2 && bCut[ jid1[2] ] >= bTh && bCut[ jid1[3]] >= bTh ) probb = 1;
      if ( n_btag < 2 && NofB == 1 && bCut[ jid1[2] ] >= bTh  ) probb = 1;
      if ( n_btag < 2 && NofB == 1 && bCut[ jid1[3] ] >= bTh  ) probb = 1;
      if ( n_btag >=0 && NofB == 0  ) probb = 0;

      // refit the W mass
      Double_t para[6] = {0,0,0,0,1,1};
      Double_t errs[6];
      Double_t probW[1] ;
      if ( type > 0 )  FitW( jvs[ jid1[0] ] , jvs[ jid1[1] ], para, errs, isJES );
      double wMass = ReFitWMass(   jvs[ jid1[0] ] , jvs[ jid1[1] ], para, probW );
      double hadTMass = ReFitHadTopMass( jvs[ jid1[0] ] , jvs[ jid1[1] ], jvs[ jid1[2] ],  para);
      double lepTMass = ReFitLepTopMass( jvs[ jid1[3] ] , nvs[ jid1[4]-1 ], mvs[0],  para);

      // determine the weighting of each entry
      double aProb = ( type > 0 ) ? probW[0] : 1 ;
      if ( n_btag >= 0 ) aProb = aProb*probb ;

      // record the result and normalization
      hTMlist.push_back( hadTMass );
      lTMlist.push_back( lepTMass );
      WMlist.push_back( wMass );
      wProb.push_back( aProb );
      norm += aProb;   
  }
  cout<<" bCount: "<<bCounter<<" / EvtCount:"<< EvtCounter  <<endl;
  if ( file != NULL ) file->Close() ;
}

// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void HadWMassFitter::TWFitter( string channelName, int type, TString DrawOpt, bool isMCMatched ){
 
  string fileName;
  string channelType = channelName.substr( 0, 2 )  ;

  vector<string> chlist;
  fitInput->GetParameters( "channel", &chlist );

  double scale = 1. ;
  if ( channelType == "tt" ) {
     vector<string> mlist;
     fitInput->GetParameters( "Signals", &mlist );
     for ( size_t j = 0; j< mlist.size(); j++){
         if ( channelName == mlist[j] ) fileName =  mlist[j];
     }
     scale = fitInput->NormalizeComponents( chlist[0] );
     
  } else { 
     vector<string> bglist;
     fitInput->GetParameters( "Backgrounds", &bglist );
     for ( size_t j = 0; j< bglist.size(); j++) {
         if ( channelName == bglist[j] ) {
            fileName = bglist[j] ;
            scale = fitInput->NormalizeComponents( chlist[j+1] );
         }
     }
  }

  TH2D* hM2M3 = new TH2D("hM2M3"," M3-M2 befor JES/JEC ", 15, 50, 350, 15, 0, 300 );
  TH2D* hM3M3 = new TH2D("hM3M3"," hadM3-LepM3 befor JES/JEC ", 15, 50, 350, 15, 50, 350 );

  string flagJES ;
  if (type == 0) flagJES = "_"  ;
  if (type == 1) flagJES = "JEC"  ;
  if (type == 2) flagJES = "JES"  ;

  ReFitSolution( fileName, hM2M3, hM3M3, type, isMCMatched );

  hM2M3->Scale( scale );
  hM3M3->Scale( scale );

  gStyle->SetOptStat("neirom");
  gStyle->SetStatY(0.95);
  gStyle->SetStatTextColor(1);
  gStyle->SetPalette(1);

  // plot  M2 vs M3 without JES
  c1 = new TCanvas("c1","", 800, 600);
  c1->SetGrid();
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1);
  gStyle->SetStatX(0.30);
  gStyle->SetNumberContours(5);
  hM2M3->Draw( DrawOpt );
  c1->Update();

  c1->cd(2);
  gStyle->SetStatX(0.90);
  TH1D* hM2M3_py = hM2M3->ProjectionY("hM2M3_py");
  hM2M3_py->Draw();
  c1->Update();

  c1->cd(3);
  TH1D* hM2M3_px = hM2M3->ProjectionX("hM2M3_px");
  hM2M3_px->Draw();
  c1->Update();

  TString plotname = hfolder + fileName + flagJES + "M2M3.gif" ;
  if ( isMCMatched ) plotname = hfolder + fileName + flagJES + "M2M3_MC.gif" ;
  c1->Print( plotname );

  // plot  M3(lep) vs M3(had) without JES
  c2 = new TCanvas("c2","", 800, 600);
  c2->SetGrid();
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->Divide(2,2);
  c2->cd(1);
  gStyle->SetStatX(0.90);
  gStyle->SetNumberContours(5);
  hM3M3->Draw( DrawOpt );
  c2->Update();

  c2->cd(2);
  TH1D* hM3M3_py = hM3M3->ProjectionY("hM3M3_py");
  hM3M3_py->Draw();
  c2->Update();

  c2->cd(3);
  TH1D* hM3M3_px = hM2M3->ProjectionX("hM3M3_px");
  hM3M3_px->Draw();
  c2->Update();

  TString plotname1 = hfolder + fileName + flagJES + "M3M3.gif" ;
  if ( isMCMatched ) plotname1 = hfolder + fileName + flagJES + "M3M3_MC.gif" ;
  c2->Print( plotname1 );

  delete c1;
  delete c2;

}

