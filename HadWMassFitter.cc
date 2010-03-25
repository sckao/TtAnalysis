#include "HadWMassFitter.h"
#include "WFormat.h"

HadWMassFitter::HadWMassFitter( double massL, double massH ) {

  fitInput = new MassAnaInput( "had", massL, massH );

  mL = massL;
  mH = massH;

  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "bThreshold", &bTh);
  fitInput->GetParameters( "n_btag", &n_btag);

}

HadWMassFitter::~HadWMassFitter(){

  delete fitInput ;

}

//FitVectors * fVec;
TMinuit minuit(6) ;

void HadWMassFitter::WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag) {
//static void WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag) {

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

void HadWMassFitter::FitW( TLorentzVector jetv1, TLorentzVector jetv2, Double_t* pars, Double_t* errs, bool isJES ) {
   
   // Initialize TMinuit via generic fitter interface with a maximum of 6 params
   FitVectors fvs  ;
   fvs.jv1 = jetv1 ;
   fvs.jv2 = jetv2 ;
   //fVec = &fvs ;

   //TMinuit minuit(6) ;

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

double HadWMassFitter::KinematicProb( int jid[] , vector<TLorentzVector> vlist, TLorentzVector mV4, TLorentzVector nV4, Double_t *par, Double_t* wtX ) {

   //cout<<" jid = ("<< jid[0]<<" , "<< jid[1]<<" , "<< jid[2]<<" , "<< jid[3]<<" ) "<<endl;
   //cout<<" vlist size = "<< vlist.size() <<endl;
   vector<TLorentzVector> nvlist = newVector( vlist[ jid[0] ], vlist[ jid[1] ], vlist[ jid[2] ], vlist[ jid[3] ] , par );

   // new hadronic W p4
   TLorentzVector nhWp4 = nvlist[0] +  nvlist[1] ;
   double sgm_nhW = Get2BodySigma( nvlist[0], nvlist[1] ) ;

   // new hadronic Top p4
   TLorentzVector nhTp4 = nvlist[0] +  nvlist[1] + nvlist[2] ;
   double sgm_nhT = Get3BodySigma( nvlist[0], nvlist[1], nvlist[2] ) ;

   // new leptonic Top p4, not used yet
   TLorentzVector nlTp4 = nvlist[3] +  mV4 + nV4 ;

   double nhW_X2 = ( nhWp4.M() -  80.4 )*( nhWp4.M() -  80.4 ) / ( sgm_nhW*sgm_nhW ) ;
   double nhT_X2 = ( nhTp4.M() - 171.2 )*( nhTp4.M() - 171.2 ) / ( sgm_nhT*sgm_nhT ) ;

   double nprob_hW = exp( -0.5*nhW_X2 ) ;
   double nprob_hT = exp( -0.5*nhT_X2 ) ;

   double prob =  nprob_hW*nprob_hT ;
   double combX =  sqrt( nhW_X2 + nhT_X2 );

   if ( wtX != NULL ) wtX[0] = sqrt( nhW_X2 );
   if ( wtX != NULL ) wtX[1] = fabs( nlTp4.M() - 171.2 );
   //cout<<" WMass= "<< nhWp4.M() <<" had W sigma = "<< sgm_nhW <<" TMass = "<< nhTp4.M() <<" had T sigma = "<< sgm_nhT <<endl;
   //cout<<" W X2 = "<< sqrt(nhW_X2) <<" W prob = "<< nprob_hW  <<"  T X2 = "<< sqrt(nhT_X2) <<" T prob = "<< nprob_hT << endl;

   //if ( sqrt(nhW_X2) > 3 && sqrt(nhT_X2 < 3) ) prob = -1 ;
   if ( sqrt(nhW_X2) > 5  ) prob = -1 ;
   if ( combX > 7 ) prob = -1 ;
   return prob ;
}

double HadWMassFitter::Get2BodySigma( TLorentzVector v1, TLorentzVector v2 ) {

   double sgm1 = 0.7*sqrt( (5.6*5.6) +  ( 1.25*1.25*v1.Pt() ) +  ( 0.033*0.033*v1.Pt()*v1.Pt() ) );
   double sgm2 = 0.7*sqrt( (5.6*5.6) +  ( 1.25*1.25*v2.Pt() ) +  ( 0.033*0.033*v2.Pt()*v2.Pt() ) );

   double A1 = ( v2.E()/v1.E() ) * sgm1*sgm1 ;
   double A2 = ( v1.E()/v2.E() ) * sgm2*sgm2 ;
   double sgm12 = 0.5 * sqrt( A1 + A2 ) ;

   return sgm12 ;
}

double HadWMassFitter::Get3BodySigma( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3 ) {

   TLorentzVector v0 = v1 + v2 ;
   double sgm0 = Get2BodySigma( v1, v2 ) ;
   double sgm3 = 0.7*sqrt( (5.6*5.6) +  ( 1.25*1.25*v3.Pt() ) +  ( 0.033*0.033*v3.Pt()*v3.Pt() ) );

   double A0 = ( v3.E()/v0.E() ) * sgm0*sgm0 ;
   double A3 = ( v0.E()/v3.E() ) * sgm3*sgm3 ;
   double sgm03 = 0.5 * sqrt( A0 + A3 ) ;

   return sgm03 ;
}

vector<TLorentzVector> HadWMassFitter::newVector( TLorentzVector v1,  TLorentzVector v2, TLorentzVector v3, TLorentzVector v4,Double_t *par ){

   double cos1 = v1.Px()/v1.Pt() ;
   double cos2 = v2.Px()/v2.Pt() ;
   double sin1 = v1.Py()/v1.Pt() ;
   double sin2 = v2.Py()/v2.Pt() ;

   double JES1 =  ( fabs( v1.Eta() ) < 1.3 ) ? par[4] : par[5] ;
   double JES2 =  ( fabs( v2.Eta() ) < 1.3 ) ? par[4] : par[5] ;
   double JES3 =  ( fabs( v3.Eta() ) < 1.3 ) ? par[4] : par[5] ;
   double JES4 =  ( fabs( v4.Eta() ) < 1.3 ) ? par[4] : par[5] ;

   double px1 = (v1.Px()+par[0]*cos1)*JES1;
   double py1 = (v1.Py()+par[0]*sin1)*JES1;
   double px2 = (v2.Px()+par[1]*cos2)*JES2;
   double py2 = (v2.Py()+par[1]*sin2)*JES2;

   double pz1 = (v1.Pz()+par[2])*JES1;
   double pz2 = (v2.Pz()+par[3])*JES2;

   double jE1 = sqrt( px1*px1 + py1*py1 + pz1*pz1 ) ;
   double jE2 = sqrt( px2*px2 + py2*py2 + pz2*pz2 ) ;

   TLorentzVector newv1( px1, py1, pz1, jE1 );
   TLorentzVector newv2( px2, py2, pz2, jE2 );
   TLorentzVector newv3 = v3 * JES3;
   TLorentzVector newv4 = v4 * JES4;

   vector<TLorentzVector> newVlist;
   newVlist.push_back( newv1 );
   newVlist.push_back( newv2 );
   newVlist.push_back( newv3 );
   newVlist.push_back( newv4 );
 
   return newVlist ;
}

Double_t HadWMassFitter::jetAngle( TLorentzVector v1, TLorentzVector v2 ){
 
   double v1v2 = v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz() ;
   double v1l  = sqrt( v1.Px()*v1.Px() + v1.Py()*v1.Py() + v1.Pz()*v1.Pz() );
   double v2l  = sqrt( v2.Px()*v2.Px() + v2.Py()*v2.Py() + v2.Pz()*v2.Pz() );
   double cosA = v1v2 / ( v1l*v2l ) ;
   double angle = acos( cosA ) ;
   return angle ;
}

double HadWMassFitter::BTagProbability( vector<double> bCut, int NofB, int jid[], bool isJES ){

   double probb = 0 ;
   if ( n_btag ==2 && NofB >= 2 && bCut[ jid[2] ] >= bTh && bCut[ jid[3]] >= bTh ) probb = 1;
   if ( n_btag ==0 && NofB >= 2 && bCut[ jid[2] ] >= bTh && bCut[ jid[3]] >= bTh ) probb = 1;
   if ( n_btag < 2 && NofB == 1 && bCut[ jid[2] ] >= bTh  ) probb = 1;
   if ( n_btag < 2 && NofB == 1 && bCut[ jid[3] ] >= bTh  ) probb = 1;
   if ( n_btag >=0 && NofB == 0  ) probb = 0 ;
   if ( n_btag == -1 && !isJES ) probb = 1 ;
   if ( n_btag == -1 && isJES ) {
      if ( NofB == 0 )   probb = 1;
      if ( NofB == 1 && bCut[ jid[2] ] >= bTh )  probb = 1;
      if ( NofB == 1 && bCut[ jid[3] ] >= bTh )  probb = 1;
      if ( NofB >= 2 && bCut[ jid[2] ] >= bTh && bCut[ jid[3]] >= bTh )  probb = 1;
      if ( NofB >= 2 && bCut[ jid[2] ] >= bTh && bCut[ jid[3]] >= bTh )  probb = 1;
   }
   return probb ;
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
void HadWMassFitter::ReFitSolution( string fileName, recoObj* wObj, int type, bool isMCMatched ) {

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

  double norm = 0;
  vector<TLorentzVector> tlist0 ;
  vector<TLorentzVector> tlist1 ;
  vector<TLorentzVector> tlist2 ;
  vector<TLorentzVector> tlist3 ;
  vector<TLorentzVector> tlist4 ;
  vector<TLorentzVector> tlist5 ;
  vector<double> wProb ;
  int N_allEvt = 0 ;
  int N_0bEvt = 0 ;
  int N_1bEvt = 0 ;
  int N_2bEvt = 0 ;
  int N_btagged = 0 ;

  for ( int j= 0 ; j< tr1->GetEntries() ; j++ ) {
      tr1->GetEntry(j);

      // reset all containers at the beginning of the event
      if ( evtId0 != evtId1 ) {
         N_allEvt++;
         if ( N_btagged == 0 ) N_0bEvt++;
         if ( N_btagged == 1 ) N_1bEvt++;
         if ( N_btagged >= 2 ) N_2bEvt++;
         N_btagged = 0 ;
         // record
         for (size_t i=0; i< wProb.size(); i++) {
             // selected best permutation
             if ( norm == 0 ) continue;
             //if ( type == 0 || isMCMatched ) {
                wObj->gethad( tlist0[i], tlist1[i], tlist2[i] );
                wObj->getlep( tlist3[i], tlist4[i], tlist5[i] );
                //wObj->Fillh( wProb[i]/norm  );
                double weighting = 1.0 / wProb.size() ;
                wObj->Fillh( weighting );
             //}
         }
         evtId0 = evtId1 ;
	 tlist0.clear();
	 tlist1.clear();
	 tlist2.clear();
	 tlist3.clear();
	 tlist4.clear();
	 tlist5.clear();
	 wProb.clear();
	 jvs.clear();
	 nvs.clear();
	 mvs.clear();
         bCut.clear();
         norm = 0;

         // accumulate objects in the event
         jvs = fitInput->GetJets( tr2, evtId1, bCut, jet_str );
         nvs = fitInput->GetNeutrinos( tr3, evtId1, neu_str );
         mvs = fitInput->GetMuons( tr4, evtId1, mu_str );
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
      N_btagged = NofB ;
      //cout <<" N of B ="<< NofB <<endl;
      // btagger switch
      if ( NofB == 0  && n_btag >= 0 ) continue;

      // b Cut  -- need more studies
      double probb = BTagProbability( bCut, NofB, jid1, isJES );

      // refit the W mass
      Double_t para[6] = {0,0,0,0,1,1};
      Double_t errs[6];
      if ( type > 0 )  { 
         FitW( jvs[ jid1[0] ] , jvs[ jid1[1] ], para, errs, isJES );
      }
      vector<TLorentzVector> newJs = newVector( jvs[ jid1[0] ], jvs[ jid1[1] ], jvs[ jid1[2] ], jvs[ jid1[3] ] , para );

      // determine the weighting of each entry
      Double_t wtX[2] = {0,0} ;
      double probWT = KinematicProb( jid1, jvs, mvs[0], nvs[ jid1[4]-1 ], para, wtX ); 
      if ( probWT == -1 && type > 0 && !isMCMatched ) continue ;

      double aProb = ( type > 0 && !isMCMatched ) ? probWT : 1 ;
      if ( n_btag >= 0 ) aProb = aProb*probb ;

      //cout<<" -------------------------------------------- "<<endl;
      //cout<<" Event "<< evtId1 <<" $$ jid("<<jid1[0]<<","<<jid1[1]<<","<<jid1[2]<<","<<jid1[3]<<","<<jid1[4]<<")" ;
      //cout<<" prob = "<< aProb <<" wX= "<< wtX[0]<<" dT ="<< wtX[1] <<endl;
      // record the result and normalization
      tlist0.push_back( newJs[0] );
      tlist1.push_back( newJs[1] );
      tlist2.push_back( newJs[2] );
      tlist3.push_back( newJs[3] );
      tlist4.push_back( nvs[ jid1[4]-1 ] );
      tlist5.push_back( mvs[0]  );
      wProb.push_back( aProb );
      norm += aProb;  
  }

  cout<<" allEvt = "<< N_allEvt <<"  0B = "<< N_0bEvt <<"  1B = "<< N_1bEvt <<" 2B = "<< N_2bEvt <<endl;
}


// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
// For ensemble test
void HadWMassFitter::ReFitSolution( string fileName, recoObj* wObj, vector< pair<int,int> >& evtlist,  int type ) {

  // get files and trees
  TFile*  file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName, "solTt",   file );
  TTree* tr2 = fitInput->GetTree( fileName, "selJet",  file );
  TTree* tr3 = fitInput->GetTree( fileName, "solNeu",  file );
  TTree* tr4 = fitInput->GetTree( fileName, "selMu",   file );

  bool isJES = ( type == 2 ) ? true : false ;

  // event solution tree
  double wm, tm , prob;
  int evtId1, entId1, entSz1, jid1[5] ;
  tr1->SetBranchAddress( "evtId" ,&evtId1 );
  tr1->SetBranchAddress( "entId" ,&entId1 );
  tr1->SetBranchAddress( "entSz" ,&entSz1 );
  tr1->SetBranchAddress( "JId"   ,&jid1   );
  tr1->SetBranchAddress( "hadWM" ,&wm   );
  tr1->SetBranchAddress( "hadTM" ,&tm   );
  tr1->SetBranchAddress( "probW" ,&prob   );

  int neu_str = 0 ;
  int jet_str = 0 ;
  int mu_str  = 0 ;
  vector<double> bCut ;
  vector<TLorentzVector> jvs;
  vector<TLorentzVector> nvs;
  vector<TLorentzVector> mvs;

  double norm = 0;
  std::vector<TLorentzVector> tlist0 ;
  std::vector<TLorentzVector> tlist1 ;
  std::vector<TLorentzVector> tlist2 ;
  std::vector<TLorentzVector> tlist3 ;
  std::vector<TLorentzVector> tlist4 ;
  std::vector<TLorentzVector> tlist5 ;
  std::vector<double> wProb ;
  int N_allEvt = 0 ;
  int N_0bEvt = 0 ;
  int N_1bEvt = 0 ;
  int N_2bEvt = 0 ;

  cout<<" -------------------------------------------- "<<endl;
  for ( size_t j= 0 ; j< evtlist.size() ; j++ ) {

      tr1->GetEntry( evtlist[j].second-1 );
      cout<<" *** entId = "<< evtlist[j].second <<" entId check = "<< entId1 <<endl;

       cout<<" Event "<< evtId1 <<" $$ jid("<<jid1[0]<<","<<jid1[1]<<","<<jid1[2]<<","<<jid1[3]<<","<<jid1[4]<<")"<<endl;
      // reset all containers at the beginning of the event
      N_allEvt++;
      tlist0.clear();
      tlist1.clear();
      tlist2.clear();
      tlist3.clear();
      tlist4.clear();
      tlist5.clear();
      wProb.clear();
      jvs.clear();
      nvs.clear();
      mvs.clear();
      bCut.clear();
      norm = 0;

      // accumulate objects in the event
      jvs = fitInput->GetJets( tr2, evtId1, bCut, jet_str );
      nvs = fitInput->GetNeutrinos( tr3, evtId1, neu_str );
      mvs = fitInput->GetMuons( tr4, evtId1, mu_str );
      cout<<" j size= "<< jvs.size() <<" nu size = "<< nvs.size() <<" mu size = "<< mvs.size() <<endl; 

      // number of b tagged
      int NofB = 0;
      for (size_t k=0; k< 4; k++) {
          if (bCut[k] > bTh) NofB++ ;
      }
      cout <<" N of B ="<< NofB <<endl;
      // statistic of btagging
      if ( NofB == 0 ) N_0bEvt++;
      if ( NofB == 1 ) N_1bEvt++;
      if ( NofB >= 2 ) N_2bEvt++;

      int headSol = evtlist[j].second - 1 ;
      int EndSol  = evtlist[j].second + entSz1 -1;
      cout<<"      size of permutation = "<< entSz1 <<" from "<< headSol << " to "<< EndSol <<endl ;
      for(int i = headSol; i < EndSol; i++ ) {
         tr1->GetEntry(i);
   
         // exclude the 5th jet
         if ( jid1[0] == 4 || jid1[1] == 4 || jid1[2] == 4 || jid1[3] == 4 ) continue;
	 // jet pt reselection
	 if ( !JetPtFilter( jvs ) ) continue;

	 // btagger switch
	 if ( NofB == 0  && n_btag >= 0 ) continue;

	 // b Cut  -- need more studies
         double probb = BTagProbability( bCut, NofB, jid1, isJES );

         // refit the W mass
         Double_t para[6] = {0,0,0,0,1,1};
         Double_t errs[6];
         cout<<" jvs size = "<< jvs.size() <<endl;
         if ( type > 0 )  { 
            FitW( jvs[ jid1[0] ] , jvs[ jid1[1] ], para, errs, isJES );
         }
         vector<TLorentzVector> newJs = newVector( jvs[ jid1[0] ], jvs[ jid1[1] ], jvs[ jid1[2] ], jvs[ jid1[3] ] , para );
         cout<<" new jvs size = "<< newJs.size() <<endl;
         // determine the weighting of each entry
         double probWT = KinematicProb( jid1, jvs, mvs[0], nvs[ jid1[4]-1 ], para ); 
         if ( probWT == -1 && type > 0 ) continue ;

	 double aProb = ( type > 0 ) ? probWT : 1. ;
	 if ( n_btag >= 0 || type > 0 ) aProb = aProb*probb ;

	 // record the result and normalization
	 if ( aProb != 0 ) {
            tlist0.push_back( newJs[0] );
	    tlist1.push_back( newJs[1] );
	    tlist2.push_back( newJs[2] );
	    tlist3.push_back( newJs[3] );
	    tlist4.push_back( nvs[ jid1[4]-1 ] );
	    tlist5.push_back( mvs[0]  );
	    wProb.push_back( aProb );
	    norm += aProb;   
         }
      }

      // record
      cout<<" accepted permutations "<< wProb.size() <<" normal =  "<< norm <<endl;
      for (size_t i=0; i< wProb.size(); i++) {
          if ( norm == 0 ) continue;
	  wObj->gethad( tlist0[i], tlist1[i], tlist2[i] );
	  wObj->getlep( tlist3[i], tlist4[i], tlist5[i] );
          // normalized the contributions from each events 
          /// btagging will have different entries from each event => normalized to 24/event
          double weighting = 1.0 / wProb.size() ;
	  if (type == 0 ) wObj->Fillh( 24*wProb[i]/norm  );
	  if ( type > 0 ) wObj->Fillh( 24*weighting );
      }
  }
  cout<<" allEvt = "<< N_allEvt <<"  0B = "<< N_0bEvt <<"  1B = "<< N_1bEvt <<" 2B = "<< N_2bEvt <<endl;

}


void HadWMassFitter::GetPermutes( int njets, vector<jlist>& jlistV ) {
 

    vector<int> pools ;
    for (int i=0; i < njets; i++) {    
        pools.push_back(i) ;
    }

    jlistV.clear() ;

    do {
       for (int i=0; i< njets-1 ; i++) {
           for (int j=i+1; j < njets; j++ ) {
               if ( pools[0] == i && pools[1] == j ) {
                  jlist theEntry ;
                  theEntry.w1 = pools[0] ;
                  theEntry.w2 = pools[1] ;
                  theEntry.bh = pools[2] ;
                  if ( pools.size() >= 4 ) theEntry.bl = pools[3] ;
                  if ( pools.size() < 4  ) theEntry.bl = 0 ;
                  jlistV.push_back( theEntry ) ;
               }
           }
       }
    } while ( next_permutation( pools.begin() ,pools.end() ) ) ;
   
}

vector<TLorentzVector> HadWMassFitter::GetLorentzVector( jlist Ls, double jpx[], double jpy[], double jpz[], double jE[] ){

    TLorentzVector wj1( jpx[Ls.w1], jpy[Ls.w1], jpz[Ls.w1], jE[Ls.w1] );
    TLorentzVector wj2( jpx[Ls.w2], jpy[Ls.w2], jpz[Ls.w2], jE[Ls.w2] );
    TLorentzVector bjh( jpx[Ls.bh], jpy[Ls.bh], jpz[Ls.bh], jE[Ls.bh] );
    TLorentzVector bjl( jpx[Ls.bl], jpy[Ls.bl], jpz[Ls.bl], jE[Ls.bl] );

    vector<TLorentzVector> tjets ;
    tjets.push_back( wj1 ) ;
    tjets.push_back( wj2 ) ;
    tjets.push_back( bjh ) ;
    tjets.push_back( bjl ) ;
    
    return tjets ;
}

// For New SolTree Format
void HadWMassFitter::ReFitSolution1( string fileName, int njets, recoObj* wObj, vector<int>* evtlist ) {

  // get files and trees
  TString treeName ;
  if ( njets == 3 ) treeName = "mu3Jets" ;
  if ( njets == 4 ) treeName = "mu4Jets" ;
  TFile*  file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName, treeName,  file );

  // event solution tree
  const static int sz = njets + 1;
  int evtId, nJ, nM, nN;
  double jpx[sz], jpy[sz], jpz[sz], jE[sz], bCut[sz] ;
  double npx[sz], npy[sz], npz[sz], nE[sz] ;
  double mpx[sz], mpy[sz], mpz[sz], mE[sz] ;
  tr1->SetBranchAddress("evtId"    ,&evtId);
  tr1->SetBranchAddress("nJ"       ,&nJ);
  tr1->SetBranchAddress("nNu"      ,&nN);
  tr1->SetBranchAddress("nMu"      ,&nM);

  tr1->SetBranchAddress("bTh"      ,&bCut);
  tr1->SetBranchAddress("jpx"      ,&jpx);
  tr1->SetBranchAddress("jpy"      ,&jpy);
  tr1->SetBranchAddress("jpz"      ,&jpz);
  tr1->SetBranchAddress("jE"       ,&jE);

  tr1->SetBranchAddress("npx"      ,&npx);
  tr1->SetBranchAddress("npy"      ,&npy);
  tr1->SetBranchAddress("npz"      ,&npz);
  tr1->SetBranchAddress("nE"       ,&nE);

  tr1->SetBranchAddress("mpx"      ,&mpx);
  tr1->SetBranchAddress("mpy"      ,&mpy);
  tr1->SetBranchAddress("mpz"      ,&mpz);
  tr1->SetBranchAddress("mE"       ,&mE);

  // clean the permutation lists and fill the permuat
  vector<jlist> jlistV ;
  GetPermutes( njets, jlistV );
  vector<double> bCutV;

  int loopEnd = tr1->GetEntries() ;
  int idx = 0;
  if ( evtlist != NULL ) loopEnd = evtlist->size() ;

  for ( int j= 0 ; j< loopEnd ; j++ ) {
     
      if ( evtlist != NULL ) idx = (*evtlist)[j] ;
      if ( evtlist == NULL ) idx = j ;
      tr1->GetEntry(idx);
       
      // number of b tagged 
      int NofB = 0;
      bCutV.clear() ;
      for (int k=0; k< nJ; k++) {
          bCutV.push_back( bCut[k] );
          if ( bCut[k] > bTh ) NofB++ ;
      }

      // set up the weighting of each permutation => make each event having the same amount of contribution
      double weighting = 1. / (2*jlistV.size()) ;
      if ( n_btag > -1 && NofB == 0 ) weighting = 1. / TMath::Factorial( njets ) ;  
      if ( n_btag > -1 && NofB == 1 ) weighting = TMath::Factorial(njets-4) / ( 2*TMath::Factorial(njets-1) );  
      if ( n_btag > -1 && NofB >= 2 ) weighting = TMath::Factorial(njets-4) / ( 2*TMath::Factorial(njets-2) ); 

      // loop over all possible permutation    
      for (size_t i=0; i< jlistV.size(); i++) {

          // get 4 momentum for jets and muon/neutrino
          // most neutrino pz has 2 solutions
          for (int k = 0; k< nN; k++) {
              vector<TLorentzVector> tjets = GetLorentzVector( jlistV[i], jpx, jpy, jpz, jE );
              TLorentzVector  muP4( mpx[0], mpy[0], mpz[0], mE[0] );
              TLorentzVector neuP4( npx[k], npy[k], npz[k], nE[k] );

              // btagging 
              int jid[4] = { jlistV[i].w1, jlistV[i].w2, jlistV[i].bh, jlistV[i].bl };
              double probb = BTagProbability( bCutV, NofB, jid, false );
              if ( n_btag > -1 && probb == 0 ) continue ;

              // recording the result
              wObj->gethad( tjets[0], tjets[1], tjets[2] );
	      wObj->getlep( tjets[3], neuP4, muP4 );
	      wObj->Fillh( weighting );
          }
      }
  }

} 

// For New SolTree Format, with JES/JEC, 
// JEC : type = 1 , JES : type = 2
void HadWMassFitter::ReFitSolution1( string fileName, int njets, recoObj* wObj, int type, vector<int>* evtlist ) {

  // get files and trees
  TString treeName ;
  if ( njets == 3 ) treeName = "mu3Jets" ;
  if ( njets == 4 ) treeName = "mu4Jets" ;
  TFile*  file = NULL ;
  TTree* tr1 = fitInput->GetTree( fileName, treeName,  file );

  bool isJES = ( type == 2 ) ? true : false ;

  // event solution tree
  const static int sz = njets + 1;
  int evtId, nJ, nM, nN;
  double jpx[sz], jpy[sz], jpz[sz], jE[sz], bCut[sz] ;
  double npx[sz], npy[sz], npz[sz], nE[sz] ;
  double mpx[sz], mpy[sz], mpz[sz], mE[sz] ;
  tr1->SetBranchAddress("evtId"    ,&evtId);
  tr1->SetBranchAddress("nJ"       ,&nJ);
  tr1->SetBranchAddress("nNu"      ,&nN);
  tr1->SetBranchAddress("nMu"      ,&nM);

  tr1->SetBranchAddress("bTh"      ,&bCut);
  tr1->SetBranchAddress("jpx"      ,&jpx);
  tr1->SetBranchAddress("jpy"      ,&jpy);
  tr1->SetBranchAddress("jpz"      ,&jpz);
  tr1->SetBranchAddress("jE"       ,&jE);

  tr1->SetBranchAddress("npx"      ,&npx);
  tr1->SetBranchAddress("npy"      ,&npy);
  tr1->SetBranchAddress("npz"      ,&npz);
  tr1->SetBranchAddress("nE"       ,&nE);

  tr1->SetBranchAddress("mpx"      ,&mpx);
  tr1->SetBranchAddress("mpy"      ,&mpy);
  tr1->SetBranchAddress("mpz"      ,&mpz);
  tr1->SetBranchAddress("mE"       ,&mE);

  // clean the permutation lists and fill the permuat
  vector<jlist> jlistV ;
  GetPermutes( njets, jlistV );
  vector<double> bCutV;
  vector<TLorentzVector> vlist0;
  vector<TLorentzVector> vlist1;
  vector<TLorentzVector> vlist2;
  vector<TLorentzVector> vlist3;
  vector<TLorentzVector> vlist4;
  vector<TLorentzVector> vlist5;

  int loopEnd = tr1->GetEntries() ;
  int idx = 0;
  if ( evtlist != NULL ) loopEnd = evtlist->size() ;

  for ( int j= 0 ; j< tr1->GetEntries() ; j++ ) {
      if ( evtlist != NULL ) idx = (*evtlist)[j] ;
      if ( evtlist == NULL ) idx = j ;
      tr1->GetEntry(idx);

      vlist0.clear();   
      vlist1.clear();   
      vlist2.clear();   
      vlist3.clear();   
      vlist4.clear();   
      vlist5.clear();   

      // number of b tagged 
      int NofB = 0;
      bCutV.clear() ;
      for (int k=0; k< nJ; k++) {
          bCutV.push_back( bCut[k] );
          if ( bCut[k] > bTh ) NofB++ ;
      }

      // loop over all possible permutation    
      for (size_t i=0; i< jlistV.size(); i++) {

          // get 4 momentum for jets and muon/neutrino
          // most neutrino pz has 2 solutions
          for (int k = 0; k< nN; k++) {
              vector<TLorentzVector> tjets = GetLorentzVector( jlistV[i], jpx, jpy, jpz, jE );
              TLorentzVector  muP4( mpx[0], mpy[0], mpz[0], mE[0] );
              TLorentzVector neuP4( npx[k], npy[k], npz[k], nE[k] );

              // btagging 
              int jid[4] = { jlistV[i].w1, jlistV[i].w2, jlistV[i].bh, jlistV[i].bl };
              double probb = BTagProbability( bCutV, NofB, jid, false );
              if ( n_btag > -1 && probb == 0 ) continue ;

              // refit the W mass if desire
	      Double_t para[6] = {0,0,0,0,1,1};
	      Double_t errs[6];
	      if ( type > 0 )  { 
                 FitW( tjets[0] , tjets[1], para, errs, isJES );
              }
              vector<TLorentzVector> newJs = newVector( tjets[0], tjets[1], tjets[2], tjets[3], para );

              // determine the weighting of each entry
	      Double_t wtX[2] = {0,0} ;
              int jid1[4] = {0,1,2,3} ;
	      double probWT = KinematicProb( jid1, tjets, muP4, neuP4, para, wtX ); 
	      if ( probWT == -1 && type > 0 ) continue ;

              // perserv the results 
              vlist0.push_back( newJs[0] );
              vlist1.push_back( newJs[1] );
              vlist2.push_back( newJs[2] );
              vlist3.push_back( newJs[3] );
              vlist4.push_back(  neuP4  );
              vlist5.push_back(   muP4  );
          }
      }
      double weighting = 1.0 /  vlist0.size() ;
      for (size_t i=0; i< vlist0.size(); i++) {
          wObj->gethad( vlist0[0], vlist1[1], vlist2[2] );
	  wObj->getlep( vlist3[3], vlist4[4], vlist5[5] );
	  wObj->Fillh( weighting );
      }
  }

}

// MCMatching information for new soltree
// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void HadWMassFitter::MCSolution( string fileName, recoObj* wObj, int type ) {

  // get files and trees
  TFile*  file = NULL ;
  TTree* tr0 = fitInput->GetTree( fileName, "mcmTt",   file );
  TTree* tr1 = fitInput->GetTree( fileName, "mu4Jets", file );

  bool isJES = ( type == 2 ) ? true : false ;

  // event solution tree
  int evtId0, jid0[5] ;
  tr0->SetBranchAddress( "evtId" ,&evtId0 );
  tr0->SetBranchAddress( "JId"   ,&jid0   );

  // event solution tree
  const static int sz = 5;
  int evtId, nJ, nM, nN;
  double jpx[sz], jpy[sz], jpz[sz], jE[sz], bCut[sz] ;
  double npx[sz], npy[sz], npz[sz], nE[sz] ;
  double mpx[sz], mpy[sz], mpz[sz], mE[sz] ;
  tr1->SetBranchAddress("evtId"    ,&evtId);
  tr1->SetBranchAddress("nJ"       ,&nJ);
  tr1->SetBranchAddress("nNu"      ,&nN);
  tr1->SetBranchAddress("nMu"      ,&nM);

  tr1->SetBranchAddress("bTh"      ,&bCut);
  tr1->SetBranchAddress("jpx"      ,&jpx);
  tr1->SetBranchAddress("jpy"      ,&jpy);
  tr1->SetBranchAddress("jpz"      ,&jpz);
  tr1->SetBranchAddress("jE"       ,&jE);

  tr1->SetBranchAddress("npx"      ,&npx);
  tr1->SetBranchAddress("npy"      ,&npy);
  tr1->SetBranchAddress("npz"      ,&npz);
  tr1->SetBranchAddress("nE"       ,&nE);

  tr1->SetBranchAddress("mpx"      ,&mpx);
  tr1->SetBranchAddress("mpy"      ,&mpy);
  tr1->SetBranchAddress("mpz"      ,&mpz);
  tr1->SetBranchAddress("mE"       ,&mE);

  vector<double> bCutV ;

  int idx = 0 ;
  for ( int j= 0 ; j< tr0->GetEntries() ; j++ ) {
      tr0->GetEntry(j);

      cout<<" MC Matched Evt = "<< evtId0 << endl; 
      for ( int i=idx; i< tr1->GetEntries() ; i++) {
          tr1->GetEntry(i);
          if ( evtId == evtId0  ) {
             idx = i ;
             
             cout<<" idx = "<< idx <<" id:"<< evtId<<" id0:"<< evtId0 << endl;
             // number of b tagged 
	     int NofB = 0;
	     bCutV.clear() ;
	     for (int k=0; k< nJ; k++) {
                 bCutV.push_back( bCut[k] );
                 if ( bCut[k] > bTh ) NofB++ ;
             }
             // get p4 for all objects
             jlist mclist;
             mclist.w1 = jid0[0] ;
             mclist.w2 = jid0[1] ;
             mclist.bh = jid0[2] ;
             mclist.bl = jid0[3] ;
             vector<TLorentzVector> tjets = GetLorentzVector( mclist, jpx, jpy, jpz, jE );
             TLorentzVector  muP4( mpx[0], mpy[0], mpz[0], mE[0] );
             TLorentzVector neuP4( npx[ jid0[4]-1 ], npy[ jid0[4]-1 ], npz[ jid0[4]-1 ], nE[ jid0[4]-1 ] );

             // btagging 
	     double probb = BTagProbability( bCutV, NofB, jid0, isJES );
	     if ( n_btag > -1 && probb == 0 ) continue ;

             // recording the result
	     wObj->gethad( tjets[0], tjets[1], tjets[2] );
	     wObj->getlep( tjets[3], neuP4, muP4 );
	     wObj->Fillh( 1. );
          }
      }
  }

}
