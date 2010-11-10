#include "HadWMassFitter.h"

HadWMassFitter::HadWMassFitter() {

  fitInput = new MassAnaInput();
  pseudoExp = new PseudoExp();

  fitInput->GetParameters( "Path",       &hfolder );
  fitInput->GetParameters( "M2M3Cuts",   &M2M3Cut );
  fitInput->GetParameters( "LepM2tCutL", &LepM2tCutL );
  fitInput->GetParameters( "dM3Cut",     &dM3Cut );
  fitInput->GetParameters( "JES",        &JES );
  fitInput->GetParameters( "JetCuts",    &jetCuts );
  fitInput->GetParameters( "MuonCuts",   &muonCuts );
  fitInput->GetParameters( "Inclusive",  &Inclusive );
  fitInput->GetParameters( "bThreshold", &bTh );
  fitInput->GetParameters( "n_btag",     &n_btag );
  fitInput->GetParameters( "JESType",    &JESType );
  inclu = ( Inclusive == "YES" ) ? true : false ;

  normMCData = true ;
}

HadWMassFitter::~HadWMassFitter(){

  delete fitInput ;
  delete pseudoExp ;

}

//FitVectors * fVec;
TMinuit minuit(6) ;

void HadWMassFitter::WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag) {

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
   //double chi2 = fabs( wp4.M() - 80.4 );
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

double HadWMassFitter::KinematicProb( vector<TLorentzVector> vlist, TLorentzVector mV4, TLorentzVector nV4, Double_t *par, Double_t* wtX ) {

   //cout<<" jid = ("<< jid[0]<<" , "<< jid[1]<<" , "<< jid[2]<<" , "<< jid[3]<<" ) "<<endl;
   //cout<<" vlist size = "<< vlist.size() <<endl;
   vector<TLorentzVector> nvlist = newVector( vlist[0], vlist[1], vlist[2], vlist[3] , par );

   // new hadronic W p4
   TLorentzVector nhWp4 = nvlist[0] +  nvlist[1] ;
   double sgm_nhW = Get2BodySigma( nvlist[0], nvlist[1] ) ;

   // new hadronic Top p4
   TLorentzVector nhTp4 = nvlist[0] +  nvlist[1] + nvlist[2] ;
   double sgm_nhT = Get3BodySigma( nvlist[0], nvlist[1], nvlist[2] ) ;

   // new leptonic Top p4, not used yet
   TLorentzVector nlTp4 = nvlist[3] +  mV4 + nV4 ;

   double nhW_X2 = ( nhWp4.M() -  80.4 )*( nhWp4.M() -  80.4 ) / ( sgm_nhW*sgm_nhW ) ;
   double nhT_X2 = ( nhTp4.M() - 172.5 )*( nhTp4.M() - 172.5 ) / ( sgm_nhT*sgm_nhT ) ;

   double nprob_hW = exp( -0.5*nhW_X2 ) ;
   double nprob_hT = exp( -0.5*nhT_X2 ) ;

   double prob =  nprob_hW*nprob_hT ;
   double combX =  sqrt( nhW_X2 + nhT_X2 );

   if ( wtX != NULL ) wtX[0] = sqrt( nhW_X2 );
   if ( wtX != NULL ) wtX[1] = fabs( nlTp4.M() - 172.5 );
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

double HadWMassFitter::jetAngle( TLorentzVector v1, TLorentzVector v2 ){
 
   double v1v2 = v1.Px()*v2.Px() + v1.Py()*v2.Py() + v1.Pz()*v2.Pz() ;
   double v1l  = sqrt( v1.Px()*v1.Px() + v1.Py()*v1.Py() + v1.Pz()*v1.Pz() );
   double v2l  = sqrt( v2.Px()*v2.Px() + v2.Py()*v2.Py() + v2.Pz()*v2.Pz() );
   double cosA = v1v2 / ( v1l*v2l ) ;
   double angle = acos( cosA ) ;
   return angle ;
}

double HadWMassFitter::minTheta_WDecay( TLorentzVector v1, TLorentzVector v2 ){

   TLorentzVector vW = v1 + v2 ;
   double pW   = sqrt( vW.Px()*vW.Px() + vW.Py()*vW.Py() + vW.Pz()*vW.Pz() );
   double theta = ( 80.4/sqrt(2) )*( pW / ( 80.4*80.4 + pW*pW ) ) ;
   double angle = 2.* atan( theta ) ;
   return angle ;
}

double HadWMassFitter::BTagProbability( double bCut[], int NofB, int jid[], bool isJES ){

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

vector<TLorentzVector> HadWMassFitter::GetLorentzVector( jlist Ls, double jpx[], double jpy[], double jpz[], double jE[] ){

    TLorentzVector wj1( jpx[Ls.w1], jpy[Ls.w1], jpz[Ls.w1], jE[Ls.w1] );
    TLorentzVector wj2( jpx[Ls.w2], jpy[Ls.w2], jpz[Ls.w2], jE[Ls.w2] );
    TLorentzVector bjh( jpx[Ls.bh], jpy[Ls.bh], jpz[Ls.bh], jE[Ls.bh] );
    TLorentzVector bjl( jpx[Ls.bl], jpy[Ls.bl], jpz[Ls.bl], jE[Ls.bl] );

    vector<TLorentzVector> tjets ;
    tjets.push_back( wj1*JES ) ;
    tjets.push_back( wj2*JES ) ;
    tjets.push_back( bjh*JES ) ;
    tjets.push_back( bjl*JES ) ;
    
    return tjets ;
}

vector<TLorentzVector> HadWMassFitter::GetLorentzVector( vector<TLorentzVector>& oblist, jlist Ls ){

    TLorentzVector wj1 = oblist[Ls.w1] ;
    TLorentzVector wj2 = oblist[Ls.w2] ;
    TLorentzVector bjh = oblist[Ls.bh] ;
    TLorentzVector bjl = oblist[Ls.bl] ;

    vector<TLorentzVector> tjets ;
    tjets.push_back( wj1*JES ) ;
    tjets.push_back( wj2*JES ) ;
    tjets.push_back( bjh*JES ) ;
    tjets.push_back( bjl*JES ) ;
    
    return tjets ;
}

vector<TLorentzVector> HadWMassFitter::NeuP4Solution( TLorentzVector muP4, TLorentzVector neuP4 ) {

    vector<TLorentzVector> neuSol ;
    double xW  = muP4.Px() + neuP4.Px() ;
    double yW  = muP4.Py() + neuP4.Py() ;
    double nuPt2 = ( neuP4.Px()*neuP4.Px() ) + ( neuP4.Py()*neuP4.Py() );

    double zl = muP4.Pz() ;
    double El = muP4.E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    double EtW  = muP4.Pt() + neuP4.E() ;
    double PtW2 = (xW*xW) + (yW*yW);
    double MtW2 = (EtW*EtW) - PtW2 ;

    if ( MtW2 < 0. ||  (B*B) < (4.*A*C) ) {
       TLorentzVector np4( neuP4.Px(), neuP4.Py(), 0, neuP4.E() );
       neuSol.push_back( np4 );
    } else {

      double nz1 = (-1.*B + sqrt(B*B - 4.*A*C )) / (2.*A) ;
      double nz2 = (-1.*B - sqrt(B*B - 4.*A*C )) / (2.*A) ;
      for (int i=1; i<3; i++ ) {
          double nz = ( i == 1) ? nz1 : nz2 ;
          double ENu2 = ( neuP4.Px()*neuP4.Px() ) + ( neuP4.Py()*neuP4.Py() ) + (nz*nz);
          //double zW = muP4.Pz() + nz ;
          //double EW = muP4.E()  + sqrt(ENu2) ;
          TLorentzVector np4 = TLorentzVector( neuP4.Px(), neuP4.Py(), nz, sqrt(ENu2) );
          neuSol.push_back( np4 ) ;
      }
    }
    return neuSol ;
}


// jets first , last one is muon
vector<TLorentzVector> HadWMassFitter::GetEventObjects( int nj, double jpx[], double jpy[], double jpz[], double jE[], double mpx, double mpy, double mpz, double mE ) {

    vector<TLorentzVector> objs ;

    for ( int i=0; i< nj; i++) {
        TLorentzVector jet_i( jpx[i], jpy[i], jpz[i], jE[i] );
        //if ( fabs(jet_i.Pt() ) < jetCuts[0] ) continue;
        //if ( fabs(jet_i.Eta()) > jetCuts[1] ) continue;
        objs.push_back( jet_i ) ;
    }

    TLorentzVector mp( mpx, mpy, mpz, mE );
  
    objs.push_back( mp ) ;
 
    
    return objs ;
}

void HadWMassFitter::ResetCuts( double m2L, double m2H, double m3L, double m3H, double lepM2tL_threshold, double dM3, bool GetDefault ) {

     // change the M2M3Cuts if desire
     if ( m2L != -1 ) M2M3Cut[0] = m2L ;
     if ( m2H != -1 ) M2M3Cut[1] = m2H ;
     if ( m3L != -1 ) M2M3Cut[2] = m3L ;
     if ( m3H != -1 ) M2M3Cut[3] = m3H ;
     // change the leptonic M2t 
     if ( lepM2tL_threshold != -1 ) LepM2tCutL = lepM2tL_threshold ;
     if ( dM3 != -1 ) dM3Cut = dM3 ;

     if ( GetDefault ) {
        fitInput->GetParameters( "M2M3Cuts",    &M2M3Cut );
	fitInput->GetParameters( "LepM2tCutL", &LepM2tCutL );
	fitInput->GetParameters( "dM3Cut", &dM3Cut );
     }
     cout<<" Using M2L = "<< M2M3Cut[0] <<" M2H= "<<M2M3Cut[1] <<" M3L= "<<M2M3Cut[2] <<" M3H= "<<M2M3Cut[3] ;
     cout<<"   LepM2tL = "<< LepM2tCutL <<"  dM3Cut = "<< dM3Cut  <<endl;
}

void HadWMassFitter::SetMCNormalization( bool normMC ){
     normMCData = normMC ;
}

void HadWMassFitter::SetMuonCuts( double pt_, double eta_, double iso_ ){

     if ( pt_  != -1 ) muonCuts[0] = pt_ ;
     if ( eta_ != -1 ) muonCuts[1] = eta_ ;
     if ( iso_ != -1 ) muonCuts[2] = iso_ ;
     
}

// For New SolTree Format, with JES/JEC, 
// JEC : type = 1 , JES : type = 2 ;  return the efficiency
// evtSplit:  0 = no splitting , negative value = the odd event number, positive value = the even event number
double HadWMassFitter::ReFitSolution( string fileName, recoObj* wObj, int nJets, double scale, vector<int>* evtlist, int evtSplit, bool smearing ) {

  // get files and trees
  TTree* tr1 = fitInput->TreeMap( fileName );

  bool isJES = ( JESType == 2 ) ? true : false ;

  // event solution tree
  // always set bigger arrary size than what has been designed originally
  double jpx[15], jpy[15], jpz[15], jE[15], bCut[15] ;
  int n90[15] ;
  double mpx[2], mpy[2], mpz[2], mE[2], X2[2];
  double mIso[3], npx[3], npy[3], npz[3], nE[3] ;
  int evtId, nJ, nM, nN;

  tr1->SetBranchAddress("evtId"    ,&evtId);
  tr1->SetBranchAddress("nJ"       ,&nJ);
  tr1->SetBranchAddress("nNu"      ,&nN);
  tr1->SetBranchAddress("nMu"      ,&nM);
  tr1->SetBranchAddress("n90"      ,&n90);

  tr1->SetBranchAddress("mpx"      ,&mpx);
  tr1->SetBranchAddress("mpy"      ,&mpy);
  tr1->SetBranchAddress("mpz"      ,&mpz);
  tr1->SetBranchAddress("mE"       ,&mE);
  tr1->SetBranchAddress("X2"       ,&X2);
  tr1->SetBranchAddress("mIso"     ,&mIso);

  tr1->SetBranchAddress("bTh"      ,&bCut);
  tr1->SetBranchAddress("jpx"      ,&jpx);
  tr1->SetBranchAddress("jpy"      ,&jpy);
  tr1->SetBranchAddress("jpz"      ,&jpz);
  tr1->SetBranchAddress("jE"       ,&jE);

  tr1->SetBranchAddress("npx"      ,&npx);
  tr1->SetBranchAddress("npy"      ,&npy);
  tr1->SetBranchAddress("npz"      ,&npz);
  tr1->SetBranchAddress("nE"       ,&nE);

  vector<TLorentzVector> vlist0;
  vector<TLorentzVector> vlist1;
  vector<TLorentzVector> vlist2;
  vector<TLorentzVector> vlist3;
  vector<TLorentzVector> vlist4;
  vector<TLorentzVector> vlist5;

  //int fk = 0;
  //int fkk = 0 ;
  //int gy = 0;
  //int gyy = 0 ;
  int idx = 0;
  int loopEnd = ( evtlist == NULL ) ? tr1->GetEntries() : evtlist->size() ;
  //cout<<" Total Evnet of "<< fileName <<" = "<< loopEnd <<endl;
  int n2sol = 0;
  for ( int j= 0 ; j< loopEnd; j++ ) {
      if ( evtlist != NULL && evtlist->size() == 0 ) break;
      idx = ( evtlist == NULL ) ? j : (*evtlist)[j] ;

      if ( evtSplit == -1 && (j%2) == 0 ) continue ;
      if ( evtSplit ==  1 && (j%2) == 1 ) continue ;
      if ( evtSplit >  1 && (j%evtSplit ) != 0 ) continue ;
      if ( evtSplit < -1 && (j%evtSplit ) != 1 ) continue ;

      tr1->GetEntry(idx);
      if ( nJ < nJets ) continue ;
      //cout<<" EventID = "<< evtId << endl;
      // muon selection
      TLorentzVector muP4( mpx[0], mpy[0], mpz[0], mE[0] );
      if ( muP4.Pt()  < muonCuts[0] ) continue ;
      if ( muP4.Eta() > muonCuts[1] ) continue ;
      if ( mIso[0]    > muonCuts[2] ) continue ;

      // jet selection
      vector<TLorentzVector> objlist;
      int NCountedJets = 0;
      for ( int i = 0; i< nJ ; i ++) {
          TLorentzVector jn( jpx[i], jpy[i], jpz[i], jE[i] );
          // sync cuts , n90 only works for calo jets
          if ( jn.Pt()          < jetCuts[0] ) continue ;
          if ( fabs( jn.Eta() ) > jetCuts[1] ) continue ;
          //if ( n90[i]           <         2  ) continue ;
          objlist.push_back( jn );
          NCountedJets++ ;
      }
      if ( NCountedJets != nJets && !inclu ) continue ;      
      if ( NCountedJets < nJets ) continue ;
      if ( NCountedJets > 6 ) continue ;
      //cout<<" N Jets of the Event = "<< NCountedJets << endl;

      objlist.push_back( muP4 );
      TLorentzVector neuP4( npx[0], npy[0], 0., sqrt( (npx[0]*npx[0]) +  (npy[0]*npy[0]) ) );
      objlist.push_back( neuP4 ) ;

      size_t obsz = objlist.size() ;

      // selection for 2 jets samples
      //if ( nJets == 2 && obsz != 5 ) continue ;
      //if ( nJets == 2 && objlist[2].Pt() > 30. ) continue ;

      if ( smearing ) pseudoExp->PhaseSmearing( objlist, j );
      vector<TLorentzVector> neuP4s = NeuP4Solution( muP4, neuP4*JES ) ;
      //cout<<" neu sol size = "<< neuP4s.size() <<endl;

      /*
      vector<TLorentzVector> neu_ck ;
      for (int k= 0; k< nN; k++) {
          TLorentzVector neuP4_tmp( npx[k], npy[k], npz[k], nE[k] );
          neu_ck.push_back( neuP4_tmp ) ;
      }
      if ( static_cast<int>(neuP4s.size()) != nN ) {
         for (size_t kk= 0; kk < neuP4s.size(); kk++ ){
              cout<<" neuP4s=( "<< neuP4s[kk].Px() <<","<< neuP4s[kk].Py()<<","<< neuP4s[kk].Pz()<<","<< neuP4s[kk].E()<<")"<<endl;
         }
         for (size_t kk= 0; kk < neu_ck.size(); kk++ ){
              cout<<" neuChk=( "<< neu_ck[kk].Px() <<","<< neu_ck[kk].Py()<<","<< neu_ck[kk].Pz()<<","<< neu_ck[kk].E()<<")"<<endl;
         }
      }
      if ( neuP4s.size() == 1 ) fk++  ;
      if ( neuP4s.size() == 2 ) fkk++ ;
      if ( npz[0] == 0 ) gy++ ;
      if ( npz[0] != 0 ) gyy++ ;
      */
     
      vlist0.clear();   
      vlist1.clear();   
      vlist2.clear();   
      vlist3.clear();   
      vlist4.clear();   
      vlist5.clear();   

      // clean the permutation lists and fill them
      vector<jlist> jlistV ;
      fitInput->GetPermutes( obsz-2 , jlistV );

      // number of b tagged 
      int NofB = 0;
      for (int k=0; k< nJ; k++) {
          if ( bCut[k] > bTh ) NofB++ ;
      }

      // loop over all possible permutation  
      int goodPermu = 0 ;
      int goodM2M3  = 0 ;
      for (size_t i=0; i< jlistV.size(); i++) {

          // simple tester to get soft jet in permutation
          //if ( jlistV[i].w1 != 2 && jlistV[i].w2 != 2 && jlistV[i].bh != 2 ) continue ;
          int jid[4] = { jlistV[i].w1, jlistV[i].w2, jlistV[i].bh, jlistV[i].bl };
          //cout<<" Permus = ( "<< jlistV[i].w1 <<","<< jlistV[i].w2 <<","<< jlistV[i].bh <<","<< jlistV[i].bl <<" )   "<<endl;
          //cout<<" BTag = ( "<< bCut[ jid[0] ] <<",   "<< bCut[ jid[1] ] <<",   "<< bCut[ jid[2] ] <<",   "<< bCut[ jid[3] ] <<" )"<<endl;
	  vector<TLorentzVector> tjets = GetLorentzVector( objlist, jlistV[i] );

	  // btagging 
	  double probb = BTagProbability( bCut, NofB, jid, isJES );
	  if ( n_btag > -1 && probb == 0 ) continue ;

          // refit the W mass if desire
	  Double_t para[6] = {0,0,0,0,1,1};
	  Double_t errs[6];
	  if ( JESType > 0 )  FitW( tjets[0] , tjets[1], para, errs, isJES );
          vector<TLorentzVector> newJs = newVector( tjets[0], tjets[1], tjets[2], tjets[3], para );

          // M2 M3 window 
          TLorentzVector vM2 = newJs[0] + newJs[1] ;
          TLorentzVector vM3 = newJs[0] + newJs[1] + newJs[2] ;
          bool goodM2 =  ( vM2.M() > M2M3Cut[0] && vM2.M() < M2M3Cut[1] ) ;
          bool goodM3 =  ( vM3.M() > M2M3Cut[2] && vM3.M() < M2M3Cut[3] ) ;
          if ( goodM2 && goodM3 ) goodM2M3++ ;

          // lepM2 Mt cut
          TLorentzVector vlepM2 = muP4 + neuP4s[0] ;
          double dphi = muP4.DeltaPhi( neuP4s[0] ) ;
          double lepMt2 = 2.*muP4.Pt()*neuP4s[0].Pt()*( 1. - cos(dphi) );
          if ( sqrt(lepMt2) < LepM2tCutL ) continue;

          for (size_t k = 0; k< neuP4s.size() ; k++) {
              // determine the weighting of each entry
              //TLorentzVector neuP4 = ( k == 0 ) ? neuP4_0 : neuP4_1 ;
	      Double_t wtX[2] = {0,0} ;
	      double probWT = KinematicProb( tjets, muP4, neuP4s[k], para, wtX ); 
	      if ( probWT == -1 && JESType > 0 ) continue ;

              // M3 M3 window
              TLorentzVector LM3 = muP4 + neuP4s[k] + newJs[3] ;
              bool goodLM3 = ( fabs( LM3.M() - vM3.M() ) < dM3Cut || obsz < 6 ) ? true : false ;
              if ( M2M3Cut[0] == 0 && M2M3Cut[1] == 999 && M2M3Cut[2] == 0 && M2M3Cut[3] == 999 ) goodLM3 = true ;
 
              if ( goodM2 && goodM3 && goodLM3 ) goodPermu++ ;
              // perserv the results 
              vlist0.push_back( newJs[0] );
              vlist1.push_back( newJs[1] );
              vlist2.push_back( newJs[2] );
              vlist3.push_back( newJs[3] );
              vlist4.push_back(  muP4    );
              vlist5.push_back( neuP4s[k] );
          }
      }
      if (  goodPermu < 1 || vlist0.size() ==  0 ) continue; 
      if ( goodM2M3 < (NCountedJets-3) ) continue; 
      n2sol++ ;

      double weighting = 1. / vlist0.size() ;
      double norm_const = ( normMCData ) ? EvtScaling( nJets, fileName ) : 1. ;
      for (size_t i=0; i< vlist0.size(); i++) {
          wObj->gethad( vlist0[i], vlist1[i], vlist2[i] );
	  wObj->getlep( vlist3[i], vlist4[i], vlist5[i] );
	  wObj->Fillh( weighting, scale*norm_const );
      }
  }
  double eff = ( loopEnd != 0 ) ? ( (n2sol*1.) / (loopEnd*1.) ) : -1 ;
  //cout<<" Passed Event of "<< fileName <<" = "<< n2sol <<"  Eff = "<< eff <<" fk = "<<fk<<" fkk= "<< fkk  ;
  //cout<<"  gy= "<<gy <<" gyy = "<<gyy <<endl;
  return eff ;

}

// MCMatching information for new soltree
// Type : 0 = original jet p4, 1 = JEC tunning , 2 = JES tunning
void HadWMassFitter::MCSolution( string fileName, recoObj* wObj ) {

  // get files and trees
  TTree* tr0 = fitInput->GetTree( fileName, "mcmTt" );
  TTree* tr1 = fitInput->GetTree( fileName, "muJets" );

  bool isJES = ( JESType == 2 ) ? true : false ;

  // event solution tree
  int evtId0, jid0[6] ;
  tr0->SetBranchAddress( "evtId" ,&evtId0 );
  tr0->SetBranchAddress( "pId"   ,&jid0   );

  // event solution tree
  const static int sz = 15;
  int evtId, nJ, nM, nN;
  double jpx[sz], jpy[sz], jpz[sz], jE[sz], bCut[sz] ;
  double npx[2], npy[2], npz[2], nE[2] ;
  double mpx[2], mpy[2], mpz[2], mE[2] ;
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


  int idx = 0 ;
  for ( int j= 0 ; j< tr0->GetEntries() ; j++ ) {
      tr0->GetEntry(j);

      //cout<<" MC Matched Evt = "<< evtId0 << endl; 
      for ( int i=idx; i< tr1->GetEntries() ; i++) {
          tr1->GetEntry(i);
          if ( evtId == evtId0  ) {
             idx = i ;
       
             //cout<<" idx = "<< idx <<" id:"<< evtId<<" id0:"<< evtId0 << endl;
             //cout<<" matchedID = ("<<jid0[0]<<","<<jid0[1]<<","<<jid0[2]<<","<<jid0[3]<<","<<jid0[4]<<","<<jid0[5]<<" )"<<endl;
             // number of b tagged 
	     int NofB = 0;
	     for (int k=0; k< nJ; k++) {
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
             
             // refit the W mass if desire
	     Double_t para[6] = {0,0,0,0,1,1};
	     Double_t errs[6];
	     if ( JESType > 0 ) FitW( tjets[0] , tjets[1], para, errs, isJES );
             
             vector<TLorentzVector> ntjets = newVector( tjets[0], tjets[1], tjets[2], tjets[3], para );

             // btagging 
	     double probb = BTagProbability( bCut, NofB, jid0, isJES );
	     if ( n_btag > -1 && probb == 0 ) continue ;

             // recording the result
	     wObj->gethad( ntjets[0], ntjets[1], ntjets[2] );
	     wObj->getlep( ntjets[3], muP4, neuP4 );
	     wObj->Fillh( 1., 1. );
          }
      }
  }

}

double HadWMassFitter::EvtScaling( int NJets, string fileName ){

    vector<double> qSet;
    fitInput->GetParameters( "qcdNorm", &qSet );
    vector<double> vSet;
    fitInput->GetParameters(  "vjNorm", &vSet );

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

