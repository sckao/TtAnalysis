#include "LepTopMassFitter.h"

LepTopMassFitter::LepTopMassFitter() {

  fitInput = new MassAnaInput();
  pseudoExp = new PseudoExp();

  fitInput->GetParameters( "Path", &hfolder );
  fitInput->GetParameters( "bThreshold", &bTh);
  fitInput->GetParameters( "n_btag", &n_btag);
  fitInput->GetParameters( "n_Jets", &n_Jets);
  fitInput->GetParameters( "IsMES", &IsMES );

}

LepTopMassFitter::~LepTopMassFitter(){

  delete fitInput ;
  delete pseudoExp ;

}

TMinuit lepTMinuit(2);

// v1 = b jet p4,  v2 = w p4
void LepTopMassFitter::FitLepTop( TLorentzVector mP4, TLorentzVector nP4, TLorentzVector bP4, Double_t* pars, Double_t* errs, bool isMES ) {

   // Initialize TMinuit via generic fitter interface with a maximum of 6 params
   LepFitVectors fvs  ;
   fvs.mp4 = mP4 ;
   fvs.np4 = nP4 ;
   fvs.bp4 = bP4 ;
   lepTMinuit.SetObjectFit( &fvs );

   pars[0] = 0.;   // dPhi
   pars[1] = 1.;   // MET Scale
   
   lepTMinuit.SetFCN( LTMFCN );

   // dPhi(0.09) ~ 5.1 degree
   // correction for dPhi
   lepTMinuit.DefineParameter(0, "dPhi", 0,   0.01, -0.05,  0.05 );
   // Re-Scale MET Energy
   lepTMinuit.DefineParameter(1, "MES",  1,   0.01,   0.95,  1.05 );
   // -1  quiet (also suppresse all warnings), 0: normal, 1: verbose
   lepTMinuit.SetPrintLevel(-1);
   //	1 for Chi2 , 0.5 for negative log likelihood
   lepTMinuit.SetErrorDef(0.5);
   
   if ( !isMES ) {
      lepTMinuit.FixParameter(1);
   } 
   /*
   minuit.Migrad();
   minuit.Release(1);
   */
   lepTMinuit.Migrad(); 
   lepTMinuit.GetParameter(0, pars[0], errs[0]);
   lepTMinuit.GetParameter(1, pars[1], errs[1]);
}


// par0 : dPhi , par1 : MET Scale
void LepTopMassFitter::LTMFCN(Int_t &npar, Double_t *, Double_t &tChi2, Double_t *par, Int_t iflag)
{

    LepFitVectors * fVec = (LepFitVectors*) lepTMinuit.GetObjectFit();
    TLorentzVector mp4 = fVec->mp4 ; 
    TLorentzVector np4 = fVec->np4 ;

    double MES =  par[1] ;
    double phi = atan2( np4.Py() , np4.Px() ) ;

    double px1 = ( np4.Pt()*cos( phi + par[0]) ) * MES;
    double py1 = ( np4.Pt()*sin( phi + par[0]) ) * MES;
    double nE1 = np4.E() * MES ;

    TLorentzVector bp4 = fVec->bp4 * MES;
    TLorentzVector np4a( px1, py1, 0, nE1 );

    //  ====== solving neutrino pz =======
    vector<TLorentzVector> nPs;

    double xW  = mp4.Px() + np4a.Px() ;
    double yW  = mp4.Py() + np4a.Py() ;
    double nuPt2 = ( np4a.Px()*np4a.Px() ) + ( np4a.Py()*np4a.Py() );

    double zl = mp4.Pz() ;
    double El = mp4.E()  ;

    double D = (80.4*80.4) - (El*El) - nuPt2 + (xW*xW) + (yW*yW) + (zl*zl) ;

    double A = 4.*( (zl*zl) - (El*El) ) ;
    double B = 4.*D*zl ;
    double C = (D*D) - (4.*El*El*nuPt2) ;

    double EtW  = mp4.Pt() + np4a.E() ;
    double PtW2 = (xW*xW) + (yW*yW);
    double MtW2 = (EtW*EtW) - PtW2 ;

    if ( MtW2 < 0. ||  (B*B) < (4.*A*C) ) {
       nPs.push_back( np4 );
    } else {

      double nz1 = (-1.*B + sqrt(B*B - 4.*A*C )) / (2.*A) ;
      double nz2 = (-1.*B - sqrt(B*B - 4.*A*C )) / (2.*A) ;
      for (int i=1; i<3; i++ ) {
          double nz = ( i == 1) ? nz1 : nz2 ;
	  double ENu2 = ( np4a.Px()*np4a.Px() ) + ( np4a.Py()*np4a.Py() ) + (nz*nz);
	  double zW = mp4.Pz() + nz ;
	  double EW = mp4.E()  + sqrt(ENu2) ;
	  TLorentzVector neu_p4 = TLorentzVector( np4a.Px(), np4a.Py(), nz, sqrt(ENu2) );
          nPs.push_back( neu_p4 ) ;
      }

    }
    // =================================================

    double dtmass = 999;
    for (size_t i =0; i< nPs.size(); i++ ) {
        double tpx = mp4.Px() + nPs[i].Px() + bp4.Px() ;
        double tpy = mp4.Py() + nPs[i].Py() + bp4.Py() ;
        double tpz = mp4.Pz() + nPs[i].Pz() + bp4.Pz() ;
        double tE  = mp4.E()  + nPs[i].E()  + bp4.E() ;
        TLorentzVector tp4( tpx, tpy, tpz, tE );
        double dM = fabs( tp4.M() - 172.5 ) ;
        if ( dM < dtmass && nPs.size() > 1 ) {
           dtmass = dM ;
        }
    }

    tChi2 =  dtmass*dtmass ;
}

// v1 = b jet p4,  v2 = w p4
vector<TLorentzVector> LepTopMassFitter::ReFitNeutrino( TLorentzVector v1, TLorentzVector mp4, Double_t *par ) {

   double MES =  par[1] ;
   double phi = atan2( v1.Py() , v1.Px() ) ;

   double px1 = ( v1.Pt()*cos( phi + par[0]) ) * MES;
   double py1 = ( v1.Pt()*sin( phi + par[0]) ) * MES;
   double nE1 = v1.E() * MES ;

   TLorentzVector np4a( px1, py1, 0, nE1 );
   vector<TLorentzVector> neuSolutions = NeuP4Solution( mp4, np4a );
    
   return neuSolutions ;
}


vector<TLorentzVector> LepTopMassFitter::NeuP4Solution( TLorentzVector muP4, TLorentzVector neuP4 ) {

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
	  double zW = muP4.Pz() + nz ;
	  double EW = muP4.E()  + sqrt(ENu2) ;
	  TLorentzVector np4 = TLorentzVector( neuP4.Px(), neuP4.Py(), nz, sqrt(ENu2) );
          neuSol.push_back( np4 ) ;
      }

    }

    return neuSol ;
}

vector<TLorentzVector> LepTopMassFitter::GetLorentzVector( jlist Ls, double jpx[], double jpy[], double jpz[], double jE[] ){

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

vector<TLorentzVector> LepTopMassFitter::GetLorentzVector( vector<TLorentzVector>& oblist, jlist Ls ){

    TLorentzVector wj1 = oblist[Ls.w1] ;
    TLorentzVector wj2 = oblist[Ls.w2] ;
    TLorentzVector bjh = oblist[Ls.bh] ;
    TLorentzVector bjl = oblist[Ls.bl] ;

    vector<TLorentzVector> tjets ;
    tjets.push_back( wj1 ) ;
    tjets.push_back( wj2 ) ;
    tjets.push_back( bjh ) ;
    tjets.push_back( bjl ) ;

    return tjets ;
}

// jets first , last one is muon
vector<TLorentzVector> LepTopMassFitter::GetEventObjects( int nj, double jpx[], double jpy[], double jpz[], double jE[], double mpx, double mpy, double mpz, double mE ) {

    vector<TLorentzVector> objs ;

    for ( int i=0; i< nj; i++) {
        TLorentzVector jet_i( jpx[i], jpy[i], jpz[i], jE[i] );
        objs.push_back( jet_i ) ;
    }
    TLorentzVector mp( mpx, mpy, mpz, mE );
    objs.push_back( mp ) ;

    return objs ;
}

double LepTopMassFitter::BTagProbability( double bCut[], int NofB, int jid[] ){

   double probb = 0 ;
   if ( n_btag ==2 && NofB >= 2 && bCut[ jid[2] ] >= bTh && bCut[ jid[3]] >= bTh ) probb = 1;
   if ( n_btag ==0 && NofB >= 2 && bCut[ jid[2] ] >= bTh && bCut[ jid[3]] >= bTh ) probb = 1;

   if ( n_btag < 2 && NofB == 1 && bCut[ jid[2] ] >= bTh  ) probb = 1;
   if ( n_btag < 2 && NofB == 1 && bCut[ jid[3] ] >= bTh  ) probb = 1;

   if ( n_btag >=0 && NofB == 0  ) probb = 0 ;
   if ( n_btag == -1 ) probb = 1 ;

   return probb ;
}


void LepTopMassFitter::ReFitLepTopSolution( string fileName, recoObj* wObj, double scale, vector<int>* evtlist, bool smearing, TTree* theTree ) {

  // get files and trees
  TTree* tr1 ;
  if ( theTree ==  NULL ) {
     tr1 = fitInput->GetTree( fileName, "muJets" );
  } else {
     tr1 = theTree ;
  }

  bool isMES = ( IsMES == "ON" ) ? true : false ;

  // event solution tree
  const static int sz = n_Jets + 1;
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

  // clean the permutation lists and fill the permuat
  vector<TLorentzVector> vlist0;
  vector<TLorentzVector> vlist1;
  vector<TLorentzVector> vlist2;
  vector<TLorentzVector> vlist3;
  vector<TLorentzVector> vlist4;
  vector<TLorentzVector> vlist5;

  int loopEnd = tr1->GetEntries() ;
  int idx = 0;
  if ( evtlist != NULL ) loopEnd = evtlist->size() ;

  for ( int j= 0 ; j< loopEnd; j++ ) {
      if ( evtlist != NULL && evtlist->size() == 0 ) break;
      if ( evtlist != NULL ) idx = (*evtlist)[j] ;
      if ( evtlist == NULL ) idx = j ;
      tr1->GetEntry(idx);

      // 0 ~ 3(4) jets, 4(5): muon, 5(6):MET
      vector<TLorentzVector> objlist = GetEventObjects( nJ, jpx, jpy, jpz, jE, mpx[0], mpy[0], mpz[0], mE[0] );
      if ( smearing ) pseudoExp->PhaseSmearing( objlist, j );

      vlist0.clear();
      vlist1.clear();
      vlist2.clear();
      vlist3.clear();
      vlist4.clear();
      vlist5.clear();

      // getting permutations list
      if ( nJ < n_Jets ) continue;
      vector<jlist> jlistV ;
      fitInput->GetPermutes( nJ, jlistV );

      int NofB = 0;
      for (int k=0; k< nJ; k++) {
          if ( bCut[k] > bTh ) NofB++ ;
      }

      // reject events without neutrino pz solution
      //cout<<" --- "<<endl;
      for (size_t i=0; i< jlistV.size(); i++) {
 
           int jid[4] = { jlistV[i].w1, jlistV[i].w2, jlistV[i].bh, jlistV[i].bl };
           //vector<TLorentzVector> tjets = GetLorentzVector( jlistV[i], jpx, jpy, jpz, jE );
           vector<TLorentzVector> tjets = GetLorentzVector( objlist, jlistV[i] );
           TLorentzVector  muP4( mpx[0], mpy[0], mpz[0], mE[0] );
           TLorentzVector neuP4( npx[0], npy[0], npz[0], nE[0] );
           size_t obsz = objlist.size() ;
           if ( smearing ) {
              muP4  = objlist[obsz-2] ;
              neuP4 = objlist[obsz-1] ;
           }
           //cout<<" Permus = ( "<< jid[0] <<","<< jid[1] <<","<<jid[2] <<","<<jid[3] <<" ) " <<endl;

           // btagging 
           double probb = BTagProbability( bCut, NofB, jid );
           if ( probb == 0 ) continue ;

           // MET adjusting 
           Double_t para[2] = {0,1};
           Double_t errs[2];
           if ( isMES ) FitLepTop( muP4, neuP4, tjets[3], para, errs, true );
           vector<TLorentzVector> neuSols = ReFitNeutrino( neuP4, muP4, para );
           //cout<<" neu sol size = " << neuSols.size() ;
           //cout<<"  p0 = "<< para[0] <<" p1 = "<< para[1] <<endl;
           //cout<<" n0 ( "<< neuP4.Px() <<", "<< neuP4.Py() <<", "<< neuP4.Pz() <<", "<< neuP4.E() <<" )"<<endl ;
           if ( neuSols.size() != 2 ) continue;


           TLorentzVector hadTP  = tjets[0] +  tjets[1]  + tjets[2] ;
           TLorentzVector lepTP0 = tjets[3] + neuSols[0] + muP4  ;
           TLorentzVector lepTP1 = tjets[3] + neuSols[1] + muP4  ;

           // ensure both dM(Had-Lep) < 30 GeV
           double dMHL0 = fabs( lepTP0.M() - hadTP.M() );
           double dMHL1 = fabs( lepTP1.M() - hadTP.M() );
           if ( dMHL0 > 30. && dMHL1 > 30. ) continue; 

           // choose the neutrino pz by looking at leptonic top mass
           double dML0 = ( fabs(lepTP0.M() - 172.5) > 30. || fabs(neuSols[0].Pz()) > 500. ) ? 999 : fabs(lepTP0.M() - 172.5) ;
           double dML1 = ( fabs(lepTP1.M() - 172.5) > 30. || fabs(neuSols[1].Pz()) > 500. ) ? 999 : fabs(lepTP1.M() - 172.5) ;
           if ( dML0 == 999. && dML1 == 999. ) continue;
           TLorentzVector theNeuSol = ( dML0 < dML1 ) ? neuSols[0] : neuSols[1] ;

           // exclude the bad pz solutions
           /*
           for (size_t j =0; j< neuSols.size(); j++) {
               if ( fabs( neuSols[j].Pz() ) > 500. ) neuSols.erase( neuSols.begin()+j );
           }
           if ( neuSols.size() == 0 ) continue;
           if ( neuSols.size() == 1 ) theNeuSol = neuSols[0] ;
           */

           vlist0.push_back( tjets[0] );
	   vlist1.push_back( tjets[1] );
	   vlist2.push_back( tjets[2] );
	   vlist3.push_back( tjets[3] );
	   vlist4.push_back( muP4  );
	   vlist5.push_back( theNeuSol  );

      }
      double weighting = 1. / vlist0.size() ;
      //double weighting = scale  ;
      for (size_t i=0; i< vlist0.size(); i++) {
           wObj->gethad( vlist0[i], vlist1[i], vlist2[i] );
           wObj->getlep( vlist3[i], vlist4[i], vlist5[i] );
           wObj->Fillh( weighting, scale );
      }
  }

}

void LepTopMassFitter::LepTopMCSolution( string fileName,  recoObj* wObj ) {

  // get files and trees
  TTree* tr0 = fitInput->GetTree( fileName, "mcmTt" );
  TTree* tr1 = fitInput->GetTree( fileName, "muJets" );

  bool isMES = ( IsMES == "ON" ) ? true : false ;

  // event solution tree
  int evtId0, jid0[6] ;
  tr0->SetBranchAddress( "evtId" ,&evtId0 );
  tr0->SetBranchAddress( "pId"   ,&jid0   );

  // event solution tree
  const static int sz = 5;
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

  vector<double> bCutV ;

  int idx = 0 ;
  for ( int j= 0 ; j< tr0->GetEntries() ; j++ ) {
      tr0->GetEntry(j);

      for ( int i=idx; i< tr1->GetEntries() ; i++) {
          tr1->GetEntry(i);
          if ( evtId == evtId0  ) {
             idx = i ;

             int NofB = 0;
             bCutV.clear() ;
             for (int k=0; k< nJ; k++) {
                 bCutV.push_back( bCut[k] );
                 if ( bCut[k] > bTh ) NofB++ ;
             }

             jlist mclist;
             mclist.w1 = jid0[0] ;
             mclist.w2 = jid0[1] ;
             mclist.bh = jid0[2] ;
             mclist.bl = jid0[3] ;
             vector<TLorentzVector> tjets = GetLorentzVector( mclist, jpx, jpy, jpz, jE );
             TLorentzVector  muP4( mpx[0], mpy[0], mpz[0], mE[0] );
             TLorentzVector neuP4( npx[ jid0[4]-1 ], npy[ jid0[4]-1 ], npz[ jid0[4]-1 ], nE[ jid0[4]-1 ] );

             Double_t para[2] = {0,1};
             Double_t errs[2];
             if ( isMES ) FitLepTop( muP4, neuP4, tjets[3], para, errs, true );
             vector<TLorentzVector> neuSols = ReFitNeutrino( neuP4, muP4, para );
             double dR = 99. ;
             int nIdx = -1 ;
             for (size_t k =0; k< neuSols.size(); k++)  {
                 double dR1 = neuP4.DeltaR( neuSols[k] );
                 if ( dR1 < dR ) { 
                    dR = dR1 ;
                    nIdx = k ;
                 }
             }

             //double probb = BTagProbability( bCutV, NofB, jid0, isJES );
             //if ( n_btag > -1 && probb == 0 ) continue ;

             wObj->gethad( tjets[0], tjets[1], tjets[2] );
             wObj->getlep( tjets[3], muP4, neuSols[nIdx]);
             wObj->Fillh( 1., 1. );
          }
      }
  }

}

