#ifndef HadWMassFitter_H
#define HadWMassFitter_H

#include <TObject.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMinuit.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "MassAnaInput.h"
#include "WFormat.h"
#include "PseudoExp.h"


class FitVectors : public TObject {

public:

   FitVectors(){}
   ~FitVectors(){}

   TLorentzVector jv1;
   TLorentzVector jv2;

   //ClassDef(FitVectors, 1);
} ;


class HadWMassFitter : public TObject {

private:

   vector<double> jetCuts;
   vector<double> muonCuts;
   vector<double> M2M3Cut;
   double LepM2tCutL ;
   double dM3Cut ;
   double JES ;

   double sigma1;
   double sigma2;
   double zigma1;
   double zigma2;

   double wmass;

   double bTh;
   int    n_btag;
   int    JESType;

   string Inclusive ;
   bool inclu ;
   bool normMCData ;

   MassAnaInput*    fitInput;
   PseudoExp*       pseudoExp;

   TCanvas* c1;
   TCanvas* c2;

   string hfolder;

public:

   HadWMassFitter();     
   ~HadWMassFitter();     

   void FitW( TLorentzVector jetv1, TLorentzVector jetv2, Double_t* pars, Double_t* errs, bool isJES );

   static void WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag);

   double ReFitWMass( TLorentzVector v1, TLorentzVector v2, Double_t *par, Double_t *ProbW );

   double ReFitHadTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par );

   double ReFitLepTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par );

   double KinematicProb( vector<TLorentzVector> nvlist, TLorentzVector mV4, TLorentzVector nV4, Double_t *par, Double_t* wtX = NULL );

   double Get2BodySigma( TLorentzVector v1, TLorentzVector v2 );

   double Get3BodySigma( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3 );

   double BTagProbability( double bCut[], int NofB, int jid[], bool isJES = false );

   vector<TLorentzVector> newVector( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, TLorentzVector v4, Double_t *par );

   double jetAngle( TLorentzVector v1, TLorentzVector v2 );
   double minTheta_WDecay( TLorentzVector v1, TLorentzVector v2 );

   vector<TLorentzVector> NeuP4Solution( TLorentzVector muP4, TLorentzVector neuP4 ) ;

   vector<TLorentzVector> GetLorentzVector( jlist Ls, double jpx[], double jpy[], double jpz[], double jE[] );

   vector<TLorentzVector> GetLorentzVector( vector<TLorentzVector>& oblist, jlist Ls );

   vector<TLorentzVector> GetEventObjects( int nj, double jpx[], double jpy[], double jpz[], double jE[], double mpx, double mpy, double mpz, double mE );

   void ResetCuts( double m2L = -1, double m2H = -1, double m3L = -1, double m3H = -1, double lepM2tL = -1, double lepM2tH = -1, bool GetDefault = false );

   void SetMCNormalization( bool normMC ) ;
   void SetMuonCuts( double pt_, double eta_, double iso_ );

   // methods for old soltree
   double ReFitSolution( string mName, recoObj* wObj, int nJets, double scale = 1., vector<int>* evtlist = NULL, int evtSplit = 0, bool smearing = false );

   void MCSolution( string fileName, recoObj* wObj );

   double EvtScaling( int NJets, string fileName );
   //ClassDef(HadWMassFitter, 1);

};


//#if !defined(__CINT__)
//    ClassImp(HadWMassFitter);
#endif

