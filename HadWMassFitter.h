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
#include "MassAna.h"
#include "MassFitFunction.h"
#include "WFormat.h"

struct jlist {
    int w1 ;
    int w2 ;
    int bh ;
    int bl ;
};


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

   double mL;
   double mH;

   double sigma1;
   double sigma2;
   double zigma1;
   double zigma2;

   double wmass;

   double bTh;
   int    n_btag;

   MassAnaInput*    fitInput;
   MassAna*         fitTools;
   MassFitFunction* fitFunc;
 
   TCanvas* c1;
   TCanvas* c2;

   string hfolder;

public:

   HadWMassFitter( double massL, double massH );     
   ~HadWMassFitter();     

   void FitW( TLorentzVector jetv1, TLorentzVector jetv2, Double_t* pars, Double_t* errs, bool isJES );

   static void WMFCN(Int_t &npar, Double_t *, Double_t &wChi2, Double_t *par, Int_t iflag);

   double ReFitWMass( TLorentzVector v1, TLorentzVector v2, Double_t *par, Double_t *ProbW );

   double ReFitHadTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par );

   double ReFitLepTopMass( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, Double_t *par );

   double KinematicProb( int jid[] , vector<TLorentzVector> nvlist, TLorentzVector mV4, TLorentzVector nV4, Double_t *par, Double_t* wtX = NULL );

   double Get2BodySigma( TLorentzVector v1, TLorentzVector v2 );

   double Get3BodySigma( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3 );

   double BTagProbability( vector<double> bCut, int NofB, int jid[], bool isJES = false );

   vector<TLorentzVector> newVector( TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, TLorentzVector v4, Double_t *par );

   bool JetPtFilter( vector<TLorentzVector> jpv );

   Double_t jetAngle( TLorentzVector v1, TLorentzVector v2 );

   // methods for old soltree
   void ReFitSolution( string mName, recoObj* wObj, int type, bool isMCMatched = false );
   void ReFitSolution( string mName, recoObj* wObj, vector< pair<int,int> >& evtlist, int type );

   void ReFitSolution1( string mName, int njets, recoObj* wObj, vector<int>* evtlist = NULL );

   void ReFitSolution1( string mName, int njets, recoObj* wObj, int type , vector<int>* evtlist = NULL );

   void MCSolution( string fileName, recoObj* wObj, int type );

   void GetPermutes( int njets, vector<jlist>& jlistV );

   vector<TLorentzVector> GetLorentzVector( jlist Ls, double jpx[], double jpy[], double jpz[], double jE[] );

   //ClassDef(HadWMassFitter, 1);

};


//#if !defined(__CINT__)
//    ClassImp(HadWMassFitter);
#endif

