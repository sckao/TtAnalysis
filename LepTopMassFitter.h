#ifndef LepTopMassFitter_H
#define LepTopMassFitter_H

#include <TObject.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "MassAnaInput.h"
#include "WFormat.h"
#include "PseudoExp.h"

class LepFitVectors : public TObject {

public:

   LepFitVectors(){}
   ~LepFitVectors(){}

   TLorentzVector mp4;
   TLorentzVector np4;
   TLorentzVector bp4;

  // ClassDef(FitVectors, 1);
};

class LepTopMassFitter : public TObject {

private:


   double tmass;

   double bTh;
   int    n_btag;
   int    n_Jets;
   string IsMES ;
   string hfolder;

   MassAnaInput*    fitInput;
   PseudoExp*       pseudoExp;


public:

   LepTopMassFitter();     
   ~LepTopMassFitter();     
 
   void FitLepTop( TLorentzVector mP4, TLorentzVector nP4, TLorentzVector bP4, Double_t* pars, Double_t* errs, bool isJES );

   static void LTMFCN(Int_t &npar, Double_t *, Double_t &tChi2, Double_t *par, Int_t iflag);

   vector<TLorentzVector> GetLorentzVector( jlist Ls, double jpx[], double jpy[], double jpz[], double jE[] );

   double BTagProbability( double bCut[], int NofB, int jid[] );

   vector<TLorentzVector> GetEventObjects( int nj, double jpx[], double jpy[], double jpz[], double jE[], double mpx, double mpy, double mpz, double mE );

   vector<TLorentzVector> GetLorentzVector( vector<TLorentzVector>& oblist, jlist Ls );

   vector<TLorentzVector> ReFitNeutrino( TLorentzVector v1, TLorentzVector muP4, Double_t *par );

   vector<TLorentzVector> NeuP4Solution( TLorentzVector muP4, TLorentzVector neuP4 );

   void ReFitLepTopSolution( string fileName, recoObj* wObj, double scale = 1., vector<int>* evtlist = NULL, bool smearing = false, TTree* theTree = NULL );

   void LepTopMCSolution( string fileName, recoObj* wObj );

   //ClassDef(LepTopMassFitter, 1);

};


//#if !defined(__CINT__)
//    ClassImp(LepTopMassFitter);
#endif

