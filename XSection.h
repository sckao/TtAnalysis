#ifndef XSection_H
#define XSection_H

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

#include "MassFitFunction.h"
#include "MassAnaInput.h"
#include "HadWMassFitter.h"
#include "PseudoExp.h"
#include "BgEstimation.h"
#include "ObjectInfo.h"


//#include "LepTopMassFitter.h"

class XSection : public TObject {

private:

   MassAnaInput*      fitInput;
   HadWMassFitter*    wmfitter;
   PseudoExp*         pseudoExp;
   BgEstimation*      bgEst;
   ObjectInfo*        objInfo ;

   string hfolder;
   TString theFolder;
   bool smearing ;
   vector<double> inputMean ;

   vector<string> File4J;
   vector<string> File2J;
   vector<string> File0J;
   vector<string> dataFile;

public:

   XSection();     
   ~XSection();     
 
   vector<double> Ratio42( int statIdx = 0 );

   void MethodTest() ;

   vector<double> GetXSection( vector<double>& R42, double n2J, double n4J ) ;
   void EnsembleTest( int nRun, int randomSeed );
   vector<double> RunEnsembles( int nJets, int treeSize, double pMean, int nRun, int randomSeed, string fileName, TH1D* hGen = NULL );

   void MethodTest1() ;
   void BgClosureTest( int nX, int nY, bool inclX = false, bool inclY = false ) ;
   vector<double> CutEff();
   vector<double> CrossSection( double nData_4J, vector<double>& nData_2j, vector<double>& R42, vector<double>& EffCut );
   void EnsembleTest1( int nRun, int randomSeed );
   vector< vector<double> >RunEnsembles1( int nJets, int treeSize, double pMean, int nRun, int randomSeed, string fileName, TH1D* hGen = NULL );

   vector<double> ErrAxB( double A, double s_A, double B, double s_B );

   void RealBackground( int nX, int nY, bool inclX = false, bool inclY = false ) ;

   //ClassDef(XSection, 1);

};

//#if !defined(__CINT__)
//    ClassImp(XSection);
#endif

