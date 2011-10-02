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
#include <TGraphAsymmErrors.h>
#include <TGraphPainter.h>
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
   string plotType ;

   bool smearing ;
   bool mcNorm ;
   vector<double> inputMean4J ;
   vector<double> inputMean2J ;

   vector<string> File4J;
   vector<string> File2J;
   vector<string> File0J;
   vector<string> dataFile;

public:

   XSection();     
   ~XSection();     
 

   void BgClosureTest( int nX, int nY, bool inclX = false, bool inclY = false ) ;

   vector<double> CutEff( int nj = 4 );
   void CutEff( string cfgFile, int nj = 4 );

   void MethodTest1( int nj = 4 ) ;
   void MethodTest2( string cfgFile, int nj ) ;
   void MethodTest3( int nj = 4 ) ;

   vector<double> CrossSection( double nData_4J, vector<double>& nData_2j, vector<double>& R42, vector<double>& EffCut );
   vector<double> CrossSection( double nData_4J, double nData_2j, vector<double>& R42, vector<double>& EffCut );

   void EnsembleTest1( int nRun, int randomSeed );

   vector<double> RunEnsembles( int nJets, int treeSize, double pMean, int nRun, int randomSeed, string fileName, TH1D* hGen = NULL );
   vector<double> RunEnsembles1( int nJets, int treeSize, double pMean, int nRun, int randomSeed, string fileName, TH1D* hGen = NULL );

   vector<double> ErrAxB( double A, double s_A, double B, double s_B );

   void RealDataAnalysis( int nX = 4, int nY = 2) ;
   void DataForSystematic( string cfgFile, int nX = 4 ) ;

   //ClassDef(XSection, 1);

};

//#if !defined(__CINT__)
//    ClassImp(XSection);
#endif

