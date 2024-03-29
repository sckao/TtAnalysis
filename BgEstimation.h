#ifndef BgEstimation_H
#define BgEstimation_H

#include <TObject.h>
#include <TMinuit.h>
#include <TFractionFitter.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <iomanip>

#include "MassAnaInput.h"
#include "MassFitFunction.h"
#include "ObjectInfo.h"
#include "WFormat.h"

class BgEstimation : public TObject {

private:

   double mL;
   double mH;

   // for hadronic permutation

   MassAnaInput*    fitInput;
   ObjectInfo*      objInfo;

   string hfolder;
   string plotType ;
   bool smearing ;

public:

   BgEstimation();     
   ~BgEstimation();     

   vector<double> RatioXY(int njX, int njY, vector<string>& fNames, int index, bool doPlot = false, bool inclX = false, bool inclY = false, bool normMC = false );
   vector<double> RatioXY( vector<string>& fNameX, vector<string>& fNameY, int index = -1., bool includeTt = false, bool doPlot = false ) ;

   vector<double> RatioXY( vector<string>& fNames, vector<double>& scaleX, vector<double>& scaleY, int njX, int njY, bool inclX = true, bool inclY = false );

   vector<double> RatioX2( vector<string>& MCFiles, vector<string>& DataFiles, vector<double>& scaleV2, vector<double>& scaleV1 ) ;
   vector<double> RatioX2( vector<string>& MCFiles, vector<string>& DataFiles ) ;
   vector<double> VQFraction( vector<string>& fNames, vector<double>& scaleV ) ;
   vector<double> VQFraction( vector<string>& fNames ) ;
   vector<double> AlphaSFromR21( vector<string>& fNames, vector<double>& scaleV2, vector<double>& scaleV1 ) ;

   void RatioForSystematic( string cfgFile ) ;

   vector<double> BgEstimate( vector<double>& R42, vector<double>& nData_2j );
   vector<double> BgEstimate( vector<double>& R42, double nData_2j );

   void LepM2Plotter( vector<TH1D*>& h1Ds, TString Tag  );

   void RatioPlotter( vector<double>& nW1, vector<double>& nW2, TString Tag, double scale = 1.  );

   //void RatioForSystem( int nx, int ny, bool inclX = false, bool inclY = false );

   vector<double> MtFitter( string& DataName, vector<string>& fNames, int phaseIdx = 0 );
   vector<double> JacobianFitter( TH1D* hData, vector<TH1D*>& hTemplates, int phaseIdx = 0, int fbin1 = 1, int fbin2 = 20);   

   vector<double> MeasureScale( string& dataFiles, vector<string>& fNames );

   double MeasureScale2D( string& dataFiles, vector<string>& fNames, double MtCut = 40, double METCut = 25, bool doQCD = true, bool isVjNorm = false );
   double MeasureScale1D( string& dataFiles, vector<string>& fNames, int hID, double MtCut = 40, double METCut = 25, bool doQCD = true, bool isVjNorm = false, int njets = 1 );
   void MeasureScale1D( string& dataFiles, string cfgFile, int hID, double MtCut = 40, double METCut = 25, bool doQCD = true, bool isVjNorm = false, int njets = 1 );
   double MeasureScale1D( vector<string>& fNames, int hID, double MtCut = 40, double METCut = 25, bool doQCD = true, bool isVjNorm = false );

   double Chi2Normalization( TH2D* hData, TH2D* hMC, int Bx1, int Bx2, int By1, int By2, double s1, double s2, bool doFit = false );
   double Chi2Normalization( TH1D* hData, TH1D* hMC, int Bx1, int Bx2, double s1, double s2, TString plotname, bool doFit = false  );

   void RatioScan( vector<string>& DataNames, vector<string>& fNames,  bool normMC );
   //ClassDef(BgEstimation, 1);
  
};

//#if !defined(__CINT__)
//    ClassImp(BgEstimation);
#endif

