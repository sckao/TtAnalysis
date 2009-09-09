#include "TObject.h"
#include "MassFitFunction.h"
#include "MassFitInput.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TH2.h>
#include <TH1.h>
#include <TH1D.h>
#include <THStack.h>
#include <TF1.h>
#include <TFormula.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFractionFitter.h>

class TemplateMassFit : public TObject {

private:

   //FILE* logfile;
   FILE* parfile;
   FILE* errfile;
   FILE* Bgpara;
   FILE* Sgpara;

   TString channel;
   TString cname;
   TString hfolder;
   TString ptype;

   TString plot1;
   TString plot2;
   TString plot3;
   TString plot4;
   TString plot5;
   TString plot6;
   TString plot7;
   TString plot8;
   TString plot9;
 
   TCanvas* c1;
   TCanvas* c2;
   TCanvas* c3;
   TCanvas* c4;
   TCanvas* c5;
   TCanvas* c6;
   TCanvas* c7;
   TCanvas* c8;
   TCanvas* c9;

   MassFitFunction* fitFunc;
   MassFitInput*    fitInput;

public:

   TemplateMassFit( TString channel, int NBTag );     
   ~TemplateMassFit();

   void MultiTemplatesFitting( TString mName, int rbin, int lowBound, int upBound );
   void TemplateFitting( TString mName, int rbin, int lowBound, int upBound );
   double TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *bPred);
   double TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *tPred, Double_t *bPred );

   void MoreCombinedFitting( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp, int NBTag = -1 );
   void CombinedFitting( TString mName, int rbin, int lowBound, int upBound, int NBTag );

   double Chi2Test( TString mName, TH1D* theData, int lowBound, int upBound, int nPar, int NBTag = -1, Double_t* statErr=NULL, int rbin = 10 ); 
   double getChi2( TH1D* theData,  TF1* theFunc, TF1* fS, TF1* fB, TF1* fW, double lowBound , double upBound, int nPar );

   void FitSignal( TString mName, int rbin, int lowBound, int upBound, Double_t *para = NULL, Double_t *perr = NULL );
   void FitBackground( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp, Double_t *para=NULL, Double_t *perr=NULL );
   void GetAllCoeff( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp );
   void TemplateDrawer( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp );
   void SetFitParameters( double mass, Double_t* para, int nPara, int NBTag, int rbin=10 );

   vector<TString> FillMassAssumption( int npoints );
   double MassDigi( TString mString );  

   void SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, Double_t *pars, Double_t *perr = NULL );
   void ScaleTemplates( double factor, TH1D* tmp, int B1, int B2 );
 
   ClassDef(TemplateMassFit, 1);

};

#if !defined(__CINT__)
    ClassImp(TemplateMassFit);
#endif

