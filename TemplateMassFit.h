#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>

class TemplateMassFit : public TObject {

private:

   FILE* logfile;
   FILE* parfile;

   TFile* file0;
   TFile* file1;
   TFile* file2;
   TFile* file3;
   TH2D*  ttmass;
   TH2D*  ttMC;
   TH1D*  ttmass_sg;

   TFile* thefile;

   TString aname;
   TString channel;
   TString hname;
   TString sname;
   TString cname;
   TString fname;
   TString hfolder;

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

public:

   TemplateMassFit( TString aname, TString channel );     

   void MultiTemplatesFitting( int rbin, int lowBound, int upBound );
   void TemplateFitting( int rbin, int lowBound, int upBound );
   double TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *bPred);
   double TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *tPred, Double_t *bPred );

   void MoreCombinedFitting( int rbin, int lowBound, int upBound );
   void CombinedFitting( int rbin, int lowBound, int upBound );

   double Chi2Test(  TH1D* theData, int lowBound, int upBound, int nPar ); 
   double getChi2( TH1D* theData,  TF1* theFunc, TF1* fS, TF1* fB, TF1* fW, double lowBound , double upBound, int nPar );

   void FitTester( TString mName, int rbin, int lowBound, int upBound );
   void TemplateDrawer( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp );
   void SetFitParameters( double mass, Double_t* para, int nPara );

   vector<TString> FillMassAssumption( int npoints );
   double MassDigi( TString mString );  

   void ScaleTemplates( double factor, TH1D* tmp, int B1, int B2 );
   void ScaleTemplates( double lumi, double nEvents, int channel, TH1D* tmp );
   void SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, Double_t *pars, Double_t *perr = NULL );

   void combineBG( TString mName, TH1D* allbg, int rbin );
   void getSignal( TH1D* h_Sg, int rbin, TString mName );
   void getBackground( TH1D* h_Bg, int type, int rbin, TString mName );
   void getData(TH1D* h_data, int rbin );
   void getFakeData(TH1D* h_data, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, int rbin );
   void getFakeData(TH1D* h_data, int rbin );

   ClassDef(TemplateMassFit, 1);

};

#if !defined(__CINT__)
    ClassImp(TemplateMassFit);
#endif

