#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>

class TemplateMassFit : public TObject {

private:

   TFile* file0;
   TFile* file1;
   TFile* file2;
   TFile* file3;
   TH2D*  ttmass;
   TH2D*  ttMC;
   TH1D*  ttmass_sg;
   TFile* file161;
   TFile* file166;
   TFile* file171;
   TFile* file176;
   TFile* file181;

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
 
   TCanvas* c1;
   TCanvas* c2;
   TCanvas* c3;
   TCanvas* c4;
   TCanvas* c5;
   TCanvas* c6;

public:

   TemplateMassFit( TString aname, TString channel );     

   void showAll(int rbin);
   void TemplateFitting( int rbin, int lowBound, int upBound );
   void MultiTemplatesFitting( int rbin, int lowBound, int upBound );
   void MoreCombinedFitting( int rbin, int lowBound, int upBound );
   void CombinedFitting( int rbin, int lowBound, int upBound );
   void ScaleTemplates( double factor, TH1D* tmp, int B1, int B2 );
   void SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, Double_t *pars );
   double TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *bPred);
   double TemplateTest( TString mName, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sPred, Double_t *tPred, Double_t *bPred );
   double Chi2Test( TString mName1, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sgpar, Double_t *bgpar);
   double Chi2Test( TString mName1, TH1D* theData, int rbin, int lowBound, int upBound, Double_t *sgpar, Double_t *tbpar, Double_t *bgpar );
   void TemplateDrawer( TString mName, int rbin, int lowBound, int upBound, Bool_t *comp );

   void combineBG( TString mName, TH1D* allbg, int rbin );
   void getSignal( TH1D* h_Sg, int rbin, TString mName );
   void getBackground( TH1D* h_Bg, int type, int rbin, TString mName );
   void getData(TH1D* h_data, int rbin );
   void getFakeData(TH1D* h_data, THStack* ttstk, TH1D* dth0, TH1D* dth1, TH1D* dth2, TH1D* dth3, int rbin );
   void getFakeData(TH1D* h_data, int rbin );
  
   Double_t fitG( Double_t* x, Double_t* par);
   Double_t fitLG( Double_t* x, Double_t* par);
   Double_t fitGS( Double_t* x, Double_t* par);
   Double_t fitTL( Double_t* x, Double_t* par);
   Double_t fitLD( Double_t* x, Double_t* par);
   Double_t fitSG( Double_t* x, Double_t* par);
   Double_t fitBW( Double_t* x, Double_t* par);
   Double_t fitData( Double_t* x, Double_t* par);
   Double_t ConvBWGS( Double_t* x, Double_t* par);
   Double_t ConvSGGS( Double_t* x, Double_t* par);

   ClassDef(TemplateMassFit, 1);

};

#if !defined(__CINT__)
    ClassImp(TemplateMassFit);
#endif

