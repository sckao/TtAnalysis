#ifndef __MassTemplate__
#define __MassTemplate__

#include "../root/include/TObject.h"

class MassTemplate : public TObject {

private:

   TFile* file0;
   TFile* file1;
   TFile* file2;
   TFile* file3;
   TH2D*  ttmass;
   TH2D*  ttMC;
   TH1D*  ttmass_sg;

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

   MassTemplate() ={ aname = "zero" ; channel = "Had"; };     
   MassTemplate( TString aname, TString channel );     

 
   void showAll(int rbin);
   void TemplateFitting( int rbin, int lowBound, int upBound );
   void TemplateFitting2( int rbin, int lowBound, int upBound );
   void CombinedFitting( int rbin, int lowBound, int upBound );
   void ScaleTemplates( double factor, TH1D* tmp, int B1, int B2 );

   void combineBG( TH1D* allbg, int rbin );
   void getSignal( TH1D* h_Sg, int rbin );
   void getBackground( TH1D* h_Bg, int type, int rbin  );
   void getData(TH1D* h_data, int rbin );
   void getFakeData(TH1D* h_data, THStack* ttstk, int rbin );
   void getFakeData(TH1D* h_data, int rbin );
   void SmoothTemplate( int type, TH1D* ds, TH1D* ds1, int lowBound, int upBound, int npar );
  
   Double_t fitG( Double_t* x, Double_t* par);
   Double_t fitGS( Double_t* x, Double_t* par);
   Double_t fitBW( Double_t* x, Double_t* par);
   Double_t fitData( Double_t* x, Double_t* par);
   

   ClassDef(MassTemplate, 1);

};
#endif
