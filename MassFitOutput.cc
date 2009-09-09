#include "MassFitOutput.h"

MassFitOutput::MassFitOutput( TString channel, int NBTag ){
 

  fitfunc  = new MassFitFunction();
  fitter   = new TemplateMassFit( channel, NBTag ); 
  fitInput = new MassFitInput( channel, NBTag );
  fitInput->Initialize( &hfolder ); 
  ptype = ".gif";
  ch_name = channel;
 

}

MassFitOutput::~MassFitOutput(){

  delete fitfunc;
  delete fitter; 
  delete fitInput;

  delete pfile;
  delete efile;

}


void MassFitOutput::CoeffCalib( int rbin, int lowBound, int upBound, Bool_t* comp ) {


     fitter->GetAllCoeff("161", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("163", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("166", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("168", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("171", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("173", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("176", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("178", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("181", rbin, lowBound, upBound, comp) ;
     fitter->GetAllCoeff("183", rbin, lowBound, upBound, comp) ;

     FILE *ffile = fopen(hfolder+"extrf.log","a");

     TCanvas *c1 = new TCanvas("c1","", 700, 600);
     c1->SetFillColor(10);
     c1->SetFillColor(10);
     gPad->SetGridx();
     gPad->SetGridy();
     gStyle->SetOptStat(kTRUE);
     gStyle->SetOptFit(111);
     c1->cd();

     Double_t inputM[10] = { 161.2, 163.7, 166.2, 168.7, 171.2, 173.7, 176.2, 178.7, 181.2, 183.7 };
     Double_t pa[10] = {0.};
     Double_t er[10] = {0.};
     Double_t ex[10] = {0.};

     TF1 *func = new TF1("func",MassFitFunction::fitPoly, 150,190, 2);

     TString plot, title_plot;
     float f1, e1;
     for ( int k= 0; k<12; k++) {

         pfile = fopen(hfolder+"paraf.log","r");
         efile = fopen(hfolder+"perrf.log","r");
         title_plot = GiveParTitle( k );
         plot = title_plot+"_Scan.gif" ;

         for ( int i=0 ; i < 10;  i++) {
             for ( int j=0 ; j< 12; j++) {
                 int fo = fscanf(pfile, "%f" , &f1 );
                 int fe = fscanf(efile, "%f" , &e1 );
                 if ( fo != 1 || fe != 1) cout<<" reading error!!  par:"<<fo<<"  pErr:"<<fe<<endl;
                 if ( j != k )  continue;
                 pa[i] = f1;
                 er[i] = e1;
             }
         }
          
         TGraphErrors* hMass = new TGraphErrors( 10, inputM, pa, ex, er );
	 hMass->SetMarkerColor(4);
	 hMass->SetMarkerStyle(20);
	 hMass->SetLineWidth(2);
	 hMass->SetLineColor(4);

         hMass->SetTitle(title_plot);
         hMass->Draw("AP");
	 double a0 = hMass->GetMean();
	 double a1 = 1.;
	 hMass->Fit("func", "R", "sames", 150., 185);
	 a0 = func->GetParameter(0);
	 a1 = func->GetParameter(1);
	 cout<<" a0= "<<a0 <<"   a1= "<<a1 <<endl;
	 fprintf(ffile," %d  %.3f  %.3f \n", k, a0, a1 );

         c1->Update();
         c1->Print(hfolder+plot);
         delete hMass;
         fclose(pfile);
         fclose(efile);
     }
     
     fclose(ffile);      
 
     delete c1;
     delete func;
}

void MassFitOutput::MassCalib( int rbin, int lowBound, int upBound, Bool_t *comp, int NBTag, int NPara) {

     TString hName = "MassReco_"+ch_name;


     if ( NPara == 9 ) {
        fitter->CombinedFitting( "161", rbin, lowBound, upBound, NBTag ) ; 
	fitter->CombinedFitting( "166", rbin, lowBound, upBound, NBTag ) ; 
	fitter->CombinedFitting( "171", rbin, lowBound, upBound, NBTag ) ; 
	fitter->CombinedFitting( "176", rbin, lowBound, upBound, NBTag ) ; 
	fitter->CombinedFitting( "181", rbin, lowBound, upBound, NBTag ) ; 
     }

     if ( NPara == 12 ) {
        fitter->MoreCombinedFitting( "161", rbin, lowBound, upBound, comp, NBTag ) ; 
	fitter->MoreCombinedFitting( "166", rbin, lowBound, upBound, comp, NBTag ) ; 
	fitter->MoreCombinedFitting( "171", rbin, lowBound, upBound, comp, NBTag ) ; 
	fitter->MoreCombinedFitting( "176", rbin, lowBound, upBound, comp, NBTag ) ; 
	fitter->MoreCombinedFitting( "181", rbin, lowBound, upBound, comp, NBTag ) ; 
     }

     FILE* fitlog = fopen(hfolder+"Outputf.log","r");

     float inputM[5]  = { 161.2, 166.2, 171.2, 176.2, 181.2 };
     float recoM[5]   = { 161.2, 166.2, 171.2, 176.2, 181.2 };
     float statErr[5] = { 0. };
     float xErr[5]    = { 0. };

     float ipm, rcm, ste;
     int fo1, fo2, fo3 ;
     for ( int i=0; i< 5; i++) {
         fo1 = fscanf(fitlog, "%f", &ipm );
         fo2 = fscanf(fitlog, "%f", &rcm );
         fo3 = fscanf(fitlog, "%f", &ste );
         inputM[i] = ipm ;
         recoM[i]  = rcm ;
         statErr[i]= ste ;
         if ( fo1 != 1 ) cout<<" Input Mass Reading Error => "<<fo1<<endl;
         if ( fo2 != 1 ) cout<<" Reco Mass Reading Error => " <<fo2<<endl;
         if ( fo3 != 1 ) cout<<" Stat Error Reading Error => "<<fo3<<endl;
     }

     TCanvas *c1 = new TCanvas("c1","", 700, 600);
     c1->SetFillColor(10);
     c1->SetFillColor(10);
     gPad->SetGridx();
     gPad->SetGridy();
     c1->cd();

     gStyle->SetOptStat(kTRUE);
     gStyle->SetOptFit(111);

     TGraphErrors* hMass = new TGraphErrors( 5, inputM, recoM, xErr, statErr );
     hMass->SetMarkerColor(4);
     hMass->SetMarkerStyle(20);
     hMass->SetLineWidth(2);
     hMass->SetLineColor(4);
     hMass->SetTitle("Mass Reconstruction");
     hMass->Draw("AP");

     TF1* func0 = new TF1("func0", MassFitFunction::fitPoly , 150, 190, 2);
     hMass->Fit( func0, "RQ","sames", inputM[0]-5, inputM[4]+5 );

     c1->Update();
     c1->Print(hfolder+hName+ptype);

     fclose(fitlog);

     delete c1;
     delete hMass;
     delete func0;
}

TString MassFitOutput::GiveParTitle( int id ) {

    TString hTitle[12] = { "Gaus_P0", "Gaus_P1", "Gaus_P2", "LogNorm_P3", "LogNorm_P4", "LogNorm_P5",
                           "LandauTt_P6", "LandauTt_P7", "LandauBG_8", "LandauBG_P9", "LandauTt_Gaus_P10", "LandauBG_Gaus_P11" };

    return hTitle[id];
}

