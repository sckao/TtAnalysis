#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

// define the fitting function
Double_t fitf(Double_t *x, Double_t *par) {

	 Double_t A0 = par[0] + 0.01 ;
         Double_t A1 = (x[0])/ A0 ;
	 Double_t A2 = pow( (x[0]-par[3]), par[1]) ;
         Double_t fitval =  par[2] / ( exp(A1) * A2 );
         return fitval;
}

Double_t fitG(Double_t *x, Double_t *par) {

	Double_t chi = (x[0]-par[1])/par[2] ;
	Double_t A3  = -0.5*chi*chi ;
	Double_t fitgau = par[0]*exp( A3 ) ;
	return fitgau ;
}

Double_t fitCb(Double_t *x, Double_t *par) {

        return fitf(x,par) + fitG(x, &par[3] );
}

void M3() {

 TFile *file  = TFile::Open("ttj_100pb_Et20a.root");
 TFile *file1 = TFile::Open("wjets_10pb_Et20a.root");
 TFile *file2 = TFile::Open("qcd_10pb_Et20a.root");
 TString hfolder = "wjets_test";
 //TString hfolder = "tt_test";

 TString name10 = "m3_j3";

 TString plot3 = "m3_j3.gif";

 m3_j3   = (TH1F *) file->Get("Jets/"+name10); 
 m3_j3a  = (TH1F *) file1->Get("Jets/"+name10); 
 m3_j3b  = (TH1F *) file2->Get("Jets/"+name10); 

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 TCanvas *c3 = new TCanvas("c3","");
 c3->SetFillColor(10);
 c3->SetFillColor(10);
 c3->Divide(2,2);

 c3->cd(1);
 m3_j3->SetAxisRange(1,150,"Y");
 m3_j3->DrawCopy();

 c3->cd(2);
 m3_j3a->SetAxisRange(1,150,"Y");
 m3_j3a->DrawCopy();

 c3->cd(3);
 m3_j3b->SetAxisRange(1,150,"Y");
 m3_j3b->DrawCopy();

 c3->Update();
 c3->Print(plot3);

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}


