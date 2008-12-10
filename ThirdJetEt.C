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

void ThirdJetEt() {

 TFile *file  = TFile::Open("ttj_100pb_Et20.root");
 TFile *file1 = TFile::Open("wjets_10pb_Et20.root");
 TFile *file2 = TFile::Open("qcd_10pb_Et20.root");
 TString hfolder = "wjets_test";
 //TString hfolder = "tt_test";

 TString name9 = "thirdJetEt";
 //TString name9 = "thirdCalEt";
 TString name10 = "m3_j3";
 TString name3 = "MET_dPhi";

 TString plot1 = "3rdJetEt.gif";
 TString plot2 = "3rdJetEt_all.gif";
 TString plot3 = "m3_j3.gif";
 TString plot4 = "met_dphi.gif";

 thirdJetEt   = (TH1F *) file->Get("Jets/"+name9); 
 thirdJetEt1  = (TH1F *) file1->Get("Jets/"+name9); 
 thirdJetEt2  = (TH1F *) file2->Get("Jets/"+name9); 

 m3_j3   = (TH1F *) file->Get("Jets/"+name10); 
 m3_j3a  = (TH1F *) file1->Get("Jets/"+name10); 
 m3_j3b  = (TH1F *) file2->Get("Jets/"+name10); 

 met_df   = (TH1F *) file->Get("METs/"+name3); 
 met_dfa  = (TH1F *) file1->Get("METs/"+name3); 
 met_dfb  = (TH1F *) file2->Get("METs/"+name3); 

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);
 
 TCanvas *c1 = new TCanvas("c1","");
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->SetLogy();

 float nSoup[100] ={0.0} ;
 float nBg[100] ={0.0} ;
 float et[100] = {0.0};
 for (int i=0; i<100; i++) {
    float nWj = thirdJetEt1->GetBinContent(i+1);
    float nTt = thirdJetEt->GetBinContent(i+1);
    nBg[i] = nWj ;
    nSoup[i] = nWj+(nTt/10.) ;
    et[i]= i;
 }	 
 
 //thirdEt = new TMultiGraph();

 WjetEt = new TGraph(100, et, nBg);
 WjetEt->SetMaximum(2000.);
 WjetEt->SetMinimum(1.);
 WjetEt->SetMarkerSize(0.5);
 WjetEt->SetMarkerColor(4);
 WjetEt->SetFillColor(7);
 WjetEt->SetMarkerStyle(21);

 TF1 *func0 = new TF1("fit0",fitf,16,60,4);
 func0->SetParLimits(0,0.1,10000.);
 func0->SetParLimits(1,0.1,10000.);
 func0->SetParLimits(2,0.1,1000000.);
 func0->SetParLimits(3,0.1,10.);
 func0->SetLineColor(4);
 gStyle->SetOptStat(kTRUE);
 gStyle->SetStatX(0.90);
 gStyle->SetOptFit(0111);  
 WjetEt->Fit(fit0,"R","",16,30);
 WjetEt->Draw("APF");
 
 SoupEt = new TGraph(100, et, nSoup);
 SoupEt->SetMarkerSize(0.5);
 SoupEt->SetMarkerColor(2);
 SoupEt->SetMarkerStyle(20);
 
 TF1 *func0a = new TF1("fit0a",fitf,16,60,4);
 func0a->SetParLimits(0,0.1,10000.);
 func0a->SetParLimits(1,0.1,10000.);
 func0a->SetParLimits(2,0.1,1000000.);
 func0a->SetParLimits(3,0.1,10.);
 func0a->SetLineColor(2);
 SoupEt->Fit(fit0a,"R","",16,30);
 // doesn't work for setting stats box !!!
 //TPave *s = (TPave*) SoupEt->GetListOfFunctions()->FindObject("stats");
 //s->SetX1NDC(0.8);
 SoupEt->Draw("SAMEP");

 c1->Update();
 c1->Print(plot1);

 TCanvas *c2 = new TCanvas("c2","");
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->SetLogy();

 thirdJetEt->SetAxisRange(-0.5,150.5,"X");
 thirdJetEt->SetAxisRange(1,20000,"Y");
 thirdJetEt->SetLineWidth(2);
 thirdJetEt->DrawCopy();
 thirdJetEt1->SetLineWidth(1);
 thirdJetEt1->SetLineColor(2);
 thirdJetEt1->Scale(10);
 thirdJetEt1->DrawCopy("same");
 thirdJetEt2->SetLineColor(4);
 thirdJetEt2->Scale(10);
 thirdJetEt2->DrawCopy("same");

 c2->Update();
 c2->Print(plot2);

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

 TCanvas *c4 = new TCanvas("c4","");
 c4->SetFillColor(10);
 c4->SetFillColor(10);
 c4->Divide(2,2);

 c4->cd(1);
 met_df->SetAxisRange(1,150,"X");
 met_df->DrawCopy();

 c4->cd(2);
 met_dfa->SetAxisRange(1,150,"X");
 met_dfa->DrawCopy();

 c4->cd(3);
 met_dfb->SetAxisRange(1,150,"X");
 met_dfb->DrawCopy();

 c4->Update();
 c4->Print(plot4);


 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}

