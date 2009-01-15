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

Double_t fitP(Double_t *x, Double_t *par) {

      Double_t fitV = par[0] + (par[1]*x[0]) + (par[2]*x[0]*x[0]) 
	            + (par[3]*x[0]*x[0]*x[0]) + (par[4]*x[0]*x[0]*x[0]*x[0]) ;
      return fitV;
}

Double_t fitCb(Double_t *x, Double_t *par) {

        return fitf(x,par) + fitG(x, &par[3] );
}

void ThirdJetEt() {

 TFile *file  = TFile::Open("ttj_100pb_Et20a.root");
 TFile *file1 = TFile::Open("wjets_10pb_Et20a.root");
 TFile *file2 = TFile::Open("qcd_10pb_Et20a.root");
 TString hfolder = "wjets_test";

 TString name1 = "thirdJetEt";
 //TString name1 = "thirdCalEt";

 TString plot1 = "3rdJetEt.gif";
 TString plot2 = "3rdJetEt_all.gif";
 TString plot3 = "Tt3rdJetEt.gif";
 TString plot4 = "wj3rdJetEt.gif";

 thirdJetEt   = (TH1F *) file->Get("Jets/"+name1); 
 thirdJetEt1  = (TH1F *) file1->Get("Jets/"+name1); 
 thirdJetEt2  = (TH1F *) file2->Get("Jets/"+name1); 

 thirdJetEt1->Scale(10.);
 thirdJetEt2->Scale(10.);

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 //gStyle->SetOptStat("nimou");
 gStyle->SetOptStat(kTRUE);
 gStyle->SetOptFit(0111);  
 //gStyle->SetStatX(0.90);
 
 
 TCanvas *c2 = new TCanvas("c2","");
 c2->cd();
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->SetLogy();

 thirdJetEt->SetAxisRange(-0.5,150.5,"X");
 thirdJetEt->SetAxisRange(1,20000,"Y");
 thirdJetEt->SetLineWidth(2);
 thirdJetEt->DrawCopy();
 thirdJetEt1->SetLineWidth(1);
 thirdJetEt1->SetLineColor(4);
 thirdJetEt1->DrawCopy("same");
 thirdJetEt2->SetLineColor(8);
 thirdJetEt2->DrawCopy("same");

 c2->Update();
 c2->Print(plot2);

 TCanvas *c1 = new TCanvas("c1","",800, 600);
 c1->cd();
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->SetLogy();

 thirdJetEt1->SetLineColor(4);
 thirdJetEt1->SetLineWidth(2);
 thirdJetEt1->SetAxisRange(-0.5,150.5,"X");
 thirdJetEt1->SetAxisRange(1,20000,"Y");
 thirdJetEt1->DrawCopy();

 TF1 *func0 = new TF1("fit0",fitf,20,30,4);
 func0->SetParLimits(0,0.1,10000.);
 func0->SetParLimits(1,0.1,10000.);
 func0->SetParLimits(2,0.1,1000000.);
 func0->SetParLimits(3,0.1,10.);
 func0->SetLineColor(4);
 gStyle->SetStatX(0.90);
 thirdJetEt1->Fit(fit0,"R","",20,60);
 thirdJetEt1->Fit(fit0,"R","",20,30);

 double p0 = func0->GetParameter(0);
 double p1 = func0->GetParameter(1);
 double p2 = func0->GetParameter(2);
 double p3 = func0->GetParameter(3);

 TF1 *func0a = new TF1("func0a",fitf,20,60,4);
 func0a->FixParameter(0, p0);
 func0a->FixParameter(1, p1);
 func0a->FixParameter(2, p2);
 func0a->FixParameter(3, p3);
 func0a->SetLineStyle(2);
 func0a->SetLineColor(4);
 func0a->Draw("same");

 TH1F* SoupEt = new TH1F("SoupEt","",500, 0., 500.);
 SoupEt->Add(thirdJetEt, 1.);
 SoupEt->Add(thirdJetEt1,1.);
 SoupEt->SetLineColor(2);
 SoupEt->DrawCopy("SAME");
 
 TF1 *func1 = new TF1("fit1",fitf,20,30,4);
 func1->SetParLimits(0,0.1,10000.);
 func1->SetParLimits(1,0.1,10000.);
 func1->SetParLimits(2,0.1,1000000.);
 func1->SetParLimits(3,0.1,10.);
 func1->SetLineColor(2);
 func1->SetLineWidth(1);
 gStyle->SetStatX(0.99);
 SoupEt->Fit(fit1,"R","same",20,30);

 //TPaveStats *ss = (TPaveStats*) SoupEt->GetListOfFunctions()->FindObject("stats");
 //ss->SetX1NDC(0.7);
 //SoupEt->FindObject("stats")->SetX1NDC(0.7);
 double t[4] ={0.0};
 double t[0] = func1->GetParameter(0);
 double t[1] = func1->GetParameter(1);
 double t[2] = func1->GetParameter(2);
 double t[3] = func1->GetParameter(3);

 TF1 *func1a = new TF1("func1a",fitf,30,60,4);
 func1a->FixParameter(0, t[0] );
 func1a->FixParameter(1, t[1] );
 func1a->FixParameter(2, t[2] );
 func1a->FixParameter(3, t[3] );
 func1a->SetLineStyle(2);
 func1a->SetLineColor(2);
 func1a->Draw("same");
 
 // doesn't work for setting stats box; fucking sucks!!!
 //TPave *s = (TPave*) SoupEt->GetListOfFunctions()->FindObject("stats");
 //s->SetX1NDC(0.8);

 c1->Update();
 c1->Print(plot1);

 TH1F* rtt = new TH1F("rtt","",500,0.,500.);
 double x1 =30.5;
 for (int i=30; i< 80; i++) {
     double Nfc = func1a->Eval( x1, 0 , 0);
     double NSoup = SoupEt->GetBinContent(i);
     double Ntt = NSoup - Nfc ;
     rtt->Fill(x1,Ntt);
     x1 = x1 + 1. ;
 }

 TCanvas *c3 = new TCanvas("c3","");
 c3->cd();
 c3->SetFillColor(10);
 c3->SetFillColor(10);
 c3->SetLogy();

 rtt->SetAxisRange(-0.5,150.5,"X");
 rtt->SetAxisRange(1,20000,"Y");
 TF1 *func2 = new TF1("fit2",fitP,30,80,5);
 rtt->Fit(fit2,"R","",30,80);

 c3->Update();
 c3->Print(plot3);

 TCanvas *c4 = new TCanvas("c4","",800, 600);
 c4->cd();
 c4->SetFillColor(10);
 c4->SetFillColor(10);
 c4->SetLogy();

 double x2 =0.5;
 double xx[150]={0.0};
 double yy[150]={0.0};
 for (int i=1; i< 150; i++) {
     double Nfc   = func2->Eval( x2, 0 , 0);
     double NSoup = SoupEt->GetBinContent(i);
     double Nwj = NSoup - Nfc ;
     int k = i-1 ;
     if (Nwj < 0. ) {
	Nwj = 0;
     }
     if ( x2 > 15. && x2 < 80. ) {
        xx[k] = x2;
        yy[k] = Nwj;
     } else {
        xx[k] = x2;
        yy[k] = NSoup;
     }	     
     x2 = x2 + 1. ;
 }

 rBg = new TGraph(150, xx, yy);
 rBg->SetMaximum(20000.0);
 rBg->SetMinimum(1.0);
 rBg->SetMarkerSize(0.5);
 rBg->SetMarkerColor(6);
 rBg->SetMarkerStyle(21);
 rBg->Draw("AP");

 TF1 *func3 = new TF1("fit3", fitf, 20, 30, 4);
 func3->SetParLimits(0,0.1,10000.);
 func3->SetParLimits(1,0.1,10000.);
 func3->SetParLimits(2,0.1,1000000.);
 func3->SetParLimits(3,0.1,10.);
 func3->SetLineColor(1);
 rBg->Fit(fit3,"R","",20,30);
 func0a->Draw("same");

 double g0 = func3->GetParameter(0);
 double g1 = func3->GetParameter(1);
 double g2 = func3->GetParameter(2);
 double g3 = func3->GetParameter(3);

 TF1 *func3a = new TF1("func3a",fitf,20,60,4);
 func3a->FixParameter(0, g0);
 func3a->FixParameter(1, g1);
 func3a->FixParameter(2, g2);
 func3a->FixParameter(3, g3);
 func3a->SetLineStyle(2);
 func3a->SetLineColor(2);
 func3a->Draw("same");

 //rBg->Draw("SAMEPL");

 c4->Update();
 c4->Print(plot4);

 //cout<<"  p0:"<<p0<<"  p1:"<<p1<<"  p2:"<<p2<<"  p3:"<<p3<<endl;

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}

