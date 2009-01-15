#include <TFile.h>
#include <iostream>
#include <memory>
//#include <stdio.h>
//#include <fstream>

void NJetsPtCut() {

 TFile *file  = TFile::Open("ttj_100pb_Et20b.root");
 TFile *file1 = TFile::Open("wjets_10pb_Et20b.root");
 TFile *file2 = TFile::Open("qcd_10pb_Et20.root");
 TString hfolder = "TtBackground";

 TString plot1 = "njets.gif";

 h0jet20 = (TH1F *) file->Get("Jets/hNjet20");
 h0jet30 = (TH1F *) file->Get("Jets/hNjet30");
 h0jet40 = (TH1F *) file->Get("Jets/hNjet40");

 h1jet20 = (TH1F *) file1->Get("Jets/hNjet20");
 h1jet30 = (TH1F *) file1->Get("Jets/hNjet30");
 h1jet40 = (TH1F *) file1->Get("Jets/hNjet40");

 h2jet20 = (TH1F *) file2->Get("Jets/hNjet20");
 h2jet30 = (TH1F *) file2->Get("Jets/hNjet30");
 h2jet40 = (TH1F *) file2->Get("Jets/hNjet40");

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 //gStyle->SetOptStat("nimou");
 TCanvas *c1 = new TCanvas("c1","",200,8,900,700);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);
 c1->cd(1);
 c1_1->SetLogy();

 h0jet20->SetTitle(" NJets w/ pt > 20 GeV ");
 h0jet20->SetAxisRange(-0.5,16.5,"X");
 h0jet20->SetAxisRange(1,10000000,"Y");
 h0jet20->SetLineWidth(2);
 h0jet20->DrawCopy();

 h1jet20->SetLineColor(kRed);
 h1jet20->Scale(10.);
 h1jet20->DrawCopy("same");
 h2jet20->SetLineColor(kBlue);
 h2jet20->Scale(10.);
 h2jet20->DrawCopy("same");

 c1->cd(2);
 c1_2->SetLogy();

 h0jet30->SetTitle(" NJets w/ pt > 30 GeV ");
 h0jet30->SetAxisRange(-0.5,16.5,"X");
 h0jet30->SetAxisRange(1,10000000,"Y");
 h0jet30->SetLineWidth(2);
 h0jet30->DrawCopy();

 h1jet30->SetLineColor(kRed);
 h1jet30->Scale(10.);
 h1jet30->DrawCopy("same");
 h2jet30->SetLineColor(kBlue);
 h2jet30->Scale(10.);
 h2jet30->DrawCopy("same");

 c1->cd(3);
 c1_3->SetLogy();

 h0jet40->SetTitle(" NJets w/ pt > 40 GeV ");
 h0jet40->SetAxisRange(-0.5,16.5,"X");
 h0jet40->SetAxisRange(1,10000000,"Y");
 h0jet40->SetLineWidth(2);
 h0jet40->DrawCopy();

 h1jet40->SetLineColor(kRed);
 h1jet40->Scale(10.);
 h1jet40->DrawCopy("same");
 h2jet40->SetLineColor(kBlue);
 h2jet40->Scale(10.);
 h2jet40->DrawCopy("same");

 c1->cd(4);
 c1_4->SetLogy();
 
 TH1F *NJetbg = new TH1F("NJetbg","", 65, -0.5, 64.5);
 NJetbg->Add( h1jet20, h2jet20, 1., 1.);
 TH1F *NJetall = new TH1F("NJetall"," ", 65, -0.5, 64.5);
 NJetall->Add( h0jet20, NJetbg, 1., 1.);

 NJetall->SetTitle(" TTJets+WJets+QCD w/ pt > 20 GeV ");
 NJetall->SetAxisRange(-0.5,16.5,"X");
 NJetall->SetAxisRange( 1,10000000,"Y");
 NJetall->DrawCopy();

 c1->Print(plot1);
 c1->Update();

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}
