#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

// define the fitting function
void ThirdJetEt2() {

 TFile *file  = TFile::Open("ttj_2Jskim.root");
 TFile *file1 = TFile::Open("wjets_2Jskim.root");
 TFile *file2 = TFile::Open("qcd_2Jskim.root");
 TString hfolder = "tt_test";

 TString name1 = "thirdJetEt";
 //TString name1 = "thirdCalEt";

 TString plot1 = "3rdJetEt_all.gif";

 thirdJetEt   = (TH1F *) file->Get("Jets/"+name1); 
 thirdJetEt1  = (TH1F *) file1->Get("Jets/"+name1); 
 thirdJetEt2  = (TH1F *) file2->Get("Jets/"+name1); 

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 gStyle->SetOptStat("nimou");
 gStyle->SetOptStat(kTRUE);
 //gStyle->SetOptFit(0111);  
 gStyle->SetStatY(0.90); 
 gStyle->SetStatTextColor(1);
 
 TCanvas *c1 = new TCanvas("c1","");
 c1->cd();
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->SetLogy();

 thirdJetEt->SetAxisRange(-0.5,150.5,"X");
 thirdJetEt->SetAxisRange(1,500,"Y");
 thirdJetEt->SetLineWidth(2);
 thirdJetEt->DrawCopy();
 c1->Update();

 //TPaveStats *st1 = (TPaveStats*) c1->GetPrimitive("stats");
 //st1->SetY1NDC(0.64); //new x start position
 //st1->SetY2NDC(0.8);  //new x end position
 
 gStyle->SetStatY(0.72); 
 gStyle->SetStatTextColor(2);
 thirdJetEt1->SetLineWidth(1);
 thirdJetEt1->SetLineColor(2);
 thirdJetEt1->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatY(0.54); 
 gStyle->SetStatTextColor(4);
 thirdJetEt2->SetLineColor(4);
 thirdJetEt2->DrawCopy("sames");

 c1->Update();
 c1->Print(plot1);

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}

