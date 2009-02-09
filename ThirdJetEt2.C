#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

// define the fitting function
void ThirdJetEt2() {

 TFile *file  = TFile::Open("ttj_fall08_4j25_3j30.root");
 TFile *file1 = TFile::Open("wjets_fall08_4j25_3j30.root");
 TFile *file2 = TFile::Open("wjets_fall08_4j25_3j30.root");
 TString hfolder = "tt_fall08";

 TString name1 = "Jet1Et";
 TString name2 = "Jet2Et";
 TString name3 = "Jet3Et";
 TString name4 = "Jet4Et";
 //TString name1 = "thirdJetEt";
 //TString name1 = "thirdCalEt";

 int qScale = 2.01 ;
 int rbin = 2;

 TString plot1 = "3rdJetEt_all.gif";

 Jet1Et   = (TH1F *) file->Get("Jets/"+name1); 
 Jet1Et1  = (TH1F *) file1->Get("Jets/"+name1); 
 Jet1Et2  = (TH1F *) file2->Get("Jets/"+name1); 

 Jet2Et   = (TH1F *) file->Get("Jets/"+name2); 
 Jet2Et1  = (TH1F *) file1->Get("Jets/"+name2); 
 Jet2Et2  = (TH1F *) file2->Get("Jets/"+name2); 

 Jet3Et   = (TH1F *) file->Get("Jets/"+name3); 
 Jet3Et1  = (TH1F *) file1->Get("Jets/"+name3); 
 Jet3Et2  = (TH1F *) file2->Get("Jets/"+name3); 

 Jet4Et   = (TH1F *) file->Get("Jets/"+name4); 
 Jet4Et1  = (TH1F *) file1->Get("Jets/"+name4); 
 Jet4Et2  = (TH1F *) file2->Get("Jets/"+name4); 

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 gStyle->SetOptStat("nimou");
 gStyle->SetOptStat(kTRUE);
 gStyle->SetStatY(0.92); 
 //gStyle->SetOptFit(0111);  
 TCanvas *c1 = new TCanvas("c1","",900,720);
 c1->cd();
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 c1->cd(1);
 c1_1->SetLogy();
 c1_1->SetGridx();
 gStyle->SetStatX(0.90); 
 gStyle->SetStatTextColor(1);
 
 Jet1Et->Rebin(rbin);
 Jet1Et->SetAxisRange(9.5,150,"X");
 Jet1Et->SetAxisRange(1,700,"Y");
 Jet1Et->SetLineWidth(2);
 Jet1Et->DrawCopy();
 c1->Update();

 //TPaveStats *st1 = (TPaveStats*) c1->GetPrimitive("stats");
 //st1->SetY1NDC(0.64); //new x start position
 //st1->SetY2NDC(0.8);  //new x end position
 
 gStyle->SetStatX(0.70); 
 gStyle->SetStatTextColor(2);
 Jet1Et1->SetLineWidth(1);
 Jet1Et1->SetLineColor(2);
 Jet1Et1->Rebin(rbin);
 Jet1Et1->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.50); 
 gStyle->SetStatTextColor(4);
 Jet1Et2->SetLineColor(4);
 Jet1Et2->Scale(qScale);
 Jet1Et2->Rebin(rbin);
 Jet1Et2->DrawCopy("sames");

 c1->Update();

 c1->cd(2);
 c1_2->SetLogy();
 c1_2->SetGridx();
 gStyle->SetStatX(0.90); 
 gStyle->SetStatTextColor(1);
 
 Jet2Et->Rebin(rbin);
 Jet2Et->SetAxisRange(9.5,150,"X");
 Jet2Et->SetAxisRange(1,700,"Y");
 Jet2Et->SetLineWidth(2);
 Jet2Et->DrawCopy();
 c1->Update();

 //TPaveStats *st1 = (TPaveStats*) c1->GetPrimitive("stats");
 //st1->SetY1NDC(0.64); //new x start position
 //st1->SetY2NDC(0.8);  //new x end position
 
 gStyle->SetStatX(0.70); 
 gStyle->SetStatTextColor(2);
 Jet2Et1->SetLineWidth(1);
 Jet2Et1->SetLineColor(2);
 Jet2Et1->Rebin(rbin);
 Jet2Et1->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.50); 
 gStyle->SetStatTextColor(4);
 Jet2Et2->SetLineColor(4);
 Jet2Et2->Scale(qScale);
 Jet2Et2->Rebin(rbin);
 Jet2Et2->DrawCopy("sames");

 c1->Update();

 c1->cd(3);
 c1_3->SetLogy();
 c1_3->SetGridx();
 gStyle->SetStatX(0.90); 
 gStyle->SetStatTextColor(1);
 
 Jet3Et->Rebin(rbin);
 Jet3Et->SetAxisRange(9.5,150,"X");
 Jet3Et->SetAxisRange(1,700,"Y");
 Jet3Et->SetLineWidth(2);
 Jet3Et->DrawCopy();
 c1->Update();

 //TPaveStats *st1 = (TPaveStats*) c1->GetPrimitive("stats");
 //st1->SetY1NDC(0.64); //new x start position
 //st1->SetY2NDC(0.8);  //new x end position
 
 gStyle->SetStatX(0.70); 
 gStyle->SetStatTextColor(2);
 Jet3Et1->SetLineWidth(1);
 Jet3Et1->SetLineColor(2);
 Jet3Et1->Rebin(rbin);
 Jet3Et1->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.50); 
 gStyle->SetStatTextColor(4);
 Jet3Et2->SetLineColor(4);
 Jet3Et2->Rebin(rbin);
 Jet3Et2->Scale(qScale);
 Jet3Et2->DrawCopy("sames");

 c1->Update();

 c1->cd(4);
 c1_4->SetLogy();
 c1_4->SetGridx();
 gStyle->SetStatX(0.90); 
 gStyle->SetStatTextColor(1);
 
 Jet4Et->Rebin(rbin);
 Jet4Et->SetAxisRange(9.5,150,"X");
 Jet4Et->SetAxisRange(1,1000,"Y");
 Jet4Et->SetLineWidth(2);
 Jet4Et->DrawCopy();
 c1->Update();

 //TPaveStats *st1 = (TPaveStats*) c1->GetPrimitive("stats");
 //st1->SetY1NDC(0.64); //new x start position
 //st1->SetY2NDC(0.8);  //new x end position
 
 gStyle->SetStatX(0.70); 
 gStyle->SetStatTextColor(2);
 Jet4Et1->SetLineWidth(1);
 Jet4Et1->SetLineColor(2);
 Jet4Et1->Rebin(rbin);
 Jet4Et1->DrawCopy("sames");
 c1->Update();

 gStyle->SetStatX(0.50); 
 gStyle->SetStatTextColor(4);
 Jet4Et2->SetLineColor(4);
 Jet4Et2->Rebin(rbin);
 Jet4Et2->Scale(qScale);
 Jet4Et2->DrawCopy("sames");

 c1->Update();
 c1->Print(plot1);

 gSystem->cd("../");
 file2->Close();
 file1->Close();
 file->Close();

}

