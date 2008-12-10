#include <vector>
#include <stdio.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

void TtMass() {
 TFile *file  = TFile::Open("ttj_100pb_Et20.root");
 TFile *file1 = TFile::Open("wjets_10pb_Et20.root");
 TFile *file2 = TFile::Open("qcd_10pb_Et20.root");
 TString hfolder = "tt_test";

 TString name1 = "hRecott";

 TString plot1 = "TtMass.gif";
 TString plot2 = "TtMassFit.gif";

 ttmass   = (TH1F *) file->Get("Tops/"+name1); 
 ttmass1  = (TH1F *) file1->Get("Tops/"+name1); 
 ttmass2  = (TH1F *) file2->Get("Tops/"+name1);

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);
 

 TCanvas *c1 = new TCanvas("c1","");
 c1->SetFillColor(10);
 c1->SetFillColor(10);

 THStack *ttall = new THStack("ttall", " Tt+Wjets+QCD ");

 ttmass1->SetFillColor(2);
 ttmass1->Rebin(2);
 ttmass1->Scale(10);
 ttall->Add(ttmass1);

 ttmass2->SetFillColor(4);
 ttmass2->Rebin(2);
 ttmass2->Scale(10);
 ttall->Add(ttmass2);

 //ttmass->Fit("gaus","R","",120,220);
 ttmass->SetFillColor(3);
 ttmass->Rebin(2);
 ttall->Add(ttmass);

 ttall->Draw();

 c1->Update();
 c1->Print(plot1);

 
 TCanvas *c2 = new TCanvas("c2","");
 c2->SetFillColor(10);
 c2->SetFillColor(10);

 TH1F *htt = ttall->GetHistogram() ;
 htt->SetLineColor(1);
 htt->Fit("gaus","R","",120,220);
 htt->DrawCopy();

 c2->Update();
 c2->Print(plot2);

 gSystem->cd("../");

 file2->Close();
 file1->Close();
 file->Close();
}
