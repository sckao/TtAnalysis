void TtPt() {

 TFile *file  = TFile::Open("ttj_fall08_3j30.root");
 TFile *file1 = TFile::Open("wjets_fall08_3j30.root");
 TFile *file2 = TFile::Open("qcd_fall08_3j30.root");
 TString hfolder = "tt_fall08";

 TString plot1 = "SysPt.gif";

 tMCPt    = (TH2F *) file->Get("Tops/PtMCtt"); 
 tRecPt1  = (TH2F *) file->Get("Tops/PtRecott1"); 
 tRecPt2  = (TH2F *) file->Get("Tops/PtRecott2"); 
 tRecPt3  = (TH2F *) file->Get("Tops/PtRecott3"); 
 tRecPt4  = (TH2F *) file->Get("Tops/PtRecott4"); 
 
 wRecPt1  = (TH2F *) file1->Get("Tops/PtRecott1"); 
 wRecPt2  = (TH2F *) file1->Get("Tops/PtRecott2"); 
 wRecPt3  = (TH2F *) file1->Get("Tops/PtRecott3"); 
 wRecPt4  = (TH2F *) file1->Get("Tops/PtRecott4"); 

 qRecPt1  = (TH2F *) file2->Get("Tops/PtRecott1"); 
 qRecPt2  = (TH2F *) file2->Get("Tops/PtRecott2"); 
 qRecPt3  = (TH2F *) file2->Get("Tops/PtRecott3"); 
 qRecPt4  = (TH2F *) file2->Get("Tops/PtRecott4"); 
 

 // re-bin value for top mass plot
 int rbin = 4;
 int nbin = 200./rbin ;

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 TCanvas *c1 = new TCanvas("c1","", 1000, 800);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);

 c1->cd(1);

 tMCPt->Rebin(rbin);
 tRecPt1->Rebin(rbin);
 tRecPt2->Rebin(rbin);
 tRecPt3->Rebin(rbin);
 tRecPt4->Rebin(rbin);

 TH1F * tRecPt = new TH1F("tMCPt0","", nbin, 0., 400.);
 TH1F * tRecPt = new TH1F("tRecPt","", nbin, 0., 400.);

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.77); 
 gStyle->SetStatTextColor(1);

 tMCPt->SetLineWidth(2);
 tMCPt->SetLineColor(1);
 tMCPt->SetTitle("Pt of Tt system in Tt+Jets Events");
 tMCPt->Draw();
 tMCPt0->Add(tMCPt);
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(2);
 tRecPt1->SetLineColor(2);
 tRecPt1->Draw("SAMES");
 tRecPt->Add(tRecPt1);
 c1->Update();

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(4);
 tRecPt2->SetLineColor(4);
 tRecPt2->Draw("SAMES");
 tRecPt->Add(tRecPt2);
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(8);
 tRecPt3->SetLineColor(8);
 tRecPt3->Draw("SAMES");
 tRecPt->Add(tRecPt3);
 c1->Update();

 gStyle->SetStatY(0.51); 
 gStyle->SetStatTextColor(5);
 tRecPt4->SetLineColor(5);
 tRecPt4->Draw("SAMES");
 tRecPt->Add(tRecPt4);
 c1->Update();

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(1);

 c1->cd(2);

 wRecPt1->Rebin(rbin);
 wRecPt2->Rebin(rbin);
 wRecPt3->Rebin(rbin);
 wRecPt4->Rebin(rbin);

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.77); 
 gStyle->SetStatTextColor(2);
 wRecPt1->SetLineColor(2);
 wRecPt1->SetTitle("Pt of Tt system in W+Jets Events ");
 wRecPt1->Draw();
 c1->Update();

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(4);
 wRecPt2->SetLineColor(4);
 wRecPt2->Draw("SAMES");
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(8);
 wRecPt3->SetLineColor(8);
 wRecPt3->Draw("SAMES");
 c1->Update();

 gStyle->SetStatY(0.51); 
 gStyle->SetStatTextColor(5);
 wRecPt4->SetLineColor(5);
 wRecPt4->Draw("SAMES");
 c1->Update();

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(1);

 
 c1->cd(3);

 qRecPt1->Rebin(rbin);
 qRecPt2->Rebin(rbin);
 qRecPt3->Rebin(rbin);
 qRecPt4->Rebin(rbin);

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.77); 
 gStyle->SetStatTextColor(2);
 qRecPt1->SetLineColor(2);
 qRecPt1->SetTitle("Pt of Tt system in QCD Events ");
 qRecPt1->Draw();
 c1->Update();

 gStyle->SetStatY(0.99); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(4);
 qRecPt2->SetLineColor(4);
 qRecPt2->Draw("SAMES");
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(8);
 qRecPt3->SetLineColor(8);
 qRecPt3->Draw("SAMES");
 c1->Update();

 gStyle->SetStatY(0.51); 
 gStyle->SetStatTextColor(5);
 qRecPt4->SetLineColor(5);
 qRecPt4->Draw("SAMES");
 c1->Update();

 c1->cd(4);
 

 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(1);

 double nMC = tMCPt0->Integral(1,nbin);
 cout<<" nMC = "<<nMC << endl; 
 tMCPt0->SetLineWidth(1);
 tMCPt0->SetLineColor(1);
 tMCPt0->DrawCopy();
 c1->Update();

 gStyle->SetStatY(0.75); 
 gStyle->SetStatTextColor(2);
 double nRec = tRecPt->Integral(1,nbin);
 double norV = nMC/nRec ;
 cout<<" nRec = "<<nRec << endl; 
 cout<<" norv = " << norV <<endl;
 tRecPt->Scale( norV ) ;
 tRecPt->SetLineColor(2);
 tRecPt->DrawCopy("SAMES");
 c1->Update();

 gStyle->SetStatY(0.95); 
 gStyle->SetStatX(0.99); 
 gStyle->SetStatTextColor(1);
 
 c1->cd();

 c1->Update();
 c1->Print(plot1);

 // Reset gStyle
 gStyle->SetStatY(0.95); 
 gStyle->SetStatTextColor(1);

 gSystem->cd("../");

 file2->Close();
 file1->Close();
 file->Close();

}
