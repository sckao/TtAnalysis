void CutEff() {

 TFile *file  = TFile::Open("ttj_fall08.root");
 TFile *file1 = TFile::Open("wjets_fall08.root");
 TFile *file2 = TFile::Open("qcd_fall08.root");
 

 TString plot1 = "ObjSelectionEff.gif";
 TString plot2 = "EvtSolutionEff.gif";

 tMuEff  = (TH2F *) file->Get("Muons/allMu_isoMu"); 
 wMuEff  = (TH2F *) file1->Get("Muons/allMu_isoMu"); 
 qMuEff  = (TH2F *) file2->Get("Muons/allMu_isoMu"); 

 tJetEff  = (TH2F *) file->Get("Jets/allJ_selJ"); 
 wJetEff  = (TH2F *) file1->Get("Jets/allJ_selJ"); 
 qJetEff  = (TH2F *) file2->Get("Jets/allJ_selJ"); 

 tMuJEff  = (TH2F *) file->Get("Ele/isoEleCut"); 
 wMuJEff  = (TH2F *) file1->Get("Ele/isoEleCut"); 
 qMuJEff  = (TH2F *) file2->Get("Ele/isoEleCut"); 

 tLepWEff  = (TH1F *) file->Get("Tops/hRecolepW"); 
 wLepWEff  = (TH1F *) file1->Get("Tops/hRecolepW"); 
 qLepWEff  = (TH1F *) file2->Get("Tops/hRecolepW"); 

 thadWEff  = (TH1F *) file->Get("Tops/hRecohadW"); 
 whadWEff  = (TH1F *) file1->Get("Tops/hRecohadW"); 
 qhadWEff  = (TH1F *) file2->Get("Tops/hRecohadW"); 

 tTtEff    = (TH2F *) file->Get("Tops/hRecott0"); 
 wTtEff    = (TH2F *) file1->Get("Tops/hRecott0"); 
 qTtEff    = (TH2F *) file2->Get("Tops/hRecott0"); 

 // 1. Objects Selection Efficiency
 TString hfolder = "tt_fall08";
 int rbin = 4;

 gSystem->mkdir(hfolder);
 gSystem->cd(hfolder);

 TCanvas *c1 = new TCanvas("c1","", 900, 750);
 c1->SetFillColor(10);
 c1->SetFillColor(10);
 c1->Divide(2,2);
 gStyle->SetOptStat(1001101);
 c1->cd(1);
 c1_1->SetLogy();

 cout<<" ========================================= "<<endl;
 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 double nTt = tMuEff->Integral();
 tMuEff->ProjectionX("tMuEff_px",2,2,"");
 tMuEff_px->SetAxisRange(1,30000,"Y");
 tMuEff_px->DrawCopy();
 double nMuTt = tMuEff_px->Integral();
 double sMuTt_eff = nMuTt/nTt ;
 cout<<"  --- Isolated Muon Selection --- "<<endl;
 cout<<"Tt n_Evt: "<<nTt<<"         n_1Mu: "<<nMuTt<<"  muEff="<< sMuTt_eff <<endl;
 c1->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 double nWj = wMuEff->Integral();
 wMuEff->ProjectionX("wMuEff_px",2,2,"");
 wMuEff_px->SetLineColor(2);
 wMuEff_px->DrawCopy("sames");
 double nMuWj = wMuEff_px->Integral();
 double sMuWj_eff = nMuWj/nWj ;
 cout<<"Wj n_Evt: "<<nWj<<"        n_1Mu: "<<nMuWj<<" muEff="<< sMuWj_eff <<endl;
 c1->Update(); 

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 double nQm = qMuEff->Integral();
 qMuEff->ProjectionX("qMuEff_px",2,2,"");
 qMuEff_px->SetLineColor(4);
 qMuEff_px->DrawCopy("sames");
 double nMuQm = qMuEff_px->Integral();
 double sMuQm_eff = nMuQm/nQm ;
 cout<<"Qm n_Evt: "<<nQm<<"  n_1Mu: "<<nMuQm<<"  muEff="<< sMuQm_eff <<endl;
 c1->Update(); 

 cout<<" ----------------------------------------- "<<endl;
 c1->cd(2);
 c1_2->SetLogy();
 c1_2->SetGridx();

 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 tJetEff->ProjectionY("tJetEff_py",1,21,"");
 tJetEff_py->SetAxisRange(1,1000000,"Y");
 tJetEff_py->DrawCopy();
 double nJetTt = tJetEff_py->Integral(5,21);
 double sJetTt_eff = nJetTt/nTt ;
 cout<<"  --- At least 4 Good Jet Selection --- "<<endl;
 cout<<"Tt n_Evt: "<<nTt<<"         n_Jet: "<<nJetTt<<"   jetEff="<< sJetTt_eff <<endl;
 c1->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 wJetEff->ProjectionY("wJetEff_py",1,21,"");
 wJetEff_py->SetLineColor(2);
 wJetEff_py->DrawCopy("sames");
 double nJetWj = wJetEff_py->Integral(5,21);
 double sJetWj_eff = nJetWj/nWj ;
 cout<<"Wj n_Evt: "<<nWj<<"        n_Jet: "<<nJetWj<<"   jetEff="<< sJetWj_eff <<endl;
 c1->Update(); 

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 qJetEff->ProjectionY("qJetEff_py",1,21,"");
 qJetEff_py->SetLineColor(4);
 qJetEff_py->DrawCopy("sames");
 double nJetQm = qJetEff_py->Integral(5,21);
 double sJetQm_eff = nJetQm/nQm ;
 cout<<"Qm n_Evt: "<<nQm<<"  n_Jet: "<<nJetQm<<" jetEff="<< sJetQm_eff <<endl;
 c1->Update(); 
 cout<<" ----------------------------------------- "<<endl;

 c1->cd(3);
 c1_3->SetLogy();
 c1_3->SetGridx();

 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 tMuJEff->ProjectionY("tMuJEff_py",1,21,"");
 tMuJEff_py->SetAxisRange(1,100000,"Y");
 tMuJEff_py->DrawCopy();
 double nMuJTt = tMuJEff_py->Integral(5,21);
 double sMuJTt_eff = nMuJTt/nTt ;
 cout<<"  --- IsoMuon + 4GoodJet Selection --- "<<endl;
 cout<<"Tt n_Evt: "<<nTt<<"         n_MuJ: "<<nMuJTt<<" combinedEff="<< sMuJTt_eff <<endl;
 c1->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 wMuJEff->ProjectionY("wMuJEff_py",1,21,"");
 wMuJEff_py->SetLineColor(2);
 wMuJEff_py->DrawCopy("sames");
 double nMuJWj = wMuJEff_py->Integral(5,21);
 double sMuJWj_eff = nMuJWj/nWj ;
 cout<<"Wj n_Evt: "<<nWj<<"        n_MuJ: "<<nMuJWj<<" combinedEff="<< sMuJWj_eff <<endl;
 c1->Update(); 

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 qMuJEff->ProjectionY("qMuJEff_py",1,21,"");
 qMuJEff_py->SetLineColor(4);
 qMuJEff_py->DrawCopy("sames");
 double nMuJQm = qMuJEff_py->Integral(5,21);
 double sMuJQm_eff = nMuJQm/nQm ;
 cout<<"Qm n_Evt: "<<nQm<<"  n_MuJ: "<<nMuJQm<<"  combinedEff="<< sMuJQm_eff <<endl;
 c1->Update(); 
 c1->Print(plot1);
 cout<<" ========================================= "<<endl;


 // 2.Event Solution Efficiency 
 TCanvas *c2 = new TCanvas("c2","", 900, 750);
 c2->SetFillColor(10);
 c2->SetFillColor(10);
 c2->Divide(2,2);
 gStyle->SetOptStat(1001101);
 c2->cd(1);

 // Leptonic W solution
 // using nMuTt,nTt, nMuWj,nWj, nMuQm, nQm
 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 cout<<"  --- Leptonic W Solution --- "<<endl;
 tLepWEff->Rebin(rbin);
 tLepWEff->DrawCopy();
 double tNlepW = tLepWEff->Integral();
 double tlepW_eff  = tNlepW/nMuTt ;
 double tlepW_eff0 = tNlepW/nTt ;
 cout<<"Tt n_lepW: "<<tNlepW<<"  lepW/nMuon= "<< tlepW_eff <<" lepW/nTt= "<< tlepW_eff0 <<endl;
 c2->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 wLepWEff->SetLineColor(2);
 wLepWEff->Rebin(rbin);
 wLepWEff->DrawCopy("sames");
 double wNlepW = wLepWEff->Integral();
 double wlepW_eff = wNlepW/nMuWj ;
 double wlepW_eff0 = wNlepW/nWj ;
 cout<<"Wj n_lepW: "<<wNlepW<<"  lepW/nMuon= "<< wlepW_eff <<" lepW/nTt= "<< wlepW_eff0 <<endl;
 c2->Update();

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 qLepWEff->SetLineColor(4);
 qLepWEff->Rebin(rbin);
 qLepWEff->DrawCopy("sames");
 double qNlepW     = qLepWEff->Integral();
 double qlepW_eff  = qNlepW/nMuQm ;
 double qlepW_eff0 = qNlepW/nQm ;
 cout<<"Qm n_lepW: "<<qNlepW<<"  lepW/nMuon= "<< qlepW_eff <<" lepW/nTt= "<< qlepW_eff0 <<endl;
 c2->Update();

 // Hadronic W solution 
 // using nJetTt,nWj
 c2->cd(2);

 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 cout<<"  --- Hadronic W Solution --- "<<endl;
 thadWEff->Rebin(rbin);
 thadWEff->DrawCopy();
 double tNhadW = thadWEff->Integral();
 double thadW_eff  = tNhadW/nJetTt ;
 double thadW_eff0 = tNhadW/nTt ;
 cout<<"Tt n_hadW: "<<tNhadW<<"   hadW/nJet= "<< thadW_eff <<"  hadW/nTt= "<< thadW_eff0 <<endl;
 c2->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 whadWEff->SetLineColor(2);
 whadWEff->Rebin(rbin);
 whadWEff->DrawCopy("sames");
 double wNhadW = whadWEff->Integral();
 double whadW_eff = wNhadW/nJetWj ;
 double whadW_eff0 = wNhadW/nWj ;
 cout<<"Wj n_hadW: "<<wNhadW<<"   hadW/nJet= "<< whadW_eff <<"  hadW/nTt= "<< whadW_eff0 <<endl;
 c2->Update();

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 qhadWEff->SetLineColor(4);
 thadWEff->Rebin(rbin);
 qhadWEff->DrawCopy("sames");
 double qNhadW     = qhadWEff->Integral();
 double qhadW_eff  = qNhadW/nJetQm ;
 double qhadW_eff0 = qNhadW/nQm ;
 cout<<"Qm n_hadW: "<<qNhadW<<"    hadW/nJet= "<< qhadW_eff <<"  hadW/nTt= "<< qhadW_eff0 <<endl;
 c2->Update();

 cout<<" ========================================= "<<endl;

 c2->cd(3);

 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 cout<<"  --- Tt Event Solution Selection --- "<<endl;
 cout<<" "<<endl;
 tTtEff->ProjectionX("tTtEff_px",-1,-1,"");
 tTtEff_px->SetLineColor(1);
 tTtEff_px->Rebin(rbin);
 tTtEff_px->DrawCopy();
 double nTtEvt = tTtEff_px->Integral();
 double sTt_eff  = nTtEvt/nTt ;
 double sTt_eff0 = nTtEvt/nMuJTt ;
 cout<<"Tt n_Evt: "<<nTtEvt <<"  Solution Eff:"<<sTt_eff0<<"  EvtEff= "<< sTt_eff <<endl;
 c2->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 wTtEff->ProjectionX("wTtEff_px",-1,-1,"");
 wTtEff_px->SetLineColor(2);
 wTtEff_px->Rebin(rbin);
 wTtEff_px->DrawCopy("sames");
 double nWjEvt = wTtEff_px->Integral();
 double sWj_eff  = nWjEvt/nWj ;
 double sWj_eff0 = nWjEvt/nMuJWj ;
 cout<<"Wj n_Evt: "<<nWjEvt <<"  Solution Eff:"<<sWj_eff0<<"    EvtEff= "<< sWj_eff <<endl;
 c2->Update();

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 qTtEff->ProjectionX("qTtEff_px",-1,-1,"");
 qTtEff_px->SetLineColor(4);
 qTtEff_px->Rebin(rbin);
 qTtEff_px->DrawCopy("sames");
 double nQmEvt = qTtEff_px->Integral();
 double sQm_eff  = nQmEvt/nQm ;
 double sQm_eff0 = nQmEvt/nMuJQm ;
 cout<<"Qm n_Evt: "<<nQmEvt <<"   Solution Eff:"<<sQm_eff0<<"  EvtEff= "<< sQm_eff <<endl;
 cout<<" ========================================= "<<endl;
 c2->Update();

 c2->cd(4);

 gStyle->SetStatY(0.99);
 gStyle->SetStatTextColor(1);
 tTtEff->ProjectionY("tTtEff_py",-1,-1,"");
 tTtEff_py->SetLineColor(1);
 tTtEff_py->Rebin(rbin);
 tTtEff_py->DrawCopy();
 c2->Update();

 gStyle->SetStatY(0.75);
 gStyle->SetStatTextColor(2);
 wTtEff->ProjectionY("wTtEff_py",-1,-1,"");
 wTtEff_py->SetLineColor(2);
 wTtEff_py->Rebin(rbin);
 wTtEff_py->DrawCopy("sames");
 c2->Update();

 gStyle->SetStatY(0.51);
 gStyle->SetStatTextColor(4);
 qTtEff->ProjectionY("qTtEff_py",-1,-1,"");
 qTtEff_py->SetLineColor(4);
 qTtEff_py->Rebin(rbin);
 qTtEff_py->DrawCopy("sames");

 c2->Update();
 c2->Print(plot2);

 gSystem->cd("../");
 gStyle->SetOptStat(1111);
 gStyle->SetStatY(0.98);
 gStyle->SetStatTextColor(1);

 file2->Close();
 file1->Close();
 file->Close();

}
