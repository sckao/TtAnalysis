#include <stdio.h>
#include <stdlib.h>
void GetCoefficient( TString idx, const Int_t nM ) {

    int id =  atoi( idx );

    TString plot = "p"+idx+"Scan.gif" ;
    TString hfolder = "AllTags0";
    TString title_plot = GiveTitle( id );

    FILE *pfile = fopen(hfolder+"/paraf.log","r");
    FILE *efile = fopen(hfolder+"/perrf.log","r");
    FILE *ffile = fopen(hfolder+"/extrf.log","a");

    Double_t pa[nM] ;
    Double_t er[nM] ;
    Double_t ex[nM] ;
    Double_t pr[nM] ;
    Double_t inputM[nM] ;

    double inputM5[5]   = { 161.2, 166.2, 171.2, 176.2, 181.2 };
    double inputM9[9]   = { 161.2, 163.7, 166.2, 168.7, 171.2, 173.7, 176.2, 178.7, 181.2 };
    double inputM10[10] = { 161.2, 163.7, 166.2, 168.7, 171.2, 173.7, 176.2, 178.7, 181.2, 183.7 };

    float f1;
    float e1;
    for ( int i=0 ; i < nM;  i++) {
        for ( int j=0 ; j< 12; j++) {

            fscanf(pfile, "%f" , &f1 );
            fscanf(efile, "%f" , &e1 );
            if ( j == id )  {
               pa[i] = f1;
               er[i] = e1;
               ex[i] = 0.;
            }            
        }
    }

    if ( nM == 5 ) FillPlotArray( inputM5, inputM, nM);
    if ( nM == 9 ) FillPlotArray( inputM9, inputM, nM);
    if ( nM == 10 ) FillPlotArray( inputM10, inputM, nM);

    gSystem->mkdir(hfolder);
    gSystem->cd(hfolder);

    TCanvas *c1 = new TCanvas("c1","", 700, 600);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    gPad->SetGridx();
    gPad->SetGridy();
    c1->cd();

    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(111);

    hMass = new TGraphErrors( nM, inputM, pa, ex, er );   
    hMass->SetMarkerColor(4);
    hMass->SetMarkerStyle(20);
    hMass->SetLineWidth(2);
    hMass->SetLineColor(4);
    hMass->SetTitle(title_plot);
    hMass->Draw("AP");

    double a0 = hMass->GetMean();
    double a1 = 1.;
    int fstatus = 999;
    //hMass->LeastSquareLinearFit(nM, a0, a1, fstatus, 150., 185);
    //cout<<" a0= "<<a0 <<"   a1= "<<a1<<"  fstatus:"<<fstatus<<endl;
    TF1 *func = new TF1("fitf",fitf,150,190, 2);
    hMass->Fit("fitf", "R0", "", 150., 185);
    a0 = func->GetParameter(0);
    a1 = func->GetParameter(1);
    cout<<" a0= "<<a0 <<"   a1= "<<a1 <<endl;
    fprintf(ffile," %d  %.3f  %.3f \n", id, a0, a1 );


    TF1 *fith = new TF1("fith",fitf,150,185,2);
    fith->FixParameter(0, a0);
    fith->FixParameter(1, a1);
    fith->SetLineStyle(2);
    //fith->SetLineColor(7);
    fith->SetLineColor(2);
    fith->Draw("same");

    c1->Update();

    c1->Print(plot);

    gSystem->cd("../");

    fclose(pfile);
    fclose(efile);
    fclose(ffile);
}

Double_t fitf(Double_t *x, Double_t *par) {

         Double_t fitval =  par[0] + par[1]*x[0] ;
         return fitval;
}

void FillPlotArray(Double_t* a1, Double_t* a2, const Int_t nM) {

    for ( int k=0; k<nM; k++ ) {
        a2[k] = a1[k] ;
    }
}

TString GiveTitle( int id ) {

    TString hTitle[12] = { "Gaus_P0", "Gaus_P1", "Gaus_P2", "LogNorm_P3", "LogNorm_P4", "LogNorm_P5",
                           "LandauTt_P6", "LandauTt_P7", "LandauBG_8", "LandauBG_P9", "LandauTt_Gaus", "LandauBG_Gaus" };

    return hTitle[id];
}
