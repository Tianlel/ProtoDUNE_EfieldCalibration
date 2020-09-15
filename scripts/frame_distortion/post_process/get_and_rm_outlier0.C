/* get outlier bins and remove them accordingly */

#include <fstream>

int YZ_deltaT_min = 4000;
int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm
double X0 = 358.6; // cm
double driftv = 0.156; //0.156461 cm/us 

// check for outlier boundary bins
bool check_out(int nbinsY, int nbinsZ, TH2F *h, int i, int j, int diff_min)
{
    float binval = h->GetBinContent(i+1, j+1);
    if (binval < 2000) return false;
    float up=0, down=0, right=0, left=0, tmp;
    if (i!=0) 
    {
        tmp = h->GetBinContent(i, j+1);
        if (tmp > 2000 && abs(tmp-binval)>=diff_min) up = 1;
    }
    if (i!=nbinsZ-1) 
    {
        tmp = h->GetBinContent(i+2, j+1);
        if (tmp > 2000 && abs(tmp-binval)>=diff_min) down = 1;
    }
    if (j!=0) 
    {
        tmp = h->GetBinContent(i+1, j);
        if (tmp > 2000 && abs(tmp-binval)>=diff_min) left = 1;
    }
    if (j!=nbinsY-1) 
    {
        tmp = h->GetBinContent(i+1, j+2);
        if (tmp > 2000 && abs(tmp-binval)>=diff_min) right = 1;
    }
   
    if (up+down+right+left > 0) return true;
    else return false;
}

float get_distortion(float deltaT, int neg)
{
    if (deltaT == 0) return 0;
    switch(neg)
    {
        case 0: 
            return 0.15-(X0-deltaT*driftv/2);
        case 1:
            return (X0-deltaT*driftv/2)-0.15;  
    }   
    return 0;
}

void get_and_rm_outlier0()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    string OUTPUT_TREE_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    TString save_name = TString::Format("%supdated_deltaT_med_TH2", OUTPUT_PATH.c_str());
    TString save_name_tree = TString::Format("%s", OUTPUT_TREE_PATH.c_str()); 

    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/ALLEVENTS_without_thermal_unordered_dT.root", "READ");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int Tbinsize = 2;
    int nbinsY = 30, nbinsZ = 36;
    int nentries_min = 1000;

    // selection range
    int deltaT_min = 4580;
    int deltaT_max = 4700;

    TH2F *h = new TH2F("dT_rmBoundary_pos","dT after removing bad boundary bins (beam left); Z (cm); Y(cm); dT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *h_neg = new TH2F("dT_rmBoundary_neg","dT after removing bad boundary bins (beam right); Z (cm); Y(cm); dT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    TH2F *deltaT_YZ_h2, *deltaT_YZ_h2_neg;
    f->GetObject("deltaT_YZ_h2_err0", deltaT_YZ_h2);
    f->GetObject("deltaT_YZ_h2_err0_neg", deltaT_YZ_h2_neg);

    float med, err;
    // beam left
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            med = deltaT_YZ_h2->GetBinContent(i+1,j+1);
            err = deltaT_YZ_h2->GetBinError(i+1,j+1);

            if (check_out(nbinsY, nbinsZ, deltaT_YZ_h2, i, j, 10))
            {
                // rm outliers
                if (i==1 && med < 4600) continue;
                if (j==29 && med < 4591) continue;
                if (j==0) continue;
                if (j==1 && med < 4595) continue;    
            }
            h->SetBinContent(i+1,j+1,med);
            h->SetBinError(i+1,j+1,err);
        }
    }

    // beam right
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            med = deltaT_YZ_h2_neg->GetBinContent(i+1,j+1);
            err = deltaT_YZ_h2_neg->GetBinError(i+1,j+1);
            if (check_out(nbinsY, nbinsZ, deltaT_YZ_h2_neg, i, j, 10))
            {
                if (i==1) continue;
                if (i==34 && med > 4600) continue;
                if (j==29) continue;
                if (j==28 && med < 4599) continue;
                if (j==1 && (med < 4594 || med > 4600)) continue; 
            } 
            h_neg->SetBinContent(i+1,j+1,med);
            h_neg->SetBinError(i+1,j+1,err);
            cout<<"err="<<err<<" seterr="<<h_neg->GetBinError(i+1,j+1,err)<<endl;
        }
    }
/*
    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    h->SetMarkerSize(0.7);
    h->Draw("COLZ TEXT");
    h->GetZaxis()->SetRangeUser(4575,4620);
    c1->SaveAs("TH2_dT_unordered_dT_without_thermal_rm_outlier.png");
    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    h_neg->SetMarkerSize(0.7);
    h_neg->Draw("COLZ TEXT");
    h_neg->GetZaxis()->SetRangeUser(4575,4620);
    c2->SaveAs("TH2_dT_unordered_dT_without_thermal_rm_outlier_neg.png");
 */

    TFile *file = new TFile("TH2_dT_unordered_dT_without_thermal_rm_outlier.root", "recreate");
    h->Write();
    h_neg->Write();
    file->Write();
    file->Close();

}
