#include <fstream>

/* For processing the root file with 2D histograms */

int YZ_deltaT_min = 4000;
int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm

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

double X0 = 358.6; // cm
double driftv = 0.16; //cm/us

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

void plot_med_deltaT_01()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    string OUTPUT_TREE_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    TString save_name = TString::Format("%supdated_deltaT_med_TH2", OUTPUT_PATH.c_str());
    TString save_name_tree = TString::Format("%s", OUTPUT_TREE_PATH.c_str()); 

    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/updateTH2.root", "READ");

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
    int deltaT_margin = 10; // ticks

    TH2F *out = new TH2F("out","frame distortion (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *out_neg = new TH2F("clone_neg","frame distortion (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    TH2F *deltaT_YZ_h2, *deltaT_YZ_h2_neg;
    f->GetObject("deltaT_YZ_h2_err0", deltaT_YZ_h2);
    f->GetObject("deltaT_YZ_h2_err0_neg", deltaT_YZ_h2_neg);
    float med, err;
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            med = deltaT_YZ_h2->GetBinContent(i+1,j+1);
            err = deltaT_YZ_h2->GetBinError(i+1,j+1);

            if (check_out(nbinsY, nbinsZ, deltaT_YZ_h2, i, j, 10))
            {
                if (j==1 || j==29) continue;
            }
            out->SetBinContent(i,j,get_distortion(med,0));
        }
    }
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            med = deltaT_YZ_h2_neg->GetBinContent(i+1,j+1);
            err = deltaT_YZ_h2_neg->GetBinError(i+1,j+1);
            if (check_out(nbinsY, nbinsZ, deltaT_YZ_h2_neg, i, j, 15))
            {
                if (med<4590) continue;
            } 
            float dis = get_distortion(med,1);
            if (med!=0) 
            {
                cout<<"bin:"<<i<<" "<<j<<endl;
                out_neg->SetBinContent(i,j,dis);
            }
        }
    }
/*
    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    out->SetMarkerSize(0.7);
    out->Draw("COLZ TEXT");
    out->GetZaxis()->SetRangeUser(7,11);
    c1->SaveAs("frameDisplacementTH2.png");

    deltaT_YZ_h2->SetMarkerSize(0.7);
    deltaT_YZ_h2->Draw("COLZ TEXT");
    deltaT_YZ_h2->GetZaxis()->SetRangeUser(4575,4620);
    //c1->SaveAs(save_name+".png");
*/
    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    out_neg->SetMarkerSize(0.7);
    out_neg->Draw("COLZ TEXT");
    out_neg->GetZaxis()->SetRangeUser(-11,-7);
    c2->SaveAs("frameDisplacementTH2_neg.png");
}
