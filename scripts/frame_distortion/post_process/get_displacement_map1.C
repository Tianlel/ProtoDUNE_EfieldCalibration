#include <fstream>

/* For processing the root file with deltaT 2D histograms */

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

void get_displacement_map1()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    string OUTPUT_TREE_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    TString save_name = TString::Format("%supdated_deltaT_med_TH2", OUTPUT_PATH.c_str());
    TString save_name_tree = TString::Format("%s", OUTPUT_TREE_PATH.c_str()); 

    TFile *f = TFile::Open("TH2_dT_unordered_dT_without_thermal_rm_outlier.root", "READ");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int Tbinsize = 2;
    int nbinsY = 30, nbinsZ = 36;
    int nentries_min = 1000;

    TH2F *out = new TH2F("frameXoffset_pos","frame Xoffset (beam left); Z (cm); Y(cm); Xoffset (cm)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *out_neg = new TH2F("frameXoffset_neg","frame Xoffset (beam right); Z (cm); Y(cm); Xoffset (cm)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    TH2F *deltaT_YZ_h2, *deltaT_YZ_h2_neg;
    f->GetObject("dT_rmBoundary_pos", deltaT_YZ_h2);
    f->GetObject("dT_rmBoundary_neg", deltaT_YZ_h2_neg);

    float med, err;
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            med = deltaT_YZ_h2->GetBinContent(i+1,j+1);
            err = deltaT_YZ_h2->GetBinError(i+1,j+1);
            if (med!=0)
                out->SetBinContent(i,j,get_distortion(med,0));
        }
    }
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            med = deltaT_YZ_h2_neg->GetBinContent(i+1,j+1);
            err = deltaT_YZ_h2_neg->GetBinError(i+1,j+1);
            if (med!=0) 
                out_neg->SetBinContent(i,j,get_distortion(med,1));
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    out->SetMarkerSize(0.7);
    out->Draw("COLZ TEXT");
    out->GetZaxis()->SetRangeUser(-2,2);
  //  c1->SaveAs("frameDisplacementTH2_unordered_with_thermal.png");
  //  deltaT_YZ_h2->SetMarkerSize(0.7);
//    deltaT_YZ_h2->Draw("COLZ TEXT");
 //   deltaT_YZ_h2->GetZaxis()->SetRangeUser(4575,4620);
    //c1->SaveAs(save_name+".png");

    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    out_neg->SetMarkerSize(0.7);
    out_neg->Draw("COLZ TEXT");
    out_neg->GetZaxis()->SetRangeUser(-2.5,1.5);
//    c2->SaveAs("frameDisplacementTH2_unordered_with_thermal_neg.png");

    TFile *file = new TFile("frameXoffset_unordered_dT_without_thermal.root", "recreate");
    out->Write();
    out_neg->Write();
    file->Write();
    file->Close();

}
