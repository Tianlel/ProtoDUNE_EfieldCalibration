/* add beam left and beam right */

#include <fstream>
#include "post_cm_fun.h"

/* Plot dT vs Z & dT vs Y */ 
int YZ_deltaT_binsize = 2; // ticks
int YZ_deltaT_min = 4000; // ticks
int YZ_deltaT_max = 5000; // ticks

int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm

void add_two_sides()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    string OUTPUT_TREE_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";

    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/updateTH2.root", "READ");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int Tbinsize = 2;
    int nbinsY = 30, nbinsZ = 36;
    int deltaT_YZ_hists_num = nbinsY * nbinsZ;

    int nentries_min = 1000;

    // selection range
    int deltaT_min = 4580;
    int deltaT_max = 4700;
    int deltaT_margin = 10; // ticks

    TH2F *out = new TH2F("frameXoffset_pos","frame Xoffset (beam left); Z (cm); Y(cm); Xoffset (cm)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *out_neg = new TH2F("frameXoffset_neg","frame Xoffset (beam right); Z (cm); Y(cm); Xoffset (cm)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *sum = new TH2F("sum","deltaT sum; Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax); 

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
            out->SetBinContent(i,j,med);
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
            if (med!=0) 
            {
                out_neg->SetBinContent(i,j,med);
            }
        }
    }

    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            float pos = out->GetBinContent(i+1,j+1);
            float neg = out_neg->GetBinContent(i+1,j+1);
            if (pos!=0 && neg!=0) sum->SetBinContent(i+1,j+1, pos+neg);

        }
    }
    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    sum->SetMarkerSize(0.7);
    sum->Draw("COLZ TEXT");
    sum->GetZaxis()->SetRangeUser(9160,9220);

 /*
    //c1->SaveAs("frameDisplacementTH2.png");
    deltaT_YZ_h2->SetMarkerSize(0.7);
    deltaT_YZ_h2->Draw("COLZ TEXT");
    deltaT_YZ_h2->GetZaxis()->SetRangeUser(4575,4620);
    //c1->SaveAs(save_name+".png");
    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    out_neg->SetMarkerSize(0.7);
    out_neg->Draw("COLZ TEXT");
    out_neg->GetZaxis()->SetRangeUser(-2.5,1.5);
    c2->SaveAs("frameDisplacementTH2_neg.png");

*/
}
