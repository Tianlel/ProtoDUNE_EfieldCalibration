/* get deltaT vs Y && deltaT vs Z */

#include <fstream>
#include "post_cm_fun.h"

void create_n_histsY(int n, TH1F *hists_pos[n], TH1F *hists_neg[n],
                    char name[], char hist_title[], char x_unit[], char y_unit[],
                    int x0, int x1, int nbinsX)
{   
    int ybin = 0, zbin = 0, y_range_min = 0, z_range_min = 0;
    for (int i=0; i<n; i++)
    {   
        hists_pos[i] = new TH1F(Form("%s_pos_%d", name, i),
                            Form("%s (Beam Left), Y = %d cm; %s; %s", hist_title, 10+20*i,
                                 x_unit, y_unit), nbinsX, x0, x1);
        hists_neg[i] = new TH1F(Form("%s_neg_%d", name, i),
                            Form("%s (Beam Right), Y = %d cm; %s; %s", hist_title, 10+20*i,
                                 x_unit, y_unit), nbinsX, x0, x1);
    }
}


void create_n_histsZ(int n, TH1F *hists_pos[n], TH1F *hists_neg[n],
                    char name[], char hist_title[], char x_unit[], char y_unit[],
                    int x0, int x1, int nbinsX)
{    
    int ybin = 0, zbin = 0, y_range_min = 0, z_range_min = 0;
    for (int i=0; i<n; i++)
    {
        hists_pos[i] = new TH1F(Form("%s_pos_%d", name, i),
                            Form("%s (Beam Left), Z = %d cm; %s; %s", hist_title, 20*i,
                                 x_unit, y_unit), nbinsX, x0, x1);
        hists_neg[i] = new TH1F(Form("%s_neg_%d", name, i),
                            Form("%s (Beam Right), Z = %d cm; %s; %s", hist_title, 20*i,
                                 x_unit, y_unit), nbinsX, x0, x1);
    }
}



/* Plot dT vs Z & dT vs Y */ 
int YZ_deltaT_binsize = 2; // ticks
int YZ_deltaT_min = 4000; // ticks
int YZ_deltaT_max = 5000; // ticks

int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm

void plot_med_deltaT_02()
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

    int nbinsT_YZ = (YZ_deltaT_max - YZ_deltaT_min) / YZ_deltaT_binsize;
    TH1F *deltaTZ[nbinsY], *deltaTZ_neg[nbinsY];
    char deltaTZ_name[] = "deltaTZ",
         deltaTZ_title[] = "deltaT vs Z";
    char deltaTZ_xunit[] = "Z (cm)", deltaTZ_yunit[] = "median deltaT (ticks)";
    create_n_histsY(nbinsY, deltaTZ, deltaTZ_neg,
                   deltaTZ_name, deltaTZ_title,
                   deltaTZ_xunit, deltaTZ_yunit,
                   -10, Zmax, nbinsZ);

    TH1F *deltaTY[nbinsZ], *deltaTY_neg[nbinsZ];
    char deltaTY_name[] = "deltaTY",
         deltaTY_title[] = "deltaT vs Y";
    char deltaTY_xunit[] = "Y (cm)", deltaTY_yunit[] = "median deltaT (ticks)";
    create_n_histsZ(nbinsZ, deltaTY, deltaTY_neg,
                   deltaTY_name, deltaTY_title,
                   deltaTY_xunit, deltaTY_yunit,
                   0, Ymax, nbinsY);

    float val, yval, zval;
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)   
        {
            val = out->GetBinContent(i+1,j+1);
            zval = out->GetXaxis()->GetBinCenter(i+1);
            yval = out->GetYaxis()->GetBinCenter(j+1);
            
            //cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<endl;

            deltaTY[i]->Fill(yval,val);
            deltaTZ[j]->Fill(zval,val);
            if (i==10 && j == 20) deltaTZ[i]->Draw("hist");

            val = out_neg->GetBinContent(i+1,j+1);
            zval = out_neg->GetXaxis()->GetBinCenter(i+1);
            yval = out_neg->GetYaxis()->GetBinCenter(j+1);
            
            deltaTY_neg[i]->Fill(yval,val);
            deltaTZ_neg[j]->Fill(zval,val);
        }
    }
/*
    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    out->SetMarkerSize(0.7);
    out->Draw("COLZ TEXT");
    out->GetZaxis()->SetRangeUser(-2,2);
    c1->SaveAs("frameDisplacementTH2.png");
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
    TFile *file = new TFile("deltaTvsY_deltaTvsZ.root", "recreate");
    for (int i=0; i<nbinsZ; i++) 
    {
        deltaTY[i]->Write();
        deltaTY[i]->GetYaxis()->SetRangeUser(4550,4630);
        deltaTY_neg[i]->Write();
        deltaTY_neg[i]->GetYaxis()->SetRangeUser(4550,4630);
    }
    for (int i=0; i<nbinsY; i++) 
    {
        deltaTZ[i]->Write();
        deltaTZ_neg[i]->Write();
    }

    file->Write();
}
