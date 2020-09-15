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

void plot_dTvsY_dTvsZ()
{
    int save_root_file = 1;
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    string OUTPUT_TREE_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";

    TFile *f = TFile::Open("../../distortion_maps/frameXoffset_unordered_dT_without_thermal.root", "READ");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int Tbinsize = 2;
    int nbinsY = 30, nbinsZ = 36;
    int deltaT_YZ_hists_num = nbinsY * nbinsZ;

    TH2F *h, *h_neg;
    f->GetObject("frameXoffset_pos", h);
    f->GetObject("frameXoffset_neg", h_neg);

    int nbinsT_YZ = (YZ_deltaT_max - YZ_deltaT_min) / YZ_deltaT_binsize;
    TH1F *deltaTZ[nbinsY], *deltaTZ_neg[nbinsY];
    char deltaTZ_name[] = "deltaTZ",
         deltaTZ_title[] = "Xoffset vs Z";
    char deltaTZ_xunit[] = "Z (cm)", deltaTZ_yunit[] = "displacement (cm)";
    create_n_histsY(nbinsY, deltaTZ, deltaTZ_neg,
                   deltaTZ_name, deltaTZ_title,
                   deltaTZ_xunit, deltaTZ_yunit,
                   -10, Zmax, nbinsZ);

    TH1F *deltaTY[nbinsZ], *deltaTY_neg[nbinsZ];
    char deltaTY_name[] = "deltaTY",
         deltaTY_title[] = "Xoffset vs Y";
    char deltaTY_xunit[] = "Y (cm)", deltaTY_yunit[] = "displacement (cm)";
    create_n_histsZ(nbinsZ, deltaTY, deltaTY_neg,
                   deltaTY_name, deltaTY_title,
                   deltaTY_xunit, deltaTY_yunit,
                   0, Ymax, nbinsY);

    float val, yval, zval, err;
    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)   
        {
            val = h->GetBinContent(i+1,j+1);
            err = h->GetBinError(i+1,j+1);
            zval = h->GetXaxis()->GetBinCenter(i+1);
            yval = h->GetYaxis()->GetBinCenter(j+1);
            
            deltaTY[i]->Fill(yval,val);
            deltaTY[i]->SetBinError(j+1,err);
            deltaTZ[j]->Fill(zval,val);
            deltaTZ[j]->SetBinError(i+1,err);

            val = h_neg->GetBinContent(i+1,j+1);
            err = h_neg->GetBinError(i+1,j+1);
            zval = h_neg->GetXaxis()->GetBinCenter(i+1);
            yval = h_neg->GetYaxis()->GetBinCenter(j+1);
            
            deltaTY_neg[i]->Fill(yval,val);
            deltaTY_neg[i]->SetBinError(j+1,err);
            deltaTZ_neg[j]->Fill(zval,val);
            deltaTZ_neg[j]->SetBinError(i+1,err);
        }
    }
    TCanvas *c = new TCanvas("c","c", 10000, 14000);
    c->Divide(6,6);
    TCanvas *c2 = new TCanvas("c2","c2", 10000, 14000);
    c2->Divide(5,6);
    
    TFile *file = new TFile("dXvsY_dXvsZ_unordered_wo_thermal.root", "recreate");
    for (int i=0; i<nbinsZ; i++) 
    {
        if (save_root_file) deltaTY[i]->Write();
        c->cd(i+1);
        deltaTY[i]->Draw();
        deltaTY[i]->SetLineColor(kRed);
        deltaTY[i]->GetYaxis()->SetRangeUser(-3,3);
    
        if (save_root_file) deltaTY_neg[i]->Write();
        deltaTY_neg[i]->Draw("SAME");
        deltaTY_neg[i]->SetLineColor(kBlue);
        deltaTY_neg[i]->GetYaxis()->SetRangeUser(-3,3);
    }
    for (int i=0; i<nbinsY; i++) 
    {
        if (save_root_file) deltaTZ[i]->Write();
        c2->cd(i+1);
        deltaTZ[i]->Draw();
        deltaTZ[i]->SetLineColor(kRed);
        deltaTZ[i]->GetYaxis()->SetRangeUser(-3,3);

        if (save_root_file) deltaTZ_neg[i]->Write();
        deltaTZ_neg[i]->Draw("SAME");
        deltaTZ_neg[i]->SetLineColor(kBlue);
        deltaTZ_neg[i]->GetYaxis()->SetRangeUser(-3,3);
    }

    c->SaveAs("./../../../plots/data/Xoffset_vs_Y_unordered_wo_thermal.pdf");
    c2->SaveAs("./../../../plots/data/Xoffset_vs_Z_unordered_wo_thermal.pdf");
    if (save_root_file) file->Write();
}
