/* cmp Ajin's 2D deltaT map and mine by subtraction */
#include "post_cm_fun.h"

void cmp_2D()
{
    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int Ymin = 0, Ymax = 600; // cm
    int Zmin = -10, Zmax = 710; // cm 
    int nbinsY = 30, nbinsZ = 36;
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_TREE_PATH = "./";
    TString save_name_tree = TString::Format("%s", OUTPUT_TREE_PATH.c_str());

    TH2F *d = new TH2F("diff","diff (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *d_neg = new TH2F("diff_neg","diff (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    TH2F *pos0, *neg0, *pos1, *neg1;
    TFile *fa = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/ALLEVENTS_with_thermal_unordered_dT.root", "READ");
    fa->GetObject("deltaT_YZ_h2_err0", pos0);
    fa->GetObject("deltaT_YZ_h2_err0_neg", neg0);

    TFile *fb = TFile::Open("/dune/app/users/apaudel/MC_sample_jan26/anode-anode_MC/anode-anode_aug16/data_YZoffset_error_v2_thermalcontraction_july30.root", "READ");
    fb->GetObject("deltaT_poshist", pos1);
    fb->GetObject("deltaT_neghist", neg1);    

    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            float diff = pos0->GetBinContent(i,j)-pos1->GetBinContent(i,j);
            if (diff > -100 && diff < 100)
                d->SetBinContent(i,j, diff);

            diff = neg0->GetBinContent(i,j)-neg1->GetBinContent(i,j);
            if (diff > -100 && diff < 100)
                d_neg->SetBinContent(i,j, diff);
        }
    }

    gStyle->SetPaintTextFormat(".1f");
    TCanvas *c1 = new TCanvas("c1","c1", 1000, 2000); 
    d->SetMarkerSize(0.7);
    d->Draw("COLZ TEXT");
    draw_Yframe();
    draw_Zframe();
    c1->SaveAs("method_diff.png");

    TCanvas *c2 = new TCanvas("c2","c2", 1000, 2000); 
    d_neg->SetMarkerSize(0.7);
    d_neg->Draw("COLZ TEXT");
    draw_Yframe();
    draw_Zframe();
    c2->SaveAs("method_diff_neg.png");


}
