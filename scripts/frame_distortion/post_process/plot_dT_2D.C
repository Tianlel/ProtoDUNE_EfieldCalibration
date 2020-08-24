/* get outlier bins and remove them accordingly */

#include <fstream>

int YZ_deltaT_min = 4000;
int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm
double X0 = 358.6; // cm

void plot_dT_2D()
{
    string FILE_NAME = "ALLEVENTS_without_thermal_unordered_dT.root";
    TString save_name = "TH2_dT_med_without_thermal_unordered_dT";

    string OPEN_FILE_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/";
    TFile *f = TFile::Open((OPEN_FILE_PATH+FILE_NAME).c_str(), "READ");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    TH2F *h, *h_neg;
    f->GetObject("deltaT_YZ_h2_err0", h);
    f->GetObject("deltaT_YZ_h2_err0_neg", h_neg);

    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    h->SetMarkerSize(0.7);
    h->Draw("COLZ TEXT");
    h->GetZaxis()->SetRangeUser(4575,4620);
    c1->SaveAs(save_name+".png");
 
    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    h_neg->SetMarkerSize(0.7);
    h_neg->Draw("COLZ TEXT");
    h_neg->GetZaxis()->SetRangeUser(4575,4620);
    c2->SaveAs(save_name+"_neg.png");
}
