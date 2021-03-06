void plot_TH2()
{
    int save_png = 0;


    string INPUT_FILE = "ALLEVENTS_without_thermal_unordered_dT.root";
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    
    TFile *f = TFile::Open((PATH+INPUT_FILE).c_str(), "READ");
    //TFile *froot = TFile::Open(save_name_tree + "deltaTmed_v2.root", "RECREATE");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    TH2F *deltaT_YZ_h2, *deltaT_YZ_h2_neg;
    f->GetObject("deltaT_YZ_h2_err0", deltaT_YZ_h2);
    f->GetObject("deltaT_YZ_h2_err0_neg", deltaT_YZ_h2_neg);

    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    deltaT_YZ_h2->SetMarkerSize(0.7);
    deltaT_YZ_h2->Draw("COLZ TEXT");
    deltaT_YZ_h2->GetZaxis()->SetRangeUser(4575,4620);
    if (save_png) c1->SaveAs("beam_left_deltaTTH2_v2.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    deltaT_YZ_h2_neg->SetMarkerSize(0.7);
    deltaT_YZ_h2_neg->Draw("COLZ TEXT");
    deltaT_YZ_h2_neg->GetZaxis()->SetRangeUser(4575,4620);
    if (save_png) c2->SaveAs("beam_right_deltaTTH2_v2.png");

}
