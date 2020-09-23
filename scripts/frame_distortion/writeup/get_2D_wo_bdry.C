
void get_2D_wo_bdry()
{
    string INPUT_FILE = "TH2_dT_unordered_dT_without_thermal_rm_outlier.root";
    string PATH = "../../frame_distortion/post_process/rootfiles_after_rm_bdry/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots/";
    
    TFile *f = TFile::Open((PATH+INPUT_FILE).c_str(), "READ");
    //TFile *froot = TFile::Open(save_name_tree + "deltaTmed_v2.root", "RECREATE");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    TH2F *deltaT_YZ_h2, *deltaT_YZ_h2_neg;
    f->GetObject("dT_rmBoundary_pos", deltaT_YZ_h2);
    f->GetObject("dT_rmBoundary_neg", deltaT_YZ_h2_neg);

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
    deltaT_YZ_h2->SetMarkerSize(0.7);
    deltaT_YZ_h2->SetTitle("beam left (x>0);z (cm);y (cm); dT (ticks)");
    deltaT_YZ_h2->Draw("COLZ TEXT");
    deltaT_YZ_h2->GetZaxis()->SetTitleOffset(0);
    deltaT_YZ_h2->GetZaxis()->SetRangeUser(4575,4620);
    c1->SaveAs("p06_2D_dT_med.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 1000);
    deltaT_YZ_h2_neg->SetMarkerSize(0.7);
    deltaT_YZ_h2_neg->SetTitle("beam right (x<0) ;z (cm);y (cm); dT (ticks)");
    deltaT_YZ_h2_neg->Draw("COLZ TEXT");
    deltaT_YZ_h2_neg->GetZaxis()->SetRangeUser(4575,4620);
    deltaT_YZ_h2_neg->GetZaxis()->SetTitleOffset(0);
    c2->SaveAs("p07_2D_dT_med_neg.pdf");

}
