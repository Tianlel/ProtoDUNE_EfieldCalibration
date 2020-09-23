TCanvas* getp1(TFile *f)
{
    TH1F *h1, *h2;
    f->GetObject("deltaT_neg_15_28", h1);
    f->GetObject("deltaT_neg_15_29", h2);

    h1->SetTitle("dT distribution (beam right, Z: 290-310 cm, Y: 560-580 cm); dT (ticks); Counts");
    h2->SetTitle("dT distribution (beam right, Z: 290-310 cm, Y: 580-600 cm); dT (ticks); Counts");

    TCanvas *c = new TCanvas("c","c",1200,500);
    c->Divide(2,1);

    c->cd(1);
    h1->Draw();
    h1->GetXaxis()->SetRangeUser(4200,4800);

    c->cd(2);
    h2->Draw();
    h2->GetXaxis()->SetRangeUser(4200,4800);
    return c;
}

void get_dT_distribution2()
{
    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/ALLEVENTS_without_thermal_unordered_dT.root", "READ");

    TCanvas *c = getp1(f);
    c->SaveAs("./plots/p05_dT_two_outlier_cmp.pdf");

}
