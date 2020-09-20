TCanvas* getp1(TFile *f)
{
    TH1F *h, *h_neg;
    f->GetObject("deltaT_pos_23_25", h);
    f->GetObject("deltaT_neg_23_25", h_neg);

    h->SetTitle("dT distribution (beam left, Z: 450-470 cm, Y: 500-520 cm); dT (ticks); Counts");
    h_neg->SetTitle("dT distribution (beam right, Z: 450-470 cm, Y: 500-520 cm); dT (ticks); Counts");

    TCanvas *c = new TCanvas("c","c",1200,500);
    c->Divide(2,1);

    c->cd(1);
    h->Draw();

    c->cd(2);
    h_neg->Draw();

    return c;
}

TCanvas* getp2(TFile *f)
{
    TH1F *h, *h_neg;
    f->GetObject("deltaT_pos_20_17", h);
    f->GetObject("deltaT_neg_20_17", h_neg);

    h->SetTitle("dT distribution (beam left, Z: 390-410 cm, Y: 340-360 cm); dT (ticks); Counts");
    h_neg->SetTitle("dT distribution (beam right, Z: 390-410 cm, Y: 340-360 cm); dT (ticks); Counts");

    TCanvas *c = new TCanvas("c","c",1200,500);
    c->Divide(2,1);

    c->cd(1);
    h->Draw();

    c->cd(2);
    h_neg->Draw();
    return c;
}

void get_dT_distribution()
{
    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/ALLEVENTS_without_thermal_unordered_dT.root", "READ");

    TCanvas *c = getp1(f);
    c->SaveAs("./plots/p02_dT_two_peak_feature.png");

    TCanvas *c2 = getp2(f);
    c2->SaveAs("./plots/p01_dT_one_peak_feature.png");
}
