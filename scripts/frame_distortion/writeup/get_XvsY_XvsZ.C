/* get deltaT vs Y && deltaT vs Z */

#include <fstream>

TCanvas* getp1(TFile *f)
{
    TH1F *h, *h_neg, *h_avg;
    f->GetObject("XoffsetY_pos_17", h);
    f->GetObject("XoffsetY_neg_17", h_neg);
    f->GetObject("XoffsetY_avg_17", h_avg);

    h->SetTitle("Xoffset vs Y (Z: 330-350 cm); Y (cm); Xoffset (ticks)");
    h_neg->SetTitle("Xoffset vs Y (beam right, Z: 330-350 cm); Y (cm); Xoffset (ticks)");

    TCanvas *c = new TCanvas("c","c",600,500);

    h->Draw();
    h->SetLineColor(kRed);
    h->SetLineWidth(2);
    h->GetYaxis()->SetRangeUser(-5,5);

    h_neg->Draw("SAME");
    h_neg->SetLineColor(kBlue);
    h_neg->SetLineWidth(2);
    h_neg->GetYaxis()->SetRangeUser(-5,5);

    h_avg->Draw("SAME");
    h_avg->SetLineColor(kBlack);
    h_avg->SetLineWidth(2);
    h_avg->GetYaxis()->SetRangeUser(-5,5);

    auto *l = new TLegend(0.1,0.7,0.3,0.9);
    l->AddEntry(h,"beam left");
    l->AddEntry(h_neg, "beam right");
    l->AddEntry(h_avg, "average");
    l->Draw("SAME");
    return c;
}

TCanvas* getp2(TFile *f)
{
    TH1F *h, *h_neg, *h_avg;
    f->GetObject("XoffsetZ_pos_15", h);
    f->GetObject("XoffsetZ_neg_15", h_neg);
    f->GetObject("XoffsetZ_avg_15", h_avg);

    h->SetTitle("Xoffset vs Z (Y: 300-320 cm); Z (cm); Xoffset (ticks)");
    h_neg->SetTitle("Xoffset vs Z (beam right, Y: 300-320 cm); Z (cm); Xoffset (ticks)");

    TCanvas *c = new TCanvas("c0","c0",600,500);

    h->Draw();
    h->SetLineColor(kRed);
    h->SetLineWidth(2);
    h->GetYaxis()->SetRangeUser(-5,5);

    h_neg->Draw("SAME");
    h_neg->SetLineColor(kBlue);
    h_neg->SetLineWidth(2);
    h_neg->GetYaxis()->SetRangeUser(-5,5);

    h_avg->Draw("SAME");
    h_avg->SetLineColor(kBlack);
    h_avg->SetLineWidth(2);
    h_avg->GetYaxis()->SetRangeUser(-5,5);

    auto *l = new TLegend(0.1,0.7,0.3,0.9);
    l->AddEntry(h,"beam left");
    l->AddEntry(h_neg, "beam right");
    l->AddEntry(h_avg, "average");
    l->Draw("SAME");

    return c;
}



void get_XvsY_XvsZ()
{
    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/post_process/rootfiles_dTvsY_dTvsZ/XvsY_XvsZ_unordered_wo_thermal.root", "READ");

    TCanvas *c = getp1(f);
    c->SaveAs("./plots/p12_XoffsetvsY.pdf");

    TCanvas *c2 = getp2(f);
    c2->SaveAs("./plots/p13_XoffsetvsZ.pdf");





}
