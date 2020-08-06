#include <fstream>

/* check distribution for particular bins */

void exam_TH1_distribution()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/updateTH2.root", "READ");

    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    TH1F *h; 
    f->GetObject("deltaT_neg_14_28", h);
    h->Draw();
    
}
