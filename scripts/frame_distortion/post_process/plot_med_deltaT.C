// helpers
double get_hist_med (TH1F *h)
{
        if (h->GetEntries() == 0) return 0;
        double x, q;
        q = 0.5;
        h->ComputeIntegral();
        h->GetQuantiles(1, &x, &q);
        return x;
}

void get_med_error (TH1F *h, double *Tlow, double *Thigh)
{
        int size = h->GetEntries();
        double sigma = 0.5/sqrt(size+2);
        double qlow = 0.5 - sigma;
        double qhigh = 0.5 + sigma; 
        h->ComputeIntegral();
        h->GetQuantiles(1, Tlow, &qlow);
        h->GetQuantiles(1, Thigh, &qhigh);
}

int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm

void plot_med_deltaT()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots";
    
    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/ALLEVENTS_deltaT_with_contraction_corr.root");

   TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int nbinsY = 30, nbinsZ = 36;

    int nentries_min = 1000;

    // selection range
    int deltaT_min = 4580;
    int deltaT_max = 4700;
    int deltaT_margin = 10; // ticks

    TH2F *deltaT_YZ_h2 = new TH2F("deltaT_YZ_h2","deltaT vs YZ bin (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_neg = new TH2F("deltaT_YZ_h2_neg","deltaT vs YZ bin (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    TH1F *h;

    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            //cout<<Form("deltaT_pos_%d_%d", i, j)<<endl;
            f->GetObject(Form("deltaT_pos_%d_%d", i, j), h);
            if (h->GetEntries() < nentries_min) continue;
            h->GetXaxis()->SetRange(deltaT_min, deltaT_max);
            int maxbin = h->GetMaximumBin();
            h->GetXaxis()->SetRange(maxbin-deltaT_margin, maxbin+deltaT_margin);
            double deltaT_med = get_hist_med(h);
            double Tlow = 0, Thigh = 0;
            get_med_error(h, &Tlow, &Thigh);
            //deltaT_YZ_h2->GetZaxis()->SetRangeUser(4578,4610);
            deltaT_YZ_h2->SetBinContent(i+1,j+1,deltaT_med);
            deltaT_YZ_h2->SetBinError(i+1,j+1,(Thigh-Tlow)/2);
        }
    }

    deltaT_YZ_h2->Draw("COLZ TEXT");
}
