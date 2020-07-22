// helpers
double get_hist_med (TH1F *h, int maxbin, int bin_margin, double *Tlow, double *Thigh)
{
    double cnt = 0, num = 0;
    vector<double> dTs;
    cout<<maxbin+bin_margin<<endl;
    for (int i=maxbin-bin_margin; i<=maxbin+bin_margin; i++)
    {
       // cout<<"i: "<<i<<endl; 
        cnt = h->GetBinContent(i);
        num = h->GetBinCenter(i);
        dTs.insert(dTs.end(), cnt, num);
    }
    int size = dTs.size();
    double* dT = &dTs[0];
    if (!size) return 0;

    double qout[3];
    double p[3];
    p[0]=0.5-0.5/sqrt(size+2);
    p[1]=0.5;
    p[2]=0.5+0.5/sqrt(size+2);
    TMath::Quantiles(size,3,&dT[0],qout,p,0,0,8);
    cout<<"p:"<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
    cout<<"qout: "<<qout[0]<<" "<<qout[1]<<" "<<qout[2]<<endl;
    *Tlow = qout[1]-qout[0];
    *Thigh = qout[2]-qout[1];
    return qout[1];
}

int YZ_deltaT_min = 4000;
int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm

void plot_med_deltaT()
{
    string PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/frame_distortion/";
    string OUTPUT_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/plots";
    TString save_name = TString::Format("%sdeltaT_med_TH2", OUTPUT_PATH.c_str());
    
    TFile *f = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/ALLEVENTS_deltaT_with_contraction_corr.root");

   TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);

    int Tbinsize = 2;
    int nbinsY = 30, nbinsZ = 36;

    int nentries_min = 1000;

    // selection range
    int deltaT_min = 4580;
    int deltaT_max = 4700;
    int deltaT_margin = 10; // ticks

    TH2F *deltaT_YZ_h2 = new TH2F("deltaT_YZ_h2","deltaT vs YZ bin (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_neg = new TH2F("deltaT_YZ_h2_neg","deltaT vs YZ bin (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    TH1F *h;

    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            //cout<<Form("deltaT_pos_%d_%d", i, j)<<endl;
            f->GetObject(Form("deltaT_pos_%d_%d", i, j), h);
            if (h->GetEntries() < nentries_min) continue;
            h->GetXaxis()->SetRangeUser(deltaT_min, deltaT_max);
            int maxbin = h->GetMaximumBin();
            int maxbin_T = YZ_deltaT_min + Tbinsize*(maxbin-1);
            //cout<<maxbin_T<<endl;
            double Tlow = 0, Thigh = 0;
            int deltaT_med = (int)get_hist_med(h, maxbin, deltaT_margin/Tbinsize, &Tlow, &Thigh);
            //cout<<maxbin_T<<"T_med: "<<deltaT_med<<endl;
            //cout<<"Tlow: "<<Tlow<<" Thigh:"<<Thigh<<endl;
            deltaT_YZ_h2->GetZaxis()->SetRangeUser(4578,4610);
            deltaT_YZ_h2->SetBinContent(i+1,j+1,deltaT_med);
            deltaT_YZ_h2->SetBinError(i+1,j+1,(Thigh-Tlow)/2);
    
            f->GetObject(Form("deltaT_neg_%d_%d", i, j), h);
            if (h->GetEntries() < nentries_min) continue;
            h->GetXaxis()->SetRangeUser(deltaT_min, deltaT_max);
            maxbin = h->GetMaximumBin();
            maxbin_T = YZ_deltaT_min + Tbinsize*(maxbin-1);
            deltaT_med = (int)get_hist_med(h, maxbin, deltaT_margin/Tbinsize, &Tlow, &Thigh);
            deltaT_YZ_h2_neg->GetZaxis()->SetRangeUser(4578,4610);
            deltaT_YZ_h2_neg->SetBinContent(i+1,j+1,deltaT_med);
            deltaT_YZ_h2_neg->SetBinError(i+1,j+1,(Thigh-Tlow)/2);


        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    deltaT_YZ_h2->SetMarkerSize(0.7);
    deltaT_YZ_h2->Draw("COLZ TEXT");
    c1->SaveAs(save_name+".png");

    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    deltaT_YZ_h2_neg->SetMarkerSize(0.7);
    deltaT_YZ_h2_neg->Draw("COLZ TEXT");
    c2->SaveAs(save_name+"_neg.png");
}
