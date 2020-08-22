/* Subtract frame_displacement map from Mike's bkwd XOffset map to get SCE map */

float get_x_displacement_mike(float x, float y, float z, TH3F *map)
{   
    return map->GetBinContent(map->FindBin(x,y,z));
}

void get_SCEmap()
{
    TStyle *st = new TStyle("Modern","my style");
    st->SetPadGridX(1);
    st->SetPadGridY(1);
    st->cd();
    gStyle->SetOptStat(0);
    
    int save_png = 1;
    TFile *f0 = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/distortion_maps/SCE_DataDriven_180kV_v4.root");
    TFile *f1 = TFile::Open("/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/distortion_maps/frameXoffset_unordered_dT_without_thermal.root");

    TH3F *dismap0, *dismap0_neg;
    f0->GetObject("RecoBkwd_Displacement_X_Pos", dismap0);
    f0->GetObject("RecoBkwd_Displacement_X_Neg", dismap0_neg);

    TH2F *dismap1, *dismap1_neg;
    f1->GetObject("frameXoffset_pos", dismap1);
    f1->GetObject("frameXoffset_neg", dismap1_neg);

    int Ymin = 0, Ymax = 600; // cm
    int Zmin = -10, Zmax = 710; // cm
    int binsize = 20;
    int nbinsZ = (Zmax-Zmin)/binsize;
    int nbinsY = (Ymax-Ymin)/binsize;

    TH2F *map = new TH2F("SCE_map","Xoffset (SCE - CPA displacement) (beam left); Z (cm); Y(cm); displacement (cm)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax); 
    TH2F *map_neg = new TH2F("SCE_map","Xoffset (SCE - CPA displacement) (beam right); Z (cm); Y(cm); displacement (cm)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);

    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            float y = Ymin+10+j*20, z = Zmin+10+i*20;
            if (y<50 || y>550 || z<50 || z>580) continue;
            float dis0 = get_x_displacement_mike(0.15, y, z, dismap0);
            float dis1 = dismap1->Interpolate(y,z);
            float dif = dis0 - dis1;
            map->SetBinContent(i,j,dif);
             
            // beam right
            dis0 = get_x_displacement_mike(-0.15, y, z, dismap0_neg);
            cout<<"y: "<<y<<" z: "<<z<<endl;
            dis1 = dismap1_neg->Interpolate(y,z);
            dif = dis0 - dis1;
            map_neg->SetBinContent(i,j,dif);
        }
    }

    gStyle->SetPaintTextFormat(".2f");
    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 4000);
    map->SetMarkerSize(0.7);
    map->Draw("COLZ TEXT");
    if (save_png) c1->SaveAs("SCEmap.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 4000, 4000);
    map_neg->SetMarkerSize(0.7);
    map_neg->Draw("COLZ TEXT");
    if (save_png) c2->SaveAs("SCEmap_neg.png");


}
