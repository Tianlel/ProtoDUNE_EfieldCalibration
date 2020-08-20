
/* Author: Tianle Liu 
 */

#define frame_cxx
#include "frame.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "frame_func.h"

using namespace std;
   
/* inut file */
string input_file_name = "deltaTmed_v3.root";

/* output file */
string output_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/";
string output_file_name = "getmed_with_thermal_unordered_dT.root";

/* optional variables */
Long64_t select_nentries = 0; // 0 -> use all nentries
int set_jentry = 0;
int print_debug_message = 0;

/* cut parameters */
int hits_size_min = 5;
int CA_crossing_cut = 0;
int med_min_hits_num = 5;

int YZ_deltaT_getmed_lowerbound = 4580; // find peak after this threshold
int YZ_deltaT_entries_min = 1000; // min 

/* purpose of analysis */
int get_deltaT_distribution_before_cut = 0;

void frame::Loop()
{
    string output = output_PATH + output_file_name;
    TFile * froot = TFile::Open(input_file_name.c_str(), "READ");
    TFile *file = new TFile(output.c_str(), "recreate");
    vector<float> *Tmed, *Tmed_neg;
    froot->GetObject("med_vals", Tmed);
    froot->GetObject("med_vals_neg", Tmed_neg);

    /*** vector declaration ***/
    vector<Float_t> hitpeakT_buffer, hity_buffer, hitz_buffer;
    vector<Int_t> hittpc_buffer, hitwire_buffer;
    vector<Float_t> hitpeakT_buffer_neg, hity_buffer_neg, hitz_buffer_neg;
    vector<Int_t> hittpc_buffer_neg, hitwire_buffer_neg;

    vector<Track> selected_tracks, selected_tracks_neg; 


    /****** histogram definition ******/

    int nbinsT_YZ = (YZ_deltaT_max - YZ_deltaT_min) / YZ_deltaT_binsize; 
    cout<<"number of YZ_bins in deltaT: "<< nbinsT_YZ<<endl;
    char deltaT_hist_xunit[] = "deltaT_total (ticks)", deltaT_hist_yunit[] = "number of hits";

    /* for plotting deltaT vs Zframe (to get bounds for cathode-anode crosser cut) */
    int deltaT_zframe_hists_num = frame_Zpositions.size() - 1;
    TH1F *deltaT_zframe_hists[deltaT_zframe_hists_num], *deltaT_zframe_hists_neg[deltaT_zframe_hists_num]; 
    vector<const char*> deltaT_zframe_hists_names (deltaT_zframe_hists_num, "zframe_ACcrossers_deltaT_distribution_hists"),
         deltaT_zframe_hists_titles {"deltaT_total distribution (0<=z<120)",
                                     "deltaT_total distribution (120<=z<230)",
                                     "deltaT_total distribution (230<=z<345)",
                                     "deltaT_total distribution (345<=z<465)",
                                     "deltaT_total distribution (465<=z<575)",
                                     "deltaT_total distribution (575<=z<695)"};
   create_n_hists(deltaT_zframe_hists_num, deltaT_zframe_hists, deltaT_zframe_hists_neg,
                  deltaT_zframe_hists_names, deltaT_zframe_hists_titles,
                  deltaT_hist_xunit, deltaT_hist_yunit,
                  nbinsT_YZ, YZ_deltaT_min, YZ_deltaT_max);


    /* for plotting deltaT vs YZ plane */
    int nbinsY = Ymax / Ybinsize, nbinsZ = (Zmax - Zmin) / Zbinsize;
    int histnum_tot = (nbinsY + 1) * (nbinsZ + 1);
    vector<vector<double>> peak_vals, peak_vals_neg;
    peak_vals.resize(histnum_tot); peak_vals_neg.resize(histnum_tot);
    int deltaT_YZ_hists_num = nbinsY * nbinsZ;
    if (print_debug_message) cout<<"deltaT_YZ_hists_num = "<<deltaT_YZ_hists_num<<endl;
    TH1F *deltaT_YZ_hists[deltaT_YZ_hists_num], *deltaT_YZ_hists_neg[deltaT_YZ_hists_num];
    char deltaT_YZ_hist_name[] = "deltaT",
         deltaT_YZ_hist_title[] = "deltaT distribution (YZ plane bin)";
    create_n_hists_YZ(nbinsY, deltaT_YZ_hists_num, deltaT_YZ_hists, deltaT_YZ_hists_neg,
                   deltaT_YZ_hist_name, deltaT_YZ_hist_title, 
                   deltaT_hist_xunit, deltaT_hist_yunit,
                   YZ_deltaT_min, YZ_deltaT_max, nbinsT_YZ);
    /* for plotting deltaT vs YZ plane end*/

    /* for plotting deltaT_med vs YZ plane */
    TH2F *deltaT_YZ_h2_err0 = new TH2F("deltaT_YZ_h2_err0","deltaT vs YZ bin (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_err0_neg = new TH2F("deltaT_YZ_h2_err0_neg","deltaT vs YZ bin (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_err1 = new TH2F("deltaT_YZ_h2_err1","deltaT vs YZ bin (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_err1_neg = new TH2F("deltaT_YZ_h2_err1_neg","deltaT vs YZ bin (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_err2 = new TH2F("deltaT_YZ_h2_err2","deltaT vs YZ bin (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_err2_neg = new TH2F("deltaT_YZ_h2_err2_neg","deltaT vs YZ bin (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);


    /* deltaT vs Z for each bin */
/*    TH1F *deltaTZ[deltaT_YZ_hists_num], *deltaTZ_neg[deltaT_YZ_hists_num];
    char deltaTZ_name[] = "deltaTZ",
         deltaTZ_title[] = "deltaT vs Z";
    char deltaTZ_xunit[] = "", deltaTZ_yunit[] = "number of hits";
    create_n_hists_YZ(nbinsY, deltaT_YZ_hists_num, deltaTZ, deltaTZ_neg,
                   deltaTZ_name, deltaTZ_title,
                   deltaTZ_xunit, deltaTZ_yunit,
                   YZ_deltaT_min, YZ_deltaT_max, nbinsT_YZ); */

    /****** histogram definition ******/

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    if (select_nentries != 0) nentries = select_nentries;   

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=set_jentry; jentry<nentries;jentry++) // event loop
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if (jentry % 10000 == 0)
            std::cout << jentry << "/" << nentries << std::endl;

        if(!trkhitz_wire2->size()) continue;

        if (print_debug_message) print("entering trk loop");
        /* trk loop: ith track */
        for(size_t i=0;i<trkhitz_wire2->size();i++)
        {
            if (print_debug_message) print("declare hit vector");
            vector<Hit> hits, hits_neg; 

            if (print_debug_message) print("entering hit loop");
            /* hit loop: jth hit */
            for(size_t j=0;j<trkhitz_wire2->at(i).size();j++)
            {
                if (print_debug_message) print("filling hits vectors");
                /****** select hits located on beam LEFT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_left_front || hit_tpc2->at(i)[j] == tpc_left_mid || hit_tpc2->at(i)[j] == tpc_left_back)
                {
                    int frame_tag = get_frame_tag(trkhitz_wire2->at(i)[j]);
                    hits.push_back( Hit(hit_peakT2->at(i)[j], 0, 0,
                                       // trkhity2->at(i)[j], trkhitz_wire2->at(i)[j],
                                        ythermal(trkhity2->at(i)[j]),
                                        zthermal(trkhitz_wire2->at(i)[j], hit_tpc2->at(i)[j]),
                                        9999.,hit_tpc2->at(i)[j],
                                        frame_tag));  
                }
                /****** select hits located on beam RIGHT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_right_front || hit_tpc2->at(i)[j] == tpc_right_mid || hit_tpc2->at(i)[j] == tpc_right_back)
                {
                    int frame_tag = get_frame_tag(trkhitz_wire2->at(i)[j]);
                    hits_neg.push_back( Hit(hit_peakT2->at(i)[j], 0, 0, // trkhity2->at(i)[j], trkhitz_wire2->at(i)[j],
                                        ythermal(trkhity2->at(i)[j]),
                                        zthermal(trkhitz_wire2->at(i)[j],hit_tpc2->at(i)[j]),
                                        -9999.,hit_tpc2->at(i)[j],
                                        frame_tag));
                }
            } // end of hit loop

            if (hits.size() > hits_size_min)
            {
                if (print_debug_message) print("create track object");
                // create Track object
                Track trk = Track(hits);
                trk.sort0();
                float Tmin = trk.Tmin();
                trk.set_deltaT_local(Tmin);
                int frame_tag = trk.get_cathode_frame_tag();
                if (frame_tag == -1) exit(0);

                float trk_deltaT_total = trk.deltaT_total();
                deltaT_zframe_hists[frame_tag]->Fill(trk_deltaT_total);
                
                if (CA_crossing_cut)
                {
                    if (is_CA_crosser(trk, frame_tag, 0))
                        if (print_debug_message) print("push_back selected_tracks");
                        selected_tracks.push_back(trk);
                }
                else
                {
                    int ybin = trk.cathode_hit().y / Ybinsize;
                    int zbin = (trk.cathode_hit().z - Zmin) / Zbinsize;
                    if (ybin < 0 || ybin >= nbinsY || zbin < 0 || zbin >= nbinsZ) continue;
                    if (print_debug_message) print("Fill deltaT_YZ hist");
                    if (print_debug_message) cout<<"YZhist_num is "<<YZhist_num(ybin,zbin,nbinsZ)<<endl;
                    deltaT_YZ_hists[YZhist_num(zbin,ybin,nbinsY)]->Fill(trk.deltaT_total());
                    if (print_debug_message) cout<<"Fill deltaT vector"<<endl;
                    fill_deltaT_vec(trk.deltaT_total(), Tcut_margin, zbin, ybin, Tmed, peak_vals);
                }
            }

            if (hits_neg.size() > hits_size_min)
            {
                if (print_debug_message) print("create track object (neg)");
                // create Track object
                Track trk = Track(hits_neg);
                trk.sort0();
                float Tmin = trk.Tmin();
                trk.set_deltaT_local(Tmin);
                int frame_tag = trk.get_cathode_frame_tag();
                if (frame_tag == -1) exit(0);

                float trk_deltaT_total = trk.deltaT_total();
                deltaT_zframe_hists_neg[frame_tag]->Fill(trk_deltaT_total);

                if (CA_crossing_cut)
                {   
                    if (is_CA_crosser(trk, frame_tag, 1))
                        if (print_debug_message) print("push_back selected_tracks_neg");
                        selected_tracks_neg.push_back(trk);
                }
                else 
                {
                    int ybin = trk.cathode_hit().y / Ybinsize;
                    int zbin = (trk.cathode_hit().z - Zmin) / Zbinsize;
                    if (ybin < 0 || ybin >= nbinsY || zbin < 0 || zbin >= nbinsZ) continue;
                    if (print_debug_message) print("Fill deltaT_YZ hist");
                    deltaT_YZ_hists_neg[YZhist_num(zbin,ybin,nbinsY)]->Fill(trk.deltaT_total());
                    fill_deltaT_vec(trk.deltaT_total(), Tcut_margin, zbin, ybin, Tmed_neg, peak_vals_neg);
                }
            }  

            if (print_debug_message) print("exiting trk loop");
        } // end of trk loop
    } // end of jentries loop

    // fill deltaT distribution histograms
    for (int i=0; i<selected_tracks.size(); i++)
    {
        int ybin = selected_tracks[i].cathode_hit().y / Ybinsize;
        int zbin = selected_tracks[i].cathode_hit().z / Zbinsize;
        if (ybin < 0 || ybin > nbinsY || zbin < 0 || zbin > nbinsZ) continue;
        deltaT_YZ_hists[YZhist_num(zbin,ybin,nbinsY)]->Fill(selected_tracks[i].deltaT_total());  

    }
    for (int i=0; i<selected_tracks_neg.size(); i++)
    {   
        int ybin = selected_tracks_neg[i].cathode_hit().y / Ybinsize;
        int zbin = selected_tracks_neg[i].cathode_hit().z / Zbinsize;
        if (ybin < 0 || ybin > nbinsY || zbin < 0 || zbin > nbinsZ) continue;
        deltaT_YZ_hists_neg[YZhist_num(zbin,ybin,nbinsY)]->Fill(selected_tracks_neg[i].deltaT_total());
    }   

    for (int i=0; i<nbinsZ; i++) 
    {
        for (int j=0; j<nbinsY; j++) 
        {
            int binnum = YZhist_num(i,j,nbinsY);
            if (peak_vals[binnum].size() < med_min_hits_num) continue;
            float Tmedval, Terr_up, Terr_down;
            get_vec_med_err(peak_vals[binnum], &Tmedval, &Terr_up, &Terr_down);
            if (print_debug_message) cout<<"Tmedout : "<<Tmedval<<" Tmedin :"<<Tmed->at(YZhist_num(i,j,nbinsY));
            if (Tmedval==0) continue;
            deltaT_YZ_h2_err0->SetBinContent(i+1,j+1,Tmedval);
            deltaT_YZ_h2_err0->SetBinError(i+1,j+1,(Terr_up+Terr_down)/2);
            deltaT_YZ_h2_err1->SetBinContent(i+1,j+1,Tmedval);
            deltaT_YZ_h2_err1->SetBinError(i+1,j+1,Terr_up);
            deltaT_YZ_h2_err2->SetBinContent(i+1,j+1,Tmedval);
            deltaT_YZ_h2_err2->SetBinError(i+1,j+1,Terr_down);
        }
    }

    for (int i=0; i<nbinsZ; i++)
    {
        for (int j=0; j<nbinsY; j++)
        {
            int binnum = YZhist_num(i,j,nbinsY);
            if (peak_vals_neg[binnum].size() < med_min_hits_num) continue;
            float Tmedval, Terr_up, Terr_down;
            get_vec_med_err(peak_vals_neg[binnum], &Tmedval, &Terr_up, &Terr_down);
            if (Tmedval==0) continue;
            deltaT_YZ_h2_err0_neg->SetBinContent(i+1,j+1,Tmedval);
            deltaT_YZ_h2_err0_neg->SetBinError(i+1,j+1,(Terr_up+Terr_down)/2);
            deltaT_YZ_h2_err1_neg->SetBinContent(i+1,j+1,Tmedval);
            deltaT_YZ_h2_err1_neg->SetBinError(i+1,j+1,Terr_up);
            deltaT_YZ_h2_err2_neg->SetBinContent(i+1,j+1,Tmedval);
            deltaT_YZ_h2_err2_neg->SetBinError(i+1,j+1,Terr_down);
        }
    }
       
    for (int i=0; i<deltaT_YZ_hists_num; i++)
    {
        deltaT_YZ_hists[i]->Write(); deltaT_YZ_hists_neg[i]->Write();
    }

    deltaT_YZ_h2_err0->Write(); deltaT_YZ_h2_err0_neg->Write();
    deltaT_YZ_h2_err1->Write(); deltaT_YZ_h2_err1_neg->Write();
    deltaT_YZ_h2_err2->Write(); deltaT_YZ_h2_err2_neg->Write();

    for (int i=0; i<deltaT_zframe_hists_num; i++) {
        deltaT_zframe_hists[i]->Write(); deltaT_zframe_hists_neg[i]->Write();
    }


    cout<<"selected_trk_count_pos = "<<selected_tracks.size()<<endl
        <<"selected_trk_count_neg = "<<selected_tracks_neg.size()<<endl;
    
    file->Write();
    file->Close();
} 
