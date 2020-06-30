#define Efield_cxx
#include "Efield.h"
#include <iostream>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <numeric>
#include <TSpline.h>
#include <TStyle.h>
#include <assert.h>
#include "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/scripts/common_funcs.h"

using namespace std;

/********** constant variables ***********/
float xbegin = 0;
float xend = 3586;
double Temp = 87.67; // K

float distanceX = xend - xbegin; 

/* hit peak time cut (select out cathode-anode crossers) */
float Tmin_pos = 4565;
float Tmax_pos = 4578;
float Tmin_neg = 4565;
float Tmax_neg = 4578;

/* volume cut */
float miny = 200;
float maxy = 400;
float minz = 240;
float maxz = 450;

/* trk size cut */
int trksize_min = 45; 

/* histogram parameters */
int nbinsT = 27;
double sidebinsizeT = 20;

int nbinsX = 20;

float Tbinsize_pos = (Tmax_pos-2*sidebinsizeT) / (nbinsT-2);
float Tbinsize_neg = (Tmax_neg-2*sidebinsizeT) / (nbinsT-2); 
float Xbinsize = xend / nbinsX;

int dT_hist_nbins = 5000;
int dT_hist_xmin = 0;
int dT_hist_xmax = 5000;

/* output file */
string output_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/";
string output_file_name = "3586_plots.root";

/* optional variables */
string output_tree_name = "tree_name";

Long64_t select_nentries = 0;

/* style selection */
//SetGrids(0); to be added ??????
//RmOptStat(0);

/* counters */
int selected_trk_count = 0, selected_trk_count_neg = 0;


/* cut selection */

/* operation selection */
int trk_rm_dupzHits = 1; // if a trk has two hits with the same z, remove both of them from the trk 

/******* to be added *********/

/* linear fit parameters */
int fit_hitnum_per_bin_min = 3; 

/* tpc number (LArSoft numbering) */ 
int tpc_left_front = 2, tpc_left_mid = 6, tpc_left_back = 10;
int tpc_right_front = 1, tpc_right_mid = 5, tpc_right_back = 9;
/******************* constant variables end **********************/

/********* struct definition *********/
struct Hit {
    float peakT;
    float y;
    float z;
    float x_calculated;
    int tpc;
    int wire;

    /*** constructors ***/
    Hit() : peakT( 0.0f ), y( 0.0f ), z( 0.0f ), x_calculated( 0.0f ),
            tpc( 0 ), wire( 0 ) {}
    Hit( float peakTIn, float yIn, float zIn, float x_calcIn, int tpcIn, int wireIn) :
        peakT( peakTIn ), y( yIn ), z( zIn ), x_calculated( x_calcIn ),
        tpc( tpcIn ), wire( wireIn ) {}
};
 
/********** helper functions **********/

TSpline3 *sp = create_ve_sp(Temp); // ??? converting drift velocity to Efield

void sort_hits_by_T(vector<Hit> &hits)
{   
    sort(hits.begin(), hits.end(), [&](const auto& hit1, const auto& hit2)
    {   
        return hit1.peakT < hit2.peakT;
    });
}

/* given sidebin size, nbins, and xmax, evenly divides the middle bins 
 *  * returns a list of floats that can be passed to the histogram constructor */
void SetVariableBins (float varbins[], int nbins, float side_bin_size, float xmax)
{

    float middle_bin_size = (xmax - 2*side_bin_size) / (nbins - 2);

    varbins[0] = 0;
    varbins[1] =  varbins[0] + side_bin_size;
    varbins[nbins] = xmax;
    varbins[nbins - 1] = xmax - side_bin_size;

    for (int i=0; i<nbins-3; i++)
    {
        varbins[i+2] = varbins[i+1] + middle_bin_size;
    }
    return;
}

// ?????????? need to be tested 
/* returns 1 if the track lies in the specified cut volume */
int inCutVolume (vector<Hit> hits, int miny, int maxy, int minz, int maxz)
{
    int checkz = 0, checky = 0, ret = 0;
    int size = hits.size();
    float z0 = hits[0].z, z1 = hits[size-1].z;
    float y0 = hits[0].y, y1 = hits[size-1].y;
    if (z0 > minz && z0 < maxz && z1 > minz && z1 < maxz) checkz = 1;
    if (y0 > miny && y0 < maxy && y1 > miny && y1 < maxy) checky = 1;
    ret = checkz * checky;
    return ret;
}

float get_fit_velocity(vector<Hit> hits)
{
    int hits_num = hits.size();
    // write T & x into vectors
    vector<float> hitT, hitX;
    for (int i=0; i<hits_num; i++)
    {
        hitT.push_back(hits[i].peakT);
        hitX.push_back(hits[i].x_calculated);
    }
    TGraph *vfit = new TGraph(hits_num, &hitT[0], &hitX[0]);
    TF1 *f = new TF1("linearfit", "pol1", 0, 4500); //linear fit
    vfit->Fit("linearfit","Q");
    float v = -2 * (f->GetParameter(1));
    return v;
}

// need to be tested
/* given the fitted velocities for each bin, 
 *  * return the median velocities and standard error*/
void get_Vmedian_per_bin (vector<vector<float>> drift_velocity, int nbins, vector<float> medianVs, vector<float> standard_error)
{
    for (int i=0; i<nbins; i++)
    {
        int vsize = drift_velocity[i].size();
        if (vsize == 0) 
        {
            medianVs.push_back(-1);
            standard_error.push_back(0);           
        }
        sort(drift_velocity[i].begin(), drift_velocity[i].end());
        float vmed = TMath::Median(drift_velocity[i].size(), &drift_velocity[i][0]);
        float trunc_rms = TMath::RMS(drift_velocity[i].begin()+(int)(vsize*0.1), drift_velocity[i].end()-(int)(vsize*0.1));
        float se = 0.31*trunc_rms/sqrt(0.8*vsize);
        medianVs.push_back(vmed);
        standard_error.push_back(se);
    }
}

float get_trunc_se (vector<float> Es)
{
    int binsize = Es.size();
    float trunc_rms = TMath::RMS(Es.begin()+(int)(binsize*0.1), Es.end()-(int)(binsize*0.1));
    float se = trunc_rms/sqrt(0.8*binsize);
    return se;
}

void fill_histograms (vector<vector<float>> values, vector<vector<float>> errors, vector<TH1F*> hists, vector<int> nbins)
{
    int nhists = hists.size();
    for (int i=0; i<nhists; i++)
    {
        for (int j=0; j<nbins[i]; j++)
        {
            hists[i]->SetBinContent(j+1, values[i][j]);
            hists[i]->SetBinError(j+1, errors[i][j]);
        }   
    }
}

/******************* helper functions end **********************/


/*********** main section of code for Efield calculation **********/
void Efield::Loop()
{
    string output = output_PATH + output_file_name;
    TFile *file = new TFile(output.c_str(), "recreate");
    
    /*** vector declaration ***/
    vector<Float_t> hitpeakT_buffer, hity_buffer, hitz_buffer;
    vector<Int_t> hittpc_buffer, hitwire_buffer;
    vector<Float_t> hitpeakT_buffer_neg, hity_buffer_neg, hitz_buffer_neg;
    vector<Int_t> hittpc_buffer_neg, hitwire_buffer_neg;
    vector<float> z_record, z_record_neg;
  
    vector<Float_t> trkhitT[nbinsT], trkhitX[nbinsT], trkhitT_neg[nbinsT], trkhitX_neg[nbinsT];
    vector<Float_t> trkhitTx[nbinsX], trkhitXx[nbinsX], trkhitTx_neg[nbinsX], trkhitXx_neg[nbinsX];

    vector<Float_t> model_Efield_T[nbinsT], model_Efield_T_neg[nbinsT];
    vector<Float_t> model_Efield_X[nbinsX], model_Efield_X_neg[nbinsX];

    vector<vector<float>> drift_velocityT, drift_velocityT_neg;
    drift_velocityT.resize(nbinsT); drift_velocityT_neg.resize(nbinsT);
    vector<vector<float>> drift_velocityX, drift_velocityX_neg;
    drift_velocityX.resize(nbinsX); drift_velocityX_neg.resize(nbinsX);  

    /********* vector declaration end *********/
 
    /*** histogram definition ***/
    TH1F *hist_dT_after_Tcut_volcut = new TH1F("hist_dT_after_Tcut_volcut", "peakT_max - peakT_min -- beam left;T (ticks);number of trks", dT_hist_nbins, dT_hist_xmin, dT_hist_xmax); // distribution of selected trks after hit peak time cut & volume cut 
    TH1F *hist_dT_after_Tcut_volcut_neg = new TH1F("hist_dT_after_Tcut_volcut_neg", "peakT_max - peakT_min -- beam right;ticks;number of trks", dT_hist_nbins, dT_hist_xmin, dT_hist_xmax);
    
    /* drift velociy histograms */
    float varbins[nbinsT+1], varbins_neg[nbinsT+1];
    SetVariableBins(varbins, nbinsT, sidebinsizeT, Tmax_pos);
    SetVariableBins(varbins_neg, nbinsT, sidebinsizeT, Tmax_neg);

    TH1F *hist_driftvel_med_fit_T = new TH1F("hist_driftvel_med_fit_T", "Median drift velocity beam left (fit per trk);T (ticks);drift velocity (mm/µs)", nbinsT, varbins); // median drift velocity from fitting x vs T. 
    TH1F *hist_driftvel_med_fit_T_neg = new TH1F("hist_driftvel_med_fit_T_neg", "Median drift velocity for beam right (fit per trk);T (ticks);drift velocity (mm/µs)", nbinsT, varbins_neg);
    
    TH1F *hist_driftvel_med_fit_x = new TH1F("hist_driftvel_med_fit_x", "Median drift velocity beam left (fit per trk);x (mm);drift velocity (mm/µs)", nbinsX, 0, xend); // median drift velocity from fitting x vs T. 
    TH1F *hist_driftvel_med_fit_x_neg = new TH1F("hist_driftvel_med_fit_x_neg", "Median drift velocity for beam right (fit per trk);x (mm);drift velocity (mm/µs)", nbinsX, 0, xend);


    /* Efield histograms */
    TH1F *hist_efield_med_fit_T = new TH1F("hist_efield_med_fit_T", "Median Electric Field beam left;T (ticks);Electric field (kV/cm)", nbinsT, varbins); // plotted against T
    TH1F *hist_efield_med_fit_T_neg = new TH1F("hist_efield_med_fit_T_neg", "Median Electric Field beam right;T (ticks);Electric field (kV/cm)", nbinsT, varbins_neg);

    TH1F *hist_efield_med_fit_x = new TH1F("hist_efield_med_fit_x", "Median Electric Field beam left;x (mm);Electric field (kV/cm)", nbinsX, 0, xend); // plotted against X
    TH1F *hist_efield_med_fit_x_neg = new TH1F("hist_efield_med_fit_x_neg", "Median Electric Field beam right;x (mm);Electric field (kV/cm)", nbinsX, 0, xend);

    /* Efield model output histograms */
    TH1F *hist_LAr_efield_T=new TH1F("hist_LAr_efield_T","LArSoft output Electric field beam left; T (ticks);Electric field (kV/cm)", nbinsT, varbins); // plotted against T
    TH1F *hist_LAr_efield_T_neg=new TH1F("hist_LAr_efield_T_neg","LArSoft output Electric field beam left; T (ticks);Electric field (kV/cm)", nbinsT, varbins_neg);
    
    TH1F *hist_LAr_efield_X=new TH1F("hist_LAr_efield_X","LArSoft output Electric field beam left; x (mm);Electric field (kV/cm)", nbinsT, varbins); // plotted ahainst X
    TH1F *hist_LAr_efield_X_neg=new TH1F("hist_LAr_efield_X_neg","LArSoft output Electric field beam left; T (mm);Electric field (kV/cm)", nbinsT, varbins_neg);


    /****** Get hit info from input tree ******/
    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntries();
    if (select_nentries != 0) nentries = select_nentries; 

    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
                break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        int no_femb = 0;
        for (int f1 = 0; f1 < 6; f1++)
        {
            no_femb = no_femb + Nactivefembs[f1];
        }
        if (jentry % 10000 == 0)
            std::cout << jentry << "/" << nentries << std::endl;
        if (!hit_peakT2->size())
            continue;
        /****** loop through all the trks ******/
        for (size_t i = 0; i < hit_peakT2->size(); i++)
        {
            if (!hit_peakT2->at(i).size()) continue;

            /****** clear vectors for each trk ******/
            hitpeakT_buffer.clear(); hitpeakT_buffer_neg.clear();  
            hittpc_buffer.clear(); hittpc_buffer_neg.clear();  
            hitwire_buffer.clear(); hitwire_buffer_neg.clear();  
            hitz_buffer.clear(); hitz_buffer_neg.clear();
            hity_buffer.clear(); hity_buffer_neg.clear();
            drift_velocityT.clear(); drift_velocityT_neg.clear();
            drift_velocityX.clear(); drift_velocityX_neg.clear();
            for (int k = 0; k < nbinsT; k++)
            { 
                    trkhitT[k].clear(); trkhitT_neg[k].clear();
                    trkhitX[k].clear(); trkhitX_neg[k].clear();
            }
            for (int k = 0; k < nbinsX; k++)
            {
                    trkhitTx[k].clear(); trkhitTx_neg[k].clear();
                    trkhitXx[k].clear(); trkhitXx_neg[k].clear();
            }              

            /****** loop through all hits of this current trk ******/
            for (size_t j = 0; j < hit_peakT2->at(i).size(); j++)
            {
                    /****** selet hits located on beam LEFT based on tpc number ******/
                    if (hit_tpc2->at(i)[j] == tpc_left_front || hit_tpc2->at(i)[j] == tpc_left_mid || hit_tpc2->at(i)[j] == tpc_left_back)
                    {
                        /**** remove hits with coinciding z ****/
                        int remove_hit = 0;
                        if (trk_rm_dupzHits) 
                        {
                            vector<float>::iterator it;
///////////////////////////// ??? need to define z_record etc. //////////////////
                            it = find(z_record.begin(), z_record.end(), trkhitz_wire2->at(i)[j]);
                            if (it != z_record.end()) remove_hit = 1;
                        } 
                        if (j==0 || remove_hit==0){
                            z_record.push_back(trkhitz_wire2->at(i)[j]);
                            hitpeakT_buffer.push_back(hit_peakT2->at(i)[j]);
                            hittpc_buffer.push_back(hit_tpc2->at(i)[j]);
                            hitwire_buffer.push_back(hit_wire2->at(i)[j]);
                            hity_buffer.push_back(trkhity2->at(i)[j]);
                            hitz_buffer.push_back(trkhitz_wire2->at(i)[j]);
                        }
                    } // end of if #tpc beam left#

                    /****** selet hits located on beam RIGHT based on tpc number ******/
                    if (hit_tpc2->at(i)[j] == tpc_right_front || hit_tpc2->at(i)[j] == tpc_right_mid || hit_tpc2->at(i)[j] == tpc_right_back)
                    {
                        /**** remove hits with coinciding z ****/
                        int remove_hit = 0;
                        if (trk_rm_dupzHits==1)
                        {   
                            vector<float>::iterator it;
                            it = find(z_record_neg.begin(), z_record_neg.end(), trkhitz_wire2->at(i)[j]);
                            if (it != z_record_neg.end()) remove_hit = 1;
                        }   
                        if (j==0 || remove_hit==0)
                        {
                            z_record_neg.push_back(trkhitz_wire2->at(i)[j]);
                            hitpeakT_buffer_neg.push_back(hit_peakT2->at(i)[j]);
                            hittpc_buffer_neg.push_back(hit_tpc2->at(i)[j]);
                            hitwire_buffer_neg.push_back(hit_wire2->at(i)[j]);
                            hity_buffer_neg.push_back(trkhity2->at(i)[j]);
                            hitz_buffer_neg.push_back(trkhitz_wire2->at(i)[j]);
                        }   
                    } // end of if #tpc beam right#
/////////////////////////// ????????????? //////////////////////////////////////////
 
            } // end of hit loop #j
        } // end of trk loop #i

        size_t trksize = hitpeakT_buffer.size(); // beam left
        if (trksize > trksize_min)
        {
            vector<Hit> hits;
            hits.reserve(trksize); 
            
            for (int i=0; i<trksize; i++)
            {
                hits.push_back( Hit( hitpeakT_buffer[i], hity_buffer[i], hitz_buffer[i],
                                     9999., hittpc_buffer[i], hitwire_buffer[i] ) );
            }
            sort_hits_by_T(hits);

            float Tmin = hits[0].peakT;
            float Tmax = hits[trksize-1].peakT;
            if (Tmax < Tmin) exit(0);

            // fill dT histogram after volume cut
            int incutvol = inCutVolume(hits, miny, maxy, minz, maxz);
            if (incutvol==1)
                    hist_dT_after_Tcut_volcut->Fill(Tmax - Tmin);

            // auxilliary cuts ???????? to be added

            if (incutvol==1)
            {
                selected_trk_count++;

                float z0 = hits[0].z, z1 = hits[trksize-1].z;
                float y0 = hits[0].y, y1 = hits[trksize-1].y;

                vector<Hit> hits_per_Tbin[nbinsT], hits_per_Xbin[nbinsX];
                // hits loop
                for (int hit_itr = 0; hit_itr < trksize; hit_itr++)
                {
                    float T = hits[hit_itr].peakT;
                    float y = hits[hit_itr].y;
                    float z = hits[hit_itr].z;
                    float x = xbegin + distanceX * (abs(z1 - z) / abs(z1 - z0)); // mm
                    float x_cm = x/10; 

                    float hitdeltaT = T - Tmin;
                    // update hit information
                    hits[hit_itr].peakT = hitdeltaT;
                    hits[hit_itr].x_calculated = x; 

                    int Tbin_num;
                    if (hitdeltaT < 0) continue;
                    else if (hitdeltaT > Tmax - Tmin) continue;
                    else if (hitdeltaT < sidebinsizeT) Tbin_num = 0; 
                    else if (hitdeltaT > Tmax - sidebinsizeT) Tbin_num = nbinsT - 1;
                    else Tbin_num = 1 + (int)(hitdeltaT - sidebinsizeT)/Tbinsize_pos;
                    hits_per_Tbin[Tbin_num].push_back(hits[hit_itr]);

                    model_Efield_T[Tbin_num].push_back(Ef_ex(x_cm,y,z)); // need to double check?????? 

                    int Xbin_num = x/Xbinsize;
                    if (Xbin_num >= nbinsX) continue;
                    hits_per_Xbin[Xbin_num].push_back(hits[hit_itr]);

                    model_Efield_X[Xbin_num].push_back(Ef_ex(x_cm,y,z));

                } // end of hits loop

                // loop through time bins to get local velocity (cm/s) from linear fit
                for (int Tbin_itr = 0; Tbin_itr < nbinsT; Tbin_itr++)
                {
                    if (hits_per_Tbin[Tbin_itr].size() < fit_hitnum_per_bin_min) continue;
                    sort_hits_by_T(hits_per_Tbin[Tbin_itr]); 

                    float v = get_fit_velocity(hits_per_Tbin[Tbin_itr]);
                    drift_velocityT[Tbin_itr].push_back(v);  
                } // end of time bin loop

                // loop through X bins to get local velocity (cm/s) from linear fit
                for (int Xbin_itr = 0; Xbin_itr < nbinsX; Xbin_itr++)
                {
                    if (hits_per_Xbin[Xbin_itr].size() < fit_hitnum_per_bin_min) continue;
                    sort_hits_by_T(hits_per_Xbin[Xbin_itr]);

                    float v = get_fit_velocity(hits_per_Xbin[Xbin_itr]);
                    drift_velocityX[Xbin_itr].push_back(v);
                } // end of X bin loop
            } // end of track selection
        } // end of track size cut
    
        trksize = hitpeakT_buffer_neg.size(); // beam right
        if (trksize > trksize_min)
        {
            vector<Hit> hits;
            hits.reserve(trksize); 
            
            for (int i=0; i<trksize; i++)
            {
                hits.push_back( Hit( hitpeakT_buffer_neg[i], hity_buffer_neg[i], hitz_buffer_neg[i],
                                     9999., hittpc_buffer_neg[i], hitwire_buffer_neg[i] ) );
            }
            sort_hits_by_T(hits);

            float Tmin = hits[0].peakT;
            float Tmax = hits[trksize-1].peakT;
            if (Tmax < Tmin) exit(0);

            // fill dT histogram after volume cut
            int incutvol = inCutVolume(hits, miny, maxy, minz, maxz);
            if (incutvol==1)
                    hist_dT_after_Tcut_volcut->Fill(Tmax - Tmin);

            // auxilliary cuts ???????? to be added

            if (incutvol==1)
            {
                selected_trk_count_neg++;

                float z0 = hits[0].z, z1 = hits[trksize-1].z;
                float y0 = hits[0].y, y1 = hits[trksize-1].y;

                vector<Hit> hits_per_Tbin[nbinsT], hits_per_Xbin[nbinsX];
                // hits loop
                for (int hit_itr = 0; hit_itr < trksize; hit_itr++)
                {
                    float T = hits[hit_itr].peakT;
                    float y = hits[hit_itr].y;
                    float z = hits[hit_itr].z;
                    float x = xbegin + distanceX * (abs(z1 - z) / abs(z1 - z0)); // mm
                    float x_cm = x/10; 

                    float hitdeltaT = T - Tmin;
                    // update hit information
                    hits[hit_itr].peakT = hitdeltaT;
                    hits[hit_itr].x_calculated = x; 

                    int Tbin_num;
                    if (hitdeltaT < 0) continue;
                    else if (hitdeltaT > Tmax - Tmin) continue;
                    else if (hitdeltaT < sidebinsizeT) Tbin_num = 0; 
                    else if (hitdeltaT > Tmax - sidebinsizeT) Tbin_num = nbinsT - 1;
                    else Tbin_num = 1 + (int)(hitdeltaT - sidebinsizeT)/Tbinsize_pos;
                    hits_per_Tbin[Tbin_num].push_back(hits[hit_itr]);

                    model_Efield_T[Tbin_num].push_back(Ef_ex(x_cm,y,z)); // need to double check?????? 

                    int Xbin_num = x/Xbinsize;
                    if (Xbin_num >= nbinsX) continue;
                    hits_per_Xbin[Xbin_num].push_back(hits[hit_itr]);

                    model_Efield_X[Xbin_num].push_back(Ef_ex(x_cm,y,z));

                } // end of hits loop

                // loop through time bins to get local velocity (cm/s) from linear fit
                for (int Tbin_itr = 0; Tbin_itr < nbinsT; Tbin_itr++)
                {
                    if (hits_per_Tbin[Tbin_itr].size() < fit_hitnum_per_bin_min) continue;
                    sort_hits_by_T(hits_per_Tbin[Tbin_itr]); 

                    float v = get_fit_velocity(hits_per_Tbin[Tbin_itr]);
                    drift_velocityT_neg[Tbin_itr].push_back(v);  
                } // end of time bin loop

                // loop through X bins to get local velocity (cm/s) from linear fit
                for (int Xbin_itr = 0; Xbin_itr < nbinsX; Xbin_itr++)
                {
                    if (hits_per_Xbin[Xbin_itr].size() < fit_hitnum_per_bin_min) continue;
                    sort_hits_by_T(hits_per_Xbin[Xbin_itr]);

                    float v = get_fit_velocity(hits_per_Xbin[Xbin_itr]);
                    drift_velocityX_neg[Xbin_itr].push_back(v);
                } // end of X bin loop
            } // end of track selection
        } // end of track size cut

    } // end of track loop

    /****** beam left ******/
    vector<float> medianVs, standard_error; 
    get_Vmedian_per_bin(drift_velocityT, nbinsT, medianVs, standard_error); 

    /****** beam right ******/
    vector<float> medianVs_neg, standard_error_neg;
    get_Vmedian_per_bin(drift_velocityT_neg, nbinsT, medianVs_neg, standard_error_neg);

    for (int i=0; i<nbinsT; i++)
    {
        /****** beam left ******/
        hist_driftvel_med_fit_T->SetBinContent(i+1, medianVs[i]);
        hist_driftvel_med_fit_T->SetBinError(i+1, standard_error[i]);
        
        float medianE = sp->Eval(medianVs[i]);
        hist_efield_med_fit_T->SetBinContent(i+1, sp->Eval(medianE));
        hist_efield_med_fit_T->SetBinError(i+1, 0.31*medianE/sqrt(0.8*drift_velocityT.size()));

        float model_binsize = model_Efield_T[i].size();
        float model_medianE = TMath::Median(model_binsize, &model_Efield_T[i][0]);
        float model_se = get_trunc_se(model_Efield_T[i]); 
        hist_LAr_efield_T->SetBinContent(i+1, model_medianE);
        hist_LAr_efield_T->SetBinError(i+1, model_se);

        /****** beam right ******/ 
        hist_driftvel_med_fit_T_neg->SetBinContent(i+1, medianVs_neg[i]);
        hist_driftvel_med_fit_T_neg->SetBinError(i+1, standard_error_neg[i]); 

        float medianE_neg = sp->Eval(medianVs_neg[i]);
        hist_efield_med_fit_T_neg->SetBinContent(i+1, sp->Eval(medianE_neg)); 
        hist_efield_med_fit_T_neg->SetBinError(i+1, 0.31*medianE_neg/sqrt(0.8*drift_velocityT_neg.size())); 

        model_binsize = model_Efield_T_neg[i].size();
        model_medianE = TMath::Median(model_binsize, &model_Efield_T_neg[i][0]);
        model_se = get_trunc_se(model_Efield_T_neg[i]); 
        hist_LAr_efield_T_neg->SetBinContent(i+1, model_medianE);
        hist_LAr_efield_T_neg->SetBinError(i+1, model_se);

    }
    
    vector<float> medianVs_x, standard_error_x;
    get_Vmedian_per_bin(drift_velocityX, nbinsX, medianVs_x, standard_error_x);

    vector<float> medianVs_x_neg, standard_error_x_neg;
    get_Vmedian_per_bin(drift_velocityX, nbinsX, medianVs_x_neg, standard_error_x_neg);

    for (int i=0; i<nbinsX; i++)
    {

        /****** beam left ******/
        hist_driftvel_med_fit_x->SetBinContent(i+1, medianVs_x[i]);
        hist_driftvel_med_fit_x->SetBinError(i+1, standard_error_x[i]);
        
        float medianE_x = sp->Eval(medianVs_x[i]);
        hist_efield_med_fit_x->SetBinContent(i+1, sp->Eval(medianE_x));
        hist_efield_med_fit_x->SetBinError(i+1, 0.31*medianE_x/sqrt(0.8*drift_velocityX.size()));

        //hist_LAr_efield_x->SetBinContent(i+1, model_Efield_x[i]);

        /****** beam right ******/ 
        hist_driftvel_med_fit_x_neg->SetBinContent(i+1, medianVs_x_neg[i]);
        hist_driftvel_med_fit_x_neg->SetBinError(i+1, standard_error_x_neg[i]); 

        float medianE_x_neg = sp->Eval(medianVs_x_neg[i]);
        hist_efield_med_fit_x_neg->SetBinContent(i+1, sp->Eval(medianE_x_neg)); 
        hist_efield_med_fit_x_neg->SetBinError(i+1, 0.31*medianE_x_neg/sqrt(0.8*drift_velocityX_neg.size())); 

        //hist_LAr_efield_x_neg->SetBinContent(i+1, model_Efield_x_neg[i]);
    }

    cout << "selected trk count in positive drift " << selected_trk_count << endl;
    cout << "selected trk count in negative drift " << selected_trk_count_neg << endl;

    hist_driftvel_med_fit_T->Write(); hist_driftvel_med_fit_T_neg->Write();
    hist_efield_med_fit_T->Write(); hist_efield_med_fit_T_neg->Write();
    
    hist_driftvel_med_fit_x->Write(); hist_driftvel_med_fit_x_neg->Write();
    hist_efield_med_fit_x->Write(); hist_efield_med_fit_x_neg->Write(); 
    
    file->Write();
    file->Close();
} // end of track loop

            


