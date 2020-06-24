#define Efield_cxx
#include "Efield.h"
#include <iostream>
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
#include "common_funcs.h"

using namespace std;

/********** constant variables ***********/
float xbegin = 0;
float xend = 3586;
double Temp = 87.67; // K

float distance = xend - xbegin; 

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
double sizebinsizeT = 20;

int nbinsX = 20;

float Tbinsize_pos = Tmax_pos / nbinsT;
float Tbinsize_neg = Tmax_neg / nbinsT; 
float Xbinsize = xend / nbinsX;

int dT_hist_nbins = 5000;
int dT_hist_xmin = 0;
int dT_hist_xmax = 5000;

/* output file */
string output_file_name = "3586_plots.root";

/* optional variables */
string output_tree_name = "tree_name";

Long64_t select_nentries = 0;

/* style selection */
SetGrids(0);
RmOptStat(0);


/* cut selection */

/* operation selection */
int trk_rm_dupzHits = 1; // if a trk has two hits with the same z, remove both of them from the trk 

/******* to be added *********/

/******************* constant variables end **********************/

/********* struct declaration *********/
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
 
struct Track {
    vector<Hit> hits;
}   

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

/******************* helper functions end **********************/


/*********** main section of code for Efield calculation **********/
void Efield::Loop()
{

    TFile *file = new TFile(output_file_name.c_str(), "recreate");
    
    /*** vector declaration ***/
    vector<Float_t> hitpeakT_buffer, hity_buffer, hitz_buffer;
    vector<Int_t> hittpc_buffer, hitwire_buffer;
    vector<Float_t> hitpeakT_buffer_neg, hity_buffer_neg, hitz_buffer_neg;
    vector<Int_t> hittpc_buffer_neg, hitwire_buffer_neg;
  
    vector<Float_t> trkhitT[nbinsT], trkhitX[nbinsT], trkhitT_neg[nbinsT], trkhitX_neg[nbinsT];
    vector<Float_t> trkhitTx[nbinsX], trkhitXx[nbinsX], trkhitTx_neg[nbinsX], trkhitXx_neg[nbinsX];

    vector<vector<double>> drift_velocity_fit, drift_velocity_fit_neg;
    drift_velocity_fit.resize(nbinsT); drift_velocity_fit_neg.resize(nbinsT);
    vector<vector<double>> drift_velocity_fitx, drift_velocity_fit_negx;
    drift_velocity_fitx.resize(nbinsX); drift_velocity_fit_negx.resize(nbinsX);  

    /********* vector declaration end *********/
 
    /*** histogram definition ***/
    TH1F *hist_dT_after_Tcut_volcut = new TH1F("hist_dT_after_Tcut_volcut", "peakT_max - peakT_min -- beam left;T (ticks);number of trks", dT_hist_nbins, dT_hist_xmin, dT_hist_xmax); // distribution of selected trks after hit peak time cut & volume cut 
    TH1F *hist_dT_after_Tcut_volcut_neg = new TH1F("hist_dT_after_Tcut_volcut_neg", "peakT_max - peakT_min -- beam right;ticks;number of trks", dT_hist_nbins, dT_hist_xmin, dT_hist_xmax);
    
    /* drift velociy histograms */
    float varbins[nbinsT+1], varbins_neg[nbinst+1];
    SetVariableBins(varbins, nbinsT, sidebinsizeT, Tmax_pos);
    SetVariableBins(varbins_neg, nbinsT, sidebinsizeT, Tmax_neg);

    TH1F *hist_driftvel_med_fit = new TH1F("hist_driftvel_med_fit", "Median drift velocity beam left (fit per trk);T (ticks);drift velocity (mm/µs)", nbinsT, varbins); // median drift velocity from fitting x vs T. 
    TH1F *hist_driftvel_med_fit_neg = new TH1F("hist_driftvel_med_fit_neg", "Median drift velocity for beam right (fit per trk);T (ticks);drift velocity (mm/µs)", mbinsT, varbins_neg);
    
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
            hist_driftvel_med_fit.clear(); hist_driftvel_med_fit_neg.clear();
            hitpeakT_buffer.clear(); hitpeakT_buffer_neg.clear();  
            hittpc_buffer.clear(); hittpc_buffer_neg.clear();  
            hitwire_buffer.clear(); hitwire_buffer_neg.clear();  
            hitz_buffer.clear(); hitz_buffer_neg.clear();
            hity_buffer.clear(); hity_buffer_neg.clear();
            for (int k = 0; k < nbinsT; k++)
            { 
                    trkhitT[k].clear(); trkhitT_neg[k].clear();
                    trkhitX[k].clear(); trkhitX_neg[k].clear();
            }
            for (int k = 0; k < nbinsX; k++)
            {
                    trkhitTx[k].clear(); trkhitT_negx[k].clear();
                    trkhitXx[k].clear(); trkhitX_negx[k].clear();
            }              

            /****** loop through all hits of this current trk ******/
            for (size_t j = 0; j < hit_peakT2->at(i).size(); j++)
            {
                    /****** selet hits located on beam LEFT based on tpc number ******/
                    if (hit_tpc2->at(i)[j] == 2 || hit_tpc2->at(i)[j] == 6 || hit_tpc2->at(i)[j] == 10)
                    {
                        /**** remove hits with coinciding z ****/
                        int remove_hit = 0;
                        if trk_rm_dupzHits 
                        {
                            vector<float>::iterator it;
///////////////////////////// ??? need to define z_vals etc. //////////////////
                            it = find(z_vals.begin(), z_vals.end(), trkhitz_wire2->at(i)[j]);
                            if (it != z_vals.end()) remove_hit = 1
                        } 
                        if (j==0 || remove_hit==0){
                            z_vals.push_back(trkhitz_wire2->at(i)[j]);
                            hitpeakT_buffer.push_back(hit_peakT2->at(i)[j]);
                            hittpc_buffer.push_back(hit_tpc2->at(i)[j]);
                            hitwire_buffer.push_back(hit_wire2->at(i)[j]);
                            hity_buffer.push_back(trkhity2->at(i)[j]);
                            hitz_buffer.push_back(trkhitz_wire2->at(i)[j]);
                        }
                    } // end of if #tpc beam left#

                    /****** selet hits located on beam RIGHT based on tpc number ******/
                    if (hit_tpc2->at(i)[j] == 1 || hit_tpc2->at(i)[j] == 5 || hit_tpc2->at(i)[j] == 9)
                    {
                        /**** remove hits with coinciding z ****/
                        int remove_hit = 0;
                        if trk_rm_dupzHits 
                        {   
                            vector<float>::iterator it;
                            it = find(z_vals_neg.begin(), z_vals_neg.end(), trkhitz_wire2->at(i)[j]);
                            if (it != z_vals_neg.end()) remove_hit = 1
                        }   
                        if (j==0 || remove_hit==0)
                        {
                            z_vals_neg.push_back(trkhitz_wire2->at(i)[j]);
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
        if (siz > trksize_min)
        {
            vector<Hit> hits;
            hits.reserve( siz ); 
            
            for (int i=0; i<siz; i++)
            {
                hits.push_back( Hit( hitpeakT_buffer[i], hity_buffer[i], hitz_buffer[i],
                                     9999., hittpc_buffer[i], hitwire_buffer[i] ) );
            }
            sort_hits_by_T(hits);
        }
