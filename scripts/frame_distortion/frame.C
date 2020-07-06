/* 
 * Author: Tianle Liu 
 */

#define frame_cxx
#include "frame.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;


/********** constant variables ***********/

/* output file */
string output_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/reco/frame_distortion/";
string output_file_name = "test_plots.root";

/* optional variables */
Long64_t select_nentries = 1000; // 0 -> use all nentries

/* cut parameters */
int hits_size_min = 0; 

/* histogram parameters */
int Ybinsize = 20; // cm
int Zbinsize = 20; // cm 
int Tbinsize = 50; // ticks 

int deltaT_min = 0; // ticks

/* detector parameters */
int detector_Tmax = 6000; // ticks 

int Ymin = 0, Ymax = 600; // cm
int Zmin = 0, Zmax = 700; // cm

/* frame parameters */
double frame_Zpositions[7] = {0,120,230,345,465,575,695};

/* tpc number (LArSoft numbering) */
int tpc_left_front = 2, tpc_left_mid = 6, tpc_left_back = 10;
int tpc_right_front = 1, tpc_right_mid = 5, tpc_right_back = 9;

/********** constant variables end ***********/

/********* struct definition *********/
struct Hit {
    float peakT;
    float deltaT_local;
    float deltaT_total;
    float y;
    float z;
    float x_calculated;
    int tpc;

    /*** constructors ***/
    Hit() : peakT( 0.0f ), deltaT_local( 0.0f ),
            y( 0.0f ), z( 0.0f ), x_calculated( 0.0f ),
            tpc( 0 ) {}
    Hit( float peakTIn, float deltaT_localIn, float deltaT_totalIn, float yIn, float zIn, float x_calcIn, int tpcIn) :
        peakT( peakTIn ), deltaT_total( deltaT_totalIn ), deltaT_local( deltaT_localIn ), 
            y( yIn ), z( zIn ), x_calculated( x_calcIn ), tpc( tpcIn ) {}

    /*** member functions ***/
    void print(int verbose);
};

// print 
void Hit::print(int verbose){
    if (verbose == 0)
        cout<<"x_calc = "<<x_calculated<<", y = "<<y<<", z = "<<z<<endl;

    if (verbose == 1)
        cout<<"x_calc = "<<x_calculated<<", y = "<<y<<", z = "<<z
        <<", deltaT_total = "<<deltaT_total<<endl;

    if (verbose == 5)
        cout<<"peakT = "<<peakT<<", deltaT_total = "<<deltaT_total
        <<", deltaT_local = "<<deltaT_local
        <<", x_calc = "<<x_calculated
        <<", y = "<<y<<", z = "<<z<<", tpc = "<<tpc<<endl;
}

struct Track { 
    vector<Hit> hits; 

    /*** constructor ***/
    Track(vector<Hit> hs): hits(hs) {}

    /*** member functions ***/
    // only to be used after the hits are sorted based on peakT
    float Tmin() {return hits[0].peakT;};
    int size() {return hits.size();};
    int deltaT_total() {return hits[hits.size()-1].peakT - hits[0].peakT;};
    Hit cathode_hit() {return hits[hits.size()-1];};
    void sort_by_T();
    void set_deltaT_local(float Tmin);
    void print(int verbose);
};

void Track::sort_by_T(){
    sort(hits.begin(), hits.end(), [&](const auto& hit1, const auto& hit2)
    {
        return hit1.peakT < hit2.peakT;
    });
}

void Track::set_deltaT_local(float Tmin){
    for (int i=0; i<hits.size(); i++){
        hits[i].deltaT_local = hits[i].peakT - Tmin;       
    }   
}

void Track::print(int verbose){
    for (int i=0; i<hits.size(); i++){
        hits[i].Hit::print(verbose);
    }
}

/********* struct definition end *********/

/********* helper functions *********/
// need to be tested
/* given a number n and a list of histogram parameters,
 *  * create 2*n such histograms (positive side & negative
 *   * side) (TH1F)*/
void create_n_hists(int n, TH1F *hists_pos[n], TH1F *hists_neg[n],
                    vector<char*> names, vector<char*> hist_titles, 
                    int x0, int x1, int nbins)
{
    for (int i=0; i<n; i++)
    {
        hists_pos[i] = new TH1F(Form("%s_pos_%d", names[i], i), 
                            Form("%s (Beam Left)", hist_titles[i]),
                            nbins, x0, x1);
        hists_neg[i] = new TH1F(Form("%s_neg_%d", names[i], i),
                            Form("%s (Beam Right)", hist_titles[i]),
                            nbins, x0, x1);
    }
}

// need to be tested
/* given a number n and a list of histogram parameters,
 *  * create 2*n such histograms (positive side & negative
 *   * side) (TH2F)*/
void create_n_hists(int n, TH2F *hists_pos[n], TH2F *hists_neg[n],
                    vector<char*> names, vector<char*> hist_titles, 
                    int x0, int x1, int nbinsX,
                    int y0, int y1, int nbinsY)
{
    for (int i=0; i<n; i++)
    {
        hists_pos[i] = new TH2F(Form("%s_pos_%d", names[i], i), 
                            Form("%s (Beam Left)", hist_titles[i]),
                            nbinsX, x0, x1, nbinsY, y0, y1);
        hists_neg[i] = new TH2F(Form("%s_neg_%d", names[i], i),
                            Form("%s (Beam Right)", hist_titles[i]),
                            nbinsX, x0, x1, nbinsY, y0, y1);
    }
}

void create_n_hists(int n, TH1F *hists_pos[n], TH1F *hists_neg[n],
                    char name[], char hist_title[], char x_unit[], char y_unit[],
                    int x0, int x1, int nbinsX)
{
    for (int i=0; i<n; i++)
    {
        hists_pos[i] = new TH1F(Form("%s_pos_%d", name, i),
                            Form("%s (Beam Left); %s; %s", hist_title,
                                 x_unit, y_unit), nbinsX, x0, x1);
        hists_neg[i] = new TH1F(Form("%s_neg_%d", name, i),
                            Form("%s (Beam Right); %s; %s", hist_title, 
                                 x_unit, y_unit), nbinsX, x0, x1);
    }
}

/*** probably need to be deleted later
vector<Hit> create_sorted_Hit_vector(vector<float> T, vector<float> y,
                              vector<float> z, vector<int> tpc)
{
    int hits_count = T.size();
    vector<Hit> hits;
    hits.reserve(hits_count);

    for (int i=0; i<hits_count; i++)
    {
        hits.push_back( Hit(T[i], y[i], z[i], 9999., tpc[i]) );
    }
    sort_hits_by_T(hits);

    return hits;   
}
*/

/* Given the Y ans Z bin numbers, return the corresponding 
 * distribution histogram number */
int YZhist_num(int i, int j, int nbinsY, int nbinsZ) {return i*nbinsY + j;}

/********* helper functions end *********/


void frame::Loop()
{
    string output = output_PATH + output_file_name;
    TFile *file = new TFile(output.c_str(), "recreate");

    /*** vector declaration ***/
    vector<Float_t> hitpeakT_buffer, hity_buffer, hitz_buffer;
    vector<Int_t> hittpc_buffer, hitwire_buffer;
    vector<Float_t> hitpeakT_buffer_neg, hity_buffer_neg, hitz_buffer_neg;
    vector<Int_t> hittpc_buffer_neg, hitwire_buffer_neg;

    vector<Track> selected_tracks, selected_tracks_neg; 

    /****** histogram definition ******/
    int nbinsY = Ymax / Ybinsize, nbinsZ = Zmax / Zbinsize;
    int deltaT_YZ_hists_num = nbinsY * nbinsZ;
    int nbinsT = detector_Tmax / Tbinsize; 
    TH1F *deltaT_YZ_hists[deltaT_YZ_hists_num], *deltaT_YZ_hists_neg[deltaT_YZ_hists_num];
    char deltaT_YZ_hist_name[] = "YZbin_deltaT_distribution_hist",
         deltaT_YZ_hist_title[] = "deltaT_total distribution (YZ plane bin)",
         deltaT_YZ_hist_xunit[] = "deltaT_total (ticks)",
         deltaT_YZ_hist_yunit[] = "number of hits";
    create_n_hists(deltaT_YZ_hists_num, deltaT_YZ_hists, deltaT_YZ_hists_neg,
                   deltaT_YZ_hist_name, deltaT_YZ_hist_title, 
                   deltaT_YZ_hist_xunit, deltaT_YZ_hist_yunit,
                   nbinsT, deltaT_min, detector_Tmax); 

    /****** histogram definition ******/

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    if (select_nentries != 0) nentries = select_nentries;   

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) // event loop
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if (jentry % 10000 == 0)
            std::cout << jentry << "/" << nentries << std::endl;

        if(!trkhitz_wire2->size()) continue;

        /* trk loop: ith track */
        for(size_t i=0;i<trkhitz_wire2->size();i++)
        {
            vector<Hit> hits, hits_neg; 

            /* hit loop: jth hit */
            for(size_t j=0;j<trkhitz_wire2->at(i).size();j++)
            {
                /****** select hits located on beam LEFT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_left_front || hit_tpc2->at(i)[j] == tpc_left_mid || hit_tpc2->at(i)[j] == tpc_left_back)
                {
                    hits.push_back( Hit(hit_peakT2->at(i)[j], 0, 0,
                                        trkhity2->at(i)[j],
                                        trkhitz_wire2->at(i)[j],
                                        9999.,hit_tpc2->at(i)[j]));  
                }
                /****** selet hits located on beam RIGHT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_right_front || hit_tpc2->at(i)[j] == tpc_right_mid || hit_tpc2->at(i)[j] == tpc_right_back)
                {
                    hits_neg.push_back( Hit(hit_peakT2->at(i)[j], 0, 0,
                                        trkhity2->at(i)[j],
                                        trkhitz_wire2->at(i)[j],
                                        -9999.,hit_tpc2->at(i)[j]));
                }
            } // end of hit loop

            if (hits.size() > hits_size_min)
            {
                // create Track object
                Track trk = Track(hits);
                trk.Track::sort_by_T();
                float Tmin = trk.Tmin();
                trk.Track::set_deltaT_local(Tmin);
                selected_tracks.push_back(trk);

           }  
        } // end of trk loop

        // fill deltaT distribution histograms
        for (int i=0; i<selected_tracks.size(); i++)
        {
            int ybin = selected_tracks[i].cathode_hit().y / Ybinsize;
            int zbin = selected_tracks[i].cathode_hit().z / Zbinsize;
            if (ybin < 0 || ybin > nbinsY || zbin < 0 || zbin > nbinsZ) continue;
            deltaT_YZ_hists[YZhist_num(ybin,zbin,nbinsY,nbinsZ)]->Fill(selected_tracks[i].deltaT_total());  

        }
 

    } // end of nentries loop

    for (int i=0; i<deltaT_YZ_hists_num; i++) {
        deltaT_YZ_hists[i]->Write();
    }


    cout<<"selected_trk_count_pos = "<<selected_tracks.size()<<endl
        <<"selected_trk_count_neg = "<<selected_tracks_neg.size()<<endl;
    
    file->Write();
    file->Close();
} 
