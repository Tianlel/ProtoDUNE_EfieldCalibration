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
string output_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/frame_distortion/";
string output_file_name = "test_plots.root";

/* optional variables */
Long64_t select_nentries = 0; // 0 -> use all nentries

/* frame parameters */
double frame_Zpositions[7] = {0,120,230,345,465,575,695};

/* tpc number (LArSoft numbering) */
int tpc_left_front = 2, tpc_left_mid = 6, tpc_left_back = 10;
int tpc_right_front = 1, tpc_right_mid = 5, tpc_right_back = 9;

/********** constant variables end ***********/

/********* struct definition *********/
struct Hit {
    float peakT;
    float y;
    float z;
    float x_calculated;
    int tpc;

    /*** constructors ***/
    Hit() : peakT( 0.0f ), y( 0.0f ), z( 0.0f ), x_calculated( 0.0f ),
            tpc( 0 ) {}
    Hit( float peakTIn, float yIn, float zIn, float x_calcIn, int tpcIn) :
        peakT( peakTIn ), y( yIn ), z( zIn ), x_calculated( x_calcIn ),
        tpc( tpcIn ){}
};

struct Track { vector<Hit> hits; };

    /*** constructors ***/
/********* struct definition end *********/

/********* helper functions *********/
void sort_hits_by_T(vector<Hit> &hits)
{
    sort(hits.begin(), hits.end(), [&](const auto& hit1, const auto& hit2)
    {
        return hit1.peakT < hit2.peakT;
    });
}

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

    /****** histogram definition ******/
    

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
            /* hit loop: jth hit */
            for(size_t j=0;j<trkhitz_wire2->at(i).size();j++)
            {
                /****** select hits located on beam LEFT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_left_front || hit_tpc2->at(i)[j] == tpc_left_mid || hit_tpc2->at(i)[j] == tpc_left_back)
                {
                    hitpeakT_buffer.push_back(hit_peakT2->at(i)[j]);
                    hittpc_buffer.push_back(hit_tpc2->at(i)[j]);
                    hity_buffer.push_back(trkhity2->at(i)[j]);
                    hitz_buffer.push_back(trkhitz_wire2->at(i)[j]);
                }
                /****** selet hits located on beam RIGHT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_right_front || hit_tpc2->at(i)[j] == tpc_right_mid || hit_tpc2->at(i)[j] == tpc_right_back)
                {
                    hitpeakT_buffer_neg.push_back(hit_peakT2->at(i)[j]);
                    hittpc_buffer_neg.push_back(hit_tpc2->at(i)[j]);
                    hity_buffer_neg.push_back(trkhity2->at(i)[j]);
                    hitz_buffer_neg.push_back(trkhitz_wire2->at(i)[j]);
                }
            } // end of hit loop
       
            vector<Hit> hits = create_sorted_Hit_vector(hitpeakT_buffer, hity_buffer, hitz_buffer, hittpc_buffer);

        } // end of trk loop
    } // end of nentries loop
}
