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
string output_file_name = "ALLEVENTS_deltaT_with_contraction_corr.root";

/* optional variables */
Long64_t select_nentries = 0; // 0 -> use all nentries
int set_jentry = 0;
int print_debug_message = 0;

/* cut parameters */
int hits_size_min = 5; 
int CA_crossing_cut = 0;

int YZ_deltaT_getmed_lowerbound = 4580; // find peak after this threshold
int YZ_deltaT_entries_min = 1000; // minimum number of entries to perform getmedian function

// deltaT cut for cathode-anode crosser selection
int Tcut_mid[6][2] = { {4595, 4605}, {4595, 4607}, {4593, 4608},
                       {4587, 4609}, {4601, 4594}, {4595, 4602} };
int Tcut_margin = 10; // ticks

/* histogram parameters */
int Ybinsize = 20; // cm
int Zbinsize = 20; // cm 
int Tbinsize = 50; // ticks 

int YZ_deltaT_binsize = 2; // ticks
int YZ_deltaT_min = 4000; // ticks
int YZ_deltaT_max = 5000; // ticks

/* detector parameters */
int detector_Tmax = 6000; // ticks 

int Ymin = 0, Ymax = 600; // cm
int Zmin = -10, Zmax = 710; // cm

/* frame parameters */
vector<double> frame_Zpositions{0,120,230,345,465,575,695};

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
    int frame_tag; // 0: [0,120); 1: [120,230); ...; 5: [575,695)

    /*** constructors ***/
    Hit() : peakT( 0.0f ), deltaT_local( 0.0f ),
            y( 0.0f ), z( 0.0f ), x_calculated( 0.0f ),
            tpc( 0 ), frame_tag( 0 ) {}
    Hit( float peakTIn, float deltaT_localIn, float deltaT_totalIn, float yIn, float zIn, float x_calcIn, int tpcIn, int frame_tagIn) :
        peakT( peakTIn ), deltaT_total( deltaT_totalIn ), deltaT_local( deltaT_localIn ), 
            y( yIn ), z( zIn ), x_calculated( x_calcIn ), 
            tpc( tpcIn ), frame_tag( frame_tagIn ) {}

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
    int get_cathode_frame_tag() {return hits[hits.size()-1].frame_tag;};
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
void print(string message) { cout<<message<<endl;}

/*** Applying thermal contraction to Y and Z values (by Ajib) ***/
double Zp[6]={0.575,230.112,232.635,462.172,464.695,694.232};
double Zn[6]={0.560,230.097,232.620,462.157,464.68,694.217};

double zthermal(double z, int tpcno){
    if(tpcno==1||tpcno==5||tpcno==9)
        return z+((Zn[(tpcno-1)/2]+Zn[(tpcno-1)/2+1])/2.0-z)*2.7e-3;
    if(tpcno==2||tpcno==6||tpcno==10)
        return z+((Zp[(tpcno-2)/2]+Zp[(tpcno-2)/2+1])/2.0-z)*2.7e-3;
    return -1;
}
double ythermal(double y){
  return y+(606.93-y)*2.7e-3;
}
/******/

int get_frame_tag(int z)
{
    if (frame_Zpositions[0]<=z && z<frame_Zpositions[1]) return 0;
    if (frame_Zpositions[1]<=z && z<frame_Zpositions[2]) return 1;
    if (frame_Zpositions[2]<=z && z<frame_Zpositions[3]) return 2;
    if (frame_Zpositions[3]<=z && z<frame_Zpositions[4]) return 3;
    if (frame_Zpositions[4]<=z && z<frame_Zpositions[5]) return 4;
    if (frame_Zpositions[5]<=z && z<frame_Zpositions[6]) return 5;
    return -1;
}

// side = 0 -> pos; side = 1 -> neg
bool is_CA_crosser(Track trk, int frame_tag, int side)
{ 
    int deltaT = trk.deltaT_total();
    int cut_midpoint = Tcut_mid[frame_tag][side];
    if (deltaT >= cut_midpoint - Tcut_margin &&
        deltaT <= cut_midpoint + Tcut_margin ) return true;
    else return false;   
}

/* given a non-empty histogram, return the median value;
 * given empty histogram, return 0*/
double get_hist_med (TH1F *h)
{
    if (h->GetEntries() == 0) return 0;
    double x, q;
    q = 0.5;
    h->ComputeIntegral();
    h->GetQuantiles(1, &x, &q);
    return x;
}

// need to be tested
/* given a number n and a list of histogram parameters,
 *  * create 2*n such histograms (positive side & negative
 *   * side) (TH1F)*/
void create_n_hists(int n, TH1F *hists_pos[n], TH1F *hists_neg[n],
                    vector<const char*> names, vector<const char*> hist_titles, 
                    char x_unit[], char y_unit[],
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

/* Given the Z and Y bin numbers, return the corresponding 
 *  *  * distribution histogram number */
int YZhist_num(int i, int j, int nbinsY) {return i*nbinsY + j;}

/* Given histogram number, get ybin, zbin number, the lower limit of the bin range in y and z*/
void getyzbin(int n, int *ybin, int *zbin, int nbinsY, int *y_range_min, int *z_range_min)
{
    *ybin = n%nbinsY;
    *zbin = n/nbinsY;
    *y_range_min = (*ybin)*Ybinsize + Ymin;
    *z_range_min = (*zbin)*Zbinsize + Zmin; 
}

void create_n_hists_YZ(int nbinsY, int n, TH1F *hists_pos[n], TH1F *hists_neg[n],
                    char name[], char hist_title[], char x_unit[], char y_unit[],
                    int x0, int x1, int nbinsX)
{
    int ybin = 0, zbin = 0, y_range_min = 0, z_range_min = 0;
    for (int i=0; i<n; i++)
    {
        getyzbin(i, &ybin, &zbin, nbinsY, &y_range_min, &z_range_min);
        hists_pos[i] = new TH1F(Form("%s_pos_%d_%d", name, zbin, ybin),
                            Form("%s (Beam Left), Ybin %d (%d-%d), Zbin %d (%d-%d); %s; %s", hist_title, ybin, y_range_min, y_range_min+Ybinsize, zbin, z_range_min, z_range_min+Zbinsize,
                                 x_unit, y_unit), nbinsX, x0, x1);
        hists_neg[i] = new TH1F(Form("%s_neg_%d_%d", name, zbin,ybin),
                            Form("%s (Beam Right),  Ybin %d (%d-%d), Zbin %d (%d-%d); %s; %s", hist_title, ybin, y_range_min, y_range_min+Ybinsize, zbin, z_range_min, z_range_min+Zbinsize,
                                 x_unit, y_unit), nbinsX, x0, x1);
    }
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
    TH2F *deltaT_YZ_h2 = new TH2F("deltaT_YZ_h2","deltaT vs YZ bin (beam left); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);
    TH2F *deltaT_YZ_h2_neg = new TH2F("deltaT_YZ_h2_neg","deltaT vs YZ bin (beam right); Z (cm); Y(cm); deltaT (ticks)", nbinsZ, Zmin, Zmax, nbinsY, Ymin, Ymax);


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
                                        ythermal(trkhity2->at(i)[j]),
                                        zthermal(trkhitz_wire2->at(i)[j], hit_tpc2->at(i)[j]),
                                        9999.,hit_tpc2->at(i)[j],
                                        frame_tag));  
                }
                /****** select hits located on beam RIGHT based on tpc number ******/
                if (hit_tpc2->at(i)[j] == tpc_right_front || hit_tpc2->at(i)[j] == tpc_right_mid || hit_tpc2->at(i)[j] == tpc_right_back)
                {
                    int frame_tag = get_frame_tag(trkhitz_wire2->at(i)[j]);
                    hits_neg.push_back( Hit(hit_peakT2->at(i)[j], 0, 0,
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
                trk.sort_by_T();
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
                }
            }

            if (hits_neg.size() > hits_size_min)
            {
                if (print_debug_message) print("create track object (neg)");
                // create Track object
                Track trk = Track(hits_neg);
                trk.sort_by_T();
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
            float bin_deltaT_med = (float) get_hist_med(deltaT_YZ_hists[YZhist_num(i,j,nbinsY)]);
            if (bin_deltaT_med == 0) continue;
            deltaT_YZ_h2->SetBinContent(i+1,j+1,bin_deltaT_med);

            float bin_deltaT_med_neg = (float) get_hist_med(deltaT_YZ_hists_neg[YZhist_num(i,j,nbinsY)]);
            if (bin_deltaT_med_neg == 0) continue;
            deltaT_YZ_h2_neg->SetBinContent(i+1,j+1,bin_deltaT_med_neg);
        }
    }
       
  /* for (int i=0; i<deltaT_YZ_hists_num; i++)
    {
        deltaT_YZ_hists[i]->Write(); deltaT_YZ_hists_neg[i]->Write();
    }

    deltaT_YZ_h2->Write(); deltaT_YZ_h2_neg->Write(); */

/*    for (int i=0; i<deltaT_zframe_hists_num; i++) {
        deltaT_zframe_hists[i]->Write(); deltaT_zframe_hists_neg[i]->Write();
    }
*/

    cout<<"selected_trk_count_pos = "<<selected_tracks.size()<<endl
        <<"selected_trk_count_neg = "<<selected_tracks_neg.size()<<endl;
    
    file->Write();
    file->Close();
} 
