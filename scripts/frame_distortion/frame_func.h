#ifndef _FRAME_FUNC_H
#define _FRAME_FUNC_H

/********** constant variables ***********/

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

/* tpc number (LArSoft numbering) */
int tpc_left_front = 2, tpc_left_mid = 6, tpc_left_back = 10;
int tpc_right_front = 1, tpc_right_mid = 5, tpc_right_back = 9;

/********** constant variables end ***********/


/* frame parameters */
vector<double> frame_Zpositions{0,120,230,345,465,575,695};

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
    return z+((Zn[(tpcno-1)/2]+Zn[(tpcno-1)/2+1])/2.0-z)*2.7e-3+0.63*(5-tpcno)/4.0;
  if(tpcno==2||tpcno==6||tpcno==10) 
    return z+((Zp[(tpcno-2)/2]+Zp[(tpcno-2)/2+1])/2.0-z)*2.7e-3+0.63*(6-tpcno)/4.0;
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

/* Specific helper function for generating the 1D deltaT distributions histograms on the YZ plane */
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

void fill_deltaT_vec(float deltaT, int deltaT_margin, int zbin, int ybin, vector<float> *deltaT_meds, vector<vector<double>> &peak_vals)
{
    int num = YZhist_num(zbin, ybin, (Ymax-Ymin)/Ybinsize);
    if (deltaT_meds->at(num)-deltaT_margin <= deltaT && deltaT <= deltaT_meds->at(num)+deltaT_margin)
    {
        peak_vals[num].push_back(deltaT);
    }
}

void get_vec_med_err(vector<double> peak_vals, float *med, float *err_up, float *err_down)
{
    int size = peak_vals.size();
    if (size < 5) 
    {
        *med = 0;
        return;
    }
    double sigma = 0.5/sqrt(size+2);
    double p[3] = {0.5-sigma, 0.5, 0.5+sigma}, qout[3];
    double *a = &peak_vals[0]; 
    TMath::Quantiles(peak_vals.size(),3,&a[0],qout,p,0,0,8); 
    *med = qout[1];
    *err_up = qout[2]-qout[1];
    *err_down = qout[1]-qout[0];
}
/********* helper functions end *********/

#endif
