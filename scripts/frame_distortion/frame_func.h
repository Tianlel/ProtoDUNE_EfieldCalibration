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

