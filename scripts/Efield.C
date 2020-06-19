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

/* constant variables */
#define X_BEGIN 0
#define X_END 3586
float Distance = X_END - X_BEGIN;

#define OUTPUT_FILE_NAME(str) ("%s_plots.root", str)
#define OUTPUT_FILE_TREE "t_err_vec"

using namespace std;

float X_BEGIN = 0;
float X_END = 3586;
float Distance = X_END - X_BEGIN;

void SetGraphStyle()
{
  TStyle *st = new TStyle("st1","my style");
  st->SetPadGridX(1);
  st->SetPadGridX(1);
  st->cd(); 
  //gStyle->SetOptStat(0);
}

//main section of the code for drift velocity
void Efield::Loop()
{
  TFile *file = new TFile(OUTPUT_FILE_NAME, "recreate");
  TTree *t = new TTree("t_err_vec","Tree with vectors");
  std::vector<float> e_rms_pos, e_rms_neg, e_se_pos, e_se_neg, e_gaus_pos, e_gaus_neg;
  std::vector<float> e_rms_posx, e_rms_negx, e_se_posx, e_se_negx, e_rms_neg_input_x, e_se_neg_input_x, e_rms_pos_input_x, e_se_pos_input_x;
  t->Branch("E_rms_pos",&e_rms_pos); t->Branch("E_rms_neg",&e_rms_neg);
  t->Branch("E_se_pos",&e_se_pos); t->Branch("E_se_neg",&e_se_neg);
  t->Branch("E_gaus_pos",&e_gaus_pos); t->Branch("E_gaus_neg",&e_gaus_neg);
  t->Branch("E_rms_posx",&e_rms_pos); t->Branch("E_rms_negx",&e_rms_negx);
  t->Branch("e_rms_neg_input_x",&e_rms_neg_input_x); t->Branch("e_se_neg_input_x",&e_se_neg_input_x);
  t->Branch("e_rms_pos_input_x",&e_rms_pos_input_x); t->Branch("e_se_pos_input_x",&e_se_pos_input_x);
  e_rms_pos.clear(); e_rms_neg.clear(); e_se_pos.clear(); e_se_neg.clear(); e_gaus_pos.clear(); e_gaus_neg.clear(); 
 e_se_neg_input_x.clear(); e_rms_pos_input_x.clear(); e_se_pos_input_x.clear();e_rms_neg_input_x.clear(); 
  // Set grid lines for all plots
  TStyle *st = new TStyle("st1","my style");
  st->SetPadGridX(1);
  st->SetPadGridY(1);
  st->cd();
  //gStyle->SetOptStat(0);
  
  float endbinsize = 20;
  int no_bins = 27;
  int no_bins_z = 140;
  int n_bins_x = 20;
  double xbin_size = X_END/n_bins_x; 
  int trk_count = 0;
  int trk_count_neg = 0;
  int z_angle_cnt = 0, z_angle_cnt_neg = 0;
  int num_hist = 0, num_hist_neg = 0; 
  float Tmax_neg = 4578; // Tmax for x<0
  float Tmin_neg = 4565; // Tmin for x<0
  float Tmax_pos = 4578;//4603;
  float Tmin_pos = 4565;//4583;
  double Temp = 87.67; //K
  float miny = 200;
  float maxy = 400;
  float minz = 240;
  float maxz = 450;

  float binsize = Tmax_pos / no_bins;
  float binsize_neg = Tmax_neg / no_bins;

  TSpline3 *sp = create_ve_sp(Temp);

  // general vectors for all calculation methods
  vector<Float_t> hitpeakT_buffer, hity_buffer, hitz_buffer;
  vector<Int_t> hittpc_buffer, hitwire_buffer;
  vector<Float_t> hitpeakT_buffer_neg, hity_buffer_neg, hitz_buffer_neg;
  vector<Int_t> hittpc_buffer_neg, hitwire_buffer_neg;

  // vectors for local drift velocity from fitting x vs T for each time bin (method 4)
  vector<Float_t> T_per_bin_per_trk[no_bins], x_per_bin_per_trk[no_bins], T_per_bin_per_trk_neg[no_bins], x_per_bin_per_trk_neg[no_bins];
  vector<Float_t> x_per_bin[no_bins], x_per_bin_neg[no_bins];
  vector<Float_t> T_vs_x[n_bins_x], T_vs_x_neg[n_bins_x];
  vector<vector<double>> drift_velocity_fit, drift_velocity_fit_neg;
  drift_velocity_fit.resize(no_bins); drift_velocity_fit_neg.resize(no_bins);  
  // for ploting efield vs x
  vector<Float_t> T_per_bin_per_trkx[n_bins_x], x_per_bin_per_trkx[n_bins_x], T_per_bin_per_trk_negx[n_bins_x], x_per_bin_per_trk_negx[n_bins_x];
  vector<vector<double>> drift_velocity_fitx, drift_velocity_fit_negx;
  drift_velocity_fitx.resize(n_bins_x); drift_velocity_fit_negx.resize(n_bins_x);

  // vectors for deviation in z from the assumed straight track (the YZ projection)
  vector<float> dz, dz_neg, dz_max, dz_max_neg;

  // vectors for hit removal (overlapping z)
  vector<float> z_vals, deltaz_deltaT, z_vals_neg, deltaz_deltaT_neg;

  // vectors for larsoft output
  vector<vector<double>> efieldx_pos, efieldx_neg;
  efieldx_neg.resize(no_bins), efieldx_pos.resize(no_bins);
  vector<vector<double>> efieldx_posx, efieldx_negx;
  efieldx_negx.resize(no_bins), efieldx_posx.resize(no_bins);
  
///////////////////////////defining all the histograms and graphs///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *deltaT_hist = new TH1F("deltaT_hist", "max peak_Time-min peakT in ticks-beam left;deltaT in ticks;no of tracks", 5000, 0, 5000);
  TH1F *deltaT_hist_neg = new TH1F("deltaT_hist_neg", "max peak_Time-min peakT in ticks for beam right;deltaT in ticks;no of tracks", 5000, 0, 5000);

  // deviation from straight line
  TH1D *dz_distribution = new TH1D("dz_distribution", "abs(z-z')_max distribution; d(mm); count", no_bins, 0, 200);
  TH1D *dz_distribution_neg = new TH1D("dz_distribution_neg", "abs(z-z')_max distribution; d(mm); count", no_bins, 0, 200);
  // difference between start angle and end angle
  TH1D *start_end_angle_diff_YZ_cos = new TH1D("start_end_angle_diff_YZ_cos", "angle difference between start and end of trk (YZ projection);cos #theta;count", 50, 0.90, 1);
  TH1D *start_end_angle_diff_YZ_cos_neg = new TH1D("start_end_angle_diff_YZ_cos_neg", "angle difference between start and end of trk (YZ projection) (x<0);cos #theta;count", 50, 0.90, 1);
  TH1D *start_end_angle_diff_YZ = new TH1D("start_end_angle_diff_YZ", "angle difference between start and end of trk (YZ projection); #theta(#circ);count", 50, 0, 10);
  TH1D *start_end_angle_diff_YZ_neg = new TH1D("start_end_angle_diff_YZ_neg", "angle difference between start and end of trk (YZ projection) (x<0); #theta(#circ);count", 50, 0, 10);
  // 
  TGraph *z_vs_y[10], *z_vs_y_neg[10], *z_vs_t[10], *z_vs_t_neg[10], *deltaz_deltaT_g[10];
  //TProfile *hitpeakt_vs_z[10], *hitpeakt_vs_z_neg[10];
  //TGraph *hitpeaktg_vs_z[10], *hitpeaktg_vs_z_neg[10];
  for (int i = 0; i < 10; i++)
  {
    //hitpeakt_vs_z[i] = new TProfile(Form("(x>0)track%d", i), Form("(x>0) track %d hitpeakt_vs_z", i), no_bins_z, 0, 700);
    //hitpeakt_vs_z_neg[i] = new TProfile(Form("(x<0)track%d", i), Form("(x<0) track %d hitpeakt_vs_z", i), no_bins_z, 0, 700);
    //hitpeaktg_vs_z[i] = new TGraph(); hitpeaktg_vs_z[i]->SetName(Form("g(x>0)track%d", i));
    //hitpeaktg_vs_z_neg[i] = new TGraph(); hitpeaktg_vs_z_neg[i]->SetName(Form("g(x<0)track%d", i));
    z_vs_y[i] = new TGraph(); z_vs_y[i]->SetName(Form("z_vs_y_%d", i));
    //z_vs_y_neg[i] = new TGraph(); z_vs_y_neg[i]->SetName(Form("z_vs_y_%d(x<0)", i)); 
    z_vs_t[i] = new TGraph(); z_vs_t[i]->SetName(Form("z_vs_t_%d", i));
    z_vs_t_neg[i] = new TGraph(); z_vs_t_neg[i]->SetName(Form("z_vs_t_neg%d", i));
    
    deltaz_deltaT_g[i] = new TGraph(); deltaz_deltaT_g[i]->SetName(Form("deltaz_deltaT_g_%d", i));
  }

  /* TH1D *r_squared_distribution_per_bin[no_bins], *r_squared_distribution_per_bin_neg[no_bins];
  for (int i=0; i<no_bins; i++)
  {
    r_squared_distribution_per_bin[i] = new TH1D(Form("r_squared_distribution_per_bin%d",i),Form("R squared distribution bin %d",i),50,0.9,1);
    r_squared_distribution_per_bin_neg[i] = new TH1D(Form("r_squared_distribution_per_bin%d(x<0)",i),Form("(x<0) R squared distribution bin %d",i),50,0.9,1);
  }
  */

  // method 4
  //TH1D *gausfit = new TH1D("gausfitv", "v from gaussian fit;time_bins;drift velocity in mm/us", no_bins, 0, Tmax_pos);
  float bin_edges_pos[no_bins+1]; bin_edges_pos[0]=0; bin_edges_pos[1]=endbinsize; bin_edges_pos[no_bins]=Tmax_pos; bin_edges_pos[no_bins-1]=Tmax_pos-endbinsize;
  float middle_binsize = (Tmax_pos-2*endbinsize)/(no_bins-2);
  for (int i=0; i<no_bins-3; i++)
  {
    bin_edges_pos[i+2]=bin_edges_pos[i+1]+middle_binsize;
    cout<<bin_edges_pos[i+2]<<endl;
  }
  float bin_edges_neg[no_bins+1]; bin_edges_neg[0]=0; bin_edges_neg[1]=endbinsize; bin_edges_neg[no_bins]=Tmax_neg; bin_edges_neg[no_bins-1]=Tmax_neg-endbinsize;
  float middle_binsize_neg = (Tmax_neg-2*endbinsize)/(no_bins-2);
  for (int i=0; i<no_bins-3; i++)
  {
    bin_edges_neg[i+2]=bin_edges_neg[i+1]+middle_binsize_neg;
    cout<<bin_edges_neg[i+2]<<endl;
  }

  TH1F *driftvel_med_fit = new TH1F("driftvel_med_fit", "Median drift velocity beam left (fit per trk);time_bins;drift velocity in mm/us", no_bins, bin_edges_pos);
  TH1F *driftvel_med_fit_neg = new TH1F("driftvel_med_fit_neg", "Median drift velocity for beam right (fit per trk);time_bins;drift velocity in mm/us", no_bins, bin_edges_neg);
  TH1F *efield_med_pos_fit = new TH1F("efield_med_pos_fit", "Electric Field beam left;time_bins;Electric field in kV/cm", no_bins, bin_edges_pos);
  TH1F *efield_med_neg_fit = new TH1F("efield_med_neg_fit", "Electric Field beam right;time_bins;Electric field in kV/cm", no_bins, bin_edges_neg);

  TH1F *efield_med_pos_fitx = new TH1F("xefield_med_pos_fit", "Electric Field beam left;x (mm);Electric field in kV/cm", n_bins_x,0, X_END);
  TH1F *efield_med_neg_fitx = new TH1F("xefield_med_neg_fit", "Electric Field beam right;x (mm);Electric field in kV/cm", n_bins_x, 0,X_END);

  TH1F *T_vs_xh = new TH1F("T_vs_x", "beam left;x (mm); time_bins", n_bins_x, 0, X_END);
  TH1F *T_vs_x_negh = new TH1F("T_vs_x_neg", "beam right;x (mm); time_bins", n_bins_x, 0, X_END);
  TH1D *efield_fit_trunc_mean = new TH1D("efield_fit_trunc_mean","Electric Field beam left (truncated mean);time_bins;Electric field in kV/cm", no_bins, bin_edges_pos);
  TH1D *efield_fit_trunc_mean_neg = new TH1D("efield_fit_trunc_mean_neg","Electric Field beam right (truncated mean);time_bins;Electric field in kV/cm", no_bins, bin_edges_neg);
  TH1D *x_trunc_mean = new TH1D("x_trunc_mean","measured x vs deltaT beam left (truncated mean);time_bins;x (mm)", no_bins, bin_edges_pos);
  TH1D *x_trunc_mean_neg = new TH1D("x_trunc_mean_neg","measured x vs deltaT beam right (truncated mean);time_bins;x (mm)", no_bins, bin_edges_neg);
  //TGraph *fitperbin[no_bins], *fitperbin_neg[no_bins];  
  TH1D *v_per_bin[no_bins], *v_per_bin_neg[no_bins];
  for (int i = 0; i < no_bins; i++)
  {
    v_per_bin[i] = new TH1D(Form("v_per_bin%d(fit_per_trk)", i), Form("v distribution bin %d (fit_per_trk)", i), 100, 1, 3);
    v_per_bin_neg[i] = new TH1D(Form("v_per_bin_neg%d(fit_per_trk)", i), Form("v distribution bin %d (x<0) (fit_per_trk)", i), 100, 1, 3);
  } 

  // larsoft efield
  TH1F *efieldx_pos_plot=new TH1F("efieldx_pos_plot","LArSoft input Electric field beam left;time in ticks;Electric Field",no_bins,bin_edges_pos);
  TH1F *efieldx_neg_plot=new TH1F("efieldx_neg_plot","LArSoft input Electric field beam right;time in ticks;Electric Field",no_bins,bin_edges_neg);

  TH1F *efieldx_pos_plotx=new TH1F("efieldx_pos_plotx","LArSoft input Electric field beam left;x (mm);Electric Field",n_bins_x,0,X_END);
  TH1F *efieldx_neg_plotx=new TH1F("efieldx_neg_plotx","LArSoft input Electric field beam right;x (mm);Electric Field",n_bins_x,0,X_END);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 
  if (fChain == 0)
    return;

  Long64_t nentries = fChain->GetEntries();

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
    // looping through all the tracks
    for (size_t i = 0; i < hit_peakT2->size(); i++)
    {
      if (!hit_peakT2->at(i).size())
        continue;
      drift_velocity_fit.clear();
      drift_velocity_fitx.clear();
      hitpeakT_buffer.clear();
      hittpc_buffer.clear();
      hitwire_buffer.clear();
      hitz_buffer.clear();
      hity_buffer.clear();
      hitpeakT_buffer_neg.clear();
      hittpc_buffer_neg.clear();
      hitwire_buffer_neg.clear();
      hitz_buffer_neg.clear();
      hity_buffer_neg.clear();
      for (int k = 0; k < no_bins; k++)
      {
        T_per_bin_per_trk[k].clear(); T_per_bin_per_trk_neg[k].clear();
        x_per_bin_per_trk[k].clear(); x_per_bin_per_trk_neg[k].clear();
      }
      for (int k = 0; k < n_bins_x; k++)
      {
        T_per_bin_per_trkx[k].clear(); T_per_bin_per_trk_negx[k].clear();
        x_per_bin_per_trkx[k].clear(); x_per_bin_per_trk_negx[k].clear();
      }

      dz.clear(); dz_neg.clear();
      z_vals.clear(); z_vals_neg.clear();
      deltaz_deltaT.clear(); deltaz_deltaT_neg.clear();

      /*****************filling the buffer variables beam left and beam right*********************************/
      // looping through all hits (ith track, jth hit)
      for (size_t j = 0; j < hit_peakT2->at(i).size(); j++)
      {
        if (hit_tpc2->at(i)[j] == 2 || hit_tpc2->at(i)[j] == 6 || hit_tpc2->at(i)[j] == 10)
        {
          //remove hits with coinciding z
          vector<float>::iterator it;
          it = find(z_vals.begin(), z_vals.end(), trkhitz_wire2->at(i)[j]);

          if (j==0 || it == z_vals.end()){
          z_vals.push_back(trkhitz_wire2->at(i)[j]);
          hitpeakT_buffer.push_back(hit_peakT2->at(i)[j]);
          hittpc_buffer.push_back(hit_tpc2->at(i)[j]);
          hitwire_buffer.push_back(hit_wire2->at(i)[j]);
          hity_buffer.push_back(trkhity2->at(i)[j]);
          hitz_buffer.push_back(trkhitz_wire2->at(i)[j]);
          }
        } //if tpc loop pos
        if (hit_tpc2->at(i)[j] == 1 || hit_tpc2->at(i)[j] == 5 || hit_tpc2->at(i)[j] == 9)
        {
          //remove hits with coinciding z
          vector<float>::iterator it;
          it = find(z_vals_neg.begin(), z_vals_neg.end(), trkhitz_wire2->at(i)[j]);

          if (j==0 || it == z_vals_neg.end()){
          z_vals_neg.push_back(trkhitz_wire2->at(i)[j]);
          hitpeakT_buffer_neg.push_back(hit_peakT2->at(i)[j]);
          hittpc_buffer_neg.push_back(hit_tpc2->at(i)[j]);
          hitwire_buffer_neg.push_back(hit_wire2->at(i)[j]);
          hity_buffer_neg.push_back(trkhity2->at(i)[j]);
          hitz_buffer_neg.push_back(trkhitz_wire2->at(i)[j]);
          }
        } //if tpc loop neg
      }   //hits loop

      /*************Filling positive drift or beam left*****************************/
      size_t siz = hitpeakT_buffer.size(); // total # of hits for x>0 for the ith track
      if (siz > 45)
      {
        vector<vector<Float_t> *> tmp_vec_for_sort;
        tmp_vec_for_sort.push_back(&hitpeakT_buffer);
        tmp_vec_for_sort.push_back(&hity_buffer);
        tmp_vec_for_sort.push_back(&hitz_buffer);
        sort_vecs_based_on_vec1(tmp_vec_for_sort);

        // find median of delta_z/delta_T 
        for (int k=0; k<siz-1; k++)
        {
          if (hitpeakT_buffer[k+1]-hitpeakT_buffer[k]==0) deltaz_deltaT.push_back(100.0);
          else deltaz_deltaT.push_back(abs((hitz_buffer[k+1]-hitz_buffer[k])/(hitpeakT_buffer[k+1]-hitpeakT_buffer[k])));
        }
        double median=TMath::Median(deltaz_deltaT.size(),&deltaz_deltaT[0]);
    
        float min = hitpeakT_buffer[0];       //peak_time at anode
        float max = hitpeakT_buffer[siz - 1]; //peak_time at cathode
        if (max < min)
          exit(0);

        float z0 = hitz_buffer[0];
        float z1 = hitz_buffer[siz - 1];
        float y0 = hity_buffer[0];
        float y1 = hity_buffer[siz - 1];
        if (hity_buffer[0] > miny && hity_buffer[0] < maxy && hity_buffer[siz - 1] > miny && hity_buffer[siz - 1] < maxy && hitz_buffer[0] > minz && hitz_buffer[0] < maxz && hitz_buffer[siz - 1] > minz && hitz_buffer[siz - 1] < maxz)
        deltaT_hist->Fill(max - min);
       
        // cuts
        int proj_devi_cut = 1;//YZ_projection_deviation_in_z(&hity_buffer, &hitz_buffer, &dz, 5, dz_distribution);
        int z_order_cut = 1; //z_order(&hitz_buffer); // make sure z is ordered
        int start_end_angle_diff_cut = angle_diff_start_end(trkstartcosxyz->at(i), trkendcosxyz->at(i), 0.2, start_end_angle_diff_YZ); // angle difference <0.2 degree bettween trk start and trk end
        int cut = z_order_cut*proj_devi_cut*start_end_angle_diff_cut;
        ////////////////////////////////////////////////////////
        // select anode-cathode crossing tracks (with more that 45 hits for x>0) within the minimum SCE effect volume
        if (cut == 1 && abs(z0-z1)>100 && (max - min) >= Tmin_pos && siz > 45 && (max - min) <= Tmax_pos && hity_buffer[0] > miny && hity_buffer[0] < maxy && hity_buffer[siz - 1] > miny && hity_buffer[siz - 1] < maxy && hitz_buffer[0] > minz && hitz_buffer[0] < maxz && hitz_buffer[siz - 1] > minz && hitz_buffer[siz - 1] < maxz)
        {
          //dz_max.push_back(maxdz);
          if(z0>z1) z_angle_cnt++;
          int bin = 0;
          float result = 999;
          
          trk_count++;

          // looping through hits with x>0
          for (int int1 = 0; int1 < siz; int1++)
          {
            if (int1<siz-1) 
            {
              float dzdt = abs((hitz_buffer[int1+1]-hitz_buffer[int1])/(hitpeakT_buffer[int1+1]-hitpeakT_buffer[int1]));
              //if (!(dzdt<1.5*median && dzdt>0.5*median)) continue;
            }
            float time = hitpeakT_buffer[int1];
            float z = hitz_buffer[int1];
            float y = hity_buffer[int1]; //new y
            float x = (abs(z1 - z) / abs(z1 - z0)) * Distance + X_BEGIN;
          
 
            if(num_hist<10){
              z_vs_t[num_hist]->SetPoint(z_vs_t[num_hist]->GetN(), hitpeakT_buffer[int1], hitz_buffer[int1]); 
              z_vs_y[num_hist]->SetPoint(z_vs_y[num_hist]->GetN(), hity_buffer[int1], hitz_buffer[int1]);
              //if(int1==0) myfile << "event " << event << " trkID " << TrkID->at(i) << " subrun " << subrun << "tpc " << hittpc_buffer[0] << "\n" << " trkstartx " << trkstartx->at(i) << " trkstarty " << trkstarty->at(i) << " trkstartz " << trkstartz->at(i) << "\n";
              if(int1<siz-1) deltaz_deltaT_g[num_hist]->SetPoint(deltaz_deltaT_g[num_hist]->GetN(), hitpeakT_buffer[int1], abs((hitz_buffer[int1+1]-hitz_buffer[int1])/(hitpeakT_buffer[int1+1]-hitpeakT_buffer[int1])));
            }

            float hitdeltatime = time - min;
            int bin1; 
            if (hitdeltatime < endbinsize) bin1 = 0;
            else if (hitdeltatime >= Tmax_pos - endbinsize) bin1 = no_bins-1;
            else bin1 = 1+(int)((hitdeltatime-endbinsize) / middle_binsize);
	    int x_bin = x/xbin_size;
            if(x_bin < n_bins_x) T_vs_x[x_bin].push_back(hitdeltatime);


            
            if (bin1 >= no_bins) {cout << bin1 << endl; continue;}

            T_per_bin_per_trk[bin1].push_back(hitpeakT_buffer[int1] - min);
            x_per_bin_per_trk[bin1].push_back(x);
            
            T_per_bin_per_trkx[x_bin].push_back(hitpeakT_buffer[int1] - min);
            x_per_bin_per_trkx[x_bin].push_back(x);

    
            // not relevent for E vs x
            x_per_bin[bin1].push_back(x);
             
            efieldx_pos[bin1].push_back(Ef_ex(0.1*x,y,z));
            efieldx_posx[x_bin].push_back(Ef_ex(0.1*x,y,z));

          } // hits loop

          for (int k = 0; k < no_bins; k++)
          {
            if (T_per_bin_per_trk[k].size() < 3)
              continue;
            vector<vector<Float_t> *> tmp_vec_for_sort;
            tmp_vec_for_sort.push_back(&T_per_bin_per_trk[k]);
            tmp_vec_for_sort.push_back(&x_per_bin_per_trk[k]);
            sort_vecs_based_on_vec1(tmp_vec_for_sort);

            TGraph *gfitbin = new TGraph(T_per_bin_per_trk[k].size(), &T_per_bin_per_trk[k][0], &x_per_bin_per_trk[k][0]);
            TF1 *f = new TF1("linearfit0", "pol1", 0, 4500); //linear fit
            gfitbin->Fit("linearfit0","Q");
	    if (k==10) {gfitbin->SetTitle("beam left time bin #10; deltaT (ticks); time-independent x (mm)");gfitbin->Write();} 
            float v_bin = -2 * (f->GetParameter(1));
            if (v_bin > 10 || v_bin < 0)
              continue;
            v_per_bin[k]->Fill(v_bin);
            drift_velocity_fit[k].push_back(v_bin);
          
            /* 
            // calculate R squared (loop thru all hits)
            float SEline = 0, SEx = 0;
	          float x_mean = TMath::Mean(x_per_bin_per_trk[k].begin(),x_per_bin_per_trk[k].end());
            for (int j=0; j<T_per_bin_per_trk[k].size(); j++)
            {
              float fitval = f->Eval(T_per_bin_per_trk[k][j]);
              SEline += pow((fitval-x_per_bin_per_trk[k][j]),2);
              SEx += pow((x_per_bin_per_trk[k][j]-x_mean),2);
            }
              //if (1-(SEline/SEx)<0.998) drift_velocity_fit[k].erase(drift_velocity_fit[k].end()-1);
              //r_squared_distribution_per_bin[k]->Fill(1-(SEline/SEx));
	*/  
          }

          // plot E vs x
	  for (int k = 0; k < n_bins_x; k++)
          {
            if (T_per_bin_per_trkx[k].size() < 3)
              continue;
            vector<vector<Float_t> *> tmp_vec_for_sort;
            tmp_vec_for_sort.push_back(&T_per_bin_per_trkx[k]);
            tmp_vec_for_sort.push_back(&x_per_bin_per_trkx[k]);
            sort_vecs_based_on_vec1(tmp_vec_for_sort);

            TGraph *gfitbin = new TGraph(T_per_bin_per_trkx[k].size(), &T_per_bin_per_trkx[k][0], &x_per_bin_per_trkx[k][0]);
            TF1 *f = new TF1("linearfit0", "pol1", 0, 4000); //linear fit
            gfitbin->Fit("linearfit0","Q");
            float v_bin = -2 * (f->GetParameter(1));
            if (v_bin > 10 || v_bin < 0) continue;
            drift_velocity_fitx[k].push_back(v_bin);
	  }



          if (num_hist < 10)
          {
            //hitpeakt_vs_z[num_hist]->Write();
            //hitpeaktg_vs_z[num_hist]->Write();
            z_vs_y[num_hist]->Write();
            z_vs_t[num_hist]->Write();
            deltaz_deltaT_g[num_hist]->Write();
          } 
          num_hist++;
        } // track selection loop
      }   // size not equal to 0


      /*********************Filling beam right or negative drift histograms************************************/
      size_t siz_neg = hitpeakT_buffer_neg.size();
      if (siz_neg >45)
      {
	      vector<vector<Float_t>*> tmp_vec_for_sort_neg;
	      tmp_vec_for_sort_neg.push_back(&hitpeakT_buffer_neg);
	      tmp_vec_for_sort_neg.push_back(&hity_buffer_neg);
	      tmp_vec_for_sort_neg.push_back(&hitz_buffer_neg);
	      sort_vecs_based_on_vec1(tmp_vec_for_sort_neg);

        // find median of delta_z/delta_T 
        for (int k=0; k<siz_neg-1; k++)
        {
          if (hitpeakT_buffer_neg[k+1]-hitpeakT_buffer_neg[k]==0) deltaz_deltaT_neg.push_back(100.0);
          else deltaz_deltaT_neg.push_back(abs((hitz_buffer_neg[k+1]-hitz_buffer_neg[k])/(hitpeakT_buffer_neg[k+1]-hitpeakT_buffer_neg[k])));
        }
        double median_neg=TMath::Median(deltaz_deltaT_neg.size(),&deltaz_deltaT_neg[0]);

        float min_neg = hitpeakT_buffer_neg[0];           //peak_time at anode
        float max_neg = hitpeakT_buffer_neg[siz_neg - 1]; //peak_time at cathode
        if (min_neg>max_neg) exit(0);

        float z0_neg = hitz_buffer_neg[0];
        float z1_neg = hitz_buffer_neg[siz_neg - 1];
        float y0_neg = hity_buffer_neg[0];
        float y1_neg = hity_buffer_neg[siz_neg - 1];
        
        if (hity_buffer_neg[0] > miny && hity_buffer_neg[0] < maxy && hity_buffer_neg[siz_neg - 1] > miny && hity_buffer_neg[siz_neg - 1] < maxy && hitz_buffer_neg[0] > minz && hitz_buffer_neg[0] < maxz && hitz_buffer_neg[siz_neg - 1] > minz && hitz_buffer_neg[siz_neg - 1] < maxz)
         deltaT_hist_neg->Fill(max_neg - min_neg);

	// cuts
        int proj_devi_cut_neg = 1;//YZ_projection_deviation_in_z(&hity_buffer_neg, &hitz_buffer_neg, &dz_neg, 10, dz_distribution_neg);
        int z_order_cut_neg = 1; //z_order(&hitz_buffer_neg); // make sure z is ordered
        int start_end_angle_diff_cut_neg = angle_diff_start_end(trkstartcosxyz->at(i), trkendcosxyz->at(i), 0.2, start_end_angle_diff_YZ_neg); // angle difference <0.2 degree bettween trk start and trk end
        int cut_neg = z_order_cut_neg*proj_devi_cut_neg*start_end_angle_diff_cut_neg;
        //cout<<"total cut: "<<cut_neg<<" dz cut: "<< proj_devi_cut_neg<<" start_end_angle_diff_cut_neg: "<<start_end_angle_diff_cut_neg<<endl;
        ////////////////////////////////////////////////////////
        if ((cut_neg == 1 && abs(z0_neg-z1_neg)>=100 && (max_neg - min_neg) > Tmin_neg && siz_neg > 45 && (max_neg - min_neg) < Tmax_neg) && (hity_buffer_neg[0] > miny && hity_buffer_neg[0] < maxy && hity_buffer_neg[siz_neg - 1] > miny && hity_buffer_neg[siz_neg - 1] < maxy && hitz_buffer_neg[0] > minz && hitz_buffer_neg[0] < maxz && hitz_buffer_neg[siz_neg - 1] > minz && hitz_buffer_neg[siz_neg - 1] < maxz))
        {
          if(z0_neg>z1_neg) z_angle_cnt_neg++;
          float t0_neg = min_neg;
          int bin_neg = 0;
          float result1 = 999;

          trk_count_neg++; // trk_counting number of tracks

          for (int int2 = 0; int2 < siz_neg; int2++)
          {
            if (int2< siz_neg-1)
	    {
              float dzdt = abs((hitz_buffer_neg[int2+1]-hitz_buffer_neg[int2])/(hitpeakT_buffer_neg[int2+1]-hitpeakT_buffer_neg[int2]));
              //if (!(dzdt<1.5*median_neg && dzdt>0.5*median_neg)) continue;
            }

            float time_neg = hitpeakT_buffer_neg[int2];
            float z_neg = hitz_buffer_neg[int2];
            float y_neg = hity_buffer_neg[int2]; //new
	    float x_neg = (abs(z1_neg - z_neg) / abs(z1_neg - z0_neg)) * Distance + X_BEGIN;
    
            // Fill histogram
            if (num_hist_neg < 10)
            {
             //hitpeakt_vs_z_neg[num_hist_neg]->Fill(hitz_buffer_neg[int2], hitpeakT_buffer_neg[int2] - min_neg);
             //hitpeaktg_vs_z_neg[num_hist]->SetPoint(hitpeaktg_vs_z_neg[num_hist]->GetN(), hitz_buffer_neg[int2], hitpeakT_buffer_neg[int2] - min_neg);
             //z_vs_y_neg[num_hist_neg]->SetPoint(z_vs_y_neg[num_hist_neg]->GetN(),hity_buffer_neg[int2],hitz_buffer_neg[int2]);
             z_vs_t_neg[num_hist_neg]->SetPoint(z_vs_t_neg[num_hist_neg]->GetN(),time_neg,z_neg);
            }
	   

            float hitdeltatime_neg = time_neg - min_neg;
            int bin1_neg; 
            int x_bin_neg = x_neg/xbin_size;
	    if(x_bin_neg < n_bins_x) T_vs_x_neg[x_bin_neg].push_back(hitdeltatime_neg);
  

            if (hitdeltatime_neg < endbinsize) bin1_neg = 0;
            else if (hitdeltatime_neg >= Tmax_neg - endbinsize) bin1_neg = no_bins-1;
            else bin1_neg = 1+(int)((hitdeltatime_neg-endbinsize) / middle_binsize_neg);
            
            T_per_bin_per_trk_neg[bin1_neg].push_back(hitpeakT_buffer_neg[int2]-min_neg);
	    x_per_bin_per_trk_neg[bin1_neg].push_back(x_neg);

            T_per_bin_per_trk_negx[x_bin_neg].push_back(hitpeakT_buffer_neg[int2]-min_neg);
            x_per_bin_per_trk_negx[x_bin_neg].push_back(x_neg);  

	    x_per_bin_neg[bin1_neg].push_back(x_neg);

            efieldx_neg[bin1_neg].push_back(Ef_ex(-0.1*x_neg,y_neg,z_neg));    
            efieldx_negx[x_bin_neg].push_back(Ef_ex(-0.1*x_neg,y_neg,z_neg));
          } //int2 hits loop

          for(int k=0; k<no_bins; k++)
          {
            if (T_per_bin_per_trk_neg[k].size()<3) continue;
            vector<vector<Float_t>*> tmp_vec_for_sort;
            tmp_vec_for_sort.push_back(&T_per_bin_per_trk_neg[k]);
            tmp_vec_for_sort.push_back(&x_per_bin_per_trk_neg[k]);
            sort_vecs_based_on_vec1(tmp_vec_for_sort);

            TGraph *gfitbin = new TGraph(T_per_bin_per_trk_neg[k].size(), &T_per_bin_per_trk_neg[k][0], &x_per_bin_per_trk_neg[k][0]);
            TF1* f_neg = new TF1("linearfit0_neg","pol1",0,4500); //linear fit	
            gfitbin->Fit("linearfit0_neg","Q");
            float v_bin = -2*(f_neg->GetParameter(1));
            if (v_bin>10 || v_bin<0) continue;
            v_per_bin_neg[k]->Fill(v_bin);
            drift_velocity_fit_neg[k].push_back(v_bin);
         }

          // E vs x
	  for(int k=0; k<no_bins; k++)
          {
            if (T_per_bin_per_trk_negx[k].size()<3) continue;
            vector<vector<Float_t>*> tmp_vec_for_sort;
            tmp_vec_for_sort.push_back(&T_per_bin_per_trk_negx[k]);
            tmp_vec_for_sort.push_back(&x_per_bin_per_trk_negx[k]);
            sort_vecs_based_on_vec1(tmp_vec_for_sort);

            TGraph *gfitbin = new TGraph(T_per_bin_per_trk_negx[k].size(), &T_per_bin_per_trk_negx[k][0], &x_per_bin_per_trk_negx[k][0]);
            TF1* f_neg = new TF1("linearfit0_neg","pol1",0,4000); //linear fit	
            gfitbin->Fit("linearfit0_neg","Q");
            float v_bin = -2*(f_neg->GetParameter(1));
            if (v_bin>10 || v_bin<0) continue;
            drift_velocity_fit_negx[k].push_back(v_bin);
         }


          if (num_hist_neg < 10)
          {
            //hitpeakt_vs_z_neg[num_hist_neg]->Write();
            //hitpeaktg_vs_z_neg[num_hist_neg]->Write();
            //z_vs_y_neg[num_hist_neg]->Write();
            z_vs_t_neg[num_hist_neg]->Write();
          }
          num_hist_neg++;
        } //big loop
      } //size neg (neg selected trk loop)
    } //i (trk) loop
  } //jentry

  for (int i = 0; i < no_bins; i++)
  {
    sort(drift_velocity_fit[i].begin(), drift_velocity_fit[i].end());
    int vsize = drift_velocity_fit[i].size();
    float medianv = TMath::Median(drift_velocity_fit[i].size(), &drift_velocity_fit[i][0]); // from method 4
    driftvel_med_fit->SetBinContent(i + 1, medianv);
    if (i==0 || i==no_bins-1) medianv = 1.4;
    efield_med_pos_fit->SetBinContent(i + 1, sp->Eval(medianv));
    float trunc_rms = TMath::RMS(drift_velocity_fit[i].begin()+(int)(vsize*0.1), drift_velocity_fit[i].end()-(int)(vsize*0.1));
    efield_med_pos_fit->SetBinError(i+1, 0.31*trunc_rms/sqrt(0.8*vsize));
    e_rms_pos.push_back(0.31*trunc_rms);
    e_se_pos.push_back(0.31*trunc_rms/sqrt(0.8*vsize));

    //truncated mean
    float trunc_mean_v = TMath::Mean(drift_velocity_fit[i].begin()+(int)(vsize*0.1), drift_velocity_fit[i].end()-(int)(vsize*0.1));
    efield_fit_trunc_mean->SetBinContent(i+1, sp->Eval(trunc_mean_v));

    // beam right
    sort(drift_velocity_fit_neg[i].begin(), drift_velocity_fit_neg[i].end()); 
    vsize = drift_velocity_fit_neg[i].size();
    float medianv_neg = TMath::Median(drift_velocity_fit_neg[i].size(), &drift_velocity_fit_neg[i][0]);
    if (i==0 || i==no_bins-1) medianv_neg = 1.4;
    driftvel_med_fit_neg->SetBinContent(i + 1, medianv_neg);
    efield_med_neg_fit->SetBinContent(i + 1, sp->Eval(medianv_neg));
    float trunc_rms_neg = TMath::RMS(drift_velocity_fit_neg[i].begin()+(int)(vsize*0.1), drift_velocity_fit_neg[i].end()-(int)(vsize*0.1));
    efield_med_neg_fit->SetBinError(i+1, 0.31*trunc_rms_neg);
    e_rms_neg.push_back(0.31*trunc_rms_neg);
    e_se_neg.push_back(0.31*trunc_rms_neg/sqrt(0.8*vsize));
  //  cout<<0.31*trunc_rms_neg/sqrt(0.8*vsize)<<endl;

    //truncated mean
    float trunc_mean_v_neg = TMath::Mean(drift_velocity_fit_neg[i].begin()+(int)(vsize*0.1), drift_velocity_fit_neg[i].end()-(int)(vsize*0.1));
    efield_fit_trunc_mean_neg->SetBinContent(i+1,sp->Eval(trunc_mean_v_neg));
    //    if (i==no_bins-2) efield_fit_trunc_mean_neg->SetBinContent(i + 2,0.5);

    efieldx_pos_plot->SetBinContent(i+1,TMath::Median(efieldx_pos[i].size(),&efieldx_pos[i][0]));
    efieldx_neg_plot->SetBinContent(i+1,TMath::Median(efieldx_neg[i].size(),&efieldx_neg[i][0]));
    //std::cout<<"efield median bin"<<i<<" is "<<TMath::Median(efieldx_pos[i].size(),&efieldx_pos[i][0])<<endl;

    // x vs T
    sort(x_per_bin[i].begin(), x_per_bin[i].end());
    int x_size = x_per_bin[i].size();
    float trunc_mean_x = TMath::Mean(x_per_bin[i].begin()+(int)(x_size*0.1), x_per_bin[i].end()-(int)(x_size*0.1));
    float trunc_rms_x = TMath::RMS(x_per_bin[i].begin()+(int)(x_size*0.1), x_per_bin[i].end()-(int)(x_size*0.1));
    x_trunc_mean->SetBinContent(i+1, trunc_mean_x);
    
    sort(x_per_bin_neg[i].begin(), x_per_bin_neg[i].end());
    int x_size_neg = x_per_bin_neg[i].size();
    trunc_mean_x = TMath::Mean(x_per_bin_neg[i].begin()+(int)(x_size_neg*0.1), x_per_bin_neg[i].end()-(int)(x_size_neg*0.1));
    trunc_rms_x = TMath::RMS(x_per_bin_neg[i].begin()+(int)(x_size_neg*0.1), x_per_bin_neg[i].end()-(int)(x_size_neg*0.1));
    x_trunc_mean_neg->SetBinContent(i+1, trunc_mean_x);
    
}



for(int i =0; i<n_bins_x; i++)
{
    // T vs x
    sort(T_vs_x[i].begin(), T_vs_x[i].end());
    int T_size = T_vs_x[i].size();
    float trunc_mean_T = TMath::Mean(T_vs_x[i].begin()+(int)(T_size*0.1), T_vs_x[i].end()-(int)(T_size*0.1));
    float trunc_rms_T = TMath::RMS(T_vs_x[i].begin()+(int)(T_size*0.1), T_vs_x[i].end()-(int)(T_size*0.1));
    T_vs_xh->SetBinContent(i+1, trunc_mean_T);
    T_vs_xh->SetBinError(i+1, trunc_rms_T/sqrt(T_size*0.8));
    //cout<< "trunc_rms_T: "<<trunc_rms_T/sqrt(T_size*0.8)<<endl;

    sort(T_vs_x_neg[i].begin(), T_vs_x_neg[i].end());
    T_size = T_vs_x_neg[i].size();
    trunc_mean_T = TMath::Mean(T_vs_x_neg[i].begin()+(int)(T_size*0.1), T_vs_x_neg[i].end()-(int)(T_size*0.1));
    trunc_rms_T = TMath::RMS(T_vs_x_neg[i].begin()+(int)(T_size*0.1), T_vs_x_neg[i].end()-(int)(T_size*0.1));
    T_vs_x_negh->SetBinContent(i+1, trunc_mean_T);
    T_vs_x_negh->SetBinError(i+1, trunc_rms_T/sqrt(T_size*0.8));
    //cout<< "trunc_rms_T_neg: "<<trunc_rms_T/sqrt(T_size*0.8)<<endl;

    // beam left
    sort(drift_velocity_fitx[i].begin(), drift_velocity_fitx[i].end());
    int vsize = drift_velocity_fitx[i].size();
    float medianv = TMath::Median(drift_velocity_fitx[i].size(), &drift_velocity_fitx[i][0]); // from method 4
    if (i==0 || i==no_bins-1) medianv = 1.4;
    efield_med_pos_fitx->SetBinContent(i + 1, sp->Eval(medianv));
    float trunc_rms = TMath::RMS(drift_velocity_fitx[i].begin()+(int)(vsize*0.1), drift_velocity_fitx[i].end()-(int)(vsize*0.1));
    efield_med_pos_fitx->SetBinError(i+1, 0.31*trunc_rms);
    e_rms_posx.push_back(0.31*trunc_rms);
    e_se_posx.push_back(0.31*trunc_rms/sqrt(0.8*vsize));


    // beam right
    sort(drift_velocity_fit_negx[i].begin(), drift_velocity_fit_negx[i].end()); 
    vsize = drift_velocity_fit_negx[i].size();
    float medianv_neg = TMath::Median(drift_velocity_fit_negx[i].size(), &drift_velocity_fit_negx[i][0]);
    if (i==0 || i==no_bins-1) medianv_neg = 1.4;
    efield_med_neg_fitx->SetBinContent(i + 1, sp->Eval(medianv_neg));
    float trunc_rms_neg = TMath::RMS(drift_velocity_fit_negx[i].begin()+(int)(vsize*0.1), drift_velocity_fit_negx[i].end()-(int)(vsize*0.1));
    efield_med_neg_fitx->SetBinError(i+1, 0.31*trunc_rms_neg);
    e_rms_negx.push_back(0.31*trunc_rms_neg);
    e_se_negx.push_back(0.31*trunc_rms_neg/sqrt(0.8*vsize));

    float inputrms = TMath::RMS(efieldx_posx[i].size(),&efieldx_posx[i][0]);
    int insiz = efieldx_posx[i].size();
cout<<insiz<<endl;
    efieldx_pos_plotx->SetBinContent(i+1,TMath::Median(efieldx_posx[i].size(),&efieldx_posx[i][0]));
    efieldx_pos_plotx->SetBinError(i+1, inputrms/sqrt(0.8*insiz));
    e_rms_pos_input_x.push_back(inputrms);
    e_se_pos_input_x.push_back(inputrms/sqrt(0.8*insiz));   

    inputrms = TMath::RMS(efieldx_negx[i].size(),&efieldx_negx[i][0]);
    insiz = efieldx_negx[i].size();
    efieldx_neg_plotx->SetBinContent(i+1,TMath::Median(efieldx_negx[i].size(),&efieldx_negx[i][0]));
    efieldx_neg_plotx->SetBinError(i+1,inputrms/sqrt(0.8*insiz));
    e_rms_neg_input_x.push_back(inputrms);
    e_se_neg_input_x.push_back(inputrms/sqrt(0.8*insiz));

 }

  deltaT_hist->Write(); deltaT_hist_neg->Write();

  driftvel_med_fit->Write(); driftvel_med_fit_neg->Write();
  efield_med_pos_fit->Write(); efield_med_neg_fit->Write();

  start_end_angle_diff_YZ->Write();
  dz_distribution->Write(); dz_distribution_neg->Write();

  efieldx_pos_plot->Write(); efieldx_neg_plot->Write();

  efield_med_pos_fitx->Write(); efield_med_neg_fitx->Write();
  efieldx_pos_plotx->Write(); efieldx_neg_plotx->Write();

  efield_fit_trunc_mean->Write(); efield_fit_trunc_mean_neg->Write();
 
  x_trunc_mean->Write(); x_trunc_mean_neg->Write();

  T_vs_xh->Write(); T_vs_x_negh->Write();

  

  cout << "trk_count in positive drift " << trk_count << endl;
  cout << "trk_count in negative drift " << trk_count_neg << endl;
  cout << z_angle_cnt<< " trks have z0>z1 for x>0"<<endl;
  cout << z_angle_cnt_neg<< " trks have z0>z1 for x<0"<<endl;
  
  t->Fill();
  file->Write();
  file->Close();
} //Loop


