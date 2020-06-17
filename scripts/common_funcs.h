#ifndef _HELPERS_H
#define _HELPERS_H

/* common helper functions for Efield Calibration */

#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;

double E0=0.4867; // kV/cm

///////////////defining parametric relation between Efield and temperature//////////////////////
double temperature=87.67;//K
double tshift = -87.203+temperature;
double xFit = 0.0938163-0.0052563*tshift-0.0001470*tshift*tshift;
double uFit = 5.18406+0.01448*tshift-0.003497*tshift*tshift-0.000516*tshift*tshift*tshift;
double vd;
// Icarus Parameter Set, use as default
double  P1 = -0.04640; // K^-1
double  P2 = 0.01712;  // K^-1
double  P3 = 1.88125;   // (kV/cm)^-1
double  P4 =  0.99408;    // kV/cm
double  P5 =  0.01172;   // (kV/cm)^-P6
double  P6 =  4.20214;
double  T0 =  105.749;  // K
// Walkowiak Parameter Set
double  P1W = -0.01481; // K^-1
double  P2W = -0.0075;  // K^-1
double   P3W =  0.141;   // (kV/cm)^-1
double   P4W =  12.4;    // kV/cm
double   P5W =  1.627;   // (kV/cm)^-P6
double   P6W =  0.317;
double   T0W =  90.371;  // K
double driftvelo(double efield,double temperature){
	if (efield < xFit) vd=efield*uFit;
	else if (efield<0.619) { 
		vd = ((P1*(temperature-T0)+1)
				*(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
				+P2*(temperature-T0));
	}
	else if (efield<0.699) {
		vd = 12.5*(efield-0.619)*((P1W*(temperature-T0W)+1)
				*(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
				+P2W*(temperature-T0W))+
			12.5*(0.699-efield)*((P1*(temperature-T0)+1)
					*(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
					+P2*(temperature-T0));
	}
	else {
		vd = ((P1W*(temperature-T0W)+1)
				*(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
				+P2W*(temperature-T0W));     
	}

	return vd;

}

TSpline3 *create_ve_sp(double Temp)
{
    vector<double> driftv, Ef;
    // TSpline loop
    for (int i = 0; i <= 400; i++)
    {
    Ef.push_back(0.300 + i * 0.001);
    double Efield1 = 0.300 + i * 0.001;
    driftv.push_back(driftvelo(Efield1, Temp));
    }
    TSpline3 *sp = new TSpline3("Cubic Spline", &driftv[0], &Ef[0], Ef.size(), "b2e2", 0, 0);
    return sp;
}

TH1D *createhist_vt(string histname, double min, double max, int no_bins)
{
    string axis (";deltaT in ticks;Drift velocity in mm/ux");
    TH1D *h = new TH1D(histname.c_str(), (histname + axis).c_str(), no_bins, min, max);
    return h;
}

// getting distortion correction from Ajib
TFile *dist_map=new TFile("/dune/app/users/tianlel/protoDUNE/E_field/sim/ajib_sceon_distortionmap_feb16.root");
TH3F *Zcorrection=(TH3F*)dist_map->Get("Zdist_3D"); 

double zcorrect(double x_calc, double y_measured, double z_measured)
{ 
  double zcorrection=Zcorrection->Interpolate(x_calc, y_measured, z_measured);
  return zcorrection;
}

////getting the variable Efield using data driven maps
// from Ajib
TFile *ef=new TFile("/dune/app/users/tianlel/protoDUNE/E_field/sim/SCE_DataDriven_180kV_v3.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
float tot_Ef(float xval,float yval,float zval){
  if(xval>=0){
    float ex=E0+E0*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
    float ey=E0*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
    float ez=E0*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
if(xval<0){
    float ex=E0+E0*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
    float ey=E0*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
    float ez=E0*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
    return sqrt(ex*ex+ey*ey+ez*ez);
  }
 else return E0;
}

float Ef_ex(float xval,float yval,float zval)
{
  if(xval>=0){
    float ex=E0+E0*xpos->GetBinContent(xpos->FindBin(xval,yval,zval));
    return ex; 
  }
if(xval<0){
    float ex=E0+E0*xneg->GetBinContent(xneg->FindBin(xval,yval,zval));
    return ex;
  }
 else return E0;
}

TH3F *zpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Pos");
TH3F *zneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Neg");
TH3F *zpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
TH3F *zneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Neg");

float zoffset(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_fd->GetBinContent(zpos_fd->FindBin(xval,yval,zval));
  }
if(xval<0){
  return zneg_fd->GetBinContent(zneg_fd->FindBin(xval,yval,zval));
  }
  return -1000;
}

float zoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
  }
if(xval<0){
  return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  }
  return -1000;
}


template <typename A, typename B>
void zip(
                const std::vector<A> &a,
                const std::vector<B> &b,
                std::vector<std::pair<A,B>> &zipped)
{
        for(size_t i=0; i<a.size(); ++i)
        {
                zipped.push_back(std::make_pair(a[i], b[i]));
        }
}

template <typename A, typename B>
void unzip(
                const std::vector<std::pair<A, B>> &zipped,
                std::vector<A> &a,
                std::vector<B> &b)
{
        for(size_t i=0; i<a.size(); i++)
        {
                a[i] = zipped[i].first;
                b[i] = zipped[i].second;
        }
}

template <typename A, typename B, typename C>
void zip3(
                const std::vector<A> &a,
                const std::vector<B> &b,
                const std::vector<C> &c,
               std::vector<std::tuple<A,B,C>> &zipped)
{
        for(size_t i=0; i<a.size(); ++i)
        {
                zipped.push_back(make_tuple(a[i], b[i], c[i]));
        }
}

template <typename A, typename B, typename C>
void unzip3(
                const std::vector<std::tuple<A, B, C>> &zipped,
                std::vector<A> &a,
                std::vector<B> &b,
                std::vector<C> &c)
{
        for(size_t i=0; i<a.size(); i++)
        {
                a[i] = get<0>(zipped[i]);
                b[i] = get<1>(zipped[i]);
                c[i] = get<2>(zipped[i]);
        }
}

template <typename A, typename B, typename C, typename D>
void zip4(
                const std::vector<A> &a,
                const std::vector<B> &b,
                const std::vector<C> &c,
		const std::vector<D> &d,
               std::vector<std::tuple<A,B,C,D>> &zipped)
{
        for(size_t i=0; i<a.size(); ++i)
        {
                zipped.push_back(make_tuple(a[i], b[i], c[i], d[i]));
        }
}

template <typename A, typename B, typename C, typename D>
void unzip4(
                const std::vector<std::tuple<A, B, C, D>> &zipped,
                std::vector<A> &a,
                std::vector<B> &b,
                std::vector<C> &c,
		std::vector<D> &d)
{
        for(size_t i=0; i<a.size(); i++)
        {
                a[i] = get<0>(zipped[i]);
                b[i] = get<1>(zipped[i]);
                c[i] = get<2>(zipped[i]);
		d[i] = get<3>(zipped[i]);
	}
}
template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
void zip7(
                const std::vector<A> &a,
                const std::vector<B> &b,
                const std::vector<C> &c,
		const std::vector<D> &d,
		const std::vector<E> &e,
                const std::vector<F> &f,
		const std::vector<G> &g,
               std::vector<std::tuple<A,B,C,D,E,F,G>> &zipped)
{
        for(size_t i=0; i<a.size(); ++i)
        {
                zipped.push_back(make_tuple(a[i], b[i], c[i], d[i], e[i], f[i], g[i]));
        }
}

template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
void unzip7(
                const std::vector<std::tuple<A, B, C, D, E,F,G>> &zipped,
                std::vector<A> &a,
                std::vector<B> &b,
                std::vector<C> &c,
		std::vector<D> &d,
		std::vector<E> &e,
                std::vector<F> &f,
		std::vector<G> &g)
{
        for(size_t i=0; i<a.size(); i++)
        {
                a[i] = get<0>(zipped[i]);
                b[i] = get<1>(zipped[i]);
                c[i] = get<2>(zipped[i]);
		d[i] = get<3>(zipped[i]);
		e[i] = get<4>(zipped[i]);
		f[i] = get<5>(zipped[i]);
		g[i] = get<6>(zipped[i]);
	}
}
/* 
 * Sort n vectors based on the first vector (supposedly hitpeakT)
 */
template<typename A>
void sort_vecs_based_on_vec1(vector<vector<A>*> &vecs)
{
	int num_vec = vecs.size();
	if (num_vec==2)
	{
		vector<pair<A, A>> zipped;
		zip(*vecs[0], *vecs[1], zipped);
	
		sort(begin(zipped), end(zipped), [&](const auto& a, const auto& b)
        	{
            		return get<0>(a) < get<0>(b);
       		});
       		unzip(zipped, *vecs[0], *vecs[1]);
		return;
	}
	if (num_vec==3)
	{
		vector<tuple<A, A, A>> zipped;
                zip3(*vecs[0], *vecs[1], *vecs[2], zipped);
        
                sort(begin(zipped), end(zipped), [&](const auto& a, const auto& b)
                {
                        return get<0>(a) < get<0>(b);
                });
                unzip3(zipped, *vecs[0], *vecs[1], *vecs[2]);
                return;
	}
	if (num_vec==4)
        {       
                vector<tuple<A, A, A, A>> zipped;
                zip4(*vecs[0], *vecs[1], *vecs[2], *vecs[3], zipped);
                
                sort(begin(zipped), end(zipped), [&](const auto& a, const auto& b)
                {       
                        return get<0>(a) < get<0>(b);
                });
                unzip4(zipped, *vecs[0], *vecs[1], *vecs[2], *vecs[3]);
                return;
        }
        if (num_vec==7)
        {
                vector<tuple<A, A, A, A, A, A, A>> zipped;
                zip7(*vecs[0], *vecs[1], *vecs[2], *vecs[3], *vecs[4], *vecs[5], *vecs[6], zipped);

                sort(begin(zipped), end(zipped), [&](const auto& a, const auto& b)
                {
                        return get<0>(a) < get<0>(b);
                });
                unzip7(zipped, *vecs[0], *vecs[1], *vecs[2], *vecs[3], *vecs[4], *vecs[5], *vecs[6]);
                return;
        }


	else 
	{
		fprintf(stderr, "error: right now this func can only sort up to 3 vectors together\n");
		return;
	}
}
///////////////////cuts//////////////////

// returns 1 if the given track passes the cut 
int YZ_projection_deviation_in_z(vector<float> *y_buffer, vector<float> *z_buffer, vector<float> *dz, float cut_dz)
{
    float siz = y_buffer->size();
    
    float z0 = z_buffer->at(0);
    float z1 = z_buffer->at(siz - 1);
    float y0 = y_buffer->at(0);
    float y1 = y_buffer->at(siz - 1);

    for (int i = 0; i < siz; i++)
    {
        float z = z_buffer->at(i);
        float y = y_buffer->at(i); //new y
        float z_str = (z1 - z0) * (y - y0) / (y1 - y0) + z0;
        dz->push_back(abs(z_str - z));
    }
    float maxdz = TMath::MaxElement(dz->size(), &dz->at(0));
    if (maxdz<cut_dz) return 1; // cut passes
    else return 0; // cut fails
}

int YZ_projection_deviation_in_z(vector<float> *y_buffer, vector<float> *z_buffer, vector<float> *dz, float cut_dz, TH1D *dz_distribution)
{
    float siz = y_buffer->size();

    float z0 = z_buffer->at(0);
    float z1 = z_buffer->at(siz - 1);
    float y0 = y_buffer->at(0);
    float y1 = y_buffer->at(siz - 1);

    for (int i = 0; i < siz; i++)
    {
        float z = z_buffer->at(i);
        float y = y_buffer->at(i); //new y
        float z_str = (z1 - z0) * (y - y0) / (y1 - y0) + z0;
        dz->push_back(abs(z_str - z));
    }
    float maxdz = TMath::MaxElement(dz->size(), &dz->at(0));
    dz_distribution->Fill(maxdz);
    //cout<<maxdz<<endl;
    if (maxdz<cut_dz) return 1; // cut passes
    else return 0; // cut fails
}

// returns 1 if the given z vector is ordered 
int z_order(vector<float> *z_buffer)
{
    float siz = z_buffer->size();

    float z0 = z_buffer->at(0);
    float z1 = z_buffer->at(siz-1);

    int order=1;
    if(z0>z1)
    {
        for(int k=0;k<siz-1;k++)
        {
            if(z_buffer->at(k)<=z_buffer->at(k+1)) {order=0; 
            //  cout<<"wrong order"<<endl;
            //  for (int j=0; j<siz-1; j++)
            //  { cout<<"z: "<<z_buffer->at(k)<<endl;
            //    break;
            }
        }
    }
    if(z0<z1)
    {
        for(int k=0;k<siz-1;k++)
        {
            if(z_buffer->at(k)>=z_buffer->at(k+1)) order=0;
        }
    }    
    return order; // 0 fails; 1 passes
}

//
int angle_diff_start_end(vector<float> trkstartcosxyz, vector<float> trkendcosxyz, float cut_angle_diff)
{
    float cosy_start = trkstartcosxyz[1], cosz_start = trkstartcosxyz[2];
    float cosy_end = trkendcosxyz[1], cosz_end = trkendcosxyz[2];
    float cos_thetayz_diff = (cosy_start * cosy_end + cosz_start * cosz_end) / TMath::Sqrt((cosy_start * cosy_start + cosz_start * cosz_start) * (cosy_end * cosy_end + cosz_end * cosz_end));
    //start_end_angle_diff_YZ_cos->Fill(cos_thetayz_diff);
    //start_end_angle_diff_YZ->Fill(TMath::ACos(cos_thetayz_diff));  
    if (TMath::ACos(cos_thetayz_diff)<cut_angle_diff) return 1;
    else return 0;
}

int angle_diff_start_end(vector<float> trkstartcosxyz, vector<float> trkendcosxyz, float cut_angle_diff, TH1D* start_end_angle_diff_YZ)
{
    float cosy_start = trkstartcosxyz[1], cosz_start = trkstartcosxyz[2];
    float cosy_end = trkendcosxyz[1], cosz_end = trkendcosxyz[2];
    float cos_thetayz_diff = (cosy_start * cosy_end + cosz_start * cosz_end) / TMath::Sqrt((cosy_start * cosy_start + cosz_start * cosz_start) * (cosy_end * cosy_end + cosz_end * cosz_end));
    //start_end_angle_diff_YZ_cos->Fill(cos_thetayz_diff);
    start_end_angle_diff_YZ->Fill(TMath::ACos(cos_thetayz_diff));  
    if (TMath::ACos(cos_thetayz_diff)<cut_angle_diff) return 1;
    else return 0;
}

int x_order_de(vector<float> x)
{
    int siz = x.size();
    for (int i=0; i<siz-1; i++)
    {
	    if (x[i] <= x[i+1]) return 1;
    }
    return 0;
}

#endif
