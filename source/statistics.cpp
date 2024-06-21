#include "statistics.h"
#include <iostream>

using namespace std;

//--------------------------------------------------
double minimum(const vector<double> &x, int &idx)
{
	double out=x[0];
	idx = 0;
	for(int i=1; i<x.size(); i++)
	{
		if(x[i]<out)
		{
			out = x[i];
			idx = i;
		}
	}

	return out;
}

//--------------------------------------------------
double maximum(const vector<double> &x, int &idx)
{
	double out=x[0];
	idx = 0;
	for(int i=1; i<x.size(); i++)
	{
		if(x[i]>out)
		{
			out = x[i];
			idx = i;
		}
	}
	return out;
}

//--------------------------------------------------------------
int find_index(const vector<double> &x, double x0)
{
   int index=0;
   int i=0;
   while(i<x.size() && index==0)
   {
      if(x[i]<x0 && x[i+1]>x0)
      {
         index = i;
      }
      i++;
   }

   return index;
}

//--------------------------------------------------------------
double average(const vector<double> &x, const vector<double> &t)
{
   double m=0.;
   for(int i=0; i<t.size()-1; i++)
   {
      m += (t[i+1]-t[i])*(x[i+1]+x[i])*.5;
   }
   m /= (t.back()-t[0]);

   return m;
}

//--------------------------------------------------
vector<double> resample(const vector<double> &x, const vector<double> &t, double dt)
{
	vector<double> out;
	int i=0;
	double ts=t[0];
	double t_end=t.back();
	while(ts<t_end && i<t.size()-1)
	{
		if(t[i]<=ts && ts<t[i+1])
		{
         out.push_back((x[i]*(t[i+1]-ts) + x[i+1]*(ts-t[i]))/(t[i+1]-t[i]));
      	ts += dt;
		}
		else
		{
			i++;
		}
	}
	return out;
}

//--------------------------------------------------
int crop_after_T(const vector<double> &x, const vector<double> &t, double T)
{
	double t_end = t.back();
	int i_crop = 0;
	int i=0;
	while(i<t.size() && i_crop==0)
	{
		if(t[i]>=T)
		{
			i_crop = i;
		}
		i++;
	}
	return i_crop;
}

//--------------------------------------------------
double systole(const vector<double> &x, const vector<double> &t, double T)
{
	int i_crop = crop_after_T(x,t,T);
	vector<double> x2(x.begin()+i_crop,x.end());

	int idx;
	double sys = maximum(x2,idx);

	return sys;
}

//--------------------------------------------------
double diastole(const vector<double> &x, const vector<double> &t, double T)
{
	int i_crop = crop_after_T(x,t,T);
	vector<double> x2(x.begin()+i_crop,x.end());

	int idx;
	double dia = minimum(x2,idx);
	return dia;
}

//--------------------------------------------------
double time_delay_min(const vector<double> &x, const vector<double> &y, const vector<double> &t, double T)
{
	int ix = crop_after_T(x,t,T);
	vector<double> x2(x.begin()+ix,x.end());
	int iy = crop_after_T(y,t,T);
	vector<double> y2(y.begin()+iy,y.end());
	vector<double> t2(t.begin()+ix,t.end());

	minimum(x2,ix);
	minimum(y2,iy);

	return t2[iy]-t2[ix];
}

//--------------------------------------------------
double time_delay_max(const vector<double> &x, const vector<double> &y, const vector<double> &t, double T)
{
	int ix = crop_after_T(x,t,T);
	vector<double> x2(x.begin()+ix,x.end());
	int iy = crop_after_T(y,t,T);
	vector<double> y2(y.begin()+iy,y.end());
	vector<double> t2(t.begin()+ix,t.end());

	maximum(x2,ix);
	maximum(y2,iy);

	return t2[iy]-t2[ix];
}

//--------------------------------------------------
double pearson_correlation(const vector<double> &x, const vector<double> &y)
{
  int n = x.size();
  double num=0.0,sx=0.0,sy=0.0,ax=0.0,ay=0.0;
  for(int i=0;i<n;i++)
  {
    ax += x[i];
    ay += y[i];
  }
  ax /= n; ay /= n;
  for(int i=0;i<n;i++)
    num += (x[i]-ax)*(y[i]-ay);
  for(int i=0;i<n;i++)
  {
    sx += pow(x[i]-ax,2.);
    sy += pow(y[i]-ay,2.);
  }
  sx = pow(sx,0.5);
  sy = pow(sy,0.5);
  return num/sx/sy;
}

//--------------------------------------------------
vector<double> cross_correlation(const vector<double> &x, const vector<double> &y)
{
	vector<double> out(x.size(),.0);

	for(int i=0; i<x.size(); i++)
	{
		vector<double> y2;
		for(int j=i; j<y.size(); j++)
		{
			y2.push_back(y[j]);
		}
		for(int j=0; j<i; j++)
		{
			y2.push_back(y[j]);
		}
		out[i] = pearson_correlation(x,y2);
	}

	return out;
}

//--------------------------------------------------
double time_delay_correl(const vector<double> &x, const vector<double> &y, const vector<double> &t, double T)
{
	int ix = crop_after_T(x,t,T);
	vector<double> x2(x.begin()+ix,x.end());
	int iy = crop_after_T(y,t,T);
	vector<double> y2(y.begin()+iy,y.end());
	vector<double> t2(t.begin()+ix,t.end());

	vector<double> cc = cross_correlation(x2,y2);
	int idx_max;
	double cc_max = maximum(cc,idx_max);

	return t[idx_max]-t[0];
}

//--------------------------------------------------
time_average::time_average(){}

time_average::~time_average()
{
	time.clear();
	average.clear();
	value.clear();
}

//--------------------------------------------------
void time_average::update(double tnew, double vnew, double T)
{
	if(tnew<T)
	{
		if(average.size()>0)
		{
			double Anew = (vnew+value.back())*.5*(tnew-time.back());
			double avenew = (average.back()*(time.back()-time[0]) + Anew)/(tnew-time[0]);
			average.push_back(avenew);
			area += Anew;
		}
		else
		{
			average.push_back(vnew);
			area = vnew*.5*tnew;
		}
		value.push_back(vnew);
		time.push_back(tnew);
	}
	else
	{
		double Anew = (vnew+value.back())*.5*(tnew-time.back());
		area += Anew;

		while(tnew-time[last_idx]>=T)
		{
			last_idx++;
			area -= (value[last_idx]+value[last_idx-1])*.5*(time[last_idx]-time[last_idx-1]);
		}

		double avenew = area/T;

		average.push_back(avenew);
		value.push_back(vnew);
		time.push_back(tnew);
	}
}

//--------------------------------------------------
void time_average::save_results(string file_name)
{
   FILE *out_file;
	out_file = fopen(file_name.c_str(),"w");
   for(int i=0; i<time.size(); i++)
   {
		fprintf(out_file, "%9.7e, %9.7e\n", time[i], average[i]);
   }
	fclose(out_file);
}

//--------------------------------------------------
void time_average::save_results(double dt, string file_name)
{
	vector<double> a_resample = resample(average,time,dt);

   FILE *out_file;
	out_file = fopen(file_name.c_str(),"w");
	double t=0.;
   for(int i=0; i<a_resample.size(); i++)
   {
		fprintf(out_file, "%9.7e, %9.7e\n", t, a_resample[i]);
		t += dt;
   }
	fclose(out_file);
}

//--------------------------------------------------
double average(const vector<double> &x){
	int n = x.size();

	double sum = 0;
	for (int i=0; i<n; i++){
		sum += x[i];
	}
	return sum/n;
}