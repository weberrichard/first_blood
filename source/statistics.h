/*===================================================================*\
									 	statistics
									 ---------------
 
	 first_blood
	 R. Weber
	 https://github.com/weberrichard/first_blood
\*==================================================================*/

#ifndef STATISTICS
#define STATISTICS

#include <vector>
#include <cmath>
#include <string>

using namespace std;

double minimum(const vector<double> &x, int &idx);
double maximum(const vector<double> &x, int &idx);
int find_index(const vector<double> &x, double x0);
double average(const vector<double> &x, const vector<double> &t);
vector<double> resample(const vector<double> &x, const vector<double> &t, double dt);
int crop_index(const vector<double> &x, const vector<double> &t, double T);
double systole(const vector<double> &x, const vector<double> &t, double T);
double diastole(const vector<double> &x, const vector<double> &t, double T);
double time_delay_min(const vector<double> &x, const vector<double> &y, const vector<double> &t, double T);
double time_delay_max(const vector<double> &x, const vector<double> &y, const vector<double> &t, double T);
double pearson_correlation(const vector<double> &x, const vector<double> &y);
vector<double> cross_correlation(const vector<double> &x, const vector<double> &y);
double time_delay_correl(const vector<double> &x, const vector<double> &y, const vector<double> &t, double T);

class time_average
{
public:
	time_average();
	~time_average();
	vector<double> average, value, time; // actual averaged time series
	double area=0.;
	int last_idx=0;
	void update(double tnew, double vnew, double T); // inicializalas!!!
	void save_results(string file_name);
	void save_results(double dt, string file_name);
};

#endif