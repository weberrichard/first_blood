#include <vector>
#include <iostream>

using namespace std;

vector<double> resample(vector<double> x, vector<double> t, double dt);

int main()
{
	vector<double> x{1.,3.,4.,5.,2.,3.};
	vector<double> t{1.,2.,3.,4.,5.,6.};
	double dt = 1.51;
	vector<double> x2 = resample(x,t,dt);
	for(int i=0; i<x2.size(); i++)
	{
		cout << x2[i] << endl;
	}

	cout << endl;
	return 0;
}

//--------------------------------------------------
vector<double> resample(vector<double> x, vector<double> t, double dt)
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
