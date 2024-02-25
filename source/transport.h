/*

Classes and functions responsible for the simulation of the transport of different substances are declared in this file.

*/

#ifndef TRANSPORT_H
#define TRANSPORT_H

#include<iostream>
#include<cstdio>
#include<string>
#include<vector>
#include<fstream>
using namespace std;

enum TransportType { RBC };
enum LumpedType { PerifCoronary0D, Perif0D, Heart0D };

int NX(double L,double dx, int minN);

class TransportNodeCl {//for 1D 
public:
    TransportType TType;

    TransportNodeCl(TransportType TType);

    void UpdateFi(int NodeIndex, double& fiNode);
};

void Virt1DforLum(vector<double> &fi_old, vector<double> &fi, double v, double dt, double dx, int n, double fiStartNode, double fiEndNode);


//a first_blood object gets one of this class. This handles 1D transport for the moc edges
class Transport1DCl {
public:
    TransportType TType;

    Transport1DCl(TransportType TType);

    void UpdateFi(vector<double> x, vector<double> v, vector<double>& fi_old, vector<double>& fi, double l, vector<double> t_act, double dt);
};



// Every lumped model gets one if necesarry. This handles the transport of substances in 0D
class D0Transport {//every 0D model gets one of this class
public:
    //for simple peripherals wirh 4 RLC circuits
    vector<double> fi_arteriole, fi_capillary, fi_venulare, fi_vein;
    vector<double> fi_old_arteriole, fi_old_capillary, fi_old_venulare, fi_old_vein;
    double dx_arteriole,dx_capillary, dx_venulare, dx_vein;
    double L_arteriole, L_capillary, L_venulare, L_vein;
    double A_arteriole, A_capillary, A_venulare, A_vein;//from file
    int nx_arteriole, nx_capillary, nx_venulare, nx_vein;
    LumpedType LType;

    //for the heart model

    vector<double> fi_RA, fi_RV, fi_LA, fi_LV;//right atrium, right ventricle, left atrium, left ventricle
    vector<double> fi_old_RA, fi_old_RV, fi_old_LA, fi_old_LV;

    D0Transport(LumpedType LType, vector<string> sv);

    void UpdateFi(int LumIndex, double dt);

    void UpdatePerifLumNodes(int LumNodeIndex, int LumIndex, double fiLeft, double fiRight);

};

#endif