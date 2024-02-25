#include "transport.h"

// the number of division points of the segments. dx can not be larger then the corresponding 1D segment but N has a minimal value
int NX(double L,double dx, int minN) {
    if (ceil(L / dx) > minN) {
        return ceil(L / dx);
    }
    else {
        return minN;
    }
}

//constructor---------------------------------
TransportNodeCl::TransportNodeCl(TransportType TType) : TType(TType) {};

//--------------------------------------------
void TransportNodeCl::UpdateFi(int NodeIndex, double& fiNode){

    if (nodes[NodeIndex]->is_master_node == true) {
        cout << "Call function -- in class --";
        return;
        };
    if (nodes[NodeIndex]->upstream_boundary > -1) {} // no idea about this

    else if (nodes[NodeIndex]->type_code == 1) {} //periphery

    else if (nodes[NodeIndex]->type_code == 0) { // junction

        int n1 = nodes[NodeIndex]->edge_in.size();
        int n2 = nodes[NodeIndex]->edge_out.size();

        double q_sum = 0.; //vfr sum of the incoming edges only
        double fiNodeOld = fiNode;
        fiNode = 0;

       for (int j = 0; j < n1; j++)
       {
           if (edges[nodes[NodeIndex]->edge_in[j]]->v.back() > 0.)
           {
            q_sum += edges[nodes[NodeIndex]->edge_in[j]]->v.back() * edges[nodes[NodeIndex]->edge_in[j]]->A.back();
            fiNode += edges[nodes[NodeIndex]->edge_in[j]]->v.back() * edges[nodes[NodeIndex]->edge_in[j]]->A.back() * edges[nodes[NodeIndex]->edge_in[j]]->fi.back();
           }
        }
        for (int j = 0; j < n2; j++)
        {
            if (edges[nodes[NodeIndex]->edge_out[j]]->v[0] < 0.){
            q_sum -= edges[nodes[NodeIndex]->edge_out[j]]->v[0] * edges[nodes[NodeIndex]->edge_out[j]]->A[0]; //not sure about the sign tho...
            fiNode -= edges[nodes[NodeIndex]->edge_out[j]]->v[0] * edges[nodes[NodeIndex]->edge_out[j]]->A[0] * edges[nodes[NodeIndex]->edge_out[j]]->fi[0];
            }
        }
        if (q_sum != 0) {
            fiNode /= q_sum;
        }
        else { fiNode = fiNodeOld; }//if nothing flows in it stays the old
    }
};

//-----------------------------------------
void Virt1DforLum(vector<double> &fi_old, vector<double> &fi, double v, double dt, double dx, int n, double fiStartNode, double fiEndNode) {
    vector<double> fi_tmp = fi;

    for (int i = 1; i < n - 1; i++) {

        if (v > 0) {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i] - fi_old[i - 1]);
        }
        else {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i + 1] - fi_old[i]);
        }
    }

    //BC
    if (v > 0) {

        fi[n - 1] = fi_old[n - 1] - v * dt / dx * (fi_old[n - 1] - fi_old[n - 2]);
        fi[0] = fiStartNode;
    }
    else {
        fi[n - 1] = fiEndNode; 
        fi[0] = fi_old[0] - v * dt / dx * (fi_old[1] - fi_old[0]);
    }

    fi_old = fi_tmp;

}

//-------------------------------------------
TransportNodeCl::TransportNodeCl(TransportType TType):  TType(TType) {};


//-------------------------------------------
void TransportNodeCl::UpdateFi(vector<double> v, vector<double>& fi_old, vector<double>& fi, double l, vector<double> t_act, double dt){
    int n = (int)v.size();
    double dx = l / (n - 1);
    vector<double> fi_tmp = fi; //this will be the new fi_old

    for (int i = 1; i < n - 1; i++) {
        if (v[i] > 0 ) {
            fi[i] = fi_old[i] - v[i] * dt / dx * (fi_old[i] - fi_old[i - 1]);
        }
        else {
            fi[i] = fi_old[i] - v[i] * dt / dx * (fi_old[i + 1] - fi_old[i]);
        }
    }
 
    // downstream BC
    if (v[n - 1] > 0) {
        fi[n - 1] = fi_old[n - 1] - v[n - 1] * dt / dx * (fi_old[n - 1] - fi_old[n - 2]);
    }
    else {
        fi[n - 1] = fi_old[n - 1]; // there should be fi_node right?
    }

    // upstream BC
    if (v[0] > 0) {
        fi[0] = fi0;  // there should be fi_node right?
    }
    else{
        fi[0] = fi_old[0] - v[0] * dt / dx * (fi_old[1] - fi_old[0]);
    }

    fi_old = fi_tmp;
}

//--------------------------------
D0Transport::D0Transport(LumpedType LType, vector<string> sv):LType(LType) {
    switch (LType) {
    case PerifCoronary0D:
    //no idea
        break;
    case Perif0D:

        L_arteriole = stod(sv[1],0);
        L_capillary = stod(sv[2],0);
        L_venulare = stod(sv[3],0);
        L_vein = stod(sv[4],0);
        A_arteriole = stod(sv[5],0);
        A_capillary = stod(sv[6],0);
        A_venulare = stod(sv[7],0);
        A_vein = stod(sv[8],0);

        //arteriole
        nx_arteriole = NX( L_arteriole, dx, 5); //at least 5 division points
        dx_arteriole = L_arteriole / (nx_arteriole - 1);
        fi_arteriole.assign(nx_arteriole, 0.);
        fi_old_arteriole.assign(nx_arteriole, 0.);

        //capillary
        nx_capillary = NX(L_capillary, dx, 5);
        dx_capillary = L_capillary / (nx_capillary - 1);
        fi_capillary.assign(nx_arteriole, 0.);
        fi_old_capillary.assign(nx_arteriole, 0.);

        //venulare
        nx_venulare = NX(L_venulare, dx, 5);
        dx_venulare = L_venulare / (nx_venulare - 1);
        fi_venulare.assign(nx_arteriole, 0.);
        fi_old_venulare.assign(nx_arteriole, 0.);

        //vein
        nx_vein = NX(L_vein, dx, 5);
        dx_vein = L_vein / (nx_vein - 1);
        fi_vein.assign(nx_arteriole, 0.);
        fi_old_vein.assign(nx_arteriole, 0.);


        break;
    case Heart0D:


        break;
    default:
        cout << "Unknown 0D type for transport.";
    }
}


//-----------------------------------
void D0Transport::UpdateFi(int LumIndex, double dt) {
    switch (this-> LType) {
    case PerifCoronary0D:
        //no idea
        break;
    case Perif0D:
      
        Virt1DforLum(fi_old_arteriole, fi_arteriole, lum[LumIndex]->edges[1]->vfr / A_arteriole, dt, dx_arteriole, nx_arteriole, lum[LumIndex]->nodes[1]->fi, lum[LumIndex]->nodes[2]->fi);
        Virt1DforLum(fi_old_capillary, fi_capillary, lum[LumIndex]->edges[2]->vfr / A_capillary, dt, dx_capillary, nx_capillary, lum[LumIndex]->nodes[2]->fi, lum[LumIndex]->nodes[3]->fi);
        Virt1DforLum(fi_old_venulare, fi_venulare, lum[LumIndex]->edges[3]->vfr / A_venulare, dt, dx_venulare, nx_venulare, lum[LumIndex]->nodes[3]->fi, lum[LumIndex]->nodes[4]->fi);
        Virt1DforLum(fi_old_vein, fi_vein, lum[LumIndex]->edges[4]->vfr / A_vein, dt, dx_vein, nx_vein, lum[LumIndex]->nodes[4]->fi, lum[LumIndex]->nodes[]->fi);
            
        //nodes
        // not sure if the virtual 1D or the nodes should be updated first
        UpdatePerifLumNodes(1, LumIndex, MMMMM , fi_arteriole[0]);// MMMMM is the master node's fi
        UpdatePerifLumNodes(2, LumIndex, fi_arteriole.back(), fi_capillary[0]);
        UpdatePerifLumNodes(3, LumIndex, fi_capillary.back(), fi_venulare[0]);
        UpdatePerifLumNodes(4, LumIndex, fi_venulare.back(), fi_vein[0]);

        //master node
        //node connecting heart and perifs


        break;
    case Heart0D:// that will be fun...
        //right atrium


        break;
    default:
        cout << "Unknown 0D type for transport.";
    }
    return;
}

//------------------------------------
//updates fi parameters of perif nodes
void D0Transport::UpdatePerifLumNodes(int LumNodeIndex, int LumIndex, double fiLeft, double fiRight) {
    int a, b, c;

    //1 if q flows towards the node
    if (lum[LumIndex]->edges[LumNodeIndex]->vfr < 0) { a = 1; } //Resistance edge
    if (lum[LumIndex]->edges[LumNodeIndex - 1]->vfr > 0) { b = 1; } //Inductance edge
    if (lum[LumIndex]->edges[LumNodeIndex + 4]->vfr < 0) { c = 1; } //Capacitance edge

    double qLeft, qRight, qDown;

    switch (a*4 + b*2 + c) {
    case 1:
        lum[LumIndex]->nodes[LumNodeIndex]->fi = fiRight;
        break;

    case 2:
        lum[LumIndex]->nodes[LumNodeIndex]->fi = fiLeft;
        break;

    case 3:
        qLeft = lum[LumIndex]->edges[LumNodeIndex - 1]->vfr;
        qDown = lum[LumIndex]->edges[LumNodeIndex - 4]->vfr;
        lum[LumIndex]->nodes[LumNodeIndex]->fi = (qLeft*fiLeft + qDown *fiRight)/(qDown + qLeft);
        break;

    case 4:
        lum[LumIndex]->nodes[LumNodeIndex]->fi = fiRight; //same az case 1
        break;

    case 5:
        lum[LumIndex]->nodes[LumNodeIndex]->fi = fiRight; //both has the same fi
        break;

    case 6:
        qLeft = lum[LumIndex]->edges[LumNodeIndex - 1]->vfr;
        qRight = lum[LumIndex]->edges[LumNodeIndex]->vfr;
        lum[LumIndex]->nodes[LumNodeIndex]->fi = (qLeft * fiLeft + qRight * fiRight) / (qRight + qLeft);
        break;
    }
}


/*
stuff to do in other classes:
- fi vector for all edges.
- fi variable for all nodes.
- fi variable for lumped 0D models
- saving stuff
- main fileba "RBCtransport on" kell beleírni
- get fv a sebességhez a moc_edge osztályba

kel v_left, right, x_left right

-first_blood class should have a BIG virtual node for connecting perifs to the heart left point

- prescribed boundaries?????
*/

/*
Ready:
-1D
-perif nodes except for master nodes or the connection between perif and heart
-1D nodes
-virtual 1D for perifs
- first_blood class-ban egy változó, ami eldönti vane transport

*/

