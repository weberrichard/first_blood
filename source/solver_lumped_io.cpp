#include "solver_lumped.h"

//--------------------------------------------------------------
void solver_lumped::load_model()
{
	ifstream file_in;
	string file_name = input_folder_path + '/' + name + ".csv";
	file_in.open(file_name);
	string line;
	if(file_in.is_open())
	{
		int nn=0,ne=0; // ne for edges, nn for nodes
		while(getline(file_in,line))
		{	
			// cleaning unnecessary characters
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());
			vector<string> sv = separate_line(line);

			if(sv[0] == "resistor" || sv[0] == "capacitor" || sv[0] == "inductor" || sv[0] == "voltage" || sv[0] == "diode" || sv[0] == "resistor2" || sv[0] == "valve" || sv[0] == "resistor_coronary" || sv[0] == "capacitor_coronary" || sv[0] == "current") // edges with one parameter
			{
				edges.push_back(new edge);
				edges[ne]->type = sv[0];
				edges[ne]->name = sv[1];
				edges[ne]->node_name_start = sv[2];
				edges[ne]->node_name_end = sv[3];
				edges[ne]->volume_flow_rate_initial = stod(sv[4],0);
				if(sv[0] == "resistor")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 0;
				}
				else if(sv[0] == "capacitor")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 1;
				}
				else if(sv[0] == "inductor")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 3;
				}
				else if(sv[0] == "voltage")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));					
					edges[ne]->type_code = 4;
				}
				else if(sv[0] == "diode")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 5;
				}
				else if(sv[0] == "resistor2" || sv[0] == "valve")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 6;
				}
				else if(sv[0] == "resistor_coronary")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 7;
				}
				else if(sv[0] == "capacitor_coronary")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));
					edges[ne]->type_code = 8;
				}
				else if(sv[0] == "current")
				{
					edges[ne]->parameter.push_back(stod(sv[5],0));					
					edges[ne]->type_code = 9;
				}
				ne++;
			}
			else if(sv[0] == "elastance")
			{
				edges.push_back(new edge);
				edges[ne]->type = sv[0];
				edges[ne]->name = sv[1];
				edges[ne]->node_name_start = sv[2];
				edges[ne]->node_name_end = sv[3];
				edges[ne]->volume_flow_rate_initial = stod(sv[4],0);
				if(sv.size()>5)
				{
					edges[ne]->parameter.push_back(stod(sv[5],0)); // elastance max SI 
					edges[ne]->parameter.push_back(stod(sv[6],0)); // elastance min SI
				}
				else
				{
					edges[ne]->parameter.push_back(elastance_max_nom*mmHg_to_Pa*1e6); // in SI
					edges[ne]->parameter.push_back(elastance_min_nom*mmHg_to_Pa*1e6); // in SI
				}
				edges[ne]->type_code = 2;
				ne++;
			}
			else if(sv[0] == "resistor_piecewise_constant")
			{
				edges.push_back(new edge);
				edges[ne]->type = sv[0];
				edges[ne]->name = sv[1];
				edges[ne]->node_name_start = sv[2];
				edges[ne]->node_name_end = sv[3];
				edges[ne]->volume_flow_rate_initial = stod(sv[4],0);
				edges[ne]->type_code = 10;

				if(sv.size()>8)
				{
					edges[ne]->parameter.push_back(stod(sv[5],0)); // R1 
					edges[ne]->parameter.push_back(stod(sv[6],0)); // R2
					edges[ne]->t1 = stod(sv[7],0);
					edges[ne]->t2 = stod(sv[8],0);
				}
				ne++;
			}
			else if(sv[0] == "node") // node
			{
				nodes.push_back(new node);
				nodes[nn]->type = sv[0];
				nodes[nn]->name = sv[1];
				nodes[nn]->pressure_initial = stod(sv[2],0);
				nodes[nn]->is_ground = false;
				nn++;
			}
			else if(sv[0] == "ground") // node with ground
			{
				nodes.push_back(new node);
				nodes[nn]->type = sv[0];
				nodes[nn]->name = sv[1];
				nodes[nn]->pressure_initial = stod(sv[2],0);
				nodes[nn]->is_ground = true;
				nn++;
			}
			else if(sv[0] == "myogenic") // setting myogenic control
			{
				if(sv[1] == "on")
				{
					do_myogenic = true;
				}
				if(sv.size()>3)
				{
					q_ref = stod(sv[2],0)*1e6; // ref volumetric flow rate, convert to NON SI
					p_ref = stod(sv[3],0)/mmHg_to_Pa; // ref pressure, convert to NON SI
				}
				if(sv.size()>7)
				{
					tao   = stod(sv[4],0); // time constant
					G     = stod(sv[5],0); // gain
					//sat1  = stod(sv[6],0); // saturation 1
					//sat2  = stod(sv[7],0); // saturation 2
				}
			}

			//RBC, HBsaturation, PlasmaO2C perif transport
			else if(sv[0] == "RBC_on" && sv[1] == "1"){
			//RBC
				do_lum_RBC_transport = true;
				RBClum = new D0_transport(RBC);
				}

			//HBsaturation
			else if (sv[0] == "HBsaturation_on" && sv[1] == "1"){
				do_lum_HB_sat_transport = true;
				HBsatlum = new D0_transport(HB_O2_saturation);
				}

			//PlasmaO2C
			else if (sv[0] == "PlasmaO2_on" && sv[1] == "1"){
				do_lum_PlasmaO2_transport = true;
				PlasmaO2lum = new D0_transport(C_Plasma_O2);
				}

			else if(sv[0]=="O2transport_init" && sv.size() > 3 ){
				fi_init_RBC_lum = stod(sv[1],0);
				init_PlasmaO2_lum = stod(sv[2],0);
				init_HB_sat_lum = stod(sv[3],0);
			}

			//CO2 transport stuff
			else if(sv[0]=="co2transport_init" && sv.size() > 3 ){
				fi_init_CO2_pla = stod(sv[1],0);
				fi_init_CO2_rbc = stod(sv[2],0);
				fi_init_HCO3_pla = stod(sv[3],0);
				fi_init_HCO3_rbc = stod(sv[4],0);
				fi_init_HbCO2 = stod(sv[5],0);
			}

			else if (sv[0] == "pla_CO2_on" && sv[1] == "1"){
				do_lum_pla_CO2_transport = true;
				CO2_pla_lum = new D0_transport(CO2_pla);
			}

			else if (sv[0] == "rbc_CO2_on" && sv[1] == "1"){
				do_lum_rbc_CO2_transport = true;
				CO2_rbc_lum = new D0_transport(CO2_rbc);
			}

			else if (sv[0] == "pla_HCO3_on" && sv[1] == "1"){
				do_lum_pla_HCO3_transport = true;
				HCO3_pla_lum = new D0_transport(HCO3_pla);
			}

			else if (sv[0] == "rbc_HCO3_on" && sv[1] == "1"){
				do_lum_rbc_HCO3_transport = true;
				HCO3_rbc_lum = new D0_transport(HCO3_rbc);
			}

			else if (sv[0] == "HbCO2_on" && sv[1] == "1"){
				do_lum_HbCO2_transport = true;
				HbCO2_lum = new D0_transport(HbCO2);
			}


			//metabolic response
			else if(sv[0] == "metabolic")
			{
				if(sv[1]=="on"){
					    Ct_ref = stod(sv[2],0);
    					 tao_met = stod(sv[3],0);
                   G_met = stod(sv[4],0);
                   sat1_met = stod(sv[5],0);
                   sat2_met = stod(sv[6],0);
    					 //Ct_ave;
    					 do_metabolic_res = true;
				}
			}

			//CO2 response
			else if(sv[0] == "CO2_response")
			{
					if(sv[1]=="on"){
				    CO2_ref = stod(sv[2],0);
    				tao_CO2 = stod(sv[3],0);
                    G_CO2 = stod(sv[4],0);
                    do_CO2_control = true;
				}
			}

			//virtual 1D edges
			else if(sv[0] == "v1D"){

				if(do_lum_RBC_transport&&sv.size()>9){
				RBClum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), RBC, fi_init_RBC_lum, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){RBClum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){RBClum->D0_edges.back()->is_pul_capillary = true;}
				}

				if(do_lum_HB_sat_transport&&sv.size()>9){
				HBsatlum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), HB_O2_saturation, init_HB_sat_lum, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){HBsatlum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){HBsatlum->D0_edges.back()->is_pul_capillary = true;}
				}
				
				if(do_lum_PlasmaO2_transport&&sv.size()>9){
				PlasmaO2lum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), C_Plasma_O2, init_PlasmaO2_lum, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){PlasmaO2lum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){PlasmaO2lum->D0_edges.back()->is_pul_capillary = true;}
				}

				//CO2 transport
				if(do_lum_pla_CO2_transport&&sv.size()>9){
				CO2_pla_lum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), CO2_pla, fi_init_CO2_pla, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){CO2_pla_lum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){CO2_pla_lum->D0_edges.back()->is_pul_capillary = true;}
				}

				if(do_lum_rbc_CO2_transport&&sv.size()>9){
				CO2_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), CO2_rbc, fi_init_CO2_rbc, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){CO2_rbc_lum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){CO2_rbc_lum->D0_edges.back()->is_pul_capillary = true;}
				}

				if(do_lum_pla_HCO3_transport&&sv.size()>9){
				HCO3_pla_lum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), HCO3_pla, fi_init_HCO3_pla, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){HCO3_pla_lum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){HCO3_pla_lum->D0_edges.back()->is_pul_capillary = true;}
				}

				if(do_lum_rbc_HCO3_transport&&sv.size()>9){
				HCO3_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), HCO3_rbc, fi_init_HCO3_rbc, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){HCO3_rbc_lum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){HCO3_rbc_lum->D0_edges.back()->is_pul_capillary = true;}
				}

				if(do_lum_HbCO2_transport&&sv.size()>9){
				HbCO2_lum->D0_edges.push_back(new D0_edge(sv[1], stod(sv[4],0), stod(sv[5],0), stoi(sv[6],0), HbCO2, fi_init_HbCO2, sv[2], sv[3], sv[9]));
				if(sv[7] == "1" ){HbCO2_lum->D0_edges.back()->is_per_capillary = true;}
				if(sv[8] == "1" ){HbCO2_lum->D0_edges.back()->is_pul_capillary = true;}
				}



			}

			//virtual 1D diodes
			else if(sv[0] == "t_diode"){
				if(do_lum_RBC_transport){
				RBClum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, RBC, fi_init_RBC_lum, sv[2], sv[3], sv[4]));
				RBClum->D0_edges.back()->is_diode = true;	
				}

				if(do_lum_HB_sat_transport){
				HBsatlum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HB_O2_saturation, init_HB_sat_lum, sv[2], sv[3], sv[4]));
				HBsatlum->D0_edges.back()->is_diode = true;
				}
				
				if(do_lum_PlasmaO2_transport){
				PlasmaO2lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, C_Plasma_O2, init_PlasmaO2_lum, sv[2], sv[3], sv[4]));
				PlasmaO2lum->D0_edges.back()->is_diode = true;			
				}


				//CO2 transport
				if(do_lum_pla_CO2_transport){
				CO2_pla_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, CO2_pla, fi_init_CO2_pla, sv[2], sv[3], sv[4]));
				CO2_pla_lum->D0_edges.back()->is_diode = true;			
				}

				if(do_lum_rbc_CO2_transport){
				CO2_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, CO2_rbc, fi_init_CO2_rbc, sv[2], sv[3], sv[4]));
				CO2_rbc_lum->D0_edges.back()->is_diode = true;			
				}

				if(do_lum_pla_HCO3_transport){
				HCO3_pla_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HCO3_pla, fi_init_HCO3_pla, sv[2], sv[3], sv[4]));
				HCO3_pla_lum->D0_edges.back()->is_diode = true;
				}

				if(do_lum_rbc_HCO3_transport){
				HCO3_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HCO3_rbc, fi_init_HCO3_rbc, sv[2], sv[3], sv[4]));
				HCO3_rbc_lum->D0_edges.back()->is_diode = true;
				}

				if(do_lum_HbCO2_transport){
				HbCO2_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HbCO2, fi_init_HbCO2, sv[2], sv[3], sv[4]));
				HbCO2_lum->D0_edges.back()->is_diode = true;
				}

			}

			

			else if(sv[0] == "0Dcapacitor"){
				if(do_lum_RBC_transport){
				RBClum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, RBC, fi_init_RBC_lum, sv[2], sv[3], sv[4]));
				RBClum->D0_edges.back()->is_capacitor = true;
				}

				if(do_lum_HB_sat_transport){
				HBsatlum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HB_O2_saturation, init_HB_sat_lum, sv[2], sv[3], sv[4]));
				HBsatlum->D0_edges.back()->is_capacitor = true;
				}
				
				if(do_lum_PlasmaO2_transport){
				PlasmaO2lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, C_Plasma_O2, init_PlasmaO2_lum, sv[2], sv[3], sv[4]));
				PlasmaO2lum->D0_edges.back()->is_capacitor = true;
				}

				//CO2 transport
				if(do_lum_pla_CO2_transport){
				CO2_pla_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, CO2_pla, fi_init_CO2_pla, sv[2], sv[3], sv[4]));
				CO2_pla_lum->D0_edges.back()->is_capacitor = true;			
				}

				if(do_lum_rbc_CO2_transport){
				CO2_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, CO2_rbc, fi_init_CO2_rbc, sv[2], sv[3], sv[4]));
				CO2_rbc_lum->D0_edges.back()->is_capacitor = true;			
				}

				if(do_lum_pla_HCO3_transport){
				HCO3_pla_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HCO3_pla, fi_init_HCO3_pla, sv[2], sv[3], sv[4]));
				HCO3_pla_lum->D0_edges.back()->is_capacitor = true;
				}

				if(do_lum_rbc_HCO3_transport){
				HCO3_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HCO3_rbc, fi_init_HCO3_rbc, sv[2], sv[3], sv[4]));
				HCO3_rbc_lum->D0_edges.back()->is_capacitor = true;
				}

				if(do_lum_HbCO2_transport){
				HbCO2_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HbCO2, fi_init_HbCO2, sv[2], sv[3], sv[4]));
				HbCO2_lum->D0_edges.back()->is_capacitor = true;
				}

			}


			else if(sv[0] == "0Delastance"){
				if(do_lum_RBC_transport){
				RBClum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, RBC, fi_init_RBC_lum, sv[2], sv[3], sv[4]));
				RBClum->D0_edges.back()->is_elastance = true;
				RBClum->D0_edges.back()->V0 = stod(sv[5],0);
				}

				if(do_lum_HB_sat_transport){
				HBsatlum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HB_O2_saturation, init_HB_sat_lum, sv[2], sv[3], sv[4]));
				HBsatlum->D0_edges.back()->is_elastance = true;
				HBsatlum->D0_edges.back()->V0 = stod(sv[5],0);
				}
				
				if(do_lum_PlasmaO2_transport){
				PlasmaO2lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, C_Plasma_O2, init_PlasmaO2_lum, sv[2], sv[3], sv[4]));
				PlasmaO2lum->D0_edges.back()->is_elastance = true;
				PlasmaO2lum->D0_edges.back()->V0 = stod(sv[5],0);
				}


				//CO2 transport
				if(do_lum_pla_CO2_transport){
				CO2_pla_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, CO2_pla, fi_init_CO2_pla, sv[2], sv[3], sv[4]));
				CO2_pla_lum->D0_edges.back()->is_elastance = true;		
				CO2_pla_lum->D0_edges.back()->V0 = stod(sv[5],0);
				}

				if(do_lum_rbc_CO2_transport){
				CO2_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, CO2_rbc, fi_init_CO2_rbc, sv[2], sv[3], sv[4]));
				CO2_rbc_lum->D0_edges.back()->is_elastance = true;
				CO2_rbc_lum->D0_edges.back()->V0 = stod(sv[5],0);
				}

				if(do_lum_pla_HCO3_transport){
				HCO3_pla_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HCO3_pla, fi_init_HCO3_pla, sv[2], sv[3], sv[4]));
				HCO3_pla_lum->D0_edges.back()->is_elastance = true;
				HCO3_pla_lum->D0_edges.back()->V0 = stod(sv[5],0);
				}

				if(do_lum_rbc_HCO3_transport){
				HCO3_rbc_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HCO3_rbc, fi_init_HCO3_rbc, sv[2], sv[3], sv[4]));
				HCO3_rbc_lum->D0_edges.back()->is_elastance = true;
				HCO3_rbc_lum->D0_edges.back()->V0 = stod(sv[5],0);
				}

				if(do_lum_HbCO2_transport){
				HbCO2_lum->D0_edges.push_back(new D0_edge(sv[1], 0., 0., 1, HbCO2, fi_init_HbCO2, sv[2], sv[3], sv[4]));
				HbCO2_lum->D0_edges.back()->is_elastance = true;
				HbCO2_lum->D0_edges.back()->V0 = stod(sv[5],0);
				}

			}

			else if(sv[0] == "Mmax"){Mmax = stod(sv[1],0);}


		}

		if(do_lum_PlasmaO2_transport&&do_lum_HB_sat_transport&&do_lum_RBC_transport){
			init_lum_tissueO2();
			RBClum->do_tissue_transport = true;
			HBsatlum->do_tissue_transport = true;
			PlasmaO2lum->do_tissue_transport = true;

			//if we calculate PlasmaO2, HB_sat and RBC and there is a pul or per capillary in the model
			int n_pul_cap = 0;
			int n_per_cap = 0;

			for(int zz=0; zz<HBsatlum->D0_edges.size(); zz++){
				if(HBsatlum->D0_edges[zz] ->is_pul_capillary){
					do_pul_O2_rtansport=true;
					n_pul_cap++;
					pul_cap_BH = HBsatlum->D0_edges[zz];
					}
				else if(HBsatlum->D0_edges[zz]->is_per_capillary){
					do_per_O2_rtansport=true;
					n_per_cap++;
					per_cap_BH = HBsatlum->D0_edges[zz];}
			}

			for(int zz=0; zz<PlasmaO2lum->D0_edges.size(); zz++){
				if(PlasmaO2lum->D0_edges[zz] ->is_pul_capillary){
					pul_cap_PO2 = PlasmaO2lum->D0_edges[zz];
					}
				else if(PlasmaO2lum->D0_edges[zz]->is_per_capillary){
					per_cap_PO2 = PlasmaO2lum->D0_edges[zz];}
			}

			for(int zz=0; zz<RBClum->D0_edges.size(); zz++){
				if(RBClum->D0_edges[zz] ->is_pul_capillary){
					pul_cap_RBC = RBClum->D0_edges[zz];
					}
				else if(RBClum->D0_edges[zz]->is_per_capillary){
					per_cap_RBC = RBClum->D0_edges[zz];}
			}

			if(n_pul_cap > 1 || n_per_cap > 1){
				std::cout << "! ERROR !" << endl << "Max one pulmonary or peripheral capillary in a lumped models: " << name << endl;
				exit(-1); }

		}


		//CO2 stuff
		if(do_lum_pla_CO2_transport&&do_lum_rbc_CO2_transport&&do_lum_pla_HCO3_transport&&do_lum_rbc_HCO3_transport&&do_lum_HbCO2_transport){
			init_lum_tissueCO2();
			CO2_pla_lum->do_tissue_CO2_transport = true;
			CO2_rbc_lum->do_tissue_CO2_transport = true;
			HCO3_pla_lum->do_tissue_CO2_transport = true;
			HCO3_rbc_lum->do_tissue_CO2_transport = true;
			HbCO2_lum->do_tissue_CO2_transport = true;

			//if we calculate PlasmaO2, HB_sat and RBC and there is a pul or per capillary in the model
			int n_pul_cap = 0;
			int n_per_cap = 0;

			for(int zz=0; zz<CO2_pla_lum->D0_edges.size(); zz++){
				if(CO2_pla_lum->D0_edges[zz]->is_pul_capillary){
					do_pul_CO2_transport=true;
					n_pul_cap++;
					pul_cap_CO2_pla = CO2_pla_lum->D0_edges[zz];}
				else if(CO2_pla_lum->D0_edges[zz]->is_per_capillary){
					do_per_CO2_transport=true;
					n_per_cap++;
					per_cap_CO2_pla = CO2_pla_lum->D0_edges[zz];}
			}

			for(int zz=0; zz<CO2_rbc_lum->D0_edges.size(); zz++){
				if(CO2_rbc_lum->D0_edges[zz]->is_pul_capillary){
					n_pul_cap++;
					pul_cap_CO2_rbc = CO2_rbc_lum->D0_edges[zz];}
				else if(CO2_rbc_lum->D0_edges[zz]->is_per_capillary){
					n_per_cap++;
					per_cap_CO2_rbc = CO2_rbc_lum->D0_edges[zz];}
			}

			for(int zz=0; zz<HCO3_pla_lum->D0_edges.size(); zz++){
				if(HCO3_pla_lum->D0_edges[zz]->is_pul_capillary){
					n_pul_cap++;
					pul_cap_HCO3_pla = HCO3_pla_lum->D0_edges[zz];}
				else if(HCO3_pla_lum->D0_edges[zz]->is_per_capillary){
					n_per_cap++;
					per_cap_HCO3_pla = HCO3_pla_lum->D0_edges[zz];}
			}

			for(int zz=0; zz<HCO3_rbc_lum->D0_edges.size(); zz++){
				if(HCO3_rbc_lum->D0_edges[zz]->is_pul_capillary){
					n_pul_cap++;
					pul_cap_HCO3_rbc = HCO3_rbc_lum->D0_edges[zz];}
				if(HCO3_rbc_lum->D0_edges[zz]->is_per_capillary){
					n_per_cap++;
					per_cap_HCO3_rbc = HCO3_rbc_lum->D0_edges[zz];}
			}

			for(int zz=0; zz<HbCO2_lum->D0_edges.size(); zz++){
				if(HbCO2_lum->D0_edges[zz]->is_pul_capillary){
					n_pul_cap++;
					pul_cap_HbCO2 = HbCO2_lum->D0_edges[zz];}
				if(HbCO2_lum->D0_edges[zz]->is_per_capillary){
					n_per_cap++;
					per_cap_HbCO2 = HbCO2_lum->D0_edges[zz];}
			}




		}
	}
	else
	{
		std::cout << "! ERROR !" << endl << " File is not open when calling load_system_csv() function!!! file: " << file_name << "\nExiting..." << endl;
		exit(-1);
	}

	// setting size of elements
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();

	file_in.close();

   //O2 transport parameters reading from file
   file_name = input_folder_path + '/' + "O2_parameters" + ".txt";
	file_in.open(file_name);
	if(file_in.is_open()){
		while(getline(file_in,line))
		{	
			// cleaning unnecessary characters
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());
			vector<string> sv = separate_line(line);
			
			if (sv[0] == "peripheries"){ assign_perif_O2_params(sv); }
			else if(sv[0] == "haemoglobin"){ assign_haemogobin_sat_params(sv); }
			else if(sv[0] == "pulmonary"){ assign_pulmonary_O2_params(sv); }

		}
	}
	else if(do_lum_PlasmaO2_transport&&do_lum_HB_sat_transport&&do_lum_RBC_transport){
		cout<<"O2_parameter default values"<<endl;
	}
	file_in.close();

	file_name = input_folder_path + '/' + "CO2_parameters" + ".txt";
	file_in.open(file_name);
	if(file_in.is_open()){
		while(getline(file_in,line))
		{	
			// cleaning unnecessary characters
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());
			vector<string> sv = separate_line(line);
			
			if (sv[0] == "perifco2"){ assign_perif_CO2_params(sv);}

		}
		
	}
	else if(true){
		cout<<"CO2_parameter default values"<<endl;
	}
	file_in.close();


}

//--------------------------------------------------------------
void solver_lumped::save_results()
{
	save_results(name);
}

//--------------------------------------------------------------
void solver_lumped::save_results(string folder_name)
{
	vector<string> edge_list, node_list;
	for(int i=0; i<number_of_edges; i++)
	{
		edge_list.push_back(edges[i]->name);
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		node_list.push_back(nodes[i]->name);
	}
	save_results(folder_name,edge_list,node_list);
}

//--------------------------------------------------------------
void solver_lumped::save_results(string folder_name, vector<string> edge_list, vector<string> node_list)
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

   string fn = "results/" + folder_name + "/" + name;

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
   	if(nodes[idx]->do_save_memory)
   	{		
	   	mkdir(fn.c_str(),0777);
		   string file_name = fn + "/" + nodes[idx]->name + ".txt";
	      out_file = fopen(file_name.c_str(),"w");
	      for(unsigned int j=0; j<nodes[idx]->pressure.size(); j++)
	      {
	      	double t = time[j];
	      	double p = nodes[idx]->pressure[j];
	         fprintf(out_file, "%9.7e, %9.7e\n", t, p);
	      }
	      fclose(out_file);
   	}
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
   	if(edges[idx]->do_save_memory)
   	{
	   	mkdir(fn.c_str(),0777);
		   string file_name = fn + "/" + edges[idx]->name + ".txt";
		   out_file = fopen(file_name.c_str(),"w");
		   for(unsigned int j=0; j<edges[idx]->volume_flow_rate.size(); j++)
		   {
		      double t = time[j];
		      double vfr = edges[idx]->volume_flow_rate[j];
		      fprintf(out_file, "%9.7e, %9.7e\n",t,vfr);
		   }
		   fclose(out_file);
   	}
   }

      if(do_myogenic)
   {
	   mkdir(fn.c_str(),0777);

		string file_name = fn + "/p_ave.txt";
		p_ave->save_results(file_name);

   }

   //if(do_lum_RBC_transport)
   //{
	//   mkdir(fn.c_str(),0777);
	//   RBClum->save_results(fn, time);
   //}

}

//--------------------------------------------------------------
void solver_lumped::save_results(double dt, string folder_name)
{
	vector<string> edge_list, node_list;
	for(int i=0; i<number_of_edges; i++)
	{
		edge_list.push_back(edges[i]->name);
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		node_list.push_back(nodes[i]->name);
	}
	save_results(dt,folder_name,edge_list,node_list);
}

//--------------------------------------------------------------
void solver_lumped::save_results(double dt, string folder_name, vector<string> edge_list, vector<string> node_list)
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);
   
   string fn = "results/" + folder_name + "/" + name;

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
   	if(nodes[idx]->do_save_memory)
   	{
   		mkdir(fn.c_str(),0777);
	      string file_name = fn + "/" + nodes[idx]->name + ".txt";
	      out_file = fopen(file_name.c_str(),"w");

	      int j=0;
			double ts=0.;
			double t_end=time.back();
			while(ts<t_end && j<time.size()-1)
			{
				if(time[j]<=ts && ts<time[j+1])
				{	
					double a0 = (time[j+1]-ts)/(time[j+1]-time[j]);
					double a1 = (ts-time[j])/(time[j+1]-time[j]);

			   	double p = nodes[idx]->pressure[j]*a0 + nodes[idx]->pressure[j+1]*a1;
			      fprintf(out_file, "%9.7e, %9.7e\n", ts, p);
	         	ts += dt;
				}
				else
				{
					j++;
				}
			}
	      fclose(out_file);
   	}
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
   	if(edges[idx]->do_save_memory)
   	{
   		mkdir(fn.c_str(),0777);
	      string file_name = fn + "/" + edges[idx]->name + ".txt";
	      out_file = fopen(file_name.c_str(),"w");

	      int j=0;
			double ts=0.;
			double t_end=time.back();
			while(ts<t_end && j<time.size()-1)
			{
				if(time[j]<=ts && ts<time[j+1])
				{
		         double a0 = (time[j+1]-ts)/(time[j+1]-time[j]);
					double a1 = (ts-time[j])/(time[j+1]-time[j]);

		         double vfr = edges[idx]->volume_flow_rate[j]*a0 + edges[idx]->volume_flow_rate[j+1]*a1;
		         fprintf(out_file, "%9.7e, %9.7e\n",ts,vfr);
	         	ts += dt;
				}
				else
				{
					j++;
				}
			}

	      fclose(out_file);
   	}
   }

    if(do_myogenic)
   {
	   mkdir(fn.c_str(),0777);
	   
		string file_name = fn + "/p_ave.txt";
		p_ave->save_results(dt, file_name);

   }
}

//--------------------------------------------------------------
void solver_lumped::save_model(string model_name, string folder_name)
{
   FILE *out_file;
   string file_name = folder_name + "/" + model_name + "/" + name + ".csv";
	out_file = fopen(file_name.c_str(),"w");

	// writing edges data to file
	fprintf(out_file,"data of edges\n");
	fprintf(out_file,"type, name, node start, node end, initial condition [SI], parameter [SI]\n");
	for(int i=0; i<number_of_edges; i++)
	{
		fprintf(out_file, "%s, %s, %s, %s, %6.3e",edges[i]->type.c_str(),edges[i]->name.c_str(),edges[i]->node_name_start.c_str(), edges[i]->node_name_end.c_str(), edges[i]->volume_flow_rate_initial);
		for(int j=0; j<edges[i]->parameter.size(); j++)
		{
			fprintf(out_file, ", %6.3e", edges[i]->parameter[j]);
		}
		fprintf(out_file,"\n");
	}
	fprintf(out_file,"\n");

	// writing nodes data to file
	fprintf(out_file,"data of nodes\n");
	fprintf(out_file,"type, name, initial condition [SI]\n");
	for(int i=0; i<number_of_nodes; i++)
	{
		fprintf(out_file, "%s, %s, %6.3e\n",nodes[i]->type.c_str(),nodes[i]->name.c_str(),nodes[i]->pressure_initial);
	}
	fprintf(out_file,"\n");
   fclose(out_file);
}

//--------------------------------------------------------------
void solver_lumped::save_initials(string model_name, string folder_name)
{
   string file_name = folder_name + "/" + model_name + "/init/" + name + ".csv";

	FILE* out_file;
	out_file = fopen(file_name.c_str(),"w");

	for(int i=0; i<number_of_edges; i++)
	{	
		double vfr = edges[i]->vfr/1.e6; // convert ml/s to m3/s
		fprintf(out_file, "edge, %s, %9.7e\n",edges[i]->name.c_str(),vfr);
	}

	for(int i=0; i<number_of_nodes; i++)
	{
		double p = nodes[i]->p*mmHg_to_Pa; // convert mmHg to Pa
		fprintf(out_file, "node, %s, %9.7e\n",nodes[i]->name.c_str(),p);
	}

   fclose(out_file);
}

//--------------------------------------------------------------
void solver_lumped::load_initials()
{
	ifstream file_in;
	string file_name = input_folder_path + "/init/" + name + ".csv";
	file_in.open(file_name);
	string line;
	if(file_in.is_open())
	{
		double p_conv = 1./(mmHg_to_Pa*elastance(0.));
		while(getline(file_in,line))
		{
			// clearing spaces and \n
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());
			
			// seperating the strings by comma
			vector<string> sv = separate_line(line);

			if(sv.size()>0)
			{
				if(sv[0] == "edge")
				{
					int idx = edge_id_to_index(sv[1]);
					double vfr = stod(sv[2],0); // m3/s
					edges[idx]->volume_flow_rate.clear();
					edges[idx]->volume_flow_rate.push_back(vfr);
					edges[idx]->vfr = vfr*1.e6; // to ml/s
				}
				else if(sv[0] == "node")
				{
					int idx = node_id_to_index(sv[1]);
					double p = stod(sv[2],0); // Pa
					nodes[idx]->pressure.clear();
					nodes[idx]->pressure.push_back(p);
					nodes[idx]->p = p/mmHg_to_Pa; // to mmHg
					nodes[idx]->y = p*p_conv;
				}
			}
		}
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_initials() function!!! file: " << file_name << "\nExiting..." << endl;
		exit(-1);
	}
	file_in.close();
}



bool solver_lumped::saveVectorToFile(const vector<double>& vec, const std::string& filePath) {
    // Create directories if they don't exist
    std::filesystem::path path(filePath);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream outFile(filePath);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file: " << filePath << std::endl;
        return false;
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i < vec.size() - 1) {
            outFile << "\n";
        }
    }

    outFile.close();
    //std::cout << "Vector saved to: " << filePath << std::endl;
    return true;
}

bool solver_lumped::saveVectorsToFile(const std::vector<double>& time, const std::vector<double>& vec, const std::string& filePath)
{
    // Create directories if they don't exist
    std::filesystem::path path(filePath);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream outFile(filePath);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file: " << filePath << std::endl;
        return false;
    }

    for (size_t i = 0; i < time.size(); ++i) {
        outFile << time[i] << ", " << vec[i] << "\n";
    }

    outFile.close();

   // std::cout << "Vectors saved to: " << filePath << std::endl;
    return true;
}