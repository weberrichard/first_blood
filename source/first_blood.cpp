#include "first_blood.h"

//--------------------------------------------------------------
first_blood::first_blood(string folder_name)
{
	// saving the name of the folder
	input_folder_path = folder_name;

	// getting rid of global path
	case_name = input_folder_path.substr(input_folder_path.rfind('/')+1);

	// loading all the input data from csv files
	load_ok = load_model();

	if(load_ok == true)
	{
		// setting constants in every model
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->set_constants(gravity, density, kinematic_viscosity, mmHg_to_Pa, atmospheric_pressure,poisson_coefficient, courant_number);
		}

		for(int i=0; i<number_of_lum; i++)
		{
			lum[i]->set_constants(gravity, density, kinematic_viscosity, mmHg_to_Pa, atmospheric_pressure);
		}
		// converting moc t-p to SI
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->convert_time_series();
		}
	}
}

//--------------------------------------------------------------
first_blood::~first_blood()
{
	initialization();
}

//--------------------------------------------------------------
bool first_blood::load_model()
{
	// loading the main csv
	bool load_ok = load_main_csv();
	if(load_ok == false)
	{
		return load_ok;
	}

	// setting the number of models
	number_of_moc = moc.size();
	number_of_lum = lum.size();
	number_of_nodes = nodes.size();

	// loading the csv for moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->load_model();
	}

	// loading the csv for lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->load_model();
	}
	return load_ok;
}

//--------------------------------------------------------------
bool first_blood::load_main_csv()
{
	load_ok = false;
	ifstream file_in;
	string file_path = input_folder_path + "/main.csv";
	file_in.open(file_path);
	string line;
	if(file_in.is_open())
	{
		int nm=0,nl=0,nn=0; // nm for moc, nl for lumped, nn for nodes
		while(getline(file_in,line))
		{
			// clearing spaces
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());

			// separating strings to vector by comma
			vector<string> sv = separate_line(line);

			if(sv[0] == "run") // run type, boundaries if given
			{
				run_type = sv[1];
				if(sv.size()>2)
				{
					upstream_boundary.node = sv[2];
					upstream_boundary.file_name = sv[3];
				}
				else
				{
					upstream_boundary.node = "";
					upstream_boundary.file_name = "";
				}
			}
			else if(sv[0] == "time") // time
			{
				time_end = stod(sv[1],0);
			}
			else if(sv[0] == "material") // material model
			{
				if(sv[1] == "linear")
				{
					material_type = 0;
				}
				else if(sv[1] == "Olufsen" || sv[1] == "olufsen")
				{
					material_type = 1;
				}
			}
			else if(sv[0] == "solver") // setting solver type, 0: maccormack, 1: moc
			{
				if(sv[1] == "maccormack" || sv[1] == "Maccormack" || sv[1] == "MacCormack")
				{
					solver_type = 0;
				}
				else if(sv[1] == "moc")
				{
					solver_type = 1;
				}
				else
				{
					cout << "Solver type: " << sv[1] << " is not known. Using MacCormack, continoiug..." << endl;
				}
			}
			else if(sv[0] == "moc") // 1D moc model
			{
				moc.push_back(new solver_moc(sv[1],input_folder_path));
				int k=2;
				while(k<sv.size())
				{
					moc[nm]->boundary_main_node.push_back(sv[k]);
					moc[nm]->boundary_model_node.push_back(sv[k+1]);
					k+=2;
				}
				nm++;
			}
			else if(sv[0] == "lumped" || sv[0] == "lum") // 0D lumed models
			{
				lum.push_back(new solver_lumped(sv[1],input_folder_path));
				int k=2;
				while(k<sv.size())
				{
					lum[nl]->boundary_main_node.push_back(sv[k]);
					lum[nl]->boundary_model_node.push_back(sv[k+1]);
					k+=2;
				}
				nl++;
			}
			else if(sv[0] == "node") // main nodes between models
			{
				nodes.push_back(sv[1]);
				nn++;
			}
			else if(sv[0] == "RBCtransport"){ //RBC transport
				if(sv[1] == "on"){
					do_RBC_transport = true;
					RBC_node_transport = new Transport_node(RBC); 
				}
			}
			else if(sv[0] == "HBsat_transport"){
				if(sv[1] == "on"){
					do_HBsat_transport = true;
					HB_O2_node_transport = new Transport_node(HB_O2_saturation);
				}
			}
			else if(sv[0] == "PlasmaO2C_transport"){
				if(sv[1] == "on"){
					do_Plasma_O2_transport = true;
					PlasmaO2_node_transport = new Transport_node(C_Plasma_O2);
				}
			}
			else if(sv[0]=="baroreflex"){
				if(sv[1] == "on"){do_baroreflex = true;}
				time_period = stod(sv[4],0);
				sys_moc = sv[2];
				sys_edge_name = sv[3];
			}

			//CO2 stuff
			else if(sv[0] == "pla_CO2_transport"){
				if(sv[1] == "on"){
					do_pla_CO2_transport = true;
					transport_node_CO2_pla = new Transport_node(CO2_pla);
				}
			}

			else if(sv[0] == "rbc_CO2_transport"){
				if(sv[1] == "on"){
					do_rbc_CO2_transport = true;
					transport_node_CO2_rbc = new Transport_node(CO2_rbc);
				}
			}

			else if(sv[0] == "pla_HCO3_transport"){
				if(sv[1] == "on"){
					do_pla_HCO3_transport = true;
					transport_node_HCO3_pla = new Transport_node(HCO3_pla);
				}
			}

			else if(sv[0] == "rbc_HCO3_transport"){
				if(sv[1] == "on"){
					do_rbc_HCO3_transport = true;
					transport_node_HCO3_rbc = new Transport_node(HCO3_rbc);
				}
			}

			else if(sv[0] == "HbCO2_transport"){
				if(sv[1] == "on"){
					do_HbCO2_transport = true;
					transport_node_HbCO2 = new Transport_node(HbCO2);
				}
			}




		}
		load_ok = true;
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_main_csv() function!!! file: " << file_path << "\n" << endl;
		load_ok = false;
		return load_ok;
	}

	// loading p-t if it is given in main.csv
	if(upstream_boundary.node != "")
	{
		for(int i=0; i<nodes.size(); i++)
		{
			for(int j=0; j<moc.size(); j++)
			{
				for(int k=0; k<moc[j]->nodes.size(); k++)
				{
					if(nodes[i] == moc[j]->nodes[k]->name)
					{
						moc[j]->load_time_series(upstream_boundary.file_name);
					}
				}
			}
		}
	}

	file_in.close();

	return load_ok;
}

//--------------------------------------------------------------
bool first_blood::run()
{
	bool is_run_ok;

	// initialization of the whole model with mocs and lums
	initialization();

	if(init_from_file)
	{
		load_initials();
	}

	check_valves();//checking thw numer of in and outgoing edges. 1+1 is the only valid 

	connect_0D_edges();

	init_time_periods_for_lum(time_period);

	set_sys_edge_pointer();//for baroreflex


	if(run_type == "forward") // simple forward calculation
	{
		is_run_ok = true;

		if(number_of_moc>0)
		{
			// initial timesteps
			int moc_idx=0, e_idx;
			for(int i=0; i<number_of_moc; i++)
			{
				moc[i]->timesteps();
			}	

			// finding lowest new timestep
			double t_act = lowest_new_time(moc_idx, e_idx);
			double t_old = -1.e10;


			// main cycle
			while(!is_run_end(t_act,t_old) && is_run_ok)
			{
				
				//cout << "t: " << t_act << endl;
				/*if(t_act>7.) // improve this
				{
					int idx = lum_id_to_index("heart_kim");
					heart_rate = 90.;
					lum[idx]->heart_rate = heart_rate;
					time_period = 60./heart_rate;
				}

				if(t_act>13.)
				{
					do_autoregulation = true;
				}

				if(do_autoregulation)
				{
					autoregulation();
				}*/

				// solving lowest edge inner points
				if(solver_type == 0)
				{
					moc[moc_idx]->edges[e_idx]->solve_maccormack();
				}
				else if(solver_type == 1)
				{
					moc[moc_idx]->edges[e_idx]->solve_moc();
				}

				// boundaries (inner)
				moc[moc_idx]->boundaries(e_idx, t_act);

				// boundaries with 0D if exist
				int si = moc[moc_idx]->edges[e_idx]->node_index_start;
				if(moc[moc_idx]->nodes[si]->is_master_node)
				{
					int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;
					solve_lum_newton(lum_idx, t_act);
				}

				int ei = moc[moc_idx]->edges[e_idx]->node_index_end;
				if(moc[moc_idx]->nodes[ei]->is_master_node)
				{
					int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;
					solve_lum_newton(lum_idx, t_act);
				}

				// postproc: interpolate, save
				moc[moc_idx]->edges[e_idx]->update();


				//O2 transport, plasmaO2, HBsat, RBC concantration
				O2_transport( moc_idx, si, ei, e_idx, t_act);
				

				//CO2 transport
				CO2_transport( moc_idx, si, ei, e_idx, t_act);
				

				//save transport variables for lumped models
				save_transport_var_for_lum( moc_idx, si, ei);
				

				//baroreflex
				double T_act_new = time_period;
				if(moc[moc_idx]->nodes[si]->is_master_node){
					int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;
					if(lum[lum_idx]->check_whole_period(t_act)){
						//cout<<systole(sys_edge->pressure_start, sys_edge->time, lum[lum_idx]->T_sum)<<endl;
						//cout<<period_of_first_lum<<"  "<<lum[lum_idx]->period<<endl;
						if(do_baroreflex && period_of_first_lum == lum[lum_idx]->period && period_of_first_lum > 7){
							double sys = systole(sys_edge->pressure_start, sys_edge->time, lum[lum_idx]->T_sum) - atmospheric_pressure;
							T_act_new = baroreflex(sys);
							//cout<<T_act_new<<endl;
						}
						lum[lum_idx]->update_period_time(T_act_new);
					}
				}

				if(moc[moc_idx]->nodes[ei]->is_master_node){
					int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;
					if(lum[lum_idx]->check_whole_period(t_act)){
						//cout<<systole(sys_edge->pressure_start, sys_edge->time, lum[lum_idx]->T_sum)<<endl;
						//cout<<period_of_first_lum<<"  "<<lum[lum_idx]->period<<endl;
						if(do_baroreflex && period_of_first_lum == lum[lum_idx]->period && period_of_first_lum > 7){
							double sys = systole(sys_edge->pressure_start, sys_edge->time, lum[lum_idx]->T_sum) - atmospheric_pressure;
							T_act_new = baroreflex(sys);
							//cout<<T_act_new<<endl;
						}
						lum[lum_idx]->update_period_time(T_act_new);
					}
				}

				// update the period
				if(T_sum + time_period - t_act<0.){period++; T_sum += time_period; time_period = T_act_new;}
				//tracking period of the "furthest" lumped model
				if(moc[moc_idx]->nodes[ei]->is_master_node || moc[moc_idx]->nodes[si]->is_master_node){period_of_first_lum = period;}


				//moc[moc_idx]->edges[e_idx]->update_variables();
				if(moc[moc_idx]->edges[e_idx]->do_save_memory)
				{
					moc[moc_idx]->edges[e_idx]->save_field_variables();
				}

				// get the time average values, e.g MAP
				//if(moc[moc_idx]->edges[e_idx]->ID == "A1")
				//{
				//	calculate_time_average();
				//}

				// new timestep
				is_run_ok = moc[moc_idx]->edges[e_idx]->new_timestep();

				// find new lowest timestep
				t_old = t_act;
				t_act = lowest_new_time(moc_idx, e_idx);

				// decreasing time_counter if it passed t_end
				if(t_act >= time_end)
				{
					time_counter--;
				}

			}
		}
		else // only LUMPED MODEL without any moc
		{
			if(number_of_lum==0)
			{
				cout << " There is no LUMPED model (nor MOC).\n Exiting..." << endl;
				return false;
			}

			double t_act = 0.;
			while(t_act<time_end)
			{
				t_act = lum[0]->time.back() + dt_lumped;
				solve_lum_newton(0, t_act);
			}
		}
	}

	return is_run_ok;
}

//--------------------------------------------------------------
double first_blood::lowest_new_time(int &moc_idx, int &e_idx)
{
	double t_act=1.e10;
	int idx;
	for(int i=0; i<number_of_moc; i++)
	{
		double t = moc[i]->min_time(idx);
		if(t<t_act)
		{
			t_act = t;
			moc_idx = i;
			e_idx = idx;
		}
	}

	return t_act;
}

//--------------------------------------------------------------
void first_blood::solve_lum_newton(int index, double t_act)
{
	/*
	x = [q1,q2,...qm,p1,p2,p3,...pn,y1,y2,...ye,qmoc1,...qmock] y: for elastance if present
	f = [edge1,edge2,...edgem,node1,node2,...node,elas1,elas2,...elase,char1,char2,...]
	*/
	int m = lum[index]->number_of_edges;
	int n = lum[index]->number_of_nodes;
	int l = lum[index]->number_of_elastance;
	int k = lum[index]->boundary_indices.size();

	// setting initial conditions for lumped part
	lum[index]->initialization_newton(t_act);
	// setting initial conditions for moc part
	for(int j=0; j<lum[index]->boundary_indices.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0];
		int moc_edge_index = lum[index]->boundary_indices[j][1];
		int edge_end = lum[index]->boundary_indices[j][2];
		int N = m+n+2*l+j;
		moc[moc_index]->initialization_newton(lum[index]->x,N,moc_edge_index,edge_end);
	}

	// start of newton iteration
	int i=0;
	do
	{
		lum[index]->coefficients_newton(t_act);

		for(int j=0; j<lum[index]->boundary_indices.size(); j++)
		{
			int moc_index = lum[index]->boundary_indices[j][0];
			int moc_edge_index = lum[index]->boundary_indices[j][1];
			int edge_end = lum[index]->boundary_indices[j][2];
			int lum_node_index = lum[index]->boundary_indices[j][3];
			double q = lum[index]->x(n+m+2*l+j);
			// double A = lum[index]->x(n+m+2*l+2*j+1);
			double p = lum[index]->x(m+lum_node_index)*mmHg_to_Pa;
			vector<double> v; // f_char,dchar_dp,dchard_dq
			if(edge_end==1)
			{
				v = moc[moc_index]->edges[moc_edge_index]->boundary_newton_end(q,p,t_act);
			}
			else
			{
				v = moc[moc_index]->edges[moc_edge_index]->boundary_newton_start(q,p,t_act);
			}

			// node continouity equation and jacobian
			lum[index]->f(m+lum_node_index) += edge_end*lum[index]->x(m+n+2*l+j)*1.e6; // sign(q)*q*1e6 [ml/s]
			lum[index]->Jac(m+lum_node_index,m+n+2*l+j) = edge_end*1.e6; // for node continouity // sign(q)*1e6 [1]

			// characteristic equation
			// also converting to non_SI for numerical stability
			lum[index]->f(m+n+2*l+j) = v[0];
			lum[index]->Jac(m+n+2*l+j,m+lum_node_index) = v[1]*mmHg_to_Pa; // dp
			lum[index]->Jac(m+n+2*l+j,m+n+2*l+j) = v[2]; // dQ
			// lum[index]->Jac(m+n+2*l+2*j,m+n+2*l+2*j+1) = -q/A/A; // dA

			// cross section - pressure equation
			// lum[index]->f(m+n+2*l+2*j+1) = v[1];
			// lum[index]->Jac(m+n+2*l+2*j+1,m+lum_node_index) = v[3]*mmHg_to_Pa; // dp
			// lum[index]->Jac(m+n+2*l+2*j+1,m+n+2*l+2*j) = v[5]; // dQ
			// lum[index]->Jac(m+n+2*l+2*j+1,m+n+2*l+2*j+1) = 1.; // dA
		}

		// cout << endl << i << " ITERATION " << i << endl;
		// cout << endl << lum[index]->Jac << endl;
		// cout << "x" << endl << lum[index]->x << endl;
		// cout << "f" << endl << lum[index]->f << endl;

		// actually solving the Newton's technique
		VectorXd dx = lum[index]->Jac.colPivHouseholderQr().solve(-lum[index]->f);
		lum[index]->x += dx;

		//cout << " ITERATION END " << endl;
		//cout << endl;
		//cin.get();
		
		i++;
	}
	while(lum[index]->f.norm() > 1e-5 && i<100);

	if(i>=100)
	{
		cout << "\n !!! ERROR !!! Newton's technique did NOT converge at Lum model " << lum[index]->name << endl;
		cout << endl << i << " ITERATION " << i << endl;
		cout << lum[index]->Jac << endl;
		cout << "x" << endl << lum[index]->x << endl;
		cout << "f" << endl << lum[index]->f << endl;
		exit(-1);
	}

	// substitute the results back to 0D
	lum[index]->substitute_newton(t_act);

	// substitute the results back to 1D
	for(int j=0; j<lum[index]->boundary_indices.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0];
		int moc_edge_index = lum[index]->boundary_indices[j][1];
		int edge_end = lum[index]->boundary_indices[j][2];
		int lum_node_index = lum[index]->boundary_indices[j][3];
		double p = lum[index]->x(m+lum_node_index)*mmHg_to_Pa;
		double q = lum[index]->x(m+n+2*l+j);
		moc[moc_index]->substitute_newton(moc_edge_index,edge_end,t_act,p,q);
	}
}

//--------------------------------------------------------------
/*
void first_blood::calculate_time_average()
{
	// mean arterial pressure
	string id = "A1";
	int idx = moc[0]->edge_id_to_index(id);
	double tn = moc[0]->edges[idx]->time.back();
	double vn = moc[0]->edges[idx]->pressure_start.back();
	map->update(tn,vn,time_period);

}*/

//--------------------------------------------------------------
void first_blood::autoregulation()
{	
	//double S = 1.;
	//double pset = 90.;
	//pset = pset*mmHg_to_Pa + atmospheric_pressure;
	//double factor = 1.0 + S*(map->average.back()-pset);
	//cout << "f: " << factor << endl;

	double pset = 90.;
	double factor = (map->average.back()-atmospheric_pressure)/mmHg_to_Pa / pset;

	vector<string> perif_brain{"p25","p26","p27","p28","p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p41","p42","p43","p44","p45","p46"};

	for(int i=0; i<perif_brain.size(); i++)
      {
         string id = perif_brain[i];
         int idx = lum_id_to_index(id);
         lum[idx]->edges[0]->parameter_factor = factor;
         lum[idx]->edges[1]->parameter_factor = factor;
         lum[idx]->edges[2]->parameter_factor = 1./factor;
      }
}

//--------------------------------------------------------------
void first_blood::save_time_average(string folder_name)
{
	mkdir("results",0777);
	mkdir(("results/"+folder_name).c_str(),0777);

	string file_name = folder_name + "/arterial/map.txt";
	map->save_results(file_name);

}

//--------------------------------------------------------------
void first_blood::save_time_average(double dt, string folder_name)
{
	mkdir("results",0777);
	mkdir(("results/"+folder_name).c_str(),0777);

	string file_name = folder_name + "/arterial/map.txt";
	map->save_results(dt, file_name);

}

//--------------------------------------------------------------
bool first_blood::is_run_end(double t_act, double t_old)
{
	if(is_periodic_run)
	{
		if(t_act<time_end_min)
		{
			return false;
		}
		else if(t_act>time_end_max)
		{
			time_end = time_end_max;
			return true;
		}
		else
		{
			int n=0;
			double t = t_act;
			while(t >= time_period)
			{
				t -= time_period;
				n++;
			}

 			// checking the end of a cycle
 			if((double)n*time_period>=t_old && (double)n*time_period<t_act)
 			{
				int idx = moc[0]->node_id_to_index(time_node);
				if(idx<0)
				{
					time_node = moc[0]->nodes[0]->name;
					idx = 0;
				}
				if(time_var == "P")
				{
					double val = (systole(moc[0]->nodes[idx]->pressure, moc[0]->nodes[idx]->time, (n-1.)*time_period)-atmospheric_pressure)/mmHg_to_Pa;
					if(abs((val-time_val_old)/time_val_old)<.0001)
					{
						time_end = t_act;
						return true;
					}
					else
					{
						time_val_old = val;
						return false;
					}
				}
				else
				{
					cout << "\n time_var: " << time_var << " is not valid, avaialble: P" << endl;
					time_end = t_act;
					return true;
				}
 			}
 			else
 			{
 				return false;
 			}
		}
	}
	else
	{
		if(time_counter>0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
}

//--------------------------------------------------------------
void first_blood::initialization()
{
	// fb class
	number_of_moc = moc.size();
	number_of_lum = lum.size();
	number_of_nodes = nodes.size();
	time_counter = 0;

	// setting initial conditions
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->initialization(pressure_initial,material_type);
		time_counter += moc[i]->number_of_edges;
	}
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->initialization(time_period);
		time_counter += 1;
	}

	// setting master nodes
	build_master();

	// setting number of moc in lum
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->set_newton_size();
	}

	// making sure time_node is saved
	if(is_periodic_run)
	{
		int idx = moc[0]->node_id_to_index(time_node);
		if(idx<0)
		{
			time_node = moc[0]->nodes[0]->name;
			idx = 0;
		}
		if(moc.size()>0)
		{
			vector<string> el, nl{time_node};
		   set_save_memory(moc[0]->name,"moc",el,nl);
		}
	}

	// saving these for time average vectors
	// vector<string> el{"A1","A5","A6","A15","A20"}, nl;
	// set_save_memory(moc[0]->name,"moc",el,nl);

	// time average stuff
	map = new time_average();

	//O2 transfer class
	RBC1D = new Transport_1D(TRBCType);
	HBsat1D = new Transport_1D(HB_O2_saturation);
	Plasma_O21D = new Transport_1D(C_Plasma_O2);

	//CO2 stuff
	CO2_pla_1D = new Transport_1D(CO2_pla);
    CO2_rbc_1D = new Transport_1D(CO2_rbc);
    HCO3_pla_1D = new Transport_1D(HCO3_pla);
    HCO3_rbc_1D = new Transport_1D(HCO3_rbc);
    HbCO2_1D = new Transport_1D(HbCO2);

}

//--------------------------------------------------------------
void first_blood::build_master()
{
	for(int i=0; i<nodes.size(); i++)
	{
		for(int j=0; j<number_of_lum; j++)
		{
			for(int k=0; k<lum[j]->boundary_main_node.size(); k++)
			{
				if(nodes[i] == lum[j]->boundary_main_node[k])
				{
					int kk = lum[j]->node_id_to_index(lum[j]->boundary_model_node[k]);
					lum[j]->nodes[kk]->is_master_node = true;

					for(int l=0; l<number_of_moc; l++)
					{
						for(int m=0; m<moc[l]->boundary_main_node.size(); m++)
						{
							if(nodes[i] == moc[l]->boundary_main_node[m])
							{
								int mm = moc[l]->node_id_to_index(moc[l]->boundary_model_node[m]);
								moc[l]->nodes[mm]->is_master_node = true;
								moc[l]->nodes[mm]->master_node_lum = j;

								for(int n=0; n<moc[l]->nodes[mm]->edge_in.size(); n++)
								{
									vector<int> v{l,moc[l]->nodes[mm]->edge_in[n],1,kk};
									// index of moc, index of edge in moc, start(-1)/end(+1), index of node in lumped
									lum[j]->boundary_indices.push_back(v);
								}
								for(int n=0; n<moc[l]->nodes[mm]->edge_out.size(); n++)
								{
									vector<int> v{l,moc[l]->nodes[mm]->edge_out[n],-1,kk};
									// index of moc, index of edge in moc, start(-1)/end(+1), index of node in lumped
									lum[j]->boundary_indices.push_back(v);
								}
							}
						}
					}
				}
			}
		}
	}
}

//--------------------------------------------------------------
void first_blood::save_results()
{
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

   string folder_name = case_name;

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(folder_name);
	}

	// saving time averages, e.g. map, cfr
	// save_time_average("results/" + folder_name);
}

//--------------------------------------------------------------
void first_blood::save_results(string folder_name)
{
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(folder_name);
	}

	// saving time averages, e.g. map, cfr
	// save_time_average("results/" + folder_name);
}

//--------------------------------------------------------------
void first_blood::save_results(double dt)
{
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

   string folder_name = case_name;

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(dt, folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(dt, folder_name);
	}

	// saving time averages, e.g. map, cfr
	save_time_average(dt, "results/" + folder_name);
}

//--------------------------------------------------------------
void first_blood::save_results(double dt, string folder_name)
{
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(dt, folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(dt, folder_name);
	}
}

//--------------------------------------------------------------
void first_blood::save_results(string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	if(model_type == "moc")
	{
		for(int i=0; i<moc.size(); i++)
		{
			if(model_name == moc[i]->name)
			{
				moc[i]->save_results(folder_name,edge_list,node_list);
			}
		}
	}
	else if(model_type == "lum")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name)
			{
				lum[i]->save_results(folder_name,edge_list,node_list);			
			}
		}
	}
	else if(model_type == "RBC_transport")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name && do_RBC_transport && lum[i]->RBClum->do_save_results)
			{
				lum[i]->RBClum->save_results(folder_name, lum[i]->time, model_name);		
			}
		}
	}
	else if(model_type == "HBsat_transport")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name && do_HBsat_transport && lum[i]->HBsatlum->do_save_results)
			{
				lum[i]->HBsatlum->save_results(folder_name, lum[i]->time, model_name);		
			}
		}
	}
	else if(model_type == "PlasmaO2_transport")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name && do_Plasma_O2_transport && lum[i]->PlasmaO2lum->do_save_results)
			{
				lum[i]->PlasmaO2lum->save_results(folder_name, lum[i]->time, model_name);		
			}
		}
	}

	else if(model_type == "CO2_transport++")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name )
			{
				lum[i]->CO2_pla_lum->save_results(folder_name, lum[i]->time, model_name);
				lum[i]->CO2_rbc_lum->save_results(folder_name, lum[i]->time, model_name);
				lum[i]->HCO3_pla_lum->save_results(folder_name, lum[i]->time, model_name);
				lum[i]->HCO3_rbc_lum->save_results(folder_name, lum[i]->time, model_name);
				lum[i]->HbCO2_lum->save_results(folder_name, lum[i]->time, model_name);		
			}
		}
	}
}

//--------------------------------------------------------------
void first_blood::save_results(double dt, string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
	mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	if(model_type == "moc")
	{
		for(int i=0; i<moc.size(); i++)
		{
			if(model_name == moc[i]->name)
			{
				moc[i]->save_results(dt,folder_name,edge_list,node_list);
			}
		}
	}
	else if(model_type == "lum")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name)
			{
				lum[i]->save_results(dt,folder_name,edge_list,node_list);			
			}
		}
	}
}

//--------------------------------------------------------------
void first_blood::save_model(string model_name)
{
	string folder_name = "../../models/";

	save_model(model_name, folder_name);
}

//--------------------------------------------------------------
void first_blood::save_model(string model_name, string folder_name)
{
   mkdir(folder_name.c_str(),0777);
   mkdir((folder_name+"/"+model_name).c_str(),0777);

	// saving moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_model(model_name, folder_name);
	}

	// saving lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_model(model_name, folder_name);
	}

	// saving main model TODO
	FILE *out_file;
	string file_name = folder_name + "/" + model_name + "/main.csv";
	out_file = fopen(file_name.c_str(),"w");

	fprintf(out_file, "run,forward\n");
	fprintf(out_file, "time,%6.3f\n",time_end);
	fprintf(out_file, "\n");

	fprintf(out_file, "type,name,main node,model node,main node,model node,...\n");
	for(int i=0; i<number_of_moc; i++)
	{
		fprintf(out_file, "moc,%s",moc[i]->name.c_str());
		for(int j=0; j<moc[i]->boundary_main_node.size(); j++)
		{
			fprintf(out_file, ",%s,%s",moc[i]->boundary_main_node[j].c_str(),moc[i]->boundary_model_node[j].c_str());
		}
		fprintf(out_file, "\n");
	}
	fprintf(out_file, "\n");

	for(int i=0; i<number_of_lum; i++)
	{
		fprintf(out_file, "lumped,%s",lum[i]->name.c_str());
		for(int j=0; j<lum[i]->boundary_main_node.size(); j++)
		{
			fprintf(out_file, ",%s,%s",lum[i]->boundary_main_node[j].c_str(),lum[i]->boundary_model_node[j].c_str());
		}
		fprintf(out_file, "\n");
	}
	fprintf(out_file, "\n");

	for(int i=0; i<number_of_nodes; i++)
	{
		fprintf(out_file, "node,%s\n",nodes[i].c_str());
	}
	fprintf(out_file, "\n");
   fclose(out_file);
}

//--------------------------------------------------------------
void first_blood::clear_save_memory()
{
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->clear_save_memory();
	}
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->clear_save_memory();
	}
}

//--------------------------------------------------------------
void first_blood::set_save_memory(string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
	if(model_type == "moc")
	{
		for(int i=0; i<number_of_moc; i++)
		{
			if(moc[i]->name == model_name)
			{
				moc[i]->set_save_memory(edge_list, node_list);
			}
		}
	}
	if(model_type == "lum" || model_type == "lumped")
	{
		for(int i=0; i<number_of_lum; i++)
		{
			if(lum[i]->name == model_name)
			{
				lum[i]->set_save_memory(edge_list, node_list);
			}
		}
	}

	if(model_type == "RBC_transport" && do_RBC_transport)
	{
		for(int i=0; i<number_of_lum; i++)
		{
			if(lum[i]->name == model_name)
			{
				lum[i]->RBClum->set_save_memory();
			}
		}
	}

	if(model_type == "HBsat_transport" && do_HBsat_transport)
	{
		for(int i=0; i<number_of_lum; i++)
		{
			if(lum[i]->name == model_name)
			{
				lum[i]->HBsatlum->set_save_memory();
			}
		}
	}

	if(model_type == "PlasmaO2_transport" && do_Plasma_O2_transport)
	{
		for(int i=0; i<number_of_lum; i++)
		{
			if(lum[i]->name == model_name)
			{
				lum[i]->PlasmaO2lum->set_save_memory();
			}
		}
	}

	if(model_type == "CO2_transport++" && do_pla_CO2_transport && do_rbc_CO2_transport && do_pla_HCO3_transport && do_rbc_HCO3_transport && do_HbCO2_transport)
	{
		for(int i=0; i<number_of_lum; i++)
		{
			if(lum[i]->name == model_name)
			{
				lum[i]->CO2_pla_lum->set_save_memory();
				lum[i]->CO2_rbc_lum->set_save_memory();
				lum[i]->HCO3_pla_lum->set_save_memory();
				lum[i]->HCO3_rbc_lum->set_save_memory();
				lum[i]->HbCO2_lum->set_save_memory();
			}
		}
	}

}

//--------------------------------------------------------------
int first_blood::lum_id_to_index(string lum_id)
{
	int i=0, idx=-1;
	bool got_it=false;
	while(i<number_of_lum && !got_it)
	{
		if(lum_id.compare(lum[i]->name) == 0)
		{
			got_it = true;
			idx = i;
		}
		i++;
	}
	if(idx == -1)
	{
		cout << "\n !!!WARNING!!!\n solver_moc::lum_id_to_index function\nLum model is not existing, lum_id: " << lum_id << "\n Continouing..." << endl;
	}
	return idx;
}

//--------------------------------------------------------------
void first_blood::save_initials(string model_name, string folder_name)
{
   mkdir(folder_name.c_str(),0777);
   mkdir((folder_name+"/"+model_name).c_str(),0777);
   mkdir((folder_name+"/"+model_name+"/init").c_str(),0777);

   for(int i=0; i<number_of_moc; i++)
   {
   	moc[i]->save_initials(model_name, folder_name);
   }

   for(int i=0; i<number_of_lum; i++)
   {
   	lum[i]->save_initials(model_name, folder_name);
   }
}

//--------------------------------------------------------------
void first_blood::load_initials()
{
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->load_initials();
	}

	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->load_initials();
	}
}


//transport stuff
//-------------------------------------------
void Transport_1D::update_fi(vector<double> v, vector<double>& fi, vector<double>& fi_new, double l, double dt, double fiStart, double fiEnd){
    int n = v.size();
    double dx = l / (n - 1);
    //vector<double> fi_tmp = fi_new; //this will be the new fi

    for (int i = 1; i < n - 1; i++) {
        if (v[i] > 0 ) {
            fi_new[i] = fi[i] - v[i] * dt / dx * (fi[i] - fi[i - 1]);
        }
        else {
            fi_new[i] = fi[i] - v[i] * dt / dx * (fi[i + 1] - fi[i]);
        }
    }
 
    // downstream BC
    if (v[n - 1] > 0) {
        fi_new[n - 1] = fi[n - 1] - v[n - 1] * dt / dx * (fi[n - 1] - fi[n - 2]);
    }
    else {
        fi_new[n - 1] = fiEnd; // fom node
    }

    // upstream BC
    if (v[0] > 0) {
        fi_new[0] = fiStart;  // from node
    }
    else{
        fi_new[0] = fi[0] - v[0] * dt / dx * (fi[1] - fi[0]);
    }

    fi = fi_new;
}

//-----------------------------------------------------------------
Transport_1D::Transport_1D(TransportType TType):TType(TType){};


//-----------------------------------------------------------------
void Transport_1D::prescribe_node_fi_CO2(TransportType TType, double& finode){
	switch(TType){

	case CO2_pla:
		finode = 0.02635920*1.0; //m3 O2/ m3 pla
		break;

	case CO2_rbc:
		finode = 0.02986472*1.0; //m3 O2/ m3 rbc cytoplasm
		break;

	case HCO3_pla:
		finode = 0.73471517*1.0; //m3 O2/ m3 pla
		break;

	case HCO3_rbc:
		finode = 0.08811732*1.0; //m3 O2/ m3 rbc cytoplasm
		break;

	case HbCO2:
		finode = 0.04748449*1.0; //m3 O2/ m3 rbc cytoplasm
		break;
	}

}



//transport stuff for 1D nodes
//separate class is needed for each type of transport. eg.: RBC...
//--------------------------------------------
Transport_node::Transport_node(TransportType TType) : TType(TType) {};

//--------------------------------------------
void Transport_node::update_fi(double& fiNode, moc_node* node, const vector<moc_edge*>& edges){
    if (node->is_master_node == true) {
        return;
        };

    if (node->upstream_boundary > -1) {} // no idea about this, TODO

    else if (node->type_code == 1) {} //periphery

    else if (node->type_code == 0) { // junction

        int n1 = node->edge_in.size();
        int n2 = node->edge_out.size();

        double q_sum = 0.; //vfr sum of the incoming edges only
        double q;
        double fiNodeOld = fiNode;
        fiNode = 0.;

       for (int j = 0; j < n1; j++)
       {
           if (edges[node->edge_in[j]]->get_velocity().back() > 0.)
           {
           	q = edges[node->edge_in[j]]->get_velocity().back() * edges[node->edge_in[j]]->get_area().back();
            q_sum += q;

            switch(TType){
            case HB_O2_saturation:
            	fiNode += q * edges[node->edge_in[j]]->HBsat_edge.back();
            	break;
            	
            case C_Plasma_O2:
                fiNode += q * edges[node->edge_in[j]]->PlasmaO2_edge.back();
            	break;

            case RBC:
            	fiNode += q * edges[node->edge_in[j]]->RBC_edge_fi.back();
            	break;

            //CO2 stuff
            case CO2_pla:
            	fiNode += q * edges[node->edge_in[j]]->CO2_pla_edge_fi.back();
            	break;

            case CO2_rbc:
            	fiNode += q * edges[node->edge_in[j]]->CO2_rbc_edge_fi.back();
            	break;

            case HCO3_pla:
            	fiNode += q * edges[node->edge_in[j]]->HCO3_pla_edge_fi.back();
            	break;

            case HCO3_rbc:
            	fiNode += q * edges[node->edge_in[j]]->HCO3_rbc_edge_fi.back();
            	break;

            case HbCO2:
            	fiNode += q * edges[node->edge_in[j]]->HbCO2_edge_fi.back();
            	break;
            }


           }
        }
        for (int j = 0; j < n2; j++)
        {
            if (edges[node->edge_out[j]]->get_velocity()[0] < 0.)
            {
            q = edges[node->edge_out[j]]->get_velocity()[0] * edges[node->edge_out[j]]->get_area()[0];
            q_sum -= q; //not sure about the sign tho...

            switch(TType){
            case HB_O2_saturation:
            	fiNode -= q * edges[node->edge_out[j]]->HBsat_edge[0];
                break;

            case C_Plasma_O2:
            	fiNode -= q * edges[node->edge_out[j]]->PlasmaO2_edge[0];
                break;

            case RBC:
                fiNode -= q * edges[node->edge_out[j]]->RBC_edge_fi[0];
                break;


            //CO2 stuff
            case CO2_pla:
            	fiNode -= q * edges[node->edge_out[j]]->CO2_pla_edge_fi[0];
            	break;

            case CO2_rbc:
            	fiNode -= q * edges[node->edge_out[j]]->CO2_rbc_edge_fi[0];
            	break;

            case HCO3_pla:
            	fiNode -= q * edges[node->edge_out[j]]->HCO3_pla_edge_fi[0];
            	break;

            case HCO3_rbc:
            	fiNode -= q * edges[node->edge_out[j]]->HCO3_rbc_edge_fi[0];
            	break;

            case HbCO2:
            	fiNode -= q * edges[node->edge_out[j]]->HbCO2_edge_fi[0];
            	break;
            }
            }
        }
        if (q_sum != 0.) {
            fiNode /= q_sum;
        }
        else { fiNode = fiNodeOld; }//if nothing flows in it stays the old
    }
};

//----------------------------------------------------------------

void Transport_node::update_master_fi(double& fiNode, moc_node* node, const vector<moc_edge*>& edges, solver_lumped& lum_mod, solver_moc& moc_mod){
    int n1 = node->edge_in.size();
    int n2 = node->edge_out.size();


    double fiNodeOld = fiNode;
    fiNode = 0.;

    double q_sum = 0.;
    double q;

//cout<<lum_mod.name<<"  "<<moc_mod.name<<endl;


    //outgoing edges
    for (int j = 0; j < n2; j++)
    {
        if (edges[node->edge_out[j]]->get_velocity()[0] < 0.){
   	       q = edges[node->edge_out[j]]->get_velocity()[0] * edges[node->edge_out[j]]->get_area()[0];
           q_sum -= q;

            switch(TType){
            	case HB_O2_saturation:
            	fiNode -= q * edges[node->edge_out[j]]->HBsat_edge[0];
                break;

            	case C_Plasma_O2:
            	fiNode -= q * edges[node->edge_out[j]]->PlasmaO2_edge[0];
                break;

                case RBC:
                fiNode -= q * edges[node->edge_out[j]]->RBC_edge_fi[0];
                break;

                //CO2
                case CO2_pla:
                fiNode -= q * edges[node->edge_out[j]]->CO2_pla_edge_fi[0];
                break;

                case CO2_rbc:
                fiNode -= q * edges[node->edge_out[j]]->CO2_rbc_edge_fi[0];
                break;

                case HCO3_pla:
                fiNode -= q * edges[node->edge_out[j]]->HCO3_pla_edge_fi[0];
                break;

                case HCO3_rbc:
                fiNode -= q * edges[node->edge_out[j]]->HCO3_rbc_edge_fi[0];
                break;

                case HbCO2:
                fiNode -= q * edges[node->edge_out[j]]->HbCO2_edge_fi[0];
                break;
            }


        }
    } 

    double r;
    for (int j = 0; j < n1; j++)
    {
        if (edges[node->edge_in[j]]->get_velocity().back() > 0.){
   	        q = edges[node->edge_in[j]]->get_velocity().back() * edges[node->edge_in[j]]->get_area().back();
            q_sum += q;

            switch(TType){
            	case HB_O2_saturation:
            	fiNode += q * edges[node->edge_in[j]]->HBsat_edge.back();
                break;

            	case C_Plasma_O2:
            	fiNode += q * edges[node->edge_in[j]]->PlasmaO2_edge.back();
                break;

                case RBC:
                fiNode += q * edges[node->edge_in[j]]->RBC_edge_fi.back();
                break;

                //CO2
                case CO2_pla:
                fiNode += q * edges[node->edge_in[j]]->CO2_pla_edge_fi.back();
                break;

                case CO2_rbc:
                fiNode += q * edges[node->edge_in[j]]->CO2_rbc_edge_fi.back();
                break;

                case HCO3_pla:
                fiNode += q * edges[node->edge_in[j]]->HCO3_pla_edge_fi.back();
                break;

                case HCO3_rbc:
                fiNode += q * edges[node->edge_in[j]]->HCO3_rbc_edge_fi.back();
                break;

                case HbCO2:
                fiNode += q * edges[node->edge_in[j]]->HbCO2_edge_fi.back();
                break;
            }


        }
    }
    

    //0D transport part
    //finding the node corresponding to the given moc model
    int indexof_node = -1;
    string lumpMaster;
    for(int i=0; i< moc_mod.boundary_model_node.size();i++ ){
    	//cout << node->name<<endl;
    	if(node->name == moc_mod.boundary_model_node[i] ){
    		lumpMaster = moc_mod.boundary_main_node[i] ;
    	}
    }


    for(int i=0; i< lum_mod.boundary_main_node.size() ; i++){
    	if (lum_mod.boundary_main_node[i] == lumpMaster){
    		for(int j=0; j<lum_mod.nodes.size(); j++){
    			if(lum_mod.nodes[j]->name == lum_mod.boundary_model_node[i]){
    				indexof_node = j;
    			}
    		}
    	}
	}


    switch(TType){
    case HB_O2_saturation:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_HBsat.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_HBsat[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_HBsat.size(); i++){//arteriole
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_HBsat[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;

    case C_Plasma_O2:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_PlasmaO2.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_PlasmaO2[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_PlasmaO2.size(); i++){//arteriole
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_PlasmaO2[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;

    case RBC:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_RBC.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_RBC[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_RBC.size(); i++){//arteriole
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_RBC[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;


    //CO2
    case CO2_pla:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_CO2_pla.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_CO2_pla[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_CO2_pla.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_CO2_pla[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;

	case CO2_rbc:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_CO2_rbc.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_CO2_rbc[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_CO2_rbc.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_CO2_rbc[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;


    case HCO3_pla:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_HCO3_pla.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_HCO3_pla[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_HCO3_pla.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_HCO3_pla[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;



	case HCO3_rbc:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_HCO3_rbc.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_HCO3_rbc[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_HCO3_rbc.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_HCO3_rbc[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;


	case HbCO2:
    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_in_HbCO2.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_in_HbCO2[i];
        q = is->corr_edge->vfr;
        if(q > 0){
        q_sum += q;
        fiNode += q * is->fi.back();
        }

        is++;
    }

    for(int i=0; i<lum_mod.nodes[indexof_node]->D0_edges_out_HbCO2.size(); i++){
    	D0_edge* is= lum_mod.nodes[indexof_node]->D0_edges_out_HbCO2[i] ;
        q = is->corr_edge->vfr;
        if(q < 0){
        q_sum -= q;
        fiNode -= q * is->fi.back();
        }

        is++;
    }
    break;



    }


    if (q_sum != 0.) {
       fiNode /= q_sum;
    }
    else { fiNode = fiNodeOld; }//if nothing flows in it stays the old

    //node concentration must be updated as well
    switch(TType){
       case HB_O2_saturation:
       lum_mod.nodes[indexof_node]->HBsat_0Dn = fiNode;
       //cout<<lum_mod.nodes[indexof_node]->HBsat_0Dn<<"  "<<fiNode<<endl;
       break;

       case C_Plasma_O2:
       lum_mod.nodes[indexof_node]->PlasmaO2_0Dn = fiNode;
       break;

       case RBC:
       lum_mod.nodes[indexof_node]->RBC_fi0Dn = fiNode;
       break;

       //
       case CO2_pla:
       lum_mod.nodes[indexof_node]->CO2_pla_n = fiNode;
       break;

       case CO2_rbc:
       lum_mod.nodes[indexof_node]->CO2_rbc_n = fiNode;
       break;

       case HCO3_pla:
       lum_mod.nodes[indexof_node]->HCO3_pla_n = fiNode;
       break;

       case HCO3_rbc:
       lum_mod.nodes[indexof_node]->HCO3_rbc_n = fiNode;
       break;

       case HbCO2:
       lum_mod.nodes[indexof_node]->HbCO2_n = fiNode;
       break;
   }



}


//--------------------------------------------------------------
void first_blood::check_valves(){
	for(int j=0;j<moc.size();j++){
	for(int i=0;i<moc[j]->nodes.size();i++){
		if(moc[j]->nodes[i]->is_diode){
			if(moc[j]->nodes[i]->edge_in.size()!=1 || moc[j]->nodes[i]->edge_out.size() !=1){
				cout<<"Diode definition between edges is not valid. More than one incoming or more than one outgoing edge.";
				exit(-1);
			}
		}
	}
}
}


//--------------------------------------------------------------
void first_blood::connect_0D_edges(){
	if(do_Plasma_O2_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_PlasmaO2_transport){
			lum[i]->PlasmaO2lum->connect_0D_edges(*lum[i]);
		}
	}}

	if(do_RBC_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_RBC_transport){
			lum[i]->RBClum->connect_0D_edges(*lum[i]);
		}
	}}

	if(do_HBsat_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_HB_sat_transport){
			lum[i]->HBsatlum->connect_0D_edges(*lum[i]);
		}
	}}

	//CO2
	if(do_pla_CO2_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_pla_CO2_transport){
			lum[i]->CO2_pla_lum->connect_0D_edges(*lum[i]);
		}
	}}

	if(do_rbc_CO2_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_rbc_CO2_transport){
			lum[i]->CO2_rbc_lum->connect_0D_edges(*lum[i]);
		}
	}}

	if(do_pla_HCO3_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_pla_HCO3_transport){
			lum[i]->HCO3_pla_lum->connect_0D_edges(*lum[i]);
		}
	}}

	if(do_rbc_HCO3_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_rbc_HCO3_transport){
			lum[i]->HCO3_rbc_lum->connect_0D_edges(*lum[i]);
		}
	}}

	if(do_HbCO2_transport){
	for(int i=0; i< lum.size(); i++){
		if(lum[i]->do_lum_HbCO2_transport){
			lum[i]->HbCO2_lum->connect_0D_edges(*lum[i]);
		}
	}}

}


//--------------------------------------------------------------
void first_blood::set_sys_edge_pointer(){
	for(int i=0; i<moc.size();i++){
		if(moc[i]->name == sys_moc){
			for (int j=0; j<moc[i]->edges.size();j++){
				if(moc[i]->edges[j]->ID == sys_edge_name){
					sys_edge = moc[i]->edges[j];
					sys_edge->do_save_memory = true;
				}
			}
		}
	}
}


//--------------------------------------------------------------
void first_blood::init_time_periods_for_lum(double T){

	for(int i=0; i<lum.size(); i++){
		lum[i]->T_act = T;
		lum[i]->T_sum = 0.;
		lum[i]->T_last = T;
		lum[i]->heart_rate = 60./T;
		lum[i]->time_period=T;
	}
}


//--------------------------------------------------------------
double first_blood::baroreflex(double sys){
	//sys is in Pa

	double sys_mmHg = sys/mmHg_to_Pa;
	//cout<<(sys_mmHg - x_B0)<<endl;
	//cout<<sys_mmHg<<endl<<endl;

	return L_B / (1 + exp(- k_B * (sys_mmHg - x_B0))) + b_B;
}

//--------------------------------------------------------------
void first_blood::O2_transport(int moc_idx, int si, int ei, int e_idx, double t_act){

	if(do_RBC_transport){
		//RBC transport in nodes

		RBC_node_transport->update_fi(moc[moc_idx]->nodes[si]->RBC_node_fi, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		RBC_node_transport->update_fi(moc[moc_idx]->nodes[ei]->RBC_node_fi, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);

		//moc[moc_idx]->nodes[0]->RBC_node_fi = 1.;
		if(moc[moc_idx]->nodes[si]->is_master_node){
			//RBC transport for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_RBC_transport){
				RBC_node_transport->update_master_fi(moc[moc_idx]->nodes[si]->RBC_node_fi, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->RBClum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);
			}

		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//RBC transport for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_RBC_transport){
				RBC_node_transport->update_master_fi(moc[moc_idx]->nodes[ei]->RBC_node_fi, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->RBClum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);
			}

		}

		//transport 1D update (actual edge)
		RBC1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->RBC_edge_fi, moc[moc_idx]->edges[e_idx]->RBC_edge_finew, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->RBC_node_fi, moc[moc_idx]->nodes[ei]->RBC_node_fi);
	}


	// HB saturation transport
	if(do_HBsat_transport){
		//HB saturation transport

		HB_O2_node_transport->update_fi(moc[moc_idx]->nodes[si]->HBsat_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		HB_O2_node_transport->update_fi(moc[moc_idx]->nodes[ei]->HBsat_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);

					
		if(moc[moc_idx]->nodes[si]->is_master_node){
			//HBsat transport for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_HB_sat_transport){
				HB_O2_node_transport->update_master_fi(moc[moc_idx]->nodes[si]->HBsat_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HBsatlum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);
			}

		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//HBsat transport for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_HB_sat_transport){
				HB_O2_node_transport->update_master_fi(moc[moc_idx]->nodes[ei]->HBsat_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HBsatlum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);
			}

		}

		//transport 1D update (actual edge)
		HBsat1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->HBsat_edge, moc[moc_idx]->edges[e_idx]->HBsat_edge_new, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->HBsat_node, moc[moc_idx]->nodes[ei]->HBsat_node);
	}



	// Plasma O2 concentration
	if(do_Plasma_O2_transport){
		//Plasma O2 concentration

		PlasmaO2_node_transport->update_fi(moc[moc_idx]->nodes[si]->PlasmaO2_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		PlasmaO2_node_transport->update_fi(moc[moc_idx]->nodes[ei]->PlasmaO2_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);

					
		if(moc[moc_idx]->nodes[si]->is_master_node){
			//Plasma O2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_PlasmaO2_transport){
				PlasmaO2_node_transport->update_master_fi(moc[moc_idx]->nodes[si]->PlasmaO2_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->PlasmaO2lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

				//capillary is only updated here
				lum[lum_idx]->capillary_O2_transport(moc[moc_idx]->edges[e_idx]->dt_act);
			}

		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//Plasma O2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_PlasmaO2_transport){
				PlasmaO2_node_transport->update_master_fi(moc[moc_idx]->nodes[ei]->PlasmaO2_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->PlasmaO2lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

				//capillary is only updated here
				lum[lum_idx]->capillary_O2_transport(moc[moc_idx]->edges[e_idx]->dt_act);
			}

		}

		//transport 1D update (actual edge)
		Plasma_O21D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->PlasmaO2_edge, moc[moc_idx]->edges[e_idx]->PlasmaO2_edge_new, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->PlasmaO2_node, moc[moc_idx]->nodes[ei]->PlasmaO2_node);
	}



}


//--------------------------------------------------------------
void first_blood::CO2_transport(int moc_idx, int si, int ei, int e_idx, double t_act){
	if(do_pla_CO2_transport){

		transport_node_CO2_pla->update_fi(moc[moc_idx]->nodes[si]->CO2_pla_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		transport_node_CO2_pla->update_fi(moc[moc_idx]->nodes[ei]->CO2_pla_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);


		if(moc[moc_idx]->nodes[si]->is_master_node){

			//Plasma CO2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_pla_CO2_transport){
				transport_node_CO2_pla->update_master_fi(moc[moc_idx]->nodes[si]->CO2_pla_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->CO2_pla_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

				//capillary is only updated here
				//lum[lum_idx]->CO2transport(moc[moc_idx]->edges[e_idx]->dt_act);
				lum[lum_idx]->capillary_CO2_transport(moc[moc_idx]->edges[e_idx]->dt_act);
			}
		}
					
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){

			//Plasma CO2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_pla_CO2_transport){
				//cout<<moc[moc_idx]->edges[e_idx]->dt_act<<endl;
				transport_node_CO2_pla->update_master_fi(moc[moc_idx]->nodes[ei]->CO2_pla_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->CO2_pla_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

				//capillary is only updated here
				//lum[lum_idx]->CO2transport(moc[moc_idx]->edges[e_idx]->dt_act);
				lum[lum_idx]->capillary_CO2_transport(moc[moc_idx]->edges[e_idx]->dt_act);
			}

		}

		//transport 1D update (actual edge)
		CO2_pla_1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->CO2_pla_edge_fi, moc[moc_idx]->edges[e_idx]->CO2_pla_edge_finew, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->CO2_pla_node, moc[moc_idx]->nodes[ei]->CO2_pla_node);
	}



	if(do_rbc_CO2_transport){

		transport_node_CO2_rbc->update_fi(moc[moc_idx]->nodes[si]->CO2_rbc_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		transport_node_CO2_rbc->update_fi(moc[moc_idx]->nodes[ei]->CO2_rbc_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);



		if(moc[moc_idx]->nodes[si]->is_master_node){
			//rbc cytoplasm CO2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_rbc_CO2_transport){
				transport_node_CO2_rbc->update_master_fi(moc[moc_idx]->nodes[si]->CO2_rbc_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->CO2_rbc_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//rbc cytoplasm CO2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_rbc_CO2_transport){
				transport_node_CO2_rbc->update_master_fi(moc[moc_idx]->nodes[ei]->CO2_rbc_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->CO2_rbc_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}

		//transport 1D update (actual edge)
		CO2_rbc_1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->CO2_rbc_edge_fi, moc[moc_idx]->edges[e_idx]->CO2_rbc_edge_finew, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->CO2_rbc_node, moc[moc_idx]->nodes[ei]->CO2_rbc_node);
	}
				

	if(do_pla_HCO3_transport){

		transport_node_HCO3_pla->update_fi(moc[moc_idx]->nodes[si]->HCO3_pla_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		transport_node_HCO3_pla->update_fi(moc[moc_idx]->nodes[ei]->HCO3_pla_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);


		if(moc[moc_idx]->nodes[si]->is_master_node){
			//rbc plasma HCO3 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_pla_HCO3_transport){
				transport_node_HCO3_pla->update_master_fi(moc[moc_idx]->nodes[si]->HCO3_pla_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HCO3_pla_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//rbc plasma HCO3 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_pla_HCO3_transport){
				transport_node_HCO3_pla->update_master_fi(moc[moc_idx]->nodes[ei]->HCO3_pla_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HCO3_pla_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
		HCO3_pla_1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->HCO3_pla_edge_fi, moc[moc_idx]->edges[e_idx]->HCO3_pla_edge_finew, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->HCO3_pla_node, moc[moc_idx]->nodes[ei]->HCO3_pla_node);
	}



	if(do_rbc_HCO3_transport){
		transport_node_HCO3_rbc->update_fi(moc[moc_idx]->nodes[si]->HCO3_rbc_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		transport_node_HCO3_rbc->update_fi(moc[moc_idx]->nodes[ei]->HCO3_rbc_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);


		if(moc[moc_idx]->nodes[si]->is_master_node){
			//rbc cytoplasm HCO3 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_rbc_HCO3_transport){
				transport_node_HCO3_rbc->update_master_fi(moc[moc_idx]->nodes[si]->HCO3_rbc_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HCO3_rbc_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//rbc cytoplasm HCO3 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_rbc_HCO3_transport){
				transport_node_HCO3_rbc->update_master_fi(moc[moc_idx]->nodes[ei]->HCO3_rbc_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HCO3_rbc_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
		HCO3_rbc_1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->HCO3_rbc_edge_fi, moc[moc_idx]->edges[e_idx]->HCO3_rbc_edge_finew, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->HCO3_rbc_node, moc[moc_idx]->nodes[ei]->HCO3_rbc_node);
	}


	if(do_HbCO2_transport){
		transport_node_HbCO2->update_fi(moc[moc_idx]->nodes[si]->HbCO2_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges);
		transport_node_HbCO2->update_fi(moc[moc_idx]->nodes[ei]->HbCO2_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges);


		if(moc[moc_idx]->nodes[si]->is_master_node){
			//HbCO2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;

			if (lum[lum_idx]->do_lum_HbCO2_transport){
				transport_node_HbCO2->update_master_fi(moc[moc_idx]->nodes[si]->HbCO2_node, moc[moc_idx]->nodes[si], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HbCO2_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
					
		if(moc[moc_idx]->nodes[ei]->is_master_node){
			//HbCO2 concentration for lum
			int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;

			if (lum[lum_idx]->do_lum_HbCO2_transport){
				transport_node_HbCO2->update_master_fi(moc[moc_idx]->nodes[ei]->HbCO2_node, moc[moc_idx]->nodes[ei], moc[moc_idx]->edges, *lum[lum_idx], *moc[moc_idx]);
				lum[lum_idx]->HbCO2_lum->update_fi(moc[moc_idx]->edges[e_idx]->dt_act, *lum[lum_idx], t_act);

			}
		}
		HbCO2_1D->update_fi(moc[moc_idx]->edges[e_idx]->get_velocity(), moc[moc_idx]->edges[e_idx]->HbCO2_edge_fi, moc[moc_idx]->edges[e_idx]->HbCO2_edge_finew, moc[moc_idx]->edges[e_idx]->length, moc[moc_idx]->edges[e_idx]->dt_act, moc[moc_idx]->nodes[si]->HbCO2_node, moc[moc_idx]->nodes[ei]->HbCO2_node);
	}
}

void first_blood::save_transport_var_for_lum(int moc_idx, int si, int ei){
	if(moc[moc_idx]->nodes[ei]->is_master_node){
		int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;
		solver_lumped* lum_e = lum[lum_idx];

		if (lum_e->do_lum_PlasmaO2_transport){
			if(lum_e->PlasmaO2lum->do_save_results){
				lum_e->PlasmaO2lum->save_variables();}
		}

		if (lum_e->do_lum_HB_sat_transport){
			if(lum_e->HBsatlum->do_save_results){
				lum_e->HBsatlum->save_variables();}
		}

		if (lum_e->do_lum_RBC_transport){
			if(lum_e->RBClum->do_save_results){
				lum_e->RBClum->save_variables();}
		}

		//CO2
		if(lum_e->do_lum_pla_CO2_transport){
			if(lum_e->CO2_pla_lum->do_save_results){
				lum_e->CO2_pla_lum->save_variables();}
		}

		if(lum_e->do_lum_rbc_CO2_transport){
			if(lum_e->CO2_rbc_lum->do_save_results){
				lum_e->CO2_rbc_lum->save_variables();}
		}

		if(lum_e->do_lum_pla_HCO3_transport){
			if(lum_e->HCO3_pla_lum->do_save_results){
				lum_e->HCO3_pla_lum->save_variables();}
		}

		if(lum_e->do_lum_rbc_HCO3_transport){
			if(lum_e->HCO3_rbc_lum->do_save_results){
				lum_e->HCO3_rbc_lum->save_variables();}
		}

		if(lum_e->do_lum_HbCO2_transport){
			if(lum_e->HbCO2_lum->do_save_results){
				lum_e->HbCO2_lum->save_variables();}
		}

	}

	if(moc[moc_idx]->nodes[si]->is_master_node){
		int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;
		solver_lumped* lum_s = lum[lum_idx];

		if (lum_s->do_lum_PlasmaO2_transport){
			if(lum_s->PlasmaO2lum->do_save_results){
				lum_s->PlasmaO2lum->save_variables();}
		}

		if (lum_s->do_lum_HB_sat_transport){
			if(lum_s->HBsatlum->do_save_results){
				lum_s->HBsatlum->save_variables();}
		}

		if (lum_s->do_lum_RBC_transport){
			if(lum_s->RBClum->do_save_results){
				lum_s->RBClum->save_variables();}
		}

		//CO2
		if(lum_s->do_lum_pla_CO2_transport){
			if(lum_s->CO2_pla_lum->do_save_results){
				lum_s->CO2_pla_lum->save_variables();}
		}

		if(lum_s->do_lum_rbc_CO2_transport){
			if(lum_s->CO2_rbc_lum->do_save_results){
				lum_s->CO2_rbc_lum->save_variables();}
		}

		if(lum_s->do_lum_pla_HCO3_transport){
			if(lum_s->HCO3_pla_lum->do_save_results){
				lum_s->HCO3_pla_lum->save_variables();}
		}

		if(lum_s->do_lum_rbc_HCO3_transport){
			if(lum_s->HCO3_rbc_lum->do_save_results){
				lum_s->HCO3_rbc_lum->save_variables();}
		}

		if(lum_s->do_lum_HbCO2_transport){
			if(lum_s->HbCO2_lum->do_save_results){
				lum_s->HbCO2_lum->save_variables();}
		}

	}
}

/*void first_blood::solve_lum(int index, double dt)
{
	vector<vector<double> > coefs;
	for(int j=0; j<lum[index]->boundary_indices.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0];
		int moc_edge_index = lum[index]->boundary_indices[j][1];
		int edge_end = lum[index]->boundary_indices[j][2];
		vector<double> c;
		if(edge_end==1)
		{
			c = moc[moc_index]->edges[moc_edge_index]->boundary_master_end(dt);
		}
		else
		{
			c = moc[moc_index]->edges[moc_edge_index]->boundary_master_start(dt);
		}
		coefs.push_back(c);
	}

	vector<vector<double> > qps = lum[index]->solve_one_step(dt, coefs);

	// updating master boundaries in mocs, n. 2
	for(int j=0; j<qps.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0]; // index of moc
		int moc_edge_index = lum[index]->boundary_indices[j][1]; // edge index in moc
		int edge_end = lum[index]->boundary_indices[j][2];
		int lum_node_index = lum[index]->boundary_indices[j][3]; // node index in lum

		double q = qps[j][0];
		double p = qps[j][1];

		double t_act = moc[moc_index]->edges[moc_edge_index]->time.back() + dt;
		if(edge_end == 1)
		{
			moc[moc_index]->edges[moc_edge_index]->boundary_end_variables(dt,p,q);
			int node_index = moc[moc_index]->edges[moc_edge_index]->node_index_end;
			moc[moc_index]->nodes[node_index]->boundary_variables(p,t_act);
		}
		else
		{
			moc[moc_index]->edges[moc_edge_index]->boundary_start_variables(dt,p,q);
			int node_index = moc[moc_index]->edges[moc_edge_index]->node_index_start;
			moc[moc_index]->nodes[node_index]->boundary_variables(p,t_act);
		}
	}
}*/
