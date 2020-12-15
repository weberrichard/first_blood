#include "moc_solver.h"

//--------------------------------------------------------------
moc_solver::moc_solver(string filename, double a_time_end) : first_blood(filename, a_time_end){}

moc_solver::~moc_solver(){}

//--------------------------------------------------------------
void moc_solver::initialization()
{
   for(unsigned int i=0; i<number_of_nodes; i++)
   {
      nodes[i]->initialization();
   }
   for(unsigned int i=0; i<number_of_edges; i++)
   {
      edges[i]->initialization();
   }

   // setting back the pressure_upstream interpolation index to 0
   index_upstream = 0;
}

//--------------------------------------------------------------
void moc_solver::full_solver()
{
   initialization();
}

//--------------------------------------------------------------
void moc_solver::forward_solver(string node_id)
{
   // getting the forward tree edges and nodes
   forward_tree(node_id);

   int i=0;
   while(time.back() < time_end)
   {
      // getting the real time step i.e. the minimum timestep from every edge
      double dt_real=1.e10;
      for(unsigned int j=0; j<forward_edges.size(); j++)
      {
         double dt=edges[forward_edges[j]]->new_timestep();
         if(dt<dt_real)
         {
            dt_real = dt;
         }
      }
      time.push_back(time.back() + dt_real);

      printf("i: %3i, time: %8.5f, dt: %8.5f \n", i, time[i], dt_real);

      // calculating new pressure and velocity field in inner points
      for(unsigned int j=0; j<forward_edges.size(); j++)
      {
         edges[forward_edges[j]]->solve();
      }

      // handling nodes i.e. boundaries conditions of edges
      boundaries(dt_real);

      // interpolation, updateing and saving field variables
      for(unsigned int j=0; j<forward_edges.size(); j++)
      {
         // interpolating to equidistant mesh
         edges[forward_edges[j]]->interpolate(dt_real);
         // updating every field variables
         edges[forward_edges[j]]->update_field_variables(dt_real);
         // saveing start and end points vars in time
         edges[forward_edges[j]]->save_field_variables();
      }

      i++;
   }

   // saving the final number of time steps
   number_of_timesteps = time.size();
}

//--------------------------------------------------------------
void moc_solver::boundaries(double dt)
{
   for(unsigned int i=0; i<forward_nodes.size(); i++)
   {
      if(nodes[forward_nodes[i]]->is_upstream_boundary) // handling the upstream boundary
      {
         // finding the position for linear interpolation
         int j=index_upstream;
         bool got_it = false;
         while(!got_it)
         {
            // making the inlet function periodic
            if(j > time_upstream.size()-1)
            {
               j -= time_upstream.size();
               period += 1;
            }

            double t = time.back()-period*time_upstream.back();

            if(t >= time_upstream[j] && t < time_upstream[j+1])
            {
               got_it=true;
               index_upstream = j;
            }
            else
            {
               j++;
            }
         }

         // interpolating
         double p_h = pressure_upstream[index_upstream+1]; // pressure at higher index
         double p_l = pressure_upstream[index_upstream]; // pressure at lower index
         double t_h = time_upstream[index_upstream+1]; // time at higher index
         double t_l = time_upstream[index_upstream]; // time at lower index

         double t_in = time.back()-period*time_upstream.back(); // actual time of the simulation
         double p_in = (p_h-p_l)/(t_h-t_l) * (t_in-t_l) + p_l; // actual pressure of the simulation

         double q_in=0.;
         // there might be several outgoing edge from an upstream node
         for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_out.size(); j++)
         {
            int edge_index = nodes[forward_nodes[i]]->edge_out[j];
            // calculating vp and pp values of edge
            q_in += edges[edge_index]->upstream_boundary(dt, p_in);
         }

         nodes[forward_nodes[i]]->pressure.push_back(p_in);
         nodes[forward_nodes[i]]->volume_flow_rate.push_back(q_in);
      }
      else if(nodes[forward_nodes[i]]->type_code == 1) // perifera
      {
         // there can be only one incoming edge
         int edge_index = nodes[forward_nodes[i]]->edge_in[0];
         double p_out = nodes[forward_nodes[i]]->pressure_out;
         double q_out = edges[edge_index]->boundary_periferia(dt,p_out);

         nodes[forward_nodes[i]]->pressure.push_back(p_out);
         nodes[forward_nodes[i]]->volume_flow_rate.push_back(q_out);
      }
      else if(nodes[forward_nodes[i]]->type_code == 0) // junctions
      {
         // first the node pressure is calculated

         // temp variables
         double num=0., denum=0.;

         // leakage part if there is, otherwise nb=[0,0]
         vector<double> nb = nodes[forward_nodes[i]]->boundary_coefficients();
         num += nb[0];
         denum += nb[1];

         // handling edges

         // for saving edge coefs
         vector<vector<double> > edge_coefs_in, edge_coefs_out;

         // ingoing edges
         for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_in.size(); j++)
         {
            int edge_index = nodes[forward_nodes[i]]->edge_in[j];
            vector<double> eb = edges[edge_index]->boundary_end_coefficients(dt);
            num   -= eb[1];
            denum += eb[0];
            edge_coefs_in.push_back(eb);
         }

         // outgoing edges
         for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_out.size(); j++)
         {
            int edge_index = nodes[forward_nodes[i]]->edge_out[j];
            vector<double> eb = edges[edge_index]->boundary_start_coefficients(dt);
            num   += eb[1];
            denum -= eb[0];
            edge_coefs_out.push_back(eb);
         }

         // nodal pressure
         double p_nodal = num/denum;
         nodes[forward_nodes[i]]->boundary_variables(p_nodal);

         double q_nodal = 0;
         // edge velocity and pressure
         for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_in.size(); j++)
         {
            double q = edge_coefs_in[j][0]*p_nodal + edge_coefs_in[j][1];
            q_nodal += q; 
            int edge_index = nodes[forward_nodes[i]]->edge_in[j];
            edges[edge_index]->boundary_end_variables(dt,p_nodal,q);
         }

         for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_out.size(); j++)
         {
            double q = edge_coefs_out[j][0]*p_nodal + edge_coefs_out[j][1];
            q_nodal -= q;
            int edge_index = nodes[forward_nodes[i]]->edge_out[j];
            edges[edge_index]->boundary_start_variables(dt,p_nodal,q);
         }

         nodes[forward_nodes[i]]->pressure.push_back(p_nodal);
         nodes[forward_nodes[i]]->volume_flow_rate.push_back(q_nodal);
      }
   }
}

//--------------------------------------------------------------
void moc_solver::backward_solver(string node_id)
{

}

//--------------------------------------------------------------
void moc_solver::forward_tree(string node_id)
{
   forward_edges.clear();
   forward_nodes.clear();

   // finding node index
   int node_index = node_id_to_index(node_id);
   if(node_index == -1)
   {
      cout << "\n ! ERROR ! Node ID not found: " << node_id << "\n Exiting..." << endl;
      exit(-1);
   }

   // finding edges directly connected to the node
   for(unsigned int i=0; i<number_of_edges; i++)
   {
      if(node_index == edges[i]->node_index_start)
      {
         forward_edges.push_back(i);
         forward_nodes.push_back(edges[i]->node_index_start);
         forward_nodes.push_back(edges[i]->node_index_end);
      }
   }

   vector<int> ns;
   vector<int> es;
   for(unsigned int i=0; i<number_of_edges; i++)
   {
      ns.push_back(edges[i]->node_index_start);
      es.push_back(i);
   }

   // finding the rest of the forward tree
   int i=0;
   while(i<forward_edges.size())
   {
      int j=0;
      bool no_new = true;
      while(j<ns.size() && no_new)
      {
         if(edges[forward_edges[i]]->node_index_end == ns[j])
         {
            forward_edges.push_back(es[j]);
            forward_nodes.push_back(edges[es[j]]->node_index_start);
            forward_nodes.push_back(edges[es[j]]->node_index_end);
            no_new = false;
            i = 0;

            ns.erase(ns.begin()+j);
            es.erase(es.begin()+j);
         }
         j++;
      }

      if(no_new)
      {
         i++;
      }
   }

   // FOR DEBUG
   cout << "forward_edges:" << endl;
   for(unsigned int i=0; i<forward_edges.size(); i++)
   {
      cout << forward_edges[i] << endl;
   }
   cout <<"OK"<<endl;

   cout << "forward_nodes all:" << endl;
   for(unsigned int i=0; i<forward_nodes.size(); i++)
   {
      cout << forward_nodes[i] << endl;
   }
   cout << "OKOK" << endl;

   forward_nodes = unique(forward_nodes);

   cout << "forward_nodes unique:" << endl;
   for(unsigned int i=0; i<forward_nodes.size(); i++)
   {
      cout << forward_nodes[i] << endl;
   }
   cout << "OKOK" << endl;

   cin.get();
}

//-------------------------------------------------------------- 
vector<int> moc_solver::unique(vector<int> x)
{
  vector<int> out;
  for(int i=0; i<x.size(); i++)
  {
    bool unique = true;
    for(int j=0; j<out.size(); j++)
    {
      if(x[i] == out[j])
      {
        unique = false;
        break;
      }
    }
    if(unique)
    {
      out.push_back(x[i]);
    }
  }

  return out;
}


//--------------------------------------------------------------
int moc_solver::backward_tree(string node_id)
{
   int edge_index = -1;
   int i=0;
   while(edge_index == -1 && i<number_of_edges)
   {
      if(edges[i]->node_name_start == node_id)
      {
         edge_index = i;
      }
      i++;
   }

   return edge_index;
}