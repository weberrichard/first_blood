#include "first_blood.h"

using namespace Eigen;

first_blood::first_blood(string file_name, double a_time_end)
{
   input_file_path = file_name;

   // getting rid of global path
   case_name = input_file_path.substr(input_file_path.rfind('/')+1);
   // getting rid of extension
   case_name = case_name.substr(0,case_name.length()-4);

   // some feedback to console
   cout << " [*] Case name: " << case_name << endl;

   // getting the path of the folder
   input_folder_path = input_file_path.substr(0,input_file_path.rfind('/')+1);

   // getting the file format
   string file_format = input_file_path.substr(input_file_path.length()-3,3);
   if(file_format == "csv")
   {
      load_system_csv();
   }
   else
   {
      cout << endl << "Unkown file format: " << file_format << endl << "Available file formats are: csv" << endl;
      exit(-1);
   }

   // calculating the number of nodes and edges
   number_of_edges = edges.size();
   number_of_nodes = nodes.size();

   // number of timesteps is an input data
   time_end = a_time_end;

   // Finding the indicies of nodes for the edges and vica versa
   build_system();

   // setting initial time
   time.push_back(0.);
}

//--------------------------------------------------------------
first_blood::~first_blood(){}

//--------------------------------------------------------------
void first_blood::print_input()
{  
   printf("\n EDGES \n ----- \n");
   printf("\n %8s, %8s, %8s, %8s, %6s, %6s, %6s, %3s, %9s, %9s, %9s, %9s, %9s", "type", "name", "n_start", "n_end", "d", "s", "l", "nx", "E1", "E2", "eta", "Rs", "Re");
   for(unsigned int i=0; i<number_of_edges; i++)
   {
      edges[i]->print_input();
   }
   printf("\n NODES \n ----- \n");
   printf("\n %8s, %8s, %6s, %6s, %6s, %6s", "type", "name", "0.", "p0", "R", "Ri");
   for(unsigned int i=0; i<number_of_nodes; i++)
   {
      nodes[i]->print_input();
   }
}

//--------------------------------------------------------------
void first_blood::build_system()
{
   // Clearing the in/out going edges from nodes
   for(unsigned int i=0; i<number_of_nodes; i++)
   {
      nodes[i]->edge_in.clear();
      nodes[i]->edge_out.clear();
   }

   for(unsigned int i=0; i<number_of_edges; i++)
   {
      // starting node
      int node_start = node_id_to_index(edges[i]->node_name_start);

      // ending node
      int node_end = node_id_to_index(edges[i]->node_name_end);

      // saving to edges
      edges[i]->node_index_start = node_start;
      edges[i]->node_index_end = node_end;

      // saving to nodes
      nodes[node_start]->edge_out.push_back(i);
      nodes[node_end]->edge_in.push_back(i);
   }
}

//--------------------------------------------------------------
int first_blood::node_id_to_index(string node_id)
{
  int i=0, idx=-1;
  bool got_it=false;
  while(i<number_of_nodes && !got_it)
  {
    if(node_id.compare(nodes[i]->name) == 0)
    {
      got_it = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n !!!WARNING!!!\n first_blood::node_id_to_index function\nNode is not existing, node_id: " << node_id << "\n Continouing..." << endl;
  }
  return idx;
}

//--------------------------------------------------------------
int first_blood::edge_id_to_index(string edge_id)
{
  int i=0, idx=-1;
  bool got_it=false;
  while(i<number_of_edges && !got_it)
  {
    if(edge_id.compare(edges[i]->name) == 0)
    {
      got_it = true;
      idx = i;
    }
    i++;
  }
  if(idx == -1)
  {
    cout << "\n!!!WARNING!!!\n first_blood::edge_id_to_index function\n Node is not existing, edge_id: " << edge_id << "\n Continouing..." << endl;
  }
  return idx;
}
