/*===================================================================*\
										 lumped_edge
									 ---------------                            

	 first_blood
	 R. Weber
\*==================================================================*/

#ifndef MOC_EDGE_H
#define MOC_EDGE_H

#include <string>
#include <vector>

using namespace std;

class lumped_edge
{
public:
	lumped_edge(string a_name);
	~lumped_edge();

	// name or ID of the lumped_edge
	string name;

	// name of the nodes at the beginning and at the end
	string node_name_start, node_name_end;

	// index of the nodes at the beginning and at the end
	int node_index_start, node_index_end;

	// type of edge
	string type;
	// 0: resistance, 1: capacity, 2: coil, 3: const pres, 4: diode
	int type_code;

	// coefficient of the edge
	double parameter;
};

#endif // MOC_EDGE_H
