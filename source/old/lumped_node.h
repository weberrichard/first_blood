/*===================================================================*\
										 lumped_node
									 ---------------

	 first_blood
	 R. Weber
\*==================================================================*/

#ifndef LUMPED_NODE_H
#define LUMPED_NODE_H

#include <string>
#include <vector>

using namespace std;

class lumped_node
{
public:
	lumped_node(string a_name);
	~lumped_node();

	// name of the lumped_node
	string name;

	// boundary type fo nodes
	string type;

	// contains indecies of edges going in and out
	vector<int> edge_in, edge_out;

	// containing pressure time series
	vector<double> pressure;
};

#endif // LUMPED_NODE_H