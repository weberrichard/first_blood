/*===================================================================*\
										 moc_solver
									 ---------------

  This class is derived from the first_blood class and provides
  simple forward and backward calculations for edges and nodes.
 
	 first_blood
	 R. Weber
	 git: 
\*==================================================================*/

#ifndef MOC_SOLVER_H
#define MOC_SOLVER_H

#include "first_blood.h"

class moc_solver : public first_blood
{
public:

	moc_solver(string filename);
	~moc_solver();

	// giving initial conditions
	void initialization();

	// evaluating full calculation, several forward and backward calculation
	void full_solver(string node_id, double time_end);

	// perfroming the moc_solver forward
	// from a certain node_id to the perif
	void forward_solver(string node_id, double time_end);
	// handling the boundaries: upstream and inner BC nodes
	void boundaries(double dt);

	// performing the backward calculation on edge_id elemment
	void backward_solver(string node_id);

private:
	// determining the downward tree from a node for the forward_solver, returning the indicies of edges
	void forward_tree(string node_id);
	// for containing the forward tree edges
	vector<int> forward_edges;
	// for containing the forward tree nodes
	vector<int> forward_nodes;
	vector<int> unique(vector<int> x);

	// founding the upward edge from a node, returning the edge index
	int backward_tree(string node_id);
};

#endif // MOC_SOLVER_H