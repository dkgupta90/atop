/*
 * neighbors.h
 *
 *  Created on: Oct 16, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef NEIGHBORS_H_
#define NEIGHBORS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

using namespace dealii;
namespace topopt{
	class StoreIndices{
	public:
		std::vector< std::vector<std::pair<unsigned int, unsigned int> > > stored_indices;
		std::vector< std::vector<std::pair<unsigned int, double> > > density_indices;
		StoreIndices(unsigned int i);
		StoreIndices();
	};

	class FEneighbors{
	public:
		std::vector<double> E_rhoi;
		std::vector<double> dE_rhoi;
		std::vector<unsigned int> q_point;
		std::vector<DoFHandler<2>::active_cell_iterator > cell_list;

	};
}




#endif /* NEIGHBORS_H_ */
