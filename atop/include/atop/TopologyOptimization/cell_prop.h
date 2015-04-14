/*
 * cell_prop.h
 *
 *  Created on: Dec 24, 2014
 *      Author: dkgupta
 */

#ifndef CELL_PROP_H_
#define CELL_PROP_H_

#include<iostream>
#include<vector>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include<deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>

using namespace dealii;
namespace topopt{
class CellProperties{

	//The first 'vector' of all the variables denotes the series of quadrature points

public:
	unsigned int quadrature_formula;
	std::vector<double> material_density;
	std::vector<double> xPhys;
	std::vector<double> E_values;
	std::vector<double> dE_values;
	double cell_area;
	unsigned int n_q_points;
	std::vector<std::vector< std::pair<unsigned int , unsigned int> > > neighbour_cells;
	/*First index = quadrature point
	 * Second index = neighbour index
	 * Pair index1 = neighbour cell iterator
	 * Pair index2 = neighbour quadrature point
	*/
	std::vector< std::vector<double> > neighbour_distance;
	std::vector< std::vector<double> > neighbour_weights;
	std::vector<double> dxPhys;
	void initialize_density(double value);
};
}




#endif /* CELL_PROP_H_ */
