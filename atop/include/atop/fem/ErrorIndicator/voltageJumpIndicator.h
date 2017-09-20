/*
 *
 *  Created on: Sep 20, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INCLUDE_ATOP_FEM_ERRORINDICATOR_VOLTAGEJUMPINDICATOR_H_
#define INCLUDE_ATOP_FEM_ERRORINDICATOR_VOLTAGEJUMPINDICATOR_H_

/**
 * This class calculates the jump in the voltages at the element edges and uses it to decide the elements for
 * refinement/coarsening.
 */

#include <vector>
#include<iostream>
#include <atop/fem/fem.h>
#include<atop/TopologyOptimization/cell_prop.h>



using namespace dealii;

namespace atop{
	template <int dim>
	class VoltageJumpIndicator{
	public:

		/**
		 * This vector stores the voltage values at the edges of an element
		 * thension euqals no. of faces, 4 in 2D and 6 in 3D
		 */
		std::vector<std::vector<double> > voltage_values;
		std::vector<double> *error_vector;
		std::vector<CellInfo> *cell_info_vector;

		FEM<dim> *fem;	//Pointer to the finite element object


		VoltageJumpIndicator(
				FEM<dim> &obj_fem,
				std::vector<double> &estimated_error_per_cell,
				std::vector<CellInfo> &cell_info_vector);
		void estimate();


	};

	template class VoltageJumpIndicator<2>;
}

#endif /* INCLUDE_ATOP_FEM_ERRORINDICATOR_VOLTAGEJUMPINDICATOR_H_ */
