/*
 *
 *  Created on: Aug 30, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INCLUDE_ATOP_FEM_ERRORINDICATOR_STRESSJUMPINDICATOR_H_
#define INCLUDE_ATOP_FEM_ERRORINDICATOR_STRESSJUMPINDICATOR_H_


/**
 * This class calculates the jump in the stresses at the element edges and uses it to decide the elements for
 * refinement/coarsening.
 */

#include <vector>
#include<iostream>
#include <atop/fem/fem.h>
#include<atop/TopologyOptimization/cell_prop.h>



using namespace dealii;

namespace atop{
	template <int dim>
	class StressJumpIndicator{
	public:

		/**
		 * This vector stores the stress values at the edges of an element
		 * first dimension denotes the number of stress compoenents e.g. 3 for 2D, 6 for 3D
		 * second dimension denotes the number of faces e.g. 4 for 2D, 6 for 3D
		 */
		std::vector<std::vector<double> > stress_values;
		std::vector<double> *error_vector;
		std::vector<CellInfo> *cell_info_vector;


		FEM<dim> *fem;	//Pointer to the finite element object


		StressJumpIndicator(
				FEM<dim> &obj_fem,
				std::vector<double> &estimated_error_per_cell,
				std::vector<CellInfo> &cell_info_vector);
		void estimate();


	};

	//template class StressJumpIndicator<2>;
	template class StressJumpIndicator<3>;
}




#endif /* INCLUDE_ATOP_FEM_ERRORINDICATOR_STRESSJUMPINDICATOR_H_ */
