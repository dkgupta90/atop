/*
 *
 *  Created on: Aug 14, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef DP_ADAPTIVITY_H_
#define DP_ADAPTIVITY_H_

#include <iostream>
#include <vector>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/designField.h>
#include <atop/fem/fem.h>

using namespace std;

namespace atop{
	template <int dim>
	class dpAdaptivity{
	public:

		unsigned int rigid_body_modes;
		void update_designField(std::vector<CellInfo>&,
				unsigned int,
				unsigned int);

		void correctify_p_order(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&,
				unsigned int);

		unsigned int get_system_design_bound(
				FEM<dim> &fem);

		/**
		 * This function simply returns the element level bound based on the dofs of the element
		 * and rigid body modes of that element.
		 */
		unsigned int get_design_bound(unsigned int p_order);

		/**
		 * This function calculates the element level design bound for a given shape function
		 * It checks the neighbors and deducts the constrained dofs also.
		 */
		void get_corrected_design_bound(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&,
				typename hp::DoFHandler<dim>::active_cell_iterator&);
	};

	template class dpAdaptivity<2>;
}



#endif /* DP_ADAPTIVITY_H_ */
