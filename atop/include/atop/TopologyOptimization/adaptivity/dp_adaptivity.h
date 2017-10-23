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
		unsigned int get_design_bound(unsigned int p_order,
				FEM<dim> &fem);

		/**
		 * This function calculates the element level design bound for a given shape function
		 * It checks the neighbors and deducts the constrained dofs also.
		 */
		unsigned int get_corrected_design_bound(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&,
				typename hp::DoFHandler<dim>::active_cell_iterator&);


		/**
		 * This function returns the system level bound excluding the hanging nodes, without the mesh having been updated
		 */
		unsigned int get_corrected_system_design_bound(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&);

		void project_design(
				std::vector<double>&,
				std::vector<std::vector<double> >&,
				std::vector<double>&,
				std::vector<std::vector<double> >&);

		void update_p_order_contrast(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&);

		void update_design_for_elem_bound_only(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&);

/*
 * This function reduces the contrast in the number of design variables between adjacent elements.
 * It ensures that the contrast is not more than d-factor difference of 1 (will be adapted later).
 * For example, if there are 25 design points in a certain FE, then the neighbors should have a minimum of
 * 16 design points in a 2D space setting.
 * Within this function, there is another subroutine which checks if the element is completely void or solid,
 * and then checks it neighbors for the same. If this condition is satisfied, then the p- and d-orders of that element are not corrected
 * for the reduction of contrast.
 */
		void update_design_contrast(
				FEM<dim> &obj_fem,
				std::vector<CellInfo>&,
				unsigned int);
	};



	//template class dpAdaptivity<2>;
	template class dpAdaptivity<3>;
}



#endif /* DP_ADAPTIVITY_H_ */
