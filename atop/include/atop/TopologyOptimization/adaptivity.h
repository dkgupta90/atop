/*
 *
 *  Created on: Aug 22, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef ADAPTIVITY_H_
#define ADAPTIVITY_H_

#include <atop/fem/fem.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/adaptivity/dp_adaptivity.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <string>
#include <vector>

using namespace dealii;

namespace atop{
	template <int dim>
	class Adaptivity{
	public:
		FEM<dim> *fem;
		std::vector<CellInfo> *cell_info_vector;
		std::vector<CellInfo> *density_cell_info_vector;

		unsigned int system_design_bound;
		dpAdaptivity<dim> dp_adap;

		void update(
				FEM<dim> &obj_fem
				);

		void mesh_refine_indicator(
				std::string& mesh_update_str);
		void coupled_refine_adaptive_grayness();
		void calc_refinement_res_multires();
		void update_element_design_bound();

		void update_cell_vectors(
				std::vector<CellInfo> &density_cell_info_vector,
				hp::DoFHandler<dim> &density_dof_handler,
				Triangulation<dim> &density_triangulation);

		void execute_coarsen_refine();

		void increase_decrease_p_order();

	private:

		unsigned int rigid_body_modes;
		std::vector<double> refineRes;
		std::vector<std::pair<double, unsigned int> > sortedRefineRes;

		void compute_sortedRefineRes();
		void dp_coarsening_refinement();

		/*
		 * This function is an improved version of the dp-coarsening/refinement approach. Here the number of design variables
		 * can be non-perfect squares which implies that for certain p-order, a higher number number of design variables can be used.
		 * Note that the design distribution is generated using k-means clustering method and even though there might be geometrical
		 * symmetry, the optimized might not be exactly symmetrical since the design distributions aint.
		 */
		void improved_dp_coarsening_refinement();

		void run_dp_analysis_based_refinement();	//Only analysis part of the refinement is done here
	};

	template class Adaptivity<2>;
}



#endif /* ADAPTIVITY_H_ */
