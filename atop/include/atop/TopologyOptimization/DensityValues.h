/*
 * DensityValues.h
 *
 *  Created on: Jul 14, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef DENSITYVALUES_H_
#define DENSITYVALUES_H_

#include <deal.II/base/point.h>
#include <vector>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <atop/TopologyOptimization/design_analysis.h>
#include <atop/TopologyOptimization/neighbors.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/physics/elasticity.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/fem/define_mesh.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

using namespace dealii;

namespace atop{
	template <int dim>
	class DensityField{
	public:
		double max_cell_area;
		unsigned int initial_no_cells;
		double cell_length;
		double volfrac;

		SparseMatrix<double> projection_matrix;	//For projecting the densities
		void create_neighbors(
				std::vector<CellInfo> &cell_info_vector,
				hp::FEValues<dim> &hp_fe_values,
				hp::DoFHandler<dim> &dof_handler,
				hp::DoFHandler<dim> &density_dof_handler,
				Projection &projection,
				DefineMesh<dim> &mesh
				);

		void create_neighbors(
				std::vector<CellInfo> &design_cell_info_vector,
				hp::DoFHandler<dim> &design_handler,
				Projection &projection
				);

		/**
		 * Function for applying projection/filter operation on density values
		 */
		void smoothing(
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector,
				DefineMesh<dim> &mesh
				);

		void smoothing(std::vector<CellInfo> &,
				std::vector<CellInfo> &
				);

		void update_design_vector(
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector,
				std::vector<double> &design_vector,
				unsigned int cycle,
				double volfrac,
				DefineMesh<dim> &mesh,
				Projection &projection);

		/* This function is written to update the lower and upper bounds for the
		 * design points used.
		 * For case where only density is used as design variable, this can be straightforward.
		 * For other cases, the vectors cannot be directly assigned and need to be properly iterated.
		 */

		void update_design_bounds(
				std::vector<double> &lb,
				std::vector<double> &ub,
				DefineMesh<dim> &mesh,
				Projection &projection,
				std::vector<double> &design_vector);

		void update_density_cell_info_vector(
				std::vector<CellInfo> &density_cell_info_vector,
				const std::vector<double> &design_vector);

		//overloaded function for uncoupled meshes, design variables are stored in cell_info for this case
		void update_density_cell_info_vector(
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector,
				const std::vector<double> &design_vector);

		double get_dxPhys_dx(CellInfo &cell_info,
				unsigned int q_point,
				unsigned int density_cell_itr2);

		double get_dxPhys_dx(CellInfo &cell_info,
				unsigned int q_point,
				unsigned int cell_itr2,
				unsigned int ngpt_itr);
		/*
		 * This implementation of get_dxPhys_dx is for movingDesignPoints.
		 * For every design point, it returns dim+2 values, each w.r.t 1 design variable
		 * design variables: rho, rmin, x_cor, y_cor
		 */
		void get_dxPhys_dx(
				std::vector<double> &dxPhys_dx,
				CellInfo &cell_info,
				unsigned int qpoint,
				Point<dim> qX,
				CellInfo &density_cell_info,
				unsigned int density_cell_itr2);

		double get_vol_fraction(
				std::vector<CellInfo> &cell_info_vector
				);

		unsigned int get_design_count(
				unsigned int cycle,
				DefineMesh<dim> &mesh,
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector);

		void get_xPhys_for_face(std::vector<double> &face_xPhys,
				hp::FECollection<dim> &temp_fe_coll,
				hp::QCollection<dim-1> &temp_q_coll,
				hp::DoFHandler<dim> &dofhandler,
				std::vector<CellInfo> &cell_info_vector,
				typename hp::DoFHandler<dim>::active_cell_iterator &cell,
				unsigned int face_itr);

	private:

		/**
		 * This function needs to be overloaded since there are at the moment
		 * some unknown issues in creating the template for this function.
		 * Currently, the implementation is only for dim = 2.
		 */
		void neighbor_search(
				hp::DoFHandler<2>::active_cell_iterator cell1,
				hp::DoFHandler<2>::active_cell_iterator cell,
				std::vector<hp::DoFHandler<2>::active_cell_iterator> &neighbor_iterators,
				double rmin
				);

		void calculate_weights(
				std::vector<CellInfo>  &cell_info_vector,
				unsigned int cell_itr1,
				double rmin,
				DefineMesh<dim> &mesh);
		void calculate_weights(
						std::vector<CellInfo>  &design_cell_info_vector,
						unsigned int,
						double);



	};

	template class DensityField<2>;
}


#endif /* DENSITYVALUES_H_ */
