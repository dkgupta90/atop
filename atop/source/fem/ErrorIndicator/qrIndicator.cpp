/*
 *
 *  Created on: Aug 7, 2017

 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/ErrorIndicator/qrIndicator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>
#include <atop/physics/mechanics/elastic.h>

using namespace dealii;
using namespace atop;

template <int dim>
QRIndicator<dim>::QRIndicator(
		FEM<dim> &obj_fem,
		std::vector<double> &qr_accuracy,
		double tol_accuracy,
		std::vector<unsigned int> &qr_p_value,
		std::vector<CellInfo> &cell_info_vector
		){
	this->cell_info_vector = &cell_info_vector;
	this->accuracy_vector = &qr_accuracy;
	this->proposed_p_values = &qr_p_value;
	this->fem = &obj_fem;
	this->tol_accuracy = tol_accuracy;
}


/* This function checks the J/Jstar values and proposed the p-order to be used for each finite element*/
template <int dim>
void QRIndicator<dim>::estimate(){

	// Getting the 1-time D-matrix
	ElasticTools elastic_tool;
	FullMatrix<double> D_matrix (3, 3);
	elastic_tool.get_D_plane_stress2D(D_matrix, 0.3);

	/*
	 * Below, an iteration over all the cells is performed.
	 * For each cell, first the solution and load and computed at all the dofs.
	 * Next, solution check is performed with increased or decreased value of p
	 * A decision is made on the value of p for each cell.
	 */

	hp::FEValues<dim> hp_fe_values(fem->fe_collection,
				fem->quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	unsigned int cell_itr;
	for (; cell != endc; ++cell){
		cell_itr = cell->user_index() - 1;
		unsigned int p_index = fem->elastic_data.get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		unsigned int q_index = fem->elastic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);
		hp_fe_values.reinit(cell, q_index);
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

		//Get the state solution for the current cell
		Vector<double> u_solution(dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i){
			u_solution(i) = fem->solution(local_dof_indices[i]);
		}

		// Get the stiffness matrix for this cell for the current p-value
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
		normalized_matrix = 0;
		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			normalized_matrix = 0;
			normalized_matrix = elastic_data.elem_stiffness_array[p_index][q_index][q_point];

			//NaN condition check ----------------------------------------------------------------------------------
			if ((*cell_info_vector)[cell_itr].E_values[q_point] != (*cell_info_vector)[cell_itr].E_values[q_point])
				std::cout<<q_point<<(*cell_info_vector)[cell_itr].E_values[q_point]<<std::endl;
			//------------------------------------------------------------------------------------------------------
			cell_matrix.add((*cell_info_vector)[cell_itr].E_values[q_point],
					normalized_matrix);
		}


	}
}


