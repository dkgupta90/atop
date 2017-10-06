/*
 *
 *  Created on: Aug 31, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/ErrorIndicator/stressJumpIndicator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>
#include <atop/physics/mechanics/elastic.h>


using namespace atop;
using namespace dealii;

template <int dim>
StressJumpIndicator<dim>::StressJumpIndicator(FEM<dim> &obj_fem,
		std::vector<double> &estimated_error_per_cell,
		std::vector<CellInfo> &obj_cell_info_vector){

	this->fem = &obj_fem;
	this->error_vector = &estimated_error_per_cell;
	this->cell_info_vector = &obj_cell_info_vector;

}

template <int dim>
void StressJumpIndicator<dim>::estimate(){

	//Updating face_B matrices
	fem->elastic_data.update_face_B_matrices(
			fem->fe_collection,
			fem->face_quadrature_collection,
			fem->dof_handler);

	//Getting the D-matrix
	ElasticTools<dim> elastic_tool;
	FullMatrix<double> D_matrix (3, 3);
	elastic_tool.get_D_plane_stress2D(D_matrix, 0.3);


	//std::cout<<"face_B_matrices updated "<<std::endl;

	//First the necessary parameters are declared/initialized here
	error_vector->clear();
	error_vector->resize(fem->triangulation.n_active_cells(), 0.0);

	hp::FEValues<dim> hp_fe_values(fem->fe_collection,
				fem->quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	hp::FEValues<dim> ng_hp_fe_values(fem->fe_collection,
				fem->quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	//Iterate over all the cells
	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	unsigned int cell_itr;
	for (; cell != endc; ++cell){
		cell_itr = cell->user_index() - 1;
		unsigned int p_index = fem->elastic_data.get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		unsigned int q_index = fem->elastic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);
		hp_fe_values.reinit(cell, q_index);
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

		//Get the state-field for the current cell
		Vector<double> cell_solution(dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i){
			cell_solution(i) = fem->solution(local_dof_indices[i]);
		}


		std::vector<double> face_xPhys;
		//Iterate over all the faces of the current cell
		for (unsigned int iface = 0; iface < GeometryInfo<dim>::faces_per_cell; ++iface){

			double face_stress_jump = 0.0;

			Vector<double> face_stress(3);
			Vector<double> ng_face_stress(3);
			face_stress = 0.0;
			ng_face_stress = 0.0;

			//Get the density values for all the quadrature points of the current face
			face_xPhys.clear();
			fem->density_field.get_xPhys_for_face(face_xPhys,
					fem->fe_collection,
					fem->face_quadrature_collection,
					fem->dof_handler,
					*cell_info_vector,
					cell,
					iface);


			unsigned int n_fq_points = fem->elastic_data.face_B_matrix_list[p_index][q_index][iface].size();

			if (face_xPhys.size() != n_fq_points)	std::cerr<<"Face quadrature points mismatch"<<std::endl;

			for (unsigned int fq = 0; fq < n_fq_points; ++fq){

				Vector<double> temp_cell_solution(dofs_per_cell);
				temp_cell_solution = cell_solution;
				temp_cell_solution *= (fem->elastic_data.face_JxW[p_index][q_index][iface][fq]);

				FullMatrix<double> calc_stress(3, dofs_per_cell);
				calc_stress = 0;

				D_matrix.mmult(calc_stress, fem->elastic_data.face_B_matrix_list[p_index][q_index][iface][fq]);

				calc_stress *= (*(fem->penal)).penalized_factor(face_xPhys[fq]);
				calc_stress.vmult_add(face_stress,
						temp_cell_solution);

			}

			//Check the stress value for this neighbor to this face of the cell
			//Do nothing if a Dirichlet boundary or no b.c. specified
			//Although there will be stress for the face, we store zero stresses for that phase
			//see dealii page on KellyErrorEstimator for the reasoning
			if (cell->at_boundary(iface) && cell->face(iface)->boundary_indicator() != 22){
				continue; //The condition above checks that faces with Dirichlet or no b.c. are skipped
			}

			//Something to be done for Neumann boundary

			typename hp::DoFHandler<dim>::active_cell_iterator ng_cell = cell->neighbor(iface);
			unsigned int ng_cell_itr = ng_cell->user_index() - 1;
			unsigned int ng_p_index = fem->elastic_data.get_p_index((*cell_info_vector)[ng_cell_itr].shape_function_order);
			unsigned int ng_q_index = fem->elastic_data.get_quad_index((*cell_info_vector)[ng_cell_itr].quad_rule);
			ng_hp_fe_values.reinit(ng_cell, ng_q_index);

			//Get the state-field for the current ng_cell
			const unsigned int ng_dofs_per_cell = ng_cell->get_fe().dofs_per_cell;
			Vector<double> ng_cell_solution(ng_dofs_per_cell);
			std::vector<types::global_dof_index> ng_local_dof_indices(ng_dofs_per_cell);
			ng_cell->get_dof_indices(ng_local_dof_indices);

			for (unsigned int i = 0; i < ng_dofs_per_cell; ++i){
				ng_cell_solution(i) = fem->solution(ng_local_dof_indices[i]);
			}

			ng_face_stress = 0.0;
			std::vector<double> ng_face_xPhys;

			//Iterating over all the faces to match the face
			for (unsigned int ng_iface = 0; ng_iface < GeometryInfo<dim>::faces_per_cell; ++ng_iface){
				if (ng_cell->neighbor(ng_iface) != cell)	continue;	//if faces do not match


				//Get the density values for all the quadrature points of the neighbour cell's current face
				ng_face_xPhys.clear();
				fem->density_field.get_xPhys_for_face(ng_face_xPhys,
						fem->fe_collection,
						fem->face_quadrature_collection,
						fem->dof_handler,
						*cell_info_vector,
						ng_cell,
						ng_iface);

				//if the common face found in the neighbor cell
				unsigned int ng_n_fq_points = fem->elastic_data.face_B_matrix_list[ng_p_index][ng_q_index][ng_iface].size();
				for (unsigned int ng_fq = 0; ng_fq < ng_n_fq_points; ++ng_fq){
					Vector<double> temp_ng_cell_solution(ng_dofs_per_cell);
					temp_ng_cell_solution = ng_cell_solution;
					temp_ng_cell_solution *= (fem->elastic_data.face_JxW[ng_p_index][ng_q_index][ng_iface][ng_fq]);

					FullMatrix<double> ng_calc_stress(3, ng_dofs_per_cell);
					ng_calc_stress = 0;
					D_matrix.mmult(ng_calc_stress, fem->elastic_data.face_B_matrix_list[ng_p_index][ng_q_index][ng_iface][ng_fq]);
					ng_calc_stress *= (*(fem->penal)).penalized_factor(ng_face_xPhys[ng_fq]);
					ng_calc_stress.vmult_add(ng_face_stress,
							temp_ng_cell_solution);
				}
			}

			ng_face_stress -= face_stress;
			//std::cout<<ng_face_stress(0)<<"  "<<ng_face_stress(1)<<"  "<<ng_face_stress(2)<<std::endl;
			face_stress_jump = ng_face_stress.l2_norm();
			(*error_vector)[cell_itr] += face_stress_jump;

			//std::cout<<cell_itr<<"  "<<iface<<"  "<<face_stress(0)<<"  "<<face_stress(1)<<"  "<<face_stress(2)<<std::endl;
		}

		(*error_vector)[cell_itr] /= (*cell_info_vector)[cell_itr].shape_function_order;
		//std::cout<<cell_itr<<"   "<<(*error_vector)[cell_itr]<<std::endl;

	}
}


