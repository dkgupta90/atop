/*
 *
 *  Created on: Aug 13, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/derivatives/compliance.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/base/quadrature.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_values.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/fem/fem.h>
#include <atop/math_tools/algebra/MatrixVector.h>
#include <atop/TopologyOptimization/DensityValues.h>
#include <iomanip>

using namespace dealii;
using namespace atop;

template <int dim>
void Compliance<dim>::set_input(
		hp::DoFHandler<dim> &obj_dof_handler,
		std::vector<CellInfo> &obj_cell_info_vector,
		std::vector<CellInfo> &obj_density_cell_info_vector,
		FEM<dim> &obj_fem
		){

	//Setting the dof_handler object
	this->dof_handler = &obj_dof_handler;
	this->cell_info_vector = &obj_cell_info_vector;
	this->density_cell_info_vector = &obj_density_cell_info_vector;
	this->fem = &obj_fem;
	this->elastic_data = &(obj_fem.elastic_data);
	this->density_field = &(obj_fem.density_field);

}

template <int dim>
void Compliance<dim>::compute(
		double &objective,
		std::vector<double> &obj_grad){

	std::cout<<"Calculating objective and sensitivities...."<<std::endl;

	objective = 0.0;

	/**
	 * Iterating over all the cells and including the contributions
	 */

	hp::FEValues<dim> hp_fe_values(fem->fe_collection,
			fem->quadrature_collection,
			update_values |
			update_gradients |
			update_quadrature_points |
			update_JxW_values
			);

	unsigned int cell_itr = 0;
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
			endc = dof_handler->end();

	for(; cell != endc; ++cell){
		unsigned int quadrature_rule = (*cell_info_vector)[cell_itr].quad_rule;
		QGauss<dim> quadrature_formula(quadrature_rule);

		unsigned int quad_index = elastic_data->get_quad_index(quadrature_rule);
		hp_fe_values.reinit(cell, quad_index);
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		const unsigned int n_q_points = quadrature_formula.size();
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
		cell_matrix = 0;
		std::vector<double> E_values, dE_values;
		E_values = (*cell_info_vector)[cell_itr].E_values;
		dE_values = (*cell_info_vector)[cell_itr].dE_values;

		unsigned int p_index = elastic_data->get_p_index((*cell_info_vector)[cell_itr].shape_function_order);

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			FullMatrix<double> normalized_matrix = elastic_data->elem_stiffness_array[p_index][quad_index][q_point];

			double area_factor = 1; //(*cell_info_vector)[cell_itr].cell_area/density_field->max_cell_area;//(cellprop[cell_itr].cell_area)/max_cell_area;
			cell_matrix.add((E_values[q_point])*area_factor,
					normalized_matrix);

		}


		//Extracting the nodal values of solution vector
		Vector<double> cell_array(dofs_per_cell);
		cell_array = 0;
		cell->get_dof_indices(local_dof_indices);
		for(unsigned int i = 0; i < dofs_per_cell; ++i){
			cell_array[i] = fem->solution(local_dof_indices[i]);
		}
		Vector<double> temp_array(dofs_per_cell);
		temp_array = 0;
		Matrix_Vector matvec;
		matvec.vector_matrix_multiply(
				cell_array,
				cell_matrix,
				temp_array,
				dofs_per_cell,
				dofs_per_cell);
		objective += matvec.vector_vector_inner_product(
				temp_array,
				cell_array);
			++cell_itr;
			//std::cout<<objective<<std::endl;
	}
	std::cout<<"Iteration: "<<fem->itr_count + 1<<"   Objective: "<<std::setprecision(10)<<objective<<std::setw(10)<<std::endl;


	//Calculating the sensitivities with respect to the density space design variables
	std::cout<<"Computing sensitivity response "<<std::endl;
	double time1 = clock();
	cell_itr = 0;

	//Initializing the obj_grad vector
	obj_grad.clear();
	obj_grad.resize(fem->density_field.get_design_count(
			fem->cycle,
			*(fem->mesh),
			*cell_info_vector,
			*density_cell_info_vector), 0.0);

	std::vector<std::vector<double> > temp_obj_grad((*cell_info_vector).size()); //used for uncoupled mesh, adaptivity_grayness
	//Below vector will be used for only adaptive_grayness and uncoupled mesh type
	for(unsigned int i = 0; i < cell_info_vector->size(); ++i){
		temp_obj_grad[i].resize((*cell_info_vector)[i].design_points.no_points, 0.0);
	}

	cell = dof_handler->begin_active(),
					endc = dof_handler->end();
	for(; cell != endc; ++cell){
		unsigned int quadrature_rule = (*cell_info_vector)[cell_itr].quad_rule;
		unsigned int quad_index = elastic_data->get_quad_index(quadrature_rule);
		unsigned int p_index = elastic_data->get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		QGauss<dim> quadrature_formula(quadrature_rule);
		std::vector<double> qweights = quadrature_formula.get_weights();	//Getting the quadrature weights
		hp_fe_values.reinit(cell, quad_index);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
		std::vector<Point<dim> > qpoints = fe_values.get_quadrature_points();
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		const unsigned int n_q_points = quadrature_formula.size();
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		std::vector<double> E_values, dE_values;
		dE_values.clear();
		dE_values = (*cell_info_vector)[cell_itr].dE_values;

		//Extracting the nodal values of solution vector
		Vector<double> cell_array(dofs_per_cell);
		cell_array = 0;
		cell->get_dof_indices(local_dof_indices);
		for(unsigned int i = 0; i < dofs_per_cell; ++i){
			cell_array[i] = fem->solution(local_dof_indices[i]);
		}

		/**
		 * cleaning the dcdx vector for cell_info_vector
		 * This vector stores the dc_dxPhys for each of the quadrature points
		 */
		(*cell_info_vector)[cell_itr].dc_dx.clear(); //Not sure if used

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){

			//Getting the normalized matrix corresponding to the quadrature point
			//std::cout<<p_index<<"   "<<quad_index<<std::endl;
			FullMatrix<double> normalized_matrix = elastic_data->elem_stiffness_array[p_index][quad_index][q_point];

			//Getting dE_dxPhys
			double dE_dxPhys = dE_values[q_point];

			double dobj;

			if (fem->mesh->coupling == false && fem->mesh->adaptivityType == "movingdesignpoints"){
				for(unsigned int i = 0 ; i < (*cell_info_vector)[cell_itr].neighbour_cells[q_point].size(); ++i){

					unsigned int density_cell_itr2 = (*cell_info_vector)[cell_itr].neighbour_cells[q_point][i];
					std::vector<double> dxPhys_dx(fem->mesh->design_var_per_point(), 0.0);	//vector for derivatives of xPhys w.r.t all design variables for that point
					density_field->get_dxPhys_dx(
							dxPhys_dx,
							(*cell_info_vector)[cell_itr],
							q_point,
							qpoints[q_point],
							(*density_cell_info_vector)[density_cell_itr2],
							density_cell_itr2);



					//Calculating all the sensitivities for the particular design point
					for (unsigned int k = 0; k < fem->mesh->design_var_per_point(); k++){
						if (fem->itr_count < 0){
							if (k > 0){
								dxPhys_dx[k] = 0.0;
							}
						}
						(*density_cell_info_vector)[density_cell_itr2].dxPhys[k] += dxPhys_dx[k];
						unsigned int design_index = (density_cell_itr2 * fem->mesh->design_var_per_point()) + k;
						double dEfactor = dE_dxPhys * dxPhys_dx[k];
						cell_matrix = 0;
						cell_matrix.add(dEfactor,
								normalized_matrix);

						Vector<double> temp_array(dofs_per_cell);
						temp_array = 0;
						Matrix_Vector matvec;
						matvec.vector_matrix_multiply(
								cell_array,
								cell_matrix,
								temp_array,
								dofs_per_cell,
								dofs_per_cell);
						dobj = matvec.vector_vector_inner_product(
								temp_array,
								cell_array);
						//Adding to the grad vector
						//std::cout<<dxPhys_dx[k]<<" ";
						obj_grad[design_index] -= dobj;
					}

				}
			}
			else if (fem->mesh->coupling == true){
				for(unsigned int i = 0 ; i < (*cell_info_vector)[cell_itr].neighbour_cells[q_point].size(); ++i){
					unsigned int density_cell_itr2 = (*cell_info_vector)[cell_itr].neighbour_cells[q_point][i];
					double dxPhys_dx = density_field->get_dxPhys_dx(
							(*cell_info_vector)[cell_itr],
							q_point,
							density_cell_itr2);

					//Adding the dxPhys_dx information into the density_cell_info_vector
					double area_factor = (*cell_info_vector)[cell_itr].cell_area/density_field->max_cell_area;
					(*density_cell_info_vector)[density_cell_itr2].dxPhys[0] += (dxPhys_dx * area_factor);

					double dEfactor = dE_dxPhys * dxPhys_dx;
					cell_matrix = 0.0;
					cell_matrix.add(dEfactor,
							normalized_matrix);

					Vector<double> temp_array(dofs_per_cell);
					temp_array = 0;
					Matrix_Vector matvec;
					matvec.vector_matrix_multiply(
							cell_array,
							cell_matrix,
							temp_array,
							dofs_per_cell,
							dofs_per_cell);
					dobj = matvec.vector_vector_inner_product(
							temp_array,
							cell_array);
					//Adding to the grad vector
					obj_grad[density_cell_itr2] -= dobj;
				}
			}
			else{
				for (unsigned int i = 0; i < (*cell_info_vector)[cell_itr].neighbour_points[q_point].size(); ++i){
					unsigned int cell_itr2 = (*cell_info_vector)[cell_itr].neighbour_points[q_point][i].first;	//index of neighbor cell
					unsigned int ngpt_itr = (*cell_info_vector)[cell_itr].neighbour_points[q_point][i].second;	//neighbor point index
					double dxPhys_dx = density_field->get_dxPhys_dx(
							(*cell_info_vector)[cell_itr],
							q_point,
							cell_itr2,
							ngpt_itr);

					//Adding the dxPhys_dx information into the cell containing this design point
					(*cell_info_vector)[cell_itr2].design_points.dxPhys_drho[ngpt_itr] += (qweights[q_point] * dxPhys_dx);

					//double area_factor = 1; //(*cell_info_vector)[cell_itr].cell_area/density_field->max_cell_area;
					//(*density_cell_info_vector)[density_cell_itr2].dxPhys[0] += (dxPhys_dx * area_factor);

					double dEfactor = dE_dxPhys * dxPhys_dx;
					cell_matrix = 0.0;
					cell_matrix.add(dEfactor,
							normalized_matrix);

					Vector<double> temp_array(dofs_per_cell);
					temp_array = 0;
					Matrix_Vector matvec;
					matvec.vector_matrix_multiply(
							cell_array,
							cell_matrix,
							temp_array,
							dofs_per_cell,
							dofs_per_cell);
					dobj = matvec.vector_vector_inner_product(
							temp_array,
							cell_array);
					//Adding to the grad vector
					temp_obj_grad[cell_itr2][ngpt_itr] -= dobj;
				}

			}

		}
		cell_itr++;
	}

	//Adding data to obj_grad for uncoupled meshes
	if (fem->mesh->coupling == false && fem->mesh->adaptivityType != "movingdesignpoints"){
		unsigned int k = 0;
		for(unsigned int i = 0; i < temp_obj_grad.size(); ++i){
			for (unsigned int j = 0; j < temp_obj_grad[i].size(); ++j){
				obj_grad[k] = temp_obj_grad[i][j];
				//std::cout<<k<<"    "<<obj_grad[k]<<std::endl;
				k++;
			}
		}
	}

	//std::cout<<obj_grad[0]<<"   "<<obj_grad[24]<<std::endl;
	std::cout<<"Size of sensitivity vector : "<<obj_grad.size()<<std::endl;
	double time2 = clock();
	time2 = (time2 - time1)/(double)CLOCKS_PER_SEC;
	std::cout<<"CPU time for sensitivity analysis: "<<time2<<" seconds"<<std::endl;
}

