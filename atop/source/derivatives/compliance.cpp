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
#include <fstream>

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

	Matrix_Vector matvec;
	objective = matvec.vector_vector_inner_product(
					fem->system_rhs,
					fem->solution);

	std::cout<<"Iteration: "<<fem->itr_count + 1<<"   Objective: "<<std::setprecision(10)<<objective<<std::setw(10)<<std::endl;

	//exit(0);

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


		hp_fe_values.reinit(cell, 0);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
		//std::vector<Point<dim> > qpoints = fe_values.get_quadrature_points();
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		unsigned int design_count = (*cell_info_vector)[cell_itr].connected_cell_iterators_2D.size();

		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		std::vector<double> E_values, dE_values;
		dE_values.clear();
		dE_values = (*cell_info_vector)[cell_itr].dE_values;
		assert(dE_values.size() == design_count);

		//Extracting the nodal values of solution vector
		Vector<double> cell_array(dofs_per_cell);
		cell_array = 0;
		cell->get_dof_indices(local_dof_indices);
		for(unsigned int i = 0; i < dofs_per_cell; ++i){
			cell_array[i] = fem->solution(local_dof_indices[i]);
		}

		/**
		 * cleaning the dcdx vector for cell_info_vector
		 * This vector stores the dc_dxPhys for each of the design cell centers
		 */
		(*cell_info_vector)[cell_itr].dc_dx.clear(); //Not sure if used


		for(unsigned int design_point = 0; design_point < design_count; ++design_point){

			//Getting the normalized matrix corresponding to the design_point
			//std::cout<<p_index<<"   "<<quad_index<<std::endl;
			FullMatrix<double> normalized_matrix = elastic_data->K_matrix_list[0][0][design_point];

			//Getting dE_dxPhys
			double dE_dxPhys = dE_values[design_point];

			double dobj;

			if (fem->mesh->coupling == true){

			}
			else{

				unsigned int design_itr = ((*cell_info_vector)[cell_itr].connected_cell_iterators_2D[design_point])->user_index() - 1;

				unsigned int no_neighbors = (*density_cell_info_vector)[design_itr].neighbor_cells.size();

				/*Iterating over all the neighbors*/
				for (unsigned int ii = 0; ii < no_neighbors; ++ii){
					//std::cout<<cell_itr<<"   "<<design_point<<"  "<<design_itr<<"  "<<ii<<"  "<<no_neighbors<<std::endl;

					//Getting dxPhys_dx
					double dxPhys_dx = (*density_cell_info_vector)[design_itr].neighbor_weights[ii];
					///std::cout<<"Collapsed here "<<std::endl;

					unsigned int ng_itr = (*density_cell_info_vector)[design_itr].neighbor_cells[ii];
					//std::cout<<ng_itr<<"  "<<dxPhys_dx<<std::endl;

					(*density_cell_info_vector)[ng_itr].dxPhys[0] += dxPhys_dx;

					//(*cell_info_vector)[cell_itr].pseudo_design_points.dxPhys_drho[ng_itr] = dxPhys_dx;

					//We give dx_drho = 1
					double dEfactor = dE_dxPhys * dxPhys_dx * 1;
					cell_matrix = 0.0;
					cell_matrix.add(dEfactor,
							normalized_matrix);
					//std::cout<<"Collapsed here 2"<<std::endl;

					Vector<double> temp_array(dofs_per_cell);
					temp_array = 0;
					Matrix_Vector matvec;
					matvec.vector_matrix_multiply(
							cell_array,
							cell_matrix,
							temp_array,
							dofs_per_cell,
							dofs_per_cell);
					//std::cout<<"Collapsed here 3"<<std::endl;

					dobj = matvec.vector_vector_inner_product(
							temp_array,
							cell_array);

					obj_grad[ng_itr] -= dobj;
					//std::cout<<"Collapsed here 4"<<std::endl;


/*
 *
					//Iterating over all the design points of the cell
					for (unsigned int j = 0; j < (*cell_info_vector)[cell_itr].pseudo_design_points.dx_drho[ngpt_itr].size(); ++j){

						double dx_drho = (*cell_info_vector)[cell_itr2].pseudo_design_points.dx_drho[ngpt_itr][j];
						if ((fabs(dx_drho) - 0) < 1e-12)	continue;


						(*cell_info_vector)[cell_itr2].design_points.dxPhys_drho[j] += (qweights[q_point] * dxPhys_dx * dx_drho);

						double dEfactor = dE_dxPhys * dxPhys_dx * dx_drho;
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
						temp_obj_grad[cell_itr2][j] -= dobj;

					}*/

				}

/*				for (unsigned int i = 0; i < (*cell_info_vector)[cell_itr].neighbour_points[q_point].size(); ++i){
					unsigned int cell_itr2 = (*cell_info_vector)[cell_itr].neighbour_points[q_point][i].first;	//index of neighbor cell
					unsigned int ngpt_itr = (*cell_info_vector)[cell_itr].neighbour_points[q_point][i].second;	//neighbor point index
					double dxPhys_dx = density_field->get_dxPhys_dx(
							(*cell_info_vector)[cell_itr],
							q_point,
							cell_itr2,
							ngpt_itr);
					double dxPhys_dx = 0;

					//Adding the dxPhys_dx information into the cell containing this pseudo-design point
					//(*cell_info_vector)[cell_itr2].design_points.dxPhys_drho[ngpt_itr] += (qweights[q_point] * dxPhys_dx);
					(*cell_info_vector)[cell_itr2].pseudo_design_points.dxPhys_drho[ngpt_itr] += (qweights[q_point] * dxPhys_dx);


					//Iterating over all the design points of the cell
					for (unsigned int j = 0; j < (*cell_info_vector)[cell_itr2].pseudo_design_points.dx_drho[ngpt_itr].size(); ++j){

						double dx_drho = (*cell_info_vector)[cell_itr2].pseudo_design_points.dx_drho[ngpt_itr][j];
						if ((fabs(dx_drho) - 0) < 1e-12)	continue;


						(*cell_info_vector)[cell_itr2].design_points.dxPhys_drho[j] += (qweights[q_point] * dxPhys_dx * dx_drho);

						double dEfactor = dE_dxPhys * dxPhys_dx * dx_drho;
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
						temp_obj_grad[cell_itr2][j] -= dobj;

					}
				}*/

			}

		}
		cell_itr++;
	}

	//Adding data to obj_grad for uncoupled meshes
/*	if (fem->mesh->coupling == false && fem->mesh->adaptivityType != "movingdesignpoints"){
		unsigned int k = 0;
		for(unsigned int i = 0; i < temp_obj_grad.size(); ++i){
			for (unsigned int j = 0; j < temp_obj_grad[i].size(); ++j){
				obj_grad[k] = temp_obj_grad[i][j];
				//std::cout<<k<<"    "<<obj_grad[k]<<std::endl;
				k++;
			}
		}
	}*/

	//std::cout<<obj_grad[0]<<"   "<<obj_grad[24]<<std::endl;
	std::cout<<"Size of sensitivity vector : "<<obj_grad.size()<<std::endl;
	double time2 = clock();
	time2 = (time2 - time1)/(double)CLOCKS_PER_SEC;
	std::cout<<"CPU time for sensitivity analysis: "<<time2<<" seconds"<<std::endl;

/*	for (unsigned int i = 0; i < obj_grad.size(); ++i){
		if (i%10 == 0)
			std::cout<<i<<"  "<<obj_grad[i]<<std::endl;
	}*/
	//exit(0);

	fem->timer->pause();
	//Add details in the log file
	std::ofstream wfile;
	if (fem->itr_count == 0){
		if (fem->cycle == 0){
			wfile.open("output/output.log", std::ios::out);
		}
		else{
			wfile.open("output/output.log", std::ios::app);

		}
		wfile<<"Cycle : "<<fem->cycle+1<<"\n";
		wfile<<"Total Dofs : "<<(*fem).dof_handler.n_dofs()<<"\n";
		wfile<<"Constraints : "<<(*fem).dof_constraints.n_constraints()<<"\n";
		wfile<<"Design count : "<<fem->design_vector->size()<<"\n";
		wfile<<fem->itr_count+1<<"\t"<<objective<<"\t"<<fem->volfrac<<"\n";
		wfile.close();
	}
	else{
		wfile.open("output/output.log", std::ios::app);
		wfile<<fem->itr_count+1<<"\t"<<objective<<"\t"<<fem->volfrac<<"\t"<<(*(fem->cell_info_vector))[0].projection_radius<<"\n";
		wfile.close();
	}
	fem->timer->resume();
}

