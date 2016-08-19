/*
 *
 *  Created on: Aug 14, 2016
 *      Author: Deepak K. Gupta
 *  
 */
#include <atop/TopologyOptimization/adaptivity/dp_adaptivity.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/designField.h>
#include <atop/fem/fem.h>
#include <math.h>

using namespace atop;

template <int dim>
void dpAdaptivity<dim>::update_designField(
		std::vector<CellInfo> &cell_info_vector,
		unsigned int cell_itr,
		unsigned int new_design_count){

	std::vector<std::vector<double> > old_pointX;
	std::vector<double> old_rho;

	old_pointX = cell_info_vector[cell_itr].design_points.pointX;
	old_rho = cell_info_vector[cell_itr].design_points.rho;

	//Updating the design field for this cell
	//Here, we currently have only manual update where points are assigned manually
	//new_design_count = cell_info_vector[cell_itr].design_points.no_points;
	std::cout<<"Hello : "<<new_design_count<<std::endl;
	//cell_info_vector[cell_itr].design_points.no_points = new_design_count;
	cell_info_vector[cell_itr].design_points.update_no_points(new_design_count);

	cell_info_vector[cell_itr].design_points.update_field_manual(0);	//manully assigning the points


	//temporarily assigning the density values
	for (unsigned int i = 0; i < new_design_count; ++i){
		cell_info_vector[cell_itr].design_points.rho[i] = old_rho[0];
	}
}


template <int dim>
void dpAdaptivity<dim>::correctify_p_order(
		FEM<dim> &fem,
		std::vector<CellInfo> &cell_info_vector,
		unsigned int rigid_body_modes){


	this->rigid_body_modes = rigid_body_modes;
	//Updating the element level bounds for all the cells

	//Update p-orders to satisfy the bounds
	unsigned int cell_itr = 0;
	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem.dof_handler->begin_active(),
			endc = fem.dof_handler->end();
	for (; cell != endc; ++cell){
		do{
			//Update the design bound for the cell
			get_corrected_design_bound(fem, cell_info_vector, cell);

			if (cell_info_vector[cell_itr].design_bound < cell_info_vector[cell_itr].design_points.no_points){
				cell_info_vector[cell_itr].shape_function_order++;	//increase the shape order by 1
				cell_info_vector[cell_itr].refine_coarsen_flag = 1;	//now this cell is already refined
			}

		}while(cell_info_vector[cell_itr].design_bound < cell_info_vector[cell_itr].design_points.no_points);
		std::cout<<cell_itr<<"    corrected element bound : "<<cell_info_vector[cell_itr].design_bound<<std::endl;

		cell_itr++;
	}
}

template <int dim>
unsigned int dpAdaptivity<dim>::get_design_bound(
		unsigned int p_order){
	return (pow(p_order + 1, dim) * dim - rigid_body_modes);
}

template <int dim>
void dpAdaptivity<dim>::get_corrected_design_bound(
		FEM<dim> &fem,
		std::vector<CellInfo> &cell_info_vector,
		typename hp::DoFHandler<dim>::active_cell_iterator &cell){

	//Get direct bound on the element based on dofs
	unsigned int cell_itr = cell->user_index() - 1;
	unsigned int design_bound = get_design_bound(cell_info_vector[cell_itr].shape_function_order);

	//Iterate over all the neighbour cells
	for (unsigned int iface = 0; iface < GeometryInfo<dim>::faces_per_cell; ++iface){
		unsigned int ng_shape_fn_order;
		if(cell->at_boundary(iface)) continue;
		if(cell->neighbor(iface)->active()){
			unsigned int ng_cell_itr = cell->neighbor(iface)->user_index() - 1;
			ng_shape_fn_order = cell_info_vector[ng_cell_itr].shape_function_order;
		}

		//Checking the hanging support point for the current cell
		//This needs to be updated for 3D, currently not right for 3D
		unsigned int shape_fn_order = cell_info_vector[cell_itr].shape_function_order;
		if (shape_fn_order <= ng_shape_fn_order)	continue;
		design_bound = design_bound - ((pow(shape_fn_order, dim-1) - pow(ng_shape_fn_order, dim-1))*dim);
	}
	//std::cout<<"cell_itr : "<<design_bound<<std::endl;
	cell_info_vector[cell_itr].design_bound = design_bound;
}


//Getting the system level bound
template <int dim>
unsigned int dpAdaptivity<dim>::get_system_design_bound(
		FEM<dim> &fem){
	unsigned int total_dofs = fem.dof_handler->n_dofs();
	unsigned int no_hanging_nodes = (fem.hanging_node_constraints).n_constraints();
	std::cout<<"No. of hanging nodes : "<<no_hanging_nodes<<std::endl;
	rigid_body_modes = 3;	//hard coded right now for 2d elastostatic problem
	unsigned int sys_bound = total_dofs - no_hanging_nodes - rigid_body_modes;
	std::cout<<"System level bound : "<<sys_bound<<std::endl;
	return sys_bound;

}
