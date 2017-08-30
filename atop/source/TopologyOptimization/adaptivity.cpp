/*
 *
 *  Created on: Aug 22, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/TopologyOptimization/adaptivity.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/physics/elasticity.h>
#include <atop/fem/define_mesh.h>
#include <atop/TopologyOptimization/adaptivity/dp_adaptivity.h>
#include <atop/fem/ErrorIndicator/stressJumpIndicator.h>
#include <atop/fem/ErrorIndicator/qrIndicator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>


#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/component_mask.h>

#include <time.h>


using namespace atop;
using namespace dealii;

template <int dim>
void Adaptivity<dim>::update(
		FEM<dim> &obj_fem){

	//Assigning the FEM object
	this->fem = &obj_fem;
	this->cell_info_vector = (fem->cell_info_vector);
	this->density_cell_info_vector = (fem->density_cell_info_vector);

	//Run the mesh adaptivity approach
	mesh_refine_indicator(fem->mesh->adaptivityType);

}

/**
 * This function chooses the right algorithm based on the provided
 * string to adaptively refine the FE mesh based on some density
 * criteria.
 	 */
template <int dim>
void Adaptivity<dim>::mesh_refine_indicator(
		std::string &update_str){

	if(update_str == "adaptive_grayness"){
		if(fem->mesh->coupling == true){
			coupled_refine_adaptive_grayness();
		}
		else{
			if (fem->mesh->amrType == "dp-refinement"){
				calc_refinement_res_multires();
				compute_sortedRefineRes();
			}
		}

	}
}

template <int dim>
void Adaptivity<dim>::coupled_refine_adaptive_grayness(){


	/*
	 * Note that this routine is not executed in dp-refinement
	 */

	//Getting the cycle number for refinement
	unsigned int cycle = fem->cycle;

	double rhomin = 0.0;
	double rhomax = 1.0;
	double rhomid = (rhomax + rhomin)/2;
	double alpha = 0.2;
	double beta = 0.8;

	double refine_lbound = rhomin + ((1 - alpha) * rhomid * exp(-beta * (double)(cycle+1)));
	double refine_ubound = rhomax - ((1 - alpha) * rhomid * exp(-beta * (double)(cycle+1)));
	double coarsen_lbound = rhomin + (alpha * rhomid * exp(-beta * (double)(cycle+1)));
	double coarsen_ubound = rhomax - (alpha * rhomid * exp(-beta * (double)(cycle+1)));


	typename Triangulation<dim>::active_cell_iterator cell = fem->triangulation.begin_active(),
			endc = fem->triangulation.end();
	typename Triangulation<dim>::active_cell_iterator fe_density_cell = fem->analysis_density_triangulation.begin_active(),
			fe_density_endc = fem->analysis_density_triangulation.end();

	typename Triangulation<dim>::active_cell_iterator density_cell = fem->design_triangulation.begin_active(),
			density_endc = fem->design_triangulation.end();

	unsigned int cell_itr = 0;

	for(; cell != endc; ++cell, ++fe_density_cell, ++density_cell){

		unsigned int n_qpoints = (*cell_info_vector)[cell_itr].density.size();
		double density = 0.0;

		for(unsigned int qpoint = 0; qpoint < n_qpoints; ++qpoint){
			density += (*cell_info_vector)[cell_itr].density[qpoint] *
								(*cell_info_vector)[cell_itr].density_weights[qpoint];

		}
		if(density > coarsen_ubound || density < coarsen_lbound){
			cell->set_coarsen_flag();
			fe_density_cell->set_coarsen_flag();
			density_cell->set_coarsen_flag();
		}
		else if (density > refine_lbound && density < refine_ubound) {
			cell->set_refine_flag();
			fe_density_cell->set_refine_flag();
			density_cell->set_refine_flag();
		}

		//Setting user index of the parent of cell which is selected for coarsening
		if (cell->coarsen_flag_set() && cell->level() > 0){
			cell->parent()->set_user_index(cell->user_index());
		}
		if (fe_density_cell->coarsen_flag_set() && fe_density_cell->level() > 0){
			fe_density_cell->parent()->set_user_index(fe_density_cell->user_index());
		}
		if (density_cell->coarsen_flag_set() && density_cell->level() > 0){
			density_cell->parent()->set_user_index(density_cell->user_index());
		}

		++cell_itr;
	}
}

/**
 * This refinement residual function is implemented for meshes where more than one design points are coupled per element
 * It extends the adaptive grayness proposed in Gupta et al 2016 and checks the values for all the design variables.
 * It does not use -1, 0, 1, rather it gives a residual value which is used for preferring the elements for refinement.
 */
template <int dim>
void Adaptivity<dim>::calc_refinement_res_multires(){

	std::cout<<"Computing the refinement residuals.......";
	//Getting the cycle number for refinement
	unsigned int cycle = fem->cycle;

	//Defining the bounds to characterize the refinement/coarsening zones
	double rhomin = 0.0;
	double rhomax = 1.0;
	double rhomid = (rhomax + rhomin)/2;
	double alpha = 0.2;
	double beta = 1.4;

	double refine_lbound = rhomin + ((1 - alpha) * rhomid * exp(-beta * (double)(cycle+1)));
	double refine_ubound = rhomax - ((1 - alpha) * rhomid * exp(-beta * (double)(cycle+1)));
	double coarsen_lbound = rhomin + (alpha * rhomid * exp(-beta * (double)(cycle+1)));
	double coarsen_ubound = rhomax - (alpha * rhomid * exp(-beta * (double)(cycle+1)));
/*	refine_lbound = 0.05;
	refine_ubound = 0.95;
	coarsen_lbound = refine_lbound;
	coarsen_ubound = refine_ubound;*/

	std::cout<<"Refinement bounds : "<<refine_lbound<<"   "<<refine_ubound<<std::endl;

	//Initializing the refineRes vector
	refineRes.clear();
	refineRes.resize(cell_info_vector->size(), 0.0);

	//Iterating over the cell_info_vector

	for (unsigned int cell_itr = 0; cell_itr < cell_info_vector->size(); ++cell_itr){

		double density = (*cell_info_vector)[cell_itr].cell_density;
		if(density > coarsen_ubound){
			refineRes[cell_itr] += (density - 1.0);
		}
		if(density < coarsen_lbound){
			refineRes[cell_itr] += (0.0 - density);
		}
		if(density >= refine_lbound && density <= rhomid){
			refineRes[cell_itr] += density - refine_lbound;
		}
		if(density <= refine_ubound && density > rhomid){
			refineRes[cell_itr] += refine_ubound - density;
		}
	}
	std::cout<<"DONE"<<std::endl;
}



template <int dim>
void Adaptivity<dim>::update_cell_vectors(
		std::vector<CellInfo> &density_cell_info_vector,
		hp::DoFHandler<dim> &density_dof_handler,
		Triangulation<dim> &density_triangulation){
	std::cout<<"Updating density info vector record"<<std::endl;

	if (fem->mesh->amrType == "dp-refinement"){

		return;
	}
	std::vector<CellInfo> temp_cellprop;
	temp_cellprop.clear();
	temp_cellprop = density_cell_info_vector;
	density_cell_info_vector.clear();
	density_cell_info_vector.resize(density_triangulation.n_active_cells());
	unsigned int density_cell_itr = 0;
	typename hp::DoFHandler<dim>::active_cell_iterator density_cell = density_dof_handler.begin_active(),
			density_endc = density_dof_handler.end();
	for(; density_cell != density_endc; ++density_cell){

		if(density_cell->level() > 0 && density_cell->parent()->user_index() > 0){
			density_cell_info_vector[density_cell_itr] = temp_cellprop[density_cell->parent()->user_index() - 1];
		}

		if(density_cell->user_index() > 0){
			density_cell_info_vector[density_cell_itr] = temp_cellprop[density_cell->user_index() - 1];
		}

		//Updating the cell area
		density_cell_info_vector[density_cell_itr].cell_area = density_cell->measure();

		++density_cell_itr;
	}

		//clearing temp_cellprop
		temp_cellprop.clear();

	//Clearing the user_index for all cells
	density_cell = density_dof_handler.begin_active();
	density_cell_itr = 0;
	for(; density_cell != density_endc; ++density_cell){
		if(density_cell->level() > 0 && density_cell->parent()->user_index() > 0){
			density_cell->parent()->clear_user_index();
		}
		if(density_cell->user_index() > 0){
			density_cell->clear_user_index();
		}
	}
}

//This function is used for calling the h-refinement module of deal.II or perform p-refinement
template <int dim>
void Adaptivity<dim>::execute_coarsen_refine(){

	if (fem->mesh->amrType == "h-refinement"){
		fem->triangulation.execute_coarsening_and_refinement();
		fem->analysis_density_triangulation.execute_coarsening_and_refinement();
		fem->design_triangulation.execute_coarsening_and_refinement();

		//std::cout<<"No. of cells after refinement : "<<fem->triangulation->n_active_cells()<<std::endl;
	}
	else if (fem->mesh->amrType == "dp-refinement"){

			system_design_bound = dp_adap.get_system_design_bound(*fem);
			if (dim == 2)	rigid_body_modes = 3;	//manually assigning the RBMs
			update_element_design_bound();
			//dp_coarsening_refinement();	//dp-adaptivity performed
			improved_dp_coarsening_refinement();
	}
}


// To calculate the sortedRefineRes vector
template <int dim>
void Adaptivity<dim>::compute_sortedRefineRes(){

	unsigned int len = refineRes.size();
	sortedRefineRes.resize(len, {0.0, 0});

	//Initialize sortedRefineRres
	for (unsigned int i = 0; i < len; ++i){
		sortedRefineRes[i].first = refineRes[i];
		sortedRefineRes[i].second = i;
	}

	//sorting the residuals
	for (unsigned int i = 0; i < len; ++i){
		for (unsigned int j = 0; j < len-1; ++j){
			if (sortedRefineRes[j].first > sortedRefineRes[j+1].first){
				double tempres = sortedRefineRes[j].first;
				unsigned int tempi = sortedRefineRes[j].second;
				sortedRefineRes[j].first = sortedRefineRes[j+1].first;
				sortedRefineRes[j].second = sortedRefineRes[j+1].second;
				sortedRefineRes[j+1].first = tempres;
				sortedRefineRes[j+1].second = tempi;
			}
		}
	}
}

//Updating the element level bounds
template <int dim>
void Adaptivity<dim>::update_element_design_bound(){
	if (fem->mesh->amrType == "dp-refinement"){

		unsigned int cell_itr = 0;	//Iterator for the triangulation vector
		typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
				endc = fem->dof_handler.end();
		for (; cell != endc; ++cell){
			unsigned int design_bound = (*cell_info_vector)[cell_itr].dofs_per_cell - rigid_body_modes;

			//Iterate over all the neighbour cells
			for (unsigned int iface = 0; iface < GeometryInfo<dim>::faces_per_cell; ++iface){
				unsigned int ng_shape_fn_order;
				if(cell->at_boundary(iface)) continue;
				if(cell->neighbor(iface)->active()){
					unsigned int ng_cell_itr = cell->neighbor(iface)->user_index() - 1;
					ng_shape_fn_order = (*cell_info_vector)[ng_cell_itr].shape_function_order;
				}

				//Checking the hanging support point for the current cell
				//This needs to be updated for 3D, currently not right for 3D
				unsigned int shape_fn_order = (*cell_info_vector)[cell_itr].shape_function_order;
				if (shape_fn_order <= ng_shape_fn_order)	continue;
				design_bound = design_bound - ((pow(shape_fn_order, dim-1) - pow(ng_shape_fn_order, dim-1))*dim);
			}
			//std::cout<<"cell_itr : "<<design_bound<<std::endl;
			(*cell_info_vector)[cell_itr].design_bound = design_bound;
			++cell_itr;
		}
	}
}


//Perform adaptive refinement/coarsening for dp-refinement strategy
template <int dim>
void Adaptivity<dim>::dp_coarsening_refinement(){

	run_dp_analysis_based_refinement();

   std::cout<<"Performing design based refinement "<<std::endl;
	unsigned int no_cells = sortedRefineRes.size();

/*
 * In the code below, the refinement residuals are used to perform dp-refinement
 * For -ve values, cells are considered for cooarsening
 * For +ve values, cells are considered for refinement
 * The refinement/coarsening are based on d-refinement and p is adjusted accordingly
 */

	std::vector<unsigned int > new_no_design_vector((*cell_info_vector).size(), 0);
	//Iterate over all the cells
	for (unsigned int i = 0; i < no_cells; ++i){

		unsigned int cell_itr = sortedRefineRes[i].second;
		unsigned int current_p_order = (*cell_info_vector)[cell_itr].shape_function_order;
		unsigned int current_no_design = (*cell_info_vector)[cell_itr].design_points.no_points;
		unsigned int new_p_order, new_no_design;
		//Coarsening
		if (fabs(sortedRefineRes[i].first) < 1e-14){
			new_p_order = (*cell_info_vector)[cell_itr].shape_function_order;
			new_no_design = (*cell_info_vector)[cell_itr].design_points.no_points;

			//std::cout<<"reached in this one"<<std::endl;
		}
		else{
			if (sortedRefineRes[i].first < 0){
				if (/*current_p_order == 1 && */ current_no_design == 1){
					new_no_design = (*cell_info_vector)[cell_itr].design_points.no_points;
				}
				else{
					new_no_design = pow(floor(sqrt((double)(current_no_design - 1))), 2);	//written for 2D
					if (new_no_design > 1){
						//new_no_design = pow(floor(sqrt((double)(new_no_design - 1))), 2);	//written for 2D
					}
					//(*cell_info_vector)[cell_itr].refine_coarsen_flag = -1;
				}
			}
			else if (sortedRefineRes[i].first > 0){
				new_no_design = pow(ceil(sqrt((double)(current_no_design + 1))), 2);	//written for 2D
				//new_no_design = pow(ceil(sqrt((double)(new_no_design + 1))), 2);	//written for 2D
				//(*cell_info_vector)[cell_itr].refine_coarsen_flag = 1;
			}

			unsigned int min_dofs  = new_no_design + rigid_body_modes;
			new_p_order = ceil(pow((min_dofs/(double)dim), 1/((double)dim)) - 1);
			if ((*cell_info_vector)[cell_itr].refine_coarsen_flag == 1){
				if (new_p_order < current_p_order){
					new_p_order = current_p_order;	//To make sure cells refined during analysis refinement are not coarsened here
				}
			}

			//std::cout<<"New p : "<<new_p_order<<"    new_no_design : "<<new_no_design<<std::endl;
			(*cell_info_vector)[cell_itr].shape_function_order = new_p_order;
			if (new_p_order <= 0){
				std::cerr<<"Zero shape function order found in Adaptivity class"<<std::endl;
				exit(0);
			}
			//Interpolate to the new design field for the current cell

			dp_adap.update_designField(*cell_info_vector,
					cell_itr,
					new_no_design);
		}

	}
//----------------------------------------------------------------------------------------------------------------------------------

	/*
	 * finding the maximum and minimum order of d field and p field
	 * This information is used to reduce the p- and d-contrasts in the active refinement areas
	 */


	unsigned int min_d_points = 999, max_d_points = 0, min_p_order = 999, max_p_order = 0;
	for (unsigned int i = 0; i < cell_info_vector->size(); ++i){
		if ((*cell_info_vector)[i].design_points.no_points < min_d_points)	min_d_points = (*cell_info_vector)[i].design_points.no_points;
		if ((*cell_info_vector)[i].design_points.no_points > max_d_points)	max_d_points = (*cell_info_vector)[i].design_points.no_points;
		if ((*cell_info_vector)[i].shape_function_order < min_p_order)	min_p_order = (*cell_info_vector)[i].shape_function_order;
		if ((*cell_info_vector)[i].shape_function_order > max_p_order)	max_p_order = (*cell_info_vector)[i].shape_function_order;


	}

	for (int i = round(sqrt(min_d_points)); i <= int(round(sqrt(max_d_points))-2); ++i){
		std::cout<<i<<"      Regularizing design field ..."<<std::endl;
		dp_adap.update_design_contrast(*fem, *cell_info_vector, rigid_body_modes);
	}

	for (int i = min_p_order; i <= (int)(max_p_order-2); ++i){
		std::cout<<i<<"  "<<max_p_order<<"      Repairing the shape functions ..."<<std::endl;
		if (i > 12) exit(0);
		dp_adap.correctify_p_order(*fem, *cell_info_vector, rigid_body_modes);
	}

	std::cout<<"Polishing the p-distribution for 1-level hanging "<<std::endl;
	dp_adap.update_p_order_contrast(*fem, *cell_info_vector);

	unsigned int cell_itr = 0;
	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	for (; cell != endc; ++cell){
		//std::cout<<"Updated p order : "<<(*cell_info_vector)[cell_itr].shape_function_order<<std::endl;
		unsigned int p_index = ((fem->elastic_data)).get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		cell->set_active_fe_index(p_index);
		cell_itr++;
	}
	//Correcting system-level violations
	std::cout<<"Correcting system level violations ..."<<std::endl;
	//unsigned int temp = dp_adap.get_corrected_system_design_bound(*fem, *cell_info_vector);
	std::cout<<"System level violations corrected "<<std::endl;

}

template <int dim>
void Adaptivity<dim>::improved_dp_coarsening_refinement(){

	//Analysis based refinement
	run_dp_analysis_based_refinement();
	//------------------------------------------------------------------

	std::cout<<"Initiating design based refinement module ....."<<std::endl;
/*
 * In the code below, the refinement residuals are used to perform dp-refinement
 * For -ve values, cells are considered for coarsening
 * For +ve values, cells are considered for refinement
 */

	//Iterate over all the cells
	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin(),
			endc = fem->dof_handler.end();
	unsigned int cell_itr = 0;
	//Updating the p-orders based on the refineRes values
	for (; cell != endc; ++cell){

			int current_p_order = (*cell_info_vector)[cell_itr].shape_function_order;
			int current_no_design = (*cell_info_vector)[cell_itr].design_points.no_points;
			int new_p_order = current_p_order;	//temporarily set to the current value

			if (fabs(refineRes[cell_itr]) < 1e-14){ // not flagged for refinement / coarsening
				//Here refineRes refers to the residual based on the density indicator (currently Gupta et al, 2016 indicator)
				//Do not update the d-value, therefore p-order is also not changed here.
			}
			else{
				int diff;
				int lower_design_bound = dp_adap.get_design_bound(current_p_order - 1);
				int upper_design_bound = dp_adap.get_design_bound(current_p_order);

				//Coarsening
				if (refineRes[cell_itr] < -1e-12){
					if ((*cell_info_vector)[cell_itr].refine_coarsen_flag != 1){	//if p not increased during analysis based check
						(*cell_info_vector)[cell_itr].refine_coarsen_flag = -1;	// coarsen design field, and adjust p based on bound check
						if (dim == 2 && current_p_order > 1){
							int current_dfactor = ceil(sqrt(current_no_design));
							diff = 1 + 2 * (current_dfactor);

							if ((*cell_info_vector)[cell_itr].old_shape_fn_order < current_p_order){//((current_no_design - diff) > lower_design_bound){
								new_p_order = current_p_order;
								(*cell_info_vector)[cell_itr].refine_coarsen_flag = -2;
							}
							else{
								new_p_order = current_p_order - 1;
								(*cell_info_vector)[cell_itr].refine_coarsen_flag = -2; //coarsened during dp-coarsening
							}
						}
						/*
						 * In the above check, the case with  p = 1 is ignored. For this case coarsening of p is not possible.
						 * However, if for this case, if old_shape_fn_order and shape_function_order are same, it means the cell was
						 * not coarsened for analysis check as well as for density based check. For such a case, the no. of design
						 * variables for this cell will be set to 1.
						 */

					}
				}
				else if (refineRes[cell_itr] > 1e-12){
					int current_dfactor = ceil(sqrt(current_no_design));

					//calculating ceil d_Factor
					if (dim == 2){
						diff = 1 + 2*current_dfactor;
					}
				    //int corrected_design_bound = dp_adap.get_corrected_design_bound(*fem, *cell_info_vector, cell);
					if ((*cell_info_vector)[cell_itr].refine_coarsen_flag == 1){//current_no_design + diff <= corrected_design_bound){
						new_p_order = current_p_order;
						(*cell_info_vector)[cell_itr].refine_coarsen_flag = 2;
					}
					else{
						new_p_order = (*cell_info_vector)[cell_itr].old_shape_fn_order + 1;
						(*cell_info_vector)[cell_itr].refine_coarsen_flag = 2;	// p-order increased during dp-refinement
					}
				}
			}
			(*cell_info_vector)[cell_itr].shape_function_order = new_p_order;
			if (new_p_order <= 0){
				std::cerr<<"Zero shape function order found in Adaptivity class"<<std::endl;
				exit(0);
			}

			//std::cout<<"refine flag : "<<(*cell_info_vector)[cell_itr].refine_coarsen_flag<<"  "<<current_p_order<<"    "<<new_p_order<<std::endl;
			cell_itr++;

	}
//----------------------------------------------------------------------------------------------------------------------------------

	/*
	 * finding the maximum and minimum order of d field and p field
	 * This information is used to reduce the p- and d-contrasts in the active refinement areas
	 */


	unsigned int min_p_order = 999, max_p_order = 0;
	for (unsigned int i = 0; i < cell_info_vector->size(); ++i){
		if ((*cell_info_vector)[i].shape_function_order < min_p_order)	min_p_order = (*cell_info_vector)[i].shape_function_order;
		if ((*cell_info_vector)[i].shape_function_order > max_p_order)	max_p_order = (*cell_info_vector)[i].shape_function_order;


	}

	std::cout<<"Polishing the p-distribution for 1-level hanging "<<std::endl;
	dp_adap.update_p_order_contrast(*fem, *cell_info_vector);


	//Update the design field to allow maximum number of permissible design variables as per element-bound in each element.
	dp_adap.update_design_for_elem_bound_only(
					*fem,
					*cell_info_vector);
/*
	//Correcting system-level violations
	std::cout<<"Correcting system level violations ..."<<std::endl;
	//unsigned int temp = dp_adap.get_corrected_system_design_bound(*fem, *cell_info_vector);
	std::cout<<"System level violations corrected "<<std::endl;
*/

	// Update the p-order to reduce the qr-patterns based on solution of previous cycle
	run_qr_based_refinement();

	cell_itr = 0;
	cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	for (; cell != endc; ++cell){
		//std::cout<<"Updated p order : "<<(*cell_info_vector)[cell_itr].shape_function_order<<std::endl;
		unsigned int p_index = ((fem->elastic_data)).get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		cell->set_active_fe_index(p_index);
		cell_itr++;
	}

}


//function to update the polynomial order of the cells using the analysis indicator
template <int dim>
void Adaptivity<dim>::increase_decrease_p_order(){
	unsigned int cell_itr = 0;	//Iterator for the triangulation vector
	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	for (; cell != endc; ++cell){
		(*cell_info_vector)[cell_itr].refine_coarsen_flag = 0;

		if (cell->refine_flag_set()){
			if ((*cell_info_vector)[cell_itr].shape_function_order < (*(fem->mesh)).max_el_order){
				// The above condition is put to ensure that p-values higher than max_el_order are not allowed
				(*cell_info_vector)[cell_itr].shape_function_order++;
				(*cell_info_vector)[cell_itr].refine_coarsen_flag = 1;
			}
		}

		if (cell->coarsen_flag_set()){
			if ((*cell_info_vector)[cell_itr].shape_function_order > 1){
				(*cell_info_vector)[cell_itr].shape_function_order--;
				(*cell_info_vector)[cell_itr].refine_coarsen_flag = -1;
			}
		}

		++cell_itr;
	}
}

template <int dim>
void Adaptivity<dim>::run_dp_analysis_based_refinement(){

	for (unsigned int cell_itr = 0; cell_itr < (*cell_info_vector).size(); ++cell_itr){
		//Line below saves a history for using it in qr-refinement
		(*cell_info_vector)[cell_itr].old_shape_fn_order = (*cell_info_vector)[cell_itr].shape_function_order;
		(*cell_info_vector)[cell_itr].old_design_count = (*cell_info_vector)[cell_itr].design_points.no_points;
	}
	if (fem->cycle < 2){
		//Update shape functions based on analysis error criterion
		cout<<"Performing analysis based refinement : "<<std::endl;

		std::vector<double> estimated_error_per_cell (fem->triangulation.n_active_cells());
		StressJumpIndicator<dim> modKellyObj(*fem,
				estimated_error_per_cell,
				*cell_info_vector);

		modKellyObj.estimate();


		Vector<double> error_indicator(estimated_error_per_cell.size());
		for (unsigned int i = 0; i < error_indicator.size(); ++i)
			error_indicator(i) = estimated_error_per_cell[i];

		GridRefinement::refine_and_coarsen_fixed_number (fem->triangulation,
														   error_indicator,
														   0.1, 0.05);

	   increase_decrease_p_order();	//refined/coarsened
	   std::cout<<"Analysis based refinement done "<<std::endl;
	}

	for (unsigned int i = 0; i < (*cell_info_vector).size(); ++i){
		if ((*cell_info_vector)[i].shape_function_order == 0){
			std::cerr<<"Zero shape function order found here\n";
			exit(0);
		}
	}

}

template <int dim>
void Adaptivity<dim>::run_qr_based_refinement(){

	if (fem->cycle >= 2)	return;

	std::vector<double> qr_accuracy(fem->triangulation.n_active_cells()); //to store accuracy of solution for each element
	std::vector<unsigned int> proposed_p_values(fem->triangulation.n_active_cells());	// obtained based on qr-check

	// Function calls below does the complete QR-test
	QRIndicator<dim> qr_test(*fem,
			qr_accuracy,
			0.8,
			proposed_p_values,
			*cell_info_vector
			);
	qr_test.estimate(qr_accuracy);
/*	for (unsigned int i = 0; i < qr_accuracy.size(); ++i){
		std::cout<<i<<"  "<<qr_accuracy[i]<<std::endl;
	}*/

	//Refine all cells which have not been refined and have qr_accuracy less than 0.1
	unsigned int cell_itr = 0;	//Iterator for the triangulation vector
	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	for (; cell != endc; ++cell){
		(*cell_info_vector)[cell_itr].refine_coarsen_flag = 0;
		if (qr_accuracy[cell_itr] > 0.1){
			++cell_itr;
			continue;
		}
		if ((*cell_info_vector)[cell_itr].shape_function_order <= (*cell_info_vector)[cell_itr].old_shape_fn_order){
			(*cell_info_vector)[cell_itr].shape_function_order = (*cell_info_vector)[cell_itr].old_shape_fn_order + 1;
			(*cell_info_vector)[cell_itr].refine_coarsen_flag = 3; //QR-refinement has been done in this cell
		}
		++cell_itr;
	}
}
