/*
 *
 *  Created on: Jun 18, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <string>
#include <iostream>
#include <atop/fem/fem.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/fem/boundary_values.h>
#include <atop/fem/output.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <atop/fem/define_mesh.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/full_matrix.h>


#include<vector>
#include <fstream>




using namespace atop;
using namespace dealii;

//Constructor for initialization
template<int dim>
FEM<dim>::FEM(
		std::vector<CellInfo> &obj_cell_info_vector,
		std::vector<CellInfo> &obj_density_cell_info_vector,
		DefineMesh<dim> &obj_mesh,
		std::vector<double> &obj_design_vector,
		Timer &obj_timer,
		std::string &problem_name):
		dof_handler(triangulation),
		analysis_density_handler(analysis_density_triangulation),
		design_handler(design_triangulation){

	this->cell_info_vector = &obj_cell_info_vector;
	this->density_cell_info_vector = &obj_density_cell_info_vector;
	this->mesh = &obj_mesh;
	this->timer = &obj_timer;
	this->problem_name = problem_name;

	//Create the mesh
	this->mesh->createMesh(
			triangulation,
			analysis_density_triangulation,
			design_triangulation);

	if (this->problem_name == "electrical_conduction"){
		std::cout<<" Defining a scalar problem"<<std::endl;
		// for scalar elements
		//For Lagrange element
		if (mesh->elementType == "FE_Q"){
			for (unsigned int degree = 1; degree <= mesh->max_el_order; ++degree){
				fe_collection.push_back(FE_Q<dim>(degree));
			}

		}

		//For Legendre elements (not tested yet)
		if (mesh->elementType == "FE_Q_hierarchical"){
			for (unsigned int degree = 1; degree <= mesh->max_el_order; ++degree){
				fe_collection.push_back(FE_Q_Hierarchical<dim>(degree));
			}
		}

	}

	else{
		//For Lagrange element
		if (mesh->elementType == "FE_Q"){
			for (unsigned int degree = 1; degree <= mesh->max_el_order; ++degree){
				fe_collection.push_back(FESystem<dim>(FE_Q<dim>(degree), dim));
			}

		}

		//For Legendre elements (not tested yet)
		if (mesh->elementType == "FE_Q_hierarchical"){
			for (unsigned int degree = 1; degree <= mesh->max_el_order; ++degree){
				fe_collection.push_back(FESystem<dim>(FE_Q_Hierarchical<dim>(degree), dim));
			}
		}
	}


	//Quadrature collection for FE
	for (unsigned int qrule = 1; qrule <= mesh->max_el_order+7; ++qrule){
		quadrature_collection.push_back(QGauss<dim>(qrule));
		face_quadrature_collection.push_back(QGauss<dim-1>(qrule));
	}

	//For density elements assuming they are DGQ type
	if (mesh->density_elementType == "FE_DGQ"){
		for (unsigned int degree = 1; degree <= mesh->max_density_el_order; ++degree){
			fe_analysis_density_collection.push_back(FESystem<dim>(FE_DGQ<dim>(degree), 1));
			fe_design_collection.push_back(FESystem<dim>(FE_DGQ<dim>(degree), 1));
		}
	}
	this->design_vector = &obj_design_vector;

	unsigned int cell_itr = 0;
	unsigned int p_index = elastic_data.get_p_index(mesh->initial_el_order);
	//intializing the type of element for each cell of analysis mesh
	for (typename hp::DoFHandler<dim>::active_cell_iterator cell =dof_handler.begin_active();
			cell != dof_handler.end(); ++cell){
		cell->set_active_fe_index(p_index);

		cell_itr++;

	}
	//intializing the type of element for each cell of analysis density mesh
	for (typename hp::DoFHandler<dim>::active_cell_iterator cell = analysis_density_handler.begin_active();
			cell != analysis_density_handler.end(); ++cell){
		cell->set_active_fe_index(0);
	}

	//intializing the type of element for each cell of design mesh
	for (typename hp::DoFHandler<dim>::active_cell_iterator cell = design_handler.begin_active();
			cell != design_handler.end(); ++cell){
		cell->set_active_fe_index(0);
	}
}

template <int dim>
FEM<dim>::~FEM(){

}

//Function that step-by-step solves the FE problem
template <int dim>
void FEM<dim>::analyze(){
	//std::cout<<"Entered FEM::analyze()"<<std::endl;
	//Setting up and assembling the Fe system



	itr_count++;	//Iterating the counter for the no. of iterations
	clean_trash();
	//boundary_info();
	std::cout<<"Entering FEM::setup_system()"<<std::endl;

	setup_system();
	assemble_system();
	//std::cout<<"System assembled"<<std::endl;
	solve();	//Solving the system
	//std::cout<<"System solved"<<std::endl;
	output_results();
}

//Function for setting up the FE system
template <int dim>
void FEM<dim>::setup_system(){

	std::cout<<"Started setting up the FEM system ...."<<std::endl;
	boundary_values.clear();
	//FE mesh
	dof_handler.distribute_dofs(fe_collection);
	hanging_node_constraints.clear();
	DoFTools::make_hanging_node_constraints(dof_handler,
			hanging_node_constraints);

	//Applying the boundary conditions
	boundary_info();

	add_boundary_constraints();	//For adding complex boundary related constraints
	BoundaryValues<dim> boundary_v;

	VectorTools::interpolate_boundary_values(dof_handler,
			42,
			/*ZeroFunction<dim>(),*/boundary_v,
			hanging_node_constraints);//boundary_values);

	hanging_node_constraints.close();
	  std::cout<< "No.of degrees of freedom: " << dof_handler.n_dofs() << "\n";
	std::cout<<"No. of hanging node constraints : "<<hanging_node_constraints.n_constraints()<<std::endl;
	analysis_density_handler.distribute_dofs(fe_analysis_density_collection);	//Used to add density on every node

	//Density mesh or design mesh
	design_handler.distribute_dofs(fe_design_collection);

	DynamicSparsityPattern dsp (dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler,
	                                 dsp,
	                                 hanging_node_constraints,
	                                 false);
	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);

	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
	lambda_solution.reinit(dof_handler.n_dofs());
	l_vector.reinit(dof_handler.n_dofs());
	nodal_density.reinit(analysis_density_handler.n_dofs());	//filtered densities for the output
	nodal_p_order.reinit(analysis_density_handler.n_dofs());
	cells_adjacent_per_node.reinit(analysis_density_handler.n_dofs());	//for normalizing the nodal density value
	nodal_d_count.reinit(analysis_density_handler.n_dofs());
}

/**
 * This function assembles the contributions from each element.
 * It creates the system matrix and the global RHS
 */
template <int dim>
void FEM<dim>::assemble_system(){
	std::cout<<"Assembling the system..."<<std::endl;

	//Initialize cell_info and density_cell_info objects
	if(cycle == 0 && itr_count == 0){
		reset();	//Initializes the variables and vectors
	}


/*
 * At every iteration, the content of design_vector needs to be copied into the
 * cell_info_vector.
 * For first iteration of every cycle from 2 to n, adaptive refinement might take place.
 * The updated content is in cell_info_vector, thus in this step, we do not perform this step.
 */
	if (!(cycle > 0  && itr_count == 0)){
		//Update the density_cell_info_vector
		if (mesh->coupling == false && mesh->adaptivityType != "movingdesignpoints"){
			density_field.update_density_cell_info_vector(
					*cell_info_vector,
					*density_cell_info_vector,
					*design_vector);
		}
		else{
			density_field.update_density_cell_info_vector(
					*density_cell_info_vector,
					*design_vector);
		}
	}

	//Initialize the components of every cycle
	if (itr_count == 0){
		initialize_cycle();
	}
	else{

	}

	//update the pseudo-design field
    update_pseudo_designField();

    //Manually updating the pseudo_densities
/*    if (itr_count == 0){
        for (unsigned int cell_itr = 0; cell_itr < (*cell_info_vector).size(); ++cell_itr){
        	if (cell_itr != 1){
        		for (unsigned int dpoint = 0; dpoint < (*cell_info_vector)[cell_itr].design_points.no_points; ++dpoint){
        			(*cell_info_vector)[cell_itr].design_points.rho[dpoint] = 1.0;
        			(*cell_info_vector)[cell_itr].pseudo_design_points.rho[dpoint] = 1.0;
        		}
        	}
        	else{
        		for (unsigned int dpoint = 0; dpoint < (*cell_info_vector)[cell_itr].design_points.no_points; ++dpoint){
        			(*cell_info_vector)[cell_itr].design_points.rho[dpoint] = 0.55;
        			(*cell_info_vector)[cell_itr].pseudo_design_points.rho[dpoint] = 0.55;

        		}
        	}
        }
    }*/


	//Apply smoothing operation on the density values
	density_field.smoothing(*cell_info_vector,
			*density_cell_info_vector,
			*mesh);
	std::cout<<"Smoothing done"<<std::endl;
	OutputData<dim> out_soln;
	if (fileReadFlag == true){
		out_soln.read_xPhys_from_file(*cell_info_vector,
					filefname);
	}
/*	OutputData<dim> out_soln;
	out_soln.read_xPhys_from_file(*cell_info_vector,
			"output_design/density_1_46.dat");*/



	//Updating the physics of the problem

	update_physics();
	std::cout<<"Physics updated"<<std::endl;

	//Compute cellwise material properties
	if (this->problem_name == "electrical_conduction")
		penal->update_param(linear_electrostatic->E0,
				linear_electrostatic->Emin,
				*cell_info_vector);
	else
		penal->update_param(linear_elastic->E, *cell_info_vector);
	std::cout<<"Cell parameters updated"<<std::endl;

	//Assembling the system and RHS
	assembly();

}

template<int dim>
void FEM<dim>::problemType(LinearElastic<dim> &obj_linear_elastic){
	this->linear_elastic = &obj_linear_elastic;
}

template<int dim>
void FEM<dim>::problemType(LinearElectrostatic<dim> &obj_linear_elastostatic){
	this->linear_electrostatic = &obj_linear_elastostatic;
}
template<int dim>
void FEM<dim>::projectionType(Projection &projection){
	this->projection = &projection;
}

template <int dim>
void FEM<dim>::penalization(Penalize &obj_penal){
	this->penal = &obj_penal;
}

template <int dim>
void FEM<dim>::solve(){

/*	SolverControl solver_control(50000, 1e-10);
	SolverCG<> cg(solver_control);
	PreconditionSSOR<> preconditioner;
	preconditioner.initialize(system_matrix, 1.5);
	cg.solve(system_matrix,
			solution,
			system_rhs,
			preconditioner);*/

/*	for (unsigned int i = 0; i < system_rhs.size(); ++i)
		std::cout<<i<<"  "<<system_rhs(i)<<std::endl;*/

	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (solution, system_rhs);

	if (self_adjoint == false){
		std::cout<<"Calculating lambda "<<std::endl;
		A_direct.vmult(lambda_solution, const_vector);
		lambda_solution *= -1;
	}
	hanging_node_constraints.distribute(solution);
	hanging_node_constraints.distribute(lambda_solution);

	//updating the system_rhs for the dirichlet boundaries
	//system_matrix.vmult(system_rhs, solution);
	std::cout<<"FEM system solved "<<std::endl;

/*	double sumQ = 0.0;
	for (unsigned int i = 0; i < system_rhs.size(); ++i){
		sumQ += system_rhs[i];
	}
	double Power = sumQ * 0.5;
	std::cout<<"Power  : "<<Power<<std::endl;
	exit(0);*/
}

template <int dim>
void FEM<dim>::output_results(){

	//Writing the state solution
	std::string filename = "solution-";
	std::stringstream ss;
	ss<< cycle +1<<"_"<<itr_count+1;
	filename += ss.str();
	filename += ".vtk";
	std::vector<std::string> solution_names;
	solution_names.clear();
	if (this->problem_name == "electrical_conduction"){
		solution_names.push_back("temperature");
	}
	else{
		switch(dim){
		case 1:
			solution_names.push_back("displacement");
			break;
		case 2:
			solution_names.push_back("x_displacement");
			solution_names.push_back("y_displacement");
			break;
		case 3:
			solution_names.push_back("x_displacement");
			solution_names.push_back("y_displacement");
			solution_names.push_back("z_displacement");
			break;
		default:
			Assert(false, ExcNotImplemented);
		}
	}

	OutputData<dim> out_soln;
	out_soln.write_fe_solution(filename, dof_handler,
			solution, solution_names);

	//Writing the density solution
	filename = "density-";
	filename += ss.str();
	filename += ".vtk";
	std::vector<std::string> density_names;
	density_names.clear();
	density_names.push_back("density");
	out_soln.write_fe_solution(filename, analysis_density_handler,
			nodal_density, density_names);

	//Writing the shape order
	if (itr_count == 1){
		filename = "shape-";
		filename += ss.str();
		filename += ".vtk";
		std::vector<std::string> density_names;
		density_names.clear();
		density_names.push_back("p-order");
		out_soln.write_fe_solution(filename, analysis_density_handler,
				nodal_p_order, density_names);

		filename = "d_count-";
		filename += ss.str();
		filename += ".vtk";
		density_names.clear();
		density_names.push_back("d-count");
		out_soln.write_fe_solution(filename, analysis_density_handler,
				nodal_d_count, density_names);
	}
	//Writing the design output
	filename = "design-";
	filename += ss.str();
	filename += ".dat";

	if(mesh->adaptivityType == "movingdesignpoints"){
		out_soln.write_design(filename,
				*design_vector,
				mesh->design_var_per_point());
	}
	else{
		out_soln.write_design(filename,
				dof_handler,
				*cell_info_vector);
	}

	//std::cout<<"Output files written"<<std::endl;

}
template <int dim>
void FEM<dim>::reset(){
	//initializing the sizes of current and running quad rule
	gauss_int.initialize_quadRuleVector(current_quad_rule,
			mesh->max_el_order, mesh->initial_dcount_per_el);

	gauss_int.max_el_order = mesh->max_el_order;
	running_quad_rule.clear();
	running_quad_rule.resize(mesh->max_el_order, 1);

	if (this->problem_name == "electrical_conduction")
		electrostatic_data.initialize_quadRuleVectors(current_quad_rule,
						running_quad_rule);
	else{
		elastic_data.initialize_quadRuleVectors(current_quad_rule,
				running_quad_rule);
		elastic_data.nu = linear_elastic->poisson;
	}

	//current_quad_rule  will be updated on adaptivity in quadrature
	/**
	 * Initialize the cell parameters for all the FE cells, will be updated in initialize function
	 */
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
	for(std::vector<CellInfo>::iterator cell_info_itr = cell_info_vector->begin();
			cell_info_itr != cell_info_vector->end();
			++cell_info_itr){
		(*cell_info_itr).dim = dim;
		(*cell_info_itr).shape_function_order = cell->active_fe_index() + 1;
		unsigned int quad_index = (*cell_info_itr).shape_function_order - 1;	//in the current quad vector based on p order
		(*cell_info_itr).quad_rule = current_quad_rule[quad_index];
		QGauss<dim> temp_quad(current_quad_rule[quad_index]);
		(*cell_info_itr).n_q_points = temp_quad.size();
		(*cell_info_itr).cell_area = 0.00001;
		(*cell_info_itr).density.resize((*cell_info_itr).n_q_points, 0.01);

		// Initialize the designField for uncoupled meshes
		if(mesh->coupling == false && mesh->adaptivityType == "adaptive_grayness"){
			(*cell_info_itr).design_points.no_points = mesh->initial_dcount_per_el;
			(*cell_info_itr).design_points.initialize_field(dim, mesh->initial_dcount_per_el, 1, volfrac);

			//here, we assume that the initial design mesh is uniform, so psuedo_design mesh is kept the same
			// note here that the no. of pseudodesign points can be a non-perfect-sqaure no., since this count is
			// used. the actual psuedo-design mesh will be created from local filtering.
			(*cell_info_itr).pseudo_design_points.no_points = mesh->initial_dcount_per_el;
			(*cell_info_itr).pseudo_design_points.initialize_field(dim, mesh->initial_dcount_per_el, 1, volfrac);
		}

		cell++;
	}
	max_design_points_per_cell = mesh->initial_dcount_per_el;	//initialized, will be updated at every cycle of refinement


	//Initialize the density cell parameters for all the density cells
	for(std::vector<CellInfo>::iterator density_cell_info_itr = density_cell_info_vector->begin();
			density_cell_info_itr != density_cell_info_vector->end();
			++density_cell_info_itr){

		//Setting the density quad rule to 1
		(*density_cell_info_itr).quad_rule = 1;
		(*density_cell_info_itr).n_q_points = 1; //	located at the centroid
		(*density_cell_info_itr).cell_area = 0.00001;	//area of the density cell
		(*density_cell_info_itr).density.resize((*density_cell_info_itr).n_q_points, volfrac);
		(*density_cell_info_itr).dxPhys.resize(mesh->design_var_per_point(), 0.0);
	}

	/**
	 * Link the cell_info_vector to the FE triangulation
	 * user_index is 1, 2, 3.......
	 */
	cell = dof_handler.begin_active();
	typename hp::DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();
	unsigned int cell_itr = 0;
	for(; cell != endc; ++cell){
		cell->set_user_index(cell_itr + 1);
		(*cell_info_vector)[cell_itr].cell_area = cell->measure(); //defining cell area
		++cell_itr;
	}


	/**
	 * Link the density_cell_info_vector to the density triangulation
	 * user index is 1, 2, 3, ...
	 */
	typename hp::DoFHandler<dim>::active_cell_iterator density_cell = design_handler.begin_active(),
			density_endc = design_handler.end();
	unsigned int density_cell_itr = 0;
	for(; density_cell != density_endc; ++density_cell){
		density_cell->set_user_index(density_cell_itr + 1);
		(*density_cell_info_vector)[density_cell_itr].cell_area = density_cell->measure();
		++density_cell_itr;
	}

	// These parameters are currently defined assuming that the final density is also represented on the analysis mesh
	density_field.max_cell_area = (*cell_info_vector)[0].cell_area;
	density_field.initial_no_cells = triangulation.n_active_cells();
	density_field.volfrac = volfrac;

}

template <int dim>
void FEM<dim>::initialize_cycle(){
	/**
	 * Link the cell_info_vector to the FE triangulation
	 * user_index is 1, 2, 3.......
	 */

	std::cout<<"Initializing the cycle "<<std::endl;
	initialize_pseudo_designField();
	std::cout<<"Psuedo design field initialized"<<std::endl;

	//Updating the hp_fe_values
	hp::FEValues<dim> hp_fe_values(fe_collection,
				quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	//Updating the cell area and quadrature index for each cell
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	unsigned int cell_itr = 0;
	GaussIntegration<dim> gauss_int;
	gauss_int.max_el_order = mesh->max_el_order;
	for(; cell != endc; ++cell){
		cell->set_user_index(cell_itr + 1);
		(*cell_info_vector)[cell_itr].cell_area = cell->measure(); //defining cell area
		unsigned int p_degree = (*cell_info_vector)[cell_itr].shape_function_order;
		(*cell_info_vector)[cell_itr].quad_rule = gauss_int.get_quadRule(
				(*cell_info_vector)[cell_itr].design_points.no_points,
				p_degree);
		QGauss<dim> temp_quad((*cell_info_vector)[cell_itr].quad_rule);
		(*cell_info_vector)[cell_itr].n_q_points = temp_quad.size();
		++cell_itr;
	}
	if (cycle > 0)	gauss_int.update_quadRuleVector(current_quad_rule, *cell_info_vector);
	std::cout<<"Quadrature rule updated "<<std::endl;

	/**
	 * Link the density_cell_info_vector to the density triangulation
	 * user index is 1, 2, 3, ...
	 */
	typename hp::DoFHandler<dim>::active_cell_iterator density_cell = design_handler.begin_active(),
			density_endc = design_handler.end();

	std::cout<<"Updating the design vector here"<<std::endl;
	density_field.update_design_vector(*cell_info_vector,
			*density_cell_info_vector,
			*design_vector,
			cycle,
			volfrac,
			*mesh,
			*projection);
	double time1 = clock();

	timer->pause();
	std::cout<<"Looking for neighbours..."<<std::endl;

	//Update thye filter radii
	projection->cycle = cycle;
	projection->update_projections(*cell_info_vector,
			dof_handler);
	std::cout<<"Projections updated "<<std::endl;
	density_field.create_neighbors(
			*cell_info_vector,
			hp_fe_values,
			dof_handler,
			design_handler,
			*projection,
			*mesh);

	double time2 = clock();
	time2 = (time2 - time1)/(double)CLOCKS_PER_SEC;
	std::cout<<"Neighbours' indices stored : time taken = "<<time2<<" seconds"<<std::endl;
	timer->resume();

}

template <int dim>
void FEM<dim>::initialize_pseudo_designField(){


	//Updating the max_design_points_per_cell parameter

	if (cycle >= 0){

		//Finding the resolution of the design mesh
		unsigned temp_max_design = 0;
		for (unsigned int i = 0; i < (*cell_info_vector).size(); ++i){
			if ((*cell_info_vector)[i].design_points.no_points > temp_max_design)
				temp_max_design = (*cell_info_vector)[i].design_points.no_points;
		}
		max_design_points_per_cell = pow(ceil(sqrt((double)temp_max_design) - 0.000000001), 2);	//this refers top pseudo-points
		//Updating the distribution and number of design points per cell
		for (unsigned int i = 0; i < (*cell_info_vector).size(); ++i){
			(*cell_info_vector)[i].pseudo_design_points.no_points = max_design_points_per_cell;
			(*cell_info_vector)[i].pseudo_design_points.initialize_field(dim, max_design_points_per_cell, 1, volfrac);

		}
	}


	//Updating the weights for projection from design mesh to pseudo-design mesh
	for (unsigned int i = 0; i < (*cell_info_vector).size(); ++i){
		(*cell_info_vector)[i].pseudo_design_points.update_pseudo_designWeights(max_design_points_per_cell,
				(*cell_info_vector)[i].design_points.pointX,
				(*cell_info_vector)[i].cell_area);
	}
}

template <int dim>
void FEM<dim>::update_pseudo_designField(){
	for (unsigned int i = 0; i < (*cell_info_vector).size(); ++i){
		(*cell_info_vector)[i].pseudo_design_points.update_pseudo_designField(
				(*cell_info_vector)[i].design_points.rho);
	}
}

template <int dim>
void FEM<dim>::update_physics(){

	//update the B and d matrices for linear elastic problem
	if(this->problem_name == "electrical_conduction"){
		electrostatic_data.update_normalized_matrices(fe_collection,
				quadrature_collection,
				dof_handler);
	}
	else{
		elastic_data.update_elastic_matrices(fe_collection,
				quadrature_collection,
				dof_handler);
	}
}

template <int dim>
void FEM<dim>::assembly(){

	hp::FEValues<dim> hp_fe_values(fe_collection,
				quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	hp::FEFaceValues<dim> hp_fe_face_values(fe_collection,
				face_quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	hp::FEValues<dim> hp_fe_analysis_density_values(fe_analysis_density_collection,
				quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	//Iterators for the FE mesh
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	unsigned int cell_itr = 0;

	//Iterator for density points on the FE mesh
	typename hp::DoFHandler<dim>::active_cell_iterator fe_den_cell = analysis_density_handler.begin_active(),
			fe_den_endc = analysis_density_handler.end();

	//Iterators for the density mesh
	typename hp::DoFHandler<dim>::active_cell_iterator density_cell = design_handler.begin_active(),
			density_endc = design_handler.end();


	add_point_source_to_rhs();
	add_point_stiffness_to_system();
	add_point_to_l_vector();

	// Defining the const vector
	const_vector.reinit(solution.size());
	if (problem_name == "electrical_conduction"){
		for (unsigned int i = 0; i < const_vector.size(); ++i){
			const_vector(i) = 1.0;
		}
	}
	else{
		for (unsigned int i = 0; i < const_vector.size(); ++i){
			const_vector(i) = l_vector(i);
		}
	}

	bool output_FE_density_field = true;
	std::string fe_fname = "fe_density/fe_density-";
	std::stringstream ss;
	ss<< cycle +1<<"_"<<itr_count+1;
	fe_fname += ss.str();
	fe_fname += ".dat";
	std::ofstream fout;
	fout.open(fe_fname, std::ios::out);
	fout.close();
	for (; cell != endc; ++cell){
		//std::cout<<cell_itr<<std::endl;
		//Getting the q_index for the cell
		unsigned int q_index;
		if (this->problem_name == "electrical_conduction")
			q_index = electrostatic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);
		else
			q_index = elastic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);

		hp_fe_values.reinit(cell, q_index);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		//std::cout<<"Dofs per cell = "<<dofs_per_cell<<std::endl;
		(*cell_info_vector)[cell_itr].dofs_per_cell = dofs_per_cell;
		const unsigned int density_per_fe_cell = fe_den_cell->get_fe().dofs_per_cell;

		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
		normalized_matrix = 0;
		Vector<double> cell_rhs(dofs_per_cell);
		Vector<double> cell_density(density_per_fe_cell);

		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		std::vector<types::global_dof_index> local_density_indices(density_per_fe_cell);

		cell_matrix = 0;
		cell_rhs = 0;
		cell_density = 0;
		unsigned int n_q_points = (*cell_info_vector)[cell_itr].n_q_points;
		//std::cout<<"NO. of q points : "<<n_q_points<<std::endl;
		QGauss<dim-1> face_quadrature_formula((*cell_info_vector)[cell_itr].quad_rule);
		unsigned int n_face_q_points = face_quadrature_formula.size();

		//Add source function to the right hand side
		std::vector<Vector<double>> rhs_values(n_q_points,
				Vector<double>(dim));
/*		add_source_to_rhs(fe_values.get_quadrature_points(),
				rhs_values);*/

		//Writing the quadrature points and density value for this cell into the file
		std::vector<Point<dim> > quad_points = fe_values.get_quadrature_points();
		assert(n_q_points == quad_points.size());
		fout.open(fe_fname, std::ios::app);
		Point<dim> em_center = cell->center();
		for (unsigned int i = 0; i < n_q_points; ++i){
		//if (cell_itr > 2) continue;
			fout<<quad_points[i](0)<<"\t"<<quad_points[i](1)<<"\t"<<
					(*cell_info_vector)[cell_itr].density[i]<<"\n";
			//std::cout<<quad_points[i](0)<<"   "<<quad_points[i](1)<<std::endl;
		}
		fout.close();

		//Calculating the cell_matrix
		double total_weight = 0.0; // For setting the density values at the nodes
		unsigned int p_index;
		if (this->problem_name == "electrical_conduction"){
			p_index = electrostatic_data.get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		}
		else{
			p_index = elastic_data.get_p_index((*cell_info_vector)[cell_itr].shape_function_order);
		}

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){

			normalized_matrix = 0;
			if (this->problem_name == "electrical_conduction"){
				normalized_matrix = electrostatic_data.elem_stiffness_array[p_index][q_index][q_point];
			}
			else
				normalized_matrix = elastic_data.elem_stiffness_array[p_index][q_index][q_point];

			//NaN condition check ----------------------------------------------------------------------------------
			if ((*cell_info_vector)[cell_itr].E_values[q_point] != (*cell_info_vector)[cell_itr].E_values[q_point])
				std::cout<<q_point<<(*cell_info_vector)[cell_itr].E_values[q_point]<<std::endl;
			//------------------------------------------------------------------------------------------------------

			cell_matrix.add((*cell_info_vector)[cell_itr].E_values[q_point],
					normalized_matrix);
			total_weight += fe_values.JxW(q_point);
			//std::cout<<fe_values.JxW(q_point)<<std::endl;
		}

		//std::cout<<total_weight<<std::endl;

		//Calculating cell_rhs
		for(unsigned int i = 0; i < dofs_per_cell; ++i){
			const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
			for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
				cell_rhs(i) += fe_values.shape_value(i, q_point) *
						rhs_values[q_point](component_i) *
						fe_values.JxW(q_point);
			}
		}

		//Adding distributed load on the edge
        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number){
          if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 62)){
              hp_fe_face_values.reinit (cell, face_number, q_index);
              const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();
              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                {
                  for (unsigned int i=0; i<dofs_per_cell; ++i){
                	  std::vector<double> distLoad = {0.0, 0.0};
          			const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
                    cell_rhs(i) += (distLoad[component_i] *
                                   fe_face_values.shape_value(i,q_point) *
                                    fe_face_values.JxW(q_point));
                  }
                	//std::cout<<q_point<<"  "<<fe_face_values.JxW(q_point)<<std::endl;

                }
            }

          if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 63)){
              hp_fe_face_values.reinit (cell, face_number, q_index);
              const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();
              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                {
                  for (unsigned int i=0; i<dofs_per_cell; ++i){
                	  std::vector<double> distLoad = {-1, 0.0};
          			const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
                    cell_rhs(i) += (distLoad[component_i] *
                                   fe_face_values.shape_value(i,q_point) *
                                    fe_face_values.JxW(q_point));
                  }
                	//std::cout<<q_point<<"  "<<fe_face_values.JxW(q_point)<<std::endl;

                }
            }

          if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 64)){
              hp_fe_face_values.reinit (cell, face_number, q_index);
              const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();
              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                {
                  for (unsigned int i=0; i<dofs_per_cell; ++i){
                	  std::vector<double> distLoad = {1, 0.0};
          			const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
                    cell_rhs(i) += (distLoad[component_i] *
                                   fe_face_values.shape_value(i,q_point) *
                                    fe_face_values.JxW(q_point));
                  }
                	//std::cout<<q_point<<"  "<<fe_face_values.JxW(q_point)<<std::endl;

                }
            }

        }
		cell->get_dof_indices(local_dof_indices);
		 hanging_node_constraints.distribute_local_to_global (cell_matrix,
		                                          cell_rhs,
		                                          local_dof_indices,
		                                          system_matrix,
		                                          system_rhs);



/*		for(unsigned int i = 0; i < dofs_per_cell; ++i){
			for(unsigned int j = 0; j < dofs_per_cell; ++j){
				system_matrix.add(local_dof_indices[i],
						local_dof_indices[j],
						cell_matrix(i, j));
			}

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}*/

		fe_den_cell->get_dof_indices(local_density_indices);
		for(unsigned int i = 0; i < density_per_fe_cell; ++i){
			for(unsigned int q_point = 0 ; q_point < n_q_points; ++q_point){
				cell_density(i) += 	(*cell_info_vector)[cell_itr].density[q_point] *
						(fe_values.JxW(q_point)/total_weight);
			}
			nodal_density(local_density_indices[i]) += cell_density(i);
			nodal_p_order(local_density_indices[i]) += (*cell_info_vector)[cell_itr].shape_function_order;
			nodal_d_count(local_density_indices[i]) += (*cell_info_vector)[cell_itr].design_points.no_points;

			cells_adjacent_per_node(local_density_indices[i]) += 1;
		}


		//Adding the density weights to the cell_info_vector
		(*cell_info_vector)[cell_itr].density_weights.clear();
		for(unsigned int qpoint = 0; qpoint < n_q_points; ++qpoint){
			(*cell_info_vector)[cell_itr].density_weights.push_back(fe_values.JxW(qpoint)/total_weight);
		}

		++fe_den_cell;
		++density_cell;
		++cell_itr;
	}

	//Add point-source function to the right hand side
	//add_point_source_to_rhs();


	double volfrac = density_field.get_vol_fraction(*cell_info_vector);
	std::cout<<"Volfrac : "<<volfrac<<std::endl;


	//Normalizing the nodal density vector
	for(unsigned int i = 0; i < nodal_density.size(); ++i){
		unsigned int denom = cells_adjacent_per_node(i);
		if (denom <= 0){
			std::cerr<<"ERROR!! Wrong no. of cells adjacent to node found"<<std::endl;
		}
		else{
			nodal_density(i) /= denom;
			nodal_p_order(i) /= denom;
			nodal_d_count(i) /= denom;

		}
	}


	//Constraining the hanging nodes
	//hanging_node_constraints.condense(system_matrix);
	//hanging_node_constraints.condense(system_rhs);

/*	MatrixTools::apply_boundary_values(boundary_values,
			system_matrix,
			solution,
			system_rhs);*/



}

template <int dim>
void FEM<dim>::add_source_to_rhs(
		const std::vector<Point<dim> > &points,
		std::vector<Vector<double> > &value_list){
	Assert (value_list.size() == points.size(),
			ExcDimensionMismatch(value_list.size(), points.size()));
	const unsigned int n_points = points.size();

	for (unsigned int p = 0; p < n_points; ++p){
		std::vector<double> x;
		x.resize(dim);
		for(unsigned int i = 0; i < dim; ++i)
			x[i] = points[p](i);
		x = mesh->source_fn(x);

		Assert(x.size() == dim,
				ExcDimensionMismatch(x.size(), dim));
		for(unsigned int i = 0; i < dim; ++i){
			value_list[p](i) = x[i];
		}
	}
}

template <int dim>
void FEM<dim>::add_point_source_to_rhs(){
		deallog.depth_console (2);

	unsigned int no_sources = (mesh->point_source_vector).size();
	for(unsigned int s = 0; s < no_sources; ++s){
		std::vector<double> load_point = mesh->point_source_vector[s].first;
		std::vector<double> load = mesh->point_source_vector[s].second;
		Assert (load_point.size() == dim,
				ExcDimensionMismatch(load_point.size(), dim));
		Assert (load.size() == dim,
				ExcDimensionMismatch(load.size(), dim));


		Point<dim> ldp, ld;
		for(unsigned int i = 0; i < dim; ++i){
			ldp(i) = load_point[i];
			ld(i) = load[i];
		}


		typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
				endc = dof_handler.end();
		//Iterating over all the cells
		for (; cell != endc; ++cell){
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
			cell->get_dof_indices(local_dof_indices);

			std::vector<Point<dim> > support_pts = cell->get_fe().get_unit_support_points();


			//Iterate over all the support points and check
			for (unsigned int i = 0; i < dofs_per_cell; ++i){
				unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
				//Below linear mapping from unit to real cell is assumed, might have to be corrected later
				Point<dim> temp_point = MappingQ<dim>(1).transform_unit_to_real_cell(cell, support_pts[i]);
				if (temp_point.distance(ldp) < 1e-12){
					system_rhs(local_dof_indices[i]) = ld(comp_i);
				}
			}
		}
	}
}

template <int dim>
void FEM<dim>::add_point_stiffness_to_system(){
	deallog.depth_console (2);
	std::cout<<"Adding point stiffnesses into the global stiffness matrix"<<std::endl;

	unsigned int no_sources = (mesh->point_stiffness_vector).size();
	if (no_sources == 0)	return;

	for(unsigned int s = 0; s < no_sources; ++s){
		std::vector<double> stiffness_point = mesh->point_stiffness_vector[s].first;
		std::vector<double> stiffness = mesh->point_stiffness_vector[s].second;
		Assert (load_point.size() == dim,
				ExcDimensionMismatch(stiffness_point.size(), dim));
		Assert (load.size() == dim,
				ExcDimensionMismatch(stiffness.size(), dim));

		Point<dim> sdp, sd;
		for(unsigned int i = 0; i < dim; ++i){
			sdp(i) = stiffness_point[i];
			sd(i) = stiffness[i];
		}

		typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
				endc = dof_handler.end();
		//Iterating over all the cells
		for (; cell != endc; ++cell){
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
			cell->get_dof_indices(local_dof_indices);

			std::vector<Point<dim> > support_pts = cell->get_fe().get_unit_support_points();


			//Iterate over all the support points and check
			for (unsigned int i = 0; i < dofs_per_cell; ++i){
				unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
				//Below linear mapping from unit to real cell is assumed, might have to be corrected later
				Point<dim> temp_point = MappingQ<dim>(1).transform_unit_to_real_cell(cell, support_pts[i]);
				if (temp_point.distance(sdp) < 1e-12){
					std::cout<<"Updated stiffness "<<std::endl;
					system_matrix.add(local_dof_indices[i], local_dof_indices[i], sd(comp_i));
				}
			}
		}

	}


}

template <int dim>
void FEM<dim>::add_point_to_l_vector(){

	std::cout<<"Adding point_l into the l_vector"<<std::endl;

	unsigned int no_sources = (mesh->point_l_vector).size();
	for(unsigned int s = 0; s < no_sources; ++s){
		std::vector<double> load_point = mesh->point_l_vector[s].first;
		std::vector<double> load = mesh->point_l_vector[s].second;
		Assert (load_point.size() == dim,
				ExcDimensionMismatch(load_point.size(), dim));
		Assert (load.size() == dim,
				ExcDimensionMismatch(load.size(), dim));


		Point<dim> ldp, ld;
		for(unsigned int i = 0; i < dim; ++i){
			ldp(i) = load_point[i];
			ld(i) = load[i];
		}


		typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
				endc = dof_handler.end();
		//Iterating over all the cells
		for (; cell != endc; ++cell){
			const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
			std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
			cell->get_dof_indices(local_dof_indices);

			std::vector<Point<dim> > support_pts = cell->get_fe().get_unit_support_points();


			//Iterate over all the support points and check
			for (unsigned int i = 0; i < dofs_per_cell; ++i){
				unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
				//Below linear mapping from unit to real cell is assumed, might have to be corrected later
				Point<dim> temp_point = MappingQ<dim>(1).transform_unit_to_real_cell(cell, support_pts[i]);
				if (temp_point.distance(ldp) < 1e-12){
					l_vector(local_dof_indices[i]) = ld(comp_i);
				}
			}
		}

	}
}

template <int dim>
void FEM<dim>::add_boundary_constraints(){


	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();

	//Iterating over all the cells
	for (; cell != endc; ++cell){
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);

		std::vector<Point<dim> > support_pts = cell->get_fe().get_unit_support_points();

		//Iterate over all the faces to check for boundary_id

		bool indic52 = false;

		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f){

			if (cell->face(f)->boundary_id() == 52 || cell->face(f)->boundary_id() == 53 || cell->face(f)->boundary_id() == 54){
				indic52 = true;
				break;
			}
		}

		if (indic52 == true){

			//Iterate over all the support points and check
			for (unsigned int i = 0; i < dofs_per_cell; ++i){
				if (dim == 2){
					//Below linear mapping from unit to real cell is assumed, might have to be corrected later
					Point<dim> temp_points = MappingQ<dim>(1).transform_unit_to_real_cell(cell, support_pts[i]);
					unsigned int boundary_indic = mesh->boundary_indicator({temp_points(0),
							temp_points(1)});
					if (boundary_indic == 52){
						unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
						//std::cout<<temp_points[0]<<"   "<<temp_points[1]<<"    "<<boundary_indic<<std::endl;

						if (comp_i == 1){
							hanging_node_constraints.add_line(local_dof_indices[i]);
							//std::cout<<hanging_node_constraints.n_constraints()<<std::endl;
						}

					}
					else if (boundary_indic == 53){
						unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
						//std::cout<<temp_points[0]<<"   "<<temp_points[1]<<"    "<<boundary_indic<<std::endl;

						if (comp_i == 0){
							hanging_node_constraints.add_line(local_dof_indices[i]);
							//std::cout<<hanging_node_constraints.n_constraints()<<std::endl;
						}
					}
					else if (boundary_indic == 54){
							hanging_node_constraints.add_line(local_dof_indices[i]);
					}
				}
			}
		}
	}
}

//Function for setting boundary indicators for the domain
template <int dim>
void FEM<dim>::boundary_info(){
	for (typename Triangulation<dim>::active_cell_iterator
			cell = triangulation.begin();
			cell != triangulation.end();
			++cell ){
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f){
			Point<dim> point1 = cell->face(f)->center();
			std::vector<double> x(dim);
			for(unsigned int  i = 0; i < dim; ++i){
				x[i] = point1(i);
			}
			unsigned int indic = mesh->boundary_indicator(x);
			if (indic == 9999)	continue;
			cell->face(f)->set_boundary_id(indic);

		}
	}
}

template <int dim>
void FEM<dim>::clean_trash(){
	system_matrix.clear();
	system_rhs = 0;
	nodal_density = 0;
	nodal_p_order = 0;
	nodal_d_count = 0;
	cells_adjacent_per_node = 0;
	solution = 0;

	l_vector = 0;
	lambda_solution = 0;

	if (mesh->coupling == true || mesh->adaptivityType == "movingdesignpoints"){
		unsigned int no_des_per_point = mesh->design_var_per_point();
		//cleaning contents of the storage vectors
		for(unsigned int i = 0 ; i < density_cell_info_vector->size(); ++i){
			(*density_cell_info_vector)[i].dxPhys.clear();
			(*density_cell_info_vector)[i].dxPhys.resize(no_des_per_point, 0.0);
		}

	}
	else{
		for (unsigned int i = 0; i < cell_info_vector->size(); ++i){
			(*cell_info_vector)[i].design_points.dxPhys_drho.clear();
			(*cell_info_vector)[i].design_points.dxPhys_drho.resize((*cell_info_vector)[i].design_points.no_points);

			for(unsigned int j = 0; j < (*cell_info_vector)[i].design_points.no_points; ++j){
				(*cell_info_vector)[i].design_points.dxPhys_drho[j] = 0.0;
			}

			(*cell_info_vector)[i].refine_coarsen_flag = 0;
		}
	}

	boundary_values.clear();
	std::cout<<"Trash cleaned "<<std::endl;

}












