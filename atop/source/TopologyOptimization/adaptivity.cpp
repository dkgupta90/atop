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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>

using namespace atop;
using namespace dealii;

template <int dim>
void Adaptivity<dim>::update(
		FEM<dim> &obj_fem){

	//Assigning the FEM object
	this->fem = &obj_fem;
	this->cell_info_vector = fem->cell_info_vector;
	this->density_cell_info_vector = fem->density_cell_info_vector;

	//Run the mesh adaptivity approach
	mesh_refine_indicator(fem->mesh->adaptivityType);


	//Run projection adaptivity


	//Run penalization continuation


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

	}
}

template <int dim>
void Adaptivity<dim>::coupled_refine_adaptive_grayness(){

	//Getting the cycle number for refinement
	unsigned int cycle = fem->cycle;

	double rhomin = 0.0;
	double rhomax = 1.0;
	double rhomid = (rhomax + rhomin)/2;
	double alpha = 0.2;
	double beta = 1.2;

	double refine_lbound = rhomin + ((1 - alpha) * rhomid * exp(-beta * (double)(cycle+1)));
	double refine_ubound = rhomax - ((1 - alpha) * rhomid * exp(-beta * (double)(cycle+1)));
	double coarsen_lbound = rhomin + (alpha * rhomid * exp(-beta * (double)(cycle+1)));
	double coarsen_ubound = rhomax - (alpha * rhomid * exp(-beta * (double)(cycle+1)));

	typename Triangulation<dim>::active_cell_iterator cell = fem->triangulation->begin_active(),
			endc = fem->triangulation->end();

	typename Triangulation<dim>::active_cell_iterator density_cell = fem->density_triangulation->begin_active(),
			density_endc = fem->density_triangulation->end();

	unsigned int cell_itr = 0;

	for(; cell != endc; ++cell, ++density_cell){

		unsigned int n_qpoints = (*cell_info_vector)[cell_itr].density.size();
		double density = 0.0;

		for(unsigned int qpoint = 0; qpoint < n_qpoints; ++qpoint){
			density += (*cell_info_vector)[cell_itr].density[qpoint] *
								(*cell_info_vector)[cell_itr].density_weights[qpoint];

		}
		if(density > coarsen_ubound || density < coarsen_lbound){
			cell->set_coarsen_flag();
			density_cell->set_coarsen_flag();
		}
		else if (density > refine_lbound && density < refine_ubound) {
			cell->set_refine_flag();
			density_cell->set_refine_flag();
		}
		++cell_itr;
	}
}
