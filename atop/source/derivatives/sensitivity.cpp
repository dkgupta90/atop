/*
 *
 *  Created on: Aug 13, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/derivatives/sensitivity.h>
#include <deal.II/dofs/dof_handler.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/fe/fe_system.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/fem/fem.h>
#include <atop/derivatives/compliance.h>

using namespace dealii;
using namespace atop;

template <int dim>
SensitivityAnalysis<dim>::SensitivityAnalysis(
		std::string &str){
	this->problem_name = str;
}

template <int dim>
void SensitivityAnalysis<dim>::set_input(
		DoFHandler<dim> &obj_dof_handler,
		std::vector<CellInfo> &obj_cell_info_vector,
		std::vector<CellInfo> &obj_density_cell_info_vector,
		FEM<dim> &obj_fem
		){

	//Setting the dof_handler object
	this->dof_handler = &obj_dof_handler;
	this->cell_info_vector = &obj_cell_info_vector;
	this->density_cell_info_vector = &obj_density_cell_info_vector;
	this->fem = &obj_fem;
}

template <int dim>
void SensitivityAnalysis<dim>::run(
		double &obj,
		std::vector<double> &obj_grad){
	if (problem_name == "minimum_compliance"){
		Compliance<dim> obj_compliance;
		obj_compliance.set_input(
				*dof_handler,
				*cell_info_vector,
				*density_cell_info_vector,
				*fem);
		obj_compliance.compute(obj, obj_grad);
	}
}



