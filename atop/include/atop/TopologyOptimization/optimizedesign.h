/*
 *
 *  Created on: Apr 2, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef OPTIMIZEDESIGN_H_
#define OPTIMIZEDESIGN_H_

#include <atop/fem/define_mesh.h>
#include <atop/fem/fem.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/physics/mechanics/elastic.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <string.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/derivatives/sensitivity.h>
#include <atop/TopologyOptimization/constraints/general.h>
#include <iostream>
#include <vector>

using namespace dealii;

namespace atop{
template <int dim>
	class Optimizedesign{
	public:

	Triangulation<dim> triangulation;
	Triangulation<dim> fe_density_triangulation;
	Triangulation<dim> density_triangulation;
	DoFHandler<dim> dof_handler, density_handler, density_dof_handler;
	std::vector<CellInfo> cell_info_vector;
	std::vector<CellInfo> density_cell_info_vector;
	DefineMesh<dim>* mesh;
	FEM<dim> *obj_fem;
	unsigned int cycle, no_cycles, itr_count;

	double volfrac;	//volume fraction

	//Objects for defining the physics of the problem
	LinearElastic<dim> *linear_elastic;

	//name of optimization problem
	std::string problem_name;

	//Object for providing the projection operator details
	Projection *projection;

	//Object for penalization of densities
	Penalize *penal;

	//Object of general constraints class
	GeneralConstraints<dim> vol_constraint;

	//Optimization related parameters
	std::vector<double> design_vector;
	double objective;
	std::vector<double> grad_vector;

		Optimizedesign();

		Optimizedesign(
				atop::DefineMesh<dim>&,
				atop::Penalize&,
				atop::Projection&,
				const std::string&,
				unsigned int cycles = 5);

		void problemType(LinearElastic<dim> &obj_linear_elastic);

		void optimize();
		void run_system();
		void update_design_vector(
				std::vector<double> &,
				const std::vector<double> &);

	private:
		void flush();


	};

template class Optimizedesign<2>;
}



#endif /* OPTIMIZEDESIGN_H_ */
