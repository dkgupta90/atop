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
#include <atop/fem/create_design.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/physics/electrical/electrostatic.h>
#include <deal.II/grid/tria.h>
#include <string.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/derivatives/sensitivity.h>
#include <atop/TopologyOptimization/constraints/general.h>
#include <iostream>
#include <vector>
#include <deal.II/hp/dof_handler.h>

using namespace dealii;

namespace atop{
template <int dim>
	class Optimizedesign{
	public:

	// dof_handler connects to triangulation
	// density_handler connects to fe_density_triangulation
	// density_dof_handler connects to density_triangulation

	std::vector<CellInfo> cell_info_vector;	//stores the information related to each analysis cell
	std::vector<CellInfo> density_cell_info_vector;	//stores information related to each cell on design mesh
	DefineMesh<dim>* mesh;
	FEM<dim> *obj_fem;
	unsigned int cycle, no_cycles;

	double volfrac;	//volume fraction

	//Objects for defining the physics of the problem
	LinearElastic<dim> *linear_elastic;
	LinearElectrostatic<dim> *lin_elecstat;

	//name of optimization problem
	std::string problem_name;
	bool is_problem_self_adjoint;

	//Object for providing the projection operator details
	Projection *projection;

	//Object for penalization of densities
	Penalize *penal;

	//Timers to time the optimization routine
	double start_time, end_time;

	//Object of general constraints class
	GeneralConstraints<dim> vol_constraint;

	bool temp1;
	std::string tempfname;	// this are temporarhy for qr-test
	unsigned int final_dcount_per_el;
	//Object of Timer class for performance evlation
	Timer timer;

	//Optimization related parameters
	std::vector<double> design_vector;
	double objective;
	std::vector<double> grad_vector;

	std::vector<double> lb, ub;
	unsigned int design_count, no_constraints;

		Optimizedesign();

		Optimizedesign(
				atop::DefineMesh<dim>&,
				atop::Penalize&,
				atop::Projection&,
				const std::string&,
				unsigned int cycles = 5);

		void problemType(LinearElastic<dim> &obj_linear_elastic);
		void problemType(LinearElectrostatic<dim> &obj_linear_electrostatic);

		void optimize();
		void run_system();
		void update_design_vector(
				std::vector<double> &,
				const std::vector<double> &);

	private:
		void flush();


	};

template class Optimizedesign<2>;
template class Optimizedesign<3>;
}



#endif /* OPTIMIZEDESIGN_H_ */
