/*
 *
 *  Created on: Aug 13, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include <string>
#include <iostream>
#include <deal.II/hp/dof_handler.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/fe/fe_system.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/fem/fem.h>

using namespace dealii;
namespace atop{
	template <int dim>
	class SensitivityAnalysis{
	public:
		std::string problem_name;

		//For objective value


		SensitivityAnalysis(
				std::string &);

		void set_input(
				hp::DoFHandler<dim> &dof_handler,
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector,
				FEM<dim> &fem
				);

		void run(double &obj,
				std::vector<double> &obj_grad
				);

	private:
		hp::DoFHandler<dim> *dof_handler;
		std::vector<CellInfo> *cell_info_vector;
		std::vector<CellInfo> *density_cell_info_vector;
		FEM<dim> *fem;
	};

	template class SensitivityAnalysis<2>;
}



#endif /* DERIVATIVES_H_ */
