/*
 *
 *  Created on: Sep 18, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef DERIVATIVES_AVERAGE_STATE_H_
#define DERIVATIVES_AVERAGE_STATE_H_

#include <deal.II/hp/dof_handler.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/fe/fe_system.h>
#include <atop/physics/electrical/electrostatic.h>
#include <atop/fem/fem.h>

using namespace dealii;
namespace atop{
	template <int dim>
	class VoltageAverage{
	public:
		void set_input(
				hp::DoFHandler<dim> &dof_handler,
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector,
				FEM<dim> &
				);
		void compute(
				double &obj,
				std::vector<double> &obj_grad);

	private:
		hp::DoFHandler<dim> *dof_handler;
		std::vector<CellInfo> *cell_info_vector;
		std::vector<CellInfo> *density_cell_info_vector;
		hp::FEValues<dim> *hp_fe_values;
		ElectrostaticData<dim> *electrostatic_data;
		FEM<dim> *fem;
		DensityField<dim> *density_field;




	};


	template class VoltageAverage<2>;
}

#endif /* INCLUDE_ATOP_DERIVATIVES_AVERAGE_state_H_ */
