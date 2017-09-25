/*
 *
 *  Created on: Sep 15, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef ELECTRICAL_COMPLIANCE_H_
#define ELECTRICAL_COMPLIANCE_H_


#include <deal.II/hp/dof_handler.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/fe/fe_system.h>
#include <atop/physics/electrical/electrostatic.h>
#include <atop/fem/fem.h>

using namespace dealii;
namespace atop{
	template <int dim>
	class ElectricalCompliance{
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


	template class ElectricalCompliance<2>;
}


#endif /* INCLUDE_ATOP_DERIVATIVES_ELECTRICAL_COMPLIANCE_H_ */