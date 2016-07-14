/*
 *
 *  Created on: Jun 22, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef CREATE_DESIGN_H_
#define CREATE_DESIGN_H_

#include <atop/fem/fem.h>
#include <deal.II/lac/vector.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/DensityValues.h>



using namespace dealii;

namespace atop{
	template <int dim>
	class CreateDesign{
	public:
		CreateDesign(){};

		void assemble_design(
				FEM<dim> &);

	private:
		FEM<dim> *fem;
		DensityField<dim> design_values;
		Vector<double> cells_adjacent_per_node;
		std::vector<CellInfo> design_info_vector;
	};

	template class CreateDesign<2>;
}


#endif /* CREATE_DESIGN_H_ */
