/*
 *
 *  Created on: Aug 16, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef GENERAL_H_
#define GENERAL_H_

#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/DensityValues.h>

namespace atop{
template <int dim>
	class GeneralConstraints{
	public:
		//material volume constraint
		double volumeConstraint(
				std::vector<double> &volume_grad_vector,
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &density_cell_info_vector,
				DensityField<dim> &density_field,
				unsigned int design_var_per_point);
	};

template class GeneralConstraints<2>;
}



#endif /* GENERAL_H_ */
