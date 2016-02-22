/*
 *
 *  Created on: Aug 16, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/TopologyOptimization/constraints/general.h>

using namespace atop;

template <int dim>
double GeneralConstraints<dim>::volumeConstraint(
		std::vector<double> &volume_grad_vector,
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		DensityField<dim> &density_field,
		unsigned int design_var_per_point
		){

	unsigned int no_design = density_cell_info_vector.size() * design_var_per_point;
	volume_grad_vector.resize(no_design, 0.0);

	//Calculating the total material volume fraction
	double volume = density_field.get_vol_fraction(cell_info_vector);
	//Calculating the first order gradient of volume constraint
	unsigned int k = 0;
	for(unsigned int i = 0; i < density_cell_info_vector.size(); ++i){
		for (unsigned int j = 0; j < design_var_per_point; j++){
			volume_grad_vector[k] = 0.25 * density_cell_info_vector[i].dxPhys[j] / density_field.initial_no_cells;	k++;
		}
		//std::cout<<volume_grad_vector[i]<<std::endl;
	}

	return (volume - density_field.volfrac);
}
