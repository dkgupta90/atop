/*
 *
 *  Created on: Aug 22, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef ADAPTIVITY_H_
#define ADAPTIVITY_H_

#include <atop/fem/fem.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <string>
#include <vector>

namespace atop{
	template <int dim>
	class Adaptivity{
	public:
		FEM<dim> *fem;
		std::vector<CellInfo> *cell_info_vector;
		std::vector<CellInfo> *density_cell_info_vector;
		void update(
				FEM<dim> &obj_fem
				);

		void mesh_refine_indicator(
				std::string& mesh_update_str);
		void coupled_refine_adaptive_grayness();
	};

	template class Adaptivity<2>;
}



#endif /* ADAPTIVITY_H_ */
