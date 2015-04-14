/*
 * cell_prop.cpp
 *
 *  Created on: Dec 28, 2014
 *      Author: dkgupta
 */

#include <atop/TopologyOptimization/cell_prop.h>

using namespace topopt;

void CellProperties::initialize_density(double value){
	for(unsigned int i = 0; i < n_q_points; ++i){
		material_density.push_back(value);
		xPhys.push_back(value);

	}
}




