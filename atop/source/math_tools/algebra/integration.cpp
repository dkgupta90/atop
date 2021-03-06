/*
 *
 *  Created on: Jul 21, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INTEGRATION_CPP_
#define INTEGRATION_CPP_

#include <atop/math_tools/algebra/integration.h>
#include <math.h>
#include <iostream>
#include <atop/TopologyOptimization/cell_prop.h>
#include <vector>

using namespace atop;

template <int dim>
unsigned int GaussIntegration<dim>::get_quadRule(
		unsigned int no_design,
		unsigned int p_degree){

	//calculating polynomial order of design field
	unsigned int p_design = ceil((double)(-3.0 + sqrt(1.0 + 8 * no_design))/2.0);
	unsigned int p_total = p_design + 2 * p_degree;
	unsigned int qrule = ceil((double)((p_total + 1)/2.0))+0;
	qrule = 10;
	return qrule;
}

template <int dim>
void GaussIntegration<dim>::initialize_quadRuleVector(
		std::vector<unsigned int> &quadRuleVector,
		unsigned int max_el_order,
		unsigned int initial_no_design){

	//This vector saves the quadRule required for each polynomial order with the initial design field chosen
	quadRuleVector.clear();
	quadRuleVector.resize(max_el_order);
	for (unsigned int i = 0; i < quadRuleVector.size(); ++i){
		unsigned int qrule = get_quadRule(initial_no_design, i+1);
		quadRuleVector[i] = qrule;
	}
}

/*
 * This function iterates over all the cells and compares the quad rule of each cell to the quad rule stored corresponding
 * to that polynomial order. If the former is larger, the later is accordingly updated.
 * Larger quad rules are possible because of the design resolution. Thus, each cell is checked if for certain polynomial order,
 * a higher quad rule might be needed,
 */
template <int dim>
void GaussIntegration<dim>::update_quadRuleVector(
		std::vector<unsigned int> &quadRuleVector,
		std::vector<CellInfo> &cell_info_vector){
	for (unsigned int i = 0; i < cell_info_vector.size(); ++i){
		if (cell_info_vector[i].quad_rule > quadRuleVector[cell_info_vector[i].shape_function_order - 1]){
			quadRuleVector[cell_info_vector[i].shape_function_order - 1] = cell_info_vector[i].quad_rule;
		}
	}
}

#endif /* INTEGRATION_CPP_ */
