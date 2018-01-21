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

	if (dim == 2){
		//calculating polynomial order of design field

		/*
		 * for dim = 2, there are 2 variables x and y
		 * for 2 variables, the no. of terms in the full polynomial are
		 * equal to (p+1)(p+2)/2, where p is the polynomial order
		 */
		unsigned int p_design = ceil((double)(-3.0 + sqrt(1.0 + 8 * no_design) - 1e-10)/2.0);
		unsigned int p_total = p_design + 2 * p_degree;
		unsigned int qrule = ceil((double)((p_total + 1)/2.0))+0;
		return qrule;
	}
	else if (dim == 3){
		//calculating polynomial order of design field

		/*
		 * for dim = 3, there are 3 vairables x, y and z
		 * for 3 variables, the no. of terms in the full polynomial are
		 * equal to (p+1)(p+2)(p+3)/6, where p is the polynomial order
		 */
		double t1 = 1.732050807 * sqrt((double)(243 * no_design * no_design - 1));
		double t2 = 27 * no_design;
		double t = pow(t1 + t2, 0.3333333333);
		double p_design = ceil((double)(0.48075498567 * t + (0.69336127435 / t) - 2 - 0.001));
		unsigned int p_total = p_design + 2 * p_degree;
		unsigned int qrule = ceil((double)((p_total + 1)/2.0))+0;
		//std::cout<<"Pdesign, Ptotal, qrule : "<<p_design<<" "<<p_total<<" "<<qrule<<std::endl;
		return qrule;
	}

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

		/*
		 * the condition below checks whether for certain quad_rule being used, a quad_rule of 2 orders higher exist.
		 * This is probably used during q-refinment, where the q value gets increased and problem is solved.
		 */
		if (cell_info_vector[i].shape_function_order < max_el_order){
			if (cell_info_vector[i].quad_rule + 2 > quadRuleVector[(cell_info_vector[i].shape_function_order + 1) - 1]){
				quadRuleVector[(cell_info_vector[i].shape_function_order + 1) - 1] = cell_info_vector[i].quad_rule + 2;
			}
		}

	}
}

#endif /* INTEGRATION_CPP_ */
