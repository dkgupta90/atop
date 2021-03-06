/*
 *
 *  Created on: Aug 19, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef OPTIMALITY_CRITERIA_H_
#define OPTIMALITY_CRITERIA_H_

#include <iostream>
#include <vector>

namespace atop{
	class OC{
	public:
		unsigned int no_design;
		double min_obj_change;
		std::vector<double> *lb, *ub, *design_vector;
		void *obj_data;

		OC(unsigned int);
		void set_lower_bounds(
				std::vector<double> &);
		void set_upper_bounds(
				std::vector<double> &);
		double (*obj_fn)(const std::vector<double> &x,
				std::vector<double> &grad,
				void *my_func_data);
		double (*constraint_fn)(
				const std::vector<double> &x,
				std::vector<double> &grad,
				void *my_func_data);
		void optimize(
				std::vector<double> &);


	};
}


#endif /* OPTIMALITY_CRITERIA_H_ */
