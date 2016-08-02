/*
 *
 *  Created on: Jul 21, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#include <math.h>
#include <vector>

namespace atop{
	template <int dim>
	class GaussIntegration{
	public:
		unsigned int get_quadRule(unsigned int,
				unsigned int);

		//This function is meant to initialize the current and running quad rule vectors
		void initialize_quadRuleVector(
				std::vector<unsigned int>&,
				unsigned int,
				unsigned int);
	};

	template class GaussIntegration<2>;
}


