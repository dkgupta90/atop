/*
 *
 *  Created on: Apr 20, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef DESIGNFIELD_H_
#define DESIGNFIELD_H_

#include<deal.II/base/point.h>

using namespace dealii;

namespace atop{
	class DesignField{
	public:

		unsigned int dim;	//dimension of the problem
		unsigned int design_rule;	//rule for distributing the design points
		DesignField();

		unsigned int no_points;	//no. of design points per analysis element

		std::vector<double> rho;	//store design density values
		std::vector<double> dxPhys_drho;	//to be used for the volume constraint
		std::vector<std::vector<double> > pointX;	//store location of each of these points

		void initialize_field(unsigned int,
				unsigned int,
				unsigned int,
				double);
		void update_field_manual(unsigned int);

	};
}



#endif /* DESIGNFIELD_H_ */