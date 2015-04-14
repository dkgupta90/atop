/*
 * penalization.h
 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef PENALIZATION_H_
#define PENALIZATION_H_

#include<vector>

namespace topopt{
	class Penalization{
	public:
		Penalization();
		void set_param(double E0, double Emin,
				std::vector<double> &E_values,
				std::vector<double> &dE_values,
				std::vector<double> &density_values,
				unsigned int penalization_model,
				double penal_power);
	};
}



#endif /* PENALIZATION_H_ */
