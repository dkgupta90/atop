/*
 * penalization.h
 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef PENALIZATION_H_
#define PENALIZATION_H_

#include<vector>
#include <string>
#include <atop/TopologyOptimization/cell_prop.h>

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

namespace atop{
	class Penalize{
	public:
		//Default constructor
		Penalize();

		std::string scheme;
		double factmin;	//Factor for determining the lower bound for material property
						// for e.g. Emin = factmin * E;
		double penal_power;
		//Constructor for initializing the penalization scheme
		Penalize(std::string);

		void update_param(
				double,
				std::vector<CellInfo>&);

		void update_param(
				double,
				CellInfo &);

		double penalized_factor(double xPhys);

	};
}


#endif /* PENALIZATION_H_ */
