/*
 *
 *  Created on: Aug 7, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INCLUDE_ATOP_FEM_ERRORINDICATOR_QRINDICATOR_H_
#define INCLUDE_ATOP_FEM_ERRORINDICATOR_QRINDICATOR_H_

/** This class compute the error in solution for each element and udpates the p-order accordingly
 * based on needed accuracy
 */


#include <vector>
#include<iostream>
#include <atop/fem/fem.h>
#include<atop/TopologyOptimization/cell_prop.h>

using namespace dealii;

namespace atop{
	template <int dim>
	class QRIndicator{
	public:

		std::vector<double> *accuracy_vector;
		std::vector<unsigned int> *proposed_p_values;
		std::vector<CellInfo> *cell_info_vector;
		double tol_accuracy;
		FEM<dim> *fem;

		QRIndicator(
				FEM<dim> &obj_fem,
				std::vector<double> &qr_accuracy,
				double tol_accuracy,
				std::vector<unsigned int> &qr_p_value,
				std::vector<CellInfo> &cell_info_vector);

		void estimate();

	};
}

#endif /* INCLUDE_ATOP_FEM_ERRORINDICATOR_QRINDICATOR_H_ */
