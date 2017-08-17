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
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <atop/fem/output.h>


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


		Vector<double> cells_adjacent_per_node;

		QRIndicator(
				FEM<dim> &obj_fem,
				std::vector<double> &qr_accuracy,
				double tol_accuracy,
				std::vector<unsigned int> &qr_p_value,
				std::vector<CellInfo> &cell_info_vector);

		void estimate();
		double get_Jvalue(hp::DoFHandler<2>::active_cell_iterator cell,
				Vector<double> &u_solution,
				Vector<double> &f_solution,
				unsigned int new_p);

	};
	template class  QRIndicator<2>;

}

#endif /* INCLUDE_ATOP_FEM_ERRORINDICATOR_QRINDICATOR_H_ */
