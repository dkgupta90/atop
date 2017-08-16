/*
 *
 *  Created on: Aug 12, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INCLUDE_ATOP_PHYSICS_MECHANICS_ELECTROSTATIC_H_
#define INCLUDE_ATOP_PHYSICS_MECHANICS_ELECTROSTATIC_H_

#include <string>
#include <vector>
#include <iostream>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>

using namespace dealii;

namespace atop{
	template <int dim>
	class LinearElectrostatic{
	public:
		double E0, Emin;
	};
}



#endif /* INCLUDE_ATOP_PHYSICS_MECHANICS_ELECTROSTATIC_H_ */
