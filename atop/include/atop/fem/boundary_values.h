/*
 *
 *  Created on: Aug 12, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef BOUNDARY_VALUES_H_
#define BOUNDARY_VALUES_H_

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>
using namespace dealii;

namespace atop{
	template <int dim>
	class BoundaryValues : public Function<dim>{
	public:
		BoundaryValues(): Function<dim>() {}
		virtual double value(const Point<dim> &p,
				 const unsigned int comp) const;
		virtual void value_list (const std::vector<Point<dim>> & points,
				std::vector<double> &value_list) const;
		virtual void vector_value(const Point<dim> &p,
				Vector<double> &values) const;
		virtual void vector_value_list (const std::vector<Point<dim>> & points,
				std::vector<Vector<double>> &value_list) const;

	};

	//template class BoundaryValues<2>;
	template class BoundaryValues<3>;
}

#endif /* BOUNDARY_VALUES_H_ */
