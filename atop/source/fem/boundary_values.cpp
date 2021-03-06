/*
 *
 *  Created on: Aug 12, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/boundary_values.h>
#include <deal.II/numerics/vector_tools.h>

using namespace atop;

template <int dim>
inline
void BoundaryValues<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values) const{
	Assert (values.size() == dim,
			ExcDimensionMismatch(values.size(), dim));
	Assert(dim >= 2,
			ExcNotImplemented());
	unsigned int model_problem = 3;
	double xmin, xmax, ymin, ymax;

	if (model_problem == 4){
		xmin = 0, ymin = 0, xmax = 2, ymax = 1;
		if (std::fabs(p(1) - (ymax)) < 1e-12){
			values(0) = 0;
			values(1) = 0;
		}
	}
	else if (model_problem == 5){
		//cantilever problem
		xmin = 0, ymin = 0, xmax = 1, ymax = 1;
		if (std::fabs(p(0) - (ymin)) < 1e-12){
			values(0) = 1;
			values(1) = 0;
		}
	}
	else if (model_problem == 3){
		//cantilever problem
		xmin = 0, ymin = 0, xmax = 2, ymax = 1;
		if (std::fabs(p(0) - (xmin)) < 1e-12){
			values(0) = 0;
			values(1) = 0;
		}
	}
	else if (model_problem == 2) {
		//Mitchel stress problem
		xmin = 0, ymin = 0, xmax = 1.5, ymax = 1;
		double rball = 0.02;
		Point<dim> center;
		center(0) = xmin + rball;
		center(1) = (ymax+ymin)/2;
		if (center.distance(p) <= rball){
			values(0) = 0;
			values(1) = 0;
		}
	}
	//below one has to be checked, looks to be wrong
	else if (model_problem == 1){
		//MBB problem
		xmin = 0, ymin = 0, xmax = 1, ymax = 1;
		if (sqrt(pow(p[0] - 0.5, 2) + pow(p[1] - 0.5, 2)) < 0.027){
			values(0) = 0;
			values(1) = 0;
		}
	}

}

template <int dim>
void BoundaryValues<dim>::vector_value_list(const std::vector<Point<dim>> &points,
		std::vector<Vector<double>> &value_list) const {
	Assert (value_list.size() == points.size(),
			ExcDimensionMismatch(value_list.size(), points.size()));
	const unsigned int n_points = points.size();

	for (unsigned int p = 0; p < n_points; ++p){
		BoundaryValues<dim>::vector_value(points[p],
				value_list[p]);
	}
}


