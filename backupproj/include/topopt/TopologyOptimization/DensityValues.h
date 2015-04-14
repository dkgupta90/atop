/*
 * DensityValues.h
 *
 *  Created on: Jul 14, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef DENSITYVALUES_H_
#define DENSITYVALUES_H_

#include <deal.II/base/point.h>
#include <vector>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <topopt/TopologyOptimization/design_analysis.h>
#include <topopt/TopologyOptimization/neighbors.h>
#include <topopt/TopologyOptimization/cell_prop.h>
#include <topopt/physics/elasticity.h>
#include <deal.II/dofs/dof_handler.h>

using namespace dealii;
namespace topopt{
	class DensityValues{
	public:
		double xmin, xmax, ymin, ymax;
		unsigned int x_count, y_count;
		double x_spacing, y_spacing;
		std::vector<double> rho_mesh;
		unsigned int itr_count;
		double max_projection_radius, min_projection_radius;
		double max_cell_area;
		double gamma;
		Point<2> design_point;
		unsigned int cycle;
		std::vector<std::vector<double> > density2d;
		void update_density_mesh(
				std::vector<CellProperties> &cellprop,
				std::vector<double> &density_mesh);
		void get_vol_fraction(
				std::vector<CellProperties> &cellprop,
				std::vector<double> &density_mesh,
				double &density_sum,
				double &max_cell_area,
				StoreElasticData &elastic_data);
		void create_neighbours(
				std::vector<CellProperties> &cellprop,
				FESystem<2> &fe,
				DoFHandler<2> &dof_handler,
				double rmin
				);

		void neighbor_search(DoFHandler<2>::active_cell_iterator cell1,
				DoFHandler<2>::active_cell_iterator cell,
				std::vector<DoFHandler<2>::active_cell_iterator> &neighbor_iterators,
				double rmin
				);

		void calculate_weights(std::vector<CellProperties>  &cellprop,
				unsigned int cell_itr1,
				double rmin);

		void filter(std::vector<CellProperties> &cellprop);

		double get_dxPhys_dx(CellProperties &cellp,
				unsigned int q_point,
				unsigned int cell_itr2,
				unsigned int q_point2);
	};
}



#endif /* DENSITYVALUES_H_ */
