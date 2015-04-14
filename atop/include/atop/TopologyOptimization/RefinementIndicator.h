/*
 * RefinementIndicator.h
 *
 *  Created on: Dec 7, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef REFINEMENTINDICATOR_H_
#define REFINEMENTINDICATOR_H_
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/physics/elasticity.h>
#include<vector>

using namespace dealii;
namespace topopt{
class DensityIndicator{
public:
	unsigned int cycle, no_cycles;
	void get_density_indicator(Triangulation<2> &triangulation,
			Triangulation<2> &density_triangulation,
			std::vector<CellProperties> &cellprop,
			StoreElasticData &elastic_data);
};

class FEIndicator{
public:
	void get_fe_indicator(Triangulation<2> &triangulation);
};
}


#endif /* REFINEMENTINDICATOR_H_ */
