/*
 * RefinementIndicator.cpp
 *
 *  Created on: Dec 7, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include<atop/TopologyOptimization/RefinementIndicator.h>
#include<iostream>
#include <math.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <atop/physics/elasticity.h>



using namespace topopt;
using namespace dealii;

void DensityIndicator::get_density_indicator(
		Triangulation<2> &triangulation,
		Triangulation<2> &density_triangulation,
		std::vector<CellProperties> &cellprop,
		StoreElasticData &elastic_data){
	typename Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
			endc = triangulation.end();
	typename Triangulation<2>::active_cell_iterator density_cell = density_triangulation.begin_active(),
			density_endc = density_triangulation.end();
	double rhomin = 0.0;
	double rhomax = 1.0;
	double rhomid = (rhomax + rhomin)/2;
	double alphar = (rhomax + rhomin)/2;
	double betar = 6;
	double alphac = 0.1;
	double betac = 4;

	double refine_lbound = rhomin + (alphar * exp(-betar * (double)(cycle+1)/(double)no_cycles));
	double refine_ubound = rhomax - (alphar * exp(-betar * (double)(cycle+1)/(double)no_cycles));
	double coarsen_lbound = rhomin + (alphac * exp(-betac * (double)(cycle+1)/(double)no_cycles));
	double coarsen_ubound = rhomax - (alphac * exp(-betac * (double)(cycle+1)/(double)no_cycles));
	//std::cout<<refine_lbound<<" "<<refine_ubound<<" "<<coarsen_lbound<<" "<<coarsen_ubound<<std::endl;
	unsigned int cell_itr = 0;
	for(; cell != endc; ++cell, ++density_cell){
		unsigned int n_qpoints = cellprop[cell_itr].material_density.size();
		double max_density = 1, min_density = 0;
		double total_weight = 0.0;
		double density = 0.0;
		unsigned int quad_index = elastic_data.get_quad_index(cellprop[cell_itr].quadrature_formula);
		for(unsigned int k = 0; k < elastic_data.JxW[quad_index].size(); ++k){
			total_weight += elastic_data.JxW[quad_index][k];
		}
		for(unsigned int qpoint = 0; qpoint < n_qpoints; ++qpoint){
			density += cellprop[cell_itr].xPhys[qpoint] *
								(elastic_data.JxW[quad_index][qpoint]/total_weight);

		}
		if(density > coarsen_ubound || density < coarsen_lbound){
			cell->set_coarsen_flag();
			density_cell->set_coarsen_flag();
		}
		else if (density > refine_lbound && density < refine_ubound) {
			cell->set_refine_flag();
			density_cell->set_refine_flag();
		}
		++cell_itr;
	}
	triangulation.prepare_coarsening_and_refinement();
	density_triangulation.prepare_coarsening_and_refinement();
	cell_itr = 0;
	cell = triangulation.begin_active(),
			endc = triangulation.end();
	for(; cell != endc; ++cell){
		//clear any previous user index
		cell->clear_user_index();

		//update the user indices
		if(cell->refine_flag_set()){
			cell->set_user_index(cell_itr + 1);
		}
		else if(cell->coarsen_flag_set()){
			cell->parent()->set_user_index(cell_itr + 1);
		}
		else{
			cell->set_user_index(cell_itr + 1);
		}
		++cell_itr;
	}

}

void FEIndicator::get_fe_indicator(Triangulation<2> &triangulation){

}

