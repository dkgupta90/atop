/*
 * DensityValues.cpp
 *
 *  Created on: Jul 14, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <topopt/TopologyOptimization/DensityValues.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <assert.h>
#include <deal.II/base/point.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <topopt/TopologyOptimization/neighbors.h>
#include <topopt/TopologyOptimization/cell_prop.h>
#include <math.h>
using namespace topopt;
using namespace dealii;

void DensityValues::update_density_mesh(
		std::vector<CellProperties> &cellprop,
		std::vector<double> &density_mesh){
	density_mesh.clear();
	unsigned int cell_count = cellprop.size();
	for(unsigned int i = 0; i < cell_count; ++i){
		for(unsigned int j = 0; j < cellprop[i].material_density.size(); ++j){
			density_mesh.push_back(cellprop[i].material_density[j]);
		}
	}

}

void DensityValues::get_vol_fraction(
		std::vector<CellProperties> &cellprop,
		std::vector<double> &density_mesh,
		double &density_sum,
		double &max_cell_area,
		StoreElasticData &elastic_data){
	unsigned int design_itr = 0;
	density_sum = 0.0;
	unsigned int no_cells = cellprop.size();
	for(unsigned int i = 0; i < no_cells; ++i){
		double temp_density = 0.0;
		unsigned int n_q_points = cellprop[i].material_density.size();
		unsigned int quadrature_rule = cellprop[i].quadrature_formula - 1;
		double total_weight = 0.0;
		for(unsigned int k = 0; k < elastic_data.JxW[quadrature_rule].size(); ++k){
			total_weight += elastic_data.JxW[quadrature_rule][k];
		}
		//std::cout<<total_weight<<std::endl;
		for(unsigned int  j = 0;  j  < n_q_points; ++j){
			//std::cout<<density_mesh[design_itr]<<std::endl;
			cellprop[i].material_density[j]  = density_mesh[design_itr];
			temp_density += (elastic_data.JxW[quadrature_rule][j]/total_weight) * cellprop[i].material_density[j];
			++design_itr;
		}
		double area_factor = cellprop[i].cell_area/max_cell_area;
		density_sum += temp_density * area_factor;
	}
}

void DensityValues::create_neighbours(
		std::vector<CellProperties>  &cellprop,
		FESystem<2> &fe,
		DoFHandler<2> &dof_handler,
				double rmin){
	unsigned int cell_itr1 = 0;
	typename DoFHandler<2>::active_cell_iterator cell1 = dof_handler.begin_active(),
				endc1= dof_handler.end();
	for(; cell1 != endc1; ++cell1){
		QGauss<2> quadrature_formula1(cellprop[cell_itr1].quadrature_formula);
		FEValues<2> fe_values1(fe,
				quadrature_formula1,
				update_values |
				update_gradients |
				update_quadrature_points |
				update_JxW_values
				);
		fe_values1.reinit(cell1);
		std::vector<std::pair<unsigned int, unsigned int> > block_indices;
		block_indices.clear();
		std::vector<DoFHandler<2>::active_cell_iterator> neighbor_iterators;
		neighbor_iterators.clear();
		double drmin = rmin + sqrt(cell1->measure()/2);
		neighbor_iterators.push_back(cell1);
		neighbor_search(cell1, cell1, neighbor_iterators, drmin);
		//std::cout<<"Cell : "<<cell_itr1<<"      No. of neighbors : "<<neighbor_iterators.size()<<std::endl;
		if(neighbor_iterators.size() == 0){
			std::cout<<"Strange condition : NO NEIGHBOR FOUND  for cell : "<<cell_itr1<<std::endl;
		}
		std::vector<Point<2> > qpoints1 = fe_values1.get_quadrature_points();
		cellprop[cell_itr1].neighbour_cells.resize(qpoints1.size());
		cellprop[cell_itr1].neighbour_distance.resize(qpoints1.size());
		unsigned int cell_itr2;
/*		cell_itr2 = cell_itr1;
		typename DoFHandler<2>::active_cell_iterator cell2 = dof_handler.begin_active(),
					endc2= dof_handler.end();*/
		typename DoFHandler<2>::active_cell_iterator cell2;
		//for(cell2 = cell1; cell2 != endc2; ++cell2){
		for(unsigned int ng_itr = 0;  ng_itr < neighbor_iterators.size(); ++ng_itr){
			cell2 = neighbor_iterators[ng_itr];
			cell_itr2 = cell2->user_index() - 1;
			QGauss<2> quadrature_formula2(cellprop[cell_itr2].quadrature_formula);
			FEValues<2> fe_values2(fe,
					quadrature_formula2,
					update_values |
					update_gradients |
					update_quadrature_points |
					update_JxW_values
					);
			fe_values2.reinit(cell2);
/*				double c2c_distance = cell1->center().distance(cell2->center());
			if(c2c_distance > (rmin + 0.04)){
				continue;
			}*/
			std::vector<Point<2> > qpoints2 = fe_values2.get_quadrature_points();
			cellprop[cell_itr2].neighbour_cells.resize(qpoints2.size());
			cellprop[cell_itr2].neighbour_distance.resize(qpoints2.size());
			double distance;
			//Checking the distance here
			unsigned int rminlevel;
			if(cell1->level() > 0) {
				rminlevel = cell1->level();
			}
			else{
				rminlevel = 0;
			}
			double exfactor1= pow(0.6, (rminlevel));
			double exfactor2 = pow(0.6, (rminlevel));
			double rmin1, rmin2;
			rmin1 = rmin * exfactor1;
			rmin2 = rmin * exfactor2;
			//std::cout<<rmin1<<"    "<<rmin2<<std::endl;
			for(unsigned int q_point1 = 0; q_point1 < qpoints1.size(); ++q_point1){
				for(unsigned int q_point2 = 0; q_point2 < qpoints2.size(); ++q_point2){
					distance  = 0.0;
					distance = qpoints1[q_point1].distance(qpoints2[q_point2]);
					if(distance > rmin1 && distance > rmin2){
						continue;
					}
					//std::cout<<distance<<"   "<<rmin<<std::endl;
					std::pair<unsigned int, unsigned int> temp_pair;
					if(distance < rmin1){
						temp_pair.first = cell_itr2;
						temp_pair.second = q_point2;
						cellprop[cell_itr1].neighbour_cells[q_point1].push_back(temp_pair);
						cellprop[cell_itr1].neighbour_distance[q_point1].push_back(distance);
					}
/*					if(cell_itr1 == cell_itr2) continue;
					if(distance < rmin2){
						temp_pair.first = cell1;
						temp_pair.second = q_point1;
						cellprop[cell_itr2].neighbour_cells[q_point2].push_back(temp_pair);
						cellprop[cell_itr2].neighbour_distance[q_point2].push_back(distance);
					}*/
				}
			}
			//++cell_itr2;
		}

		//std::cout<<cell_itr1<<"    "<<qpoints1.size()<<"   "<<cellprop[cell_itr1].neighbour_cells[0].size()<<std::endl;
		calculate_weights(cellprop,
				cell_itr1,
				rmin);
		++cell_itr1;
	}
}

void DensityValues::neighbor_search(DoFHandler<2>::active_cell_iterator cell1,
		DoFHandler<2>::active_cell_iterator cell,
		std::vector<DoFHandler<2>::active_cell_iterator> &neighbor_iterators,
		double rmin
		){
	for(unsigned int iface = 0; iface < GeometryInfo<2>::faces_per_cell; ++iface){
		if(cell->at_boundary(iface)) continue;
		if(cell->neighbor(iface)->active()){
			if(cell1->center().distance(cell->neighbor(iface)->center()) < rmin){
				if(find(neighbor_iterators.begin(), neighbor_iterators.end(), cell->neighbor(iface)) == neighbor_iterators.end()){
					neighbor_iterators.push_back(cell->neighbor(iface));
					neighbor_search(cell1, cell->neighbor(iface), neighbor_iterators, rmin);
				}
			}
		}
		else{
			for(unsigned int ichildface = 0; ichildface < cell->face(iface)->n_children(); ++ichildface){
				if(cell1->center().distance(cell->neighbor_child_on_subface(iface, ichildface)->center()) < rmin){
					if(find(neighbor_iterators.begin(), neighbor_iterators.end(), cell->neighbor_child_on_subface(iface, ichildface)) == neighbor_iterators.end()){
						neighbor_iterators.push_back(cell->neighbor_child_on_subface(iface, ichildface));
						neighbor_search(cell1, cell->neighbor_child_on_subface(iface, ichildface), neighbor_iterators, rmin);
					}
				}
			}
		}
	}
}

void DensityValues::calculate_weights(std::vector<CellProperties> &cellprop,
		unsigned int cell_itr1,
		double rmin){
	unsigned int n_q_points = cellprop[cell_itr1].neighbour_distance.size();
	cellprop[cell_itr1].neighbour_weights.resize(n_q_points);
	for(unsigned int qpoint = 0; qpoint < n_q_points; ++qpoint){
		cellprop[cell_itr1].neighbour_weights[qpoint].resize(cellprop[cell_itr1].neighbour_distance[qpoint].size());
		double sum_weights = 0;
		for(unsigned int i = 0 ; i < cellprop[cell_itr1].neighbour_distance[qpoint].size(); ++i){
			double temp1 = rmin - cellprop[cell_itr1].neighbour_distance[qpoint][i];
			double area_factor = 1;//cellprop[cell_itr1].cell_area/max_cell_area;
			temp1 = temp1*area_factor;
			cellprop[cell_itr1].neighbour_weights[qpoint][i] = temp1;
			sum_weights += temp1;
		}
		//std::cout<<sum_weights<<std::endl;

		double xPhys = 0.0;

		for(unsigned int i = 0; i < cellprop[cell_itr1].neighbour_weights[qpoint].size(); ++i){
			if(sum_weights == 0){
				std::cerr<<"DensityValues::calculate_weights(..) : Divide by zero exception"<<std::endl;
			}
			else{
				cellprop[cell_itr1].neighbour_weights[qpoint][i] /= sum_weights;
			}
		}
	}
}

void DensityValues::filter(std::vector<CellProperties> &cellprop){
	unsigned int no_cells = cellprop.size();
	for(unsigned int cell_itr = 0 ; cell_itr < no_cells; ++cell_itr){
		for(unsigned int qpoint = 0 ; qpoint < cellprop[cell_itr].neighbour_cells.size(); ++qpoint){
			double xPhys = 0.0;
			unsigned int cell_itr2, qpoint2;
			for(unsigned int i = 0; i < cellprop[cell_itr].neighbour_weights[qpoint].size(); ++i){
				cell_itr2= cellprop[cell_itr].neighbour_cells[qpoint][i].first; //Denotes the cell iterator for the neighbor cell
				qpoint2 = cellprop[cell_itr].neighbour_cells[qpoint][i].second;  //denotes the q_point for the neighbor cell
				xPhys += cellprop[cell_itr].neighbour_weights[qpoint][i] * cellprop[cell_itr2].material_density[qpoint2];
			}
			cellprop[cell_itr].xPhys[qpoint] = xPhys;
		}
	}

}

double DensityValues::get_dxPhys_dx(CellProperties &cellp,
		unsigned int q_point,
		unsigned int cell_itr2,
		unsigned int q_point2){
	unsigned int n_cell, n_qpoint;
	for(unsigned int i = 0 ; i < cellp.neighbour_cells[q_point2].size(); ++i){
		n_cell = cellp.neighbour_cells[q_point2][i].first;
		n_qpoint = cellp.neighbour_cells[q_point2][i].second;
		if(n_cell == cell_itr2 && n_qpoint == q_point){
			double output = cellp.neighbour_weights[q_point2][i];
			return output;
		}
	}
}



