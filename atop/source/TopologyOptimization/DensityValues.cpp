/*
 * DensityValues.cpp
 *
 *  Created on: Jul 14, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <atop/TopologyOptimization/DensityValues.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <assert.h>
#include <deal.II/base/point.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/grid/grid_tools.h>
#include <atop/TopologyOptimization/neighbors.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <assert.h>

#include <stdlib.h>

#include <math.h>
using namespace topopt;
using namespace dealii;
using namespace atop;

//---------------------DENSITYFIELD CLASS MODULES-------------------------------------------------------------------

template <int dim>
void DensityField<dim>::create_neighbors(
		std::vector<CellInfo> &cell_info_vector,
		hp::FEValues<dim> &hp_fe_values,
		hp::DoFHandler<dim> &dof_handler,
		hp::DoFHandler<dim> &density_dof_handler,
		Projection &projection,
		DefineMesh<dim> &mesh
		){


	/*
	 * Iterate over all the cells to check for neighbors
	 * iterated over the cells in triangulation
	 */
	unsigned int cell_itr1 = 0;
	typename hp::DoFHandler<dim>::active_cell_iterator cell1 = dof_handler.begin_active(),
				endc1= dof_handler.end();
	for(; cell1 != endc1; ++cell1){

		//Just a random check, can be deleted
		if(cell_info_vector[cell_itr1].quad_rule == 0){
			exit(0);
		}

		hp_fe_values.reinit(cell1,
				cell_info_vector[cell_itr1].quad_rule - 1);	//to relate to quad_index
		const FEValues<dim> &fe_values1 = hp_fe_values.get_present_fe_values();

		//Stores cell iterators for all neighbors of the current cell
		std::vector<hp::DoFHandler<2>::active_cell_iterator> neighbor_iterators;
		neighbor_iterators.clear();


		//Computing the cell specific filter radius
		double proj_radius = cell_info_vector[cell_itr1].projection_radius;
		double drmin = proj_radius + sqrt(cell1->measure()); //added term is the distance from center of square element to the corner


		//Computing the cell specific filter radius
		//double proj_radius = projection.radius * pow(projection.gamma, (double)(cell1->level()));
		//double drmin = proj_radius + sqrt(cell1->measure()/2); //added term is the distance from center of square element to the corner


		//If the two meshes are decoupled, and design points are distributed w.r.t analysis mesh
		if (mesh.coupling == false && mesh.adaptivityType == "adaptive_grayness"){

			//Computing the cell specific filter radius
			//proj_radius = cell_info_vector[cell_itr1].projection_radius;
			//drmin = proj_radius + sqrt(cell1->measure()/2); //added term is the distance from center of square element to the corner

			//The following function gets the neighbors of the current cell lying within a distance of drmin
			neighbor_iterators.push_back(cell1);
			neighbor_search(cell1, cell1, neighbor_iterators, drmin);

			//std::cout<<"Cell : "<<cell_itr1<<"      No. of neighbors : "<<neighbor_iterators.size()<<std::endl;
			if(neighbor_iterators.size() == 0){
				std::cout<<"Strange condition : NO NEIGHBOR FOUND  for cell : "<<cell_itr1<<std::endl;
			}
			std::vector<Point<dim> > qpoints1 = fe_values1.get_quadrature_points();
			cell_info_vector[cell_itr1].neighbour_points.clear();
			cell_info_vector[cell_itr1].neighbour_distance.clear();
			cell_info_vector[cell_itr1].neighbour_cell_area_fraction.clear();
			cell_info_vector[cell_itr1].neighbour_points.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_distance.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_cell_area_fraction.resize(qpoints1.size());

			//Defining the virtual cell area
			cell_info_vector[cell_itr1].cell_area_fraction = (1.0/((double)cell_info_vector[cell_itr1].pseudo_design_points.no_points));


			unsigned int cell_itr2;
			typename hp::DoFHandler<2>::active_cell_iterator cell2;
			//Iterate over all neighboring cells to check distance with Gauss points
			for(unsigned int ng_itr = 0;  ng_itr < neighbor_iterators.size(); ++ng_itr){
				cell2 = neighbor_iterators[ng_itr];
				cell_itr2 = cell2->user_index() - 1;

				double distance;

				double rmin1;
				rmin1 = proj_radius;

				for(unsigned int q_point1 = 0; q_point1 < qpoints1.size(); ++q_point1){
					//Iterating over all the psuedo-design points of cell 2
					unsigned int ng_no_points = cell_info_vector[cell_itr2].pseudo_design_points.no_points;
					for (unsigned int ngpt_itr = 0; ngpt_itr < ng_no_points; ++ngpt_itr){

							Point<dim> point2;
							//converting vector to point coordinates
							Point<dim> centroid = cell2->center();	//getting the centre for scaling the points
							double side_length = pow(cell2->measure(), 1.0/dim);
							for(unsigned int dimi = 0; dimi < dim; ++dimi){
								point2(dimi) = centroid(dimi) +
								(cell_info_vector[cell_itr2].pseudo_design_points.pointX[ngpt_itr][dimi]) * (side_length/2.0);
							}

						distance = 0.0;
						distance = qpoints1[q_point1].distance(point2);

						if(distance > rmin1){
							continue;
						}

						//std::cout<<"And I reached here"<<std::endl;
						//Adding the point to the neighbour vector
						cell_info_vector[cell_itr1].neighbour_points[q_point1].push_back(
								std::make_pair(cell_itr2, ngpt_itr));

						//Adding the respective distance
						cell_info_vector[cell_itr1].neighbour_distance[q_point1].push_back(distance);

						//Adding the virtual area fraction
						cell_info_vector[cell_itr1].neighbour_cell_area_fraction[q_point1].push_back(1.0/((double)cell_info_vector[cell_itr2].pseudo_design_points.no_points));


					}
				}
			}

			//computing the respective weights
			calculate_weights(cell_info_vector,
					cell_itr1,
					proj_radius,
					mesh);
			++cell_itr1;
		}
		else{
/*			//Indentifying the density cell which contains the centroid of this FE cell
			typename hp::DoFHandler<dim>::active_cell_iterator density_cell1;
			//density_cell1 = GridTools::find_active_cell_around_point(density_dof_handler, cell1->center());
			density_cell1 = cell1;

			//The following function gets the neighbors of the current cell lying within a distance of drmin
			neighbor_iterators.push_back(density_cell1);
			neighbor_search(density_cell1, density_cell1, neighbor_iterators, drmin);

			//std::cout<<"Cell : "<<cell_itr1<<"      No. of neighbors : "<<neighbor_iterators.size()<<std::endl;
			if(neighbor_iterators.size() == 0){
				std::cout<<"Strange condition : NO NEIGHBOR FOUND  for cell : "<<cell_itr1<<std::endl;
			}
			std::vector<Point<dim> > qpoints1 = fe_values1.get_quadrature_points();
			cell_info_vector[cell_itr1].neighbour_cells.clear();
			cell_info_vector[cell_itr1].neighbour_distance.clear();
			cell_info_vector[cell_itr1].neighbour_cell_area.clear();
			cell_info_vector[cell_itr1].neighbour_cell_area_fraction.clear();

			cell_info_vector[cell_itr1].neighbour_cells.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_distance.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_cell_area.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_cell_area_fraction.resize(qpoints1.size());


			unsigned int density_cell_itr2;
			typename hp::DoFHandler<2>::active_cell_iterator density_cell2;
			//Iterate over all neighboring cells to check distance with Gauss points
			for(unsigned int ng_itr = 0;  ng_itr < neighbor_iterators.size(); ++ng_itr){

				density_cell2 = neighbor_iterators[ng_itr];
				density_cell_itr2 = density_cell2->user_index() - 1;
				if (fabs(cell_info_vector[density_cell_itr2].cell_density - 100) < 1e-12){
					//This one is only for the final output design in dp-adaptivity
					continue;
				}

				double distance;

				double exfactor1= 1;
				double rmin1;

				rmin1 = proj_radius * exfactor1;

				for(unsigned int q_point1 = 0; q_point1 < qpoints1.size(); ++q_point1){
						distance  = 0.0;
						distance = qpoints1[q_point1].distance(density_cell2->center());
						if(distance > rmin1){
							continue;
						}
						cell_info_vector[cell_itr1].neighbour_cells[q_point1].push_back(density_cell_itr2);
						cell_info_vector[cell_itr1].neighbour_distance[q_point1].push_back(distance);
						cell_info_vector[cell_itr1].neighbour_cell_area[q_point1].push_back(density_cell2->measure());
						cell_info_vector[cell_itr1].neighbour_cell_area_fraction[q_point1].push_back(cell_info_vector[density_cell_itr2].cell_area_fraction);
				}
			}


			//std::cout<<cell_itr1<<"    "<<qpoints1.size()<<"   "<<cellprop[cell_itr1].neighbour_cells[0].size()<<std::endl;
			calculate_weights(cell_info_vector,
					cell_itr1,
					proj_radius,
					mesh);
			++cell_itr1;*/
		}
	}
}

//---------------------DENSITYFIELD CLASS MODULES-------------------------------------------------------------------

template <int dim>
void DensityField<dim>::create_neighbors(
		std::vector<CellInfo> &design_cell_info_vector,
		hp::DoFHandler<dim> &design_handler,
		Projection &projection
		){
	/*
	 * Iterate over all the design cells to find the neighbors
	 */

	unsigned int design_itr1 = 0;
	double rmin1 = 1.00000001 * projection.true_radius;	//True radius based on the resolution of design mesh

	typename hp::DoFHandler<dim>::active_cell_iterator design_cell1 = design_handler.begin_active(),
				design_endc1= design_handler.end();
	for(; design_cell1 != design_endc1; ++design_cell1){

		//Stores cell iterators for all neighbors of the current cell
		std::vector<hp::DoFHandler<2>::active_cell_iterator> neighbor_iterators;
		neighbor_iterators.clear();

		//Getting the neighbors within the mentioned rmin1
		neighbor_iterators.push_back(design_cell1);
		neighbor_search(design_cell1, design_cell1, neighbor_iterators, rmin1);

		assert(neighbor_iterators.size() > 0);	/*Check that atleast one neighbor exists*/

		//Check all cell iterators and save the necessary information
		design_cell_info_vector[design_itr1].neighbor_cells.clear();
		design_cell_info_vector[design_itr1].neighbor_distance.clear();
		design_cell_info_vector[design_itr1].neighbor_weights.clear();
		typename hp::DoFHandler<dim>::active_cell_iterator design_cell2;
		for (unsigned int ng_itr = 0; ng_itr < neighbor_iterators.size(); ++ng_itr){

			/*Adding the neighbor design cell iterator index*/
			design_cell2 = neighbor_iterators[ng_itr];
			unsigned int design_itr2 = design_cell2->user_index()-1;
			design_cell_info_vector[design_itr1].neighbor_cells.push_back(design_itr2);

			/*Adding the center-center distance with the current neighbor design cell*/
			double distance = (design_cell1->center()).distance(design_cell2->center());
			design_cell_info_vector[design_itr1].neighbor_distance.push_back(distance);
		}

		/*computing the respective weights for smoothing*/
		calculate_weights(design_cell_info_vector,
				design_itr1,
				projection.true_radius);

		++design_itr1;
	}
}

template <int dim>
void DensityField<dim>::neighbor_search(hp::DoFHandler<2>::active_cell_iterator cell1,
		hp::DoFHandler<2>::active_cell_iterator cell,
		std::vector<hp::DoFHandler<2>::active_cell_iterator> &neighbor_iterators,
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

/*
 * This function calculates weights for the case where the design points cannot move in the domain
 * Currently implemented for only coupled mesh, however, it can be extended to decoupled mesh as well.
 */
template <int dim>
void DensityField<dim>::calculate_weights(std::vector<CellInfo> &cell_info_vector,
		unsigned int cell_itr1,
		double rmin,
		DefineMesh<dim> &mesh){
	unsigned int n_q_points = cell_info_vector[cell_itr1].neighbour_distance.size();
	cell_info_vector[cell_itr1].neighbour_weights.resize(n_q_points);
	for(unsigned int qpoint = 0; qpoint < n_q_points; ++qpoint){
		cell_info_vector[cell_itr1].neighbour_weights[qpoint].resize(cell_info_vector[cell_itr1].neighbour_distance[qpoint].size());
		double sum_weights = 0;
		for(unsigned int i = 0 ; i < cell_info_vector[cell_itr1].neighbour_distance[qpoint].size(); ++i){
			double temp1 = rmin - cell_info_vector[cell_itr1].neighbour_distance[qpoint][i];

			//check for mesh type
			double area_factor = 1;
			if(mesh.coupling == true){
				//area_factor = cell_info_vector[cell_itr1].neighbour_cell_area[qpoint][i]/max_cell_area;
				area_factor = cell_info_vector[cell_itr1].neighbour_cell_area_fraction[qpoint][i];
			}
			else{
				area_factor = cell_info_vector[cell_itr1].neighbour_cell_area_fraction[qpoint][i];
				//std::cout<<area_factor<<std::endl;

			}
			temp1 = temp1*area_factor;
			cell_info_vector[cell_itr1].neighbour_weights[qpoint][i] = temp1;
			sum_weights += temp1;
		}

		for(unsigned int i = 0; i < cell_info_vector[cell_itr1].neighbour_weights[qpoint].size(); ++i){
			if(sum_weights == 0){
				std::cerr<<"atop::DensityField:calculate_weights(..) : Divide by zero exception"<<std::endl;
			}
			else{
				//Corresponds to dx
				cell_info_vector[cell_itr1].neighbour_weights[qpoint][i] /= sum_weights;
			}
		}
	}
}

template <int dim>
void DensityField<dim>::calculate_weights(std::vector<CellInfo> &design_cell_info_vector,
		unsigned int design_itr1,
		double rmin){

	unsigned int no_neighbors = design_cell_info_vector[design_itr1].neighbor_cells.size();
	design_cell_info_vector[design_itr1].neighbor_weights.resize(no_neighbors);

	/*Calculating the weights*/
	double sum_weights = 0.0;
	for(unsigned int i = 0; i < no_neighbors; ++i){
		double temp1 = rmin - design_cell_info_vector[design_itr1].neighbor_distance[i];

		design_cell_info_vector[design_itr1].neighbor_weights[i] = temp1;
		sum_weights += temp1;
	}

	/*Normalizing the weights*/
	for (unsigned int i = 0; i < no_neighbors; ++i){
		assert(sum_weights > 0.0);	//total sum of weights should be none zero
		design_cell_info_vector[design_itr1].neighbor_weights[i] /= sum_weights;
	}
}


/*//--------------------------------------------------------------------------------------------------------------
template <int dim>
void DensityField<dim>::smoothing(
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		DefineMesh<dim> &mesh
		){
	unsigned int no_cells = cell_info_vector.size();

	if (mesh.coupling == true){
		for(unsigned int cell_itr = 0 ; cell_itr < no_cells; ++cell_itr){
			for(unsigned int qpoint = 0 ; qpoint < cell_info_vector[cell_itr].neighbour_cells.size(); ++qpoint){
				double xPhys = 0.0;
				unsigned int density_cell_itr2;
				for(unsigned int i = 0; i < cell_info_vector[cell_itr].neighbour_weights[qpoint].size(); ++i){
					double tempi = density_cell_info_vector[density_cell_itr2].density[0];
					density_cell_itr2 = cell_info_vector[cell_itr].neighbour_cells[qpoint][i];
					xPhys += cell_info_vector[cell_itr].neighbour_weights[qpoint][i]
							  * density_cell_info_vector[density_cell_itr2].density[0];

				}
				cell_info_vector[cell_itr].density[qpoint] = xPhys;
			}
		}
	}
	else{

		for(unsigned int cell_itr = 0 ; cell_itr < no_cells; ++cell_itr){

			cell_info_vector[cell_itr].density.clear();
			cell_info_vector[cell_itr].density.resize(cell_info_vector[cell_itr].n_q_points);
			for(unsigned int qpoint = 0 ; qpoint < cell_info_vector[cell_itr].neighbour_points.size(); ++qpoint){
				double xPhys = 0.0;
				unsigned int cell_itr2, ng_pt_itr;
				for(unsigned int i = 0; i < cell_info_vector[cell_itr].neighbour_weights[qpoint].size(); ++i){
					cell_itr2 = cell_info_vector[cell_itr].neighbour_points[qpoint][i].first;
					ng_pt_itr = cell_info_vector[cell_itr].neighbour_points[qpoint][i].second;
					xPhys += cell_info_vector[cell_itr].neighbour_weights[qpoint][i]
					          * cell_info_vector[cell_itr2].pseudo_design_points.rho[ng_pt_itr];
				}
				cell_info_vector[cell_itr].density[qpoint] = xPhys;
				//std::cout<<"xPhys : "<<xPhys<<std::endl;
			}
		}
	}
}*/

template <int dim>
void DensityField<dim>::smoothing(
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &design_cell_info_vector){
	unsigned int no_design_cells = design_cell_info_vector.size();

	for (unsigned int cell_itr = 0; cell_itr < cell_info_vector.size(); ++cell_itr){
		cell_info_vector[cell_itr].density.clear();
	}


	for(unsigned int design_itr = 0 ; design_itr < no_design_cells; ++design_itr){
		double xPhys = 0.0;
		unsigned int design_cell_itr2;
		//std::cout<<design_cell_info_vector[design_itr].neighbour_weights.size()<<std::endl;
		for(unsigned int i = 0; i < design_cell_info_vector[design_itr].neighbor_weights.size(); ++i){
			design_cell_itr2 = design_cell_info_vector[design_itr].neighbor_cells[i];
			xPhys += design_cell_info_vector[design_itr].neighbor_weights[i]
					  * design_cell_info_vector[design_cell_itr2].cell_density;

		}
		design_cell_info_vector[design_itr].filtered_density = xPhys;
		//std::cout<<"filtered density "<<xPhys<<std::endl;
		unsigned int cell_itr = design_cell_info_vector[design_itr].connected_cell_iterators_2D[0]->user_index() - 1;
		cell_info_vector[cell_itr].density.push_back(xPhys);
	}
}

/*//----------------------------------------------------------------
//Only for the output design mesh...used in create_design.cpp
template <int dim>
void DensityField<dim>::smoothing(
		std::vector<CellInfo> &design_info_vector){
	unsigned int no_cells = design_info_vector.size();

	for(unsigned int cell_itr = 0 ; cell_itr < no_cells; ++cell_itr){
		for(unsigned int qpoint = 0 ; qpoint < design_info_vector[cell_itr].neighbour_cells.size(); ++qpoint){
			double xPhys = 0.0;
			unsigned int density_cell_itr2;
			for(unsigned int i = 0; i < design_info_vector[cell_itr].neighbour_weights[qpoint].size(); ++i){
				density_cell_itr2 = design_info_vector[cell_itr].neighbour_cells[qpoint][i];
				xPhys += design_info_vector[cell_itr].neighbour_weights[qpoint][i]
						  * design_info_vector[density_cell_itr2].cell_density;

			}
			design_info_vector[cell_itr].density[qpoint] = xPhys;
		}
	}
}*/

//--------------------------------------------------------------------------------------------------------------------
/**
 * This function updates the design_vector
 * Design_vector is used for for the optimization purpose
 */
template <int dim>
void DensityField<dim>::update_design_vector(
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		std::vector<double> &design_vector,
		unsigned int cycle,
		double volfrac,
		DefineMesh<dim> &mesh,
		Projection &projection){


	unsigned int cell_count = density_cell_info_vector.size();

	//Update the design vector

	/*
	 * Below is the case of a coupled mesh where the same element is used for analysis as well design.
	 * For clarity purpose two different meshes are used. However, a one-one correlation exists between
	 * elements of the 2 meshes.
	 */
	if (mesh.coupling == true){
		design_vector.clear();
		if (cycle == 0){
			design_vector.resize(cell_count, volfrac);



		}
		else{
			for(unsigned int i = 0; i < cell_count; ++i){
				for(unsigned int j = 0; j < density_cell_info_vector[i].density.size(); ++j){
					design_vector.push_back(density_cell_info_vector[i].density[j]);
				}
			}
		}
	}
	else{
		if (cycle == 0){
			unsigned int no_design_points = design_vector.size();

			design_vector.clear();
			design_vector.resize(no_design_points, volfrac);
/*
			for (unsigned int i = 0; i < design_vector.size(); ++i){
				if ((i < 32)){
					design_vector[i] = 1.0;
				}
			}*/
/*			design_vector[63]*=1.5;
			design_vector[56]*=1.5;*/

		}
		else{
			design_vector.clear();
			design_vector.resize(density_cell_info_vector.size());
			for (unsigned int i = 0; i < cell_info_vector.size(); ++i){
				for (unsigned int j = 0; j < cell_info_vector[i].design_points.no_points; j++){
					unsigned int design_itr = (cell_info_vector[i].connected_cell_iterators_2D[j])->user_index() - 1;
					design_vector[design_itr] = cell_info_vector[i].design_points.rho[j];
				}
			}
		}
	}
}

//--------------------------------------------------------------------------------------------------------------------------
/*
 * This function is used to update the values of the design bounds for different types of adaptive methods.
 * For the case of only density design bounds, it will be straightforward.
 */
template <int dim>
void DensityField<dim>::update_design_bounds(
		std::vector<double> &lb,
		std::vector<double> &ub,
		DefineMesh<dim> &mesh,
		Projection &projection,
		std::vector<double> &design_vector){

	if (mesh.coupling == true){

		//in this case the number of design variables is equal to the number of cells in the
		//density_trinagulation mesh
		unsigned int design_count = mesh.design_triangulation->n_active_cells();
		lb.resize(design_count, 0.0);
		ub.resize(design_count, 1.0);
	}

	else{
		unsigned int design_count = design_vector.size();
		lb.resize(design_count, 0.0);
		ub.resize(design_count, 1.0);


/*		for (unsigned int i = 0; i < lb.size(); ++i){
			if ((i%24 <= 7) || (i%24 >=16)){
				lb[i] = 1.0 - 1e-6;
				ub[i] = 1.0;
			}
		}*/




	}

}

//------------------------------------------------------------------------------------------------------
/*template <int dim>
double DensityField<dim>::get_dxPhys_dx(CellInfo &cell_info,
		unsigned int q_point,
		unsigned int density_cell_itr2){
	unsigned int n_cell;
	for(unsigned int i = 0 ; i < cell_info.neighbour_cells[q_point].size(); ++i){
		n_cell = cell_info.neighbour_cells[q_point][i];
		if(n_cell == density_cell_itr2){
			double output = cell_info.neighbour_weights[q_point][i];
			return output;
		}
	}
	return 0;
}*/

/*template <int dim>
double DensityField<dim>::get_dxPhys_dx(CellInfo &cell_info,
		unsigned int q_point,
		unsigned int cell_itr2,
		unsigned int ngpt_itr){
	unsigned int n_cell, ng_pt;
	for(unsigned int i = 0; i < cell_info.neighbour_points[q_point].size(); ++i){
		n_cell = cell_info.neighbour_points[q_point][i].first;
		ng_pt = cell_info.neighbour_points[q_point][i].second;
		if(n_cell == cell_itr2 && ng_pt == ngpt_itr){
			double output = cell_info.neighbour_weights[q_point][i];
			return output;
		}

	}
}*/


template <int dim>
double DensityField<dim>::get_vol_fraction(
		std::vector<CellInfo> &design_cell_info_vector
		){
	double volume = 0.0;
	for(unsigned int i = 0; i < design_cell_info_vector.size(); ++i){

		volume += (design_cell_info_vector[i].filtered_density);
	}
	volume /= design_cell_info_vector.size();
	return volume;
}

template <int dim>
void DensityField<dim>::update_density_cell_info_vector(
		std::vector<CellInfo> &density_cell_info_vector,
		const std::vector<double> &design_vector){
	/*
	 * For the case of coupled mesh, for now, the size of density_cell_info_vector and design_vector will be same.
	 * This check will be used to check the type of adaptivity and accordingly update the vector
	 */

	if (density_cell_info_vector.size() == design_vector.size()){
		//This is the case of the coupled mesh approach and the design points do not move.
		unsigned int cell_count = density_cell_info_vector.size();
		unsigned int design_count = 0;
		for(unsigned int i = 0; i < cell_count; ++i){
			for(unsigned int j = 0; j < density_cell_info_vector[i].density.size(); ++j){
				density_cell_info_vector[i].density[j] = design_vector[design_count];
				++design_count;
			}
		}
	}
	else if (((dim + 2) * density_cell_info_vector.size()) == design_vector.size()){
		//This is the case for moving design points, here, every design point is associated with dim+2 design variables
		unsigned int k = 0;	//iterate over the design variables
		for (unsigned int i = 0 ; i < density_cell_info_vector.size(); ++i){
			density_cell_info_vector[i].density[0] = design_vector[k];	k++;	//updated design point density
			density_cell_info_vector[i].projection_fact = design_vector[k];	k++;	//updated projection radius factor

			//updating the density point location
			for (unsigned int j = 0; j < dim; j++){
				density_cell_info_vector[i].pointX[j] = design_vector[k];	k++;
			}
		}

		for (unsigned int i = 0; i < density_cell_info_vector.size(); ++i){
			std::cout<<density_cell_info_vector[i].density[0]<<" "<<
					density_cell_info_vector[i].projection_fact<<" "<<
					density_cell_info_vector[i].pointX[0]<<" "<<
					density_cell_info_vector[i].pointX[1]<<" "<<std::endl;
		}
	}


}

//Overloaded function for uncoupled meshes
template <int dim>
void DensityField<dim>::update_density_cell_info_vector(
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		const std::vector<double> &design_vector){

	//This is the case where the meshes are decoupled, but excludes the case of 'movingdesignpoints'

	//std::cout<<"Size of the design vector = "<<design_vector.size()<<std::endl;
	for (unsigned int i = 0; i < cell_info_vector.size(); ++i){
		cell_info_vector[i].design_points.rho.clear();
	}

	//Iterating over the design vector
	for (unsigned int design_itr = 0; design_itr < design_vector.size(); ++design_itr){
		unsigned int cell_itr = density_cell_info_vector[design_itr].connected_cell_iterators_2D[0]->user_index()-1;

		cell_info_vector[cell_itr].design_points.rho.push_back(design_vector[design_itr]);

	}
}


//return size of the design_vector
template <int dim>
unsigned int DensityField<dim>::get_design_count(
		unsigned int cycle,
		DefineMesh<dim> &mesh,
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector){

	if (cycle == 0){
		if (mesh.coupling == false){
			return (cell_info_vector.size() * mesh.initial_dcount_per_el);
		}
		else{
			return (density_cell_info_vector.size() * mesh.design_var_per_point());
		}
	}
	else{
		if (mesh.coupling == false){
			//Calculating the no. of design variables from cell_info_vector
			unsigned int no_design_count = 0;
			for(unsigned int i = 0; i < cell_info_vector.size(); ++i){
				no_design_count += cell_info_vector[i].design_points.no_points;
			}
			return no_design_count;
		}
	}

	return 0;
}

/*
 * This function is meant to calculate the projected density values at all the quadrature points of a certain face
 */
template <int dim>
void DensityField<dim>::get_xPhys_for_face(std::vector<double> &face_xPhys,
		hp::FECollection<dim> &temp_fe_coll,
		hp::QCollection<dim-1> &temp_fq_coll,
		hp::DoFHandler<dim> &dofhandler,
		std::vector<CellInfo> &cell_info_vector,
		typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		unsigned int face_itr){

	hp::FEFaceValues<dim> hp_fe_face_values(temp_fe_coll,
			temp_fq_coll,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	unsigned int cell_itr = cell->user_index() - 1;
	unsigned int q_index = cell_info_vector[cell_itr].quad_rule - 1;
	hp_fe_face_values.reinit(cell, face_itr, q_index);
	const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();
	const std::vector<Point<dim> > face_q_points = fe_face_values.get_quadrature_points();

	double proj_radius = cell_info_vector[cell_itr].projection_radius;

	//Getting all the neighbor cells for the current cell
	double drmin = proj_radius + sqrt(cell->measure()); //added term is the distance from center of square element to the corner

	std::vector<hp::DoFHandler<2>::active_cell_iterator> neighbor_iterators;
	neighbor_iterators.clear();
	//The following function gets the neighbors of the current cell lying within a distance of drmin
	neighbor_iterators.push_back(cell);
	neighbor_search(cell, cell, neighbor_iterators, drmin);

	//std::cout<<"Cell : "<<cell_itr1<<"      No. of neighbors : "<<neighbor_iterators.size()<<std::endl;
	if(neighbor_iterators.size() == 0){
		std::cout<<"Strange condition : NO NEIGHBOR FOUND  for cell : "<<cell_itr<<"at face : "<<face_itr<<std::endl;
	}

	face_xPhys.clear();
	face_xPhys.resize(face_q_points.size());

	double sum_weighted_densities, sum_weights;

	//Iterating over all the quadrature points of the face
	for (unsigned int i = 0; i < face_xPhys.size(); ++i){

		sum_weighted_densities = 0.0;
		sum_weights = 0.0;

		Point<dim> point1 = face_q_points[i];

		//Iterate over all the cells
		for (unsigned int j = 0; j < neighbor_iterators.size(); ++j){
			unsigned int ng_cell_itr = neighbor_iterators[j]->user_index() - 1;

			//Iterate over all the pseudo-design points
			unsigned int no_points = cell_info_vector[ng_cell_itr].pseudo_design_points.no_points;
			for (unsigned int pt = 0; pt < no_points; ++pt){

				Point<dim> point2;
				//converting vector to point coordinates
				Point<dim> centroid = neighbor_iterators[j]->center();	//getting the centre for scaling the points
				double side_length = pow(neighbor_iterators[j]->measure(), 1.0/dim);
				for(unsigned int dimi = 0; dimi < dim; ++dimi){
					point2(dimi) = centroid(dimi) +
					(cell_info_vector[ng_cell_itr].pseudo_design_points.pointX[pt][dimi]) * (side_length/2.0);
				}

				double distance = point1.distance(point2);
				if (distance > proj_radius)	continue;

				double diff_len = proj_radius - distance;
				double rho = cell_info_vector[ng_cell_itr].pseudo_design_points.rho[pt];

				sum_weighted_densities += (diff_len * rho);
				sum_weights += diff_len;

			}
		}

		if (fabs(sum_weights) < 1e-12)	std::cerr<<"Zero weights error ..."<<std::endl;
		sum_weighted_densities /= sum_weights;
		face_xPhys[i] = sum_weighted_densities;
	}
}
