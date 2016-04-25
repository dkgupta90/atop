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
#include <deal.II/fe/fe_values.h>
#include <atop/TopologyOptimization/neighbors.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <stdlib.h>

#include <math.h>
using namespace topopt;
using namespace dealii;
using namespace atop;

//---------------------DENSITYFIELD CLASS MODULES-------------------------------------------------------------------

template <int dim>
void DensityField<dim>::create_neighbors(
		std::vector<CellInfo> &cell_info_vector,
		FESystem<dim> &fe,
		FESystem<dim> &density_fe,
		DoFHandler<dim> &dof_handler,
		DoFHandler<dim> &density_dof_handler,
		Projection &projection,
		DefineMesh<dim> &mesh
		){

	/*
	 * Iterate over all the cells to check for neighbors
	 * iterated over the cells in triangulation
	 */
	unsigned int cell_itr1 = 0;
	typename DoFHandler<dim>::active_cell_iterator cell1 = dof_handler.begin_active(),
				endc1= dof_handler.end();
	for(; cell1 != endc1; ++cell1){

		//Just a random check, can be deleted
		if(cell_info_vector[cell_itr1].quad_rule == 0){
			exit(0);
		}
		QGauss<dim> quadrature_formula1(cell_info_vector[cell_itr1].quad_rule);
		FEValues<dim> fe_values1(fe,
				quadrature_formula1,
				update_values |
				update_gradients |
				update_quadrature_points |
				update_JxW_values
				);
		fe_values1.reinit(cell1);

		//Stores cell iterators for all neighbors of the current cell
		std::vector<DoFHandler<2>::active_cell_iterator> neighbor_iterators;
		neighbor_iterators.clear();

		//Computing the cell specific filter radius
		double proj_radius = projection.radius * pow(projection.gamma, (double)(cell1->level()));
		double drmin = proj_radius + sqrt(cell1->measure()/2); //added term is the distance from center of square element to the corner

		//If the two meshes are decoupled, and design points are distributed w.r.t analysis mesh
		if (mesh.coupling == false && mesh.adaptivityType == "adaptive_grayness"){

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
			cell_info_vector[cell_itr1].neighbour_points.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_distance.resize(qpoints1.size());

			unsigned int cell_itr2;
			typename DoFHandler<2>::active_cell_iterator cell2;
			//Iterate over all neighboring cells to check distance with Gauss points
			for(unsigned int ng_itr = 0;  ng_itr < neighbor_iterators.size(); ++ng_itr){
				cell2 = neighbor_iterators[ng_itr];
				cell_itr2 = cell2->user_index() - 1;

				double distance;

				double exfactor1= 1;	//no idea whether this one is needed or not
				double rmin1;
				rmin1 = proj_radius * exfactor1;
				for(unsigned int q_point1 = 0; q_point1 < qpoints1.size(); ++q_point1){

					//Iterating over all the density points of cell 2
					unsigned int ng_no_points = cell_info_vector[cell_itr2].design_points.no_points;
					for (unsigned int ngpt_itr = 0; ngpt_itr < ng_no_points; ++ngpt_itr){
							Point<dim> point2;
							//converting vector to point coordinates
							Point<dim> centroid = cell2->center();	//getting the centre for scaling the points
							double side_length = pow(cell2->measure(), 1.0/dim);
							for(unsigned int dimi = 0; dimi < dim; ++dimi){
								point2(dimi) = centroid(dimi) +
								(cell_info_vector[cell_itr2].design_points.pointX[ngpt_itr][dimi]) * (side_length/2.0);
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
			//Indentifying the density cell which contains the centroid of this FE cell
			typename DoFHandler<dim>::active_cell_iterator density_cell1;
			density_cell1 = GridTools::find_active_cell_around_point(density_dof_handler, cell1->center());

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
			cell_info_vector[cell_itr1].neighbour_cells.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_distance.resize(qpoints1.size());
			cell_info_vector[cell_itr1].neighbour_cell_area.resize(qpoints1.size());

			unsigned int density_cell_itr2;
			typename DoFHandler<2>::active_cell_iterator density_cell2;
			//Iterate over all neighboring cells to check distance with Gauss points
			for(unsigned int ng_itr = 0;  ng_itr < neighbor_iterators.size(); ++ng_itr){
				density_cell2 = neighbor_iterators[ng_itr];
				density_cell_itr2 = density_cell2->user_index() - 1;

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
				}
			}


			//std::cout<<cell_itr1<<"    "<<qpoints1.size()<<"   "<<cellprop[cell_itr1].neighbour_cells[0].size()<<std::endl;
			calculate_weights(cell_info_vector,
					cell_itr1,
					proj_radius,
					mesh);
			++cell_itr1;
		}
	}
}

//---------------------------------------------------------------------------------------------------------------
/*
 * This overloaded function is implemented for decoupled mesh with moving design points
 * density_triangulation should not be used in this instance
 *
 */
template <int dim>
void DensityField<dim>::create_neighbors(
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		FESystem<dim> &fe,
		DoFHandler<dim> &dof_handler
		){

	projection_matrix.clear();	//This is the regularization operator

	/*
	 * Clearing the cell_info_vector and initializing it
	 */
	for (unsigned int i = 0; i < cell_info_vector.size(); ++i){
		unsigned int n_qpoints = cell_info_vector[i].n_q_points;
		cell_info_vector[i].neighbour_cells.clear();
		cell_info_vector[i].neighbour_distance.clear();
		cell_info_vector[i].neighbour_cell_area.clear();
		cell_info_vector[i].neighbour_cells.resize(n_qpoints);
		cell_info_vector[i].neighbour_distance.resize(n_qpoints);
		cell_info_vector[i].neighbour_cell_area.resize(n_qpoints);
		cell_info_vector[i].sum_weights.resize(n_qpoints);
	}

	/*
	 * Calculating the projection radius for the current design point
	 * Here, we fix the radius w.r.t to the cell size of the initial triangulation
	 * However, this needs to changed in future to be adapted w.r.t refinement in the triangulation
	 */
	typename DoFHandler<dim>::active_cell_iterator cell1 = dof_handler.begin_active();
	double cell_len = sqrt(cell1->measure());
	this->cell_length = cell_len;	//for use in dxPhys_dx();
	/*
	 * Iterate over all the density points and checks the neighbor cells in triangulation
	 * iterated over the CellInfo vector for the design points
	 */
	unsigned int no_design_points = density_cell_info_vector.size();
	for (unsigned int i = 0; i < density_cell_info_vector.size(); i++){

		//For storing the cell iterators for all neighbor cells in triangulation
		std::vector<DoFHandler<2>::active_cell_iterator> neighbor_iterators;
		neighbor_iterators.clear();

		//Computing the projection radius of density point
		//double rmin = density_cell_info_vector[i].projection_fact * cell_len;
		double rmin = density_cell_info_vector[i].projection_fact * cell_len;	//scaled cut-off radius

		//Identifying the cell in the triangulation around the current density point
		Point<dim> des_pt;
		for (unsigned int k = 0; k < dim; k++)	des_pt[k] = density_cell_info_vector[i].pointX[k];

		cell1 = GridTools::find_active_cell_around_point(dof_handler, des_pt);
		double drmin = rmin + sqrt(cell1->measure()/2);	//for making sure all the cells are covered

		//Getting the neighbors of current cell within a distance of rmin
		neighbor_iterators.push_back(cell1);
		neighbor_search(cell1, cell1, neighbor_iterators, drmin);
		//std::cout<<"Density point: "<<i<<"  No. of neighbors : "<<neighbor_iterators.size()<<std::endl;

		//clearing the memory for saving information of the neighbors
		density_cell_info_vector[i].neighbour_cells.clear();
		density_cell_info_vector[i].neighbour_distance.clear();
		density_cell_info_vector[i].neighbour_cell_area.clear();

		//initializing the size of above vectors
		density_cell_info_vector[i].neighbour_cells.resize(1);
		density_cell_info_vector[i].neighbour_distance.resize(1);
		density_cell_info_vector[i].neighbour_cell_area.resize(1);	//1 since only one design point stored in one density_cell

		//Iterating over all the found out neighbor cells
		typename DoFHandler<2>::active_cell_iterator cell2;
		unsigned int cell_itr2;
		for (unsigned int ng_itr = 0; ng_itr < neighbor_iterators.size(); ++ng_itr){
			cell2 = neighbor_iterators[ng_itr];
			cell_itr2 = cell2->user_index() - 1;

			QGauss<dim> quadrature_formula1(cell_info_vector[cell_itr2].quad_rule);
			FEValues<dim> fe_values2(fe,
					quadrature_formula1,
					update_values |
					update_gradients |
					update_quadrature_points |
					update_JxW_values
					);
			fe_values2.reinit(cell2);
			std::vector<Point<dim> > qpoints2 = fe_values2.get_quadrature_points();	//Getting the quad points for the current neighbor cell

			double distance;
			//Iterating over all the Gauss points of the neighbor cell
			for (unsigned int qpoint2 = 0; qpoint2 < qpoints2.size(); ++qpoint2){
				distance = des_pt.distance(qpoints2[qpoint2]);	//distance of design point from the quad point
				if (distance > rmin){
					continue;
				}
				cell_info_vector[cell_itr2].neighbour_cells[qpoint2].push_back(i);
				cell_info_vector[cell_itr2].neighbour_distance[qpoint2].push_back(distance);
			}


		}
	}

	//Update the weights
	calculate_weights(cell_info_vector,
			density_cell_info_vector,
			cell_len);
}

template <int dim>
void DensityField<dim>::neighbor_search(DoFHandler<2>::active_cell_iterator cell1,
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
				area_factor = cell_info_vector[cell_itr1].neighbour_cell_area[qpoint][i]/max_cell_area;
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

/*
 * This implementation for calculating the weights for the 'movingDesignPoints' approach
 * Unlike the above one, here the weights w.r.t. whole triangulation is calculated at once.
 * It is NOT called again and again, but just ONCE for one iteration of optimization.
 */
template <int dim>
void DensityField<dim>::calculate_weights(std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		double cell_len){

	//Iterating over all the cells of the triangulation
	for (unsigned int cell_itr = 0; cell_itr < cell_info_vector.size(); ++cell_itr){
		unsigned int n_q_points = cell_info_vector[cell_itr].n_q_points;
		cell_info_vector[cell_itr].neighbour_weights.resize(n_q_points);

		//Iterating over all the Gauss points of the current cell
		for (unsigned int qpoint = 0; qpoint < n_q_points; ++qpoint){

			//setting the size of the empty weights vector for the current Gauss point
			cell_info_vector[cell_itr].neighbour_weights[qpoint].resize(cell_info_vector[cell_itr].neighbour_distance[qpoint].size());
			double sum_weights = 0.0;	//to store the total sum for using for division later

			//Iterating over all the neighbor DESIGN points for the current Gauss point
			for (unsigned int i = 0; i < cell_info_vector[cell_itr].neighbour_cells[qpoint].size(); ++i){
				unsigned int density_itr = cell_info_vector[cell_itr].neighbour_cells[qpoint][i];	//index of the neighbor
/*				double gamma = (density_cell_info_vector[i].projection_fact - 2.0) / 5.0;			//equivalent to the variance in gaussian function
				//double rmin = density_cell_info_vector[density_itr].projection_fact * cell_len;		//rmin for the current design point
				double temp1 = cell_info_vector[cell_itr].neighbour_distance[qpoint][i];
				double a = 1 / sqrt(2*3.14159*gamma);
				double b = exp(-(temp1 * temp1)/(2 * gamma * cell_len * cell_len));
				cell_info_vector[cell_itr].neighbour_weights[qpoint][i] = a * b;
				sum_weights += (a * b);*/

/*				//Below is the linear cone projection
				double rmin = density_cell_info_vector[density_itr].projection_fact * cell_len;		//rmin for the current design point
				double temp1 = rmin - cell_info_vector[cell_itr].neighbour_distance[qpoint][i];
				cell_info_vector[cell_itr].neighbour_weights[qpoint][i] = temp1;
				sum_weights += temp1;*/

				//Below is the cubic spline
				double rmin = density_cell_info_vector[density_itr].projection_fact * cell_len;		//rmin for the current design point
				double temp1 = cell_info_vector[cell_itr].neighbour_distance[qpoint][i];
				double r = temp1/rmin;
				double weight = 0.0;
				if (r <= 0.5)	weight = (2.0/3.0) - 4*r*r + 4*r*r*r;
				else if(r >0.5 && r <= 1.0)		weight = (4.0/3.0) - 4*r + 4*r*r - (4.0/3.0)*r*r*r;
				cell_info_vector[cell_itr].neighbour_weights[qpoint][i] = weight;
				sum_weights += weight;
			}
			//Normalizing the weights by the sum total of the weights for any Gauss point
			for (unsigned int i = 0; i < cell_info_vector[cell_itr].neighbour_weights[qpoint].size(); ++i){
				if (sum_weights == 0){
					std::cerr<<"atop::DensityField:calculate_weights(..) : Divide by zero exception"<<std::endl;
				}
				else{
					cell_info_vector[cell_itr].neighbour_weights[qpoint][i] /= sum_weights;
				}
			}
			cell_info_vector[cell_itr].sum_weights[qpoint] = sum_weights;
		}
	}
}

//--------------------------------------------------------------------------------------------------------------
template <int dim>
void DensityField<dim>::smoothing(
		std::vector<CellInfo> &cell_info_vector,
		std::vector<CellInfo> &density_cell_info_vector,
		DefineMesh<dim> &mesh
		){
	unsigned int no_cells = cell_info_vector.size();

	if (mesh.coupling == true || mesh.adaptivityType == "movingdesignpoints"){
		for(unsigned int cell_itr = 0 ; cell_itr < no_cells; ++cell_itr){
			for(unsigned int qpoint = 0 ; qpoint < cell_info_vector[cell_itr].neighbour_cells.size(); ++qpoint){
				double xPhys = 0.0;
				unsigned int density_cell_itr2;
				for(unsigned int i = 0; i < cell_info_vector[cell_itr].neighbour_weights[qpoint].size(); ++i){
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
			for(unsigned int qpoint = 0 ; qpoint < cell_info_vector[cell_itr].neighbour_points.size(); ++qpoint){
				double xPhys = 0.0;
				unsigned int cell_itr2, ng_pt_itr;
				for(unsigned int i = 0; i < cell_info_vector[cell_itr].neighbour_weights[qpoint].size(); ++i){
					cell_itr2 = cell_info_vector[cell_itr].neighbour_points[qpoint][i].first;
					ng_pt_itr = cell_info_vector[cell_itr].neighbour_points[qpoint][i].second;
					xPhys += cell_info_vector[cell_itr].neighbour_weights[qpoint][i]
					          * cell_info_vector[cell_itr2].design_points.rho[ng_pt_itr];
				}
				cell_info_vector[cell_itr].density[qpoint] = xPhys;
				//std::cout<<"xPhys : "<<xPhys<<std::endl;
			}
		}
	}


}

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
	 * This is a case where the design points are decoupled and they actually do not form a mesh.
	 * For convenience it is implemented as a separate mesh, however, there is only one point per mesh.
	 * No quadrature sort of implementation exists for the design mesh for such a case.
	 * In the worst implementation, q = 1 point might be assumed for quadrature, however, that should be avoided.
	 */
	if (mesh.coupling == false && mesh.adaptivityType == "movingdesignpoints"){
		design_vector.clear();
		if (cycle == 0){
			design_vector.resize(cell_count * mesh.design_var_per_point());

			typename Triangulation<dim>::active_cell_iterator cell = mesh.triangulation->begin_active(),
						endc = mesh.triangulation->end();

			for (unsigned int i = 0; i < design_vector.size();){
				double randno = (double)(rand() % 100 + 1);
				randno /= 1000.0;
				design_vector[i] = volfrac - randno; //adding the density for the current design point
				i++;
				design_vector[i] = projection.fact; //Adding the projection radius factor
				i++;
				// initializing the dim no. of coordinates for the location of the design point
				for (unsigned int j = 0; j < dim; j++){
					design_vector[i] = cell->center()[j];
/*					if (design_vector[i] > 1.5){
						design_vector[i] -= 1.5;
					}
					if (design_vector[i] > 1.0){
						design_vector[i] -= 1.0;
					}
					if (design_vector[i] > 0.5){
						design_vector[i] -= 0.5;
					}*/
					i++;
				}
				cell++;

			}
		}
		else{
			design_vector.clear();
			unsigned int k = 0; //iterate over the whole design vector
			for (unsigned int i = 0; i < cell_count; ++i){
				design_vector.push_back(density_cell_info_vector[i].density[0]);	//adding the density value
				design_vector.push_back(density_cell_info_vector[i].projection_radius);		//adding rmin value for the cell
				for(unsigned int j = 0; j < density_cell_info_vector[i].pointX.size(); ++j){
					//adding all the components of dim coordinates for the location
					design_vector.push_back(density_cell_info_vector[i].pointX[j]);
				}
			}

		}
	}
	/*
	 * Below is the case of a coupled mesh where the same element is used for analysis as well design.
	 * For clarity purpose two different meshes are used. However, a one-one correlation exists between
	 * elements of the 2 meshes.
	 */
	else if (mesh.coupling == true){
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
		}
		else{
			design_vector.clear();
			for (unsigned int i = 0; i < cell_info_vector.size(); ++i){
				for (unsigned int j = 0; j < cell_info_vector[i].design_points.no_points; j++){
					design_vector.push_back(cell_info_vector[i].design_points.rho[j]);
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

	if (mesh.coupling == false && mesh.adaptivityType == "movingdesignpoints"){
		unsigned int no_cells = mesh.triangulation->n_active_cells();	// No of design points in the domain
		//This case assumes only one design points exists with respect to every finite element in the start of MTO
		unsigned int design_count = no_cells * mesh.design_var_per_point();	//Total number of design variables for optimization

		lb.resize(design_count);
		ub.resize(design_count);	//setting the sizes of the bound vectors

		unsigned int k = 0;	//iterating over all design variables

		for (unsigned int i = 0; i < no_cells; ++i){

			lb[k] = 0.25;	//density value
			ub[k] = 1;	//density value
			k++;
			lb[k] = projection.minFact;	//minimum projection factor (elem size)
			ub[k] = projection.maxFact;	//maximum projection factor (elem size)
			k++;

			// updating the location bounds for the dim dimensions
			for (unsigned int j = 0; j < dim; ++j){
				lb[k] = mesh.coordinates[j][0];
				ub[k] = mesh.coordinates[j][1];
				k++;
			}
		}

	}
	else if (mesh.coupling == true){

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
	}

}

//------------------------------------------------------------------------------------------------------
template <int dim>
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
}

template <int dim>
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
}


//---------------------------------------------------------------------------------------------------------
/*
 * This implementation of get_dxPhys_dx is for movingDesignPoints.
 * For every design point, it returns dim+2 values, each w.r.t 1 design variable
 * design variables: rho, rmin, x_cor, y_cor
 */
template <int dim>
void DensityField<dim>::get_dxPhys_dx(
		std::vector<double> &dxPhys_dx,
		CellInfo &cell_info,
		unsigned int qpoint,
		Point<dim> qX,
		CellInfo &density_cell_info,
		unsigned int density_cell_itr2){

	//searching the design point in the neighbor list
	unsigned int design_index;
	for (unsigned int i = 0; i < cell_info.neighbour_cells[qpoint].size(); ++i){
		design_index = cell_info.neighbour_cells[qpoint][i];
		if (design_index != density_cell_itr2)	continue;

		//When the design index matches, the derivatives are calculated below
		double dxPhys_dH = (density_cell_info.density[0] - cell_info.density[qpoint])/cell_info.sum_weights[qpoint];

/*
		double gamma = (density_cell_info.projection_fact - 2.0) / 5.0;	//gamma for the current density point
		double distance = cell_info.neighbour_distance[qpoint][i];
		double dH_dgamma = (cell_info.neighbour_weights[qpoint][i] * cell_info.sum_weights[qpoint])/(2 * gamma) * ((distance * distance) / (gamma * cell_length * cell_length) - 1.0);
		double dgamma_dprojFact = 1.0 / 5.0;
		double dH_dD = (cell_info.neighbour_weights[qpoint][i] * cell_info.sum_weights[qpoint]) * (-distance / (gamma * cell_length * cell_length));	// D is the distance between the gauss point and the design point
		dxPhys_dx[0] = cell_info.neighbour_weights[qpoint][i];	//w.r.t design density
		dxPhys_dx[1] = dxPhys_dH * dH_dgamma * dgamma_dprojFact;	//w.r.t. projection factor
*/


/*
		//For linear projection
		double distance = cell_info.neighbour_distance[qpoint][i];
		double dH_dproj = cell_length;
		double dH_dD = -1.0;	// D is the distance between the gauss point and the design point
		dxPhys_dx[0] = cell_info.neighbour_weights[qpoint][i];	//w.r.t design density
		dxPhys_dx[1] = dxPhys_dH * dH_dproj;	//w.r.t. projection factor
*/


		//For cubic spline

		double rmin = density_cell_info.projection_fact * cell_length;		//rmin for the current design point
		double temp1 = cell_info.neighbour_distance[qpoint][i];
		double r = temp1/rmin;
		double dH_dr = 0.0;
		if (r <= 0.5)	dH_dr = -8.0*r + 12.0*r*r;
		else if(r > 0.5 && r <= 1.0)		dH_dr = -4.0 + 8.0*r - 4.0*r*r;
		double dr_drmin = -temp1 / (rmin * rmin);
		double dr_dD = 1.0/rmin;
		double dH_dD = dH_dr * dr_dD;
		double drmin_dproj = cell_length;
		dxPhys_dx[0] = cell_info.neighbour_weights[qpoint][i];	//w.r.t design density
		dxPhys_dx[1] = dxPhys_dH * dH_dr * dr_drmin * drmin_dproj;	//w.r.t. projection factor

		for (unsigned int k = 0; k < dim; k++){
			double dD_ddimk = -(qX[k] - density_cell_info.pointX[k])/cell_info.neighbour_distance[qpoint][i];
			dxPhys_dx[k+2] = dxPhys_dH * dH_dD * dD_ddimk;	//adding the sens w.r.t dim_k
		}
/*		if (density_cell_info.pointX[0] == 0.1 && density_cell_info.pointX[1] == 0.1){
			std::cout<<"Left bottom : "<<dxPhys_dx[0]<<" "<<dxPhys_dx[1]<<" "<<dxPhys_dx[2]<<" "<<dxPhys_dx[3]<<std::endl;
		}

		if (density_cell_info.pointX[0] == 0.1 && density_cell_info.pointX[1] == 0.9){
			std::cout<<"Left top : "<<dxPhys_dx[0]<<" "<<dxPhys_dx[1]<<" "<<dxPhys_dx[2]<<" "<<dxPhys_dx[3]<<std::endl;
		}*/
	}

}

//---------------------------------------------------------------------------------------------------------
template <int dim>
double DensityField<dim>::get_vol_fraction(
		std::vector<CellInfo> &cell_info_vector
		){
	double volume = 0.0;
	for(unsigned int i = 0; i < cell_info_vector.size(); ++i){

		cell_info_vector[i].cell_density = 0.0;
		if(cell_info_vector[i].density_weights.size() == 0)
			std::cout<<"ERRORRRRRRRR! no quad point here "<<std::endl;
		for(unsigned int qpoint = 0; qpoint < cell_info_vector[i].density_weights.size(); ++qpoint){
			cell_info_vector[i].cell_density += cell_info_vector[i].density_weights[qpoint] * cell_info_vector[i].density[qpoint];
		}
		double area_fraction = cell_info_vector[i].cell_area / max_cell_area;
		volume += (cell_info_vector[i].cell_density * area_fraction);
	}
	volume /= initial_no_cells;
	//std::cout<<"Volume fraction: "<<volume<<std::endl;
	return volume;
}

template <int dim>
void DensityField<dim>::update_density_cell_info_vector(
		std::vector<CellInfo> &density_cell_info_vector,
		const std::vector<double> &design_vector){
	//std::cout<<density_cell_info_vector.size()<<"..........."<<design_vector.size()<<std::endl;

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

	unsigned int k = 0;	//iterator over the design points
	for (unsigned int i = 0; i < cell_info_vector.size(); ++i){
		for (unsigned int j = 0; j < cell_info_vector[i].design_points.no_points; ++j){
			cell_info_vector[i].design_points.rho[j] = design_vector[k];
			k++;
		}
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
		if (mesh.coupling == false && mesh.adaptivityType != "movingdesignpoints"){
			return (cell_info_vector.size());
		}
		else{
			return (density_cell_info_vector.size() * mesh.design_var_per_point());
		}
	}
	else{
		if (mesh.coupling == false && mesh.adaptivityType != "movingdesignpoints"){
			//Calculating the no. of design variables from cell_info_vector
			unsigned int no_design_count = 0;
			for(unsigned int i = 0; i < cell_info_vector.size(); ++i){
				no_design_count += cell_info_vector[i].design_points.no_points;
			}
			return no_design_count;
		}

		else{
			return (density_cell_info_vector.size() * mesh.design_var_per_point());
		}
	}

	return 0;
}
