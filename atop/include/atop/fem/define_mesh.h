/*
 *
 *  Created on: Mar 15, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef DEFINE_MESH_H_
#define DEFINE_MESH_H_

#include <iostream>
#include <vector>
#include <string>
#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <atop/TopologyOptimization/cell_prop.h>

using namespace dealii;

namespace atop{
//using namespace dealii;
	template <int dim>
	class DefineMesh{

	public:
		DefineMesh() = default;
		DefineMesh& operator=(DefineMesh&&) = default;
		//Default constructor

		//Mesh type: coupled means cells of both meshes overlap
		bool coupling;

		//int dim; //Dimensions of the mesh
		unsigned int dofs_per_node;

		//Defines the type of mesh e.g. subdivided hyper rectangle
		std::string meshType;

		//defines the type of element e.g. scalar lagrange element
		std::string elementType;
		std::string density_elementType;

		std::string adaptivityType;
		std::string amrType;

		//defines the element order
		unsigned int initial_el_order, initial_density_el_order;
		unsigned int max_el_order;
		unsigned int max_density_el_order;

		//initial no. of design points per elemetnt
		unsigned int initial_dcount_per_el, max_dcount_per_el;

		//Saves the coordinates for physical dimensions of the mesh
		std::vector<std::vector<double> > coordinates;

		//Provides number of subdivisions along each physical dimension of the mesh for FE and density
		std::vector<unsigned int> subdivisions;
		std::vector<unsigned int> density_subdivisions;

		//Saves the path of the input mesh (if supplied, else empty)
		std::string input_file;

		//Holds all coordinates and magnitudes of all point sources
		std::vector< std::pair< std::vector<double> , std::vector<double> > > point_source_vector;
		std::vector< std::pair< std::vector<double> , std::vector<double> > > point_stiffness_vector;
		std::vector< std::pair< std::vector<double> , std::vector<double> > > point_l_vector;


		//Constructor for getting only number of physical dimensions
		DefineMesh(int);

		//Pointer to the source function supplied by the user
		std::vector<double> (*source_fn)(std::vector<double>);

		//Pointer to the boundary indicator function supplied by the user
		unsigned int (*boundary_indicator)(std::vector<double>);

		//Triangulation objects
		Triangulation<dim>* triangulation;
		Triangulation<dim>* analysis_density_triangulation;

		Triangulation<dim>* design_triangulation;

		//DOF Handler object
		//hp::DoFHandler<dim> dof_handler;

		//Function for creating the mesh
		void createMesh(
				Triangulation<dim> &triangulation,
				Triangulation<dim> &analysis_density_triangulation,
				Triangulation<dim> &design_triangulation);

		void boundary_info();

		unsigned int design_var_per_point();

		void update_outputDesignMesh(Triangulation<dim> &design_triangulation,
				unsigned int);

		void update_analysis_design_connections(hp::DoFHandler<dim>&,
				hp::DoFHandler<dim>&,
				std::vector<CellInfo>&,
				std::vector<CellInfo>&);

	};
	template class DefineMesh<2>;

}

#endif /* DEFINE_MESH_H_ */
