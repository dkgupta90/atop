/*
 *
 *  Created on: Mar 15, 2015
 *      Author: Deepak K. Gupta
 *  
 */
#include <string>
#include <atop/fem/define_mesh.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>

using namespace atop;
using namespace dealii;


//Constructor for getting only number of physical dimensions
template <int dim>
DefineMesh<dim>::DefineMesh(int obj_dofs_per_node){
	dofs_per_node = obj_dofs_per_node;
	input_file = "";
	source_fn = NULL;
	point_source_vector.clear();
	subdivisions.clear();
	density_subdivisions.clear();
}


template <int dim>
void DefineMesh<dim>::createMesh(
		Triangulation<dim> &triangulation,
		Triangulation<dim> &fe_density_triangulation,
		Triangulation<dim> &density_triangulation){

	/* Note that below we assume the moving design points method to be the non-coupled way,
	 * however, there can be non moving design points based methods as sell which are non-coupled.
	 */

	//For the coupled meshes
	if (this->coupling == true){
		this->triangulation = &triangulation;
		this->fe_density_triangulation = &fe_density_triangulation;
		this->density_triangulation = &density_triangulation;
		if (meshType == "subdivided_hyper_rectangle"){
			Point<dim> point1, point2;
			if (dim == 2){
				point1 = Point<2>(coordinates[0][0], coordinates[1][0]);
				point2 = Point<2>(coordinates[0][1], coordinates[1][1]);
			}
			GridGenerator::subdivided_hyper_rectangle(
					triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			GridGenerator::subdivided_hyper_rectangle(
					fe_density_triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			GridGenerator::subdivided_hyper_rectangle(
					density_triangulation,
					density_subdivisions,
					point1,
					point2,
					false
					);
			std::cout<<"Active cells in FE mesh: "<<this->triangulation->n_active_cells()<<std::endl;
			std::cout<<"Active cells in density mesh: "<<this->density_triangulation->n_active_cells()<<std::endl;
			boundary_info();
		}
	}
	else{
		this->triangulation = &triangulation;
		this->fe_density_triangulation = &fe_density_triangulation;
		this->density_triangulation = &density_triangulation;
		if (meshType == "subdivided_hyper_rectangle"){
			Point<dim> point1, point2;
			if (dim == 2){
				point1 = Point<2>(coordinates[0][0], coordinates[1][0]);
				point2 = Point<2>(coordinates[0][1], coordinates[1][1]);
			}
			GridGenerator::subdivided_hyper_rectangle(
					triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			GridGenerator::subdivided_hyper_rectangle(
					fe_density_triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			std::cout<<"Active cells in FE mesh: "<<this->triangulation->n_active_cells()<<std::endl;
			boundary_info();
		}
	}
}

template <int dim>
void DefineMesh<dim>::boundary_info(){

	typename Triangulation<dim>::active_cell_iterator
	cell = triangulation->begin(),
	endc = triangulation->end();
	for (; cell != endc; ++cell ){
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f){
			Point<dim> midpoint = cell->face(f)->center();
				std::vector<double> X;
				X.clear();
				for(unsigned int i = 0; i <dim; i++){
					X.push_back(midpoint(i));
				}
				unsigned int indicator = boundary_indicator(X);
				cell->face(f)->set_boundary_indicator(indicator);
		}
	}
}

template <int dim>
unsigned int DefineMesh<dim>::set_no_of_design_parameters(){
	std::transform(adaptivityType.begin(), adaptivityType.end(),adaptivityType.begin(), ::tolower);
	if (this->coupling == false && this->adaptivityType == "movingdesignpoints"){
		return ((dim + 2) * triangulation->n_active_cells());
	}
	else{
		return (density_triangulation->n_active_cells());
	}
	return 0;
}

