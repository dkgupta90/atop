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
DefineMesh<dim>::DefineMesh(unsigned int obj_dofs_per_node){
	dofs_per_node = obj_dofs_per_node;
	input_file = "";
	source_fn = NULL;
	point_source_vector.clear();
	subdivisions.clear();
	density_subdivisions.clear();
}


template <int dim>
void DefineMesh<dim>::createMesh(
		Triangulation<dim> &obj_triangulation,
		Triangulation<dim> &obj_analysis_density_triangulation,
		Triangulation<dim> &obj_design_triangulation){

	/* Note that below we assume the moving design points method to be the non-coupled way,
	 * however, there can be non moving design points based methods as sell which are non-coupled.
	 */

	//For the coupled meshes
	if (this->coupling == true){
		this->triangulation = &obj_triangulation;
		this->analysis_density_triangulation = &obj_analysis_density_triangulation;
		this->design_triangulation = &obj_design_triangulation;
		if (meshType == "subdivided_hyper_rectangle"){
			Point<dim> point1, point2;
			if (dim == 2){
				point1 = Point<dim>(coordinates[0][0], coordinates[1][0]);
				point2 = Point<dim>(coordinates[0][1], coordinates[1][1]);
			}
			else if (dim == 3){

				point1 = Point<dim>(coordinates[0][0], coordinates[1][0], coordinates[2][0]);
				point2 = Point<dim>(coordinates[0][1], coordinates[1][1], coordinates[2][1]);
			}
			GridGenerator::subdivided_hyper_rectangle(
					*triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			GridGenerator::subdivided_hyper_rectangle(
					*analysis_density_triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			GridGenerator::subdivided_hyper_rectangle(
					*design_triangulation,
					density_subdivisions,
					point1,
					point2,
					false
					);
			std::cout<<"Active cells in FE mesh: "<<this->triangulation->n_active_cells()<<std::endl;
			std::cout<<"Active cells in density mesh: "<<this->design_triangulation->n_active_cells()<<std::endl;
			boundary_info();
		}
	}
	else{
		this->triangulation = &obj_triangulation;
		this->analysis_density_triangulation = &obj_analysis_density_triangulation;
		this->design_triangulation = &obj_design_triangulation;
		if (meshType == "subdivided_hyper_rectangle"){
			Point<dim> point1, point2;
			if (dim == 2){
				point1 = Point<dim>(coordinates[0][0], coordinates[1][0]);
				point2 = Point<dim>(coordinates[0][1], coordinates[1][1]);
			}
			else if (dim == 3){

				point1 = Point<dim>(coordinates[0][0], coordinates[1][0], coordinates[2][0]);
				point2 = Point<dim>(coordinates[0][1], coordinates[1][1], coordinates[2][1]);
			}
			GridGenerator::subdivided_hyper_rectangle(
					obj_triangulation,
					subdivisions,
					point1,
					point2,
					false
					);
			GridGenerator::subdivided_hyper_rectangle(
					*analysis_density_triangulation,
					subdivisions,
					point1,
					point2,
					false
					);

			//For decoupled mesh, this one can have a different resolution
			GridGenerator::subdivided_hyper_rectangle(
					*design_triangulation,
					density_subdivisions,
					point1,
					point2,
					false
					);
			std::cout<<"Active cells in FE mesh: "<<this->triangulation->n_active_cells()<<std::endl;
			std::cout<<"Active cells in density mesh: "<<this->design_triangulation->n_active_cells()<<std::endl;

			//boundary_info();
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
unsigned int DefineMesh<dim>::design_var_per_point(){
	std::transform(adaptivityType.begin(), adaptivityType.end(),adaptivityType.begin(), ::tolower);
	if (this->coupling == false && this->adaptivityType == "movingdesignpoints"){
		return (dim + 2);
	}
	else{
		return 1;
	}
	return 0;
}

template <int dim>
void DefineMesh<dim>::update_outputDesignMesh(Triangulation<dim> &design_triangulation,
		unsigned int design_per_dim){

	Point<dim> point1, point2;

	if (dim == 2){

		point1 = Point<dim>(coordinates[0][0], coordinates[1][0]);
		point2 = Point<dim>(coordinates[0][1], coordinates[1][1]);
	}
	else if (dim == 3){

		point1 = Point<dim>(coordinates[0][0], coordinates[1][0], coordinates[2][0]);
		point2 = Point<dim>(coordinates[0][1], coordinates[1][1], coordinates[2][1]);
	}
	design_triangulation.clear();
	GridGenerator::subdivided_hyper_rectangle(
			design_triangulation,
			{design_per_dim*subdivisions[0], design_per_dim*subdivisions[1], design_per_dim*subdivisions[2]},
			point1,
			point2,
			false
			);
}
