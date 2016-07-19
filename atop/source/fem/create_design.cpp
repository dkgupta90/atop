/*
 *
 *  Created on: Jun 22, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <atop/fem/create_design.h>
#include <atop/fem/fem.h>
#include <deal.II/grid/grid_tools.h>
#include <math.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/DensityValues.h>
#include <atop/fem/output.h>



using namespace atop;

template<int dim>
void CreateDesign<dim>::assemble_design(
		FEM<dim> &obj_fem){

	//Necessary initializations
	this->fem = &obj_fem;

	std::cout<<"Collecting design point values from the FE domain............."<<std::endl;

	//Define the temporary design_info_vector
	design_info_vector.clear();
	design_info_vector.resize(fem->design_triangulation->n_active_cells());

	//Initializing the density vector
	for (unsigned int i = 0; i < design_info_vector.size(); ++i){
		design_info_vector[i].quad_rule = 1;
		design_info_vector[i].density.resize(1);
	}

	//Iterator for the design and analysis meshes
	typename DoFHandler<dim>::active_cell_iterator design_cell;
	typename DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler->begin_active(),
			endc = fem->dof_handler->end();

	//Iterating over the cell_info_vector
	for (unsigned int cell_itr = 0; cell_itr < fem->cell_info_vector->size(); ++cell_itr){
		std::cout<<"Analysis cell : "<<cell_itr<<std::endl;
		//Iterating over all the design points for current element
		for(unsigned int j = 0; j < (*(fem->cell_info_vector))[cell_itr].design_points.no_points; j++){
			//Creating the point
			Point<dim> pointX;
				for (unsigned int k = 0; k < dim; k++){
					pointX(k) = (*(fem->cell_info_vector))[cell_itr].design_points.pointX[j][k];	//in local coordinates
			}

			//Converting to global coordinates
			Point<dim> point2;
			//converting vector to point coordinates
			Point<dim> centroid = cell->center();	//getting the centre for scaling the points
			double side_length = pow(cell->measure(), 1.0/dim);
			for(unsigned int dimi = 0; dimi < dim; ++dimi){
				point2(dimi) = centroid(dimi) +
				pointX(dimi) * (side_length/2.0);
			}
			//Find surrounding design cell

			design_cell = GridTools::find_active_cell_around_point(*(fem->design_handler), point2);

			//Assigning the design value to the respective location in the design_info_vector
			unsigned int cell_index = design_cell->index();
			design_info_vector[cell_index].cell_density = (*(fem->cell_info_vector))[cell_itr].design_points.rho[j];

		}
		++cell;
	}

	std::cout<<"Finding neighbors for design field....";

	//Adjusting the coupling parameter of the mesh
	fem->mesh->coupling = true;
	design_values.create_neighbors(
			design_info_vector,
			*(fem->fe_design),
			*(fem->fe_design),
			*(fem->design_handler),
			*(fem->design_handler),
			*(fem->projection),
			*(fem->mesh));

	std::cout<<"DONE"<<std::endl;

	std::cout<<"Smoothing the design densities.....";
	design_values.smoothing(design_info_vector);
	std::cout<<"DONE"<<std::endl;


	const unsigned int density_per_design_cell = fem->fe_design->dofs_per_cell;	//dofs per design cell
	std::vector<types::global_dof_index> local_design_indices(density_per_design_cell);
	Vector<double> cell_density(density_per_design_cell);
	cells_adjacent_per_node = 0;
	cells_adjacent_per_node.reinit(fem->design_handler->n_dofs());	//for normalizing the nodal density value


	Vector<double> nodal_design_density;	//For creating the decoupled design mesh
	nodal_design_density.reinit(fem->design_handler->n_dofs());	//	decoupled design output


	//Iterators for the design mesh
	design_cell = fem->design_handler->begin_active();
	typename DoFHandler<dim>::active_cell_iterator design_endc = fem->design_handler->end();

	for (; design_cell != design_endc; ++design_cell){

		unsigned int cell_itr  = design_cell->index();

		cell_density = 0;
		design_cell->get_dof_indices(local_design_indices);
		for(unsigned int i = 0; i < density_per_design_cell; ++i){
			unsigned int n_q_points = 1;
			for(unsigned int q_point = 0 ; q_point < n_q_points; ++q_point){
				cell_density(i) += 	design_info_vector[cell_itr].density[0];
			}
			nodal_design_density(local_design_indices[i]) += cell_density(i);
			cells_adjacent_per_node(local_design_indices[i]) += 1;
		}
	}


	//Normalizing the nodal design density vector
	for(unsigned int i = 0; i < nodal_design_density.size(); ++i){
		unsigned int denom = cells_adjacent_per_node(i);
		if (denom <= 0){
			std::cerr<<"ERROR!! Wrong no. of cells adjacent to node found"<<std::endl;
		}
		else{
			nodal_design_density(i) /= denom;
		}
	}

	//Writing the design output
	std::string filename = "design-";
	std::stringstream ss;
	ss<< fem->cycle +1<<"_"<<fem->itr_count+1;
	filename += ss.str();
	filename += ".vtk";
	std::vector<std::string> density_names;
	density_names.clear();
	density_names.push_back("density");
	OutputData<dim> out_soln;
	out_soln.write_fe_solution(filename, *(fem->design_handler),
			nodal_design_density, density_names);

	nodal_design_density = 0;

	//re-adjusting the tuned parameters
	fem->mesh->coupling = false;
	std::cout<<"Design field created and saved.."<<std::endl;
}
