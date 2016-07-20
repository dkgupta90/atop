/*
 *
 *  Created on: Aug 12, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/output.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <string>
#include <iostream>
#include <fstream>
#include <deal.II/hp/dof_handler.h>

using namespace atop;
using namespace dealii;


template <int dim>
void OutputData<dim>::write_fe_solution(
		std::string &filename,
		hp::DoFHandler<dim> &dof_handler,
		Vector<double> &solution,
		std::vector<std::string> &solution_names
		){

	std::ofstream output(("output/" + filename).c_str());
	DataOut<dim, hp::DoFHandler<dim> > data_out;
	data_out.attach_dof_handler(dof_handler);

	data_out.add_data_vector(solution, solution_names);
	data_out.build_patches();
	data_out.write_vtk(output);
}

template <int dim>
void OutputData<dim>::write_design(
		std::string &filename,
		std::vector<double> &design_vector,
		unsigned int design_var_per_point){

		std::ofstream wfile;
		wfile.open("output_design/" + filename, std::ios::out);
		//Writing the design data
		for (unsigned int i = 0; i < design_vector.size(); ){
			for (unsigned int j = 0; j < design_var_per_point; j++){
				wfile<<design_vector[i]<<"\t";
				i++;
			}
			wfile<<"\n";
		}
		wfile.close();
}

template <int dim>
void OutputData<dim>::write_design(
				std::string &filename,
				hp::DoFHandler<dim> &dof_handler,
				std::vector<CellInfo> &cell_info_vector){
	std::ofstream wfile;
	wfile.open("output_design/" + filename, std::ios::out);

	//Writing the design data
	unsigned int cell_itr = 0;
	typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
				endc= dof_handler.end();
	for(; cell != endc; ++cell){

		//Iterating over the density points
		for (unsigned int i = 0; i < cell_info_vector[cell_itr].design_points.no_points; ++i){
			double density = cell_info_vector[cell_itr].design_points.rho[i];
			wfile<<density<<"\t";
			//converting vector to point coordinates
			Point<dim> centroid = cell->center();	//getting the centre for scaling the points
			double side_length = pow(cell->measure(), 1.0/dim);
			for(unsigned int dimi = 0; dimi < dim; ++dimi){
				double pointi = centroid(dimi) +
				(cell_info_vector[cell_itr].design_points.pointX[i][dimi]) * (side_length/2.0);
				wfile<<pointi<<"\t";
			}
			wfile<<"\n";
		}
		++cell_itr;
	}
}

