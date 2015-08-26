/*
 *
 *  Created on: Aug 12, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/output.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <string>
#include <fstream>

using namespace atop;
using namespace dealii;


template <int dim>
void OutputData<dim>::write_fe_solution(
		std::string &filename,
		DoFHandler<dim> &dof_handler,
		Vector<double> &solution,
		std::vector<std::string> &solution_names
		){

	std::ofstream output(filename.c_str());
	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);

	data_out.add_data_vector(solution, solution_names);
	data_out.build_patches();
	data_out.write_vtk(output);
}



