/*
 *
 *  Created on: Aug 12, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef OUTPUTDATA_H_
#define OUTPUTDATA_H_

#include <deal.II/numerics/data_out.h>
#include <string>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/hp/dof_handler.h>


using namespace dealii;

namespace atop{
	template <int dim>
	class OutputData{
	public:
		void write_fe_solution(
				std::string &filename,
				hp::DoFHandler<dim> &dof_handler,
				Vector<double> &solution,
				std::vector<std::string> &solution_names);
		void write_design(
				std::string &,
				std::vector<double> &,
				unsigned int);
		void write_design(
				std::string &,
				hp::DoFHandler<dim> &,
				std::vector<CellInfo> &);

		/*
		 * This function is used to write the filtered density field in a text format.
		 * This function is implemented since right now density field cannot be read from
		 * vtk files for post-processing, post-evaluation.
		 * Currently, the density is iteratively read from design_handler field and stored.
		 * The first line in the file denotes the number of cells
		 * The later lines denote the density values
		 * Each line has only one value corresponding to constant density of the corresponding
		 * cell in the design_triangulation
		 */
		void write_density(std::vector<CellInfo> &design_info_vector,
				unsigned int cycle,
				unsigned int itr_count);

		/*
		 * This function reads density values from saved file
		 * The number of cells in the analysis mesh should be same as the number
		 * of density entries in the file.
		 * Currently, it only reads for cases, where the density is constant inside the analysis element
		 * This is the case with which the dp-adaptive and FCM cases will be compared for performance
		 */
		void read_xPhys_from_file(std::vector<CellInfo> &,
				const std::string &);
	};

	template class OutputData<2>;
}



#endif /* OUTPUTDATA_H_ */
