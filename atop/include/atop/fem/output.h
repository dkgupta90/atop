/*
 *
 *  Created on: Aug 12, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef OUTPUTDATA_H_
#define OUTPUTDATA_H_

#include <deal.II/numerics/data_out.h>
#include <string.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>


using namespace dealii;

namespace atop{
	template <int dim>
	class OutputData{
	public:
		void write_fe_solution(
				std::string &filename,
				DoFHandler<dim> &dof_handler,
				Vector<double> &solution,
				std::vector<std::string> &solution_names);
	};

	template class OutputData<2>;
}



#endif /* OUTPUTDATA_H_ */
