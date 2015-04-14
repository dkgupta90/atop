/*
 * elasticity.h
 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef ELASTICITY_H_
#define ELASTICITY_H_

#include <vector>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include<deal.II/fe/fe_values.h>

using namespace dealii;
namespace topopt{
	class ElasticTools{
	public:
		void get_lambda_mu(std::vector<double> &E_values,
				double nu,
				std::vector<double> &lambda_values,
				std::vector<double> &mu_values);
		void get_D_plane_stress2D(FullMatrix<double> &D_matrix,
				double nu);
		void get_B_matrix_2D(std::vector<FullMatrix<double> > &B_matrix_vector,
				std::vector<double> &JxW,
				unsigned int quadrature_rule,
				FESystem<2> &fe,
				DoFHandler<2> &dofhandler);
		void get_normalized_matrix(FullMatrix<double> &D_matrix,
				std::vector<FullMatrix<double> > &B_matrix_vector,
				std::vector<double> &JxW,
				std::vector<FullMatrix<double> > &elem_stiffness_array);
		void display_matrix(FullMatrix<double> &mat);
	};

class StoreElasticData{
public:
	//StoreElasticData();
	unsigned int current_quad_rule, running_quad_rule;
	double nu;
	std::vector<std::vector<FullMatrix<double> > > B_matrix_list;
	std::vector<std::vector<double> > JxW;
	FullMatrix<double> D_matrix;
	std::vector<std::vector<FullMatrix<double> > >  elem_stiffness_array;
	void update_elastic_matrices(FESystem<2> &fe,
			DoFHandler<2> &dofhandler);
	void check_linker();
	unsigned int get_quad_index(unsigned int quad_rule);
};

}





#endif /* ELASTICITY_H_ */
