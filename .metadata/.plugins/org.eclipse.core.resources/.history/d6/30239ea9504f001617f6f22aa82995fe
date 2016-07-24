/*
 *
 *  Created on: Apr 2, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef ELASTIC_H_
#define ELASTIC_H_

#include <string>
#include <vector>
#include <iostream>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/dof_handler.h>

using namespace dealii;

namespace atop{


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
				hp::DoFHandler<2> &dofhandler);
		void get_normalized_matrix(FullMatrix<double> &D_matrix,
				std::vector<FullMatrix<double> > &B_matrix_vector,
				std::vector<double> &JxW,
				std::vector<FullMatrix<double> > &elem_stiffness_array);
		void display_matrix(FullMatrix<double> &mat);
	};

	class ElasticData{
	public:
		double nu; //Poisson coefficient
		unsigned int current_quad_rule, running_quad_rule;
		std::vector<std::vector<FullMatrix<double> > > B_matrix_list;
		std::vector<std::vector<double> > JxW;
		FullMatrix<double> D_matrix;
		std::vector<std::vector<FullMatrix<double> > >  elem_stiffness_array;
		void update_elastic_matrices(FESystem<2> &fe,
				hp::DoFHandler<2> &dofhandler);
		void check_linker();
		unsigned int get_quad_index(unsigned int quad_rule);
	};

	template<int dim>
	class LinearElastic{
	public:
		double E, poisson; //Define properties of the material
		std::string planarType;  //For 2D problems
		ElasticData obj_elas_data;
		LinearElastic();

	};

	template class LinearElastic<2>;

}


#endif /* ELASTIC_H_ */
