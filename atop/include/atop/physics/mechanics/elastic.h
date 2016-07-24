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
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

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
				unsigned int p_index,
				unsigned int q_index,
				hp::FECollection<2>&,
				hp::QCollection<2>&,
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
		std::vector<unsigned int> current_quadRuleVector, running_quadRuleVector;
		hp::FECollection<2> fe_collection;
		hp::QCollection<2> quadrature_collection;

		std::vector<std::vector<std::vector<FullMatrix<double> > > > B_matrix_list;
		std::vector<std::vector<std::vector<double> > > JxW;
		FullMatrix<double> D_matrix;
		std::vector<std::vector<std::vector<FullMatrix<double> > > > elem_stiffness_array;
		void initialize_quadRuleVectors(std::vector<unsigned int>&,
				std::vector<unsigned int>&);
		void update_elastic_matrices(hp::FECollection<2> &temp_fe_coll,
				hp::QCollection<2> &temp_q_coll,
				hp::DoFHandler<2> &dofhandler);
		void check_linker();
		unsigned int get_quad_index(unsigned int quad_rule);
		unsigned int get_p_index(unsigned int);
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
