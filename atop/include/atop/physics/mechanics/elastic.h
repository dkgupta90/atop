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
#include <deal.II/hp/fe_values.h>
#include <atop/TopologyOptimization/cell_prop.h>


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
		void get_K_matrix_2D(std::vector<FullMatrix<double> > &B_matrix_vector,
				std::vector<double> &JxW,
				FullMatrix<double> &D_matrix,
				unsigned int p_order,
				hp::FECollection<2>&,
				hp::QCollection<2>&,
				hp::DoFHandler<2> &dofhandler,
				hp::DoFHandler<2> &design_handler,
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &design_cell_info_vector
				);
		void get_point_B_matrix_2D(FullMatrix<double> &B_matrix,
				double &JxW,
				typename hp::DoFHandler<2>::active_cell_iterator&,
				hp::FEValues<2>&,
				unsigned int,
				unsigned int);
		void get_face_B_matrices_2D(std::vector<std::vector<FullMatrix<double> > > &B_matrix_vector,
				std::vector<std::vector<double> > &JxW,
						unsigned int p_index,
						unsigned int q_index,
						hp::FECollection<2> &fe_collection,
						hp::QCollection<1> &face_quadrature_collection,
						hp::DoFHandler<2> &dofhandler);

		//Updates for every cell
		void get_face_B_matrix_2D(std::vector<std::vector<FullMatrix<double> > > &B_matrix_vector,
				std::vector<std::vector<double> > &JxW,
						unsigned int q_index,
						hp::DoFHandler<2>::active_cell_iterator &cell,
						hp::FECollection<2> &fe_collection,
						hp::QCollection<1> &face_quadrature_collection);


		void get_normalized_matrix(FullMatrix<double> &D_matrix,
				std::vector<FullMatrix<double> > &B_matrix_vector,
				std::vector<double> &JxW,
				std::vector<FullMatrix<double> > &K_matrix_vector);
		void display_matrix(FullMatrix<double> &mat);
	};

	class ElasticData{
	public:
		double nu; //Poisson coefficient
		std::vector<unsigned int> *current_quadRuleVector, *running_quadRuleVector;
		hp::FECollection<2> *fe_collection;
		hp::QCollection<2> *quadrature_collection;
		hp::QCollection<1> *face_quadrature_collection;

		std::vector<std::vector<std::vector<FullMatrix<double> > > > K_matrix_list;
		std::vector<std::vector<std::vector<double> > > JxW;
		FullMatrix<double> D_matrix;
		std::vector<std::vector<std::vector<FullMatrix<double> > > > elem_stiffness_array;
		std::vector<std::vector<std::vector<std::vector<FullMatrix<double> > > > > face_B_matrix_list;
		std::vector<std::vector<std::vector<std::vector<double> > > > face_JxW;

		void initialize_quadRuleVectors(std::vector<unsigned int>&,
				std::vector<unsigned int>&);
		void update_elastic_matrices(hp::FECollection<2> &temp_fe_coll,
				hp::QCollection<2> &temp_q_coll,
				hp::DoFHandler<2> &dofhandler,
				hp::DoFHandler<2> &design_handler,
				std::vector<CellInfo> &cell_info_vector,
				std::vector<CellInfo> &design_cell_info_vector,
				unsigned int);
		void update_face_B_matrices(hp::FECollection<2> &temp_fe_coll,
				hp::QCollection<1> &temp_face_q_coll,
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
