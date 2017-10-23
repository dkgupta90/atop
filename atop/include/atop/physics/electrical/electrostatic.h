/*
 *
 *  Created on: Aug 12, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INCLUDE_ATOP_PHYSICS_MECHANICS_ELECTROSTATIC_H_
#define INCLUDE_ATOP_PHYSICS_MECHANICS_ELECTROSTATIC_H_

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

using namespace dealii;

namespace atop{
	template <int dim>
	class ElectrostaticTools{
	public:
		void get_normalized_matrix(unsigned int p_index,
				unsigned int q_index,
				hp::FECollection<dim> &fe_collection,
				hp::QCollection<dim> &quadrature_collection,
				hp::DoFHandler<dim> &dofhandler,
				std::vector<FullMatrix<double> > &elem_stiffness_array);
		void get_face_B_matrices_2D(std::vector<std::vector<FullMatrix<double> > > &B_matrix_vector,
				std::vector<std::vector<double> > &JxW,
						unsigned int p_index,
						unsigned int q_index,
						hp::FECollection<dim> &fe_collection,
						hp::QCollection<dim-1> &face_quadrature_collection,
						hp::DoFHandler<dim> &dofhandler);
		void display_matrix(FullMatrix<double> &mat);

	};


	template <int dim>
	class ElectrostaticData{
	public:
		std::vector<unsigned int> *current_quadRuleVector, *running_quadRuleVector;
		hp::FECollection<dim> *fe_collection;
		hp::QCollection<dim> *quadrature_collection;
		hp::QCollection<dim-1> *face_quadrature_collection;

		std::vector<std::vector<std::vector<FullMatrix<double> > > > B_matrix_list;
		std::vector<std::vector<std::vector<double> > > JxW;
		FullMatrix<double> D_matrix;
		std::vector<std::vector<std::vector<FullMatrix<double> > > > elem_stiffness_array;
		std::vector<std::vector<std::vector<std::vector<FullMatrix<double> > > > > face_B_matrix_list;
		std::vector<std::vector<std::vector<std::vector<double> > > > face_JxW;

		void initialize_quadRuleVectors(std::vector<unsigned int>&,
				std::vector<unsigned int>&);
		void update_normalized_matrices(hp::FECollection<dim> &temp_fe_coll,
				hp::QCollection<dim> &temp_q_coll,
				hp::DoFHandler<dim> &dofhandler);
		void check_linker();
		unsigned int get_quad_index(unsigned int quad_rule);
		unsigned int get_p_index(unsigned int);
	};


	template <int dim>
	class LinearElectrostatic{
	public:
		double E0, Emin;
		LinearElectrostatic();
	};

	//template class LinearElectrostatic<2>;
	//template class ElectrostaticData<2>;
	//template class ElectrostaticTools<2>;
	template class LinearElectrostatic<3>;
	template class ElectrostaticData<3>;
	template class ElectrostaticTools<3>;
}



#endif /* INCLUDE_ATOP_PHYSICS_MECHANICS_ELECTROSTATIC_H_ */
