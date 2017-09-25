 /*
 *  Created on: Sep 6, 2017
 *      Author: Deepak K. Gupta
 *
 */

#include <atop/physics/electrical/electrostatic.h>
#include <vector>
#include<iostream>
#include<deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include<deal.II/hp/fe_values.h>

using namespace atop;
using namespace dealii;

template <int dim>
LinearElectrostatic<dim>::LinearElectrostatic(){

}


template <int dim>
void ElectrostaticTools<dim>::get_normalized_matrix(unsigned int p_index,
		unsigned int q_index,
		hp::FECollection<dim> &fe_collection,
		hp::QCollection<dim> &quadrature_collection,
		hp::DoFHandler<dim> &dofhandler,
		std::vector<FullMatrix<double> > &elem_stiffness_array){

	elem_stiffness_array.clear(); // to save all point stiffness matrices for certain p-value and q-rule
	hp::FEValues<dim> hp_fe_values(fe_collection,
			quadrature_collection,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	typename hp::DoFHandler<dim>::active_cell_iterator cell = dofhandler.begin_active();

	unsigned int real_p_index = cell->active_fe_index();	// for saving back after the computation

	cell->set_active_fe_index(p_index);
	hp_fe_values.reinit(cell, q_index);

	const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

	unsigned int n_q_points = fe_values.n_quadrature_points;

	FullMatrix<double> K_matrix(dofs_per_cell, dofs_per_cell);

	for (unsigned int qpoint = 0; qpoint < n_q_points; ++qpoint){
		K_matrix = 0;

		for (unsigned int i = 0; i < dofs_per_cell; ++i){
			for (unsigned int j = 0; j < dofs_per_cell; ++j){
				K_matrix(i, j) += (fe_values.shape_grad (i, qpoint) *
									fe_values.shape_grad (j, qpoint) *
									fe_values.JxW(qpoint));
			}
		}

		elem_stiffness_array.push_back(K_matrix);
	}
	//reverting to the actual fe index
	cell->set_active_fe_index(real_p_index);

}

//The outer dimension of the vectors below denotes the number of faces
template <int dim>
void ElectrostaticTools<dim>::get_face_B_matrices_2D(std::vector<std::vector<FullMatrix<double> > > &B_matrix_vector,
		std::vector<std::vector<double> > &JxW,
				unsigned int p_index,
				unsigned int q_index,
				hp::FECollection<dim> &fe_collection,
				hp::QCollection<dim-1> &face_quadrature_collection,
				hp::DoFHandler<dim> &dofhandler){

	hp::FEFaceValues<dim> hp_fe_face_values(fe_collection,
			face_quadrature_collection,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	typename hp::DoFHandler<dim>::active_cell_iterator cell = dofhandler.begin_active();

	unsigned int real_p_index = cell->active_fe_index();	//saving back after the computation

	cell->set_active_fe_index(p_index);

	//Iterate over all the faces
	B_matrix_vector.clear();
	B_matrix_vector.resize(GeometryInfo<dim>::faces_per_cell);
	JxW.clear();
	JxW.resize(B_matrix_vector.size());

	std::vector<Point<dim> > support_pts = cell->get_fe().get_unit_support_points();

	for (unsigned int iface = 0; iface < GeometryInfo<dim>::faces_per_cell; ++iface){
		hp_fe_face_values.reinit(cell, iface, q_index);

		const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values();
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		unsigned int n_q_points = fe_face_values.n_quadrature_points;


		FullMatrix<double> B_matrix(dofs_per_cell, dofs_per_cell);

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			B_matrix = 0;
			for (unsigned int k = 0; k < 3; ++k){
				unsigned int k0_itr = 0;
				for(unsigned int i = 0; i < dofs_per_cell; ++i){
					if (k == 2){
						unsigned int t1 = 0, t2;
						unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
						if (comp_i == 0)	t1 = 1;
						//Now we need to find dof with same coordinate as ith dof and opposite comp_i
						for (unsigned int j = 0; j < dofs_per_cell; ++j){
							if (fabs(support_pts[i].distance(support_pts[j])) > 1e-10)
								continue;
							unsigned int comp_j = cell->get_fe().system_to_component_index(j).first;
							if (comp_j == t1){
								t2 = j;
								break;
							}

						}
						B_matrix(k, i) = B_matrix(t1, t2);
					}
					else{
						unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;

							B_matrix(comp_i, k0_itr) = fe_face_values.shape_grad(i, q_point)[comp_i];
							k0_itr++;
					}
				}
			}
			B_matrix_vector[iface].push_back(B_matrix);
			JxW[iface].push_back(fe_face_values.JxW(q_point));
		}
	}
	//reverting to the actual fe index
	cell->set_active_fe_index(real_p_index);
}


template <int dim>
void ElectrostaticTools<dim>::display_matrix(FullMatrix<double> &mat){
	unsigned int cols = mat.n_cols();
	unsigned int rows = mat.n_rows();
	for(unsigned int i = 0; i < rows; ++i){
		for(unsigned int j = 0; j < cols; ++j){
			std::cout<<mat(i, j)<<"  ";
		}
		std::cout<<std::endl;
	}
}

template <int dim>
void ElectrostaticData<dim>::update_normalized_matrices(hp::FECollection<dim> &temp_fe_coll,
		hp::QCollection<dim> &temp_q_coll,
		hp::DoFHandler<dim> &dofhandler){

	this->fe_collection = &temp_fe_coll;
	this->quadrature_collection = &temp_q_coll;

	ElectrostaticTools<dim> electrostatic_tool;

	unsigned int max_p_degree = running_quadRuleVector->size();
	//std::cout<<"Max p degree : "<<max_p_degree<<std::endl;
	elem_stiffness_array.resize(max_p_degree);	//referred as stiffness matrix, even for electrical conduction problems

	for (unsigned int degree = 1; degree <= max_p_degree; ++degree){
		//Updating the sizes based on the new current quad rules
		unsigned int p_index = degree - 1;

		if ((*running_quadRuleVector)[p_index] > (*current_quadRuleVector)[p_index])
			continue;
		elem_stiffness_array[p_index].resize((*current_quadRuleVector)[p_index]);
		for (unsigned int i = (*running_quadRuleVector)[p_index]; i <= (*current_quadRuleVector)[p_index]; ++i){
			unsigned int q_index = i - 1;

			electrostatic_tool.get_normalized_matrix(p_index,
					q_index,
					*fe_collection,
					*quadrature_collection,
					dofhandler,
					elem_stiffness_array[p_index][q_index]
					);
		}
		(*running_quadRuleVector)[p_index] = (*current_quadRuleVector)[p_index] + 1;
	}
	//check_linker();
}


template <int dim>
unsigned int ElectrostaticData<dim>::get_quad_index(unsigned int quad_rule){
	return (quad_rule - 1);
}

template <int dim>
unsigned int ElectrostaticData<dim>::get_p_index(unsigned int p_order){
	return (p_order - 1);
}

template <int dim>
void ElectrostaticData<dim>::check_linker(){
	std::vector<FullMatrix<double> > KEquads;
	KEquads = elem_stiffness_array[0][1];
	std::cout<<"Reached here"<<std::endl;
	FullMatrix<double> output(4, 4);
	output = 0;
	for(unsigned int i = 0 ; i < KEquads.size(); ++i){
		output.add(KEquads[i], 1);
	}
	ElectrostaticTools<dim> elecs;
	elecs.display_matrix(output);
}

template <int dim>
void ElectrostaticData<dim>::initialize_quadRuleVectors(std::vector<unsigned int> &temp_current,
		std::vector<unsigned int> &temp_running){
	this->running_quadRuleVector = &temp_running;
	this->current_quadRuleVector = &temp_current;
}





