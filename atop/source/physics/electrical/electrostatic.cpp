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
									fe_values.shape_grad (i, qpoint) *
									fe_values.JxW(qpoint));
			}
		}

		elem_stiffness_array.push_back(K_matrix);
	}
	//reverting to the actual fe index
	cell->set_active_fe_index(real_p_index);

}

template <int dim>
void ElectrostaticData<dim>::update_normalized_matrices(hp::FECollection<dim> &temp_fe_coll,
		hp::QCollection<dim> &temp_q_coll,
		hp::DoFHandler<dim> &dofhandler){

	std::cout<<"Entered here "<<std::endl;
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
}

template <int dim>
unsigned int ElectrostaticData<dim>::get_quad_index(unsigned int quad_rule){
	return (quad_rule - 1);
}

template <int dim>
unsigned int ElectrostaticData<dim>::get_p_index(unsigned int p_order){
	return (p_order - 1);
}


/*void ElasticData::check_linker(){
	std::vector<FullMatrix<double> > KEquads;
	KEquads = elem_stiffness_array[0][1];
	std::vector<double> JxWquads = JxW[0][1];
	FullMatrix<double> output(8, 8);
	output = 0;
	for(unsigned int i = 0 ; i < KEquads.size(); ++i){
		output.add(KEquads[i], 1);
	}
	ElasticTools els;
	els.display_matrix(output);
}*/

template <int dim>
void ElectrostaticData<dim>::initialize_quadRuleVectors(std::vector<unsigned int> &temp_current,
		std::vector<unsigned int> &temp_running){
	this->running_quadRuleVector = &temp_running;
	this->current_quadRuleVector = &temp_current;
}





