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

void ElasticTools::get_lambda_mu(std::vector<double> &E_values,
		double nu,
		std::vector<double> &lambda_values,
		std::vector<double> &mu_values){
	unsigned int n_points = E_values.size();
	for(unsigned int i = 0; i < n_points; i++){
		lambda_values[i] = (E_values[i] * nu) / ((1 + nu) * (1 - (2 * nu)));
		//lambda_values[i] = (E_values[i] * nu) / (1 - (nu * nu));
		double K = (1 - (2 * nu))/(2 * nu);
		mu_values[i] = lambda_values[i] * K;
		//std::cout<<lambda_values[i]<<" "<<mu_values[i]<<std::endl;
	}
}

void ElasticTools::display_matrix(FullMatrix<double> &mat){
	unsigned int cols = mat.n_cols();
	unsigned int rows = mat.n_rows();
	for(unsigned int i = 0; i < rows; ++i){
		for(unsigned int j = 0; j < cols; ++j){
			std::cout<<mat(i, j)<<"  ";
		}
		std::cout<<std::endl;
	}
}
void ElasticTools::get_D_plane_stress2D(FullMatrix<double> &D_matrix,
		double nu){
	double k = 1/(1 - (nu*nu));
	D_matrix(0, 0) = k;
	D_matrix(0, 1) = nu*k;
	D_matrix(0, 2) = 0;
	D_matrix(1, 0) = nu*k;
	D_matrix(1, 1) = k;
	D_matrix(1, 2) = 0;
	D_matrix(2, 0) = 0;
	D_matrix(2, 1) = 0;
	D_matrix(2, 2) = ((1 - nu)*k)/2;
	//display_matrix(D_matrix);
}

void ElasticTools::get_B_matrix_2D(std::vector<FullMatrix<double> > &B_matrix_vector,
		std::vector<double> &JxW,
				unsigned int p_index,
				unsigned int q_index,
				hp::FECollection<2> &fe_collection,
				hp::QCollection<2> &quadrature_collection,
				hp::DoFHandler<2> &dofhandler){

	hp::FEValues<2> hp_fe_values(fe_collection,
			quadrature_collection,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	typename hp::DoFHandler<2>::active_cell_iterator cell = dofhandler.begin_active();

	unsigned int real_p_index = cell->active_fe_index();	//saving back after the computation

	cell->set_active_fe_index(p_index);
	hp_fe_values.reinit(cell, q_index);

	const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values();

	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

	unsigned int n_q_points = fe_values.n_quadrature_points;

	FullMatrix<double> B_matrix(3, dofs_per_cell);

	std::vector<Point<2> > support_pts = cell->get_fe().get_unit_support_points();

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

						B_matrix(comp_i, k0_itr) = fe_values.shape_grad(i, q_point)[comp_i];
						k0_itr++;
				}
			}
		}
		B_matrix_vector.push_back(B_matrix);
		JxW.push_back(fe_values.JxW(q_point));
	}

	//reverting to the actual fe index
	cell->set_active_fe_index(real_p_index);
}

void ElasticTools::get_point_B_matrix_2D(FullMatrix<double> &B_matrix,
		double &JxW,
		typename hp::DoFHandler<2>::active_cell_iterator &cell,
		hp::FEValues<2> &hp_fe_values,
		unsigned int q_index,
		unsigned int q_point){

	hp_fe_values.reinit(cell, q_index);
	const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values();

	//std::cout<<"No. of quad points : "<<(fe_values.get_quadrature_points()).size()<<std::endl;

	std::vector<Point<2> > support_pts = cell->get_fe().get_unit_support_points();
	//unsigned int n_q_points = fe_values.n_quadrature_points;
	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	//unsigned int n_q_points = fe_values.n_quadrature_points;
/*	for(unsigned int i = 0; i < support_pts.size(); ++i)
		std::cout<<support_pts[i](0)<<"    "<<support_pts[i](01)<<std::endl;*/
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

					B_matrix(comp_i, k0_itr) = fe_values.shape_grad(i, q_point)[comp_i];
					k0_itr++;
			}
		}
	}
	JxW = fe_values.JxW(q_point);

	//B_matrix.print(std::cout);
}

//The outer dimension of the vectors below denotes the number of faces
void ElasticTools::get_face_B_matrices_2D(std::vector<std::vector<FullMatrix<double> > > &B_matrix_vector,
		std::vector<std::vector<double> > &JxW,
				unsigned int p_index,
				unsigned int q_index,
				hp::FECollection<2> &fe_collection,
				hp::QCollection<1> &face_quadrature_collection,
				hp::DoFHandler<2> &dofhandler){

	hp::FEFaceValues<2> hp_fe_face_values(fe_collection,
			face_quadrature_collection,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	typename hp::DoFHandler<2>::active_cell_iterator cell = dofhandler.begin_active();

	unsigned int real_p_index = cell->active_fe_index();	//saving back after the computation

	cell->set_active_fe_index(p_index);

	//Iterate over all the faces
	B_matrix_vector.clear();
	B_matrix_vector.resize(GeometryInfo<2>::faces_per_cell);
	JxW.clear();
	JxW.resize(B_matrix_vector.size());

	std::vector<Point<2> > support_pts = cell->get_fe().get_unit_support_points();

	for (unsigned int iface = 0; iface < GeometryInfo<2>::faces_per_cell; ++iface){
		hp_fe_face_values.reinit(cell, iface, q_index);

		const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values();
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		unsigned int n_q_points = fe_face_values.n_quadrature_points;
		FullMatrix<double> B_matrix(3, dofs_per_cell);

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

/*		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			B_matrix = 0;
			for (unsigned int k = 0; k < 3; ++k){
				unsigned int k0_itr = 0, k1_itr = 1;
				for(unsigned int i = 0; i < dofs_per_cell; ++i){
					if (k == 2){
						int t1, t2;
						if (i % 2 == 0){
							t1 = 1;
							t2 = i + 1;
						}
						else{
							t1 = 0;
							t2 = i - 1;
						}
						B_matrix(k, i) = B_matrix(t1, t2);
					}
					else{
						unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
						if (comp_i == 0 && k == 0){
							B_matrix(0, k0_itr) = fe_face_values.shape_grad(i, q_point)[k];
							k0_itr += 2;
						}
						else if (comp_i == 1 && k == 1){
							B_matrix(1, k1_itr) = fe_face_values.shape_grad(i, q_point)[k];
							k1_itr += 2;
						}

					}
				}
			}
			B_matrix_vector[iface].push_back(B_matrix);

			JxW[iface].push_back(fe_face_values.JxW(q_point));

		}*/
	}
	//reverting to the actual fe index
	cell->set_active_fe_index(real_p_index);
}


//The outer dimension of the vectors below denotes the number of faces
void ElasticTools::get_face_B_matrix_2D(std::vector<std::vector<FullMatrix<double> > > &B_matrix_vector,
		std::vector<std::vector<double> > &JxW,
				unsigned int q_index,
				hp::DoFHandler<2>::active_cell_iterator &cell,
				hp::FECollection<2> &fe_collection,
				hp::QCollection<1> &face_quadrature_collection){

	hp::FEFaceValues<2> hp_fe_face_values(fe_collection,
			face_quadrature_collection,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);


	//Iterate over all the faces
	B_matrix_vector.clear();
	B_matrix_vector.resize(GeometryInfo<2>::faces_per_cell);
	JxW.clear();
	JxW.resize(B_matrix_vector.size());

	for (unsigned int iface = 0; iface < GeometryInfo<2>::faces_per_cell; ++iface){
		hp_fe_face_values.reinit(cell, iface, q_index);

		const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values();
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		unsigned int n_q_points = fe_face_values.n_quadrature_points;
		FullMatrix<double> B_matrix(3, dofs_per_cell);

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			B_matrix = 0;
			for (unsigned int k = 0; k < 3; ++k){
				unsigned int k0_itr = 0, k1_itr = 1;
				for(unsigned int i = 0; i < dofs_per_cell; ++i){
					if (k == 2){
						int t1, t2;
						if (i % 2 == 0){
							t1 = 1;
							t2 = i + 1;
						}
						else{
							t1 = 0;
							t2 = i - 1;
 						}
						B_matrix(k, i) = B_matrix(t1, t2);
					}
					else{
						unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
						if (comp_i == 0 && k == 0){
							B_matrix(0, k0_itr) = fe_face_values.shape_grad(i, q_point)[k];
							k0_itr += 2;
						}
						else if (comp_i == 1 && k == 1){
							B_matrix(1, k1_itr) = fe_face_values.shape_grad(i, q_point)[k];
							k1_itr += 2;
						}

					}
				}
			}
			B_matrix_vector[iface].push_back(B_matrix);

			JxW[iface].push_back(fe_face_values.JxW(q_point));

		}


	}
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


void ElasticData::update_face_B_matrices(hp::FECollection<2> &temp_fe_coll,
		hp::QCollection<1> &temp_face_q_coll,
		hp::DoFHandler<2> &dofhandler){

	this->fe_collection = &temp_fe_coll;
	this->face_quadrature_collection = &temp_face_q_coll;

	ElasticTools elastic_tool;

	unsigned int max_p_degree = running_quadRuleVector->size();
	//std::cout<<"Max p degree : "<<max_p_degree<<std::endl;
	face_B_matrix_list.resize(max_p_degree);
	face_JxW.resize(max_p_degree);

	for (unsigned int degree = 1; degree <= max_p_degree; ++degree){
		//Updating the sizes based on the new current quad rules
		unsigned int p_index = degree - 1;

		unsigned int new_q_len = B_matrix_list[p_index].size();

		if (new_q_len <= face_B_matrix_list[p_index].size())
			continue;

		unsigned int old_q_len = face_B_matrix_list[p_index].size();

		face_B_matrix_list[p_index].resize(new_q_len);
		face_JxW[p_index].resize(new_q_len);

		for (unsigned int i = old_q_len + 1; i <= new_q_len; ++i){
			unsigned int q_index = i - 1;
			face_B_matrix_list[p_index][q_index].clear();
			face_JxW[p_index][q_index].clear();
			elastic_tool.get_face_B_matrices_2D(face_B_matrix_list[p_index][q_index], face_JxW[p_index][q_index],
					p_index, q_index,
					*fe_collection, *face_quadrature_collection, dofhandler);

		}
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

void ElasticData::initialize_quadRuleVectors(std::vector<unsigned int> &temp_current,
		std::vector<unsigned int> &temp_running){
	this->running_quadRuleVector = &temp_running;
	this->current_quadRuleVector = &temp_current;
}




