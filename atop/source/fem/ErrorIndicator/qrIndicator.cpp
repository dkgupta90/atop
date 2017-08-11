/*
 *
 *  Created on: Aug 7, 2017

 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/ErrorIndicator/qrIndicator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/math_tools/algebra/MatrixVector.h>



using namespace dealii;
using namespace atop;

template <int dim>
QRIndicator<dim>::QRIndicator(
		FEM<dim> &obj_fem,
		std::vector<double> &qr_accuracy,
		double tol_accuracy,
		std::vector<unsigned int> &qr_p_value,
		std::vector<CellInfo> &cell_info_vector
		){
	this->cell_info_vector = &cell_info_vector;
	this->accuracy_vector = &qr_accuracy;
	this->proposed_p_values = &qr_p_value;
	this->fem = &obj_fem;
	this->tol_accuracy = tol_accuracy;
}


/* This function checks the J/Jstar values and proposed the p-order to be used for each finite element*/
template <int dim>
void QRIndicator<dim>::estimate(){

/*	// Getting the 1-time D-matrix
	ElasticTools elastic_tool;
	FullMatrix<double> D_matrix (3, 3);
	elastic_tool.get_D_plane_stress2D(D_matrix, 0.3);*/


	hp::FEValues<dim> hp_fe_values(fem->fe_collection,
				fem->quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);


	// for storing the output of qr distribution
	typename hp::DoFHandler<dim>::active_cell_iterator analysis_density_cell = fem->analysis_density_handler.begin_active(),
			analysis_density_endc = fem->analysis_density_handler.end();
	hp::FEValues<dim> hp_fe_analysis_density_values(fem->fe_analysis_density_collection,
				fem->quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	cells_adjacent_per_node = 0;
	cells_adjacent_per_node.reinit(fem->analysis_density_handler.n_dofs());	//for normalizing the nodal density value

	Vector<double> nodal_qrValue;	//For creating the decoupled design mesh
	nodal_qrValue.reinit(fem->analysis_density_handler.n_dofs());	//	decoupled design output

	/*
	 * Below, an iteration over all the cells is performed.
	 * For each cell, first the solution and load are computed at all the dofs.
	 * Next, solution check is performed with increased or decreased value of p
	 * A decision is made on the value of p for each cell.
	 */

	typename hp::DoFHandler<dim>::active_cell_iterator cell = fem->dof_handler.begin_active(),
			endc = fem->dof_handler.end();
	unsigned int cell_itr;
	for (; cell != endc; ++cell){
		cell_itr = cell->user_index() - 1;
		unsigned int p_index = fem->elastic_data.get_p_index((*cell_info_vector)[cell_itr].old_shape_fn_order);
		unsigned int q_index = fem->elastic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);
		hp_fe_values.reinit(cell, q_index);
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

		//Get the state solution for the current cell
		Vector<double> u_solution(dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i){
			u_solution(i) = fem->solution(local_dof_indices[i]);
			//std::cout<<u_solution(i)<<std::endl;
		}

		// Get the stiffness matrix for this cell for the current p-value
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
		unsigned int n_q_points = (*cell_info_vector)[cell_itr].n_q_points;
		normalized_matrix = 0;
		cell_matrix = 0;
		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			normalized_matrix = 0;
			normalized_matrix = (*fem).elastic_data.elem_stiffness_array[p_index][q_index][q_point];
			//NaN condition check ----------------------------------------------------------------------------------
			if ((*cell_info_vector)[cell_itr].E_values[q_point] != (*cell_info_vector)[cell_itr].E_values[q_point])
				std::cout<<q_point<<(*cell_info_vector)[cell_itr].E_values[q_point]<<std::endl;
			//------------------------------------------------------------------------------------------------------
			cell_matrix.add((*cell_info_vector)[cell_itr].E_values[q_point],
					normalized_matrix);
		}

		// Computing J value for the current cell
		Vector<double> temp_array(dofs_per_cell);
		temp_array = 0;
		Matrix_Vector matvec;
		matvec.vector_matrix_multiply(
				u_solution,
				cell_matrix,
				temp_array,
				dofs_per_cell,
				dofs_per_cell);
		double Jvalue = matvec.vector_vector_inner_product(
				temp_array,
				u_solution);
		/*
		 * Here, we start with looking for higher values of p and then lower values
		 */
		unsigned int current_p_order = (*cell_info_vector)[cell_itr].old_shape_fn_order;
		unsigned int max_p = current_p_order + 3;

		// Getting the solution for lower values of p
		unsigned int new_p = current_p_order + 4;
		double sum_JJstar = 0.0;
		while (new_p >= 2){
			// Get the Jvalue for this value of p
			double Jstar = get_Jvalue(cell, u_solution, new_p);
			sum_JJstar += (Jvalue/Jstar);
			std::cout<<cell_itr<<"  "<<new_p<<"  "<<Jvalue/Jstar<<std::endl;
			new_p-=1;
		}

		sum_JJstar /= 5;

		// Adding to the qrValue vector
		const unsigned int density_per_design_cell = analysis_density_cell->get_fe().dofs_per_cell;

		std::vector<types::global_dof_index> local_design_indices(density_per_design_cell);
		Vector<double> cell_qr(density_per_design_cell);

		cell_qr = 0;
		analysis_density_cell->get_dof_indices(local_design_indices);

		for(unsigned int i = 0; i < density_per_design_cell; ++i){
			unsigned int n_q_points = 1;
			for(unsigned int q_point = 0 ; q_point < n_q_points; ++q_point){
				cell_qr(i) += 	sum_JJstar;

			}
			nodal_qrValue(local_design_indices[i]) += cell_qr(i);
			cells_adjacent_per_node(local_design_indices[i]) += 1;

		}
		++analysis_density_cell;
		//std::cout<<q_index<<"   "<<(*cell_info_vector)[cell_itr].old_shape_fn_order<<"   "<<n_q_points<<"  "<<Jvalue<<std::endl;
	}

	//Normalizing the nodal qr vector
	for(unsigned int i = 0; i < nodal_qrValue.size(); ++i){
		unsigned int denom = cells_adjacent_per_node(i);
		if (denom <= 0){
			std::cerr<<"ERROR!! Wrong no. of cells adjacent to node found"<<std::endl;
		}
		else{
			nodal_qrValue(i) /= denom;
			nodal_qrValue(i) /= denom;
			nodal_qrValue(i) /= denom;

		}
	}

	//writing the output qrvalue file
	std::string filename = "qrValue-";
	std::stringstream ss;
	ss<< fem->cycle +1;
	filename += ss.str();
	filename += ".vtk";
	std::vector<std::string> property_names;
	property_names.clear();
	property_names.push_back("p-order");
	OutputData<dim> out_soln;
	out_soln.write_fe_solution(filename, fem->analysis_density_handler,
			nodal_qrValue, property_names);
	std::cout<<"Output QR file created ..."<<std::endl;


}

/*
 * This function takes a p-value as an input and a cell iterator and returns the local strain energy value for this cell
 * It needs to interpolate the displacement field for this p-value from the old_shape_fn_order and construct the stiffness
 * matrix as well for this order of the shape function.
 */
template <int dim>
double QRIndicator<dim>::get_Jvalue(hp::DoFHandler<2>::active_cell_iterator cell,
		Vector<double> &u_solution,
		unsigned int new_p){

	unsigned int cell_itr = cell->user_index() - 1;
	unsigned int p_index = fem->elastic_data.get_p_index(new_p);

	//Create a triangulation of one element
	Triangulation<dim> temp_tria;
	GridGenerator::hyper_cube(temp_tria);
	DoFHandler<dim> temp_dofh(temp_tria);
	FESystem<dim> fe(FE_Q<dim>(new_p), dim);

	temp_dofh.clear();
	temp_dofh.distribute_dofs(fe);
	//Getting the quad rule
	unsigned int qrule = fem->gauss_int.get_quadRule((*cell_info_vector)[cell_itr].pseudo_design_points.no_points,
																		new_p);
	QGauss<dim> quad_formula(qrule);
	FEValues<dim> fe_values (fe, quad_formula,
	                             update_values   | update_gradients |
	                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points = quad_formula.size();
    //std::cout<<dofs_per_cell<<"  "<<n_q_points<<std::endl;

    // construct the interpolated solution for this element
	typename hp::DoFHandler<dim>::active_cell_iterator new_cell = temp_dofh.begin_active();
	fe_values.reinit(new_cell);

	//Get the coordinates for all the support points
	std::vector<Point<dim> > support_pts = fe.get_unit_support_points();
	Vector<double> new_solution(dofs_per_cell);

    std::vector<Vector<double> > temp_solution(support_pts.size(), Vector<double>(dim));



	// Create artifical quad for solution interpolation
	Quadrature<dim> temp_quad_formula(support_pts);

	hp::QCollection<dim> temp_qcollection;
	temp_qcollection.push_back(temp_quad_formula);
	hp::FEValues<dim> hp_fe_values(fem->fe_collection,
				temp_qcollection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);

	hp_fe_values.reinit(cell, 0);
	const FEValues<dim> &old_fe_values = hp_fe_values.get_present_fe_values();

	std::vector<Point<dim> > quad_points = old_fe_values.get_quadrature_points();

/*	//priting out the gauss points used for evaluation
	std::cout<<"Evaluation gauss points : "<<std::endl;
	for (unsigned int i = 0; i < quad_points.size(); ++i){
		std::cout<<quad_points[i](0)<<"   "<<quad_points[i](1)<<std::endl;
	}*/
	old_fe_values.get_function_values((*fem).solution, temp_solution);

	//Iterate over all the support points and check
	for (unsigned int i = 0; i < dofs_per_cell; ++i){
		unsigned int comp_i = fe.system_to_component_index(i).first;
		if (comp_i == 0)
			new_solution(i) = temp_solution[i](0);
		else
			new_solution(i) = temp_solution[i](1);
	}

/*	if (new_p == 4){
		// print the original solution
		std::cout<<"Original solution : "<<std::endl;
		for (unsigned int i = 0; i < u_solution.size(); ++i){
			std::cout<<u_solution(i)<<std::endl;
		}

		// print the new solution
		std::cout<<"New solution : "<<std::endl;
		for (unsigned int i = 0; i < new_solution.size(); ++i){
			std::cout<<new_solution(i)<<std::endl;
		}

		exit(0);
	}*/
	//for (unsigned int i = 0; i < new_solution.size(); ++i)	std::cout<<new_solution(i)<<std::endl;

/*
 * Next, the stiffness matrix needs to be calculated for this element
 * This requires calculating the stifness matrices at all the integration points of new_cell and
 * calculating the E_values at each of these points. The calculation of E_Values requires doing the
 * projection again on the pseudo-design.
 */
	CellInfo new_cell_info;	// to save details related to new cell
	//Find the neighbors of the current cell for smoothing
	(*fem).density_field.find_neighbors(cell,
			new_cell,
			fe_values,
			new_cell_info,
			*(cell_info_vector));
	//exit(0);
	// Apply smoothing
	(*fem).density_field.smoothing(new_cell_info,
			*cell_info_vector);

	// Update the Evalues
	(*fem).penal->update_param((*fem).linear_elastic->E, new_cell_info);

/*	if (new_p == 3 && cell_itr == 177){
		//printing the actual E_values
		std::cout<<"Original E values"<<std::endl;
		for (unsigned int i = 0; i < (*cell_info_vector)[cell_itr].E_values.size(); ++i)
			std::cout<<(*cell_info_vector)[cell_itr].E_values[i]<<std::endl;

		std::cout<<"Current E values"<<std::endl;
		for (unsigned int i = 0; i < new_cell_info.E_values.size(); ++i)
			std::cout<<new_cell_info.E_values[i]<<std::endl;

		exit(0);
	}*/

	// Calculating the objective value
	// Get the stiffness matrix for this cell for the current p-value
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
	normalized_matrix = 0;
	cell_matrix = 0;
	unsigned int q_index = fem->elastic_data.get_quad_index(qrule);
	for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
		normalized_matrix = 0;
		normalized_matrix = (*fem).elastic_data.elem_stiffness_array[p_index][q_index][q_point];
		//NaN condition check ----------------------------------------------------------------------------------
		if (new_cell_info.E_values[q_point] != new_cell_info.E_values[q_point])
			std::cout<<cell_itr<<"  "<<new_p<<"  "<<q_point<<"   "<<n_q_points<<"  "<<new_cell_info.E_values[q_point]<<std::endl;
		//------------------------------------------------------------------------------------------------------
		cell_matrix.add(new_cell_info.E_values[q_point],
				normalized_matrix);
	}

	// Computing J value for the current cell
	Vector<double> temp_array(dofs_per_cell);
	temp_array = 0;
	Matrix_Vector matvec;
	matvec.vector_matrix_multiply(
			new_solution,
			cell_matrix,
			temp_array,
			dofs_per_cell,
			dofs_per_cell);
	double Jstar = matvec.vector_vector_inner_product(
			temp_array,
			new_solution);

	return Jstar;
}
