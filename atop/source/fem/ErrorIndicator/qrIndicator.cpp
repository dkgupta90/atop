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
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/hp/fe_values.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/math_tools/algebra/MatrixVector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>





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
void QRIndicator<dim>::scalar_estimate(std::vector<double> &qr_error){

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
		/*for (unsigned int i = 1; i <= 35; ++i)	++cell;*/
		cell_itr = cell->user_index() - 1;
		unsigned int p_index = fem->electrostatic_data.get_p_index((*cell_info_vector)[cell_itr].old_shape_fn_order);
		unsigned int qrule = (*cell_info_vector)[cell_itr].quad_rule;
		unsigned int no_design_pts = (*cell_info_vector)[cell_itr].old_design_count;
		unsigned int q_index = fem->electrostatic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);
		hp_fe_values.reinit(cell, q_index);
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

		//if (fem->cycle == 1 && cell_itr == 0)	exit(0);
		//Get the state solution for the current cell
		Vector<double> u_solution(dofs_per_cell);	// for the current displacement solution
		Vector<double> f_solution(dofs_per_cell);	// for the current force solution
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i){
			u_solution(i) = fem->solution(local_dof_indices[i]);
		}

		// Get the stiffness matrix for this cell for the current p-value
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
		FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
		unsigned int n_q_points = (*cell_info_vector)[cell_itr].n_q_points;
		normalized_matrix = 0;
		cell_matrix = 0;

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
			normalized_matrix = 0;
			normalized_matrix = (*fem).electrostatic_data.elem_stiffness_array[p_index][q_index][q_point];

			//NaN condition check ----------------------------------------------------------------------------------
			if ((*cell_info_vector)[cell_itr].E_values[q_point] != (*cell_info_vector)[cell_itr].E_values[q_point])
				std::cout<<q_point<<(*cell_info_vector)[cell_itr].E_values[q_point]<<std::endl;
			//------------------------------------------------------------------------------------------------------
			cell_matrix.add((*cell_info_vector)[cell_itr].E_values[q_point],
					normalized_matrix);
		}

		// Computing J value for the current cell
/*		Vector<double> temp_array(dofs_per_cell);
		temp_array = 0;*/
		Matrix_Vector matvec;
		matvec.vector_matrix_multiply(
				u_solution,
				cell_matrix,
				f_solution,
				dofs_per_cell,
				dofs_per_cell);
		Vector<double> const_vector(u_solution.size());
		for (unsigned int i = 0; i < const_vector.size(); ++i){
			const_vector(i) = 1.0;
		}
		double Jvalue = matvec.vector_vector_inner_product(
				const_vector,
				u_solution);
		//std::cout<<"Jvalue : "<<Jvalue<<std::endl;

		/*
		 * Here, we start with looking for higher values of p and then lower values
		 */
		unsigned int current_p_order = (*cell_info_vector)[cell_itr].old_shape_fn_order;
		//unsigned int max_p = current_p_order + 3;

		// Getting the solution for lower values of p
		unsigned int new_p = current_p_order + 1;
		double sum_JJstar = 0.0;
		while (new_p >= (*cell_info_vector)[cell_itr].old_shape_fn_order + 1){
			// Get the Jvalue for this value of p
			double Jstar = get_scalar_Jvalue(cell, u_solution, f_solution, new_p);
			sum_JJstar += (Jvalue/Jstar);
			//std::cout<<cell_itr<<"  "<<new_p<<"  "<<Jvalue<<"  "<<Jstar<<"  "<<Jstar/Jvalue<<std::endl;
			new_p-=1;
		}
		//std::cout<<sum_JJstar<<std::endl;
		//exit(0);
		//sum_JJstar /= 5;
		qr_error[cell_itr] = sum_JJstar;	//assigning the error
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
double QRIndicator<dim>::get_scalar_Jvalue(hp::DoFHandler<2>::active_cell_iterator cell,
		Vector<double> &u_solution,
		Vector<double> &f_solution,
		unsigned int new_p){

	unsigned int cell_itr = cell->user_index() - 1;
	unsigned int p_index = fem->electrostatic_data.get_p_index(new_p);

	//Create a triangulation of one element
	double cell_length = sqrt(cell->measure());
	Triangulation<dim> temp_tria, temp_tria2;
	GridGenerator::hyper_cube(temp_tria, 0, cell_length);
	GridGenerator::hyper_cube(temp_tria2, 0, cell_length);
	DoFHandler<dim> temp_dofh(temp_tria);
	DoFHandler<dim> temp_dofh2(temp_tria2);
	FESystem<dim> fe(FE_Q<dim>(new_p), dim);
	FESystem<dim> fe2(FE_Q<dim>((*cell_info_vector)[cell_itr].old_shape_fn_order), dim);

	//Getting the quad rule
	unsigned int qrule = fem->gauss_int.get_quadRule((*cell_info_vector)[cell_itr].old_design_count,
																		new_p);
	QGauss<dim> quad_formula(qrule);

	// fe_values declared below is for a new op value
	FEValues<dim> fe_values (fe, quad_formula,
	                             update_values   | update_gradients |
	                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points = quad_formula.size();
	Vector<double> new_solution(dofs_per_cell);	// for the displacement field on new support points
	Vector<double> new_f_solution(dofs_per_cell);	// for the force field on new support points
	Vector<double> actual_solution(dofs_per_cell);

	SparseMatrix<double> system_matrix;
	Vector<double> system_rhs;
	// Compute the solution for new_p
	actual_solution = 0;

	temp_dofh.clear();
	temp_dofh.distribute_dofs(fe);
	temp_dofh2.clear();
	temp_dofh2.distribute_dofs(fe2);

    // construct the interpolated solution for this element
	typename DoFHandler<dim>::active_cell_iterator new_cell = temp_dofh.begin_active();
	typename DoFHandler<dim>::active_cell_iterator new_cell2 = temp_dofh2.begin_active();

	fe_values.reinit(new_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
	new_cell->get_dof_indices(local_dof_indices);

	//Get the coordinates for all the support points
	std::vector<Point<dim> > support_pts = fe.get_unit_support_points();
    std::vector<Vector<double> > temp_u_solution(support_pts.size(), Vector<double>(dim));
    std::vector<Vector<double> > temp_f_solution(support_pts.size(), Vector<double>(dim));

	// Create artifical quad for solution interpolation
	Quadrature<dim> quad_formula2(support_pts);
	FEValues<dim> fe_values2(fe2, quad_formula2,
								update_values | update_gradients |
								update_quadrature_points | update_JxW_values);
	fe_values2.reinit(new_cell2);
	std::vector<Point<dim> > quad_points = fe_values2.get_quadrature_points();
	fe_values2.get_function_values(u_solution, temp_u_solution);	//Getting displacement at new support points
	fe_values2.get_function_values(f_solution, temp_f_solution);	//Getting force field at new support points

	//Iterate over all the support points and check
	for (unsigned int i = 0; i < dofs_per_cell; ++i){
		unsigned int comp_i = fe.system_to_component_index(i).first;
		if (comp_i == 0){
			new_solution(i) = temp_u_solution[i](0);
			new_f_solution(i) = temp_f_solution[i](0);
		}
		else
		{
			new_solution(i) = temp_u_solution[i](1);
			new_f_solution(i) = temp_f_solution[i](1);
		}

		// Setting internal nodes to 0
		if (fabs(support_pts[i](0) - 0.5)<(0.5-1e-10) && fabs(support_pts[i](1) - 0.5)<(0.5-1e-10)){
			new_f_solution(i) = 0.0;	// case where support point is internal
		}
	}

	ConstraintMatrix constraints;
	constraints.clear();
	constraints.add_line(local_dof_indices[0]);	//since local and global are same at 1-element level
	constraints.add_line(local_dof_indices[1]);
	constraints.set_inhomogeneity(0, new_solution(local_dof_indices[0]));
	constraints.set_inhomogeneity(1, new_solution(local_dof_indices[1]));

	constraints.close();
	DynamicSparsityPattern dsp (temp_dofh.n_dofs());
	DoFTools::make_sparsity_pattern (temp_dofh,
	                                 dsp,
	                                 constraints,
	                                 false);
	SparsityPattern sparsity_pattern;
	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	system_rhs.reinit(temp_dofh.n_dofs());

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
	(*fem).penal->update_param((*fem).linear_electrostatic->E0,
			(*fem).linear_electrostatic->Emin,
			new_cell_info);
/*	for (unsigned int i = 0; i < new_cell_info.E_values.size(); ++i){
		std::cout<<i<<"    "<<new_cell_info.E_values[i]<<std::endl;
	}*/

	// Calculating the objective value
	// Get the stiffness matrix for this cell for the current p-value
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
	normalized_matrix = 0;
	cell_matrix = 0;
	unsigned int q_index = fem->electrostatic_data.get_quad_index(qrule);
	//std::cout<<"No. of quad points : "<<n_q_points<<std::endl;
	for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
		normalized_matrix = 0;
		normalized_matrix = (*fem).electrostatic_data.elem_stiffness_array[p_index][q_index][q_point];
		//NaN condition check ----------------------------------------------------------------------------------
		if (new_cell_info.E_values[q_point] != new_cell_info.E_values[q_point])
			std::cout<<cell_itr<<"  "<<new_p<<"  "<<q_point<<"   "<<n_q_points<<"  "<<new_cell_info.E_values[q_point]<<std::endl;
		//------------------------------------------------------------------------------------------------------
		cell_matrix.add(new_cell_info.E_values[q_point],
				normalized_matrix);
	}
	//cell_matrix.print(std::cout);
	constraints.distribute_local_to_global (cell_matrix,
	                                          new_f_solution,
	                                          local_dof_indices,
	                                          system_matrix,
	                                          system_rhs);

/*	std::cout<<"cell_rhs    system_rhs     actual_solution"<<std::endl;
	for (unsigned int i = 0; i < new_f_solution.size(); ++i)
		std::cout<<new_f_solution(i)<<"   "<<system_rhs(i)<<"   "<<actual_solution(i)<<std::endl;*/
	//std::cout<<"Reached here "<<std::endl;
	//system_matrix.print(std::cout);
	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (actual_solution, system_rhs);
	constraints.distribute(actual_solution);

/*	std::cout<<"Solution for new p: "<<std::endl;
	for (unsigned int i = 0; i < actual_solution.size(); ++i){
		std::cout<<new_solution(i)<<"   "<<actual_solution(i)<<std::endl;
	}*/
	//cell_matrix.print(std::cout);
	// Computing J value for the current cell
	Vector<double> temp_array(dofs_per_cell);
	temp_array = 0;
	Matrix_Vector matvec;
	matvec.vector_matrix_multiply(
			actual_solution,
			cell_matrix,
			temp_array,
			dofs_per_cell,
			dofs_per_cell);
	Vector<double> const_vector(u_solution.size());
	for (unsigned int i = 0; i < const_vector.size(); ++i){
		const_vector(i) = 1.0;
	}
	double Jstar = matvec.vector_vector_inner_product(
			const_vector,
			u_solution);

	return Jstar;
}



/* This function checks the J/Jstar values and proposed the p-order to be used for each finite element*/
template <int dim>
void QRIndicator<dim>::estimate(std::vector<double> &qr_error){

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
		/*for (unsigned int i = 1; i <= 35; ++i)	++cell;*/
		cell_itr = cell->user_index() - 1;
		unsigned int p_index = fem->elastic_data.get_p_index((*cell_info_vector)[cell_itr].old_shape_fn_order);
		unsigned int qrule = (*cell_info_vector)[cell_itr].quad_rule;
		unsigned int no_design_pts = (*cell_info_vector)[cell_itr].old_design_count;
		unsigned int q_index = fem->elastic_data.get_quad_index((*cell_info_vector)[cell_itr].quad_rule);
		hp_fe_values.reinit(cell, q_index);
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

		//if (fem->cycle == 1 && cell_itr == 0)	exit(0);
		//Get the state solution for the current cell
		Vector<double> u_solution(dofs_per_cell);	// for the current displacement solution
		Vector<double> f_solution(dofs_per_cell);	// for the current force solution
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < dofs_per_cell; ++i){
			u_solution(i) = fem->solution(local_dof_indices[i]);
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
/*		Vector<double> temp_array(dofs_per_cell);
		temp_array = 0;*/
		Matrix_Vector matvec;
		matvec.vector_matrix_multiply(
				u_solution,
				cell_matrix,
				f_solution,
				dofs_per_cell,
				dofs_per_cell);
		double Jvalue = matvec.vector_vector_inner_product(
				f_solution,
				u_solution);
		//std::cout<<"Jvalue : "<<Jvalue<<std::endl;

		/*
		 * Here, we start with looking for higher values of p and then lower values
		 */
		unsigned int current_p_order = (*cell_info_vector)[cell_itr].old_shape_fn_order;
		//unsigned int max_p = current_p_order + 3;

		// Getting the solution for lower values of p
		unsigned int new_p = current_p_order + 1;
		double sum_JJstar = 0.0;
		while (new_p >= (*cell_info_vector)[cell_itr].old_shape_fn_order + 1){
			// Get the Jvalue for this value of p
			double Jstar = get_Jvalue(cell, u_solution, f_solution, new_p);
			sum_JJstar += (Jvalue/Jstar);
			//std::cout<<cell_itr<<"  "<<new_p<<"  "<<Jvalue<<"  "<<Jstar<<"  "<<Jstar/Jvalue<<std::endl;
			new_p-=1;
		}
		//std::cout<<sum_JJstar<<std::endl;
		//exit(0);
		//sum_JJstar /= 5;
		qr_error[cell_itr] = sum_JJstar;	//assigning the error
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
		Vector<double> &f_solution,
		unsigned int new_p){

	unsigned int cell_itr = cell->user_index() - 1;
	unsigned int p_index = fem->elastic_data.get_p_index(new_p);

	//Create a triangulation of one element
	double cell_length = sqrt(cell->measure());
	Triangulation<dim> temp_tria, temp_tria2;
	GridGenerator::hyper_cube(temp_tria, 0, cell_length);
	GridGenerator::hyper_cube(temp_tria2, 0, cell_length);
	DoFHandler<dim> temp_dofh(temp_tria);
	DoFHandler<dim> temp_dofh2(temp_tria2);
	FESystem<dim> fe(FE_Q<dim>(new_p), dim);
	FESystem<dim> fe2(FE_Q<dim>((*cell_info_vector)[cell_itr].old_shape_fn_order), dim);

	//Getting the quad rule
	unsigned int qrule = fem->gauss_int.get_quadRule((*cell_info_vector)[cell_itr].old_design_count,
																		new_p);
	QGauss<dim> quad_formula(qrule);

	// fe_values declared below is for a new op value
	FEValues<dim> fe_values (fe, quad_formula,
	                             update_values   | update_gradients |
	                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points = quad_formula.size();
	Vector<double> new_solution(dofs_per_cell);	// for the displacement field on new support points
	Vector<double> new_f_solution(dofs_per_cell);	// for the force field on new support points
	Vector<double> actual_solution(dofs_per_cell);

	SparseMatrix<double> system_matrix;
	Vector<double> system_rhs;
	// Compute the solution for new_p
	actual_solution = 0;

	temp_dofh.clear();
	temp_dofh.distribute_dofs(fe);
	temp_dofh2.clear();
	temp_dofh2.distribute_dofs(fe2);

    // construct the interpolated solution for this element
	typename DoFHandler<dim>::active_cell_iterator new_cell = temp_dofh.begin_active();
	typename DoFHandler<dim>::active_cell_iterator new_cell2 = temp_dofh2.begin_active();

	fe_values.reinit(new_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
	new_cell->get_dof_indices(local_dof_indices);

	//Get the coordinates for all the support points
	std::vector<Point<dim> > support_pts = fe.get_unit_support_points();
    std::vector<Vector<double> > temp_u_solution(support_pts.size(), Vector<double>(dim));
    std::vector<Vector<double> > temp_f_solution(support_pts.size(), Vector<double>(dim));

	// Create artifical quad for solution interpolation
	Quadrature<dim> quad_formula2(support_pts);
	FEValues<dim> fe_values2(fe2, quad_formula2,
								update_values | update_gradients |
								update_quadrature_points | update_JxW_values);
	fe_values2.reinit(new_cell2);
	std::vector<Point<dim> > quad_points = fe_values2.get_quadrature_points();
	fe_values2.get_function_values(u_solution, temp_u_solution);	//Getting displacement at new support points
	fe_values2.get_function_values(f_solution, temp_f_solution);	//Getting force field at new support points

	//Iterate over all the support points and check
	for (unsigned int i = 0; i < dofs_per_cell; ++i){
		unsigned int comp_i = fe.system_to_component_index(i).first;
		if (comp_i == 0){
			new_solution(i) = temp_u_solution[i](0);
			new_f_solution(i) = temp_f_solution[i](0);
		}
		else
		{
			new_solution(i) = temp_u_solution[i](1);
			new_f_solution(i) = temp_f_solution[i](1);
		}

		// Setting internal nodes to 0
		if (fabs(support_pts[i](0) - 0.5)<(0.5-1e-10) && fabs(support_pts[i](1) - 0.5)<(0.5-1e-10)){
			new_f_solution(i) = 0.0;	// case where support point is internal
		}
	}

/*	if (cell_itr == 35){
		// print the original solution
		std::cout<<"Original f solution : "<<std::endl;
		for (unsigned int i = 0; i < f_solution.size(); ++i){
			std::cout<<f_solution(i)<<std::endl;
		}

		// print the new solution
		std::cout<<"New solution : "<<std::endl;
		for (unsigned int i = 0; i < new_f_solution.size(); ++i){
			std::cout<<new_solution(i)<<"   "<<new_f_solution(i)<<std::endl;
		}
	}*/

	ConstraintMatrix constraints;
	constraints.clear();
	constraints.add_line(local_dof_indices[0]);	//since local and global are same at 1-element level
	constraints.add_line(local_dof_indices[1]);
	constraints.add_line(local_dof_indices[2]);
	constraints.add_line(local_dof_indices[3]);
	constraints.set_inhomogeneity(0, new_solution(local_dof_indices[0]));
	constraints.set_inhomogeneity(1, new_solution(local_dof_indices[1]));
	constraints.set_inhomogeneity(2, new_solution(local_dof_indices[2]));
	constraints.set_inhomogeneity(3, new_solution(local_dof_indices[3]));

	constraints.close();
	DynamicSparsityPattern dsp (temp_dofh.n_dofs());
	DoFTools::make_sparsity_pattern (temp_dofh,
	                                 dsp,
	                                 constraints,
	                                 false);
	SparsityPattern sparsity_pattern;
	sparsity_pattern.copy_from(dsp);
	system_matrix.reinit(sparsity_pattern);
	system_rhs.reinit(temp_dofh.n_dofs());

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
/*	for (unsigned int i = 0; i < new_cell_info.E_values.size(); ++i){
		std::cout<<i<<"    "<<new_cell_info.E_values[i]<<std::endl;
	}*/

	// Calculating the objective value
	// Get the stiffness matrix for this cell for the current p-value
	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
	normalized_matrix = 0;
	cell_matrix = 0;
	unsigned int q_index = fem->elastic_data.get_quad_index(qrule);
	//std::cout<<"No. of quad points : "<<n_q_points<<std::endl;
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
	//cell_matrix.print(std::cout);
	constraints.distribute_local_to_global (cell_matrix,
	                                          new_f_solution,
	                                          local_dof_indices,
	                                          system_matrix,
	                                          system_rhs);

/*	std::cout<<"cell_rhs    system_rhs     actual_solution"<<std::endl;
	for (unsigned int i = 0; i < new_f_solution.size(); ++i)
		std::cout<<new_f_solution(i)<<"   "<<system_rhs(i)<<"   "<<actual_solution(i)<<std::endl;*/
	//std::cout<<"Reached here "<<std::endl;
	//system_matrix.print(std::cout);
	SparseDirectUMFPACK  A_direct;
	A_direct.initialize(system_matrix);
	A_direct.vmult (actual_solution, system_rhs);
	constraints.distribute(actual_solution);

/*	std::cout<<"Solution for new p: "<<std::endl;
	for (unsigned int i = 0; i < actual_solution.size(); ++i){
		std::cout<<new_solution(i)<<"   "<<actual_solution(i)<<std::endl;
	}*/
	//cell_matrix.print(std::cout);
	// Computing J value for the current cell
	Vector<double> temp_array(dofs_per_cell);
	temp_array = 0;
	Matrix_Vector matvec;
	matvec.vector_matrix_multiply(
			actual_solution,
			cell_matrix,
			temp_array,
			dofs_per_cell,
			dofs_per_cell);
	double Jstar = matvec.vector_vector_inner_product(
			temp_array,
			actual_solution);

	return Jstar;
}
