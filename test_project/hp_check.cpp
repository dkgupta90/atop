/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */
#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <vector>
using namespace dealii;
class Step3
{
public:
  Step3 ();
  void run ();
private:
  void make_grid ();
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results () const;

  void set_boundary_indicator();

  void add_point_source_to_rhs();

  void get_D_plane_stress2D(FullMatrix<double>&,
  		double);
  void get_point_B_matrix_2D(FullMatrix<double>&,
  		double &JxW,
  		typename hp::DoFHandler<2>::active_cell_iterator&,
  		hp::FEValues<2>&,
  		unsigned int,
  		unsigned int);

  Triangulation<2>     triangulation;
  hp::FECollection<2> fe_collection;
  hp::QCollection<2>     quadrature_collection;
  hp::DoFHandler<2>        dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;
  const unsigned int max_degree;
  ConstraintMatrix hanging_node_constraints;
  std::map<types::global_dof_index, double> boundary_values;

};


Step3::Step3 ()
  :
  dof_handler (triangulation),  max_degree(5)
{
	  {
	    for (unsigned int degree=1; degree<=max_degree; ++degree)
	      {
	        fe_collection.push_back(FESystem<2>(FE_Q<2>(degree), 2));
	        quadrature_collection.push_back (QGauss<2>(degree+1));
	      }
	  }
}
void Step3::make_grid ()
{
	std::vector<unsigned int> repetitions;
	repetitions.clear();
	repetitions.push_back(4);
	repetitions.push_back(2);
	Point<2> point1(0, 0);
	Point<2> point2(2, 1);
  GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions, point1, point2);
  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
}
void Step3::setup_system ()
{

	//Assigning FE for each cell

	unsigned int cell_itr = 0;
	for (typename hp::DoFHandler<2>::active_cell_iterator cell =dof_handler.begin_active();
			cell != dof_handler.end(); ++cell){
		cell->set_active_fe_index(3);	//p-degree-1
		if (cell_itr == 5) 		cell->set_active_fe_index(3);
		cell_itr++;
	}

	boundary_values.clear();
  dof_handler.distribute_dofs (fe_collection);
  hanging_node_constraints.clear();
	DoFTools::make_hanging_node_constraints(dof_handler,
			hanging_node_constraints);

	set_boundary_indicator();	//define the Dirichlet boundary
	//Applying the boundary conditions
	VectorTools::interpolate_boundary_values(dof_handler,
			42,
			ZeroFunction<2>(),
			hanging_node_constraints);

	hanging_node_constraints.close();

  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

	std::cout<<"No. of hanging node constraints : "<<hanging_node_constraints.n_constraints()<<std::endl;


  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, dsp, hanging_node_constraints, false);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit (sparsity_pattern);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}
void Step3::assemble_system ()
{

	//Updating the hp_fe_values
	hp::FEValues<2> hp_fe_values(fe_collection,
				quadrature_collection,
				update_values | update_gradients |
				update_quadrature_points | update_JxW_values);


	unsigned int cell_itr = 0;

  hp::DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
	  std::cout<<"Cell : "<<cell_itr<<std::endl;
		hp_fe_values.reinit(cell);
		const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values();

		  const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;

		  const unsigned int   n_q_points    = (fe_values.get_quadrature_points()).size();
		  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
		  Vector<double>       cell_rhs (dofs_per_cell);
		  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

		  FullMatrix<double> normalized_matrix(dofs_per_cell, dofs_per_cell);
          cell_matrix = 0;
          cell_rhs = 0;

		FullMatrix<double> D_matrix(3, 3);
		FullMatrix<double> B_matrix(3, dofs_per_cell);
		double JxW;
		get_D_plane_stress2D(D_matrix, 0.3);

/*        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int
            component_i = cell->get_fe().system_to_component_index(i).first;
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const unsigned int
                component_j = cell->get_fe().system_to_component_index(j).first;
                for (unsigned int q_point=0; q_point<n_q_points;
                     ++q_point)
                  {
                    cell_matrix(i,j)
                    +=
                      (
                        (fe_values.shape_grad(i,q_point)[component_i] *
                         fe_values.shape_grad(j,q_point)[component_j] *
                         0.053)
                        +
                        (fe_values.shape_grad(i,q_point)[component_j] *
                         fe_values.shape_grad(j,q_point)[component_i] *
                         0.035)
                        +
                        ((component_i == component_j) ?
                         (fe_values.shape_grad(i,q_point) *
                          fe_values.shape_grad(j,q_point) *
                          0.035)  :
                         0)
                      )
                      *
                      fe_values.JxW(q_point);
                  }
              }
          }*/

		for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){

			normalized_matrix = 0;

			get_point_B_matrix_2D(B_matrix,  JxW, cell, hp_fe_values, cell->active_fe_index(), q_point);

			//normalized_matrix = elastic_data.elem_stiffness_array[p_index][q_index][q_point];
			normalized_matrix.triple_product(D_matrix,
					B_matrix,
					B_matrix,
					true,
					false,
					JxW);

			cell_matrix.add(0.091125,
					normalized_matrix);
		}

/*		//Calculating cell_rhs
		for(unsigned int i = 0; i < dofs_per_cell; ++i){
			const unsigned int component_i = cell->get_fe().system_to_component_index(i).first;
			for(unsigned int q_point = 0; q_point < n_q_points; ++q_point){
				double body_load = 0;
				if (component_i == 1 && cell_itr%4 == 3) body_load = 0.5;
				cell_rhs(i) += fe_values.shape_value(i, q_point) *
						body_load * //rhs_values[q_point](component_i) *
						fe_values.JxW(q_point);
			}
		}*/


	  add_point_source_to_rhs();
      cell->get_dof_indices (local_dof_indices);
		 hanging_node_constraints.distribute_local_to_global (cell_matrix,
		                                          cell_rhs,
		                                          local_dof_indices,
		                                          system_matrix,
		                                          system_rhs);
		 cell_itr++;
    }
  	  //system_rhs.print(std::cout);
/*
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);*/

}
void Step3::solve ()
{
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              solver (solver_control);
  solver.solve (system_matrix, solution, system_rhs,
                PreconditionIdentity());

	hanging_node_constraints.distribute(solution);

}
void Step3::output_results () const
{

    DataOut<2,hp::DoFHandler<2> > data_out;
    data_out.attach_dof_handler (dof_handler);
	std::vector<std::string> solution_names;
	solution_names.clear();
	solution_names.push_back("x_displacement");
	solution_names.push_back("y_displacement");
    data_out.add_data_vector (solution, solution_names);
    solution_names.clear();
	solution_names.push_back("x_load");
	solution_names.push_back("y_load");
    data_out.add_data_vector (system_rhs, solution_names);

    data_out.build_patches ();
    const std::string filename = "solution.vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);

}

void Step3::run ()
{
  make_grid ();
  std::cout<<"Grid made "<<std::endl;
  setup_system ();
  std::cout<<"system set up"<<std::endl;
  assemble_system ();
  std::cout<<"System assembled"<<std::endl;
  solve ();
  std::cout<<"System solved"<<std::endl;
  output_results ();
  std::cout<<"Results stored"<<std::endl;
}
int main ()
{
  deallog.depth_console (2);
  Step3 laplace_problem;
  laplace_problem.run ();
  return 0;
}


//For defining the boundary indicator
void Step3::set_boundary_indicator(){
	for (typename Triangulation<2>::active_cell_iterator
			cell = triangulation.begin();
			cell != triangulation.end();
			++cell ){
		for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f){
			Point<2> point1 = cell->face(f)->center();

			if (fabs(point1(0) - 0) < 1e-12)
				cell->face(f)->set_boundary_id(42);
		}
	}
}

//Adding point load
void Step3::add_point_source_to_rhs(){
	deallog.depth_console (2);

	typename hp::DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
			endc = dof_handler.end();
	for (; cell != endc; ++cell){
		const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
		std::vector<types::global_dof_index> global_dof_indices(dofs_per_cell);
		cell->get_dof_indices(global_dof_indices);
		for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v){
			if (fabs(cell->vertex(v)[0] - 2.0) < 1e-12 && fabs(cell->vertex(v)[1] - 0.5) < 1e-12){
				system_rhs(cell->vertex_dof_index(v-1, 1)) = 0;
				system_rhs(cell->vertex_dof_index(v-1, 2)) = 1;
				//system_rhs[global_dof_indices[2*v]] = 0;
				//system_rhs[global_dof_indices[2*v + 1]] = 1;
			}
		}
	}

/*	hp::MappingCollection<2> mapping;
	mapping.push_back(MappingQ<2>(1));
	VectorTools::create_point_source_vector(mapping, dof_handler, {2.0, 0.5}, {0, 1}, system_rhs);*/
}

//Additional functions related to the physics of the problem
void Step3::get_point_B_matrix_2D(FullMatrix<double> &B_matrix,
		double &JxW,
		typename hp::DoFHandler<2>::active_cell_iterator &cell,
		hp::FEValues<2> &hp_fe_values,
		unsigned int q_index,
		unsigned int q_point){

	hp_fe_values.reinit(cell, q_index);
	const FEValues<2> &fe_values = hp_fe_values.get_present_fe_values();
	//std::cout<<"No. of quad points : "<<(fe_values.get_quadrature_points()).size()<<std::endl;


	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

	std::vector<Point<2> > support_pts = cell->get_fe().get_unit_support_points();
	//unsigned int n_q_points = fe_values.n_quadrature_points;
	for(unsigned int i = 0; i < support_pts.size(); ++i)
		std::cout<<support_pts[i](0)<<"    "<<support_pts[i](01)<<std::endl;
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

//constitutive matrix
void Step3::get_D_plane_stress2D(FullMatrix<double> &D_matrix,
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
