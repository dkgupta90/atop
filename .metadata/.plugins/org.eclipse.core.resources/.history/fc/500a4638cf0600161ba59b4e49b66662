/*
 *
 *  Created on: Jun 18, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef FEM_H_
#define FEM_H_


#include <iostream>
#include <vector>
#include <string>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <atop/fem/define_mesh.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/fem/define_mesh.h>
#include <atop/TopologyOptimization/DensityValues.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/fem/boundary_values.h>


using namespace dealii;

namespace atop{
	template <int dim>
	class FEM{
	public:
		FEM() = default;
		FEM& operator=(FEM&&) = default;

		// Mesh object
		DefineMesh<dim> *mesh;

		//Design vector for optimization purposes
		std::vector<double> *design_vector;

		//DOFHandler objects
		DoFHandler<dim> *dof_handler, *analysis_density_handler, *design_handler;

		//Triangulation objects
		Triangulation<dim> *triangulation, *analysis_density_triangulation, *design_triangulation;

		//Objects for cell and density cell properties
		std::vector<CellInfo> *cell_info_vector;
		std::vector<CellInfo> *density_cell_info_vector;

		//Object storing parameters related to the physics of the problem
		ElasticData elastic_data;
		LinearElastic<dim> *linear_elastic;
		//FESystem objects
		FESystem<dim> *fe, *fe_analysis_density, *fe_design;

		//Projection object for defining the regularization properties
		Projection *projection;

		//Penalization object for penalizing the density values
		Penalize *penal;

		//DensityField object containing all functions related to density distribution
		DensityField<dim> density_field;
		//Some comment to be added
		ConstraintMatrix hanging_node_constraints;
		SparsityPattern sparsity_pattern;
		SparseMatrix<double> system_matrix;
		Vector<double> solution;
		Vector<double> system_rhs;
		Vector<double> nodal_density;


		unsigned int cycle, itr_count;

		unsigned int current_quad_rule; //For integrating over an element
		unsigned int running_quad_rule; //Used for quad related adaptivity

		double volfrac;	//maximum permissible volume fraction

		//Constructor for initializing the dof_handlers
		FEM(
		Triangulation<dim>&,
		Triangulation<dim>&,
		Triangulation<dim>&,
		DoFHandler<dim>&,
		DoFHandler<dim>&,
		DoFHandler<dim>&,
		std::vector<CellInfo>&,
		std::vector<CellInfo>&,
		DefineMesh<dim>&,
		std::vector<double>&);

		//Solving the FE problem
		void analyze();

		//For the physics of the problem
		void problemType(LinearElastic<dim>&);

		//For the regularizarion type
		void projectionType(Projection&);

		//FOr the penalization
		void penalization(Penalize&);

		//Destructor for the class
		~FEM();

	private:
		Vector<double> cells_adjacent_per_node;
		double projection_radius;
		void clean_trash();
		void boundary_info();
		void setup_system();
		void assemble_system();
		void solve();
		void output_results();
		void reset();
		void initialize_cycle();
		void update_physics();
		void assembly();
		void add_source_to_rhs(
				const std::vector<Point<dim> > &,
				std::vector<Vector<double> > &);
		void add_point_source_to_rhs();

	};
	template class FEM<2>;
}


#endif /* FEM_H_ */
