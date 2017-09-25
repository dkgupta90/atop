/*
 *
 *  Created on: Aug 12, 2017
 *      Author: Deepak K. Gupta
 *  
 */

/*
 * This examples aims at a 2D as well as 3D implementation of electrical conduction problems
 * A linear case will be considered for simplicity.
 * However, in future it is aimed to extend it to the nonlinear solar cell problem
 * This is a scalar problem with only one dof per node/support point
 */

#include <atop/fem/define_mesh.h>
#include <atop/physics/electrical/electrostatic.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/TopologyOptimization/optimizedesign.h>
#include <atop/fem/fem.h>
#include <string>
#include <vector>
#include <math.h>
#include <ctime>

using namespace atop;

//Source term for current source distribution in the entire domain
std::vector<double> source_function(std::vector<double> X){
	//This function returns the value of the source function in the whole domain
	if (X.size() == 2){
		std::vector<double> output_source = {310.0};	//photoillumination current
		return output_source;
	}
	else if (X.size() == 3){
		std::vector<double> output_source = {0.0, 0.0, 0.0};
		return output_source;
	}
}

unsigned int get_boundary_indicator(std::vector<double> X){
	//This function defines the boundary indicators

	if (fabs(X[0] - 0) < 1e-12 && fabs(X[1] - 0) < 0.001)
		return 42;
	else
		return 9999;

}

unsigned int get_boundary_indicator_dist(std::vector<double> X){
	//This function defines the boundary indicators

	if (X.size() == 2){	// for two-dimensions
		if (fabs(X[1] - 1) < 1e-12)
			return 62;	//for adding distributed electrical load
		else if (fabs(X[0] - 0) < 1e-12 && fabs(X[1] - 1.0))	//fixed Dirichlet b.c.
			return 42;
		else
			return 9999;
	}
	else if (X.size() == 3){ // for three-dimensional problems

	}


}

int main(){

	deallog.depth_console (2);
	//Define the mesh
	DefineMesh<2> mesh(2);
	mesh.coupling = false;
	mesh.elementType = "FE_Q";
	mesh.density_elementType = "FE_DGQ";
	mesh.adaptivityType = "adaptive_grayness";
	mesh.amrType = "dp-refinement";


	Projection filter("density_filter",
			"dp-refinement", 0.0001875, 1.0);

	//Define the penalization scheme
	Penalize penal("SIMP");
	penal.factmin = 1e-3;
	penal.penal_power = 3.0;

	//Define the physics of the problem
	LinearElectrostatic<2> material1;
	material1.E0 = 100;	//electrode conductivity
	material1.Emin = 0.02;	//TCo conductivity as in Gupta et al 2015, SMO

	//Define the optimization parameters
	Optimizedesign<2> opt(mesh, penal, filter, "OC", 5);
	opt.problem_name = "electrical_conduction";
	//opt.problem_name = "compliant_mechanism";
	opt.is_problem_self_adjoint = false;
	opt.problemType(material1);
	opt.volfrac = 0.3; //Maximum permissible volume fraction

	//Initializing the compulsory variables
	mesh.point_stiffness_vector.clear();
	mesh.point_l_vector.clear();

	//Parameters for defining the test cases for dp-refinement
	//std::string test_problem = "compliant_mechanism2D";
	std::string test_problem = "elec_cond2D";
	unsigned int dim = 2;

	if (dim == 2){
		if (test_problem == "elec_cond2D"){
			mesh.coordinates = {{0, 0.015}, {0, 0.0075}};
			mesh.subdivisions = {80, 40};
			mesh.meshType = "subdivided_hyper_rectangle";

			mesh.initial_el_order = 2;
			mesh.initial_density_el_order = 1;
			mesh.max_el_order = 7;
			mesh.max_density_el_order = 1;
			mesh.initial_dcount_per_el = 16;
			mesh.max_dcount_per_el = 81;
			unsigned int d_per_line = round(sqrt(mesh.initial_dcount_per_el));
			mesh.density_subdivisions = {d_per_line*mesh.subdivisions[0], d_per_line*mesh.subdivisions[1]};


			mesh.source_fn = source_function;

			//Define loads
			std::string loadType = "pointLoad";
			if (loadType == "pointLoad"){
				mesh.boundary_indicator = get_boundary_indicator;
				//Define point force
				std::vector<double> point = {0.0, 0.0};
				std::vector<double> source = {0, 0.0};
				mesh.point_source_vector.push_back(std::make_pair(point, source)); //make pairs and push
				   //empty dist load
			}
			else if (loadType == "distLoad"){

				mesh.boundary_indicator = get_boundary_indicator_dist;
				mesh.point_source_vector.clear();	//no point load
			}
		}
	}
	else if (dim ==3){

	}

	clock_t begin = clock();
	opt.start_time = double (begin)/CLOCKS_PER_SEC;
	opt.temp1 = false;
	opt.tempfname = "";

	opt.optimize();
	clock_t end = clock();
	double elapsed_secs = double (end - begin)/CLOCKS_PER_SEC;
	std::cout<<"Optimization completed......Computing time : "<<elapsed_secs<<std::endl;
	double objMTO = opt.objective;
	opt.temp1 = true;
	std::string filename = "output_design/density_";
	std::stringstream ss;
	ss<<opt.no_cycles<<"_"<<opt.obj_fem->itr_count+1;
	filename += ss.str();
	filename += ".dat";
	opt.tempfname = filename;
	mesh.initial_el_order = 3;
	unsigned int d_per_line = round(sqrt(opt.final_dcount_per_el));
	std::cout<<d_per_line<<std::endl;
	mesh.subdivisions = {d_per_line * 20, d_per_line * 10};
	filter.radius /= d_per_line;
	mesh.initial_dcount_per_el = 1;
	mesh.density_subdivisions = {mesh.initial_dcount_per_el*mesh.subdivisions[0], mesh.initial_dcount_per_el*mesh.subdivisions[1]};
	opt.no_cycles = 1;
	opt.optimize();
	double objTO = opt.objective;
	std::cout<<"Solution accuracy : "<<objMTO/objTO<<std::endl;
	std::cout<<"SUCCESS.... : "<<std::endl;
}




