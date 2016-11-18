/*
 *
 *  Created on: Apr 19, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/define_mesh.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/TopologyOptimization/projection.h>
#include <atop/TopologyOptimization/optimizedesign.h>
#include <atop/fem/fem.h>
#include <string>
#include <vector>
#include <math.h>
#include <ctime>

using namespace atop;


//Empty source term in compliance minimization
std::vector<double> source_function(std::vector<double> X){
	//This function returns the value of the source function in the whole domain
	if (X.size() == 2){
		std::vector<double> output_source = {0.0, 0.0};
		return output_source;
	}
	else if (X.size() == 3){
		std::vector<double> output_source = {0.0, 0.0, 0.0};
		return output_source;
	}
}

unsigned int get_boundary_indicator(std::vector<double> X){
	//This function defines the boundary indicators

	if (fabs(X[0] - 0) < 1e-12)
		return 42;
	else
		return 9999;

}

unsigned int get_boundary_indicator_dist(std::vector<double> X){
	//This function defines the boundary indicators

	if (fabs(X[0] - 0) < 1e-12)	//fixed Dirichlet b.c.
		return 42;
	if (fabs(X[1] - 1) < 1e-12)
		return 62;	//Rolling Dirichlet b.c.
	else
		return 9999;

}

//Boundary indicator function for a compliant force inverter problem
unsigned int get_boundary_indicator_force_inv(std::vector<double> X){
	//This function defines the boundary indicators

	if (fabs(X[0] - 0) < 1e-12 && fabs(X[1] - 0.2) < 1e-12)
		return 42;
	else if (fabs(X[1] - 1) < 1e-12)
		return 62;	//Neumann boundary for distributed load
	else
		return 9999;

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
			"dp-refinement", 0.25, 1.0);

	//Define the penalization scheme
	Penalize penal("SIMP");
	penal.factmin = 1e-9;
	penal.penal_power = 3.0;

	//Define the physics of the problem
	LinearElastic<2> material1;
	material1.E = 1.0;
	material1.poisson = 0.3;
	material1.planarType = "planar_stress";

	//Define the optimization parameters
	Optimizedesign<2> opt(mesh, penal, filter, "OC", 2);
	opt.problem_name = "minimum_compliance";
	opt.problemType(material1);
	opt.volfrac = 0.45; //Maximum permissible volume fraction

	//Parameters for defining the test cases for dp-refinement
	std::string test_problem = "cantilever2D";
	unsigned int dim = 2;

	if (dim == 2){
		if (test_problem == "cantilever2D"){
			mesh.coordinates = {{0, 2}, {0, 1}};
			mesh.subdivisions = {16, 8};
			mesh.meshType = "subdivided_hyper_rectangle";

			mesh.initial_el_order = 2;
			mesh.initial_density_el_order = 1;
			mesh.max_el_order = 11;
			mesh.max_density_el_order = 1;
			mesh.initial_dcount_per_el = 9;
			mesh.max_dcount_per_el = 64;
			unsigned int d_per_line = round(sqrt(mesh.initial_dcount_per_el));
			mesh.density_subdivisions = {d_per_line*mesh.subdivisions[0], d_per_line*mesh.subdivisions[1]};


			mesh.source_fn = source_function;

			//Define loads
			std::string loadType = "distLoad";
			if (loadType == "pointLoad"){
				mesh.boundary_indicator = get_boundary_indicator;
				//Define point force
				std::vector<double> point = {2.0, 0.5};
				std::vector<double> source = {0, 1.0};
				mesh.point_source_vector.push_back(std::make_pair(point, source)); //make pairs and push
				   //empty dist load
			}
			else if (loadType == "distLoad"){

				mesh.boundary_indicator = get_boundary_indicator_dist;
				mesh.point_source_vector.clear();	//no point load
			}
		}

		else if (test_problem == "compliant_mechanism2D"){
			mesh.coordinates = {{0, 2}, {0, 1}};
			mesh.subdivisions = {40, 20};
			mesh.meshType = "subdivided_hyper_rectangle";

			mesh.initial_el_order = 3;
			mesh.initial_density_el_order = 1;
			mesh.max_el_order = 11;
			mesh.max_density_el_order = 1;
			mesh.initial_dcount_per_el = 16;
			unsigned int d_per_line = round(sqrt(mesh.initial_dcount_per_el));
			mesh.density_subdivisions = {d_per_line*mesh.subdivisions[0], d_per_line*mesh.subdivisions[1]};


			mesh.source_fn = source_function;

			//Define loads
			std::string loadType = "distLoad";
			if (loadType == "pointLoad"){
				mesh.boundary_indicator = get_boundary_indicator;
				//Define point force
				std::vector<double> point = {2.0, 0.5};
				std::vector<double> source = {0, 1.0};
				mesh.point_source_vector.push_back(std::make_pair(point, source)); //make pairs and push
				   //empty dist load
			}
			else if (loadType == "distLoad"){

				mesh.boundary_indicator = get_boundary_indicator_force_inv;
				mesh.point_source_vector.clear();	//no point load
			}
		}
	}
	else if (dim ==3){

	}

	clock_t begin = clock();
	opt.start_time = double (begin)/CLOCKS_PER_SEC;
	opt.optimize();
	clock_t end = clock();
	double elapsed_secs = double (end - begin)/CLOCKS_PER_SEC;
	std::cout<<"SUCCESS......Computing time : "<<elapsed_secs<<std::endl;
}
