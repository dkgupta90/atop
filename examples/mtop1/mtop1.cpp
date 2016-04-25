/*
 *
 *  Created on: Feb 8, 2016
 *      Author: Deepak K. Gupta
 *  
 */

/*
 *
 *  Created on: Apr 2, 2015
 *      Author: Deepak K. Gupta
 *
 */

/*In this step, we test the moving design points based method.
 * The algorithm starts with 1 design point per element and these points are allowed to move in space.
 * The radius of projection of each of these points is adapted.
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

using namespace atop;

std::vector<double> source_function(std::vector<double> X){
	//This function returns the value of the source function in the whole domain
	std::vector<double> output_source = {0.0, 0.0};
	return output_source;
}

unsigned int get_boundary_indicator(std::vector<double> X){
	//This function defines the boundary indicators

	if (fabs(X[0] - 0) < 1e-12)
		return 42;
	else
		return 0;

}

int main(){
	using namespace atop;

	//Define the mesh
	DefineMesh<2> mesh(2);
	mesh.coordinates = {{0, 2}, {0, 1}};
	mesh.subdivisions = {40, 20};
	mesh.density_subdivisions = {40, 20};
	mesh.coupling = false;
	mesh.source_fn = source_function;
	mesh.boundary_indicator = get_boundary_indicator;
	mesh.meshType = "subdivided_hyper_rectangle";
	mesh.elementType = "FE_Q";
	//mesh.density_elementType = "FE_DGQ";
	mesh.el_order = 1;
	//mesh.density_el_order = 1;
	mesh.adaptivityType = "movingDesignPoints";

	//Define point force
	std::vector<double> point = {2.0, 0.5};
	std::vector<double> source = {0, 1.0};
	mesh.point_source_vector.push_back(std::make_pair(point, source)); //make pairs and push

	//Define the FE related parameters

	//Define the physics of the problem
	LinearElastic<2> material1;
	material1.E = 1.0;
	material1.poisson = 0.3;
	material1.planarType = "planar_stress";

	//Define the penalization scheme
	Penalize penal("SIMP");
	penal.factmin = 1e-9;
	penal.penal_power = 1.5;

	//Define the projection scheme
	Projection filter("density_filter",
		3, 2, 10);	//Implementation the Gaussian filter, so the lowest is quite high

	//Define the optimization parameters
	Optimizedesign<2> opt(mesh, penal, filter, "MMA", 1);
	opt.problem_name = "minimum_compliance";
	opt.problemType(material1);
	opt.volfrac = 0.45; //Maximum permissible volume fraction
	opt.optimize();
	std::cout<<"SUCCESS ATTAINED"<<std::endl;
}



