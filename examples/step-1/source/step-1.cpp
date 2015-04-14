/*
 *
 *  Created on: Apr 2, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/fem/define_mesh.h>
#include <atop/physics/mechanics/elastic.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/TopologyOptimization/optimizedesign.h>
#include <string>
#include <vector>

std::vector<double> source_function(std::vector<double> X){
	//This function returns the value of the source function in the whole domain
	unsigned int dim = 2;
	std::vector<double> output_source = {0.0, 0.0};
	return output_source;
}

int main(){
	using namespace atop;

	//Define the finite element mesh
	DefineMesh mesh(2, 2);
	mesh.coordinates = {{0, 2}, {0, 1}};
	mesh.subdivisions = {20, 10};
	mesh.source_fn = source_function;

	//Define point force
	std::vector<double> point = {2.0, 0.5};
	std::vector<double> source = {0, 1.0};
	mesh.point_source_vector.push_back(std::make_pair(point, source)); //make pairs and push

	mesh.showdata();
	//Define the physics of the problem
	LinearElastic material1;
	material1.E = 1.0;
	material1.poisson = 0.3;
	material1.planarType = "planar_stress";

	//Define the penalization scheme
	Penalize penal("SIMP");
	penal.factmin = 1e-3;
	penal.penal_power = 3.0;

	//Define the optimization parameters
	Optimizedesign opt();

	std::cout<<"SUCCESS"<<std::endl;
}
