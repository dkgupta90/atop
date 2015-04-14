/*
 *
 *  Created on: Mar 15, 2015
 *      Author: Deepak K. Gupta
 *  
 */
#include <string>
#include <atop/fem/define_mesh.h>

using namespace atop;

//Constructor for getting only number of physical dimensions
DefineMesh::DefineMesh(int dim, int dofs_per_node){
	this->dim = dim;
	this->dofs_per_node = dofs_per_node;
	input_file = "";
	source_fn = NULL;
	point_source_vector.clear();
}

void DefineMesh::showdata(){
	source_fn({2.0, 2.0});
}
