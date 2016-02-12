/*
 *
 *  Created on: Aug 5, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/TopologyOptimization/projection.h>
#include <atop/fem/define_mesh.h>
#include <math.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

using namespace atop;

Projection::Projection(
		std::string type,
		double r,
		double g){
	projection_type = type;
	radius = r;
	gamma = g;
}

Projection::Projection(
		std::string type,
		double r){
	projection_type = type;
	radius = r;
	gamma = 1.0;
}

Projection::~Projection(){

}

Projection::Projection(
		std::string type,
		double r,
		double minR,
		double maxR){
	projection_type = type;
	this->fact = r;	//factor for current radius
	this->minFact = minR;	//factor for minimum projection
	this->maxFact = maxR;	//factor for maximum projection
}
