/*
 *
 *  Created on: Aug 5, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/TopologyOptimization/projection.h>

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

