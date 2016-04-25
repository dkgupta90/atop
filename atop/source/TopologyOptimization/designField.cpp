/*
 *
 *  Created on: Apr 20, 2016
 *      Author: Deepak K. Gupta
 *  
 */


#include <atop/TopologyOptimization/designField.h>

using namespace atop;

DesignField::DesignField(){

}



void DesignField::initialize_field(
		unsigned int dim,
		unsigned int no_points,
		unsigned int design_rule,
		double volfrac){

	//no_points defines the number of points to be put inside the element
	//design_rule defines the rule for distributing/adding points in the domain
	//At every p-refinement, the new design points should not disturb the position of previous points.
	this->dim = dim;
	this->design_rule = design_rule;

	//assigning the no. of points for the density and point vectors
	rho.resize(no_points, volfrac);
	dxPhys_drho.resize(no_points);
	pointX.resize(no_points);
	for(unsigned int i = 0; i < pointX.size(); ++i){
		pointX[i].resize(dim);
	}

	//for manually assigning the points
	if (design_rule == 1){
		update_field_manual(0);
	}
}

//This version is called only during initialization, overloaded forms will follow in future
void DesignField::update_field_manual(unsigned int ctr){

	//ctr !=0 means that this is not the initialization cycle,
	//which would have meant new rho will be interpolated from existing ones

	unsigned int dim = pointX[0].size();	//getting the no. of dimensions of the problem
	unsigned int no_points = pointX.size();
	if (dim == 2){
		for (unsigned int i = 0; i < no_points; ++i){
			if (i == 0){
				pointX[i][0] = 0.0;
				pointX[i][1] = 0.0;
			}
			else if (i == 1){
				pointX[i][0] = -0.5;
				pointX[i][1] = -0.5;
			}
			else if (i == 2){
				pointX[i][0] = 0.5;
				pointX[i][1] = -0.5;
			}
			else if (i == 3){
				pointX[i][0] = 0.5;
				pointX[i][1] = 0.5;
			}
			else if (i == 4){
				pointX[i][0] = -0.5;
				pointX[i][1] = 0.5;
			}
		}
	}



}
