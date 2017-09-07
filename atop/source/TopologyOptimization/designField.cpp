/*
 *
 *  Created on: Apr 20, 2016
 *      Author: Deepak K. Gupta
 *  
 */


#include <atop/TopologyOptimization/designField.h>
#include <math.h>
#include <iostream>
#include <fstream>

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
		if (isPerfectSquare(no_points) == true){
			unsigned int dxcount = (unsigned int)(round(sqrt(no_points)));
			double dx = 2.0/(dxcount);

			for (unsigned int i = 0; i < dxcount; i++){
				for (unsigned int j = 0; j < dxcount; j++){
					pointX[i*dxcount + j][0] = -1 + (j+0.5)*dx;
					pointX[i*dxcount + j][1] = -1 + (i+0.5)*dx;
				}
			}
		}
		else{
			std::ifstream rfile;
			rfile.open("designField/designField.dat", std::ios::in);

			if (!rfile.is_open()){
				std::cerr<<"designField.dat not found\n";
			}

			std::string temps;
			for (unsigned int i = 1; i <= no_points-1; i++)
			        std::getline(rfile, temps);
			double k;
			rfile>>k;	//reading the no.of points

			if (fabs((double)no_points - k) > 1e-12){

				std::cerr<<"Mismatch in reading the no. of points in file : "<<no_points<<"  "<<k<<"\n";
				exit(0);
			}
			for (unsigned int pt = 0; pt < no_points; ++pt){
				rfile>>pointX[pt][0];
				rfile>>pointX[pt][1];
				//std::cout<<pt<<" "<<pointX[pt][0]<<"   "<<pointX[pt][1]<<std::endl;
			}

		}
	}
}

void DesignField::update_no_points(unsigned int no_design_points){
	this->no_points = no_design_points;
	unsigned int dim = 2;//pointX[0].size();
	rho.clear();
	pointX.clear();
	rho.resize(no_points);
	pointX.resize(no_points);

	for (unsigned int i = 0; i < rho.size(); ++i){
		pointX[i].clear();
		pointX[i].resize(dim, 0.);
	}

}

void DesignField::update_pseudo_designWeights(unsigned int max_design_points_per_cell,
				std::vector<std::vector<double> > &dp_PointX,
				double cell_area){
	unsigned int dim = pointX[0].size();


	if (dim == 2){
		//Calculating the pseudo-filter radius
		unsigned int max_d_factor = round(sqrt(max_design_points_per_cell));
		unsigned int d_factor = floor(sqrt(dp_PointX.size()));

		//Below, 2.0 is used as the length of the cell, since all the design points defined within the cell
		//assume that the cell is 2X2 in length and are relative to it with center at 0,0.
		double rmin = (((double)2.0)/(d_factor*sqrt(2.0))) * 1.05;	//5% tolerance added
		/*
		 * This is a 2 dimensional vector with the first dimension iterating over the number of pseudo-design points
		 * and the second dimension iterating over the number of actual non-uniformly distributed design points
		 */
		dx_drho.clear();
		dx_drho.resize(max_design_points_per_cell);

		for(unsigned int j = 0; j < max_design_points_per_cell; j++){
			dx_drho[j].clear();
			dx_drho[j].resize(dp_PointX.size(), 0.0);
			double sum_weights = 0.0;

			for (unsigned int k = 0; k < dp_PointX.size(); ++k){
				double distance = 0.0;
				distance = pow(pointX[j][0] - dp_PointX[k][0], 2);
				distance += pow(pointX[j][1] - dp_PointX[k][1], 2);
				distance = sqrt(distance);
				if (rmin <= distance)	continue;

				dx_drho[j][k] = fabs(rmin - distance);
				sum_weights+= fabs(rmin - distance);
			}

			//Dividing by the weights
			if (sum_weights == 0.0){
				std::cerr<<"DesignField::update_pseudo_designWeights : Zero sum_weights error"<<std::endl;
				exit(0);
			}
			for (unsigned int k = 0; k < dp_PointX.size(); ++k){
				//Note that these derivatives can also be used to sum the value of any point in the pseudo-design mesh
				dx_drho[j][k] /= sum_weights;
			}

		}
	}


}

void DesignField::update_pseudo_designField(
		std::vector<double> &dp_rho){

	for(unsigned int j = 0; j < rho.size(); j++){

		if (dx_drho[j].size() != dp_rho.size()){
			std::cerr<<"DesignField::update_psuedo_designField - Dimension Mismatch "<<std::endl;
			exit(0);
		}

		rho[j] = 0.0;
		for (unsigned int k = 0; k < dp_rho.size(); ++k){
			rho[j] += (dx_drho[j][k] * dp_rho[k]);
		}
	}
}

bool DesignField::isPerfectSquare(unsigned int no_design_points){
	double sqroot = sqrt((double)no_design_points);
	if (fabs(sqroot - (double)(floor(sqroot))) < 1e-12){
		return true;
	}
	else{
		return false;
	}
}
