/*
 * penalization.cpp

 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <atop/TopologyOptimization/penalization.h>
#include<vector>
#include <math.h>
#include <iostream>
using namespace topopt;

using namespace atop;

Penalization::Penalization(){}

void Penalization::set_param(double E0, double Emin,
		std::vector<double> &E_values,
		std::vector<double> &dE_values,
		std::vector<double> &density_values,
		unsigned int penalization_model,
		double penal_power){
	E_values.resize(density_values.size());
	dE_values.resize(density_values.size());
	unsigned int n_points = density_values.size();
	for (int i = 0; i < n_points; i++){
		E_values[i] = (Emin + (E0 - Emin)*(pow(density_values[i], penal_power)));
		dE_values[i] = penal_power * ((E0 - Emin)*(pow(density_values[i], penal_power-1)));

	}
}


Penalize::Penalize(std::string scheme){
	//This function initializes the penalization scheme
	this->scheme = scheme;
}

