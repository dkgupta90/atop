/*
 * penalization.cpp

 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <atop/TopologyOptimization/penalization.h>
#include <atop/TopologyOptimization/cell_prop.h>
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

Penalize::Penalize(){}

Penalize::Penalize(std::string scheme){
	//This function initializes the penalization scheme
	this->scheme = scheme;
}

void Penalize::update_param(
		double E0,
		std::vector<CellInfo> &cell_info_vector){
	double Emin = E0 * factmin;
	unsigned int no_cells = cell_info_vector.size();

	//Iterating over every cell
	for(unsigned int i = 0; i < no_cells; ++i){
		unsigned int design_count =cell_info_vector[i].density.size();
		cell_info_vector[i].E_values.resize(design_count);
		cell_info_vector[i].dE_values.resize(design_count);

		for(unsigned int q_point = 0; q_point < design_count; ++q_point){
			double density = cell_info_vector[i].density[q_point];
			double Evalue, dEvalue;
			if (scheme == "SIMP"){
				Evalue = (Emin + (E0 - Emin)*(pow(density, penal_power)));
				dEvalue = penal_power * ((E0 - Emin)*(pow(density, penal_power-1)));
			}
			else if (scheme == "RAMP"){
					double denom = 1 + (penal_power * (1 - density));
					Evalue = Emin + (density / denom) * (E0 - Emin);
					dEvalue = ((1 + penal_power)/(denom * denom)) * (E0 - Emin);
			}

			cell_info_vector[i].E_values[q_point] = Evalue;
			cell_info_vector[i].dE_values[q_point] = dEvalue;
		}
	}
}

double Penalize::penalized_factor(double xPhys){

	double E0 = 1.0;
	double Emin = E0 * factmin;

	if (scheme == "SIMP"){
		double Evalue = (Emin + (E0 - Emin)*(pow(xPhys, penal_power)));
		return Evalue;
	}
	else if (scheme == "RAMP"){
			double denom = 1 + (penal_power * (1 - xPhys));
			double Evalue = Emin + (xPhys / denom) * (E0 - Emin);
			return Evalue;
	}
}

