/*
 *
 *  Created on: Aug 19, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/optimizer/optimality_criteria.h>
#include <atop/TopologyOptimization/optimizedesign.h>
#include <math.h>
#include <cstdlib>

using namespace atop;

OC::OC(unsigned int n){
	no_design = n;
}

void OC::set_lower_bounds(
		std::vector<double> &vec1){
	this->lb = &vec1;
}

void OC::set_upper_bounds(
		std::vector<double> &vec1){
	this->ub = &vec1;
}

void OC::optimize(
		std::vector<double> &vec1){
	this->design_vector = &vec1;

	double old_objective, objective = 9999999999;
	Optimizedesign<2> *opt_design2d = static_cast<Optimizedesign<2>*>(obj_data);



	do{
		old_objective = objective;
		std::cout<<"Iteration : "<<opt_design2d->obj_fem->itr_count + 2<<std::endl;
		std::cout<<"OC method initiated"<<std::endl;
		std::vector<double> obj_grad(no_design), old_design_vector;
		old_design_vector.clear();
		old_design_vector = *design_vector;
		obj_grad.clear();

		//Calculating the objective and derivative
		objective = obj_fn(
				*design_vector,
				obj_grad,
				obj_data);
		double current_volfrac = 0.0;

/*		for (unsigned int i = 0; i < obj_grad.size(); ++i){
			std::cout<<obj_grad[i]<<std::endl;
		}*/


		//Calculating volume derivatives
		std::vector<double> vol_grad(obj_grad.size(), 0.0);
		current_volfrac = opt_design2d->vol_constraint.volumeConstraint(
				vol_grad,
				opt_design2d->cell_info_vector,
				opt_design2d->density_cell_info_vector,
				opt_design2d->obj_fem->density_field,
				*(opt_design2d->mesh)
				);




		double l1 = 1e-5, l2 = 100000000, move = 0.08;

/*		if (opt_design2d->obj_fem->itr_count == 0){
			l1 = l2;

			//updating the density_cell_info_vector
			if (opt_design2d->mesh->coupling == false){
				opt_design2d->obj_fem->density_field.update_density_cell_info_vector(
						opt_design2d->cell_info_vector,
						opt_design2d->density_cell_info_vector,
						*design_vector);
			}
			else{
				opt_design2d->obj_fem->density_field.update_density_cell_info_vector(
						opt_design2d->density_cell_info_vector,
						*design_vector);
			}

			//update the pseudo-design field
			opt_design2d->obj_fem->update_pseudo_designField();
			opt_design2d->obj_fem->add_density_to_design_cell_info_vector();

			//Applying smoothing
			opt_design2d->obj_fem->density_field.smoothing(
					opt_design2d->cell_info_vector,
					opt_design2d->density_cell_info_vector,
					*(opt_design2d->mesh));
			opt_design2d->obj_fem->density_field.smoothing(
					opt_design2d->cell_info_vector,
					opt_design2d->density_cell_info_vector);

			//Computing current volume fraction
			//Note that the elastic_data object makes it hardcoded, it needs to be removed
			current_volfrac = opt_design2d->obj_fem->density_field.get_vol_fraction(
					opt_design2d->density_cell_info_vector);
		}*/



		while ((l2 - l1) > 1e-4){
			double lmid = 0.5*(l1  + l2);
			for (unsigned int i = 0; i < design_vector->size(); ++i){

				//std::cout<<obj_grad[i]<<"   "<<vol_grad[i]<<std::endl;
				double octemp1;
				double old_density = old_design_vector[i];
				if (obj_grad[i] > 0){
					std::cout<<"Error located "<<obj_grad[i]<<std::endl;
					obj_grad[i] = 0.0;
				}

				//For compliant mechanism
				if (opt_design2d->problem_name == "compliant_mechanism" && obj_grad[i] > 0){
					obj_grad[i] = -1e-10;
					octemp1 = old_density * (pow(-(obj_grad[i]/(lmid * vol_grad[i])), 0.3));

				}
				else{
					octemp1 = old_density * (sqrt(-(obj_grad[i]/(lmid * vol_grad[i]))));

				}
				//std::cout<<octemp1<<std::endl;
				if ((old_density + move) < octemp1){
					octemp1 = old_density + move;
				}
				if (octemp1 > (*ub)[i]){
					octemp1 = (*ub)[i];
				}
				if ((old_density - move) > octemp1){
					octemp1 = old_density - move;
				}
				if (octemp1 < (*lb)[i] + 1e-9){
					octemp1 = (*lb)[i] + 1e-9;
				}

				(*design_vector)[i] = octemp1;

			}


			//updating the density_cell_info_vector
			if (opt_design2d->mesh->coupling == false){
				opt_design2d->obj_fem->density_field.update_density_cell_info_vector(
						opt_design2d->cell_info_vector,
						opt_design2d->density_cell_info_vector,
						*design_vector);
			}
			else{
				opt_design2d->obj_fem->density_field.update_density_cell_info_vector(
						opt_design2d->density_cell_info_vector,
						*design_vector);
			}

			//update the pseudo-design field
			opt_design2d->obj_fem->update_pseudo_designField();
			opt_design2d->obj_fem->add_density_to_design_cell_info_vector();

			//Applying smoothing
/*			opt_design2d->obj_fem->density_field.smoothing(
					opt_design2d->cell_info_vector,
					opt_design2d->density_cell_info_vector,
					*(opt_design2d->mesh));*/
			opt_design2d->obj_fem->density_field.smoothing(
					opt_design2d->cell_info_vector,
					opt_design2d->density_cell_info_vector);

			//Computing current volume fraction
			//Note that the elastic_data object makes it hardcoded, it needs to be removed
			current_volfrac = opt_design2d->obj_fem->density_field.get_vol_fraction(
					opt_design2d->density_cell_info_vector);

			//std::cout<<(density_sum/(x_count*y_count))<<std::endl;
			if (current_volfrac < 0 || current_volfrac > 1){
				std::cerr<<"ERROR!! Current_volfrac : "<<current_volfrac;
				exit(0);

			}
			if((current_volfrac - opt_design2d->volfrac) > 0){
				l1 = lmid;
			}
			else{
				l2 = lmid;
			}
		}
		//std::cout<<density_sum<<std::endl;
		std::cout<<"Volfrac: "<<current_volfrac<<std::endl;
		//exit(0);
	}while(fabs(old_objective - objective) > min_obj_change && opt_design2d->obj_fem->itr_count < 1500);


}
