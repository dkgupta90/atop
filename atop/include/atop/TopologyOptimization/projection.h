/*
 *
 *  Created on: Aug 5, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef PROJECTION_H_
#define PROJECTION_H_

#include <iostream>
#include <string>
#include <atop/fem/define_mesh.h>
#include <atop/TopologyOptimization/cell_prop.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>


namespace atop{
	class Projection{
	public:
		std::string projection_type;
		std::string adaptivityType;
		double radius, true_radius;
		double gamma;
		double fact, minFact, maxFact;	//factors for calculating projection radius w.r.t element size
		unsigned int cycle;
		//Above it is assumed that elements are square or regular hexahedrons for now.

		/**
		 * Constructor to choose the type of projection operator
		 * Allows to choose the radius of projection and minimum feature size
		 */
		Projection(std::string,
				double);

		/**
		 * Last parameter of this constructor denotes 'gamma'.
		 * This is used to adapt the size of the projection radius based on the size of the cell
		 */
		Projection(std::string,
				double,
				double);

		/*
		 * The constructor below is being written for decoupled mesh with type movingDesginPoints
		 * The inputs here are the factors which multiply with element length (square) to determine radius
		 */
		Projection(
				std::string,
				double,
				double,
				double);

		/*
		 * This constructor is currently specifically tailored for dp-refinement and the projection length is the
		 * distance between the design points
		 */
		Projection(std::string,
				std::string,
				double);

		Projection(std::string,
				std::string,
				double,
				double);

		void update_projections(std::vector<CellInfo> &cell_info_vector,
				hp::DoFHandler<2> &dof_handler);

		void update_projection(std::vector<CellInfo> &cell_info_vector);

		~Projection();



	};
}



#endif /* PROJECTION_H_ */
