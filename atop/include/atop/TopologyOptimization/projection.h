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

namespace atop{
	class Projection{
	public:
		std::string projection_type;
		double radius;
		double gamma;
		double fact, minFact, maxFact;	//factors for calculating projection radius w.r.t element size
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

		~Projection();

	};
}



#endif /* PROJECTION_H_ */
