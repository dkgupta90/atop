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

namespace atop{
	class Projection{
	public:
		std::string projection_type;
		double radius;
		double gamma;

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
		~Projection();

	};
}



#endif /* PROJECTION_H_ */
