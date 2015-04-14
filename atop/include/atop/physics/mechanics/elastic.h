/*
 *
 *  Created on: Apr 2, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef ELASTIC_H_
#define ELASTIC_H_

#include <string>
//using namespace dealii;

namespace atop{
	class LinearElastic{
	public:
		unsigned int E, poisson; //Define properties of the material
		std::string planarType;  //For 2D problems
	};
}




#endif /* ELASTIC_H_ */
