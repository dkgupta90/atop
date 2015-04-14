/*
 *
 *  Created on: Apr 2, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef OPTIMIZEDESIGN_H_
#define OPTIMIZEDESIGN_H_

#include <atop/fem/define_mesh.h>
#include <atop/TopologyOptimization/penalization.h>
#include <atop/physics/elasticity.h>


namespace atop{
	class Optimizedesign{
	public:
		DefineMesh mesh;
		Penalize penal;
	};
}



#endif /* OPTIMIZEDESIGN_H_ */
