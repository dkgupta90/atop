/*
 * neighbors.cpp
 *
 *  Created on: Oct 16, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <topopt/TopologyOptimization/neighbors.h>
#include <iostream>
#include <vector>

using namespace topopt;

StoreIndices::StoreIndices(unsigned int i){
	stored_indices.resize(i);
	density_indices.resize(i);
}

StoreIndices::StoreIndices(){

}
