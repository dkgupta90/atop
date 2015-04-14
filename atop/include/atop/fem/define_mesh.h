/*
 *
 *  Created on: Mar 15, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef DEFINE_MESH_H_
#define DEFINE_MESH_H_

#include <iostream>
#include <vector>
#include <string>

namespace atop{
//using namespace dealii;
	class DefineMesh{

	public:
		int dim; //Dimensions of the mesh
		unsigned int dofs_per_node;

		//Saves the coordinates for physical dimensions of the mesh
		std::vector<std::vector<double> > coordinates;

		//Provides number of subdivisions along each physical dimension of the mesh
		std::vector<unsigned int> subdivisions;

		//Saves the path of the input mesh (if supplied, else empty)
		std::string input_file;

		//Holds all coordinates and magnitudes of all point sources
		std::vector< std::pair< std::vector<double> , std::vector<double> > > point_source_vector;

		//Constructor for getting only number of physical dimensions
		DefineMesh(int, int);

		//Pointer to the source function supplied by the user
		std::vector<double> (*source_fn)(std::vector<double>);

		//Pointer to the function for defining the boundaries


		void showdata();
	};

}


#endif /* DEFINE_MESH_H_ */
