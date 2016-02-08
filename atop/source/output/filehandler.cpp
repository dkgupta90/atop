/*
 *
 *  Created on: Nov 8, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/output/filehandler.h>
#include <fstream>

using namespace atop;

template <int dim>
void FileHandler<dim>::addToFile(
		std::string filename,
		unsigned int cycle,
		unsigned int iter,
		unsigned int no_dofs,
		double output_value,
		double iter_time,
		double current_time
		){

	//Creating/Opening the file to add the data
	std::ofstream myfile1;
	if(cycle == 1 && iter == 1){
		myfile1.open(filename.c_str(), std::ofstream::out);
	}
	else{
		myfile1.open(filename.c_str(), std::ofstream::app);
	}

	//Adding the information to the opened file
	myfile1<<cycle<<"\t"<<iter<<"\t"<<no_dofs<<"\t"<<output_value<<"\t"<<current_time<<"\t"<<iter_time<<"\n";

	//closing the file
	myfile1.close();
}