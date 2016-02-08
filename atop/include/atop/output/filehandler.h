/*
 *
 *  Created on: Nov 8, 2015
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

#include <iostream>
#include <string>

namespace atop{
	template <int dim>
	class FileHandler{
	public:

		void addToFile(
				std::string,
				unsigned int,
				unsigned int,
				unsigned int,
				double,
				double,
				double);

	};

	template class FileHandler<2>;
}




#endif /* FILEHANDLER_H_ */
