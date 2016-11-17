/*
 * design_analysis
 *
 *  Created on: Aug 29, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */
#include <iostream>
#include <vector>

#ifndef DESIGN_ANALYSIS_H_
#define DESIGN_ANALYSIS_H_

namespace topopt{
	class Design_Analysis{
	public:
		std::vector< std::vector< std::pair<unsigned int, double> > > drhoi_ddn;
		Design_Analysis(int q_points);
	};
}



#endif /* DESIGN_ANALYSIS_H_ */
