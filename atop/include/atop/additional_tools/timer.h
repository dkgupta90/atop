/*
 *
 *  Created on: Nov 7, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#ifndef INCLUDE_ATOP_ADDITIONAL_TOOLS_TIMER_H_
#define INCLUDE_ATOP_ADDITIONAL_TOOLS_TIMER_H_

#include <ctime>
/*
 * This class contains functions for timing selectively different parts of a program
 */

namespace atop{

	class Timer{
	public:


		clock_t begin, resume_clock;

		double totalTime;
		Timer();
		void pause();
		void resume();
		void stop();
	};
}



#endif /* INCLUDE_ATOP_ADDITIONAL_TOOLS_TIMER_H_ */
