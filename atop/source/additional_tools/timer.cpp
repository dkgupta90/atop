/*
 *
 *  Created on: Nov 7, 2016
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/additional_tools/timer.h>
#include <ctime>
#include <iostream>

using namespace atop;

Timer::Timer(){

	this->begin = clock();
	resume_clock = begin;
	this->totalTime = 0.0;
}

void Timer::pause(){
	clock_t cur_clock = clock();
	this->totalTime += ((double)(cur_clock - resume_clock)/CLOCKS_PER_SEC);

}

void Timer::resume(){
	this->resume_clock = clock();
}

void Timer::stop(){
	clock_t cur_clock = clock();
	this->totalTime += ((double)(cur_clock - resume_clock)/CLOCKS_PER_SEC);
	std::cout<<"Timer::algorithm time : "<<totalTime<<std::endl;

	double actualTime = (double)(cur_clock - this->begin)/CLOCKS_PER_SEC;
	std::cout<<"Timer::total time : "<<actualTime<<std::endl;

}
