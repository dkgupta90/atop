/*
 * bitmap.h
 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 *
 * Bitmap class contains functions for generating bitmap images from arrays.
 * It can deal with several different colormaps.
 */

#ifndef BITMAP_H_
#define BITMAP_H_

#include <string>
#include <vector>
namespace topopt{
	class Bitmap{
	public:
		std::string colormap;
		double **scaled_density;
		unsigned int scaled_x_count;
		unsigned int scaled_y_count;
		Bitmap();
		void set_colormap(std::string str);
		void create_2d_to_bitmap(
				double **data,
				std::string fname,
				unsigned int x_count,
				unsigned int y_count);

		void int_upscale_image(
				double **data,
				unsigned int x_count,
				unsigned int y_count,
				unsigned int factor);
	};
}



#endif /* BITMAP_H_ */
