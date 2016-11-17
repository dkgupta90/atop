/*
 * bitmap.cpp
 *
 *  Created on: Jul 15, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <topopt/data_format/bitmap/bitmap_image.hpp>
#include <topopt/data_format/bitmap.h>
#include <cmath>
#include <string>

using namespace topopt;
#include <string>
#include <vector>

//Constructor for setting the default value of colormap to grayscale
Bitmap::Bitmap(){
	colormap = "gray";
}

//This functiopn converts the colormap value to lowercase
void Bitmap::set_colormap(std::string str){
	colormap = str;
	int i =0;
	while(colormap[i]){
		colormap[i] = std::tolower((char)colormap[i]);
		++i;
	}
}

void Bitmap::int_upscale_image(
		double **data,
		unsigned int x_count,
		unsigned int y_count,
		unsigned int factor){

	scaled_x_count = x_count * factor;
	scaled_y_count = y_count * factor;
	scaled_density = (double**)malloc(scaled_x_count*sizeof(double*));
	for(unsigned int i = 0; i < scaled_x_count; ++i){
		scaled_density[i] = (double*)malloc(scaled_y_count*sizeof(double));
	}

	for(unsigned int i = 0; i < scaled_x_count; ++i){
		for(unsigned int j = 0; j < scaled_y_count; ++j){
			unsigned int irem = (int)(i / factor);
			unsigned int jrem = (int)(j / factor);
			scaled_density[i][j] = data[irem][jrem];
		}
	}
}

//Function that defines the grid and fills the pixels with the user defined colormap
void Bitmap::create_2d_to_bitmap(
		double **data,
		std::string fname,
		unsigned int x_count,
		unsigned int y_count){
	bitmap_image image(scaled_x_count, scaled_y_count);
	for(unsigned int i = 0; i < scaled_x_count; ++i){
		for(unsigned int j = 0; j < scaled_y_count; j++){
			unsigned int temp_col = (int)floor((scaled_density[i][j]) * 1000);
			if (temp_col == 1000) --temp_col;
			temp_col %= 1000;
			rgb_store col;

			if (colormap.compare("jet") == 0){
				col = jet_colormap[temp_col];
			}
			else if(colormap.compare("hot") == 0){
				col = hot_colormap[temp_col];
			}
			else if(colormap.compare("hsv") == 0){
				col = hsv_colormap[temp_col];
			}
			else if (colormap.compare("gray") == 0){
				col = gray_colormap[temp_col];
			}
			else{
				col = gray_colormap[temp_col];
			}
			image.set_pixel(i, j, col.red, col.green, col.blue);
		}
	}
	image.save_image(fname);
}


