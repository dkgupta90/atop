/*
 * MatrixVector.cpp
 *
 *  Created on: Jul 27, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <atop/math_tools/algebra/MatrixVector.h>

using namespace dealii;
using namespace topopt;

void MatrixVector::vector_matrix_multiply(
		Vector<double> &array,
		FullMatrix<double> &matrix,
		Vector<double> &output,
		int i_length,
		int j_length){
	int array_length = array.size();
	if (array_length != i_length)
		std::cerr<<"Compatibility Issue!! Vector and matrix not compatible for multiplication"<<std::endl;
	for(int j = 0; j < j_length; ++j){
		double sum = 0.0;
		for(int i = 0; i < i_length; ++i){
			sum += (array(i) * matrix(i, j));
		}
		output(j) = sum;
	}
}

double MatrixVector::vector_vector_inner_product(
		Vector<double> &vector1,
		Vector<double> &vector2){
	if (vector1.size() != vector2.size()){
		std::cerr<<"Vector dimensions do not match for inner product, hence exiting"<<std::endl;
		exit(0);
	}
	double sum = 0.0;
	for (int i = 0; i < vector1.size(); ++i){
		sum += (vector1(i)*vector2(i));
	}

	return sum;
}


