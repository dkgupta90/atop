/*
 * MatrixVector.h
 *
 *  Created on: Jul 27, 2014 at Precision & Microsystems Engineering group, TU Delft
 *      Author: Deepak K. Gupta
 */

#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_


#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>



namespace topopt{
using namespace dealii;
	class MatrixVector{
	public:
		void vector_matrix_multiply(
				Vector<double> &array,
				FullMatrix<double> &matrix,
				Vector<double> &output,
				int i_length,
				int j_length);

		double vector_vector_inner_product(
				Vector<double> &vector1,
				Vector<double> &vector2
				);
	};
}



#endif /* MATRIXVECTOR_H_ */
