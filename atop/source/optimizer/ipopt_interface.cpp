/*
 *
 *  Created on: Jan 10, 2017
 *      Author: Deepak K. Gupta
 *  
 */

#include <atop/optimizer/ipopt_interface.hpp>
#include <atop/TopologyOptimization/optimizedesign.h>
#include <cassert>
#include <iostream>

using namespace Ipopt;
using namespace atop;

// constructor
IpOpt_IF::IpOpt_IF(Optimizedesign<2> &obj_opt){
	std::cout<<"Optimization object assigned "<<std::endl;
	this->opt_design2D = &obj_opt;	//Initializing the TO object to IpOpt
}

//destructor
IpOpt_IF::~IpOpt_IF(){
	delete this;
}

// returns the size of the problem
bool IpOpt_IF::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	//No. of design variables
	n = opt_design2D->design_count;

    // Only one equality constraint for now: volume constraint
    m = opt_design2D->no_constraints;

    // in this example the jacobian is dense and contains 8 nonzeros
    nnz_jac_g = n * m;

    // the hessian is also dense and has 16 total nonzeros, but we
    // only need the lower left corner (since it is symmetric)
    //nnz_h_lag = 10;	//We use the quasi-newton and do not provide the Hessian

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool IpOpt_IF::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == opt_design2D->design_count);
  assert(m == opt_design2D->no_constraints);

  // the variables have lower and upper bounds of be set
  for (Index i=0; i < n; i++) {
    x_l[i] = (*opt_design2D).lb[i];
    x_u[i] = (*opt_design2D).ub[i];
  }

  //Set the upper and lower bounds for the constraint

  // the first constraint g1 has a lower bound of 25
  g_l[0] = 0.0;
  g_u[0] = opt_design2D->volfrac;	//volume constraint

  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  //g_u[0] = 2e19;

  return true;
}


// returns the initial point for the problem
bool IpOpt_IF::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  for (Index i=0; i < n; i++) {
    x[i] = opt_design2D->volfrac;
  }
  return true;
}


// returns the value of the objective function
bool IpOpt_IF::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	assert(n == opt_design2D->design_count);

	opt_design2D->obj_fem->cycle = opt_design2D->cycle;
	std::vector<double> x_value(opt_design2D->design_count);

	  for (Index i=0; i < n; i++) {
	    x_value[i] = (double)x[i];
	  }
	opt_design2D->update_design_vector(opt_design2D->design_vector, x_value);
	opt_design2D->run_system();
	grad_obj = (opt_design2D->grad_vector);
	obj_value = opt_design2D->objective;

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IpOpt_IF::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == opt_design2D->design_count);
  for (Index i=0; i < n; i++) {
    grad_f[i] = grad_obj[i];
  }


  return true;
}


// return the value of the constraints: g(x)
bool IpOpt_IF::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == opt_design2D->design_count);
  assert(m == opt_design2D->no_constraints);

	g[0] = opt_design2D->vol_constraint.volumeConstraint(
			grad_cons,
			opt_design2D->cell_info_vector,
			opt_design2D->density_cell_info_vector,
			opt_design2D->obj_fem->density_field,
			*(opt_design2D->mesh)
			);

  return true;
}

// return the structure or values of the jacobian
bool IpOpt_IF::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    // this particular jacobian is dense
	  for (Index i=0; i < n; i++) {
	    iRow[i] = 0;
	    jCol[i] = i;
	  }
  }
  else {
    // return the values of the jacobian of the constraints

	  for (Index i=0; i < n; i++) {
	    values[i] = grad_cons[i];
	  }
  }

  return true;
}

//return the structure or values of the hessian
bool IpOpt_IF::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx=0;
    for (Index row = 0; row < 4; row++) {
      for (Index col = 0; col <= row; col++) {
        iRow[idx] = row;
        jCol[idx] = col;
        idx++;
      }
    }

    assert(idx == nele_hess);
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    // fill the objective portion
    values[0] = obj_factor * (2*x[3]); // 0,0

    values[1] = obj_factor * (x[3]);   // 1,0
    values[2] = 0.;                    // 1,1

    values[3] = obj_factor * (x[3]);   // 2,0
    values[4] = 0.;                    // 2,1
    values[5] = 0.;                    // 2,2

    values[6] = obj_factor * (2*x[0] + x[1] + x[2]); // 3,0
    values[7] = obj_factor * (x[0]);                 // 3,1
    values[8] = obj_factor * (x[0]);                 // 3,2
    values[9] = 0.;                                  // 3,3


    // add the portion for the first constraint
    values[1] += lambda[0] * (x[2] * x[3]); // 1,0

    values[3] += lambda[0] * (x[1] * x[3]); // 2,0
    values[4] += lambda[0] * (x[0] * x[3]); // 2,1

    values[6] += lambda[0] * (x[1] * x[2]); // 3,0
    values[7] += lambda[0] * (x[0] * x[2]); // 3,1
    values[8] += lambda[0] * (x[0] * x[1]); // 3,2

    // add the portion for the second constraint
    values[0] += lambda[1] * 2; // 0,0

    values[2] += lambda[1] * 2; // 1,1

    values[5] += lambda[1] * 2; // 2,2

    values[9] += lambda[1] * 2; // 3,3
  }

  return false;
}


void IpOpt_IF::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i<n; i++) {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
}

