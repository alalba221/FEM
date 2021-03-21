#include <iostream>
#include <Eigen/Dense>
#include "PoissonSolver1D.h"
#include "GaussInteger.h"


int main()
{
	PoissonSolver1D solver(1.0/4,BasisType::_1D_linear, BasisType::_1D_linear);
	solver.solve();	

}