#include <iostream>
#include <Eigen/Dense>
#include "PoissonSolver1D.h"
#include "PoissonSolver2D.h"
#include "GaussianIntegerTriangle.h"
double Func(std::array<double, 2> p) {
	return 1.0 / (sqrt(p[0] + p[1]) * (1 + p[0] + p[1]) * (1 + p[0] + p[1]));
}

int main()
{
	//1D
	//PoissonSolver1D solver(1.0/8,BasisType::_1D_quadratic, BasisType::_1D_quadratic);
	// solver.solve();	

	 //2D
	PoissonSolver2D solver({ { 0,1 } }, { 0,1 }, { 1.0 / 2, 1.0 / 2 }, BasisType::_2D_Lagrange_quadratic, BasisType::_2D_Lagrange_quadratic);

	solver.solve();

	/// <summary>
	///  test GaussianInteferTriangle
	/// </summary>
	/// <returns></returns>

	//GaussianIntegerTriangle* gs = new GaussianIntegerTriangle(9);
	////std::array< std::array<double, 2>, 3 > bnd = {  {0,0},{1,0},{01}  };
	//gs->set_bound({ { {5,5},{2,5},{2,2} } });
	//double r = gs->integral(Func);

	//std::cout << r << std::endl;
}
