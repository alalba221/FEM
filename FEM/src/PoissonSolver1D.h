#pragma once
#include"global.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "GaussInteger.h"
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>

class PoissonSolver1D
{
public:
	PoissonSolver1D(double h_partition, BasisType trial_basis_type, BasisType  test_basis_type);
	~PoissonSolver1D() {};
private:
	GaussInteger* gs;

	Eigen::MatrixXd P, T, Pb_test, Tb_test, Pb_trial, Tb_trial;
	int number_of_local_basis_test, number_of_local_basis_trial;

	//Eigen::MatrixXd A;
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b;
	BasisType trial_basis_type, test_basis_type;
	
	std::vector<BoundaryNode>* boundary_nodes;
	/// [comment]
	/// N_FE:The N for the FE, not the partition.
	/// N_partition : The N for the partition, not the FE basis functions.
	/// [/comment]
	int N_fe, N_partition;
	double h_partition, h_fe;
	
	///Mesh information for partition
	void generate_P_T_1D(double left, double right, double h_partition, BasisType basis_type);

	// Mesh information for finite element basis functions.
	void generate_Pb_Tb_1D(double left, double right, double h_partition, BasisType trial_basis_type, BasisType test_basis_type);

	/// <summary>
	/// 三个函数要参数化
	/// </summary>
	/// <param name="x"></param>
	/// <returns></returns>
	inline double coefficient_function(double x) {
		return exp(x);
	}
	inline double function_f(double x) {
		return -exp(x)*(cos(x) - 2 * sin(x) - x*cos(x) - x*sin(x));
	}

	inline double function_g(double x) {
		if (x == 0)
			return 0;
		else if (x == 1)
			return cos(1);
	}
	inline double exact_solution(double x) {
		return x * cos(x);
	}
	double exact_solution_derivative(double x) {
		return cos(x) - x*sin(x);
	}
		

		
	// x: the coordinate of the point where we want to evaluate the local FE basis function.
	// basis_type : the type of the FE.
	// basis_type = 101 : 1D linear FE.
	// basis_index : the index of basis function to specify which basis function we want to use.
	// derivative_degree : the derivative degree of the FE basis function.
	double local_basis_1D(double x, const std::vector<double>& vertices, BasisType basis_type, int basis_index, int derivative_degree);
	
	
	void assemble_matrix_from_1D_integral(int number_of_partition, int test_derivative_degree, int trial_derivative_degree,
		int number_of_trial_local_basis,
		int number_of_test_local_basis);

	void assemble_vector_from_1D_integral(int number_of_partition, int test_derivative_degree,
		int number_of_test_local_basis);

	void generate_boundary_nodes_1D();

	void treat_Dirichlet_boundary_1D();
	// 在以veretices为界的mesh上的局部FE function
	double local_FE_function_1D(double x, int mesh, const Eigen::VectorXd & uh, BasisType basis_type, int derivative_degree){
		double result = 0;
		std::vector<double>vertices(2);
		vertices[0] = P(0, T(0, mesh));
		vertices[1] = P(0, T(1, mesh));
		double lowerbound = std::min(vertices[0], vertices[1]);
		double upperbound = std::max(vertices[0], vertices[1]);
		vertices[0] = lowerbound;
		vertices[1] = upperbound;
		for (int k = 0; k < this->number_of_local_basis_test; k++) {
			result += uh(Tb_trial(k, mesh)) * local_basis_1D(x,vertices,basis_type,k,derivative_degree);
		}
		return result;
	}
	double FE_solution_error_1D(const Eigen::VectorXd& uh, const std::function<double(double)>& exact_solution, int derivative_degree, BasisType test_basis_type) {
		double result = 0;
		for (int n = 0; n < N_partition; n++) {
			std::vector<double>vertices(2);
			vertices[0] = P(0, T(0, n));
			vertices[1] = P(0, T(1, n));
			double lowerbound = std::min(vertices[0], vertices[1]);
			double upperbound = std::max(vertices[0], vertices[1]);
			vertices[0] = lowerbound;
			vertices[1] = upperbound;

			gs->set_bound(lowerbound, upperbound);
			auto integral_Func = [this,&exact_solution,&n,&uh,&test_basis_type,&derivative_degree]
			(double x)->double {return (exact_solution(x)-local_FE_function_1D(x,n,uh,test_basis_type,derivative_degree))* 
				(exact_solution(x) - local_FE_function_1D(x, n, uh, test_basis_type, derivative_degree)); };
			result += gs->integral(integral_Func);
		}
		return result;
	}
	
public:	
	void solve() {
		generate_P_T_1D(0, 1, this->h_partition, BasisType::_1D_linear);
		generate_Pb_Tb_1D(0, 1, this->h_partition, this->trial_basis_type,this->test_basis_type);
		assemble_matrix_from_1D_integral(N_partition, 1, 1, this->number_of_local_basis_trial, this->number_of_local_basis_test);
		// std::cout << A << std::endl;
		assemble_vector_from_1D_integral(N_partition, 0, this->number_of_local_basis_test);
		//std::cout << b << std::endl;
		generate_boundary_nodes_1D();
		treat_Dirichlet_boundary_1D();
		Eigen::VectorXd x;
		//std::cout << x << std::endl;

		// solve Ax = b
		Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		if (solver.info() != Eigen::Success) {
			// decomposition failed
			std::cout << "decomposition failed";
			return;
		}
		x = solver.solve(b);
		if (solver.info() != Eigen::Success) {
			// solving failed
			std::cout << "solving failed";
			return;
		}
		//std::function<double(double)> a = std::bind(&exact_solution, this, std::placeholders::_1);
		double L2error = FE_solution_error_1D(x, std::bind(&PoissonSolver1D:: exact_solution, this, std::placeholders::_1), 0, test_basis_type);
		std::cout << sqrt(L2error) << std::endl;

		double H1error = FE_solution_error_1D(x, std::bind(&PoissonSolver1D::exact_solution_derivative, this, std::placeholders::_1), 1, test_basis_type);
		std::cout << sqrt(H1error) << std::endl;
	}
};

