#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "GaussInteger.h"
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>
enum class  BoundaryConditionType{
	Dirichlet
};
enum class BasisType {
	_1D_linear,
	_1D_quadratic
};

struct BoundaryNode {
	BoundaryConditionType bc_type;
	int global_index = -1;
	BoundaryNode(BoundaryConditionType type, int index):
		bc_type(type), global_index(index){
	}
};

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
	// x: the coordinate of the point where we want to evaluate the local FE basis function.
	// basis_type : the type of the FE.
	// basis_type = 101 : 1D linear FE.
	// basis_index : the index of basis function to specify which basis function we want to use.
	// derivative_degree : the derivative degree of the FE basis function.
	double local_basis_1D(double x, const std::vector<double>& vertices, BasisType basis_type, int basis_index, int derivative_degree);
	
	
	void assemble_matrix_from_1D_integral(int number_of_element, int test_derivative_degree, int trial_derivative_degree,
		int number_of_trial_local_basis,
		int number_of_test_local_basis);

	void assemble_vector_from_1D_integral(int number_of_elements, int test_derivative_degree,
		int number_of_test_local_basis);


	void generate_boundary_nodes_1D();

	void treat_Dirichlet_boundary_1D();
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
		std::cout << x << std::endl;
	}
};

