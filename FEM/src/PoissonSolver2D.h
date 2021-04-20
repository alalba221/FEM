#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "GaussInteger.h"
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU>
#include "global.h"
#include <array>
class PoissonSolver2D
{
public:
	PoissonSolver2D(const std::array<double,2>& x_range, const std::array<double,2>& y_range, 
		const std::array<double,2>& h_partition, 
		BasisType trial_basis_type, BasisType  test_basis_type);
	
	~PoissonSolver2D() {};
private:
	GaussInteger* gs;

	Eigen::MatrixXd P_partition, T_partition, Pb_test, Tb_test, Pb_trial, Tb_trial;
	int number_of_local_basis_test, number_of_local_basis_trial;

	//Eigen::MatrixXd A;
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b;
	BasisType trial_basis_type, test_basis_type;

	// boundary edges is along with mesh, boundary nodes are along with FE not partition
	std::vector<BoundaryEdge>* boundary_edges;
	std::vector<BoundaryNode>* boundary_nodes;
	/// [comment]
	/// N_FE:The N for the FE, not the partition.
	/// N_partition : The N for the partition, not the FE basis functions.
	/// [/comment]
	std::array<double,2> X_range, Y_range;
	std::array<int,2> N_basis, N_partition;
	std::array<double,2> h_partition, h_basis;
	int tnp;
	///Mesh information for partition
	void generate_P_T_2D(const std::array<double, 2>& x_range, const std::array<double, 2>& y_range, const std::array<double,2>& h_partition);

	// Mesh information for finite element basis functions.
	void generate_Pb_Tb_2D(const std::array<double, 2>& x_range, const std::array<double, 2>& y_range, const std::array<double, 2>& h_partition, BasisType trial_basis_type, BasisType test_basis_type);

	void generate_boundary_edges(const std::array<int, 2>& n_partition);
	void generate_boundary_nodes(const std::array<int, 2>& n_basis);

	// This is the reference FE basis function on triangle ABC where A = (0, 0), B = (1, 0) and C = (0, 1).
	double triangular_reference_basis(double x, double y, BasisType basis_type, int basis_index, int derivative_degree_x, int derivative_degree_y);
public:
	void solve() {
		/// Mesh information
		generate_P_T_2D(this->X_range, this->Y_range, this->h_partition);
		generate_boundary_edges(this->N_partition);
		for (auto n : *boundary_edges) {
			std::cout << n.mesh_index<<" "<<n.first_node_global_index<<" "<<n.second_node_global_index << std::endl;
		}
		/// FE info
		generate_Pb_Tb_2D(this->X_range, this->Y_range, this->h_partition, BasisType::_2D_Lagrange_quadratic, BasisType::_2D_Lagrange_quadratic);
		generate_boundary_nodes(this->N_basis);
		std::cout << Tb_test << std::endl;
		for (auto n : *boundary_nodes) {
			std::cout << n.global_index << std::endl;
		}
	}
};



