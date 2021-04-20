#include "PoissonSolver2D.h"

PoissonSolver2D::PoissonSolver2D(const std::array<double, 2>& x_range, const std::array<double, 2>& y_range,
	const std::array<double, 2>& h_partition,
	BasisType trial_basis_type, BasisType  test_basis_type)
	:X_range(x_range),Y_range(y_range),h_partition(h_partition), trial_basis_type(trial_basis_type), test_basis_type(test_basis_type)
{
	N_partition[0] = (X_range[1] - X_range[0]) / h_partition[0];
	N_partition[1] = (Y_range[1] - Y_range[0]) / h_partition[1];

	
	if (test_basis_type == BasisType::_2D_Lagrange_quadratic && trial_basis_type == BasisType::_2D_Lagrange_quadratic) {
		N_basis[0] = N_partition[0] * 2;
		N_basis[1] = N_partition[1] * 2;

		h_basis[0] = h_partition[0] / 2;
		h_basis[1] = h_partition[1] / 2;
	}else if (test_basis_type == BasisType::_2D_linear && trial_basis_type == BasisType::_2D_linear) {
		
		N_basis = N_partition;
		h_basis = h_partition;	
	}

	this->boundary_edges = new std::vector<BoundaryEdge>;
	this->boundary_nodes = new std::vector<BoundaryNode>;

}

void PoissonSolver2D::generate_P_T_2D(const std::array<double, 2>& x_range, const std::array<double, 2>& y_range, const std::array<double, 2>& h_partition)
{
	double left = x_range[0];
	double right = x_range[1];
	double bottom = y_range[0];
	double top = y_range[1];
	double h1 = h_partition[0];
	double h2 = h_partition[1];
	int tnp = (N_partition[0] + 1) * (N_partition[1] + 1);
	P_partition = Eigen::MatrixXd::Zero(2, tnp);
	T_partition = Eigen::MatrixXd::Zero(3, 2*N_partition[0]*N_partition[1]);

	for (int j = 0; j < tnp; j++) {
		int rn = j % (N_partition[1] + 1);
		int cn = j / (N_partition[1] + 1);
		
		P_partition(0, j) = left + cn * h1;
		P_partition(1, j) = bottom + rn * h2;
	}
	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(N_partition[0]+1,N_partition[1]+1);
	for (int i = 0; i < Q.rows(); i++) {
		for (int j = 0; j < Q.cols(); j++) {
			Q(i, j) = j * (N_partition[1] + 1) + i;
		}
	}

	for (int n = 0; n < (N_partition[0]*N_partition[1]); n++) {
		
		int re = n % N_partition[1];
		int ce = n / N_partition[1];
		T_partition(0, 2 * n) = Q(re, ce);
		T_partition(1, 2 * n) = Q(re, ce+1);
		T_partition(2, 2 * n) = Q(re+1, ce);

		T_partition(0, 2 * n+1) = Q(re+1, ce);
		T_partition(1, 2 * n+1) = Q(re, ce + 1);
		T_partition(2, 2 * n+1) = Q(re + 1, ce+1);
	}
	
}

void PoissonSolver2D::generate_Pb_Tb_2D(const std::array<double, 2>& x_range, const std::array<double, 2>& y_range, const std::array<double, 2>& h_partition, BasisType trial_basis_type, BasisType test_basis_type)
{
	double left = x_range[0];
	double right = x_range[1];
	double bottom = y_range[0];
	double top = y_range[1];
	
	
	if (test_basis_type == BasisType::_2D_linear && trial_basis_type == BasisType::_2D_linear) {

		Pb_test = P_partition;
		Tb_test = T_partition;
		Pb_trial = P_partition;
		Tb_trial = T_partition;

	}
	else if (test_basis_type == BasisType::_2D_Lagrange_quadratic && trial_basis_type == BasisType::_2D_Lagrange_quadratic) {
		
		double h1 = h_basis[0];
		double h2 = h_basis[1];
		int tnp = (N_basis[0] + 1) * (N_basis[1] + 1);

		Pb_test = Eigen::MatrixXd::Zero(2, tnp);
		// Tb_test = Eigen::MatrixXd::Zero(3, N_basis); Tb should be the size of N_basis or N_partition?
		Tb_test = Eigen::MatrixXd::Zero(6, 2*N_partition[0]* N_partition[1]);

		Pb_trial = Eigen::MatrixXd::Zero(2, tnp);;
		Tb_trial = Eigen::MatrixXd::Zero(6, 2*N_partition[0] * N_partition[1]);

		for (int j = 0; j < tnp; j++) {
			int rn = j % (N_basis[1] + 1);
			int cn = j / (N_basis[1] + 1);

			Pb_test(0, j) = left + cn * h1;
			Pb_test(1, j) = bottom + rn * h2;

			Pb_trial(0, j) = left + cn * h1;
			Pb_trial(1, j) = bottom + rn * h2;
		}

		Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(N_basis[0] + 1, N_basis[1] + 1);
		for (int i = 0; i < Q.rows(); i++) {
			for (int j = 0; j < Q.cols(); j++) {
				Q(i, j) = j * (N_basis[1] + 1) + i;
			}
		}
		
		for (int n = 0; n < (N_partition[0] * N_partition[1]); n++) {

			int re = n % N_partition[1];
			int ce = n / N_partition[1];

			Tb_test(0, 2 * n) = Q(2 * re, 2 * ce);
			Tb_test(1, 2 * n) = Q(2 * re, 2 * ce + 2);
			Tb_test(2, 2 * n) = Q(2 * re + 2, 2 * ce);
			Tb_test(3, 2 * n) = Q(2 * re, 2 * ce + 1);
			Tb_test(4, 2 * n) = Q(2 * re + 1, 2 * ce + 1);
			Tb_test(5, 2 * n) = Q(2 * re + 1, 2 * ce);


			Tb_test(0, 2 * n + 1) = Q(2 * re + 2, 2 * ce);
			Tb_test(1, 2 * n + 1) = Q(2 * re, 2 * ce + 2);
			Tb_test(2, 2 * n + 1) = Q(2 * re + 2, 2 * ce + 2);
			Tb_test(3, 2 * n + 1) = Q(2 * re + 1, 2 * ce + 1);
			Tb_test(4, 2 * n + 1) = Q(2 * re + 1, 2 * ce + 2);
			Tb_test(5, 2 * n + 1) = Q(2 * re + 2, 2 * ce + 1);
		}
		Tb_trial = Tb_test;
	}
}
// suppose all Dirichlet edges
void PoissonSolver2D::generate_boundary_edges(const std::array<int, 2>& n_partition)
{
	//Information matrix for boundary edges. It uses the index of partition, not the index of FE.
	int nbe = (n_partition[0] + n_partition[1]) * 2;
	
	
	
	for (int k = 0; k < nbe; k++) {
		int mesh = -1;
		int node1 = -1;
		int node2 = -1;
		//	bottom edge
		// row_idx +  col idx * 2 * each mesh per col
		if (k < n_partition[0]) {
			mesh = 0+ k * 2 * n_partition[0];
			node1 = T_partition(0, mesh);
			node2 = T_partition(1, mesh);
		}
		//	right edge
		else if (k >= n_partition[0] && k < n_partition[0] + n_partition[1]) {
		
			// 2*row_idx+1 +  col idx * 2 * each mesh per col
			mesh = 2* (k- n_partition[0]) +1 + 2*(N_partition[0]-1)*N_partition[1];
			node1 = T_partition(1, mesh);
			node2 = T_partition(2, mesh);
		}
		//	top edge
		else if (k >= N_partition[0] + N_partition[1] && k < 2 * N_partition[0] + N_partition[1]) {
			
			// 2*row_idx+1 +  col idx * 2 * each mesh per col
			mesh = 2 * (N_partition[1] - 1) + 1 + (2 * N_partition[0] + N_partition[1] - 1 - k) * 2 * N_partition[1];
			node1 = T_partition(2, mesh);
			node2 = T_partition(0, mesh);
		}
		//	left edge
		else {			
			// 2*row_idx +  col idx * 2 * each mesh per col
			mesh = 2 * (2 * N_partition[1] + 2 * N_partition[0] - 1 - k) + 2 * 0 * N_partition[1];
			node1 = T_partition(2, mesh);
			node2 = T_partition(0, mesh);
		}
		
		//BoundaryEdge edge(BoundaryConditionType::Dirichlet, mesh, node1, node2);
		boundary_edges->push_back(BoundaryEdge(BoundaryConditionType::Dirichlet, mesh, node1, node2));
	}
	
}

void PoissonSolver2D::generate_boundary_nodes(const std::array<int, 2>& n_basis)
{
	int nbn = 2 * (n_basis[0] + n_basis[1]);
	int node = -1;

	for (int k = 0; k < nbn; k++) {
		// bottom
		if (k < n_basis[0]) {
			// row_index + col_index* #nodes per col
			node = 0 + (k) * (n_basis[1] + 1);
		}
		else if (k >= n_basis[0] && k < n_basis[0] + n_basis[1]) {
			node = (k - n_basis[0]) + n_basis[0] * (n_basis[1] + 1);
		}
		else if (k >= n_basis[0] + n_basis[1] && k < 2 * n_basis[0] + n_basis[1]) {
			node = n_basis[1] + (2 * n_basis[0] + n_basis[1] - k) * (n_basis[1] + 1);
		}
		else {
			node = (2 * (n_basis[0] + n_basis[1]) - k) + 0 * (n_basis[1] + 1);
		}
		boundary_nodes->push_back(BoundaryNode(BoundaryConditionType::Dirichlet, node));
	}
}

double PoissonSolver2D::triangular_reference_basis(double x, double y, BasisType basis_type, int basis_index, int derivative_degree_x, int derivative_degree_y)
{
	return 0.0;
}
