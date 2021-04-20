#include "PoissonSolver1D.h"

PoissonSolver1D::PoissonSolver1D(double h_partition, BasisType trial_basis_type, BasisType  test_basis_type)
	:h_partition(h_partition),trial_basis_type(trial_basis_type),test_basis_type(test_basis_type) {
	
	boundary_nodes = new std::vector<BoundaryNode>;
	gs = new GaussInteger(4);


	if (trial_basis_type == BasisType::_1D_linear) {
		number_of_local_basis_trial = 2;
	}
	else if (trial_basis_type == BasisType::_1D_quadratic) {
		number_of_local_basis_trial = 3;
		h_fe = h_partition / 2;
		
	}
	if (test_basis_type == BasisType::_1D_linear) {
		number_of_local_basis_test = 2;
	}
	else if (test_basis_type == BasisType::_1D_quadratic) {
		number_of_local_basis_test = 3;
		h_fe = h_partition / 2;
	}

};
void PoissonSolver1D::generate_P_T_1D(double left, double right, double h_partition, BasisType basis_type)
{
	double h = h_partition;
	if (basis_type == BasisType::_1D_linear) {
		
		N_partition = (right - left) / h;
		int N = (right - left) / h;
		P = Eigen::MatrixXd::Zero(1, N_partition + 1);
		T = Eigen::MatrixXd::Zero(2, N_partition);

		for (int i = 0; i < P.cols(); i++) {
			P(0, i) = left + i * h;
		}

		for (int i = 0; i < T.cols(); i++) {
			T(0, i) = i;
			T(1, i) = i + 1;
		}
	}	
}

void PoissonSolver1D::generate_Pb_Tb_1D(double left, double right, double h_partition, BasisType trial_basis_type, BasisType test_basis_type)
{
	if (test_basis_type == BasisType::_1D_linear && trial_basis_type==BasisType::_1D_linear) {
		
		N_partition = (right - left) / h_partition;
		N_fe = N_partition;
		
		Pb_test = P;
		Tb_test = T;
		Pb_trial = P;
		Tb_trial = T;
		
	}
	else if (test_basis_type == BasisType::_1D_quadratic && trial_basis_type == BasisType::_1D_quadratic) {
		N_partition = (right - left) / h_partition;
		N_fe = N_partition*2;
		
		Pb_test = Eigen::MatrixXd::Zero(1, N_fe + 1);
		// Tb_test = Eigen::MatrixXd::Zero(3, N_fe); Tb should be the size of N_fe or N_partition?
		Tb_test = Eigen::MatrixXd::Zero(3, N_partition);
		
		Pb_trial = Eigen::MatrixXd::Zero(1, N_fe + 1);
		// Tb_trial = Eigen::MatrixXd::Zero(3, N_fe);
		Tb_trial = Eigen::MatrixXd::Zero(3, N_partition);

		for (int i = 0; i < N_fe+1; i++) {
			Pb_test(0, i) = left + i * h_fe;
			Pb_trial(0, i) = left + i * h_fe;
		}
		// left right  mid
		for (int i = 0; i < N_partition; i++) {
			Tb_test(0, i) = 2*i;
			Tb_test(1, i) = 2*i + 2;
			Tb_test(2, i) = 2*i + 1;
			
			Tb_trial(0, i) = 2*i;
			Tb_trial(1, i) = 2*i + 2;
			Tb_trial(2, i) = 2*i + 1;			
		}
	}
}
// vertices stoire the bound on the mesh not FE!
double PoissonSolver1D::local_basis_1D(double x, const std::vector<double>& vertices, BasisType basis_type, int local_basis_index, int derivative_degree)
{
	double result =0;
	if (basis_type == BasisType::_1D_linear) {
		if (derivative_degree == 0) {
			if (local_basis_index == 0) {
				result = (vertices[1] - x) / (vertices[1] - vertices[0]);
			}
			else if (local_basis_index == 1) {
				result = (x - vertices[0]) / (vertices[1] - vertices[0]);
			}
		}
		else if (derivative_degree == 1) {
			if (local_basis_index == 0) {
				result = 1.0 / (vertices[0] - vertices[1]);
			}
			else if (local_basis_index == 1) {
				result = 1.0 / (vertices[1] - vertices[0]);
			}
		}
		else if (derivative_degree >= 2) {
			result = 0;
		}
	}
	else if (basis_type == BasisType::_1D_quadratic) {
		
		double x_bar = (x - vertices[0]) / (vertices[1] - vertices[0]);
		double dx_bar = 1.0 / h_partition;
		if (derivative_degree == 0) {
			if (local_basis_index == 0) {
				result = 2 * x_bar * x_bar - 3 * x_bar + 1;
			}
			else if (local_basis_index == 1) {
				result = 2 * x_bar * x_bar - x_bar;
			}
			else if (local_basis_index == 2) {
				result = -4 * x_bar * x_bar + 4 * x_bar;
			}
		}
		else if (derivative_degree == 1) {
			if (local_basis_index == 0) {
				result = 4*x_bar*dx_bar-3*dx_bar;
			}
			else if (local_basis_index == 1) {
				result = 4 * x_bar * dx_bar - dx_bar;
			}
			else if (local_basis_index == 2) {
				result = -8 * x_bar * dx_bar + 4 * dx_bar;
			}
		}
		else if (derivative_degree == 2) {
			if (local_basis_index == 0) {
				result = 4 * dx_bar * dx_bar;
			}
			else if (local_basis_index == 1) {
				result = 4 * dx_bar * dx_bar;
			}
			else if (local_basis_index == 2) {
				result = -8 * dx_bar * dx_bar;
			}	
		}
		else if (derivative_degree >= 3) {
			result = 0;
		}		
	}
	return result;
}



void PoissonSolver1D::assemble_matrix_from_1D_integral(int number_of_partition, int test_derivative_degree, int trial_derivative_degree,
	int number_of_trial_local_basis,
	int number_of_test_local_basis)
{
	// A should be the size of # of test rows x # of trial cols
	A = Eigen::SparseMatrix<double>(N_fe + 1, N_fe + 1);
	for (int n = 0; n < number_of_partition; n++) {
		std::vector<double>vertices(2);
		vertices[0] = P(0, T(0, n));
		vertices[1] = P(0, T(1, n));
		double lowerbound = std::min(vertices[0], vertices[1]);
		double upperbound = std::max(vertices[0], vertices[1]);
		vertices[0] = lowerbound;
		vertices[1] = upperbound;
		
		
		gs->set_bound(lowerbound, upperbound);
		for (int alpha = 0; alpha < number_of_trial_local_basis; alpha++) {
			for (int beta = 0; beta < number_of_test_local_basis; beta++) {
				
				auto trial_basis = [this, &vertices, alpha, trial_derivative_degree](double x) {return local_basis_1D(x, vertices, trial_basis_type, alpha, trial_derivative_degree); };
				auto test_basis = [this, &vertices, beta, test_derivative_degree](double x) {return local_basis_1D(x, vertices, test_basis_type, beta, test_derivative_degree); };
				auto finalFunc = [this, &trial_basis, &test_basis](double x)->double {return coefficient_function(x) * trial_basis(x) * test_basis(x); };
				double temp = gs->integral(finalFunc);
				
				/*if (Tb_test(beta, n) == Tb_trial(alpha, n)) {
					A.coeffRef(Tb_test(beta, n), Tb_trial(alpha, n)) += temp;
				}
				else
					A.coeffRef(Tb_test(beta, n), Tb_trial(alpha, n)) = temp;*/
				// mesh 边界需要+=， 其他元素直接=, 因为mesh边界点对应的基函数位于两个相邻的mesh上
				A.coeffRef(Tb_test(beta, n), Tb_trial(alpha, n)) += temp;
			}
		}
	}
}

void PoissonSolver1D::assemble_vector_from_1D_integral(int number_of_partition, int test_derivative_degree
	, int number_of_test_local_basis)
{
	b = Eigen::VectorXd::Zero(N_fe+1);

	for (int n = 0; n < number_of_partition; n++) {
		std::vector<double>vertices(2);
		vertices[0] = P(0, T(0, n));
		vertices[1] = P(0, T(1, n));
		double lowerbound = std::min(vertices[0], vertices[1]);
		double upperbound = std::max(vertices[0], vertices[1]);
		vertices[0] = lowerbound;
		vertices[1] = upperbound;

		gs->set_bound(lowerbound, upperbound);
		for (int beta = 0; beta < number_of_test_local_basis; beta++) {
	
			auto test_basis = [this, &vertices, beta, test_derivative_degree](double x) {return local_basis_1D(x, vertices, test_basis_type, beta, test_derivative_degree); };
			auto finalFunc = [this, &test_basis](double x)->double {return function_f(x)  * test_basis(x); };
			double temp = gs->integral(finalFunc);
			b(Tb_test(beta, n)) += temp;
		}
		
	}
}
// It uses the index of FE, not the index of partition.
void PoissonSolver1D::generate_boundary_nodes_1D()
{
	boundary_nodes->push_back(BoundaryNode(BoundaryConditionType::Dirichlet, 0));
	boundary_nodes->push_back(BoundaryNode(BoundaryConditionType::Dirichlet, N_fe));
}

void PoissonSolver1D::treat_Dirichlet_boundary_1D()
{
	for (int i = 0; i < boundary_nodes->size(); i++) {
		if (boundary_nodes->at(i).bc_type == BoundaryConditionType::Dirichlet) {
			int index = boundary_nodes->at(i).global_index;
			A.row(index)*=0;
			A.coeffRef(index, index) = 1;
			b(index) = function_g(Pb_test(0, index));
		}
	}
}
