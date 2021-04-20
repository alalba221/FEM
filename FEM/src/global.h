#pragma once
enum class  BoundaryConditionType {
	Dirichlet,
	Neumann,
	Robins
};
enum class BasisType {
	_1D_linear,
	_1D_quadratic,
	_2D_linear,
	_2D_Lagrange_quadratic
};

struct BoundaryNode {
	BoundaryConditionType bc_type;
	int global_index = -1;
	BoundaryNode(BoundaryConditionType type, int index) :
		bc_type(type), global_index(index) {
	}
};

struct BoundaryEdge {
	BoundaryConditionType bc_type;
	int mesh_index = -1;
	int first_node_global_index = -1;
	int second_node_global_index = -1;
	BoundaryEdge(BoundaryConditionType type, int index, int node1, int node2) :
		bc_type(type), mesh_index(index), first_node_global_index(node1),second_node_global_index(node2) {
	}
};
