#pragma once
#include "stdafx.h"
#include <vector>

class Polygon
{
public:
	friend static void add_figure_point(const double alpha1, const double alpha2, const double poly_edge_angle, double &alpha, std::vector<double> &vec);
	Polygon();
	Polygon(uint32_t NVerts, uint32_t NMesh);
	~Polygon();
private:
	uint32_t n_vertices;
	uint32_t mesh_size;
	double step;
	std::vector<double> poly;
	std::vector<double> on_mesh;
	void create_from_regular_poly(uint16_t n_angle);
	void project_to_mesh();
};

class Solver
{
public:
	Solver() { }
	~Solver() { }

private:

};