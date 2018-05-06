#pragma once
#include "stdafx.h"
#include <vector>

class Polygon
{
public:
	friend static void add_figure_point(const double alpha1, const double alpha2, const double poly_edge_angle, double &alpha, Polygon &pol);
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
};

class Solver
{
public:
	Solver() : mesh_size(1), step(2. * M_PI / mesh_size){ f.resize(mesh_size); }
	Solver(uint32_t N) : mesh_size(N), step(2. * M_PI / mesh_size) { f.resize(mesh_size); }
	~Solver() { f.clear(); };

	void fill_polygon(int n_angles);
private:
	uint32_t mesh_size;
	double step;
	std::vector<double> f;
	std::vector<double> g;
};