#pragma once
#include "stdafx.h"
#include <vector>

class Polygon
{
public:
	friend static void add_figure_point(const double alpha1, const double alpha2, const double r1, const double r2, const double poly_edge_angle, double &alpha, std::vector<double> &vec);
	Polygon();
	Polygon(uint32_t NVerts, uint32_t NMesh);
	~Polygon();

	const uint32_t get_mesh_size() { return mesh_size; };
	std::vector<double> get_poly_on_mesh() { return on_mesh; }
	void create_from_regular_poly(uint16_t n_angle);
	void project_to_mesh(uint32_t rot);
private:
	uint32_t n_vertices;
	uint32_t mesh_size;
	double step;
	std::vector<double> poly;
	std::vector<double> on_mesh;
};

struct func
{
	Polygon f;
	Polygon g;
};

class Solver
{
public:
	Solver(func &poly_pare);
	~Solver();
	void add_task(func &poly_pare);
	void solve();
private:
	std::vector<func> tasks;
};