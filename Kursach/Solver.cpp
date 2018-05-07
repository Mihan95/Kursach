#include "stdafx.h"
#include "Solver.h"

Polygon::Polygon()
{
	n_vertices = mesh_size = 1;
	step = 2. * M_PI;

	poly.reserve(n_vertices);
	on_mesh.reserve(mesh_size);
}

Polygon::Polygon(uint32_t NVerts, uint32_t NMesh)
{
	n_vertices = NVerts;
	mesh_size = NMesh;
	step = 2. * M_PI / mesh_size;
	
	poly.reserve(n_vertices);
	on_mesh.reserve(mesh_size);
}

Polygon::~Polygon()
{
	poly.clear();
	on_mesh.clear();
}

static void add_figure_point(const double alpha1, const double alpha2, const double poly_edge_angle, double &alpha, std::vector<double> &vec)
{
	double r = sin(alpha1 - alpha2) / (sin(alpha - alpha1) - sin(alpha - alpha2));
	vec.push_back(r);
	alpha += poly_edge_angle;
}

void Polygon::create_from_regular_poly(uint16_t n_angle)
{
	const double edge_angle = 2. * M_PI / n_angle;
	const double poly_edge_angle = 2. * M_PI / n_vertices;

	double alpha1 = 0.;
	double alpha2 = edge_angle;
	double alpha = poly_edge_angle;

	poly.push_back(1.);
	for (uint32_t i = 1; i < n_vertices; i++)
	{
		if (alpha > alpha2)
		{
			alpha1 = alpha2;
			alpha2 += edge_angle;
		}
		add_figure_point(alpha1, alpha2, poly_edge_angle, alpha, poly);
	}

}

void Polygon::project_to_mesh()
{
	double alpha1 = 0;
	double alpha2 = poly[0];
	double alpha = step;

	on_mesh.push_back(1.);
	uint16_t poly_edge = 1;
	for (uint32_t i = 1; i < mesh_size; i++)
	{
		if (alpha > alpha2)
		{
			alpha1 = alpha2;
			alpha2 = poly[poly_edge++];
		}
		add_figure_point(alpha1, alpha2, step, alpha, on_mesh);
	}

}

void Solver::fill_polygon(int n_angles)
{

}