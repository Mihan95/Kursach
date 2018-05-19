#include "stdafx.h"
#include "Solver.h"
#include <array>
#include <memory>

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

static void add_figure_point(const double alpha1, const double alpha2, const double r1, const double r2, const double poly_edge_angle, double &alpha, std::vector<double> &vec)
{
	double r = r1 * r2 * sin(alpha1 - alpha2) / (r1 * sin(alpha - alpha1) - r2 * sin(alpha - alpha2));
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
		add_figure_point(alpha1, alpha2, 1., 1., poly_edge_angle, alpha, poly);
	}

}

void Polygon::project_to_mesh(uint32_t rot)
{
	const double poly_edge_angle = 2. * M_PI / n_vertices;

	double alpha1 = 0.;
	double alpha2 = poly_edge_angle;
	double alpha = 0.;

	std::rotate(poly.rbegin(), poly.rbegin() + rot, poly.rend());

	double r1 = poly[0];
	double r2 = poly[1];

	uint16_t poly_edge_count;
	for (uint32_t i = 0; i < mesh_size; i++)
	{
		if (alpha > alpha2)
		{
			alpha1 = alpha2;
			alpha2 += poly_edge_angle;
			r1 = r2;
			r2 = poly[poly_edge_count++];
		}
		add_figure_point(alpha1, alpha2, r1, r2, step, alpha, on_mesh);
	}

}

Solver::Solver(func &poly_pare)
{
	tasks.push_back(poly_pare);
}

Solver::~Solver()
{
	tasks.clear();
}

void Solver::add_task(func &poly_pare)
{
	tasks.push_back(poly_pare);
}

void Solver::solve()
{
	for (auto task : tasks)
	{
		const uint32_t N = task.f.get_mesh_size();

		double *f_in  = (double* )fftw_malloc(sizeof(double)*N);
		double *f_out = (double* )fftw_malloc(sizeof(double)*N);
		double *g_in  = (double* )fftw_malloc(sizeof(double)*N);
		double *g_out = (double* )fftw_malloc(sizeof(double)*N);

		fftw_plan f_my_plan, g_my_plan;
		f_my_plan = fftw_plan_r2r_1d(N, f_in, f_out, FFTW_R2HC, FFTW_ESTIMATE);
		g_my_plan = fftw_plan_r2r_1d(N, g_in, g_out, FFTW_R2HC, FFTW_ESTIMATE);

		std::copy(task.f.get_poly_on_mesh().begin(), task.f.get_poly_on_mesh().end(),  f_in);
		std::copy(task.f.get_poly_on_mesh().begin(), task.f.get_poly_on_mesh().end(), f_out);

		std::copy(task.g.get_poly_on_mesh().begin(), task.g.get_poly_on_mesh().end(), g_in);
		std::copy(task.g.get_poly_on_mesh().begin(), task.g.get_poly_on_mesh().end(), g_out);
		
		fftw_execute(f_my_plan);
		fftw_execute(g_my_plan);

		f_my_plan = fftw_plan_r2r_1d(N, f_in, f_out, FFTW_HC2R, FFTW_ESTIMATE);

		/*сгруппируем полученные коэффициенты*/
		memset(f_in, 0, N * sizeof(double));
		for (uint32_t i = 0; i < N; i++)
		{
			f_in[0] += f_out[i] * f_out[i];
			f_in[0] += g_out[i] * g_out[i];
		}
		f_in[0] -= 2. * (f_out[0] * g_out[0]);
		f_in[0] *= 8. * M_PI / N;

		for (uint32_t i = 1; i < N/2; i++)
		{
			f_in[i] = f_out[i] * g_out[i] + f_out[N - i] * g_out[N - i];
			f_in[i] *= -16. * M_PI / N;
		}
		f_in[N / 2] = f_out[N/2] * g_out[N/2];
		f_in[N / 2] *= -16. * M_PI / N;

		for (uint32_t i = 1; i < N / 2; i++)
		{
			f_in[N - i] = f_out[i] * g_out[N-i] - f_out[N - i] * g_out[i];
			f_in[N - i] *= -16. * M_PI / N;
		}

		uint32_t rot_angle = std::min_element(f_out, f_out + N) - f_out;

		/*найдем норму наложения с новым сдвигом*/


		fftw_execute(f_my_plan);

		fftw_destroy_plan(f_my_plan);
		fftw_destroy_plan(g_my_plan);
		fftw_free(f_in);
		fftw_free(f_out);
		fftw_free(g_in);
		fftw_free(g_out);
	}
}