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
	double r = r1 * r2 * sin(alpha2 - alpha1) / (r1 * sin(alpha - alpha1) + r2 * sin(alpha2 - alpha));
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

	double r1 = poly[0];
	double r2 = poly[1];

	uint16_t poly_edge_count = 0;
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

	std::rotate(on_mesh.rbegin(), on_mesh.rbegin() + rot, on_mesh.rend());

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
		double *temp  = (double*)fftw_malloc(sizeof(double)*N);

		fftw_plan f_my_plan, g_my_plan, my_plan;
		f_my_plan = fftw_plan_r2r_1d(N, f_in, f_out, FFTW_R2HC, FFTW_ESTIMATE);
		g_my_plan = fftw_plan_r2r_1d(N, g_in, g_out, FFTW_R2HC, FFTW_ESTIMATE);

		for (int i = 0; i < N; i++)
		{
			f_in[i] = task.f.get_poly_on_mesh()[i];
			g_in[i] = task.g.get_poly_on_mesh()[i];
		}
		
		fftw_execute(f_my_plan);
		fftw_execute(g_my_plan);

		for (int i = 0; i < N; i++)
		{
			temp[i] = f_out[i];
		}

		my_plan = fftw_plan_r2r_1d(N, f_in, f_out, FFTW_HC2R, FFTW_ESTIMATE);
		memset(f_in, 0, N * sizeof(double));

		for (int i = 0; i < N; i++)
		{
			f_out[i] = temp[i];
		}
		/*сгруппируем полученные коэффициенты*/
		for (int i = 0; i <= N / 2; i++)
		{
			f_out[i] *= 2.;
			g_out[i] *= 2.;
		}
		for (int i = N / 2 + 1; i < N; i++)
		{
			f_out[i] *= -2.;
			g_out[i] *= -2.;
		}

		for (uint32_t i = 0; i < N; i++)
		{
			f_in[0] += f_out[i] * f_out[i];
			f_in[0] += g_out[i] * g_out[i];
		}
		f_in[0] -= 2. * (f_out[0] * g_out[0]);
		f_in[0] *= 2. * M_PI / N;
		f_in[0] /= (N*0.5);

		for (uint32_t i = 1; i < N / 2; i++)
		{
			f_in[i] = f_out[i] * g_out[i] + f_out[N - i] * g_out[N - i];
			f_in[i] *= -4. * M_PI / N;
			f_in[i] /= (N*0.5);
		}
		f_in[N / 2] = f_out[N / 2] * g_out[N / 2];
		f_in[N / 2] *= -4. * M_PI / N;
		f_in[N / 2] /= (N*0.5);

		for (uint32_t i = 1; i < N / 2; i++)
		{
			f_in[N - i] = f_out[i] * g_out[N - i] - f_out[N - i] * g_out[i];
			f_in[N - i] *= 4. * M_PI / N;
			f_in[N - i] /= (N*0.5);
		}
		/*for (int i = 0; i <= N / 2; i++)
		{
			f_in[i] /= 2.;
		}
		for (int i = N / 2 + 1; i < N; i++)
		{
			f_in[i] /= -2.;
		}*/

		memset(f_out, 0, N * sizeof(double));		
		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k <= N / 2; k++)
			{
				f_out[i] += (f_in[k] * cos(2.*M_PI * k * i/N));
			}
			for (int k = N-1; k >= N / 2+1; k--)
			{
				f_out[i] += (f_in[k] * sin(2.*M_PI * k * i/N));
			}
		}

		//fftw_execute(my_plan);

		double *min_ptr =  std::min_element(f_out, f_out + N);
		uint64_t rot_angle = min_ptr - f_out;

		std::cout << "Angle = " << rot_angle << std::endl;
		std::cout << "Difference = " << *min_ptr << std::endl;

		fftw_destroy_plan(f_my_plan);
		fftw_destroy_plan(g_my_plan);
		fftw_destroy_plan(my_plan);
		fftw_free(f_in);
		fftw_free(f_out);
		fftw_free(g_in);
		fftw_free(g_out);
	}
}