#include "stdafx.h"
#include "Solver.h"
#include <array>
#include <memory>

#define USE_FFTW

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

void Solver::solve1()
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

static void group_coeffs(uint32_t N, double *h_out, double *g_out, double *h_in)
{
	const double  mul =  2.;
	const double _mul = -2.;
	const int border = N / 2 + 1;

	std::transform(h_out, h_out + border, h_out, std::bind2nd(std::multiplies<double>(), mul));
	std::transform(g_out, g_out + border, g_out, std::bind2nd(std::multiplies<double>(), mul));

	std::transform(h_out + border, h_out + N, h_out + border, std::bind2nd(std::multiplies<double>(), _mul));
	std::transform(g_out + border, g_out + N, g_out + border, std::bind2nd(std::multiplies<double>(), _mul));

	h_in[0] = h_out[0] * g_out[0];
	for (uint32_t i = 1; i < N / 2; i++)
	{
		h_in[i] = h_out[i] * g_out[i] + h_out[N - i] * g_out[N - i];
	}
	h_in[N / 2] = h_out[N / 2] * g_out[N / 2];

	for (uint32_t i = 1; i < N / 2; i++)
	{
		h_in[N - i] = h_out[N - i] * g_out[i] - h_out[i] * g_out[N - i];
	}
}

static void free(fftw_plan &plan1, fftw_plan &plan2, fftw_plan &plan3, double *buff1, double *buff2, double *buff3, double *buff4, double *buff5)
{
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
	fftw_destroy_plan(plan3);
	fftw_free(buff1);
	fftw_free(buff2);
	fftw_free(buff3);
	fftw_free(buff4);
	fftw_free(buff5);
}

void Solver::solve2()
{
	for (auto task : tasks)
	{
		const uint32_t N = task.f.get_mesh_size();

		double *h_in  = (double*)fftw_malloc(sizeof(double)*N);
		double *h_out = (double*)fftw_malloc(sizeof(double)*N);
		double *g_in  = (double*)fftw_malloc(sizeof(double)*N);
		double *g_out = (double*)fftw_malloc(sizeof(double)*N);
		double *w_out = (double*)fftw_malloc(sizeof(double)*N);

		fftw_plan my_plan_j, my_plan_s, my_plan;
		my_plan_j = fftw_plan_r2r_1d(N, h_in, h_out, FFTW_R2HC, FFTW_ESTIMATE);
		my_plan_s = fftw_plan_r2r_1d(N, g_in, g_out, FFTW_R2HC, FFTW_ESTIMATE);

		for (uint32_t i = 0; i < N; i++)
		{
			h_in[i] = task.f.get_poly_on_mesh()[i] * task.w[i];
			g_in[i] = task.g.get_poly_on_mesh()[i];
		}

		fftw_execute(my_plan_j);
		fftw_execute(my_plan_s);
		
		group_coeffs(N, h_out, g_out, h_in);

		double *e_in  = g_in;   g_in = nullptr;
		double *e_out = g_out; g_out = nullptr;
		double *w_in  = h_out; h_out = nullptr;

		my_plan_j = fftw_plan_r2r_1d(N, w_in, w_out, FFTW_R2HC, FFTW_ESTIMATE);
		my_plan_s = fftw_plan_r2r_1d(N, e_in, e_out, FFTW_R2HC, FFTW_ESTIMATE);

		for (int i = 0; i < N; i++)
		{
			w_in[i] = task.w[i];
			e_in[i] = task.g.get_poly_on_mesh()[i] * task.g.get_poly_on_mesh()[i];
		}

		fftw_execute(my_plan_j);
		fftw_execute(my_plan_s);

		double *ew_in = w_in;  w_in = nullptr;
		double *ew_out = e_in; e_in = nullptr;

		my_plan = fftw_plan_r2r_1d(N, ew_in, ew_out, FFTW_HC2R, FFTW_ESTIMATE);

		/*сгруппируем полученные коэффициенты*/
		group_coeffs(N, w_out, e_out, ew_in);
		const double N_1 = 1. / N;
		for (int i = 0; i < N; i++)
		{
			ew_in[i] -= 2.* h_in[i];
		}
		std::for_each(ew_in, ew_in + N, [&](double &n) { n *= (2. * N_1 * N_1 * N_1); });

#ifdef USE_FFTW
		fftw_execute(my_plan);
#else
		memset(w_out, 0, N * sizeof(double));
		for (int i = 0; i < N; i++)
		{
			for (int k = 0; k <= N / 2; k++)
			{
				w_out[i] += (ew_in[k] * cos(2.*M_PI * k * i / N));
			}
			for (int k = N - 1; k >= N / 2 + 1; k--)
			{
				w_out[i] += (ew_in[k] * sin(2.*M_PI * k * i / N));
			}
		}

#endif

		double *min_ptr = std::min_element(ew_out, ew_out + N);
		uint64_t rot_angle = min_ptr - ew_out;
		std::cout << "Angle = " << rot_angle << std::endl;

		g_in = ew_in;
		for (int i = 0; i < N; i++)
		{
			g_in[i] = task.g.get_poly_on_mesh()[i];
		}
		std::rotate(g_in, g_in + rot_angle, g_in + N);
		double norm = 0;
		for (int i = 0; i < N; i++)
		{
			norm += ((task.f.get_poly_on_mesh()[i] - g_in[i]) * (task.f.get_poly_on_mesh()[i] - g_in[i]) * task.w[i]);
		}
		std::cout << "Difference = " << norm << std::endl;

		free(my_plan, my_plan_s, my_plan_j, ew_in, ew_out, e_out, h_in, w_out);

	}
}

bool Solver::is_refinement()
{
	for (auto task : tasks)
	{
		const double g_angle = 2. * M_PI / task.g.get_mesh_size();
		const double f_angle = 2. * M_PI / task.f.get_mesh_size();
		double alpha1 = g_angle;
		double alpha  = f_angle;
		double sum_f   = 0.;
		double sum_g   = 0.;
		double sum_g_1 = 0.;

		for (auto f_i : task.f.get_poly_on_mesh())
		{
			if (alpha < alpha1)
			{
				double g_i   = sin(g_angle) / (sin((alpha - f_angle) - (alpha1 - g_angle)) + sin(alpha1 - (alpha - f_angle)));
				double g_i_1 = sin(g_angle) / (sin(alpha - (alpha1 - g_angle)) + sin(alpha1 - alpha));
				
				double diff_g_i = g_i_1 - g_i;
				sum_g   += g_i * diff_g_i;
				sum_g_1 += g_i_1 * diff_g_i;
				sum_f   += f_i * diff_g_i;
				alpha += f_angle;
			}
			else
			{
				alpha1 += g_angle;
				alpha  += f_angle;
			}
		}
		double t = (sum_f - sum_g_1) / (sum_g - sum_g_1);

		/*посчитаем нормы*/
		alpha1 = g_angle;
		alpha  = f_angle;
		double norm_before = 0.;
		double norm_after  = 0.;
		for (auto f_i : task.f.get_poly_on_mesh())
		{
			if (alpha < alpha1)
			{
				double g_i = sin(g_angle) / (sin((alpha - f_angle) - (alpha1 - g_angle)) + sin(alpha1 - (alpha - f_angle)));
				double g_i_1 = sin(g_angle) / (sin(alpha - (alpha1 - g_angle)) + sin(alpha1 - alpha));
				norm_before += (f_i - g_i) * (f_i - g_i);
				norm_after  += (f_i - t * g_i - (1. - t) * g_i_1) * (f_i - t * g_i - (1. - t) * g_i_1);
				alpha += f_angle;
			}
			else
			{
				alpha1 += g_angle;
				alpha  += f_angle;
			}
		}

		norm_before = sqrt(norm_before);
		norm_after  = sqrt(norm_after);
		norm_before /= task.f.get_mesh_size();
		norm_after  /= task.f.get_mesh_size();
		std::cout << "Norm before: " << norm_before << std::endl;
		std::cout << "Norm after: "  << norm_after  << std::endl;
	}
	return true;

}