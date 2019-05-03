#include "stdafx.h"
#include "Solver.h"
#include <array>
#include <memory>
#include "poly34.h"
#include "gsl_min.h"
#include <gsl/gsl_errno.h>
#include "gsl_multimin.h"
#include <fstream>

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

static double shift_point(double r1, double r2, double cos_, double t)
{
	double t_1 = t - 1.;
	return sqrt(t_1*t_1 * r1*r1 + t * t*r2*r2 - 2.*t*t_1*r1*r2*cos_);
}

static double shift_point2(const double r1, const double r2, const double t, const uint32_t N)
{
	const double M_2PI = 2. * M_PI;
	return (r1*r2*sin(M_2PI / N) / (r2 * sin(M_2PI * (1 - t) / N) + r1 * sin(M_2PI * t / N)));
}

void Polygon::shift_poly_on_mesh(double t)
{
	const uint32_t N = on_mesh.size();

	std::vector<double> new_on_mesh(N);
	//new_on_mesh.resize(N);
	{
		new_on_mesh[0] = shift_point2(on_mesh[0], on_mesh[N-1], t, N);
	}
	for (uint32_t i = 1; i < N; i++)
	{
		new_on_mesh[i] = shift_point2(on_mesh[i], on_mesh[i - 1], t, N);
	}
	on_mesh = new_on_mesh;
}

//void Polygon::shift_poly_on_mesh_op(double t)
//{
//	const uint32_t N = on_mesh.size();
//	const double cos_ = cos(2. * M_PI / N);
//
//	std::vector<double> tmp_on_mesh(N);
//	std::vector<double> diff_on_mesh(N);
//
//	for (uint32_t i = 0; i < N-1; i++)
//	{
//		tmp_on_mesh[i] = shift_point2(new_on_mesh[i], new_on_mesh[i + 1], t, N);
//		diff_on_mesh[i] = fabs(on_mesh[i] - tmp_on_mesh[i]);
//	}
//	{
//		tmp_on_mesh[N-1] = shift_point2(new_on_mesh[N - 1], new_on_mesh[0], t, N);
//		diff_on_mesh[N - 1] = fabs(on_mesh[N - 1] - tmp_on_mesh[N - 1]);
//	}
//	on_mesh = new_on_mesh;
//}

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

void Polygon::read_polygon_on_mesh(std::string input_name)
{
	std::ifstream in(input_name);
	if (in.is_open())
	{
		uint32_t n;
		in >> n;
		on_mesh.clear();
		poly.clear();
		on_mesh.reserve(n);
		poly.reserve(n);
		double x, y, r;

		uint32_t i = 0;
		while (in >> x >> y)
		{
			i++;
			r = sqrt(x * x + y * y);
			on_mesh.push_back(r);
			poly.push_back(r);
		}
		if (n != i)
			std::cout << "n = " << n << " != " << i << " = i\n";
	}
	else
		std::cout << input_name.c_str() << " open error\n";
	in.close();
	//for (uint32_t i = 0; i < 10; i++)
	//	std::cout << on_mesh[i] << std::endl;
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

		for (uint32_t i = 0; i < N; i++)
		{
			f_in[i] = task.f.get_poly_on_mesh()[i];
			g_in[i] = task.g.get_poly_on_mesh()[i];
		}
		
		fftw_execute(f_my_plan);
		fftw_execute(g_my_plan);

		for (uint32_t i = 0; i < N; i++)
		{
			temp[i] = f_out[i];
		}

		my_plan = fftw_plan_r2r_1d(N, f_in, f_out, FFTW_HC2R, FFTW_ESTIMATE);
		memset(f_in, 0, N * sizeof(double));

		for (uint32_t i = 0; i < N; i++)
		{
			f_out[i] = temp[i];
		}
		/*сгруппируем полученные коэффициенты*/
		for (uint32_t i = 0; i <= N / 2; i++)
		{
			f_out[i] *= 2.;
			g_out[i] *= 2.;
		}
		for (uint32_t i = N / 2 + 1; i < N; i++)
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
		for (uint32_t i = 0; i < N; i++)
		{
			for (uint32_t k = 0; k <= N / 2; k++)
			{
				f_out[i] += (f_in[k] * cos(2.*M_PI * k * i/N));
			}
			for (uint32_t k = N-1; k >= N / 2+1; k--)
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

		for (uint32_t i = 0; i < N; i++)
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
		for (uint32_t i = 0; i < N; i++)
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
		for (uint32_t i = 0; i < N; i++)
		{
			g_in[i] = task.g.get_poly_on_mesh()[i];
		}
		std::rotate(g_in, g_in + rot_angle, g_in + N);
		double norm = 0;
		for (uint32_t i = 0; i < N; i++)
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
		const double sin_g_angle = sin(g_angle);
		double alpha1 = g_angle;
		double alpha  = f_angle;
		double sum_f   = 0.;
		double sum_g   = 0.;
		double sum_g_1 = 0.;

		std::vector<double> gi;
		std::vector<double> gi_1;

		for (auto f_i : task.f.get_poly_on_mesh())
		{
			if (alpha < alpha1)
			{  /* Для оптимизации НАЧАЛО */
				const double angle_gi1_1 = alpha - (alpha1 - g_angle);
				const double angle_gi1 = angle_gi1_1 - f_angle;

				const double angle_gi1_2 = alpha1 - alpha;
				const double angle_gi2 = angle_gi1_2 + f_angle;
				/* КОНЕЦ для оптимизации */

				gi.push_back( sin_g_angle / ( sin(angle_gi1) + sin(angle_gi2) ) );
				gi_1.push_back( sin_g_angle / ( sin(angle_gi1_1) + sin(angle_gi1_2) ) );

				double diff_g_i = gi_1.back() - gi.back();
				sum_g   += gi.back() * diff_g_i;
				sum_g_1 += gi_1.back() * diff_g_i;
				sum_f   += f_i * diff_g_i;
				alpha += f_angle;
			}
			else
			{
				alpha1 += g_angle;
				alpha  += f_angle;
			}
		}
		const double t = (sum_f - sum_g_1) / (sum_g - sum_g_1);

		/*посчитаем нормы*/
		alpha1 = g_angle;
		alpha  = f_angle;
		double norm_before = 0.;
		double norm_after  = 0.;
		for (uint32_t i = 0; i < task.f.get_mesh_size(); i++)
		{
			if (alpha < alpha1)
			{  
				const auto fi = task.f.get_poly_on_mesh()[i];
				/* Для оптимизации НАЧАЛО */
				const double fi_gi = fi - gi[i];
				const double fi_gi_t = fi - t * gi[i] - (1. - t) * gi_1[i];
				/* КОНЕЦ для оптимизации */
				
				norm_before += fi_gi * fi_gi;
				norm_after  += fi_gi_t * fi_gi_t;
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
		gi.clear();
		gi_1.clear();
	}
	return true;

}

//s2;
static double compute_angle_quad_functional(const std::vector<double> & g_coo, const std::vector<double> & f_coo, double t)
{
	const double M_2PI = 2. * M_PI;
	const uint32_t N = f_coo.size();

	double sum = 0.;
	double tmp_sum;
	for (uint32_t i = 0; i < N-1; i++)
	{
		tmp_sum = (t*t - 2.*t + 1) * g_coo[i] * g_coo[i] + t * t * g_coo[i + 1] * g_coo[i + 1]
			- 2. * (t*t - t) * g_coo[i] * g_coo[i + 1] * cos(M_2PI / N) - f_coo[i] * f_coo[i];
		sum = sum + (tmp_sum * tmp_sum);
	}
	{
		tmp_sum = (t*t - 2.*t + 1) * g_coo[N - 1] * g_coo[N - 1] + t * t * g_coo[0] * g_coo[0]
			- 2. * (t*t - t) * g_coo[N - 1] * g_coo[0] * cos(M_2PI / N) - f_coo[N - 1] * f_coo[N - 1];
		sum = sum + (tmp_sum * tmp_sum);
	}

	return sum;
}

//s1;
struct my_f_params { std::vector<double> f_coo; std::vector<double> g_coo; };

double compute_angle_functional(const gsl_vector *v, void * p)
{
	struct my_f_params * params = (struct my_f_params *) p;
	std::vector<double> f_coo = params->f_coo;
	std::vector<double> g_coo = params->g_coo;
	const double M_2PI = 2. * M_PI;
	const uint32_t N = f_coo.size();

	double sum = 0.;
	double tmp_sum;
	double t = gsl_vector_get(v, 0);
	for (uint32_t i = 0; i < N - 1; i++)
	{
		tmp_sum = g_coo[i] * g_coo[i+1] * sin(M_2PI / N)
			/ (g_coo[i+1] * sin(M_2PI * (1 - t) / N) + g_coo[i] * sin(M_2PI * t / N))
			- f_coo[i];
		sum = sum + (tmp_sum * tmp_sum);
	}
	{
		tmp_sum = g_coo[N - 1] * g_coo[0] * sin(M_2PI / N)
			/ (g_coo[0] * sin(M_2PI * (1 - t) / N) + g_coo[N - 1] * sin(M_2PI * t / N))
			- f_coo[N - 1];
		sum = sum + (tmp_sum * tmp_sum);
	}

	return sum;
}

void compute_derivative_angle_functional(const gsl_vector *v, void *p, gsl_vector *df)
{
	struct my_f_params * params = (struct my_f_params *) p;
	std::vector<double> f_coo = params->f_coo;
	std::vector<double> g_coo = params->g_coo;
	const uint32_t N = f_coo.size();
	const double M_2PI = 2. * M_PI;
	const double M_2PI_N = M_2PI / N;

	const double t = gsl_vector_get(v, 0);
	const double _1_t = 1. - t;

	double a(0.), b(0.), b_der(0);
	double sum(0.), tmp_sum(0.);
	for (uint32_t i = 0; i < N - 1; i++)
	{
		b = g_coo[i+1] * sin(M_2PI_N * _1_t) + g_coo[i] * sin(M_2PI_N * t);

		if (fabs(b) < 1.e-15)
		{
			gsl_vector_set(df, 0, GSL_NAN);
			printf("b == 0\n");
			return;
		}

		a = g_coo[i] * g_coo[i+1] * sin(M_2PI_N);
		b_der = M_2PI_N * (g_coo[i] * cos(M_2PI_N * t) - g_coo[i+1] * cos(M_2PI_N * (_1_t)));

		tmp_sum = -2. * a * b_der * (a / b - f_coo[i]) / (b*b);
		sum = sum + tmp_sum;
	}
	{
		b = g_coo[0] * sin(M_2PI_N * _1_t) + g_coo[N - 1] * sin(M_2PI_N * t);

		if (fabs(b) < 1.e-15)
		{
			gsl_vector_set(df, 0, GSL_NAN);
			printf("b == 0\n");
			return;
		}

		a = g_coo[N - 1] * g_coo[0] * sin(M_2PI_N);
		b_der = M_2PI_N * (g_coo[N - 1] * cos(M_2PI_N * t) - g_coo[0] * cos(M_2PI_N * (_1_t)));

		tmp_sum = -2. * a * b_der * (a / b - f_coo[N - 1]) / (b*b);
		sum = sum + tmp_sum;
	}

	gsl_vector_set(df, 0, sum);
}

void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
	*f = compute_angle_functional(x, params);
	compute_derivative_angle_functional(x, params, df);
}

struct EngVal
{
	double t;
	double s2;
	EngVal(double par, double s2_par) : t(par), s2(s2_par) {}
	EngVal() : t(0.), s2(0.) {}
};

uint32_t Solver::edit_angle_quad(double &ret_rot_angle)
{

	const double M_2PI = 2. * M_PI;
	//Construct equation coeffs
	for (auto task : tasks)
	{
		const Polygon & f = task.f;
		const Polygon & g = task.g;
		const std::vector<double> & g_coo = g.get_poly_on_mesh();
		const std::vector<double> & f_coo = f.get_poly_on_mesh();

		const uint32_t N = task.f.get_mesh_size();
		//TODO: calculate rot_angle first

		/////////////////////////////////////////////
		double border_value_0(0.), border_value_1(0.);
		for (uint32_t i = 0; i < N; i++)
		{
			double tmp_sum = g_coo[i] * g_coo[i] - f_coo[i] * f_coo[i];
			border_value_0 = border_value_0 + tmp_sum * tmp_sum;
		}
		for (uint32_t i = 0; i < N-1; i++)
		{
			double tmp_sum = g_coo[i+1] * g_coo[i+1] - f_coo[i] * f_coo[i];
			border_value_1 = border_value_1 + tmp_sum * tmp_sum;
		}
		{
			double tmp_sum = g_coo[0] * g_coo[0] - f_coo[N - 1] * f_coo[N - 1];
			border_value_1 = border_value_1 + tmp_sum * tmp_sum;
		}
		//(t, s2(t))
		EngVal border_min(0., 0.);
		if (border_value_0 < border_value_1)
		{			
			border_min.t = 0.;
			border_min.s2 = border_value_0;
		}
		else
		{
			border_min.t = 1.;
			border_min.s2 = border_value_1;			
		}

		
		//////////////////////////////////////////////

		double *a = new double[N];
		double a_coeff(0.), b_coeff(0.), c_coeff(0.), d_coeff(0.);
		double tmp_sum_A(0.), tmp_sum_B(0.), tmp_sum_C(0.);
		//a, b, c and A, B, C are absolutely independent 
		double gg = 0.;
		double ggcos = 0.;
		for (uint32_t i = 0; i < N-1; i++)
		{
			gg = g_coo[i] * g_coo[i];
			ggcos = g_coo[i] * g_coo[i + 1] * cos(M_2PI / N);

			tmp_sum_A = gg + g_coo[i + 1] * g_coo[i + 1] - 2 * ggcos;

			tmp_sum_B = ggcos - gg;

			tmp_sum_C = gg - f_coo[i] * f_coo[i];

			//Construct equation coeffs
			a_coeff = a_coeff + (tmp_sum_A * tmp_sum_A);
			b_coeff = b_coeff + (tmp_sum_A * tmp_sum_B);
			c_coeff = c_coeff + (2. * tmp_sum_B + tmp_sum_C * tmp_sum_A);
			d_coeff = d_coeff + (tmp_sum_B * tmp_sum_C);
		}
		{
			gg = g_coo[N - 1] * g_coo[N - 1];
			ggcos = g_coo[N - 1] * g_coo[0] * cos(M_2PI / N);

			tmp_sum_A = gg + g_coo[0] * g_coo[0] - 2 * ggcos;

			tmp_sum_B = ggcos - gg;

			tmp_sum_C = gg - f_coo[N - 1] * f_coo[N - 1];

			//Construct equation coeffs
			a_coeff = a_coeff + (tmp_sum_A * tmp_sum_A);
			b_coeff = b_coeff + (tmp_sum_A * tmp_sum_B);
			c_coeff = c_coeff + (2. * tmp_sum_B + tmp_sum_C * tmp_sum_A);
			d_coeff = d_coeff + (tmp_sum_B * tmp_sum_C);
		}
		b_coeff *= 3.;

		if(fabs(a_coeff) < 1.e-15 && fabs(b_coeff) < 1.e-15 && fabs(c_coeff) < 1.e-15)
		{
			std::cout << "Solver::edit_angle_quad(): wrong a_coeff, b_coeff, c_coeff\n";
			return -1;
		}
		
		if (fabs(a_coeff) < 1.e-15)
		{
			if (fabs(b_coeff) < 1.e-15)
			{
				//This is linear equation
				ret_rot_angle = -d_coeff / c_coeff;
				double value_0 = compute_angle_quad_functional(g_coo, f_coo, ret_rot_angle);
				ret_rot_angle = value_0 < border_min.s2 ? ret_rot_angle : border_min.t;
				if (ret_rot_angle < 0. || ret_rot_angle > 1.)
					ret_rot_angle = border_min.t;
				return 1;
			}

			//This is quadratic equation
			c_coeff /= b_coeff;
			d_coeff /= b_coeff;
			double * x = new double[2];
			int32_t ret_val = SolveP2(x, c_coeff, d_coeff);
			if (ret_val == 0)
				std::cout << "Solver::edit_angle_quad() : no real root in quadratic equation\n";

			EngVal value_x_0(x[0], compute_angle_quad_functional(g_coo, f_coo, x[0]));
			EngVal value_x_1(x[1], compute_angle_quad_functional(g_coo, f_coo, x[1]));

			std::vector<EngVal> vals2comp;
			vals2comp.push_back(border_min);
			if (!(value_x_0.t < 0. || value_x_0.t > 1.))
				vals2comp.push_back(value_x_0);
			if (!(value_x_1.t < 0. || value_x_1.t > 1.))
				vals2comp.push_back(value_x_1);

			double min_s = vals2comp[0].s2;
			ret_rot_angle = vals2comp[0].t;
			for (auto value : vals2comp)
			{
				if (value.s2 < min_s)
				{
					min_s = value.s2;
					ret_rot_angle = value.t;
				}
			}
			return 1;
		}

		//This is cubic equation
		b_coeff /= a_coeff;
		c_coeff /= a_coeff;
		d_coeff /= a_coeff;
		double * x = new double[3];
		int32_t ret_val = SolveP3(x, b_coeff, c_coeff, d_coeff);
		if (ret_val == 3)
		{
			EngVal value_x_0(x[0], compute_angle_quad_functional(g_coo, f_coo, x[0]));
			EngVal value_x_1(x[1], compute_angle_quad_functional(g_coo, f_coo, x[1]));
			EngVal value_x_2(x[2], compute_angle_quad_functional(g_coo, f_coo, x[2]));


			std::vector<EngVal> vals2comp;
			vals2comp.push_back(border_min);
			if (!(value_x_0.t < 0. || value_x_0.t > 1.))
				vals2comp.push_back(value_x_0);
			if (!(value_x_1.t < 0. || value_x_1.t > 1.))
				vals2comp.push_back(value_x_1);
			if (!(value_x_2.t < 0. || value_x_2.t > 1.))
				vals2comp.push_back(value_x_2);

			double min_s = vals2comp[0].s2;
			ret_rot_angle = vals2comp[0].t;
			for (auto value : vals2comp)
			{
				if (value.s2 < min_s)
				{
					min_s = value.s2;
					ret_rot_angle = value.t;
				}
			}
		}
		else
		{
			EngVal value_x_0(x[0], compute_angle_quad_functional(g_coo, f_coo, x[0]));

			std::vector<EngVal> vals2comp;
			vals2comp.push_back(border_min);
			if (!(value_x_0.t < 0. || value_x_0.t > 1.))
				vals2comp.push_back(value_x_0);

			double min_s = vals2comp[0].s2;
			ret_rot_angle = vals2comp[0].t;
			for (auto value : vals2comp)
			{
				if (value.s2 < min_s)
				{
					min_s = value.s2;
					ret_rot_angle = value.t;
				}
			}
		}

	}
	return 1;
}

uint32_t Solver::edit_angle(double &ret_rot_angle)
{
	for (auto task : tasks)
	{
		int status;
		int iter = 0, max_iter = 1000;
		const gsl_multimin_fdfminimizer_type  *T;
		gsl_multimin_fdfminimizer  *s;

		gsl_multimin_function_fdf F;
		struct my_f_params params = { task.f.get_poly_on_mesh(), task.g.get_poly_on_mesh() };
		F.n = 1;
		F.f = &compute_angle_functional;
		F.df = &compute_derivative_angle_functional;
		F.fdf = &fdf;
		F.params = &params;

		gsl_vector *x = gsl_vector_alloc(1);
		gsl_vector_set(x, 0, 0.5);

		T = gsl_multimin_fdfminimizer_vector_bfgs2;
		s = gsl_multimin_fdfminimizer_alloc(T,1);
		gsl_multimin_fdfminimizer_set(s, &F, x, 0.5, 1.e-1);

		do
		{
			iter++;
			status = gsl_multimin_fdfminimizer_iterate(s);

			if (status)
				break;

			status = gsl_multimin_test_gradient(s->gradient, 1.e-15);


		} while (status == GSL_CONTINUE && iter < max_iter);

		printf("iter = %d\n", iter);

		gsl_vector *_shift = gsl_vector_alloc(1);
		gsl_vector_set(_shift, 0, ret_rot_angle);
		double s1_shift = compute_angle_functional(_shift, &params);

		x = gsl_multimin_fdfminimizer_x(s);
		ret_rot_angle = gsl_vector_get(x, 0);
		gsl_vector *_0 = gsl_vector_alloc(1);
		gsl_vector_set(_0, 0, 0.);
		gsl_vector *_1 = gsl_vector_alloc(1);
		gsl_vector_set(_1, 0, 1.);
		

		double s1_0 = compute_angle_functional(_0, &params);
		double s1_1 = compute_angle_functional(_1, &params);
		double s1_t = compute_angle_functional(x, &params);
		printf("\n%17s | %17s | %17s | %17s\n", "s1(0)", "s1(t)", "s1(1)", "s1(shift)");
		printf("%.15f | %.15f | %.15f | %.15f\n", s1_0, s1_t, s1_1, s1_shift);

		gsl_multimin_fdfminimizer_free(s);
		return status;
	}
	return 1;
}

double Solver::edge2edge_cut_angles_mismatching(double t)
{
	auto & f = tasks[0].f.get_poly_on_mesh();
	auto & g = tasks[0].g.get_poly_on_mesh();

	const uint32_t N = f.size();
	const uint32_t N_ang = tasks[0].f.get_angles_number();
	const double cos_ = cos(2. * M_PI / N);

	if (N % N_ang != 0 || N < N_ang)
		return -1.;

	const uint32_t verts2poly_edge = N / N_ang;
	double tmp(0.), measure(0.);
	{
		tmp = shift_point2(g[N-1], g[0], t, N);
		tmp = tmp - f[0];
		measure = measure + tmp * tmp;
	}
	for (uint32_t i = verts2poly_edge; i < N; i += verts2poly_edge)
	{
		tmp = shift_point2(g[i-1], g[i], t, N);
		tmp = tmp - f[i];
		measure = measure + tmp * tmp;
	}
	return measure;
}
