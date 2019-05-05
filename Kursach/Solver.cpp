#include "stdafx.h"
#include "Solver.h"
#include <array>
#include <memory>
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

static double shift_point(const double r1, const double r2, const double t, const uint32_t N)
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
		new_on_mesh[0] = shift_point(on_mesh[0], on_mesh[N-1], t, N);
	}
	for (uint32_t i = 1; i < N; i++)
	{
		new_on_mesh[i] = shift_point(on_mesh[i], on_mesh[i - 1], t, N);
	}
	on_mesh = new_on_mesh;
}

void Polygon::shift_poly_on_mesh_op(double t)
{
	const uint32_t N = on_mesh.size();
	const double cos_ = cos(2. * M_PI / N);

	std::vector<double> tmp_on_mesh(N);

	for (uint32_t i = 0; i < N-1; i++)
	{
		tmp_on_mesh[i] = shift_point(on_mesh[i], on_mesh[i + 1], t, N);
	}
	{
		tmp_on_mesh[N-1] = shift_point(on_mesh[N - 1], on_mesh[0], t, N);
	}
	on_mesh = tmp_on_mesh;
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
		/*����������� ���������� ������������*/
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

		/*����������� ���������� ������������*/
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
	std::vector<double> one_point_norm(N);
	uint32_t n = N / 4;
	for (uint32_t i = 0; i < N - 1; i++)
	{
		//	if (i == 0 || i == n || i == 2 * n || i == 3 * n)
		//		continue;
		tmp_sum = g_coo[i] * g_coo[i+1] * sin(M_2PI / N)
			/ (g_coo[i+1] * sin(M_2PI * (1 - t) / N) + g_coo[i] * sin(M_2PI * t / N))
			- f_coo[i];
		sum = sum + (tmp_sum * tmp_sum);
		one_point_norm[i] = tmp_sum * tmp_sum;
	}
	{
		tmp_sum = g_coo[N - 1] * g_coo[0] * sin(M_2PI / N)
			/ (g_coo[0] * sin(M_2PI * (1 - t) / N) + g_coo[N - 1] * sin(M_2PI * t / N))
			- f_coo[N - 1];
		sum = sum + (tmp_sum * tmp_sum);
		one_point_norm[N-1] = tmp_sum * tmp_sum;
	}

	/***********************************************************/

	double average_norm = sum / N;
	std::vector<double> angle_points_norm;
	for (uint32_t i = 0; i < N - 1; i++)
	{
		if (i == 0 || i == n || i == 2 * n || i == 3 * n)
		{
			tmp_sum = g_coo[i] * g_coo[i + 1] * sin(M_2PI / N)
				/ (g_coo[i + 1] * sin(M_2PI * (1 - t) / N) + g_coo[i] * sin(M_2PI * t / N))
				- f_coo[i];
			sum = sum + (tmp_sum * tmp_sum);
			angle_points_norm.push_back( tmp_sum * tmp_sum );
		}
	}
	std::sort(one_point_norm.begin(), one_point_norm.end(), std::greater<double>());

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
	uint32_t n = N / 4;
	for (uint32_t i = 0; i < N - 1; i++)
	{
		//if (i == 0 || i == n || i == 2 * n || i == 3 * n)
		//	continue;
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
