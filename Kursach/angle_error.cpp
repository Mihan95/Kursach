#include "stdafx.h"
#include "angle_error.h"

using Point = std::pair<double, double>;

Point lineLineIntersection(Point A, Point B, Point C, Point D)
{
	// Line AB represented as a1x + b1y = c1 
	double a1 = B.second - A.second;
	double b1 = A.first - B.first;
	double c1 = a1 * (A.first) + b1 * (A.second);

	// Line CD represented as a2x + b2y = c2 
	double a2 = D.second - C.second;
	double b2 = C.first - D.first;
	double c2 = a2 * (C.first) + b2 * (C.second);

	double determinant = a1 * b2 - a2 * b1;

	if (fabs(determinant) < EPS)
	{
		// The lines are parallel. This is simplified 
		// by returning a pair of DBL_MAX 
		return std::make_pair(DBL_MAX, DBL_MAX);
	}
	else
	{
		double x = (b2*c1 - b1 * c2) / determinant;
		double y = (a1*c2 - a2 * c1) / determinant;
		return std::make_pair(x, y);
	}
}

double find_r_A_angle(double r, uint32_t i, Point A, uint32_t N)
{
	const double M_2PI_i_N = 2. * M_PI * i / N;
	double t = (cos(M_2PI_i_N) * A.first + sin(M_2PI_i_N) * A.second)
		/ sqrt(A.first*A.first + A.second*A.second);
	if (t < -1.) t = -1.;
	else if (t > 1.) t = 1.;
	return acos(t);
}

bool is_point_on_segment(Point A, Point B, Point P)
{
	double a1 = B.second - A.second;
	double b1 = A.first - B.first;
	double c1 = a1 * (A.first) + b1 * (A.second);

	double min_x, min_y, max_x, max_y;
	min_x = std::min(A.first, B.first);
	max_x = std::max(A.first, B.first);
	min_y = std::min(A.second, B.second);
	max_y = std::max(A.second, B.second);

	if (fabs(a1 * P.first + b1 * P.second + c1) < EPS)
	{
		if (P.first >= min_x && P.first <= max_x
			&& P.second >= min_y && P.second <= max_y)
			return true;
		if (fabs(P.first - min_x) < EPS && fabs(P.second - min_y) < EPS
			|| fabs(P.first - max_x) < EPS && fabs(P.second - min_y) < EPS
			|| fabs(P.first - min_x) < EPS && fabs(P.second - max_y) < EPS
			|| fabs(P.first - max_x) < EPS && fabs(P.second - max_y) < EPS)
			return true;
	}
	return false;
}

double compute_polar_line(double r1, double r2, double alpha, double t)
{
	return r1 * r2 * sin(alpha) / (r2 * sin(alpha*(1.-t)) + r1 * sin(alpha * t));
}

struct my_params { double r1; double r2; double alpha; double rho; };

double compute_polar_line_p(const gsl_vector *v, void * p)
{
	struct my_params * params = (struct my_params *) p;
	double r1 = params->r1;
	double r2 = params->r2;
	double alpha = params->alpha;
	double rho = params->rho;
	double t = gsl_vector_get(v, 0);

	if (t < 0. || t > 1.)
		return 0.;

	double r = compute_polar_line(r1, r2, alpha, t);
	double r_rho = fabs(r - rho);
	double ret_val = -r_rho * r_rho;
	return ret_val;
}

void compute_derivative_polar_line_p(const gsl_vector *v, void *p, gsl_vector *df)
{
	struct my_params * params = (struct my_params *) p;
	double r1 = params->r1;
	double r2 = params->r2;
	double alpha = params->alpha;
	double rho = params->rho;
	double t = gsl_vector_get(v, 0);

	const double _1_t = 1. - t;

	if(t < 0. || t > 1.)
	{
		gsl_vector_set(df, 0, GSL_NAN);
		printf("t < 0. || t > 1. == 0\n");
		return;
	}

	double b = r2 * sin(alpha * _1_t) + r1 * sin(alpha * t);

	if (fabs(b) < EPS)
	{
		gsl_vector_set(df, 0, GSL_NAN);
		//printf("b == 0\n");
		return;
	}

	double a = r1 * r2 * sin(alpha);
	double b_der = alpha * (r1 * cos(alpha * t) - r2 * cos(alpha * (_1_t)));

	double der = -2. * a * b_der * (a / b - rho) / (b*b);
	gsl_vector_set(df, 0, der);
}

void fdf_line(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
	*f = compute_polar_line_p(x, params);
	compute_derivative_polar_line_p(x, params, df);
}

double line_min(double r1, double r2, double alpha, double rho)
{
	int status;
	int iter = 0, max_iter = 1000;
	const gsl_multimin_fdfminimizer_type  *T;
	gsl_multimin_fdfminimizer  *s;

	gsl_multimin_function_fdf F;
	struct my_params params = { r1, r2, alpha, rho };
	F.n = 1;
	F.f = &compute_polar_line_p;
	F.df = &compute_derivative_polar_line_p;
	F.fdf = &fdf_line;
	F.params = &params;
	
	gsl_vector *x = gsl_vector_alloc(1);
	gsl_vector_set(x, 0, 0.5);

	T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc(T, 1);
	gsl_multimin_fdfminimizer_set(s, &F, x, 0.5, 1.e-1);

	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 1.e-15);


	} while (status == GSL_CONTINUE && iter < max_iter);

	x = gsl_multimin_fdfminimizer_x(s);
	gsl_vector *_0 = gsl_vector_alloc(1);
	gsl_vector_set(_0, 0, 0.);
	gsl_vector *_1 = gsl_vector_alloc(1);
	gsl_vector_set(_1, 0, 1.);

	double min_0 = compute_polar_line_p(_0, &params);
	double min_1 = compute_polar_line_p(_1, &params);
	double min_t = compute_polar_line_p(x, &params);

	return std::min({ min_0, min_1, min_t });
}

void print_(std::vector<double> on_mesh)
{
	std::ofstream fout;
	fout.open("output.txt");
	int N = on_mesh.size();
	double M_2PI_N = 2. * M_PI / N;
	for (uint32_t i = 0; i < N; i++)
	{
		fout << on_mesh[i] * cos(M_2PI_N * i) << "   " << on_mesh[i] * sin(M_2PI_N * i) << std::endl;
	}

	fout.close();
}

double max_angle_norm(const Polygon & f_pol, const Polygon & g_pol)
{
	const auto & f = f_pol.get_poly_on_mesh();
	const auto & g = g_pol.get_poly_on_mesh();
	const uint32_t N = f_pol.get_mesh_size();
	const double M_2PI_N = 2. * M_PI / N;
	double global_max = 0;

	auto single_angle_norm = [&](uint32_t i0, uint32_t i1, uint32_t i2, uint32_t i3)
	{
		Point A(g[i0] * cos(M_2PI_N * (i0)), g[i0] * sin(M_2PI_N * (i0))),
			 B(g[i1] * cos(M_2PI_N * (i1)), g[i1] * sin(M_2PI_N * (i1))),
			 C(g[i2] * cos(M_2PI_N * (i2)), g[i2] * sin(M_2PI_N * (i2))),
			 D(g[i3] * cos(M_2PI_N * (i3)), g[i3] * sin(M_2PI_N * (i3)));			  

		Point P = lineLineIntersection(A, B, C, D);
		if (P.first == DBL_MAX || P.second == DBL_MAX) {
			//std::cout << "Parallel lines\n";
			return;
		}
		if (is_point_on_segment(A, B, P) || is_point_on_segment(C, D, P)) {
			std::cout << "Intersected segments\n";
			return;
		}
		Point g_p(g[i1] * cos(M_2PI_N * i1), g[i1] * sin(M_2PI_N * i1));
		double angle = find_r_A_angle(g[i1], i1, P, N);
		if (angle >= 2. * M_PI)
			angle -= 2. * M_PI;
		if (angle < 0.)
			angle += 2. * M_PI;
		if (angle < EPS)
			return;
		double r_P = sqrt(P.first*P.first + P.second*P.second);
		double min1 = line_min(g[i1], r_P, angle, f[i1]);
		double min2 = line_min(r_P, g[i2], fabs(M_2PI_N - angle), f[i1]);
		double min = std::min(min1, min2);
		double max = -min;
		if (max > global_max)
			global_max = max;

	};

	for (uint32_t i = 1; i < N - 2; i++) {
		single_angle_norm(i - 1, i, i + 1, i + 2);
	}
	
	single_angle_norm(N - 3, N - 2, N - 1, 0);
	single_angle_norm(N - 2, N - 1, 0, 1);
	single_angle_norm(N - 1, 0, 1, 2);

	return global_max;
}