#include <stdio.h>
#include <stdlib.h>
#include "stdafx.h"
#include "Solver.h"
#include "angle_error.h"

#define TEST_REF
//TO DO: добавить умножение на константы и на N во все нормы
int main(int argc, char* argv[])
{
	const uint32_t N = 128;
	const double shift = 0.25;

#ifdef TEST_REF

	Polygon f(4, N);
	f.create_from_regular_poly(4);
	f.project_to_mesh(0);

	Polygon g(4, N);
	g.create_from_regular_poly(4);
	g.project_to_mesh(0);
	g.shift_poly_on_mesh(shift);
	//g.shift_poly_on_mesh_op(shift);



	func task;
	task.f = f;
	task.g = g;
	task.w.resize(N);
	std::fill(task.w.begin(), task.w.end(), 1.);
	//std::fill(task.w.begin(), task.w.begin() + 3*N/4+1, 1.);
	//std::fill(task.w.begin() + 3 * N / 4 + 1, task.w.end(), 0.);

	Solver slv(task);
	double t(shift);
	slv.edit_angle(t);
	//g.shift_poly_on_mesh_op(shift);
	std::cout.precision(15);
	std::cout << "\nminimum: " << std::fixed << t << std::endl;
	std::cout << "origin_shift - finded_shift = " << std::fixed << /*2.*M_PI**/fabs(t-shift)/*/N */<< std::endl;

	std::cout << "Angle norm = " << max_angle_norm(f, g) << std::endl;
#else
	Polygon f(4, N);
	f.create_from_regular_poly(4);
	f.project_to_mesh(0);

	Polygon g(4, N);
	g.create_from_regular_poly(4);
	g.project_to_mesh(15);

	func task;
	task.f = f;
	task.g = g;
	task.w.resize(N);
	std::fill(task.w.begin(), task.w.end(), 1.);
	//std::fill(task.w.begin(), task.w.begin() + 3*N/4+1, 1.);
	//std::fill(task.w.begin() + 3 * N / 4 + 1, task.w.end(), 0.);

	Solver slv(task);
	slv.solve2();
#endif
	system("pause");
}