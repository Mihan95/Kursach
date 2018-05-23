#include <stdio.h>
#include <stdlib.h>
#include "stdafx.h"
#include "Solver.h"

#define TEST_REF
//TO DO: добавить умножение на константы и на N во все нормы
int main(int argc, char* argv[])
{
	const uint32_t N = 128;

#ifdef TEST_REF

	Polygon f(4, N);
	f.create_from_regular_poly(4);
	f.project_to_mesh(31);

	Polygon g(4, 4);
	g.create_from_regular_poly(4);
	g.project_to_mesh(0);

	func task;
	task.f = f;
	task.g = g;
	task.w.resize(N);
	std::fill(task.w.begin(), task.w.end(), 1.);
	//std::fill(task.w.begin(), task.w.begin() + 3*N/4+1, 1.);
	//std::fill(task.w.begin() + 3 * N / 4 + 1, task.w.end(), 0.);

	Solver slv(task);
	slv.is_refinement();

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