#include <stdio.h>
#include <stdlib.h>
#include "stdafx.h"
#include "Solver.h"

int main(int argc, char* argv[])
{
	Polygon f(4, 16);
	f.create_from_regular_poly(4);
	f.project_to_mesh(0);

	Polygon g(4, 16);
	g.create_from_regular_poly(4);
	g.project_to_mesh(0);

	func task;
	task.f = f;
	task.g = g;

	Solver slv(task);
	slv.solve();

	system("pause");
}