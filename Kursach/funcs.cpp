#include "stdafx.h"

double pho(int j, int N_2)
{
	if (j != 0 && j != N_2)
		return 1.;
	else return 0.5;
}

double coeff_a(int k, int N, double *r)
{
	/*что делать с ро?*/
	double coeff_a = 0.;
	double inv_N = 1. / N;
	double N_2 = N / 2;
	for (int j = 0; j < N; j++)
	{
		coeff_a += r[j] * cos(2. * M_PI * k * j * inv_N);
	}
	coeff_a *= pho(k, N_2);
	return coeff_a;
}

double coeff_b(int k, int N, double *r)
{
	double coeff_b = 0.;
	double inv_N = 1. / N;
	for (int j = 1; j < N; j++)
	{
		coeff_b += r[j] * sin(2. * M_PI * k * j * inv_N);
	}
	return coeff_b;
}