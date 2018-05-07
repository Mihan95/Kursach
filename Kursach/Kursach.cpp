#include <stdio.h>
#include <stdlib.h>
#include "stdafx.h"
#include <complex>
#include <time.h>
#include <fftw3.h>

inline double function1(double, double*);

int main(int argc, char* argv[])
{
	double *xpoints;
	
	xpoints = (double*)malloc(Npoints * sizeof(double));

	for (i = 0; i < Npoints; i++)
	{
		xpoints[i] = function1(t, coef);
	}

	fftw_complex *in, *out;
	fftw_plan plan;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Npoints);			//pay attention
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Npoints);		//pay attention
	plan = fftw_plan_dft_1d(Npoints, in, out, FFTW_FORWARD, FFTW_ESTIMATE); 	//Here we set which kind of transformation we want to perform

	for (i = 0; i < Npoints; i++)
	{
		*in[i] = xpoints[i];
	}
	fftw_execute(plan); //Execution of FFT

	fftw_destroy_plan(plan);	 //Free memory
	fftw_free(in);			 //Free memory
	fftw_free(out);			 //Free memory
	free(xpoints);
	return 0;
}