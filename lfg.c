#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "lfg.h"

typedef struct {

	double		seed[55];
	int		ra, rb;
}
lfg_t;

static lfg_t		lfg;

static uint32_t
lfg_lcgu(uint32_t rseed)
{
	/* Linear Congruential generator.
	 * */

	return rseed * 17317U + 1U;
}

void lfg_start(int seed)
{
	uint32_t	lcgu;
	int		i;

	lcgu = lfg_lcgu(seed);
	lcgu = lfg_lcgu(lcgu);

	for (i = 0; i < 55; ++i) {

		lfg.seed[i] = (double) (lcgu = lfg_lcgu(lcgu)) / 4294967296.;
	}

	lfg.ra = 0;
	lfg.rb = 31;
}

double lfg_rand()
{
	double		x, a, b;

	/* Lagged Fibonacci generator.
	 * */

	a = lfg.seed[lfg.ra];
	b = lfg.seed[lfg.rb];

	x = (a < b) ? a - b + 1. : a - b - 1.;

	lfg.seed[lfg.ra] = x;

	lfg.ra = (lfg.ra < 54) ? lfg.ra + 1 : 0;
	lfg.rb = (lfg.rb < 54) ? lfg.rb + 1 : 0;

	return x;
}

double lfg_gauss()
{
	double		x, s;

	/* Box-Muller transform.
	 * */

	do {
		s = lfg_rand();
		x = lfg_rand();

		s = s * s + x * x;

		if (s > 0. && s < 1.)
			break;
	}
	while (1);

	x *= sqrt(- 2. * log(s) / s);

	return x;
}

