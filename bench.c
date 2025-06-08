#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lse.h"
#include "lfg.h"

#define LSE_EPSF		5E-7

#define LSE_printf(s)		fprintf(stderr, "%s in %s:%i\n", (s), __FILE__, __LINE__)
#define LSE_assert(x)		if ((x) == 0) { LSE_printf(#x); exit(-1); }

#define LSE_assert_relative(x, r)	LSE_assert(fabs((x) - (r)) < LSE_EPSF * (r))
#define LSE_assert_absolute(x, r)	LSE_assert(fabs((x) - (r)) < LSE_EPSF)
#define LSE_assert_approx(x, r, a)	LSE_assert(fabs((x) - (r)) < LSE_EPSF * (a))

void lse_basic_test()
{
	lse_t		ls, lb;
	lse_float_t	v[LSE_FULL_MAX];

	lse_construct(&ls, 1, 4, 1);

	v[0] = (lse_float_t) 0.;
	v[1] = (lse_float_t) - 3.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.;
	v[4] = (lse_float_t) 2.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) 0.;
	v[2] = (lse_float_t) - 1.;
	v[3] = (lse_float_t) 1.;
	v[4] = (lse_float_t) 3.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) - 7.;
	v[1] = (lse_float_t) 2.;
	v[2] = (lse_float_t) - 1.;
	v[3] = (lse_float_t) 0.;
	v[4] = (lse_float_t) - 6.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) - 1.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.;
	v[4] = (lse_float_t) 7.;

	lse_insert(&ls, v);
	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 2.);
	LSE_assert_relative(ls.sol.m[2], 3.);
	LSE_assert_relative(ls.sol.m[3], 4.);

	lse_construct(&ls, 1, 4, 1);
	lse_nostd(&ls);

	v[0] = (lse_float_t) 0.;
	v[1] = (lse_float_t) - 3.E+30;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.E+30;
	v[4] = (lse_float_t) 2.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 2.E+30;
	v[1] = (lse_float_t) 0.;
	v[2] = (lse_float_t) - 1.E+30;
	v[3] = (lse_float_t) 1.E+30;
	v[4] = (lse_float_t) 3.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) - 7.E+30;
	v[1] = (lse_float_t) 2.E+30;
	v[2] = (lse_float_t) - 1.E+30;
	v[3] = (lse_float_t) 0.;
	v[4] = (lse_float_t) - 6.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.E+30;
	v[1] = (lse_float_t) - 1.E+30;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.E+30;
	v[4] = (lse_float_t) 7.;

	lse_insert(&ls, v);
	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.E-30);
	LSE_assert_relative(ls.sol.m[1], 2.E-30);
	LSE_assert_relative(ls.sol.m[2], 3.E-30);
	LSE_assert_relative(ls.sol.m[3], 4.E-30);

	lse_construct(&ls, 2, 1, 4);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) 1.;
	v[2] = (lse_float_t) - 2.;
	v[3] = (lse_float_t) 4.;
	v[4] = (lse_float_t) 1.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) - 2.;
	v[2] = (lse_float_t) 4.;
	v[3] = (lse_float_t) 1.;
	v[4] = (lse_float_t) - 2.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) 4.;
	v[2] = (lse_float_t) 1.;
	v[3] = (lse_float_t) - 2.;
	v[4] = (lse_float_t) 4.;

	lse_insert(&ls, v);

	lse_solve(&ls);
	lse_std(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 1.);
	LSE_assert_relative(ls.sol.m[2], 1.);
	LSE_assert_relative(ls.sol.m[3], 1.);

	LSE_assert_relative(ls.std.m[0], 3.);
	LSE_assert_relative(ls.std.m[1], 3.);
	LSE_assert_relative(ls.std.m[2], 3.);
	LSE_assert_relative(ls.std.m[3], 3.);

	lse_construct(&ls, 2, 1, 1);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) 1.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) - 2.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) 4.;

	lse_insert(&ls, v);

	lse_solve(&ls);
	lse_std(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.std.m[0], 3.);

	lse_construct(&ls, 3, 3, 1);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) - 3.;
	v[2] = (lse_float_t) 2.;
	v[3] = (lse_float_t) 1.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) 0.;
	v[2] = (lse_float_t) - 1.;
	v[3] = (lse_float_t) - 1.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) - 7.;
	v[1] = (lse_float_t) 2.;
	v[2] = (lse_float_t) - 1.;
	v[3] = (lse_float_t) - 6.;

	lse_insert(&ls, v);
	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 2.);
	LSE_assert_relative(ls.sol.m[2], 3.);

	v[0] = (lse_float_t) - 1.;
	v[1] = (lse_float_t) 1.;
	v[2] = (lse_float_t) 1.;
	v[3] = (lse_float_t) 4.;

	lse_insert(&ls, v);
	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 2.);
	LSE_assert_relative(ls.sol.m[2], 3.);

	lse_construct(&ls, 4, 3, 1);

	v[0] = (lse_float_t) 3.;
	v[1] = (lse_float_t) 1.;
	v[2] = (lse_float_t) 1.;
	v[3] = (lse_float_t) 5.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) 2.;
	v[2] = (lse_float_t) 2.;
	v[3] = (lse_float_t) 6.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 3.;
	v[1] = (lse_float_t) 4.;
	v[2] = (lse_float_t) 2.;
	v[3] = (lse_float_t) 9.;

	lse_insert(&ls, v);

	lse_esv(&ls, 2);

	printf("\ncase  ---- basic ridge ----\n");

	printf("esv  % .8lE % .8lE\n", ls.esv.min, ls.esv.max);

	lse_ridge(&ls, ls.n_len_of_x * ls.esv.max * LSE_EPSF);
	lse_solve(&ls);

	printf("sol  % .8lE % .8lE % .8lE\n", ls.sol.m[0],
			ls.sol.m[1], ls.sol.m[2]);

	LSE_assert_approx(ls.sol.m[0], 1., 10.);
	LSE_assert_approx(ls.sol.m[1], 1., 10.);
	LSE_assert_approx(ls.sol.m[2], 1., 10.);

	lse_construct(&ls, 1, 3, 1);

	v[0] = (lse_float_t) 0.;
	v[1] = (lse_float_t) 4.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 4.;

	lse_insert(&ls, v);
	lse_forget(&ls, (lse_float_t) .5);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) 0.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.;

	lse_insert(&ls, v);
	lse_forget(&ls, (lse_float_t) .5);

	v[0] = (lse_float_t) 0.;
	v[1] = (lse_float_t) 0.;
	v[2] = (lse_float_t) 1.;
	v[3] = (lse_float_t) 1.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) 2.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 1.;

	lse_insert(&ls, v);
	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1. / 3.);
	LSE_assert_relative(ls.sol.m[1], 1. / 3.);
	LSE_assert_relative(ls.sol.m[2], 1.);

	lse_std(&ls);

	LSE_assert_relative(ls.std.m[0], 1. / sqrt(3.));

	lse_construct(&ls, 3, 3, 1);
	lse_nostd(&ls);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) -1.;
	v[2] = (lse_float_t) 1.;
	v[3] = (lse_float_t) 3.;

	lse_insert(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) 1.;
	v[2] = (lse_float_t) 2.;
	v[3] = (lse_float_t) 9.;

	lse_insert(&ls, v);

	lse_construct(&lb, 2, 3, 1);
	lse_nostd(&lb);

	v[0] = (lse_float_t) -1.;
	v[1] = (lse_float_t) -1.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) -3.;

	lse_insert(&lb, v);
	lse_merge(&ls, &lb);

	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 2.);
	LSE_assert_relative(ls.sol.m[2], 3.);
}

void lse_large_test(int n_full)
{
	lse_t		*ls;
	lse_float_t	*v;

	int		len, n, i;

	len = lse_getsize(LSE_CASCADE_MAX, n_full);

	ls = (lse_t *) malloc(len);
	v = (lse_float_t *) malloc(n_full * sizeof(lse_float_t));

	lse_construct(ls, LSE_CASCADE_MAX, n_full - 1, 1);
	lse_nostd(ls);

	for (n = 0; n < n_full - 1; ++n) {

		for (i = 0; i < n; ++i)
			v[i] = (lse_float_t) 1.;

		v[n] = (lse_float_t) 2.;

		for (i = n + 1; i < n_full; ++i)
			v[i] = (lse_float_t) 1.;

		lse_insert(ls, v);
	}

	lse_solve(ls);

	printf("\ncase  ---- large (%i) ----\n", n_full);

	printf("sol  % .8lE % .8lE % .8lE % .8lE\n",
			ls->sol.m[0], ls->sol.m[1],
			ls->sol.m[2], ls->sol.m[3]);

	for (i = 0; i < ls->n_len_of_x; ++i) {

		LSE_assert_approx(ls->sol.m[i], 1. / (double) n_full, 10.);
	}

	lse_esv(ls, 4);

	printf("esv  % .8lE % .8lE\n", ls->esv.min, ls->esv.max);

	LSE_assert_approx(ls->esv.max, n_full, 100. * n_full);
	LSE_assert_approx(ls->esv.min, 1., 100.);

	free(ls);
	free(v);
}

void lse_advanced_test()
{
	lse_t		ls;
	lse_float_t	v[LSE_FULL_MAX];

	int		i;

	lse_construct(&ls, LSE_CASCADE_MAX, 1, 2);
	lse_nostd(&ls);

	for (i = 0; i < 100000; ++i) {

		v[0] = (lse_float_t) 1.;
		v[1] = (lse_float_t) (i % 50);
		v[2] = (lse_float_t) - (i % 20);

		lse_insert(&ls, v);

		v[0] = (lse_float_t) 1.;
		v[1] = (lse_float_t) - (i % 50);
		v[2] = (lse_float_t) (i % 20);

		lse_insert(&ls, v);
	}

	lse_solve(&ls);

	printf("\ncase  ---- advanced ----\n");

	LSE_assert_absolute(ls.sol.m[0], 0.);
	LSE_assert_absolute(ls.sol.m[1], 0.);

	lse_construct(&ls, LSE_CASCADE_MAX, 1, 4);

	for (i = 0; i < 100000; ++i) {

		v[0] = (lse_float_t) 1.;
		v[1] = (lse_float_t) 1.;
		v[2] = (lse_float_t) 1.;
		v[3] = (lse_float_t) 1.;
		v[4] = (lse_float_t) 1.;

		lse_insert(&ls, v);
	}

	lse_solve(&ls);
	lse_std(&ls);

	printf("sol  % .8lE % .8lE % .8lE % .8lE\n",
			ls.sol.m[0], ls.sol.m[1],
			ls.sol.m[2], ls.sol.m[3]);

	printf("std  % .8lE % .8lE % .8lE % .8lE\n",
			ls.std.m[0], ls.std.m[1],
			ls.std.m[2], ls.std.m[3]);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 1.);
	LSE_assert_relative(ls.sol.m[2], 1.);
	LSE_assert_relative(ls.sol.m[3], 1.);

	LSE_assert_absolute(ls.std.m[0], 0.);
	LSE_assert_absolute(ls.std.m[1], 0.);
	LSE_assert_absolute(ls.std.m[2], 0.);
	LSE_assert_absolute(ls.std.m[3], 0.);
}

void lse_mean_std_test()
{
	lse_t		ls;
	lse_float_t	v[LSE_FULL_MAX];

	double		x[4];
	int		i;

	lfg_start(24);

	lse_construct(&ls, LSE_CASCADE_MAX, 1, 4);

	for (i = 0; i < 100000; ++i) {

		x[0] = 11. + lfg_gauss() * 2.;
		x[1] = - 3. + lfg_gauss() * 1.;
		x[2] = 5. + lfg_gauss() * 7.;
		x[3] = - 4. + lfg_gauss() * 5.;

		v[0] = (lse_float_t) 1.;
		v[1] = (lse_float_t) x[0];
		v[2] = (lse_float_t) x[1];
		v[3] = (lse_float_t) x[2];
		v[4] = (lse_float_t) x[3];

		lse_insert(&ls, v);
	}

	lse_solve(&ls);
	lse_std(&ls);

	printf("\ncase  ---- mean std ----\n");

	printf("sol  % .8lE % .8lE % .8lE % .8lE\n",
			ls.sol.m[0], ls.sol.m[1],
			ls.sol.m[2], ls.sol.m[3]);

	printf("std  % .8lE % .8lE % .8lE % .8lE\n",
			ls.std.m[0], ls.std.m[1],
			ls.std.m[2], ls.std.m[3]);

	LSE_assert_approx(ls.sol.m[0], 11.,  0.1 / LSE_EPSF);
	LSE_assert_approx(ls.sol.m[1], - 3., 0.1 / LSE_EPSF);
	LSE_assert_approx(ls.sol.m[2], 5.,   0.1 / LSE_EPSF);
	LSE_assert_approx(ls.sol.m[3], - 4., 0.1 / LSE_EPSF);

	LSE_assert_approx(ls.std.m[0], 2., 0.1 / LSE_EPSF);
	LSE_assert_approx(ls.std.m[1], 1., 0.1 / LSE_EPSF);
	LSE_assert_approx(ls.std.m[2], 7., 0.1 / LSE_EPSF);
	LSE_assert_approx(ls.std.m[3], 5., 0.1 / LSE_EPSF);
}

void lse_random_test()
{
	lse_t		ls;
	lse_float_t	v[LSE_FULL_MAX];

	double		x[5], z[5], b[25], rel[25];
	int		n, i;

	lfg_start(27);

	for (i = 0; i < 25; ++i) {

		b[i] = lfg_gauss();
	}

	lse_construct(&ls, LSE_CASCADE_MAX, 5, 5);

	for (n = 0; n < 5000000; ++n) {

		x[0] = lfg_urand() * 1.;
		x[1] = lfg_urand() * 1.;
		x[2] = lfg_urand() * 1.;
		x[3] = lfg_urand() * 1.;
		x[4] = lfg_urand() * 1.;

		z[0] = 0.;
		z[1] = 0.;
		z[2] = 0.;
		z[3] = 0.;
		z[4] = 0.;

		for (i = 0; i < 5; ++i) {

			z[0] += x[i] * b[i + 0];
			z[1] += x[i] * b[i + 5];
			z[2] += x[i] * b[i + 10];
			z[3] += x[i] * b[i + 15];
			z[4] += x[i] * b[i + 20];
		}

		v[0] = (lse_float_t) x[0];
		v[1] = (lse_float_t) x[1];
		v[2] = (lse_float_t) x[2];
		v[3] = (lse_float_t) x[3];
		v[4] = (lse_float_t) x[4];

		v[5] = (lse_float_t) z[0];
		v[6] = (lse_float_t) z[1];
		v[7] = (lse_float_t) z[2];
		v[8] = (lse_float_t) z[3];
		v[9] = (lse_float_t) z[4];

		lse_insert(&ls, v);
	}

	lse_solve(&ls);

	for (i = 0; i < 25; ++i) {

		rel[i] = fabs((ls.sol.m[i] - b[i]) / b[i]);
	}

	printf("\ncase  ---- random ----\n");

	printf("rel  % .8lE % .8lE % .8lE % .8lE % .8lE\n",
			rel[0], rel[1], rel[2], rel[3], rel[4]);
	printf("rel  % .8lE % .8lE % .8lE % .8lE % .8lE\n",
			rel[5], rel[6], rel[7], rel[8], rel[9]);
	printf("rel  % .8lE % .8lE % .8lE % .8lE % .8lE\n",
			rel[10], rel[11], rel[12], rel[13], rel[14]);
	printf("rel  % .8lE % .8lE % .8lE % .8lE % .8lE\n",
			rel[15], rel[16], rel[17], rel[18], rel[19]);
	printf("rel  % .8lE % .8lE % .8lE % .8lE % .8lE\n",
			rel[20], rel[21], rel[22], rel[23], rel[24]);

	lse_std(&ls);

	printf("std  % .8lE % .8lE % .8lE % .8lE % .8lE\n",
			ls.std.m[0], ls.std.m[1], ls.std.m[2],
			ls.std.m[3], ls.std.m[4]);

	lse_esv(&ls, 2);

	printf("cond % .8lE\n\n", ls.esv.max / ls.esv.min);
}

int main(int argc, char *argv[])
{
	lse_basic_test();
	lse_large_test(50);
	lse_large_test(200);
	lse_advanced_test();
	lse_mean_std_test();
	lse_random_test();

	return 0;
}

