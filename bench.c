#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "lse.h"
#include "lfg.h"

#define LSE_EPSF		4E-7

#define LSE_printf(s)		fprintf(stderr, "%s in %s:%i\n", (s), __FILE__, __LINE__)
#define LSE_assert(x)		if ((x) == 0) { LSE_printf(#x); exit(-1); }

#define LSE_assert_relative(x, r)	LSE_assert(fabs((x) - (r)) < LSE_EPSF * (r))
#define LSE_assert_absolute(x, r)	LSE_assert(fabs((x) - (r)) < LSE_EPSF)
#define LSE_assert_approx(x, r, a)	LSE_assert(fabs((x) - (r)) < LSE_EPSF * (a))

void lse_basic_test()
{
	lse_t		ls;
	lse_float_t	v[LSE_FULL_MAX];

	/* ------------------- TEST ------------------- */

	lse_construct(&ls, 1, 4, 1);

	v[0] = (lse_float_t) 0.;
	v[1] = (lse_float_t) - 3.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.;
	v[4] = (lse_float_t) 2.;

	lse_reduce(&ls, v);

	v[0] = (lse_float_t) 2.;
	v[1] = (lse_float_t) 0.;
	v[2] = (lse_float_t) - 1.;
	v[3] = (lse_float_t) 1.;
	v[4] = (lse_float_t) 3.;

	lse_reduce(&ls, v);

	v[0] = (lse_float_t) - 7.;
	v[1] = (lse_float_t) 2.;
	v[2] = (lse_float_t) - 1.;
	v[3] = (lse_float_t) 0.;
	v[4] = (lse_float_t) - 6.;

	lse_reduce(&ls, v);

	v[0] = (lse_float_t) 1.;
	v[1] = (lse_float_t) - 1.;
	v[2] = (lse_float_t) 0.;
	v[3] = (lse_float_t) 2.;
	v[4] = (lse_float_t) 7.;

	lse_reduce(&ls, v);
	lse_solve(&ls);

	LSE_assert_relative(ls.sol.m[0], 1.);
	LSE_assert_relative(ls.sol.m[1], 2.);
	LSE_assert_relative(ls.sol.m[2], 3.);
	LSE_assert_relative(ls.sol.m[3], 4.);

	/* ------------------- TEST ------------------- */

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

	/* ------------------- TEST ------------------- */

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

	/* ------------------- TEST ------------------- */

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

	/* ------------------- TEST ------------------- */

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

	lse_cond(&ls, 4);

	lse_ridge(&ls, ls.n_len_of_x * ls.l_max * LSE_EPSF);
	lse_solve(&ls);

	LSE_assert_approx(ls.sol.m[0], 1., 10.);
	LSE_assert_approx(ls.sol.m[1], 1., 10.);
	LSE_assert_approx(ls.sol.m[2], 1., 10.);

	/* ------------------- TEST ------------------- */

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

	for (n = 0; n < n_full - 1; ++n) {

		for (i = 0; i < n; ++i)
			v[i] = (lse_float_t) 1.;

		v[n] = (lse_float_t) 2.;

		for (i = n + 1; i < n_full; ++i)
			v[i] = (lse_float_t) 1.;

		lse_insert(ls, v);
	}

	lse_cond(ls, 4);

	LSE_assert_approx(ls->l_max, n_full, 100. * n_full);
	LSE_assert_approx(ls->l_min, 1., 100.);

	lse_solve(ls);

	for (i = 0; i < ls->n_len_of_x; ++i) {

		LSE_assert_approx(ls->sol.m[i], 1. / (double) n_full, 10.);
	}

	free(ls);
	free(v);
}

void lse_advanced_test()
{
	lse_t		ls;
	lse_float_t	v[LSE_FULL_MAX];

	int		i;

	/* ------------------- TEST ------------------- */

	lse_construct(&ls, LSE_CASCADE_MAX, 1, 2);

	for (i = 0; i < 1000; ++i) {

		v[0] = (lse_float_t) 1.;
		v[1] = (lse_float_t) - i;
		v[2] = (lse_float_t) i;

		lse_insert(&ls, v);

		v[0] = (lse_float_t) 1.;
		v[1] = (lse_float_t) i;
		v[2] = (lse_float_t) - i;

		lse_insert(&ls, v);
	}

	lse_solve(&ls);

	LSE_assert_approx(ls.sol.m[0], 0., 10.);
	LSE_assert_approx(ls.sol.m[1], 0., 10.);

	/* ------------------- TEST ------------------- */

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

	lfg_start(20);

	lse_construct(&ls, LSE_CASCADE_MAX, 1, 4);

	for (i = 0; i < 1000000; ++i) {

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

	LSE_assert_approx(ls.sol.m[0], 11., .1 / LSE_EPSF);
	LSE_assert_approx(ls.sol.m[1], - 3., .1 / LSE_EPSF);
	LSE_assert_approx(ls.sol.m[2], 5., .1 / LSE_EPSF);
	LSE_assert_approx(ls.sol.m[3], - 4., .1 / LSE_EPSF);

	LSE_assert_approx(ls.std.m[0], 2., .1 / LSE_EPSF);
	LSE_assert_approx(ls.std.m[1], 1., .1 / LSE_EPSF);
	LSE_assert_approx(ls.std.m[2], 7., .1 / LSE_EPSF);
	LSE_assert_approx(ls.std.m[3], 5., .1 / LSE_EPSF);
}

void lse_bench()
{
	lse_t		ls;
	lse_float_t	v[LSE_FULL_MAX];

	double		x[3], z[3], b[9], rel[9];
	int		i;

	lfg_start(27);

	b[0] = lfg_gauss();
	b[1] = lfg_gauss();
	b[2] = lfg_gauss();
	b[3] = lfg_gauss();
	b[4] = lfg_gauss();
	b[5] = lfg_gauss();
	b[6] = lfg_gauss();
	b[7] = lfg_gauss();
	b[8] = lfg_gauss();

	lse_construct(&ls, LSE_CASCADE_MAX, 3, 3);

	for (i = 0; i < 5000000; ++i) {

		x[0] = 0. + lfg_gauss() * 1.;
		x[1] = 0. + lfg_gauss() * 1.;
		x[2] = 0. + lfg_gauss() * 1.;

		z[0] = x[0] * b[0] + x[1] * b[1] + x[2] * b[2] + lfg_gauss() * 0.;
		z[1] = x[0] * b[3] + x[1] * b[4] + x[2] * b[5] + lfg_gauss() * 0.;
		z[2] = x[0] * b[6] + x[1] * b[7] + x[2] * b[8] + lfg_gauss() * 0.;

		v[0] = (lse_float_t) x[0];
		v[1] = (lse_float_t) x[1];
		v[2] = (lse_float_t) x[2];

		v[3] = (lse_float_t) z[0];
		v[4] = (lse_float_t) z[1];
		v[5] = (lse_float_t) z[2];

		lse_insert(&ls, v);
	}

	lse_solve(&ls);

	rel[0] = fabs((ls.sol.m[0] - b[0]) / b[0]);
	rel[1] = fabs((ls.sol.m[1] - b[1]) / b[1]);
	rel[2] = fabs((ls.sol.m[2] - b[2]) / b[2]);
	rel[3] = fabs((ls.sol.m[3] - b[3]) / b[3]);
	rel[4] = fabs((ls.sol.m[4] - b[4]) / b[4]);
	rel[5] = fabs((ls.sol.m[5] - b[5]) / b[5]);
	rel[6] = fabs((ls.sol.m[6] - b[6]) / b[6]);
	rel[7] = fabs((ls.sol.m[7] - b[7]) / b[7]);
	rel[8] = fabs((ls.sol.m[8] - b[8]) / b[8]);

	printf("rel  % .8lE % .8lE % .8lE\n", rel[0], rel[3], rel[6]);
	printf("rel  % .8lE % .8lE % .8lE\n", rel[1], rel[4], rel[7]);
	printf("rel  % .8lE % .8lE % .8lE\n", rel[2], rel[5], rel[8]);

	lse_std(&ls);
	lse_cond(&ls, 4);

	printf("std  % .8lE % .8lE % .8lE\n", ls.std.m[0], ls.std.m[1], ls.std.m[2]);
	printf("cond % .8lE\n", ls.l_max / ls.l_min);
}

int main(int argc, char *argv[])
{
	lse_basic_test();
	lse_large_test(50);
	lse_large_test(200);
	lse_advanced_test();
	lse_mean_std_test();

	lse_bench();

	return 0;
}

