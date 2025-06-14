#include "lse.h"

#if LSE_FAST_GIVENS != 0
/* Define the maixmal allowed scale of the fast transformation. The input data
 * range is reduced by this number. Large number allows us to do scaling rarely.
 * */
#define LSE_DMAX 		((lse_float_t) 1048576)
#endif /* LSE_FAST_GIVENS */

#define LSE_RM_TOP(ls)		((ls)->rm + (ls)->n_cascades - 1)

/* Define built-in branch prediction functions.
 * */
#define likely(x)		__builtin_expect(!!(x), 1)
#define unlikely(x)		__builtin_expect(!!(x), 0)

/* Define which external math functions to use in LSE.
 * */
#define lse_fabsf(x)		__builtin_fabsf(x)
#define lse_sqrtf(x)		__builtin_sqrtf(x)

static void
#if LSE_FAST_GIVENS != 0
lse_qrupdate(lse_t *ls, lse_upper_t *rm, lse_float_t *xz, lse_float_t d0, int nz)
#else /* LSE_FAST_GIVENS */
lse_qrupdate(lse_t *ls, lse_upper_t *rm, lse_float_t *xz, int nz)
#endif
{
	lse_float_t	*m = rm->m;
#if LSE_FAST_GIVENS != 0
	lse_float_t	*d = rm->d;
#endif /* LSE_FAST_GIVENS */

	lse_float_t	x0, xi, alpa, beta;
	int		n, i, j;

#if LSE_FAST_GIVENS != 0
	lse_float_t	di;
#endif /* LSE_FAST_GIVENS */

	n = (rm->rows < rm->keep) ? rm->rows : rm->keep;

	/* Do we have leading zeros?
	 * */
	if (unlikely(nz > 0)) {

		m += nz * rm->len - nz * (nz - 1) / 2;
	}

	for (i = nz; i < n; ++i) {

		m += - i;

		x0 = - xz[i];
		xi = m[i];

		if (x0 != (lse_float_t) 0) {

#if LSE_FAST_GIVENS != 0
			di = d[i];

			/* We build the fast Givens transformation.
			 * */
			alpa = x0 * di;
			beta = xi * d0;

			if (x0 * alpa < xi * beta) {

				beta = - alpa / beta;
				alpa = x0 / xi;

				m[i] = m[i] + beta * xz[i];

				for (j = i + 1; j < rm->len; ++j) {

					xi = m[j] + beta * xz[j];
					x0 = alpa * m[j] + xz[j];

					xz[j] = x0;
					 m[j] = xi;
				}

				x0 = (lse_float_t) 1 - alpa * beta;

				d[i] = di * x0;
				  d0 = d0 * x0;
			}
			else {
				beta = - beta / alpa;
				alpa = xi / x0;

				m[i] = beta * m[i] + xz[i];

				for (j = i + 1; j < rm->len; ++j) {

					xi = beta * m[j] + xz[j];
					x0 = m[j] + alpa * xz[j];

					xz[j] = x0;
					 m[j] = xi;
				}

				x0 = (lse_float_t) 1 - alpa * beta;

				d[i] = d0 * x0;
				  d0 = di * x0;
			}

			/* Keep diagonal is in allowed range.
			 * */
			if (unlikely(d[i] > LSE_DMAX * LSE_DMAX)) {

				alpa = (lse_float_t) 1 / LSE_DMAX;

				for (j = i; j < rm->len; ++j) {

					x0 = m[j];
					m[j] = x0 * alpa;
				}

				d[i] *= (alpa * alpa);
			}

			if (unlikely(d0 > LSE_DMAX * LSE_DMAX)) {

				alpa = (lse_float_t) 1 / LSE_DMAX;

				for (j = i + 1; j < rm->len; ++j) {

					x0 = xz[j];
					xz[j] = x0 * alpa;
				}

				d0 *= (alpa * alpa);
			}
#else /* LSE_FAST_GIVENS */

			/* WARNING: We use naive hypot implementation as it is
			 * the fastest one and quite ulp-accurate.
			 */
			alpa = lse_sqrtf(x0 * x0 + xi * xi);
			beta = (lse_float_t) 1 / alpa;

			m[i] = alpa;

			/* We build the orthogonal transformation.
			 * */
			alpa = x0 * beta;
			beta = xi * beta;

			for (j = i + 1; j < rm->len; ++j) {

				xi = beta * m[j] - alpa * xz[j];
				x0 = alpa * m[j] + beta * xz[j];

				xz[j] = x0;
				 m[j] = xi;
			}
#endif /* LSE_FAST_GIVENS */
		}

		m += rm->len;
	}

	if (unlikely(n < rm->rows)) {

		m += - n;

		if (rm->lazy != 0) {

			/* We merge the retained row-vector into the upper
			 * cascade matrix before copying the new content.
			 * */
#if LSE_FAST_GIVENS != 0
			lse_qrupdate(ls, rm + 1, m, d[n], n);
#else /* LSE_FAST_GIVENS */
			lse_qrupdate(ls, rm + 1, m, n);
#endif
		}

		/* Copy the tail content.
		 * */
		for (i = n; i < rm->len; ++i)
			m[i] = xz[i];

#if LSE_FAST_GIVENS != 0
		d[n] = d0;
#endif /* LSE_FAST_GIVENS */
	}

	rm->keep += 1;

	if (unlikely(		rm->keep >= ls->n_threshold
				&& ls->n_cascades >= 2)) {

		if (rm < LSE_RM_TOP(ls)) {

			/* Mark the cascade matrix content as lazily merged.
			 * */
			rm->keep = 0;
			rm->lazy = 1;
		}
		else {
			/* Update the threshold value based on amount of data
			 * rows in top cascade.
			 * */
			ls->n_threshold = (rm->keep > ls->n_threshold)
				? rm->keep : ls->n_threshold;
		}
	}
}

static void
lse_qrmerge(lse_t *ls, lse_upper_t *rm, lse_upper_t *um)
{
	lse_float_t	*m = um->m;
#if LSE_FAST_GIVENS != 0
	lse_float_t	*d = um->d;
#endif /* LSE_FAST_GIVENS */

	int		n0, i;

	n0 = (um->lazy != 0) ? um->rows
		: (um->rows < um->keep) ? um->rows : um->keep;

	for (i = 0; i < n0; ++i) {

		m += - i;

		/* We extract one by one the row-vectors from cascade
		 * matrix and merge them into the upper cascade matrix.
		 * */
#if LSE_FAST_GIVENS != 0
		lse_qrupdate(ls, rm, m, d[i], i);
#else /* LSE_FAST_GIVENS */
		lse_qrupdate(ls, rm, m, i);
#endif

		m += um->len;
	}

	um->keep = 0;
	um->lazy = 0;
}

static void
lse_qrfinal(lse_t *ls)
{
	lse_upper_t	*rm = LSE_RM_TOP(ls);

	int		len, nul, i;

	for (i = 0; i < ls->n_cascades - 1; ++i) {

		/* We merge all cascades into the top \rm matrix.
		 * */
		lse_qrmerge(ls, &ls->rm[i + 1], &ls->rm[i]);
	}

	if (unlikely(rm->keep < rm->rows)) {

		/* Zero out uninitialized tail content.
		 * */
		len = rm->keep * rm->len - rm->keep * (rm->keep - 1) / 2;
		nul = rm->rows * rm->len - rm->rows * (rm->rows - 1) / 2;

		for (i = len; i < nul; ++i)
			rm->m[i] = (lse_float_t) 0;

#if LSE_FAST_GIVENS != 0
		for (i = rm->keep; i < rm->rows; ++i)
			rm->d[i] = (lse_float_t) 1;
#endif /* LSE_FAST_GIVENS */

		rm->keep = rm->rows;
	}
}

static void
lse_qrstep(lse_t *ls, lse_upper_t *um, lse_upper_t *im, lse_float_t *u)
{
	lse_float_t	*mq = im->m;
#if LSE_FAST_GIVENS != 0
	lse_float_t	*d = um->d;
#endif /* LSE_FAST_GIVENS */
	lse_float_t	*m;

	int		i, j;

	um->keep = 0;
	um->lazy = 0;

	/* Here we transpose the input matrix \im and bring it to the
	 * upper-triangular form again and store into \um.
	 * */
	for (i = 0; i < um->rows; ++i) {

		m = mq;

		for (j = 0; j < i + 1; ++j) {

			u[j] = m[0];
			m += im->len - (j + 1);
		}

		for (j = i + 1; j < um->len; ++j)
			u[j] = (lse_float_t) 0;

#if LSE_FAST_GIVENS != 0
		lse_qrupdate(ls, um, u, d[i], 0);
#else /* LSE_FAST_GIVENS */
		lse_qrupdate(ls, um, u, 0);
#endif

		mq += 1;
	}
}

int lse_getsize(int n_cascades, int n_full)
{
	int		n_lse, n_vm;

	n_lse = sizeof(lse_t) - sizeof(((lse_t *) 0)->vm);

	n_vm = n_cascades * n_full * (n_full + 1) / 2

#if LSE_FAST_GIVENS != 0
		  + n_cascades * n_full
#endif /* LSE_FAST_GIVENS */

		  + n_full * n_full / 4 + n_full / 2 + 1;

	return n_lse + sizeof(lse_float_t) * n_vm;
}

void lse_construct(lse_t *ls, int n_cascades, int n_len_of_x, int n_len_of_z)
{
	lse_float_t	*vm = ls->vm;

	int		n_full, i;

	ls->n_cascades = n_cascades;
	ls->n_len_of_x = n_len_of_x;
	ls->n_len_of_z = n_len_of_z;

	n_full = n_len_of_x + n_len_of_z;

	ls->n_threshold = n_full * 2;
	ls->n_total = 0;

	for (i = 0; i < ls->n_cascades; ++i) {

		ls->rm[i].len = n_full;
		ls->rm[i].rows = n_full;
		ls->rm[i].keep = 0;
		ls->rm[i].lazy = 0;
		ls->rm[i].m = vm;

		vm += n_full * (n_full + 1) / 2;

#if LSE_FAST_GIVENS != 0
		ls->rm[i].d = vm;
		vm += n_full;
#endif /* LSE_FAST_GIVENS */
	}

	ls->sol.len = ls->n_len_of_x * ls->n_len_of_z;
	ls->sol.m = vm;

	ls->std.len = ls->n_len_of_z;
	ls->std.m = vm + ls->sol.len;

	ls->esv.max = (lse_float_t) 0;
	ls->esv.min = (lse_float_t) 0;
}

void lse_nostd(lse_t *ls)
{
	int		n_full, i;

	/* Do not update lower triangle block.
	 * */
	n_full = ls->n_len_of_x;

	ls->n_threshold = n_full * 2;

	for (i = 0; i < ls->n_cascades; ++i) {

		ls->rm[i].rows = n_full;
	}
}

void lse_insert(lse_t *ls, lse_float_t *xz)
{
#if LSE_FAST_GIVENS != 0
	lse_qrupdate(ls, ls->rm, xz, (lse_float_t) 1, 0);
#else /* LSE_FAST_GIVENS */
	lse_qrupdate(ls, ls->rm, xz, 0);
#endif

	ls->n_total += 1;
}

void lse_ridge(lse_t *ls, lse_float_t la)
{
	lse_float_t	*xz = ls->sol.m;

	int		i, j;

	/* Add bias using the unity matrix multiplied by \la.
	 * */
	for (i = 0; i < ls->n_len_of_x; ++i) {

		xz[i] = la;

		for (j = i + 1; j < ls->rm[0].len; ++j)
			xz[j] = (lse_float_t) 0;

#if LSE_FAST_GIVENS != 0
		lse_qrupdate(ls, ls->rm, xz, (lse_float_t) 1, i);
#else /* LSE_FAST_GIVENS */
		lse_qrupdate(ls, ls->rm, xz, i);
#endif
	}
}

void lse_forget(lse_t *ls, lse_float_t la)
{
	lse_upper_t	*rm;

	int		n0, len, i, j;

	for (i = 0; i < ls->n_cascades; ++i) {

		rm = &ls->rm[i];

		n0 = (rm->lazy != 0) ? rm->rows
			: (rm->rows < rm->keep) ? rm->rows : rm->keep;

		if (n0 != 0) {

			len = n0 * rm->len - n0 * (n0 - 1) / 2;

			/* We just scale \rm matrices with factor \la.
			 * */
			for (j = 0; j < len; ++j)
				rm->m[j] *= la;
		}
	}
}

void lse_merge(lse_t *ls, lse_t *lb)
{
	lse_upper_t	*um = LSE_RM_TOP(lb);

	lse_float_t	*m = um->m;
#if LSE_FAST_GIVENS != 0
	lse_float_t	*d = um->d;
#endif /* LSE_FAST_GIVENS */

	int		i;

	lse_qrfinal(lb);

	for (i = 0; i < um->rows; ++i) {

		m += - i;

		/* TODO: Here you could shuffle \um columns.
		 * */

		/* We extract one by one the row-vectors from \lb instance and
		 * merge them into the \ls instance.
		 * */
#if LSE_FAST_GIVENS != 0
		lse_qrupdate(ls, ls->rm, m, d[i], i);
#else /* LSE_FAST_GIVENS */
		lse_qrupdate(ls, ls->rm, m, i);
#endif

		m += um->len;
	}
}

void lse_solve(lse_t *ls)
{
	lse_upper_t	*rm = LSE_RM_TOP(ls);

	lse_float_t	*sol = ls->sol.m;
	lse_float_t	*mq, *m, u;

	int		n, i, j;

	lse_qrfinal(ls);

	mq = rm->m + (ls->n_len_of_x - 1) * rm->len
		- ls->n_len_of_x * (ls->n_len_of_x - 1) / 2;

	/* We calculate solution \b with backward substitution.
	 * */
	for (n = 0; n < ls->n_len_of_z; ++n) {

		m = mq;

		for (i = ls->n_len_of_x - 1; i >= 0; --i) {

			u = (lse_float_t) 0;

			for (j = i + 1; j < ls->n_len_of_x; ++j)
				u += sol[j] * m[j];

			sol[i] = (m[ls->n_len_of_x + n] - u) / m[i];

			m += i - rm->len;
		}

		sol += ls->n_len_of_x;
	}
}

void lse_std(lse_t *ls)
{
	lse_upper_t	*rm = LSE_RM_TOP(ls);

	lse_float_t	*std = ls->std.m;
	lse_float_t	*mq, *m, u, ratio;

#if LSE_FAST_GIVENS != 0
	lse_float_t	*d = rm->d + ls->n_len_of_x;
#endif /* LSE_FAST_GIVENS */

	int		i, j;

	lse_qrfinal(ls);

	mq = rm->m + ls->n_len_of_x * rm->len
		- ls->n_len_of_x * (ls->n_len_of_x - 1) / 2;

	ratio = (lse_float_t) 1 / (lse_float_t) (ls->n_total - 1);

	/* We calculate l2 norm over \rm columns.
	 * */
	for (i = 0; i < ls->n_len_of_z; ++i) {

		m = mq;

#if LSE_FAST_GIVENS != 0
		u = m[0] * m[0] / d[0];
#else /* LSE_FAST_GIVENS */
		u = m[0] * m[0];
#endif

		for (j = 1; j < i + 1; ++j) {

			m += rm->len - (ls->n_len_of_x + j);

#if LSE_FAST_GIVENS != 0
			u += m[0] * m[0] / d[j];
#else /* LSE_FAST_GIVENS */
			u += m[0] * m[0];
#endif
		}

		std[i] = lse_sqrtf(u * ratio);

		mq += 1;
	}
}

void lse_esv(lse_t *ls, int n_approx)
{
	lse_upper_t	um, im, *rm = LSE_RM_TOP(ls);
	lse_float_t	*m, u;

	int		len, i;

	lse_qrfinal(ls);

	len = ls->n_len_of_x * (ls->n_len_of_x + 1) + ls->n_len_of_x * 3;

	if (ls->rm[0].m + len <= rm->m) {

		/* We allocate temporal \um matrices instead of \rm
		 * cascades that are empty after merge.
		 * */
		m = ls->rm[0].m;
	}
	else {
		/* WARNING: We allocate temporal \um matrices in tail
		 * of LSE memory instead of \b and so on.
		 * */
		m = ls->sol.m;
	}

	um.len = ls->n_len_of_x;
	um.rows = ls->n_len_of_x;
	um.m = m;

	m += ls->n_len_of_x * (ls->n_len_of_x + 1) / 2;

#if LSE_FAST_GIVENS != 0
	um.d = m;
	m += ls->n_len_of_x;
#endif /* LSE_FAST_GIVENS */

	im.len = ls->n_len_of_x;
	im.rows = ls->n_len_of_x;
	im.m = m;

	m += ls->n_len_of_x * (ls->n_len_of_x + 1) / 2;

#if LSE_FAST_GIVENS != 0
	im.d = m;
	m += ls->n_len_of_x;
#endif /* LSE_FAST_GIVENS */

#if LSE_FAST_GIVENS != 0
	for (i = 0; i < ls->n_len_of_x; ++i) {

		um.d[i] = (lse_float_t) 1;
		im.d[i] = rm->d[i];
	}
#endif /* LSE_FAST_GIVENS */

	/* First step of QR algorithm.
	 * */
	lse_qrstep(ls, &um, rm, m);

	for (i = 1; i < n_approx; ++i) {

		/* Swap the matrices content.
		 * */
		{ lse_upper_t qm = um; um = im; im = qm; }

		/* We run the reduced form of QR algorithm. With each
		 * iteration off-diagonal elements tend to zero so the
		 * diagonal approaches singular values.
		 * */
		lse_qrstep(ls, &um, &im, m);
	}

	m = um.m;

	/* We are looking for the largest and smallest diagonal elements of \um.
	 * */
	for (i = 0; i < ls->n_len_of_x; ++i) {

#if LSE_FAST_GIVENS != 0
		u = lse_fabsf(m[0] / lse_sqrtf(um.d[i] * im.d[i]));
#else /* LSE_FAST_GIVENS */
		u = lse_fabsf(m[0]);
#endif

		if (likely(i != 0)) {

			ls->esv.max = (ls->esv.max < u) ? u : ls->esv.max;
			ls->esv.min = (ls->esv.min > u) ? u : ls->esv.min;
		}
		else {
			ls->esv.max = u;
			ls->esv.min = u;
		}

		m += um.len - i;
	}
}

