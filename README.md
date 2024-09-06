# Least Squares Estimate

LSE is a linear least squares solver library written in C.

## Brief

LSE is aimed to be used in signal processing on small embedded systems with
single precision FPU. The main goals are low precision losses on large datasets
and low memory footprint.

Typical usage scenario is about ~10 unknown parameters and over ~1kk of
measurements.

## Features

* Ordinary least squares (LS) solver using QR transformation.
* Low precision losses on large datasets.
* Able to solve rank deficient ill-posed problem using l2 regularization.
* Standard deviation (STD) estimate.
* Largest and smallest singular values approximation.
* Static and dynamic memory allocation.

## Approach

We will give the formulas in matlab-like syntax to make them readable in plain
text and monospace font. Let us define linear model equation.

	x(i) * b = z(i)

Where \x(i) and \z(i) is some known dataset of row-vectors, \b is unknown
matrix to be found. We need to find \b that minimizes the sum of the squared
residuals over all measurements.

	                         2
	sum || x(i) * b - z(i) ||  -->  min
	 i

We use a well-known method based on QR decomposition.

	X = [x(1); x(2); x(3); ...]
	Z = [z(1); z(2); z(3); ...]

	[Q,rm] = qr([X Z], 0)

	     [ RX  S  ]
	rm = [     RZ ]

	b = RX \ S

In this formulation we do not need the matrix \Q at all. And the matrix \rm can
be obtained by the sequential QR updating.

	rm = qrupdate([rm; [x(i) z(i)]])

Thus we do not need to store all dataset in a tall skinny matrix. But that is
not all. On large dataset we loss precision more and more due to roundoff
errors at each QR update. To overcome this we propose the cascading update.

Keep aggregating data in \rm(0) until number of aggregated rows reaches some
threshold number \n. After that merge \rm(0) into the next level matrix \rm(1).

	rm(1) = qrupdate([rm(1); rm(0)])

Then reset \rm(0) and continue to aggregate data in \rm(0). When number of
aggregated rows in \rm(i) reaches some threshold number \n we merge \rm(i) into
the next level matrix \rm(i+1).

The number of cascades is selected depending on desired error compensation and
available amount of memory. The threshold number \n increases each time when
top-level \rm(i) matrix is updated to produce balanced cascades.

To get the final solution we merge all cascade \rm(i) matrices into the one.

	rm = qrupdate([rm(0); rm(1); rm(2); ...])

We use fast Givens transformation to implement QR update procedure. Thus the
square root function is not used to find the solution.

## Example

See the API description in LSE header file and testbench code.

```
	lse_t       ls;
	
	lse_construct(&ls, LSE_CASCADE_MAX, 2, 1);
	
	for (...) {
	
		v[0] = (lse_float_t) x[0];
		v[1] = (lse_float_t) x[1];
		v[2] = (lse_float_t) z[0];
	
		lse_insert(&ls, v);
	}
	
	lse_solve(&ls);
	
	printf("sol %.4f %.4f\n", ls.sol.m[0], ls.sol.m[1]);
```

## Status

We are using LSE in motor control software on Cortex-M4F embedded processor.

