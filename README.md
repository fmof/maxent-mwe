# GSL LBFGS-2 Interaction with a Simple Maxent Model

This repository documents an issue with (a portion of) GSL's multimin library: specifically, lbfgs2 does not properly optimize a convex function, even though it returns successfully (GSL_SUCCESS). 
The convex function is the (negative) log-likelihood of a maximum entropy model on a simple test case.

This has been tested with GSL 1.16.

## How-to

### Getting the code

```
$ git clone https://github.com/fmof/maxent-mwe
$ cd maxent-mwe
```

### Dependencies

* An active internet connection (to download the correct version of GSL).
* None, except for the standard math library.
* A C99 compiler (tested with GCC 4.8.2).

GSL (version 1.16) is pulled down and installed locally to a subdirectory of this repository. 
All GSL linking is done statically.

### Compile:

* `make all`
This will:
  * Perform a local install to `$(pwd)/gsl_1.16_build`.
  * Build the executable with `maxent_mwe` with `-std=c99 -O2 -g`.

#### Changing options
* Turn on verbose output by adding `-DVERBOSE` to `CFLAGS`.

## What The Program Does
The program `maxent_mwe` optimizes a very simple [log-linear (/maximum entropy/logistic regression) model](http://en.wikipedia.org/wiki/Log-linear_model).
Given a support `S`, a vector `f` of `F` separate features, and an `F`-length vector of real-valued weights θ, this defines a probability distribution

```
p(x) = 1/Z * exp(θ \cdot f(x))
```

where `Z` is the partition function, `sum_{y \in S} exp(θ \cdot f(y))`.

This program sets `S = {0, 1}` and defines two binary features `f_0(x)` and `f_1(x)`, by `f_i(x) = (i == x)`. 
Therefore, `p(x) = 1/Z * exp(θ_x)`, and `Z = exp(θ_0) + exp(θ_1)`.

The program observes a simple dataset, where "0" appears twice and "1" appears once. 
It fits \hat{θ} to these data (maximum likelihood estimation) via the [GSL BFGS-2](https://www.gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-with-Derivatives.html) optimization method. 
The objective function `L(θ)` is the (negative) log-likelihood of the observed data. 
Using `c(x)` to indicate the number of times `x` is observed in the data, we can write

```
L(θ) = -( c(0) * log p(0 | θ) + c(1) * log p(1 | θ) )
     = -( 2 * log p(0 | θ) + 1 * log p(1 | θ) ).
```

The gradient `G(θ)` of `L` follows the "observed minus expected" formulation of maxent models:

```
G(θ) = < -c(0) + 3 * p(0 | θ),
         -c(1) + 3 * p(1 | θ) >
     = < -2 + 3 * p(0 | θ),
         -1 + 3 * p(1 | θ) >.
```

These functions are implemented in `src/maxent_mwe.c#neg_ll_eval` (`L(θ)`) and `src/maxent_mwe.c#neg_ll_grad` (`G(θ)`).

## Output: Actual vs. Expected

The program should properly fit a model to the provided data (described above). 
The first three lines output demonstrate some tests for computing the partition function, the objective and the gradient. 
The output **should** be:

```
$ ./maxent_mwe 
(GOOD) Z(3.200000, -2.000000) = 24.667865, should be 24.667865 (tolerance 1.000000e-06)
(GOOD) nLL(3.200000, -2.000000) = 5.216504, should be 5.216504 (tolerance 1.000000e-06)
(GOOD) nGrad(3.200000, -2.000000) = (0.983541, -0.983541), should be (0.983541, -0.983541) (tolerance 1.000000e-06)
(GOOD) LBFGS2 optimization status = 0, should be 0
(GOOD) point[0] = 0.346574, should be 0.346574 (tolerance = 1.000000e-06)
(GOOD) point[0] = -0.346574, should be -0.346574 (tolerance = 1.000000e-06)
(GOOD) Ending status agrees with optimization result
```

However, it is actually

```
$ ./maxent_mwe 
(GOOD) Z(3.200000, -2.000000) = 24.667865, should be 24.667865 (tolerance 1.000000e-06)
(GOOD) nLL(3.200000, -2.000000) = 5.216504, should be 5.216504 (tolerance 1.000000e-06)
(GOOD) nGrad(3.200000, -2.000000) = (0.983541, -0.983541), should be (0.983541, -0.983541) (tolerance 1.000000e-06)
(GOOD) LBFGS2 optimization status = 0, should be 0
(BAD) point[0] = 0.000000, should be 0.346574 (tolerance = 1.000000e-06)
(BAD) point[0] = 0.000000, should be -0.346574 (tolerance = 1.000000e-06)
(BAD) Ending status does not agree with optimization result
```

This indicates an error with the optimization function. 

By compiling with verbosity on (`CFLAGS+=-DVERBOSE`), you can see intermediate function and gradient evaluations. 
Specifically, intermediate `θ` values are tried that minimize the objective. 
This suggests an issue with the line search. 
Nevertheless, a poorly terminating line search (and its containing caller) should not return a "success."

```
    p = (0.0000000000, 0.0000000000) ==eval==> 2.079442
    grad|_(0.0000000000, 0.0000000000) ==grad==> (-0.500000, 0.500000)
    grad|_(0.0000000000, 0.0000000000) ==grad==> (-0.500000, 0.500000)
    p = (0.0070710678, -0.0070710678) ==eval==> 2.072445
    grad|_(0.0070710678, -0.0070710678) ==grad==> (-0.489394, 0.489394)
    grad|_(0.0070710678, -0.0070710678) ==grad==> (-0.489394, 0.489394)
    p = (0.0707106781, -0.0707106781) ==eval==> 2.016225
    grad|_(0.0707106781, -0.0707106781) ==grad==> (-0.394110, 0.394110)
    grad|_(0.0707106781, -0.0707106781) ==grad==> (-0.394110, 0.394110)
    p = (0.3373587446, -0.3373587446) ==eval==> 1.909656
    grad|_(0.3373587446, -0.3373587446) ==grad==> (-0.012324, 0.012324)
    grad|_(0.3373587446, -0.3373587446) ==grad==> (-0.012324, 0.012324)
    ...
```
