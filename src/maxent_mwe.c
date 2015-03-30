#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_log.h"
#include "stdio.h"

typedef struct {
  int* counts;
  int total_count;
} UnigramClosure;

double compute_partition(const gsl_vector* log_weights) {
  double partition = 0.0;
  const size_t data_size = log_weights->size;
  for(size_t i = 0; i < data_size; ++i) {
    partition += gsl_sf_exp(gsl_vector_get(log_weights, i));
  }
  return partition;
}

double neg_ll_eval(const gsl_vector* log_weights, void *fparams) {
  UnigramClosure* closure = (UnigramClosure*)fparams;
  const double log_partition = gsl_sf_log(compute_partition(log_weights));
  double ll = 0.0;
  const size_t data_size = log_weights->size;
  for(size_t i = 0; i < data_size; ++i) {
    ll += ((double)closure->counts[i]) * (gsl_vector_get(log_weights,i) - log_partition);
  }
#ifdef VERBOSE
  printf("    p = (%f, %f) ==eval==> %f\n",
	 gsl_vector_get(log_weights, 0),
	 gsl_vector_get(log_weights, 1),
	 -ll);
#endif
  return -ll;
}
void neg_ll_grad(const gsl_vector* log_weights, void *fparams, gsl_vector* gsl_grad) {
  UnigramClosure* closure = (UnigramClosure*)fparams;
  const double partition = compute_partition(log_weights);
  const size_t data_size = gsl_grad->size;
  const double N = (double)(closure->total_count);
  for(size_t i = 0; i < data_size; ++i) {
    double prob = gsl_sf_exp(gsl_vector_get(log_weights, i)) / partition;
    double val = closure->counts[i] - N*prob;
    gsl_vector_set(gsl_grad, i, -val);
#ifdef VERBOSE
    printf("    grad|_(%f, %f) ==grad==> (%f, %f)\n",
	   gsl_vector_get(log_weights,0),
	   gsl_vector_get(log_weights,1),
	   gsl_vector_get(gsl_grad,0),
	   gsl_vector_get(gsl_grad,1));
#endif
  }
}
void neg_ll_eg(const gsl_vector* trial_weights, void *fparams,
	       double* f, gsl_vector* grad) {
  *f = neg_ll_eval(trial_weights, fparams);
  neg_ll_grad(trial_weights, fparams, grad);
}

gsl_multimin_function_fdf* get_fdf(UnigramClosure* params, const size_t size) {
  gsl_multimin_function_fdf *objective = malloc(sizeof(gsl_multimin_function_fdf));
  objective->n   = size;
  objective->f   = &neg_ll_eval;
  objective->df  = &neg_ll_grad;
  objective->fdf = &neg_ll_eg;
  objective->params = (void*)params;
  return objective;
}

// return a zero-vector of type gsl_vector* with dimension == size_
gsl_vector* get_optimization_initial_point(int size_) {
  gsl_vector *vec = gsl_vector_alloc(size_);
  for(int i = 0; i < size_; ++i)
    gsl_vector_set(vec, i, 0.0);
  return vec;
}

UnigramClosure* get_model() {
  UnigramClosure* closure = malloc(sizeof(UnigramClosure));
  int *counts = malloc(2 * sizeof(int));
  counts[0] = 2;
  counts[1] = 1;
  closure->counts = counts;
  closure->total_count = 0;
  for(size_t i = 0; i < 2; ++i) {
    closure->total_count += counts[i];
  }
  return closure;
}

int main() {
  UnigramClosure* closure = get_model();
  size_t iter = 0;
  int status;
  const int arity = 2;

  gsl_multimin_function_fdf* my_func = get_fdf(closure, arity);

  gsl_vector *point = get_optimization_initial_point(arity);

  char *is_good[] = {"BAD", "GOOD"};

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
  gsl_multimin_fdfminimizer *minimizer = gsl_multimin_fdfminimizer_alloc (T, arity);
  gsl_multimin_fdfminimizer_set(minimizer, my_func, point, 0.01, 1e-4);

  do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(minimizer);
      if (status)
        break;
      status = gsl_multimin_test_gradient(minimizer->gradient, 1e-3);
  } while (status == GSL_CONTINUE && iter < 100);

  printf("(%s) LBFGS2 optimization status = %d, should be %d\n", is_good[status == GSL_SUCCESS], status, GSL_SUCCESS);
  double expected_point[] = {0.346574, -0.346574};
  double eps = 1E-6;
  int good = 1;
  for(int i = 0; i < arity; ++i) {
    int within = fabs(gsl_vector_get(point, i) - expected_point[i]) < eps;
    good &= within;
    printf("(%s) point[0] = %.6f, should be %.6f (tolerance = %e)\n", 
	   is_good[within], 
	   gsl_vector_get(point, i),
	   expected_point[i], eps);
  }
  if(good) {
    printf("(%s) Ending status agrees with optimization result\n", is_good[good]);
  } else {
    printf("(%s) Ending status does not agree with optimization result\n", is_good[good]);
  }

  free(my_func);
  free(closure->counts);
  free(closure);
  gsl_multimin_fdfminimizer_free(minimizer);
  gsl_vector_free(point);
  return !good;
}
