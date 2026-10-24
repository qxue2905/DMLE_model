#ifndef SOLVER_H
#define SOLVER_H

int solve_tridiagonal(const double *a, const double *b, const double *c,
                      const double *d, double *x, int n);

#endif /* SOLVER_H */