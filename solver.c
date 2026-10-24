#include "solver.h"
#include <math.h>
#include <stdlib.h>

int solve_tridiagonal(const double *a, const double *b, const double *c,
                      const double *d, double *x, int n)
{
    if (n <= 0)
        return -1;

    if (n == 1)
    {
        if (fabs(b[0]) < 1e-20)
            return -1;
        x[0] = d[0] / b[0];
        return 0;
    }

    double *cp = (double *)malloc(n * sizeof(double));
    double *dp = (double *)malloc(n * sizeof(double));

    if (!cp || !dp)
    {
        free(cp);
        free(dp);
        return -1;
    }

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; i++)
    {
        double denom = b[i] - a[i] * cp[i - 1];
        if (fabs(denom) < 1e-20)
        {
            free(cp);
            free(dp);
            return -1;
        }
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    x[n - 1] = dp[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    free(cp);
    free(dp);

    return 0;
}