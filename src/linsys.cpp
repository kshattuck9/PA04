#include "linsys.h"
#include <cmath>

void gauss_solve(int n, double* a, double* b)
{
    for (int k = 0; k < n - 1; ++k) {
        for (int i = k + 1; i < n; ++i) {
            double m = a[i*n+k] / a[k*n+k];
            a[i*n+k] = 0.0;
            for (int j = k + 1; j < n; ++j)
                a[i*n+j] -= m * a[k*n+j];
            b[i] -= m * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j)
            b[i] -= a[i*n+j] * b[j];
        b[i] /= a[i*n+i];
    }
}

void lu_factor(int n, double* a, int* p)
{
    for (int i = 0; i < n; ++i)
        p[i] = i;
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(a[j*n+i]) > std::abs(a[pivot*n+i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            for (int j = 0; j < n; ++j) {
                double temp = a[i*n+j];
                a[i*n+j] = a[pivot*n+j];
                a[pivot*n+j] = temp;
            }
            int temp = p[i];
            p[i] = p[pivot];
            p[pivot] = temp;
        }

        for (int j = i + 1; j < n; ++j) {
            a[j*n+i] /= a[i*n+i];
            for (int k = i + 1; k < n; ++k) {
                a[j*n+k] -= a[j*n+i] * a[i*n+k];
            }
        }
    }
    // TODO: select each pivot, exchange complete rows and p entries,
    // then store multipliers below the diagonal and update the trailing block.
    // Remove this placeholder when implementing the routine.
}

void forward_substitution(int n, const double* lu, const double* b, double* x)
{
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
        for (int j = 0; j < i; ++j)
            x[i] -= lu[i*n+j] * x[j];
        x[i] /= lu[i*n+i];
    }
    // TODO: traverse the rows in increasing order; L has unit diagonal.
}

void back_substitution(int n, const double* lu, const double* b, double* x)
{
    for (int i = n-1; i >= 0; --i){
        x[i] = b[i];
        for (int j = i+1; j < n; ++j)
            x[i] -= lu[i*n+j] * x[j];
        x[i] /= lu[i*n+i];
    }
    // TODO: traverse the rows in decreasing order.
}

void lu_solve(int n, const double* lu, const int* p, const double* b, double* x)
{
    double* y = new double[n];
    for (int i = 0; i < n; ++i)
        y[i] = b[p[i]];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j)
            y[i] -= lu[i*n+j] * y[j];
    }
    for (int i = n - 1; i >= 0; --i)
    {
        for (int j = i + 1; j < n; ++j)
            y[i] -= lu[i*n+j] * x[j];

        x[i] = y[i] / lu[i*n+i];
    }
    delete[] y;
    // TODO: permute the right-hand side, then use the two triangular solves.
}
