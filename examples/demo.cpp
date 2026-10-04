#include "linsys.h"
#include <iostream>

int main()
{
    // The two-row-exchange example in PA04 Written.
    double a[] = {0, 2, 1,  1, 1, 0,  2, 0, 1};
    double b[] = {7, 3, 5};
    double c[] = {3, 2, 3};
    double x[3] = {};
    int p[3];

    lu_factor(3, a, p);
    lu_solve(3, a, p, b, x);
    std::cout << "First solution (expected 1 2 3): ";
    for (int i = 0; i < 3; ++i)
        std::cout << x[i] << ' ';
    std::cout << '\n';

    lu_solve(3, a, p, c, x);
    std::cout << "Second solution (expected 1 1 1): ";
    for (int i = 0; i < 3; ++i)
        std::cout << x[i] << ' ';
    std::cout << '\n';

    //adding another case for practice 
    double d[] = {2, 1, 1, 4, 3, 3, 8, 7, 9};
    double e[] = {4, 10, 24};
    double y[3] = {};
    int q[3];

    lu_factor(3, d, q);
    lu_solve(3, d, q, e, y);

    std::cout << "Third solution (expected 1 1 1): ";
    for (int i = 0; i < 3; ++i)
        std::cout << y[i] << ' ';
    std::cout << '\n';

    // Case 4 testing a 4x4 matrix
double f[] = {
    1, 2, 0, 1,
    2, 1, 1, 0,
    0, 1, 2, 1,
    1, 0, 1, 2
};

double g[] = {9, 7, 12, 12};
double z[4] = {};
int r[4];

lu_factor(4, f, r);
lu_solve(4, f, r, g, z);

std::cout << "Fourth solution (expected 1 2 3 4): ";
for (int i = 0; i < 4; ++i)
    std::cout << z[i] << ' ';
std::cout << '\n';
}
