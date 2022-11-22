
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include "excerpt.h"

typedef float fp_t;
const double PI = acos(-1.0); //exact calc of pi

using namespace std;

//TRIGONOMETRIC SOLUTION LAPAZ 1949
template<typename fp_t>
void trigonometric(fp_t a, fp_t b, fp_t c, fp_t *x0, fp_t *x1) {
    //init
    fp_t a, a, a, TAU;

    cout << endl << A << "x^2 + " << B << "x + " << C << endl; // ax^2+bx+c
    fp_t arg1 = 2 * sqrt(abs(A * C)) / B; //arg at tetta
    fp_t arg2 = sqrt(abs(C / A)); //arg at x
    fp_t tetta_p_0;
    fp_t tetta_p_1;
    fp_t tetta_n_0;
    fp_t tetta_n_1;
    fp_t tetta0;
    if (C > 0 and abs(arg1) <= 1) {
        TAU = asin(-arg1); //support var
        tetta_p_0 = TAU;
        tetta_p_1 = TAU + PI;

        (0 < abs(tetta_p_0) <= PI/2) ? tetta0 = tetta_p_0 : tetta0 = tetta_p_1; // 0<|tetta_p|<pi/2

        *x0 = arg2 * tan(tetta0/2);
        *x1 = arg2 * 1/tan(tetta0/2);
    } else if (C <= 0) {
        TAU = atan(arg1);
        tetta_n_0 = TAU;
        tetta_n_1 = TAU + PI/2;
        (0 < abs(tetta_n_0) <= PI/2) ? tetta0 = tetta_n_0 : tetta0 = tetta_n_1; // 0<|tetta_n|<pi/2

        *x0 = arg2 * tan(tetta0/2);
        *x1 = -arg2 * 1/tan(tetta0/2);
    } else
        cout << "Only complex roots ";
    cout << "Roots: " << x0 << "; " << x1 << endl;
}
template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t x0, x1, max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, 1e-5, -1, 1, roots, coefficients);
    trigonometric<fp_t>(coefficients[2], coefficients[1], coefficients[0], &x0, &x1);
    vector<fp_t> roots_computed = {x1, x0};
    auto result = compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, max_absolute_error, max_relative_error);
    switch (result) {
        case PR_2_INFINITE_ROOTS:
            cout << "INFINITE ROOTS";
            break;
        case PR_AT_LEAST_ONE_ROOT_IS_FAKE:
            cout << "AT LEAST ONE ROOT IS FAKE";
            break;
        case PR_AT_LEAST_ONE_ROOT_LOST:
            cout << "AT LEAST ONE ROOT LOST";
            break;
        default:
            break;
    }
    return max_absolute_error;
}

typedef float fp_t;

int main() {
    float deviation, max_deviation = 0;
    for (auto i = 0; i < 10'000; ++i) {
        deviation = testPolynomial<float>(2);
        if (deviation != std::numeric_limits<float>::infinity()) {
            cout << "deviation = " << deviation << endl;
            if (deviation > max_deviation) {
                max_deviation = deviation;
            }
        }
        else cout << "\t\tComplex roots!"<<endl;
    }
    cout<< endl<<"MAX_deviation = "<< max_deviation<<endl;

}
