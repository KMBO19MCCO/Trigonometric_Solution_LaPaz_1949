// КМБО-03-19 Успенский Артём / Кукулиев Андрей
// почта: uspenskiy.artem@yandex.ru / kukuliew@mail.ru
// статья: https://articles.adsabs.harvard.edu/pdf/1949PASP...61..222L


#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "excerpt.h"
#define MAX_DISTANCE 1e-5

using namespace std;


//TRIGONOMETRIC SOLUTION LAPAZ 1949
template<typename fp_t>
int trigonometric(vector <fp_t> coefficients, vector <fp_t> &roots) {
    //init
    static const double PI = numbers::pi_v<fp_t>;             // Точное пи
    static const double PId_2 = static_cast<fp_t>(0.5) * PI;  // PI/2
    fp_t a, b, c;

    a = coefficients[2];
    b = coefficients[1];
    c = coefficients[0];

    //cout << endl << a << "x^2 + " << b << "x + " << c << endl; // ax^2+bx+c
    fp_t arg1;    //аргумент при tetta
    fp_t arg2;    //аргумент при x
    fp_t tetta_p; //пер. при c>0
    fp_t tetta_n; //пер. при c<0
    fp_t tetta_0; //главная пер.

    int j = 0;
    int cnt_roots = 0;

    // переменные для одинаковых проверок
    fp_t c_b = c / b;
    fp_t c_a = c / a;
    fp_t a_b = a / b;

    if (!isinf(a_b) && !isinf(c_a) && !isinf(c_b) && a!=0 && b!=0) {                // Проверка отсутствия линейности
        arg1 = 2 * sqrt(abs(a * c)) / b;
        arg2 = sqrt(abs(c_a));

        if (c > 0 and abs(arg1) <= 1) {              // Проверка на вещественные корни и "С" больше 0
            tetta_p = asin(-arg1);
            tetta_0 = tetta_p;

            while (!(0 < abs(tetta_0) < PId_2)) {           // Диапазон (0,PI/2)
                // (tetta_p - j * PI) / pow((-1), j)
                fp_t arg3 = fma(-PI, j, tetta_p);
                tetta_0 = fma(-2 * (j % 2 != 0), arg3, arg3); // tetta_0 = arg3 - 2 * (j % 2 != 0) * arg3

                (tetta_0 > PId_2) ? j++ : j--;              // Если выходит справа, иначе слева
            }

            fp_t arg4 = tan(static_cast<fp_t>(0.5) * tetta_0);
            fp_t root_1 = arg2 / arg4;

            roots[0] = arg2 * arg4;

            if (arg4 != 0 && !isinf(root_1))
                roots[1] = root_1;
            else
                cnt_roots = 1;

            //cout << "Roots: " << roots[0] << "; " << roots[1] << endl;
            cnt_roots = 2;
        } else if (c < 0 and abs(arg1) <= 1) {    // Проверка на вещественные корни и "С" меньше 0
            tetta_n = atan(arg1);
            tetta_0 = tetta_n;
            while (!(0 < abs(tetta_0) < PId_2)) { // Диапазон (0,PI/2)
                // fma(-PI,j,tetta_n), (tetta_n - j * PI)
                tetta_0 = fma(-PI, j, tetta_n);

                (tetta_0 > PId_2) ? j++ : j--;    // Если выходит справа, иначе справа
            }

            fp_t arg4 = tan(static_cast<fp_t>(0.5) * tetta_0);
            fp_t root_1 = arg2 / arg4;

            roots[0] = arg2 * arg4;

            if (arg4 != 0 && !isinf(root_1))
                roots[1] = -root_1;
            else
                cnt_roots = 1;

            //cout << "Roots: " << roots[0] << "; " << roots[1] << endl;
            cnt_roots = 2;
        } else { // Если комплексные корни
            cnt_roots = 0;
            //cout << "Only complex roots " << endl;
        }
    } else if (isinf(c_a)) {    // c/a = inf - значит уравнение формально не квадратное, находим единственный корень
        roots[0] = -c_b;      // не проверяем b!=0, т.к. делаем проверку корня на бесконечность
        if (!isinf(roots[0])) cnt_roots = 1;
        else cnt_roots = 0;
    } else if (isinf(c_b)) {    // уже знаем что a != 0
        if (-c_a >= 0) {        // Проверка на вещественность
            fp_t b_0 = sqrt(-c_a);
            roots[0] = -b_0;
            roots[1] = b_0;
            cnt_roots = 2;
        } else cnt_roots = 0;
    }
    /*ax^2+bx+c=0 if a is any
    ax^2+c=0 // if b==0
    x^2 = -c/a
    x = +-sqrt(-c/a)

    ax^2+bx+c=0 if a and b == 0
    c=0
    null set
    */
    return cnt_roots;
}


template<typename fp_t>
pair <fp_t, fp_t> testPolynomial(unsigned int roots_count) {
    fp_t max_absolute_error, max_relative_error;
    vector <fp_t> roots_computed(roots_count);
    vector <fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots,
                              coefficients);
    int cnt_real_roots = trigonometric(coefficients, roots_computed);
    if (cnt_real_roots != 0 && cnt_real_roots != -1) {
        compare_roots2<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
                            max_absolute_error, max_relative_error);

    } else max_absolute_error = 0, max_relative_error = 0;
    return pair<fp_t, fp_t>(max_absolute_error, max_relative_error);
}


typedef float fp_t;

int main() {
    float max_absolut_deviation = 0;
    float max_relative_deviation = 0;
    for (auto i = 0; i < 10'000'000; ++i) {
        auto deviation = testPolynomial<float>(2);
        //cout << "ABS_DEVATION = " << deviation.first << endl;

        if (deviation.first > max_absolut_deviation) {
            max_absolut_deviation = deviation.first;
        }
        if (deviation.second > max_relative_deviation) {
            max_relative_deviation = deviation.second;
        }
    }
    cout << endl << "MAX_ABSOLUT_deviation = " << max_absolut_deviation << endl;
    cout << endl << "MAX_RELATIVE_deviation = " << max_relative_deviation << endl;

}
