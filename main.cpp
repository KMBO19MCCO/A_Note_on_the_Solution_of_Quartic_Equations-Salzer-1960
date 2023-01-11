#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>
#include "except.h"

#define PRINT true

typedef float fp_t;

using namespace std;

// ================================= ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ================================= //

template<typename fp_t>
inline bool isZero(const fp_t& x)
{
    return FP_ZERO == fpclassify(x);
}

template <typename fp_t>
inline fp_t fpFix(fp_t x)
{
    return x = abs(x) < numeric_limits<fp_t>::epsilon() ? static_cast<fp_t>(0.0L) : x;
}

template<typename fp_t>
inline fp_t fms(fp_t a, fp_t b, fp_t c, fp_t d) // a * b - c * c
{
    fp_t cd = -c * d;

    return fma(a, b, cd) - fma(c, d, cd);
}

template <typename fp_t>
inline int sgn(fp_t x)
{
    return (fp_t(0) < x) - (x < fp_t(0));
}

// ======================== ФУНКЦИИ ПОИСКА КОРНЕЙ СТЕПЕННЫХ ПОЛИНОМОВ ======================== //

template<typename fp_t>
unsigned int solveLinear(fp_t n, fp_t a, vector<fp_t>& roots) { return 0; }

template<typename fp_t>
int solveQuadratic(fp_t n, fp_t a, fp_t b, vector<fp_t>& roots) { return 0; }

/*
    Имплементация метода решения кубического уравнения — Cubic_Formula-Wolfram
    Информация о методе — https://mathworld.wolfram.com/CubicFormula.html
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveCubic(fp_t n, fp_t a, fp_t b, fp_t c, vector<fp_t>& roots)
{
    // Нормировка исходных коэффициентов, где n - старший коэффициент уравнения
    if (isZero(n) || isinf(a /= n))
        return solveQuadratic(a, b, c, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;

    // Объявление констант
    const fp_t PI = static_cast<fp_t>(M_PI);
    const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);

    // Расчетные коэффициенты
    fp_t Q = fpFix(fms(static_cast<fp_t>(3), b, a, a)) / static_cast<fp_t>(9);
    fp_t R = fpFix(fms(fms(static_cast<fp_t>(9), b, static_cast<fp_t>(2) * a, a), a, static_cast<fp_t>(27), c)) / static_cast<fp_t>(54);
    fp_t C = R * pow(abs(Q), static_cast<fp_t>(-1.5L));

    // Дискриминант
    fp_t D = fpFix(fms(Q, Q * Q, -R, R));

    // Количество вещественных корней
    unsigned int numberOfRoots = D <= 0 ? 3 : 1;

    if (D < 0) // D < 0: Все три корня вещественные и различные
    {
        fp_t fi = acos(R / pow(-Q, static_cast<fp_t>(1.5L)));

        fp_t sqrtQ = static_cast<fp_t>(2) * sqrt(-Q);

        roots =
        {
            fms(sqrtQ, cos(fi * ONE_THIRD), a, ONE_THIRD),
            fms(sqrtQ, cos(fma(static_cast<fp_t>(2), PI, fi) * ONE_THIRD), a, ONE_THIRD),
            fms(sqrtQ, cos(fma(static_cast<fp_t>(4), PI, fi) * ONE_THIRD), a, ONE_THIRD)
        };
    }
    else if (D > 0) // D > 0: Только один корень вещественный, остальные комплексные
    {
        fp_t y = Q > 0 ? static_cast<fp_t>(2) * sinh(asinh(C) * ONE_THIRD) :
                         Q < 0 && C >= 1 ? static_cast<fp_t>(2) * cosh((acosh(C) * ONE_THIRD)) :
                         Q < 0 && C <= -1 ? static_cast<fp_t>(2) * -cosh((acosh(abs(C)) * ONE_THIRD)) :
                         static_cast<fp_t>(2) * cos(acos(C) * ONE_THIRD);

        roots =
        {
            fms(sqrt(abs(Q)), y, a, ONE_THIRD)
        };
    }
    else // D = 0: Все три корня вещественные, один из которых кратен двум
    {
        roots =
        {
            fms(sgn(-R) * static_cast<fp_t>(-2), sqrt(abs(Q)), a, ONE_THIRD),
            fms(-a, ONE_THIRD, static_cast<fp_t>(sgn(R)), sqrt(abs(Q))),
            fms(-a, ONE_THIRD, static_cast<fp_t>(sgn(R)), sqrt(abs(Q)))
        };
    }

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — A_Note_on_the_Solution_of_Quartic_Equations-Salzer-1960
    Информация о методе — https://www.ams.org/journals/mcom/1960-14-071/S0025-5718-1960-0117882-6/S0025-5718-1960-0117882-6.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveQuartic(fp_t N, fp_t A, fp_t B, fp_t C, fp_t D, vector<fp_t>& roots)
{
    // Нормировка исходных коэффициентов, где N - старший коэффициент уравнения
    if (isZero(N) || isinf(A /= N))
        return solveCubic(A, B, C, D, roots);
    if (isinf(B /= N))
        return 0;
    if (isinf(C /= N))
        return 0;
    if (isinf(D /= N))
        return 0;

    // Количество вещественных корней
    unsigned int numberOfRoots = 0;

    // Вычисление коэффициентов кубического уравнения
    fp_t a = -B;
    fp_t b = fms(A, C, static_cast<fp_t>(4), D);
    fp_t c = fms(D, fms(static_cast<fp_t>(4), B, A, A), C, C);
    vector<fp_t> cubicRoots(3);

    // Решение кубического уравнения
    solveCubic(static_cast<fp_t>(1), a, b, c, cubicRoots);

    // Вещественный корень кубического уравнения
    fp_t x = cubicRoots[0];

    // Вычисление начальных расчетных коэффициентов
    fp_t m;
    fp_t n;
    fp_t mm = fpFix(fma(static_cast<fp_t>(0.25L) * A, A, -B) + x);

    if (mm > 0)
    {
        m = sqrt(mm);
        n = static_cast<fp_t>(0.25L) * fpFix(fms(A, x, static_cast<fp_t>(2), C)) / m;
    }
    else if (isZero(mm))
    {
        m = 0;
        n = sqrt(fpFix(fma(static_cast<fp_t>(0.25L) * x, x, -D)));
    }
    else // m - комплексное, следовательно уравнение не будет иметь вещественных корней
    {
        return 0;
    }

    // Вычисление расчетных коэффициентов
    fp_t alpha = fma(static_cast<fp_t>(0.5L) * A, A, -x) - B;
    fp_t beta = fms(static_cast<fp_t>(4), n, A, m);
    fp_t gamma = fpFix(alpha + beta);
    fp_t delta = fpFix(alpha - beta);

    if (gamma >= 0) // Если gamma >= 0, то уравнение имеет два либо больше вещественных корней
    {
        gamma = sqrt(gamma);

        roots[0] = fma(static_cast<fp_t>(0.5L), gamma, fms(static_cast<fp_t>(0.5L), m, static_cast<fp_t>(0.25L), A));
        roots[1] = fma(static_cast<fp_t>(-0.5L), gamma, fms(static_cast<fp_t>(0.5L), m, static_cast<fp_t>(0.25L), A));

        numberOfRoots += 2;
    }
    if (delta >= 0) // Если delta >= 0, то уравнение имеет два либо больше вещественных корней
    {
        delta = sqrt(delta);

        roots[numberOfRoots] = fma(static_cast<fp_t>(0.5L), delta, fms(static_cast<fp_t>(-0.5L), m, static_cast<fp_t>(0.25L), A));
        roots[numberOfRoots + 1] = fma(static_cast<fp_t>(-0.5L), delta, fms(static_cast<fp_t>(-0.5L), m, static_cast<fp_t>(0.25L), A));

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

// ======================== ФУНКЦИИ ТЕСТИРОВАНИЯ МЕТОДОВ ======================== //

template <typename fp_t>
void testCubicAdv(const int testCount, const fp_t dist)
{
    int P = 3; // power, total number of tests
    fp_t low = -1, high = 1; // [low, high], max distance between clustered roots
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0;
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;
    int count = 0;

    vector<fp_t> coefficients(P + 1);
    vector<fp_t> trueRoots(P);

    for (size_t i = 0; i < testCount; ++i)
    {
        vector<fp_t> foundRoots(P);

        generate_polynomial<fp_t>(P, 0, P, 0, dist, low, high, trueRoots, coefficients);
        numOfFoundRoots = solveCubic<fp_t>(coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        compare_roots<fp_t>(numOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError);

        if (isinf(absMaxError))
            cantFind += 1;
        else
        {
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;

            count += relMaxError > 1 ? 1 : 0;
        }
    }

    if (PRINT)
    {
        cout << "\n\n\t\t\tCUBIC TEST RESULTS\n\n";
        cout << "Max distance: " << dist << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "Mean absMaxError = " << absErrors / (testCount - cantFind) << endl;
        cout << "Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: " << maxAbsAllofTest << endl;
        cout << "Mean RelMaxError = " << relError / (testCount - cantFind) << endl;
        cout << "Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: " << maxRelAllofTest << endl;
        cout << "RelMaxError > 1: " << count << " times" << endl;
    }
}

template<typename fp_t>
void testQuarticAdv(const int testCount, const fp_t dist) {
    int P = 4; // power
    fp_t low = -1, high = 1; // [low, high]
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0;
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;
    int count = 0;

    vector<fp_t> coefficients(P + 1);

    for (size_t i = 0; i < testCount; ++i) {
        vector<fp_t> foundRoots(P);
        vector<fp_t> trueRoots(P);

        generate_polynomial<fp_t>(P, 0, P, 0, dist,
            low, high, trueRoots, coefficients);
        numOfFoundRoots = solveQuartic<fp_t>(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        compare_roots<fp_t>(numOfFoundRoots, 4, foundRoots, trueRoots, absMaxError, relMaxError);

        if (isinf(absMaxError))
            cantFind += 1;
        else {
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;
            count += relMaxError > 1 ? 1 : 0;
        }
    }

    if (PRINT)
    {
        cout << setprecision(10);
        cout << "\n\n\t\t\tQUARTIC TEST RESULTS\n\n";
        cout << "Max distance: " << dist << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "Mean absMaxError = " << absErrors / (testCount - cantFind) << endl;
        cout << "Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: " << maxAbsAllofTest << endl;
        cout << "Mean RelMaxError = " << relError / (testCount - cantFind) << endl;
        cout << "Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: " << maxRelAllofTest << endl;
        cout << "RelMaxError > 1: " << count << " times" << endl;
    }
}

int main()
{
    testCubicAdv(1000000, static_cast<fp_t>(1e-5));
    testQuarticAdv(1000000, static_cast<fp_t>(1e-5));

    return 0;
}
