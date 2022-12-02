#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

typedef float fp_t;

using namespace std;

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

/*
    ������������� ������ ������� ����������� ��������� � Cubic_Formula-Wolfram
    ���������� � ������ � https://mathworld.wolfram.com/CubicFormula.html
    ������ �������� � ������� ������� (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveCubic(fp_t n, fp_t a, fp_t b, fp_t c, vector<fp_t>& roots)
{
    // ���������� �������� �������������, ��� n - ������� ����������� ���������
    a /= n;
    b /= n;
    c /= n;

    // ���������� ��������
    const fp_t PI = static_cast<fp_t>(M_PI);
    const fp_t ONE_THIRD = static_cast<fp_t>(1.L / 3.L);

    // ��������� ������������
    fp_t Q = fms(static_cast<fp_t>(3), b, a, a) / static_cast<fp_t>(9);
    fp_t R = fms(fms(static_cast<fp_t>(9), b, static_cast<fp_t>(2) * a, a), a, static_cast<fp_t>(27), c) / static_cast<fp_t>(54);
    fp_t C = R * pow(abs(Q), static_cast<fp_t>(-1.5));

    // ������������
    fp_t D = fms(Q, Q * Q, -R, R);

    // ���������� ������������ ������
    unsigned int numberOfRoots = D <= 0 ? 3 : 1;

    if (D < 0) // D < 0: ��� ��� ����� ������������ � ���������
    {
        fp_t fi = acos(R / pow(-Q, static_cast<fp_t>(1.5)));

        fp_t sqrtQ = static_cast<fp_t>(2) * sqrt(-Q);

        roots =
        {
            fms(sqrtQ, cos(fi * ONE_THIRD), a, ONE_THIRD),
            fms(sqrtQ, cos(fma(static_cast<fp_t>(2), PI, fi) * ONE_THIRD), a, ONE_THIRD),
            fms(sqrtQ, cos(fma(static_cast<fp_t>(4), PI, fi) * ONE_THIRD), a, ONE_THIRD)
        };
    }
    else if (D > 0) // D > 0: ������ ���� ������ ������������, ��������� �����������
    {
        fp_t y = Q > 0            ? static_cast<fp_t>(2) * sinh(asinh(C) * ONE_THIRD) :
                 Q < 0 && C >= 1  ? static_cast<fp_t>(2) * cosh((acosh(C) * ONE_THIRD)) :
                 Q < 0 && C <= -1 ? static_cast<fp_t>(2) * -cosh((acosh(abs(C)) * ONE_THIRD)) :
                                    static_cast<fp_t>(2) * cos(acos(C) * ONE_THIRD);

        roots =
        {
            fms(sqrt(abs(Q)), y, a, ONE_THIRD)
        };
    }
    else // D = 0: ��� ��� ����� ������������, ���� �� ������� ������ ����
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
    ������������� ������ ������� ��������� ��������� ������� � A_Note_on_the_Solution_of_Quartic_Equations-Salzer-1960
    ���������� � ������ � https://www.ams.org/journals/mcom/1960-14-071/S0025-5718-1960-0117882-6/S0025-5718-1960-0117882-6.pdf
    ������ �������� � ������� ������� (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solve(fp_t N, fp_t A, fp_t B, fp_t C, fp_t D, vector<fp_t>& roots)
{
    // ���������� �������� �������������, ��� n - ������� ����������� ���������
    A /= N;
    B /= N;
    C /= N;
    D /= N;

    // ���������� ������������ ������
    unsigned int numberOfRoots = 0;

    // ���������� ������������� ����������� ���������
    fp_t a = -B;
    fp_t b = fms(A, C, static_cast<fp_t>(4), D);
    fp_t c = fms(D, fms(static_cast<fp_t>(4), B, A, A), C, C);
    vector<fp_t> cubicRoots(3);

    // ������� ����������� ���������
    solveCubic(static_cast<fp_t>(1), a, b, c, cubicRoots);

    // ������������ ������ ����������� ���������
    fp_t x = cubicRoots[0];

    // ���������� ��������� ��������� �������������
    fp_t m;
    fp_t n;
    fp_t mm = fma(static_cast<fp_t>(0.25) * A, A, -B) + x;

    if (mm > 0)
    {
        m = sqrt(mm);
        n = static_cast<fp_t>(0.25) * fms(A, x, static_cast<fp_t>(2), C) / m;
    }
    else if (mm == 0)
    {
        m = 0;
        n = sqrt(fma(static_cast<fp_t>(0.25) * x, x, -D));
    }
    else // m - �����������, ������������� ��������� �� ����� ����� ������������ ������
    {
        return 0;
    }

    // ���������� ��������� �������������
    fp_t alfa = fma(static_cast<fp_t>(0.5) * A, A, -x) - B;
    fp_t beta = fms(static_cast<fp_t>(4), n, A, m);
    fp_t gamma = alfa + beta;
    fp_t delta = alfa - beta;

    if (gamma >= 0) // ���� gamma >= 0, �� ��������� ����� ��� ���� ������ ������������ ������
    {
        gamma = sqrt(gamma);

        roots[0] = fma(static_cast<fp_t>(0.5), gamma, fms(static_cast<fp_t>(0.5), m, static_cast<fp_t>(0.25), A));
        roots[1] = fma(static_cast<fp_t>(-0.5), gamma, fms(static_cast<fp_t>(0.5), m, static_cast<fp_t>(0.25), A));

        numberOfRoots += 2;
    }
    if (delta >= 0) // ���� delta >= 0, �� ��������� ����� ��� ���� ������ ������������ ������
    {
        delta = sqrt(delta);

        roots[numberOfRoots] = fma(static_cast<fp_t>(0.5), delta, fms(static_cast<fp_t>(-0.5), m, static_cast<fp_t>(0.25), A));
        roots[numberOfRoots + 1] = fma(static_cast<fp_t>(-0.5), delta, fms(static_cast<fp_t>(-0.5), m, static_cast<fp_t>(0.25), A));

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

int main()
{
    /*
    vector<fp_t> r(4);

    fp_t a = 1;
    fp_t b = 1;
    fp_t c = 1;
    fp_t d = 1;
    fp_t e = 1;

    solve(a, b, c, d, e, r);
    */

    return 0;
}