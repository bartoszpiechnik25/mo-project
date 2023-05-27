//
// Created by barto on 5/23/2023.
//

#ifndef PROJEKTLAB11_STALE_H
#define PROJEKTLAB11_STALE_H
#include <cmath>

constexpr short T_MAX = 2;
constexpr short T_MIN = 0;
constexpr short D = 1;
constexpr double tau = 0.1;
constexpr double EPS = 1e-20;
constexpr short ITER = 200;
constexpr double H = 0.1;
const double a = 6 * sqrt(D * (tau + T_MAX));
const double X_MAX = a;
const double X_MIN = -a;
constexpr double LAMBDA_KMB = 0.4;
constexpr double LAMBDA_CN = 1.;

//delta_t (dt) dla dyskretyzacji bezpośredniej
const double DT_KMB = (0.4 * H * H) / D;
//delta_t (dt) dla dyskretyzacji pośredniej Cranka-Nicolson
const double DT_CN = (1.0 * H * H) / D;

void warunekPoczatkowy(double** matrix, const int n, const int m) {
    for (int i = 0; i < n; ++i) {
        matrix[i][0] = .0;
        matrix[i][m - 1] = .0;
    }
}

void warunekBrzegowy(double** matrix, const int m) {
    double x = X_MIN;
    for (int i = 0; i < m; ++i) {
        matrix[0][i] = ((1.0 / (2.0 * sqrt(M_PI * D * tau))) * exp(-1.0 * (x * x) / (4 * D * tau)));
        x += h;
    }
}

#endif //PROJEKTLAB11_STALE_H
