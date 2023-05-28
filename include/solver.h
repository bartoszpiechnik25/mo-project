//
// Created by barto on 5/23/2023.
//

#ifndef PROJEKTLAB11_SOLVER_H
#define PROJEKTLAB11_SOLVER_H
#include "stale.h"
#include "utils.h"
#include <functional>

struct wektory {
    double* l, * d, * u, * b, * x;
    const int size;

    explicit wektory(const int size): size(size) {
        l = alokacjaWektora<double>(size);
        d = alokacjaWektora<double>(size);
        u = alokacjaWektora<double>(size);
        b = alokacjaWektora<double>(size);
        x = alokacjaWektora<double>(size);
    }
    ~wektory() {
        usuwanieWektora<double>(l);
        usuwanieWektora<double>(d);
        usuwanieWektora<double>(u);
        usuwanieWektora<double>(b);
        usuwanieWektora<double>(x);
    }
};

constexpr double rozwiazanieAnalityczne(const double x, const double t) {
    return (1.0 / (2.0 * sqrt(M_PI * D * (t + tau)))) * exp(-1.0 * (x * x) / (4 * D * (t + tau)));
}

double** macierzRozwiazanieAnalityczne(const int n, const int m, const double h, const double dt) {
    auto matrix = alokujMacierz<double>(n, m);
    double x = X_MIN;
    double t = T_MIN;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            matrix[i][j] = rozwiazanieAnalityczne(x, t);
            x += h;
        }
        x = X_MIN;
        t += dt;
    }
    return matrix;
}

void algorytmThomasa(wektory& w, const int m) {
    int i = 2;
    for (; i < m; ++i) {
        //obliczanie wyrazu na przekątnej
        w.d[i] = w.d[i] - (w.l[i - 1] / w.d[i - 1]) * w.u[i - 1];
        //obliczanie wyrazu wolnego
        w.b[i] = w.b[i] - (w.l[i - 1] / w.d[i - 1]) * w.b[i - 1];
    }
    w.x[m - 1] = w.b[m - 1] / w.d[m - 1];

    i = m - 2;
    //obliczenie wektora rozwiązań
    for (; i >= 0; ++i)
        w.x[i] = (w.b[i] - w.u[i] * w.x[i + 1]) / w.d[i];
}

void algorytmJacobiego(wektory& w, int m) {
    // zamienic ldu na A (alokacja i wypelnianie zerami)
    double** A = new double* [m];
    for (int i = 0; i < m; ++i) {
        A[i] = new double[m];
        for (int j = 0; j < m; ++j) A[i][j] = 0.0;
    }

    // wypelnianie macierzy trzema wektorami
    A[0][0] = w.d[0];
    A[0][1] = w.u[0];
    for (int i = 1; i < m - 1; ++i) {
        A[i][i - 1] = w.l[i];
        A[i][i] = w.d[i];
        A[i][i + 1] = w.u[i];
    }
    A[m - 1][m - 1] = w.d[m - 1];
    A[m - 1][m - 2] = w.l[m - 1];
    // wektor kolejnego przyblizenia
    double* nastepnyX = alokacjaWektora < double >(m);
    for (int step = 0; step < ITER; step++) {
        for (int i = 0; i < m; i++) {
            // obliczenie wartosci funkcji dla i wiersza (dla L+U)
            nastepnyX[i] = (w.b[i] - policzWartosc(A, w.x, i, m)) / A[i][i];
        }
        // warunki koncowe
        double wartoscEstymatora = estymator(w.x, nastepnyX, m);
        double wartoscReziduum = residuum(A, nastepnyX, w.b, m);
        if (wartoscReziduum < EPS || wartoscEstymatora < EPS)
            break;
        // x następne do x
        vector_clone(nastepnyX, w.x, m);
    }
    // wektor x wynikowy
    vector_clone(nastepnyX, w.x, m);
}

void dyskretyzacjaKMB(double **A, const int n, const int m) {
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m - 1; ++j)
            A[i][j] = A[i - 1][j] + LAMBDA_KMB * (A[i - 1][j - 1] - (2 * A[i - 1][j]) + A[i - 1][j + 1]);
    }
}

void dyskretyzacjaCrankaNicolson(wektory& w, double** A, const int k, const int m) {
    for (int i = 1; i < m - 1; ++i) {
        //lower
        w.l[i] = LAMBDA_CN / 2.;
        //diagonal
        w.d[i] = -(1. + LAMBDA_CN);
        //upper
        w.u[i] = LAMBDA_CN / 2.;
        //wektor b
        w.b[i] = -(LAMBDA_CN / 2.0 * A[k - 1][i - 1] + (1.0 - LAMBDA_CN) * A[k - 1][i] + (LAMBDA_CN / 2.0) * A[k - 1][i + 1]);
    }
}

void equationSolver(double** A,
                      const int n, const int m,
                      std::function<void(wektory&, int)> solver,
                      std::function<void(wektory&, double**, int, int)> discretization) {
    wektory w(m);
    for (int k = 1; k < n; ++k) {
        w.l[0] = 0.0;
        w.d[0] = 1.0;
        w.u[0] = 0.0;
        w.b[0] = A[k - 1][0];
        discretization(w, A, k, m);
        w.l[m - 1] = 0.0;
        w.d[m - 1] = 1.0;
        w.u[m - 1] = 0.0;
        w.b[m - 1] = 0.0;
        solver(w, m);
        for (int i = 1; i < m - 1; ++i) {
            A[k][i] = w.x[i];
        }
    }
}

#endif //PROJEKTLAB11_SOLVER_H
