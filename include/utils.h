//
// Created by barto on 5/23/2023.
//

#ifndef PROJEKTLAB11_UTILS_H
#define PROJEKTLAB11_UTILS_H

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

template < typename T >
T** alokujMacierz(int n, int m) {
    T** matrix = new T * [n];
    for (int i = 0; i < n; ++i) {
        matrix[i] = new T[m];
        for (int j = 0; j < m; j++)
            matrix[i][j] = 0;
    }
    return matrix;
}

template < typename T >
T* alokacjaWektora(int n) {
    T* _vector = new T[n];
    for (int i = 0; i < n; i++) _vector[i] = 0;
    return _vector;
}

template < typename T >
void usuwanieMacierzy(T** matrix, int n, int m) {
    if (matrix) {
        for (int i = n - 1; i >= 0; --i) {
            delete[] matrix[i];
        }
        delete[] matrix;
    }
}

template < typename T >
void usuwanieWektora(T* _vector) {
    delete[] _vector;
}

//template < typename T >
//void zapiszMacierzDoPliku(T** matrix, int n, int m, std::string& filename) {
//    std::ofstream file(filename, std::ios::out);
//
//    if (file.is_open()) {
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < m; ++j) {
//                file << matrix[i][j] << ";";
//            }
//            file << "\n";
//        }
//    } else {
//        throw std::runtime_error("Nie można otworzyć pliku do zapisu");
//    }
//    file.close();
//}

template < typename T >
void zapiszWektorDoPliku(T* _vector, int n, std::string& filename) {
    std::fstream file(filename.c_str(), std::ios::out);
    if (file.is_open()) {
        for (int i = 0; i < n; ++i) {
            file << _vector[i] << std::endl;
        }
    }
    file.close();
}

template < typename T >
void zapiszDwaWektoryDoPliku(T* _vector, T* _second_vector, int n, std::string& filename) {
    std::ofstream file(filename.c_str(), std::ios::out);
    if (file.is_open()) {
        for (int i = 0; i < n; ++i) {
            file << _second_vector[i] << ", " << _vector[i] << std::endl;
        }
    }
    file.close();
}

template < typename T >
T* vector_subtract(T* vector_1, T* vector_2, const int w) {
    T* vec_buf = new T[w];
    for (int i = 0; i < w; i++) {
        vec_buf[i] = vector_1[i] - vector_2[i];
    }
    return vec_buf;
}

template < typename T >
T vector_max(T* _vector, int n) {
    T max = std::fabs(_vector[0]);
    for (int i = 1; i < n; ++i) {
        if (max < std::fabs(_vector[i])) {
            max = _vector[i];
        }
    }
    return max;
}

template < typename T >
T* vector_residuum(T** matrix, T* x, T* b, const int w) {
    T sum = T();
    T* vector_res = new double[w];
    for (int i = 0; i < w; i++, sum = 0) {
        for (int j = 0; j < w; j++) {
            sum += matrix[i][j] * x[j];
        }
        vector_res[i] = sum - b[i];
    }
    return vector_res;
}

template < typename T >
T vector_norm(T* vector, const int w) {
    T norm = std::fabs(vector[0]);
    for (int i = 0; i < w - 1; i++) {
        if (std::fabs(vector[i]) < std::fabs(vector[i + 1])) {
            norm = std::fabs(vector[i + 1]);
        }
    }
    return norm;
}

template < typename T >
T estymator(T* vector_1, T* vector_2, const int w) {
    T* vec_buf;
    vec_buf = vector_subtract(vector_1, vector_2, w);
    T error = vector_norm(vec_buf, w);
    return error;
}

template < typename T >
T residuum(T** matrix, T* x, T* b, const int w) {
    T* vec_residuum = vector_residuum(matrix, x, b, w);
    return vector_norm(vec_residuum, w);
}

template < typename T >
T policzWartosc(T** matrix, T* x, int i, const int w) {
// obliczenie wartosci funkcji dla i wiersza (dla L+U)
    T scalar = 0;
    for (int j = 0; j < w; ++j) {
        if (i != j) {
            scalar += matrix[i][j] * x[j];
        }
    }
    return scalar;
}

template < typename T >
void vector_clone(T* vector_1, T* vector_2, const int w) {
    for (int i = 0; i < w; i++)
        vector_2[i] = vector_1[i];
}
void zapisWynikowDoPliku(double** macierz_wynikow,
                         double const* x,
                         const int n,
                         const int m,
                         std::string& filename,
                         const double dt) {
    std::fstream file;
    for (int i = 0; i < n; ++i) {
        //krok czasowy
        if (i == 1 || i == 100 || i == 199) {
            file = std::fstream(filename + std::to_string(i) + ".csv",std::ios::out);
            if (file.is_open()) {
                for (int j = 0; j < m; ++j) {
                    // krok przestrzenny
                    file << x[j] << ", " << macierz_wynikow[i][j] << std::endl;
                }
                file << "\n";
            }
        }
    }
    file.close();
}

double* krokPrzestrzenny(const int m) {
    auto kroki = alokacjaWektora<double>(m);
    double x = X_MIN;
    for (int i = 0; i < m; ++i) {
        kroki[i] = x;
        x += H;
    }
    return kroki;
}
double* krokCzasowy(const double dt, const int n) {
    auto kroki = alokacjaWektora<double>(n);
    double t = T_MIN;
    for (int i = 0; i < n; ++i) {
        kroki[i] = t;
        t += dt;
    }
    return kroki;
}

//fuction to print array
template < typename T >
void printArray(T* array, int n) {
    for (int i = 0; i < n; ++i) {
        std::cout << array[i] << ", ";
    }
    std::cout << std::endl;
}

//function to save matrix to file
template < typename T >
void zapiszMacierzDoPliku(T** matrix, int n, int m, std::string& filename) {
    std::ofstream file(filename, std::ios::out);
    if (file.is_open()) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file << matrix[i][j];
                if (j != m)
                    file << ",";
            }
            file << std::endl;
        }
    }
    file.close();
}


#endif //PROJEKTLAB11_UTILS_H
