#include <iostream>
#include "stale.h"
#include "utils.h"
#include "solver.h"
#include <cassert>

using namespace std;

double** macierzCrankNicolsonThomas(int n, int m) {
    double** thomas = alokujMacierz < double >(n, m);
    warunekPoczatkowy(thomas, m);
    warunkiBrzegowe(thomas, n, m);
    equationSolver(thomas, n, m, algorytmThomasa, dyskretyzacjaCrankaNicolson);
    return thomas;
}

int main() {

    const int t_kmb = (int)((T_MAX - T_MIN) / DT_KMB) + 1 + 1;
    const int t_cn = (int)((T_MAX - T_MIN) / DT_CN) + 1 + 1;
    const int m = (int)((X_MAX - X_MIN) / H) + 1;
    const string path = "../wyniki";

    cout << "Wymiar czasu KMB: " << t_kmb << endl
        << "Wymoar czasu CN: " << t_cn << endl
        << "Wymiar przestrzenny m: " << m << endl
        << "Krok czasowy dla KMB: " << DT_KMB << endl
        << "Krok czasowy dla CN: " << DT_CN << endl
        << "Krok przestrzenny: " << H << endl;

    auto analityczna = macierzRozwiazanieAnalityczne(t_cn, m, H, DT_CN);
    auto krokiT = krokCzasowy(DT_CN, t_cn);
    auto krokiX = krokPrzestrzenny(m);

    auto p1 = path + "/analityczna_1.csv";
    auto p2 = path + "/analityczna_3_kroki.csv";
    auto p3 = path + "/x_wartosxi.csv";
    auto p4 = path + "/var_wartosxii.csv";

    printArray(krokiX, m);
    cout << endl << endl;
    printArray(krokiT, t_cn);

    zapiszMacierzDoPliku <double>(analityczna, t_cn, m, p1);
    zapisWynikowDoPliku(analityczna, krokiX, t_cn, m, p2, DT_CN);
    zapiszWektorDoPliku <double>(krokiX, m, p3);
    zapiszWektorDoPliku <double>(krokiT, t_cn, p4);
    return 0;
}
