#include <iostream>
#include <map>
#include "stale.h"
#include "utils.h"
#include "solver.h"

using namespace std;


auto macierzCrankNicolsonThomas(const int n, const int m) {
    double** thomas = alokujMacierz <double>(n, m);
    warunekPoczatkowy(thomas, m);
    warunkiBrzegowe(thomas, n, m);
    equationSolver(thomas, n, m, algorytmThomasa, dyskretyzacjaCrankaNicolson);
    return thomas;
}

auto macierzCrankNicolsonJacobi(const int n, const int m) {
    auto jacobi = alokujMacierz <double>(n, m);
    warunekPoczatkowy(jacobi, m);
    warunkiBrzegowe(jacobi, n, m);
    equationSolver(jacobi, n, m, algorytmJacobiego, dyskretyzacjaCrankaNicolson);
    return jacobi;
}

auto macierzKMB(const int n, const int m) {
    auto kmb = alokujMacierz <double>(n, m);
    warunekPoczatkowy(kmb, m);
    warunkiBrzegowe(kmb, n, m);
    dyskretyzacjaKMB(kmb, n, m);
    return kmb;
}

int main() {

    const int t_kmb = static_cast<int>((T_MAX - T_MIN) / DT_KMB) + 1 + 1;
    const int t_cn = static_cast<int>((T_MAX - T_MIN) / DT_CN) + 1 + 1;
    const int m = static_cast<int>((X_MAX - X_MIN) / H) + 1;
    const string path = "../wyniki";

    cout << "Wymiar czasu KMB: " << t_kmb << endl
         << "Wymoar czasu CN: " << t_cn << endl
         << "Wymiar przestrzenny m: " << m << endl
         << "Krok czasowy dla KMB: " << DT_KMB << endl
         << "Krok czasowy dla Cranka-Nicolson: " << DT_CN << endl
         << "Krok przestrzenny: " << H << endl;

    const auto krokiX = krokPrzestrzenny(m);
    const auto krokiT = krokCzasowy(DT_CN, t_cn);
    const auto krokiT_KMB = krokCzasowy(DT_KMB, t_kmb);
    const auto analitycznaKMB = macierzRozwiazanieAnalityczne(t_kmb, m, DT_KMB);
    zapisWynikowDoPliku(analitycznaKMB, krokiX, t_kmb, m, path + "/analitycznaKMB.csv");

    map<string, double**> functions = {
            {"analityczna", macierzRozwiazanieAnalityczne(t_cn, m, DT_CN)},
            {"KMB", macierzKMB(t_kmb, m)},
            {"Crank-NicolsonThomas", macierzCrankNicolsonThomas(t_cn, m)},
            {"Crank-NicolsonJacobi", macierzCrankNicolsonJacobi(t_cn, m)}
    };

    for(auto& pair: functions) {
        auto filename = path + "/" + pair.first;
        auto t = pair.first == "KMB" ? t_kmb : t_cn;
        auto kroki = pair.first == "KMB" ? krokiT_KMB : krokiT;
        auto analityczne = pair.first == "KMB" ? analitycznaKMB : functions["analityczna"];
        zapisWynikowDoPliku(pair.second, krokiX, t, m, filename);
        auto bledy = macierzBledu(analityczne, pair.second, t, m);
        auto maxBlad = bladMaksymalny(bledy, t, m);
        zapiszDwaWektoryDoPliku<double>(maxBlad, kroki, t,path + "/max_bledy_" + pair.first + ".csv");
    }
    zapiszWektorDoPliku(krokiT, t_cn, path + "/krokiT.csv");
    zapiszWektorDoPliku(krokiT_KMB, t_kmb, path + "/krokiT_KMB.csv");
    zapiszWektorDoPliku(krokiX, m, path + "/krokiX.csv");

    return 0;
}
