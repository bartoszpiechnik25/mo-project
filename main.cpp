#include <iostream>
#include <map>
#include <mutex>
#include <future>
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

void calculateMaxErrors(const int t_kmb, const int t_cn, const string& path) {
    double h = 0.1;

    ofstream file(path + "/max_errors.csv");
    file << "log_h,log_error_analytical,log_error_kmb,log_error_cn_thomas,log_error_cn_jacobi" << endl;

    std::mutex file_mutex;

        for(auto i = 0; i < 5; ++i) {
            cout << "Calculating for h = " << h << endl;
            const int m = static_cast<int>((X_MAX - X_MIN) / h) + 1;

            auto analityczna_future = std::async(std::launch::async, macierzRozwiazanieAnalityczne, t_cn, m, h);
            auto analitycznaKMB_future = std::async(std::launch::async, macierzRozwiazanieAnalityczne, t_kmb, m, h);
            auto kmb_future = std::async(std::launch::async, macierzKMB, t_kmb, m);
            auto cn_thomas_future = std::async(std::launch::async, macierzCrankNicolsonThomas, t_cn, m);
            auto cn_jacobi_future = std::async(std::launch::async, macierzCrankNicolsonJacobi, t_cn, m);

            const auto analityczna = analityczna_future.get();
            const auto analitycznaKMB = analitycznaKMB_future.get();
            const auto kmb = kmb_future.get();
            const auto cn_thomas = cn_thomas_future.get();
            const auto cn_jacobi = cn_jacobi_future.get();

            const auto kmb_error = macierzBledu(analitycznaKMB, kmb, t_kmb, m);
            const auto cn_thomas_error = macierzBledu(analityczna, cn_thomas, t_cn, m);
            const auto cn_jacobi_error = macierzBledu(analityczna, cn_jacobi, t_cn, m);

            const auto max_kmb_error = bladMaksymalny(kmb_error, t_kmb, m)[t_kmb - 1];
            const auto max_cn_thomas_error = bladMaksymalny(cn_thomas_error, t_cn, m)[t_cn - 1];
            const auto max_cn_jacobi_error = bladMaksymalny(cn_jacobi_error, t_cn, m)[t_cn - 1];

            std::lock_guard<std::mutex> lock(file_mutex);
            file << log10(h) << ","
                 << log10(max_kmb_error) << ","
                 << log10(max_cn_thomas_error) << ","
                 << log10(max_cn_jacobi_error) << endl;

            usuwanieMacierzy(analityczna, t_cn);
            usuwanieMacierzy(analitycznaKMB, t_kmb);
            usuwanieMacierzy(kmb, t_kmb);
            usuwanieMacierzy(cn_thomas, t_cn);
            usuwanieMacierzy(cn_jacobi, t_cn);
            h /= 2;

    }

    file.close();
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
    calculateMaxErrors(t_kmb, t_cn, path);

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
