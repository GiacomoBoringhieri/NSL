#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <string>
#include "psi.h"

using namespace std;

// Costruttore: inizializza mu, sigma e il generatore di numeri casuali
Psi::Psi(double mu_, double sigma_) {
    mu = mu_;
    sigma = sigma_;
};

// Valuta la funzione d'onda come somma di due gaussiane
double Psi::eval(double x) {
    double psi_plus = exp(-((x+mu)*(x+mu)/(2.*sigma*sigma)));
    double psi_min = exp(-((x-mu)*(x-mu)/(2.*sigma*sigma)));
    return psi_min + psi_plus;
};

// Calcola la derivata seconda della funzione d'onda
double Psi::D2_eval(double x) {
    double psi_plus = exp(-((x+mu)*(x+mu)/(2.*sigma*sigma)));
    double psi_min = exp(-((x-mu)*(x-mu)/(2.*sigma*sigma)));
    double factor_plus = -1./(sigma*sigma) + (x+mu)*(x+mu)/(pow(sigma,4.));
    double factor_min = -1./(sigma*sigma) + (x-mu)*(x-mu)/(pow(sigma,4.));
    return factor_plus * psi_plus + factor_min * psi_min;
};

// Calcola l'energia locale: termine cinetico + potenziale
double Psi::loc_energy(double x) {
    double psi_eval = eval(x);
    double D2_psi_eval = D2_eval(x);
    return - 0.5 * D2_psi_eval / psi_eval + V(x);
};

// Calcola il modulo quadro della funzione d'onda
double Psi::abs_val(double x) {
    return eval(x) * eval(x);
};
