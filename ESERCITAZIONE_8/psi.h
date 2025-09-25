#ifndef __Functions__
#define __Functions__

//using namespace std;

class Psi {
private:
    double mu;     // Parametro mu della funzione d'onda
    double sigma;  // Parametro sigma della funzione d'onda

    // Potenziale esterno
    double V(double x) {
        return x * x * x * x - 2.5 * x * x;
    };

public:
    // Costruttore
    Psi(double mu_, double sigma_);

    //Prendere e Restituire i Parametri
    void setParameters(double _mu, double _sigma){mu = _mu; sigma = _sigma;};
    double getMu() {return mu;};
    double getSigma() {return sigma;};

    // Metodi pubblici
    double eval(double x);          // Valuta la funzione d'onda
    double D2_eval(double x);       // Valuta la derivata seconda
    double loc_energy(double x);    // Calcola l'energia locale
    double abs_val(double x);       // Calcola |ψ(x)|²

};

#endif