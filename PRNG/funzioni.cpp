#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "funzioni.h"
#include "random.h"

using namespace std;

// Metodo per integrare una funzione tramite il metodo dei blocchi
void FunzioneBase::integraBlocchi(int M, int N, Random& rnd, const string& file_name) const {

    // Controllo che il numero di punti e blocchi sia positivo
    if (M <= 0 || N <= 0) {
        cerr << "Errore: M e N devono essere positivi." << endl;
        return;
    }

    // Controllo che M sia divisibile per N (ogni blocco deve avere lo stesso numero di punti)
    if (M % N != 0) {
        cerr << "Errore: M deve essere divisibile per N." << endl;
        return;
    }

    int L = M / N; // Numero di punti per blocco
    double meanTotal = 0.;            // Media cumulativa dei blocchi
    double stdDevTotal = 0.;          // Deviazione standard cumulativa
    double sumBlock;                  // Somma dei valori di un blocco
    double sumMeanBlocks = 0.;        // Somma delle medie dei blocchi
    double meanBlock = 0.;            // Media del singolo blocco
    double sumSquaredMeanBlocks = 0.; // Somma dei quadrati delle medie dei blocchi
    double invL = 1.0 / L;            // Inverso di L, per calcolare la media velocemente

    // Apertura del file di output
    ofstream outputFile(file_name);
    if (!outputFile.is_open()) {
        cerr << "Errore: impossibile aprire il file " << file_name << endl;
        return;
    }

    // Ciclo sui blocchi
    for (int i = 0; i < N; i++) {
        sumBlock = 0.;

        // Ciclo sui punti all'interno del blocco
        for (int j = 0; j < L; j++) {
            sumBlock += eval(rnd.Rannyu()); // Eval valuta la funzione in un punto casuale uniforme [0,1)
        }

        meanBlock = sumBlock * invL;          // Media del blocco
        sumMeanBlocks += meanBlock;           // Aggiornamento somma medie blocchi
        sumSquaredMeanBlocks += meanBlock * meanBlock; // Aggiornamento somma quadrati medie

        // Calcolo della media cumulativa e della deviazione standard
        meanTotal = sumMeanBlocks / (i + 1);
        stdDevTotal = (i == 0) ? 0. : sqrt((sumSquaredMeanBlocks / (i + 1) - meanTotal * meanTotal) / i);

        // Scrittura su file: numero di punti, media cumulativa, deviazione standard cumulativa
        outputFile << L * (i + 1) << " " << meanTotal << " " << stdDevTotal << endl;
    }

    outputFile.close();
}

// Metodo per integrare usando il metodo dei blocchi con una funzione di distribuzione cumulativa (CDF)
void FunzioneBase::integraBlocchiCDF(int M, int N, Random& rnd, FunzioneBase& P_r, const string& file_name) const {

    // Controllo parametri
    if (M <= 0 || N <= 0) {
        cerr << "Errore: M e N devono essere positivi." << endl;
        return;
    }

    if (M % N != 0) {
        cerr << "Errore: M deve essere divisibile per N." << endl;
        return;
    }

    int L = M / N;
    double meanTotal = 0.;
    double stdDevTotal = 0.;
    double sumBlock;
    double sumMeanBlocks = 0.;
    double meanBlock = 0.;
    double sumSquaredMeanBlocks = 0.;
    double invL = 1.0 / L;
    double x;

    ofstream outputFile(file_name);
    if (!outputFile.is_open()) {
        cerr << "Errore: impossibile aprire il file " << file_name << endl;
        return;
    }

    // Ciclo sui blocchi
    for (int i = 0; i < N; i++) {
        sumBlock = 0.;

        // Ciclo sui punti del blocco
        for (int j = 0; j < L; j++) {
            x = P_r.eval(rnd.Rannyu()); // Trasformazione tramite la CDF
            sumBlock += eval(x);         // Valutazione della funzione originale sul valore trasformato
        }

        meanBlock = sumBlock * invL;
        sumMeanBlocks += meanBlock;
        sumSquaredMeanBlocks += meanBlock * meanBlock;

        meanTotal = sumMeanBlocks / (i + 1);
        stdDevTotal = (i == 0) ? 0. : sqrt((sumSquaredMeanBlocks / (i + 1) - meanTotal * meanTotal) / i);

        outputFile << L * (i + 1) << " " << meanTotal << " " << stdDevTotal << endl;
    }

    outputFile.close();
}

// Implementazioni dei metodi eval delle varie funzioni
double Linear::eval(double x) const {
    return x; // Funzione lineare f(x) = x
}

double Cos_g::eval(double x) const {
    return (M_PI / 4.) * (cos((M_PI / 2.) * x) / (1. - x)); // Funzione coseno modificata
}

double Cos::eval(double x) const {
    return (M_PI / 2.) * (cos((M_PI / 2.) * x)); // Funzione coseno semplice
}

double var::eval(double x) const {
    return (x - 0.5) * (x - 0.5); // Varianza rispetto a 0.5
}

double P_r_2::eval(double x) const {
    return 1 - sqrt(1. - x); // Trasformazione inversa per la CDF
}

double P_r_exp::eval(double x) const {
    return -log(1 - x) / lambda; // Trasformazione inversa per distribuzione esponenziale
}

double P_r_Lorentz::eval(double x) const {
    return x0 + gamma * tan(M_PI * (x - 0.5)); // Trasformazione inversa per distribuzione di Lorentz
}

