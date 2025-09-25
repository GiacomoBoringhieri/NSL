#include <iostream>
#include <fstream>
#include <string>
#include "../../PRNG/random.h"
#include "../../PRNG/funzioni.h"

using namespace std;

// Funzione max per calcolare il payoff (utile per call e put)
double max(double x, double y){
    if(x > y || x == y) return x;
        else return y;
};

int main(int argc, char *argv[]){

    int M=100000;  // Numero totale di simulazioni Monte Carlo
    int N=100;     // Numero di blocchi per il calcolo delle medie

    Random rnd;    // Generatore di numeri casuali
    rnd.Start(1);  // Inizializzazione con seed 1

    // Parametri del contratto e modello (Black-Scholes)
    double S0 = 100.0, T = 1.0, K = 100.0, mu = 0.1, sigma = 0.25;
    double z, S;
    int L = M / N; 
    double invL = 1.0 / L;

    // File di output per call e put (valore finale)
    ofstream outputFile_call("data_call.dat");
    ofstream outputFile_put("data_put.dat");

    // Controllo apertura file
    if (!outputFile_call.is_open() || !outputFile_put.is_open()) {
        cerr << "Errore: impossibile aprire i file" << endl;
        return -1;
    }

    // Variabili per medie progressive e deviazioni standard (call e put)
    double meanTotal_call = 0., stdDevTotal_call = 0., sumBlock_call, sumMeanBlocks_call = 0., meanBlock_call = 0., sumSquaredMeanBlocks_call = 0.;
    double meanTotal_put = 0., stdDevTotal_put = 0., sumBlock_put, sumMeanBlocks_put = 0., meanBlock_put = 0., sumSquaredMeanBlocks_put = 0.;

    // === PRIMO METODO: simulazione diretta con T finale ===
    for (int i = 0; i < N; i++) {
        sumBlock_call = 0.;
        sumBlock_put = 0.;

        for (int j = 0; j < L; j++) {
            z = rnd.Gauss(0., 1.);  // Estrazione normale standard
            // Evoluzione del prezzo secondo modello log-normale
            S = S0*exp((mu - 0.5*sigma*sigma)*T + sigma*z*sqrt(T));
            // Payoff scontato
            sumBlock_call += exp(-mu*T)*max(0., S-K);
            sumBlock_put += exp(-mu*T)*max(0., K-S);
        }

        // Medie e deviazioni progressive (metodo dei blocchi)
        meanBlock_call = sumBlock_call * invL; 
        sumMeanBlocks_call += meanBlock_call;
        sumSquaredMeanBlocks_call += meanBlock_call * meanBlock_call; 

        meanBlock_put = sumBlock_put * invL; 
        sumMeanBlocks_put += meanBlock_put;
        sumSquaredMeanBlocks_put += meanBlock_put * meanBlock_put; 

        meanTotal_call = sumMeanBlocks_call / (i + 1);
        stdDevTotal_call = sqrt((sumSquaredMeanBlocks_call / (i + 1) - (meanTotal_call * meanTotal_call)) / (i + 1));

        meanTotal_put = sumMeanBlocks_put / (i + 1);
        stdDevTotal_put = sqrt((sumSquaredMeanBlocks_put / (i + 1) - (meanTotal_put * meanTotal_put)) / (i + 1));

        outputFile_call << L * (i + 1) << " " << meanTotal_call << " " << stdDevTotal_call << endl;
        outputFile_put << L * (i + 1) << " " << meanTotal_put << " " << stdDevTotal_put << endl;
    }

    outputFile_call.close();
    outputFile_put.close();

    // === SECONDO METODO: simulazione con divisione T a step ===
    // Inizializzazione variabili progressive
    meanTotal_call = meanTotal_put = 0.; 
    stdDevTotal_call = stdDevTotal_put = 0.;
    sumMeanBlocks_call = sumMeanBlocks_put = 0.;
    meanBlock_call = meanBlock_put = 0.;
    sumSquaredMeanBlocks_call = sumSquaredMeanBlocks_put = 0.;

    ofstream outputFile2_call("data_step_call.dat");
    ofstream outputFile2_put("data_step_put.dat");

    if (!outputFile2_call.is_open() || !outputFile2_put.is_open()) {
        cerr << "Errore: impossibile aprire i file" << endl;
        return -1;
    }

    for (int i = 0; i < N; i++) {
        sumBlock_call = 0.;
        sumBlock_put = 0.;

        for (int j = 0; j < L; j++) {
            S = S0; // Reset prezzo iniziale

            // Evoluzione del prezzo in 100 step temporali
            for(int l =  0; l < 100; l++){
                z = rnd.Gauss(0., 1.);
                S *= exp((mu-0.5*sigma*sigma)*0.01 + sigma*z*sqrt(0.01));
            }

            // Payoff scontato
            sumBlock_call += exp(-mu*T)*max(0., S - K);
            sumBlock_put += exp(-mu*T)*max(0., K - S);
        }

        // Medie e deviazioni progressive
        meanBlock_call = sumBlock_call * invL; 
        sumMeanBlocks_call += meanBlock_call;
        sumSquaredMeanBlocks_call += meanBlock_call * meanBlock_call; 

        meanBlock_put = sumBlock_put * invL; 
        sumMeanBlocks_put += meanBlock_put;
        sumSquaredMeanBlocks_put += meanBlock_put * meanBlock_put; 

        meanTotal_call = sumMeanBlocks_call / (i + 1);
        stdDevTotal_call = sqrt((sumSquaredMeanBlocks_call / (i + 1) - (meanTotal_call * meanTotal_call)) / (i + 1));

        meanTotal_put = sumMeanBlocks_put / (i + 1);
        stdDevTotal_put = sqrt((sumSquaredMeanBlocks_put / (i + 1) - (meanTotal_put * meanTotal_put)) / (i + 1));

        outputFile2_call << L * (i + 1) << " " << meanTotal_call << " " << stdDevTotal_call << endl;
        outputFile2_put << L * (i + 1) << " " << meanTotal_put << " " << stdDevTotal_put << endl;
    }

    outputFile2_call.close();
    outputFile2_put.close();

    return 0;
}