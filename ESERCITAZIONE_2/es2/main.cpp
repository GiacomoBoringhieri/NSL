#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "../../PRNG/random.h"
#include "../../PRNG/funzioni.h"

using namespace std;

// Funzione che fa uno step del random walk discreto a 3D (coordinate x,y,z) usando un dado a 6 facce
void Step(vector<int>& P, Random& rnd){
   double y = rnd.Rannyu();
   int pos = static_cast<int>(y*6); // Scelta della direzione casuale
   P[pos] += 1;                     // Aggiorna la componente scelta
};

// Calcola la distanza dall'origine per un random walk discreto
double Distanza(const vector<int>& P){
    double d = 0.0;
    for (int i = 0; i < 6; i += 2) {
        double diff = (P[i] - P[i + 1]);
        d += diff * diff;
    }
    return sqrt(d);
};

// Funzione che fa uno step del random walk continuo su sfera unitaria
void Step_Continuo(vector<double>& P, Random& rnd){
   double phi = rnd.Rannyu()*2.*M_PI;
   double r = rnd.Rannyu();
   double theta = acos(1. - 2.*r); // Distribuzione uniforme sulla sfera
   P[0] += sin(theta)*cos(phi);     // Aggiorna componente x
   P[1] += sin(theta)*sin(phi);     // Aggiorna componente y
   P[2] += cos(theta);              // Aggiorna componente z
};

// Calcola la distanza dall'origine per random walk continuo
double Distanza_Continuo(const vector<double>& P){
    double d = 0.0;
    for (int i = 0; i < 3; i++) {
        d += P[i]*P[i];
    }
    return sqrt(d);
};

int main(int argc, char *argv[]){

    int M=10000;   // Numero totale di random walk
    int N=100;     // Numero di blocchi per calcolo delle medie
    int I = 100;   // Numero di step per ogni random walk
    int L = M / N; // Numero di RW per blocco

    // Variabili per calcolo medie e deviazioni progressive (discreto e continuo)
    double meanTotal = 0.0, stdDevTotal = 0.0, sumBlock, sumMeanBlocks = 0.0, meanBlock = 0.0, sumSquaredMeanBlocks = 0.0;
    double invL = 1.0 / L;

    double meanTotal_c = 0.0, stdDevTotal_c = 0.0, sumBlock_c, sumMeanBlocks_c = 0.0, meanBlock_c = 0.0, sumSquaredMeanBlocks_c = 0.0;

    Random rnd;
    rnd.Start(2);          // Inizializza il generatore casuale
    rnd.SaveSeed();        // Salva stato del generatore

    // Inizializza tutti i random walk all'origine
    vector<vector<int>> P (M, vector<int>(6, 0));        // RW discreto
    vector<vector<double>> P_c (M, vector<double>(3, 0.0)); // RW continuo

    ofstream outputFile("data.dat");    // File output discreto
    ofstream outputFile_c("data_c.dat"); // File output continuo

    if (!outputFile.is_open() || !outputFile_c.is_open()) {
        cerr << "Errore: impossibile aprire i file di output" << endl;
        return 0;
    }

    // Ciclo sugli step dei random walk
    for(int i = 0; i < I; i++){

        meanBlock = meanBlock_c = 0.0;
        sumMeanBlocks = sumMeanBlocks_c = 0.0;
        sumSquaredMeanBlocks = sumSquaredMeanBlocks_c = 0.0;

        // Ciclo sui blocchi
        for(int k = 0; k < N; k++){

            sumBlock = sumBlock_c = 0.0;

            // Ciclo sui random walk all'interno del blocco
            for(int j = 0; j < L; j++){

                Step(P[k*N + j], rnd);                    // Step RW discreto
                sumBlock += Distanza(P[k*N + j]);        // Aggiorna distanza dall'origine

                Step_Continuo(P_c[k*N + j], rnd);       // Step RW continuo
                sumBlock_c += Distanza_Continuo(P_c[k*N + j]); // Aggiorna distanza dall'origine
            }

            // Calcolo medie dei blocchi
            meanBlock = sumBlock * invL; 
            sumMeanBlocks += meanBlock;
            sumSquaredMeanBlocks += meanBlock * meanBlock; 

            meanBlock_c = sumBlock_c * invL; 
            sumMeanBlocks_c += meanBlock_c;
            sumSquaredMeanBlocks_c += meanBlock_c * meanBlock_c; 
        }

        // Medie progressive e deviazioni standard (discreto)
        meanTotal = sumMeanBlocks / N;
        stdDevTotal = sqrt((sumSquaredMeanBlocks / N - (meanTotal * meanTotal)) / N);
        outputFile << i + 1 << " " << meanTotal << " " << stdDevTotal << endl;

        // Medie progressive e deviazioni standard (continuo)
        meanTotal_c = sumMeanBlocks_c / N;
        stdDevTotal_c = sqrt((sumSquaredMeanBlocks_c / N - (meanTotal_c * meanTotal_c)) / N);
        outputFile_c << i + 1 << " " << meanTotal_c << " " << stdDevTotal_c << endl;
    }
   
    outputFile.close();
    outputFile_c.close();

    return 0;
}
