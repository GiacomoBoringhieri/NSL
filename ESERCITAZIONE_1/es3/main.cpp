#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../PRNG/random.h"

using namespace std;

// Funzione che stima π usando il metodo di Buffon per un blocco di L lanci
// Se 'alternativo' è true, usa un metodo alternativo per scegliere l'angolo theta
double stima_pi_blocco(Random &rnd, int L, double L_ago, double d, bool alternativo) {
    int hit = 0; // Contatore per gli aghi che attraversano una linea

    for (int j = 0; j < L; j++) {
        double x = rnd.Rannyu() * d; // Posizione iniziale dell'ago
        double theta = alternativo ? 2.0 * acos(rnd.Rannyu()) : rnd.Rannyu() * M_PI; // Angolo dell'ago

        double x1 = x - (L_ago / 2.0) * cos(theta); // Estremità 1 dell'ago
        double x2 = x + (L_ago / 2.0) * cos(theta); // Estremità 2 dell'ago

        // Incrementa hit se l'ago attraversa una linea
        if (floor(x1 / d) != floor(x2 / d)) hit++;
    }

    if (hit == 0) return 0.0; // Evita divisione per zero
    return 2.0 * L_ago * L / (hit * d); // Formula di Buffon per stimare π
}

int main(int argc, char *argv[]) {
    const int M = 100000;  // Numero totale di lanci
    const int N = 100;     // Numero di blocchi
    const int L = M / N;   // Lanci per blocco

    const double L_ago = 1.0; // Lunghezza ago
    const double d = 2.0;     // Distanza tra le linee

    Random rnd;
    rnd.Start(1);              // Inizializzazione generatore casuale

    ofstream outputFile("data.dat");      // File per il metodo tradizionale
    ofstream outputFileAlt("_data.dat");  // File per il metodo alternativo

    if (!outputFile.is_open() || !outputFileAlt.is_open()) {
        cerr << "Errore: impossibile aprire i file di output" << endl;
        return 1;
    }

    // Variabili per medie progressive e deviazioni standard
    double sumMean = 0.0, sumSquaredMean = 0.0;
    double sumMeanAlt = 0.0, sumSquaredMeanAlt = 0.0;

    for (int i = 0; i < N; i++) {
        // Stima π per il blocco corrente (metodo tradizionale e alternativo)
        double pi_block = stima_pi_blocco(rnd, L, L_ago, d, false);
        double pi_block_alt = stima_pi_blocco(rnd, L, L_ago, d, true);

        // Aggiornamento medie progressive e deviazioni standard (tradizionale)
        sumMean += pi_block;
        sumSquaredMean += pi_block * pi_block;
        double mean = sumMean / (i + 1);
        double error = (i == 0) ? 0.0 :
            sqrt((sumSquaredMean / (i + 1) - mean * mean) / i);

        outputFile << L * (i + 1) << " " << mean << " " << error << endl;

        // Aggiornamento medie progressive e deviazioni standard (alternativo)
        sumMeanAlt += pi_block_alt;
        sumSquaredMeanAlt += pi_block_alt * pi_block_alt;
        double meanAlt = sumMeanAlt / (i + 1);
        double errorAlt = (i == 0) ? 0.0 :
            sqrt((sumSquaredMeanAlt / (i + 1) - meanAlt * meanAlt) / i);

        outputFileAlt << L * (i + 1) << " " << meanAlt << " " << errorAlt << endl;
    }

    outputFile.close();
    outputFileAlt.close();

    return 0;
}