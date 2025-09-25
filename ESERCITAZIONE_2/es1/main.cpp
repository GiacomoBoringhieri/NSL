#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../PRNG/random.h"
#include "../../PRNG/funzioni.h"

using namespace std;

int main(int argc, char *argv[]){

    int M=100000;              // Numero totale di punti per l'integrazione Monte Carlo
    int N=100;                 // Numero di blocchi
    Random rnd;                // Primo generatore di numeri casuali
    Random rnd_2;              // Secondo generatore di numeri casuali

    P_r_2 P;                   // Funzione per trasformazione con CDF
    Cos_g Cg;                   // Funzione coseno modificata
    Cos C;                      // Funzione coseno semplice

    rnd.Start(2);               // Inizializzazione del primo generatore con seed 2
    rnd_2.Start(2);             // Inizializzazione del secondo generatore con seed 2

    // Stampa 20 numeri casuali generati dal primo generatore (verifica)
    for(int i=0; i<20; i++){
       cout << rnd.Rannyu() << endl;
    }

    // Salva lo stato del primo generatore su file "seed.out"
    rnd.SaveSeed();

    // Integra la funzione Cos usando il metodo dei blocchi e salva su file
    C.integraBlocchi(M, N, rnd, "data_2_1.dat");

    // Integra la funzione Cos_g usando la trasformazione tramite CDF P_r_2 e salva su file
    Cg.integraBlocchiCDF(M, N, rnd_2, P, "data_2_1_CDF.dat");

    return 0;
}