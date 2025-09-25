#include <iostream>
#include <fstream>
#include <string>
#include "../../PRNG/random.h"
#include "../../PRNG/funzioni.h"

using namespace std;

// Funzione che calcola il valore del test chi-quadro
// confrontando le frequenze osservate con una frequenza attesa costante
double chi_quadro(const std::vector<int>& observed, double expected);

int main(int argc, char *argv[]){
    
   int M=100000;              // Numero totale di lanci per l'integrazione Monte Carlo
   int N=100;                 // Numero di blocchi per il metodo dei blocchi
   Random rnd;                // Oggetto generatore di numeri casuali
   
   rnd.Start(2);              // Inizializzazione del generatore con seed specifico
   rnd.SaveSeed();            // Salvataggio dello stato corrente del generatore su file

   Linear f;                  // Funzione lineare f(x) = x
   var sigma;                 // Funzione della varianza (x-0.5)^2

   // === PRIMO ESERCIZIO ===
   // Integra la funzione lineare a blocchi, salva i risultati su "data_1_1.dat"
   f.integraBlocchi(M, N, rnd, "data_1_1.dat");

   // === SECONDO ESERCIZIO ===
   // Integra la funzione della varianza a blocchi, salva i risultati su "data_1_2.dat"
   sigma.integraBlocchi(M, N, rnd, "data_1_2.dat");

   // === TERZO ESERCIZIO ===
   Random rnd_chi;             // Nuovo generatore per il test chi-quadro
   rnd_chi.Start(3);           // Inizializzazione con seed diverso

   M = 100;                    // Numero di sottointervalli (bins) per il chi-quadro
   int n = 10000;              // Numero di lanci per ogni test chi-quadro

   ofstream fout_chi("data_chi.dat"); // File di output dei risultati

   for (int i = 0; i < 100; i++) { // 100 esperimenti indipendenti
      vector<int> conteggi(M, 0);  // Inizializza i conteggi dei bin a zero

      for (int j = 0; j < n; j++) {
         // Determina in quale bin cade il numero casuale uniforme [0,1)
         int index = floor(rnd_chi.Rannyu() * M);
         conteggi[index]++;        // Incrementa il conteggio del bin corrispondente
      }

      // Calcola il chi-quadro rispetto alla distribuzione uniforme attesa (n/M per bin)
      // e salva su file insieme all'indice del test
      fout_chi << i+1 << "\t" << chi_quadro(conteggi, n/M) << endl;
   }

    return 0;
}

// Implementazione della funzione chi-quadro
// Somma ((osservato - atteso)^2 / atteso) su tutti i bin
double chi_quadro(const std::vector<int>& observed, double expected) {
    double chi = 0;
    for (int i = 0; i < (int)observed.size(); i++) {
        chi += pow(observed[i] - expected, 2) / expected;
    }
    return chi;
}