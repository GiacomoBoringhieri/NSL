#include <iostream>
#include <fstream>
#include <string>
#include "../../PRNG/random.h"
#include "../../PRNG/funzioni.h"

using namespace std;

int main(int argc, char *argv[]){

   Random rnd;                // Oggetto generatore di numeri casuali
   rnd.Start(1);              // Inizializzazione con seed 1

   // Vettore che indica il numero di variabili indipendenti da sommare per ciascun esperimento
   vector<int> N = {1, 2, 10, 100};

   // Parametri per le distribuzioni esponenziale e lorentziana
   double lambda {1}; // parametro della distribuzione esponenziale
   double mu {0};     // media della lorentziana
   double gamma {1};  // larghezza della lorentziana

   // ======================================================
   // Generazione dei dati
   // ======================================================

   for (int i = 0; i < (int) N.size(); i++) {

      // Nome del file che conterrà i dati relativi a N[i] variabili sommate
      string filename = "data_N" + std::to_string(N[i]) + ".dat";
      ofstream fout(filename);

      // Ogni riga del file conterrà tre colonne:
      // media di N[i] variabili uniformi, esponenziali e lorentziane
      if (fout.is_open()) {

         for (int j = 0; j < 10000; j++) { // Genero 10^4 campioni per ciascun N
            double sum_std = 0; // somma delle N[i] variabili uniformi in [0,1)
            double sum_exp = 0; // somma delle N[i] variabili esponenziali
            double sum_lor = 0; // somma delle N[i] variabili lorentziane

            // Ciclo sulle variabili da sommare
            for (int k = 0; k < N[i]; k++) {
               sum_std += rnd.Rannyu();               // Genera uniforme [0,1)
               sum_exp += rnd.Exponential(lambda);    // Genera esponenziale
               sum_lor += rnd.Lorentz(mu, gamma);     // Genera lorentziana
            }

            // Scrivo nel file le medie delle variabili sommate
            fout << sum_std / N[i] << "\t"
                 << sum_exp / N[i] << "\t"
                 << sum_lor / N[i] << endl;
         }

         fout.close(); // Chiudo il file

      } else {
         // Messaggio di errore se il file non si apre correttamente
         cerr << "Errore nell'apertura del file: " << filename << endl;
      }
   }

   return 0;
}