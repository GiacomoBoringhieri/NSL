#ifndef __CHROMOSOME_H__
#define __CHROMOSOME_H__

#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>

#include "random.h"

using namespace std;
using namespace arma;

// Il gene è la città
// OSSERVAZIONE: QUELLO CHE IO HO CHIAMATO GENE IN SENSO BIOLOGICO SAREBBE L'ALLENE
// IL GENE È: IL GENE COLORE DEGLI OCCHI
// L'ALLENE È: AZZURRO, VERDE, ...
// QUI IL GENE È LA POSIZIONE NUMERO n, L'ALLENE È L'INTERO CHE LA RICOPRE
struct Gene{
        // Coordinate Cartesiane
        double x, y;
        // Indice
        int index;
    };

// Il cromosoma è il percorso
class Chromosome{
    private:
        // Lunghezza del cromosoma
        int _genes_number;
        // Cromosoma fatto di geni
        arma::field<Gene> _genes;
        // Lunghezza del percorso, è il valore loss function di quel cromosoma
        double _lenght;

    public:

        // Costruttore vuoto
        Chromosome() {};

        // Costruttore con lunghezza e file per lettura
        Chromosome(int n, string filename);
        // Getter per l'array di geni
        field<Gene> &get_genes() { return _genes; }
        // Resituisce la lunghezza del percorso di un cromosoma, che è nel private
        const double get_lenght() const {
            return _lenght;
        }

        // OPERAZIONI

        // Controlla che ogni città venga visitata una e una sola volta, no indici ripetuti tra i geni del cromosoma
        bool check();
        // Mi valuta la loss function, È VERO CHE IL FIT È UNA COSA PROPRIA DEL CROMOSOMA MA LA MATRICE È IN PUPOLATION
        double evaluateFit();
        // Stampa su file il percorso di questo cromosoma
        void printPath(string filename);
        // Stampa su output il percorso di questo cromosoma
        void printString();
        // Muta l'individuo stesso attraverso una permutazione
        void mutate_Permutation(Random &rnd);
        // Muta l'individuo stesso attraverso uno shift di n posizioni
        void mutate_Shift(Random &rnd);
        // Permutazione con un altro blocco di città
        void mutate_Blocks(Random &rnd);
        // Inversione in cui appaiono
        void mutate_Inversion(Random &rnd);
        // Ha sensazione di quale sia il più vicino a lui e il prossimo è quello
        // void mutate_Sense(mat matrix, Random &rnd);

        // Muta con certa probabilità nei vari modi desctritti nelle altre funzioni
        void mutate(Random &rnd);


};

// METODI DI MUTAZIONE CHE FUNZIONANO: Inversion, Permutation
// METODI DI MUTAZIONE CHE NON FUNZIONANO: Shift, Blocks (da solo no problemi ma in mutate fa rogne)


#endif //__CHROMOSOME_H__