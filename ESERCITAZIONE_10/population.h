#ifndef __POPULATION_H__
#define __POPULATION_H__

#include "chromosome.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Population{
    private:

    int _rank;

    field <Chromosome> _chromosomes;
    uword _number_chromosomes;

    mat lenght_matrix;

    Random rnd_pop;

    pair<double, double> averageFitness;
    pair<double, double> averageFitnessBestHalf;

    double _prob_mutation;
    double _prob_cross;
    double _prob_elitismo;

    double _prob_select_cross;
    double _prob_select_mutation;
    double _prob_select_elitismo;

    public:

    // Costruttore vuoto
    Population() {};
    // Costruttore con cromosoma iniziale
    Population(int n, Chromosome c, int rank);

    // Getter per i cromosomi
    field<Chromosome> &get_chromosomes() { return _chromosomes; }

    // Evolve la popolazione per portarla alla generazione successiva
    void evolve(double cross, double mutation, double pass, double p_cross = 1., double p_mutation = 1.);

    // Ordina la popolazione dal migliore al peggiore
    void sort();

    // Selection Operator
    Chromosome select(double p = 1);

    // CROSSOVER
    void crossover_cut(Chromosome mother, Chromosome father, Chromosome &son, Chromosome &daugher);

    // OPERAZIONI SUL FIT

    // Return il migliore fit della popolazione
    double getBestFit();

    // Return il migliore della popolazione, in una popolazione ordinata è sempre al primo posto
    Chromosome getBest();

    // Retun della media, con errore
    pair<double, double> getAverageFit();

    // Retun della media, con errore. Della metà megliore
    pair<double, double> getAverageFitBestHalf();

    // ESPERIMENTO in cui fa tutto da solo e in base ai cambiamenti si adatta
    void evolve();

    // Funzione per stampare le probabilità
    void printProbabilities() const {
        std::cout << "=== Probabilità della Popolazione ===" << std::endl;
        std::cout << "_prob_mutation         : " << _prob_mutation << std::endl;
        std::cout << "_prob_cross            : " << _prob_cross << std::endl;
        std::cout << "_prob_elitismo         : " << _prob_elitismo << std::endl;
        std::cout << "_prob_select_cross     : " << _prob_select_cross << std::endl;
        std::cout << "_prob_select_mutation  : " << _prob_select_mutation << std::endl;
        std::cout << "_prob_select_elitismo  : " << _prob_select_elitismo << std::endl;
        std::cout << "===================================" << std::endl;
    };

    // Nuovo: Sostituisce il peggior cromosoma con un nuovo percorso
    void replaceWorstWithVector(const vector<int>& new_percorso);
    
    // Nuovo: Sostituisce il peggior cromosoma (versione semplificata per compatibilità)
    void replaceWorst(const vector<int>& new_percorso) { 
        replaceWorstWithVector(new_percorso); 
    }

};

#endif //__POPULATION_H__