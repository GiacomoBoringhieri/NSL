#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include "psi.h"
#include "simulator.h"
#include <vector>

class SimulatedAnnealing {
private:
    Psi& psi;               // Riferimento alla funzione d'onda da ottimizzare
    MCSimulator& simulator; // Simulatore MC per valutare l'energia
    double T;               // Temperatura corrente
    double mu_step;         // Passo per perturbare mu
    double sigma_step;      // Passo per perturbare sigma
    double Tmin;            // Temperatura minima
    double cooling_rate;    // Fattore di raffreddamento

    double par_ratio;

public:
    SimulatedAnnealing(Psi& psi_, MCSimulator& sim_, double T0, double dmu, double dsigma, double T_min = 1E-5, double cooling = 0.95);
    void optimize(int sa_steps, int mc_steps_per_sa, int mc_blocks_per_sa, bool ratio = false);  // Esegue l'ottimizzazione
    void saveTrajectory(const std::string& filename);  // Salva la storia dell'ottimizzazione
};

#endif