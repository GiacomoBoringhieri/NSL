#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "psi.h"
#include "simulator.h"
#include "annealing.h"
#include "random.h"

using namespace std;

int main() {

    // Parametri iniziali
    const double mu0 = 1.0;
    const double sigma0 = 0.5;
    const double delta_step = 1.5;       // Passo per Metropolis
    
    // Parametri Simulated Annealing
    const double T0 = 5.0;              // Temperatura iniziale
    const double Tmin = 1E-8;           // Temperatura minima
    const double cooling_rate = 0.95;   // Cooling Rate
    const double mu_step = 0.05;        // Passo per mu
    const double sigma_step = 0.01;     // Passo per sigma
    
    // Passi di simulazione
    const int sa_steps = 100000;        // Passi SA
    const int mc_steps_per_sa = 1000000; // Passi MC per ogni passo SA
    const int mc_blocks_per_sa = 1000;  // Blocchi MC per ogni passo SA
    const int equil_steps = 1000;       // Passi di equilibrazione

    // Inizializzazione
    Psi psi(mu0, sigma0);
    MCSimulator simulator(psi, delta_step, equil_steps);
    SimulatedAnnealing sa(psi, simulator, T0, mu_step, sigma_step, Tmin, cooling_rate);

    // Esegui ottimizzazione
    cout << "\nStarting optimization with initial parameters:" << endl;
    cout << "mu = " << mu0 << ", sigma = " << sigma0 << endl;

    cout << "Calcolo enegia con questi parametri iniziali: " << endl;
    MCResult forst_res = simulator.compute_energy(10000, 100, true, true, "OUTPUT/initial_energy.dat");
    cout << "Finito" << endl;
    
    cout << "Inizio ottimizzazione: ..." << endl;
    
    sa.optimize(sa_steps, mc_steps_per_sa, mc_blocks_per_sa, true);

    // Risultati finali
    cout << "\nOptimization completed. Final parameters:" << endl;
    cout << "mu = " << psi.getMu() << ", sigma = " << psi.getSigma() << endl;

    // Esegui una simulazione finale più accurata
    cout << "\nRunning final accurate simulation..." << endl;
    MCResult final_res = simulator.compute_energy(10000, 100, true, true);
    cout << "Final energy: " << final_res.mean_energy << " ± " << final_res.std_dev << endl;

    // Salva risultati finali
    final_res.stampa("OUTPUT/final_energy.dat");

    //Stampa su file dati per l'istogramma
    final_res.compute_histogram(50, "OUTPUT/psi_histogram.dat");

    //Mi dice la percentuale di accettati del metropolis
    cout << "Percentuale di accettati metropolis sul totale nell'ultima simulazione: " << final_res.metro_ratio * 100. << "%" << endl;

    return 0;
}