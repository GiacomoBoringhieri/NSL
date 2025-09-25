#ifndef MC_SIMULATOR_H
#define MC_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "psi.h"
#include "random.h"

struct MCResult {
    double mean_energy;
    double std_dev;
    std::vector<double> samples;
    double metro_ratio;
    std::vector<double> block_ratios; // ratio di accettazione per ogni blocco

    //Stampa su file energia e suo errore
    void stampa(std::string file_name) {
        std::ofstream outputFile(file_name);
        if (!outputFile.is_open()) {
            std::cerr << "Errore: impossibile aprire il file " << file_name << std::endl;
            return;
        }
        outputFile << mean_energy << " " << std_dev << std::endl;
        outputFile.close();
    }

    void compute_histogram(int num_bins, const std::string& output_filename) {
    if (samples.empty() || num_bins <= 0) return;

    double min_val = *std::min_element(samples.begin(), samples.end());
    double max_val = *std::max_element(samples.begin(), samples.end());

    // Se tutti i valori sono uguali, evitiamo divisione per zero
    double bin_width = (max_val - min_val) / num_bins;
    if (bin_width == 0) bin_width = 1.0;

    std::vector<int> bins(num_bins, 0);

    // Conta le occorrenze in ogni bin
    for (double val : samples) {
        int bin = std::min(static_cast<int>((val - min_val) / bin_width), num_bins - 1);
        bins[bin]++;
    }

    // Normalizzazione: frequenza / (totale campioni * larghezza bin)
    std::ofstream outfile(output_filename);
    for (int i = 0; i < num_bins; ++i) {
        double bin_center = min_val + (i + 0.5) * bin_width;
        double normalized = static_cast<double>(bins[i]) / (samples.size() * bin_width);
        outfile << bin_center << " " << normalized << "\n";
    }
    outfile.close();
}

};

class MCSimulator {
private:
    Psi& psi;                // Riferimento alla funzione d'onda trial
    double delta_step;       // Ampiezza del passo Metropolis
    int equilibration_steps; // Passi di termalizzazione iniziali
    Random rnd;

public:

    MCSimulator(Psi& psi_, double delta_, int equil_steps = 1000, int n = 1);
    MCResult compute_energy(int total_steps, int n_blocks, bool print = false, bool ratio = false, std::string filename = "OUTPUT/block_energies.dat");  // Esegue la simulazione
    bool metropolis(double& x);  // Algoritmo Metropolis
    Random& getRND() {return rnd;};
    
};

#endif