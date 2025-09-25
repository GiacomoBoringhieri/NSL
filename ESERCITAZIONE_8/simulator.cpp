#include "simulator.h"

//Costruttore
MCSimulator::MCSimulator(Psi& psi_, double delta_, int equil_steps, int n) 
    : psi(psi_), delta_step(delta_), equilibration_steps(equil_steps) {rnd.Start(n);}

// Implementazione dell'algoritmo Metropolis
bool MCSimulator::metropolis(double &x) {
    // Proposta di nuovo stato
    double x_new = x + rnd.Rannyu(- delta_step, delta_step);

    // Calcolo delle probabilità
    double p_old = psi.abs_val(x);
    double p_new = psi.abs_val(x_new);

    // Probabilità di accettazione
    double A = std::min(1., p_new/p_old);

    // Decisione di accettare o rifiutare
    if(rnd.Rannyu() < A) {x = x_new; return true;}
        else return false;

};

//Restituisce MCResult con energia e suo errore e anche il ratio Metropolis
MCResult MCSimulator::compute_energy(int total_steps, int n_blocks, bool print, bool ratio, std::string file_name) {
    MCResult result;
    int block_size = total_steps / n_blocks;
    double x = 0.0; // Posizione iniziale
    int accepted = 0;
    int rejected = 0;
    std::ofstream out_file;
    if (print) {
        out_file.open(file_name);
        if (!out_file)
            throw std::runtime_error("Impossibile aprire il file per la scrittura.");
    }

    // ====== AUTOTUNING DEL DELTA_STEP ======
    double target_ratio = 0.5;
    for (int tune = 0; tune < 10; ++tune) { // ripeto qualche iterazione di tuning
        int acc = 0, rej = 0;
        double xt = x;
        for (int i = 0; i < 1000; ++i) { // 1000 tentativi per stimare il ratio
            if (metropolis(xt)) acc++;
            else rej++;
        }
        double r = double(acc) / (acc + rej);
        if (r > target_ratio) {
            delta_step *= 1.05; // aumenta un po'
        } else {
            delta_step *= 0.95; // diminuisci un po'
        }
    }

    // Termalizzazione
    for (int i = 0; i < equilibration_steps; ++i) {
        metropolis(x);
    }

    // Simulazione principale
    std::vector<double> block_energies(n_blocks);
    double progressive_sum = 0.0;
    double progressive_sum_sq = 0.0;
    
    // Variabili per l'accettazione progressiva
    int progressive_accepted = 0;
    int progressive_rejected = 0;

    for (int block = 0; block < n_blocks; ++block) {
        double energy_sum = 0.0;
        int block_accepted = 0; // Contatore per questo blocco
        int block_rejected = 0; // Contatore per questo blocco

        for (int step = 0; step < block_size; ++step) {
            bool metro = metropolis(x);
            if(metro) {
                accepted++;
                block_accepted++;
            } else {
                rejected++;
                block_rejected++;
            }
            energy_sum += psi.loc_energy(x);
            result.samples.push_back(x);
        }

        double block_energy = energy_sum / block_size;
        block_energies[block] = block_energy;

        // Calcolo del ratio per questo blocco
        double block_ratio = 0.0;
        if (block_accepted + block_rejected > 0) {
            block_ratio = double(block_accepted) / double(block_accepted + block_rejected);
        }

        // Aggiorna i contatori progressivi
        progressive_accepted += block_accepted;
        progressive_rejected += block_rejected;

        // Calcolo dell'accettazione media progressiva
        double progressive_acceptance_ratio = 0.0;
        if (progressive_accepted + progressive_rejected > 0) {
            progressive_acceptance_ratio = double(progressive_accepted) / double(progressive_accepted + progressive_rejected);
        }

        // Progressiva per l'energia
        progressive_sum += block_energy;
        progressive_sum_sq += block_energy * block_energy;
        double mean = progressive_sum / (block + 1);
        double err = 0.0;
        if (block > 0) {
            err = std::sqrt((progressive_sum_sq / (block + 1) - mean * mean) / block);
        }

        if (print) {
            // Formato: blocco | media progressiva | errore | ratio del blocco | accettazione media progressiva
            out_file << block << '\t' << mean << '\t' << err << '\t' << block_ratio << '\t' << progressive_acceptance_ratio << '\n';
        }
    }

    if(ratio) result.metro_ratio = double(accepted) / double(accepted + rejected);

    // Statistiche finali
    double sum = 0.0, sum_sq = 0.0;
    for (double e : block_energies) {
        sum += e;
        sum_sq += e * e;
    }
    result.mean_energy = sum / n_blocks;
    result.std_dev = std::sqrt((sum_sq / n_blocks - result.mean_energy * result.mean_energy) / n_blocks);

    if (print) {
        out_file.close();
    }

    return result;
}