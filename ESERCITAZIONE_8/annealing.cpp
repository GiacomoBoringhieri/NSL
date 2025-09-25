#include "annealing.h"
#include "simulator.h"
#include <fstream>
#include <cmath>
#include <random>

SimulatedAnnealing::SimulatedAnnealing(Psi& psi_, MCSimulator& sim_, double T0, double dmu, double dsigma, double T_min, double cooling)
    : psi(psi_), simulator(sim_), T(T0), mu_step(dmu), sigma_step(dsigma), Tmin(T_min), cooling_rate(cooling) {}

void SimulatedAnnealing::optimize(int sa_steps, int mc_steps_per_sa, int mc_blocks_per_sa, bool ratio) {

    std::ofstream out("OUTPUT/sa_trajectory.dat");
    std::ofstream par_out("OUTPUT/parameters_trajectory.dat");
    out << "#step" << " " << "T" << " " << "mu" << " " << "sigma" << " " << "energy" << " " << "err_energy" <<"\n";
    par_out << "#Mu" << " " << "Sigma" <<"\n";
    double current_mu = psi.getMu();   
    double current_sigma = psi.getSigma();

    par_out << current_mu << " " << current_sigma << "\n";

    int accepted = 0;
    int rejected = 0;

    for (int step = 0; step < sa_steps && T > Tmin; ++step){
        // 1. Genera nuovi parametri candidati
        double new_mu = current_mu + simulator.getRND().Rannyu(-mu_step, mu_step);
        double new_sigma = current_sigma + simulator.getRND().Rannyu(-sigma_step, sigma_step);
        if(new_sigma <= 0.) new_sigma = current_sigma;  // Evita sigma <= 0

        // 2. Valuta l'energia con i vecchi e nuovi parametri
        psi.setParameters(current_mu, current_sigma);
        MCResult old_res = simulator.compute_energy(mc_steps_per_sa, mc_blocks_per_sa);  

        psi.setParameters(new_mu, new_sigma);
        MCResult new_res = simulator.compute_energy(mc_steps_per_sa, mc_blocks_per_sa);

        // 3. Calcola la probabilità di accettazione (Metropolis)
        double delta_E = new_res.mean_energy - old_res.mean_energy;
        double acceptance_prob = (delta_E < 0.) ? 1.0 : exp(- delta_E / T);

        if (simulator.getRND().Rannyu() < acceptance_prob) {
            accepted++;
            current_mu = new_mu;
            current_sigma = new_sigma;
            out << step << " " << T << " " << current_mu << " " << current_sigma << " " << new_res.mean_energy << " " << new_res.std_dev <<"\n";
            par_out << current_mu << " " << current_sigma  <<"\n";
        } else {
            psi.setParameters(current_mu, current_sigma);
            rejected++;
            out << step << " " << T << " " << current_mu << " " << current_sigma << " " << old_res.mean_energy << " " << old_res.std_dev <<"\n";
            par_out << current_mu << " " << current_sigma  <<"\n";
        }
        
        // 5. Raffredda il sistema
        T *= cooling_rate;
    }

    out.close();
    par_out.close();

    par_ratio = double(accepted) / (accepted + rejected);

    if(ratio){
        std::cerr << "\nPERCENTUALE DI ACCETTAZIONE NELL'OTTIMIZZAZIONE: " << par_ratio * 100. << " %" << std::endl;
    }


}


void SimulatedAnnealing::saveTrajectory(const std::string& filename) {
    // (Opzionale: implementare se si vogliono salvare più dati)
}