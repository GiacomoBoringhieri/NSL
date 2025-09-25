/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "system.h"

using namespace std;

// Funzione per salvare acceptance finale
void save_acceptance(const string &source_file, const string &target_file, double T) {
    ifstream in(source_file);
    ofstream out(target_file, ios::app); // append

    if (!in.is_open() || !out.is_open()) {
        cerr << "Errore apertura file acceptance" << endl;
        return;
    }

    string line, last_line;
    while (getline(in, line)) {
        if (!line.empty())
            last_line = line;
    }

    // estrai l'ultima colonna
    stringstream ss(last_line);
    double dummy, acceptance;
    ss >> dummy >> acceptance; // assumendo due colonne

    out << T << " " << acceptance << endl;
}

int main(int argc, char *argv[]) {
    int nconf = 1;
    System SYS;

    // Dummy initialize per leggere tipo e campo
    SYS.initialize(1.0);
    int type = SYS.getType();
    double H = SYS.get_H();
    SYS.finalize();

    if (type != 2 && type != 3) {
        cerr << "Tipo sistema non riconosciuto" << endl;
        return 1;
    }

    string base_output = (type == 2) ? "../OUTPUT/T_FUNCTIONS_METRO/" : "../OUTPUT/T_FUNCTIONS_GIBBS/";

    // Prepara lista file da azzerare
    vector<string> files_to_clear;

    if (H == 0) {
        files_to_clear = {"specific_heat_T.dat", "susceptibility_T.dat", "total_energy_T.dat"};
    } else {
        files_to_clear = {"magnetization_T.dat"};
    }

    if (type == 2) files_to_clear.push_back("acceptance.dat");

    for (auto &fname : files_to_clear)
        ofstream(base_output + fname).close();

    // Ciclo temperature
    for (double T = 0.5; T <= 2.0; T += 0.05) {
        SYS.initialize(T);
        SYS.initialize_properties();
        SYS.block_reset(0);

        for (int i = 0; i < SYS.get_nbl(); i++) {
            for (int j = 0; j < SYS.get_nsteps(); j++) {
                SYS.step();
                SYS.measure();
            }
            SYS.averages(i + 1);
            SYS.block_reset(i + 1);
        }

        SYS.finalize();

        // Salvataggio condizionale
        if (H == 0) {
            SYS.print_last_two_columns("../OUTPUT/specific_heat.dat", base_output + "specific_heat_T.dat", T);
            SYS.print_last_two_columns("../OUTPUT/susceptibility.dat", base_output + "susceptibility_T.dat", T);
            SYS.print_last_two_columns("../OUTPUT/total_energy.dat", base_output + "total_energy_T.dat", T);
        } else {
            SYS.print_last_two_columns("../OUTPUT/magnetization.dat", base_output + "magnetization_T.dat", T);
        }

        // Acceptance solo per Metropolis
        if (type == 2) {
            save_acceptance("../OUTPUT/acceptance.dat", base_output + "acceptance.dat", T);
        }
    }

    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
