#include "random.h"
#include "chromosome.h"
#include "population.h"
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <vector>
using namespace std;

int main (int argc, char *argv[]){
    int rank, size;
    
    MPI_Init(&argc, &argv); // inizializzazione MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    Random _rnd;
    _rnd.Start(1);
    
    // Parametri del problema
    int lunghezza_cromosoma = 110;
    int lunghezza_popolazione = 500;
    int number_generations = 10000;
    int n_migrazione = 200; // Ogni quante generazioni fare la migrazione
    bool auto_adattamento = true;
    bool prov_ita = true;
    string input_file = "cap_prov_ita.dat";
    string prefix = "cap_prov";
    
    if (rank == 0) {
        cout << "=== PARAMETRI SIMULAZIONE ===" << endl;
        cout << "Avvio evoluzione parallela con " << size << " processi" << endl;
        cout << "USO DEL FILE INIZIALE: " << input_file << endl;
        cout << "Numero processi MPI: " << size << endl;
        cout << "Lunghezza cromosoma: " << lunghezza_cromosoma << endl;
        cout << "Popolazione per processo: " << lunghezza_popolazione << endl;
        cout << "Generazioni: " << number_generations << endl;
        cout << "Migrazione ogni: " << n_migrazione << " generazioni" << endl;
        if (auto_adattamento) cout << "SIMULAZIONE CON AUTO-ADATTAMENTO DELLA POPOLAZIONE" << endl;
        cout << "==============================" << endl;
    }
    
    // Inizializzo il generatore random (diverso per ogni processo)
    Random rnd;
    rnd.Start(rank + 1); // Seed diverso per ogni processo
    
    // Genero le città (stesso per tutti i processi)
    rnd.city_circle(lunghezza_cromosoma);
    rnd.city_box(lunghezza_cromosoma);
    
    // Creo il cromosoma iniziale Adam usando il file scelto
    Chromosome Adam(lunghezza_cromosoma, input_file);
    
    if (rank == 0) {
        cout << "CROMOSOMA INIZIALE: \n";
        Adam.printString();
    }
    
    // Creo la popolazione iniziale (diversa per ogni processo)
    Population generation_0(lunghezza_popolazione, Adam, rank);
    generation_0.sort();

    string rank_file = "RANK/";
    
    // File di output specifici per ogni processo
    string best_half_file = rank_file + prefix + "_rank" + to_string(rank) + 
                           (auto_adattamento ? "_best_half_auto.dat" : "_best_half.dat");
    string champion_path_file = rank_file +prefix + "_rank" + to_string(rank) + 
                               (auto_adattamento ? "_best_path_final_auto.dat" : "_best_path_final.dat");
    string champions_file = rank_file + prefix + "_rank" + to_string(rank) + 
                           (auto_adattamento ? "_champions_auto.dat" : "_champions.dat");
    
    // Apro i file per salvare i risultati
    ofstream off_file(best_half_file);
    ofstream off_champion_file(champions_file);
    
    if (!off_file.is_open()) {
        cerr << "Processo " << rank << ": Unable to open file " << best_half_file << endl;
        MPI_Finalize();
        return -1;
    }
    
    // Evoluzione per generazioni con migrazione periodica
    for (int i = 0; i < number_generations; i++) {
        
        // Evoluzione locale
        if (auto_adattamento)
            generation_0.evolve();
        else
            generation_0.evolve(0.7, 0.5, 0.2);
        
        generation_0.sort();
        
        // Migrazione circolare ogni n_migrazione generazioni
        if (i > 0 && i % n_migrazione == 0 && size > 1) {
            
            // Ottengo il miglior cromosoma come vector<int>
            vector<int> best_percorso = generation_0.getBest().getAsVector();
            
            // Buffer per ricevere
            vector<int> recv_percorso(lunghezza_cromosoma);
            
            // Calcolo i ranks per lo scambio circolare
            int send_to = (rank + 1) % size;      // processo successivo
            int recv_from = (rank - 1 + size) % size;  // processo precedente
            
            // Scambio circolare
            MPI_Sendrecv(best_percorso.data(), lunghezza_cromosoma, MPI_INT, send_to, 0,
                        recv_percorso.data(), lunghezza_cromosoma, MPI_INT, recv_from, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Sostituisco il peggior cromosoma con quello ricevuto
            generation_0.replaceWorst(recv_percorso);
            generation_0.sort();
            
            if (rank == 0) {
                cout << "GEN " << i << " - Migrazione circolare completata" << endl;
            }
        }
        
        // Salvataggio e stampa ogni 100 generazioni
        if (i % 100 == 0) {
            auto fit_half = generation_0.getAverageFitBestHalf();
            
            if (rank == 0) {
                cout << "GEN " << i << " -> MEDIA: " << fit_half.first
                     << " PM " << fit_half.second << endl;
            }
            
            off_file << i << " " << fit_half.first << " " << fit_half.second << "\n";
            off_champion_file << i << " " << generation_0.getBestFit() << "\n";
        }
    }
    
    off_file.close();
    off_champion_file.close();
    
    // Stampa delle probabilità finali
    if (rank == 0) {
        generation_0.printProbabilities();
    }
    
    // Salvo il miglior percorso locale di ogni processo
    generation_0.getBest().printPath(champion_path_file);
    
    cout << "Processo " << rank << " - Evoluzione completata" << endl;
    cout << "CAMPIONE LOCALE PROCESSO " << rank << ":" << endl;
    generation_0.getBest().printString();
    cout << "Migliore fit locale: " << generation_0.getBestFit() << endl;
    
    auto fit_final = generation_0.getAverageFitBestHalf();
    cout << "MEDIA FINALE: " << fit_final.first << " PM " << fit_final.second << endl;
    
    if (generation_0.getBest().check())
        cout << "CAMPIONE VALIDO" << endl;
    
    // Sincronizzazione prima della raccolta finale
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Preparazione per l'invio del miglior cromosoma usando i nuovi metodi
    vector<int> best_chromosome_out = generation_0.getBest().getAsVector();
    double fitness_out = generation_0.getBestFit();
    
    // Raccolta dei migliori risultati sul processo 0
    if (rank == 0) {
        // Vettori per memorizzare i risultati di tutti i processi
        vector<vector<int>> all_best_chromosomes(size, vector<int>(lunghezza_cromosoma));
        vector<double> all_best_fitness(size);
        
        // Il rank 0 salva i suoi risultati
        all_best_chromosomes[0] = best_chromosome_out;
        all_best_fitness[0] = fitness_out;
        
        // Riceve i risultati dagli altri processi
        for (int l = 1; l < size; l++) {
            MPI_Recv(all_best_chromosomes[l].data(), lunghezza_cromosoma, MPI_INT, l, 0, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&all_best_fitness[l], 1, MPI_DOUBLE, l, 1, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        // Trova il miglior risultato globale (fitness minimo)
        int best_global_rank = 0;
        double best_global_fitness = all_best_fitness[0];
        
        for (int k = 1; k < size; k++) {
            if (all_best_fitness[k] < best_global_fitness) {
                best_global_fitness = all_best_fitness[k];
                best_global_rank = k;
            }
        }
        
        // Crea un cromosoma dal miglior risultato globale per salvarlo
        Chromosome global_champion(lunghezza_cromosoma, input_file);
        global_champion.setFromVector(all_best_chromosomes[best_global_rank]);
        
        // Salva il miglior percorso globale
        string global_champion_file = prefix + "_GLOBAL_CHAMPION_UNO.dat";
        global_champion.printPath(global_champion_file);
        
        // Stampa riassunto finale
        cout << "\n========================================" << endl;
        cout << "=== RISULTATI FINALI GLOBALI ===" << endl;
        cout << "========================================" << endl;
        
        for (int k = 0; k < size; k++) {
            cout << "Processo " << k << " - Miglior fitness: " << all_best_fitness[k];
            if (k == best_global_rank) cout << " <- CAMPIONE GLOBALE";
            cout << endl;
        }
        
        cout << "\nCAMPIONE GLOBALE:" << endl;
        global_champion.printString();
        cout << "Miglior fitness globale: " << best_global_fitness << endl;
        cout << "Trovato dal processo: " << best_global_rank << endl;
        
        if (global_champion.check()) {
            cout << "CAMPIONE GLOBALE VALIDO" << endl;
        }
        
        cout << "========================================" << endl;
        
    } else {
        // Gli altri processi inviano i loro risultati al rank 0
        MPI_Send(best_chromosome_out.data(), lunghezza_cromosoma, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&fitness_out, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}