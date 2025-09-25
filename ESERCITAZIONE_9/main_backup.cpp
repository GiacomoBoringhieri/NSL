#include "random.h"
#include "chromosome.h"
#include "population.h"

#include <iostream>
#include <fstream>

using namespace std;

int main (int argc, char *argv[]){

  // Lunghezza del cromosoma (numero di città)
  int lunghezza_cromosoma = 34;

  // Dimensione della popolazione
  int lunghezza_popolazione = 500;

  int numero_generazioni = 5000;

    // Imposto se usare auto-adattamento
  bool auto_adattamento = true;

  // Selezione del file di input iniziale: "box.out" o "circle.out"
  string input_file = "box.out";   // default
  if (argc > 1) input_file = argv[1]; // possibilità di passare come argomento



  // Prefisso per i file di output in base al file iniziale
  string prefix = (input_file == "box.out") ? "box" : "circle";

  bool prov_ita = false;

  if (prov_ita){
    input_file = "cap_prov_ita.dat";
    prefix = "cap_prov";
  }

  cout << "USO DEL FILE INIZIALE: " << input_file << endl;
  

  // Inizializzo il generatore random
  Random rnd;
  rnd.Start(3);

  // Genero le città su cerchio e su box
  rnd.city_circle(lunghezza_cromosoma);
  rnd.city_box(lunghezza_cromosoma);

  // Creo il cromosoma iniziale Adam usando il file scelto
  Chromosome Adam(lunghezza_cromosoma, input_file);

  // Stampo il cromosoma iniziale
  cout << "CROMOSOMA INIZIALE: \n";
  Adam.printString();

  // Creo la popolazione iniziale
  Population generation_0(lunghezza_popolazione, Adam);

  // Ordino la popolazione iniziale per fitness
  generation_0.sort();

  // Stampo il miglior fit della generazione 0
  cout << "Migliore fit della generazione 0: " << generation_0.getBestFit() << endl;

  // Salvo il percorso iniziale del miglior cromosoma
  generation_0.getBest().printPath(prefix + "_best_path_0.dat");

  // Definisco i file di output finali
  string best_half_file  = prefix + (auto_adattamento ? "_best_half_auto.dat" : "_best_half.dat");
  string champion_path_file   = prefix + (auto_adattamento ? "_best_path_final_auto.dat" : "_best_path_final.dat");
  string champions_file = prefix + (auto_adattamento ? "_champions_auto.dat" : "_champions.dat");

  if (auto_adattamento)
    cout << "SIMULAZIONE CON AUTO-ADATTAMENTO DELLA POPOLAZIONE" << endl;

  // Apro il file per salvare l'andamento della media della metà migliore
  ofstream off_file(best_half_file);
  ofstream off_champion_file(champions_file);

  if (off_file.is_open()) {
    // Evoluzione per generazioni
    for (int i = 0; i < numero_generazioni; i++) {
      if (auto_adattamento)
        generation_0.evolve();  // evolve auto-adattativo
      else
        generation_0.evolve(0.65, 0.25, 0.15, 2, 3); // evolve manuale

      // Ordino la popolazione
      generation_0.sort();

      // Salvataggio e stampa ogni 100 generazioni
      if (i % 100 == 0) {
        auto fit_half = generation_0.getAverageFitBestHalf();
        cout << "GEN " << i << " -> MEDIA: " << fit_half.first
             << " PM " << fit_half.second << endl;
        off_file << i << " " << fit_half.first << " " << fit_half.second << "\n";
        off_champion_file << i << " " << generation_0.getBestFit() << "\n";
      }
    }
  } else {
    cerr << "PROBLEM: Unable to open file " << best_half_file << endl;
  }

  off_file.close();
  off_champion_file.close();

  generation_0.printProbabilities();

  // Salvo il miglior percorso finale
  generation_0.getBest().printPath(champion_path_file);

  // Stampo informazioni sul cromosoma migliore
  cout << "CAMPIONE DELLA GENERAZIONE FINALE:\n";
  generation_0.getBest().printString();
  cout << "Migliore fit della generazione finale: " << generation_0.getBestFit() << endl;

  auto fit_final = generation_0.getAverageFitBestHalf();
  cout << "MEDIA: " << fit_final.first << " PM " << fit_final.second << endl;

  if (generation_0.getBest().check())
    cout << "CAMPIONE VALIDO" << endl;

  return 0;
}