#include "population.h"

// Un costruttore di una popolazione costruisce la generazione 0 a partire da un cromosoma mutandolo
Population :: Population(int n, Chromosome c) {
    _number_chromosomes = n;
    _chromosomes.set_size(_number_chromosomes);

    // Inizializzo il generatore random con la funzione che ho fatto
    rnd_pop.Start(3);

    // Costruisco la matrice delle distanze a partire dal primo cromosoma
    double number_genes = c.get_genes().size();
    mat matrix(number_genes, number_genes);
    for (int i = 0; i < number_genes; i++) {
        for (int j = 0; j < number_genes; j++) {
            matrix(i, j) = pow(c.get_genes()(i).x - c.get_genes()(j).x, 2) + pow(c.get_genes()(i).y - c.get_genes()(j).y, 2);
        }
    }

    lenght_matrix = matrix;

    int a = 0;
    int m = 0;
    // Muto il cromosoma c e creo gli altri cromosomi per la generazione0, ovviamente lo muto con probabilità 1
    // Guinea Pig è lì così che tutte le modifiche hanno come radice c, è sempre lui in modo diretto il paziente 0 da cui nasce la prima generazione
    for (int i = 0; i < _number_chromosomes; i++) {
        Chromosome guinea_pig = c;
        // guinea_pig.printString();
        // Lo modifico sempre 100 volte, si può chiaramente cambiare. Penso che molte meno bastino comunque
        for (int j = 0; j < 100; j++) {
            if (guinea_pig.check()) a++;
                else {m++; cout << "ERRORE";};
            guinea_pig.mutate(rnd_pop); // mutazione avviene al 100% ogni volta, quale/i avvenga/avvengano è sempre diverso
        }
        // guinea_pig.printString();
        _chromosomes(i) = guinea_pig;
    }

    // cerr << "Errori: " << m << " Giusti: " << a << endl;
    // Ordino la generazione 0
    sort();

    // Setto dei valori inzili che andò poi a cambiare dinamicamente a seconda della convergenza

    _prob_mutation = 0.35;
    _prob_cross = 0.65;
    _prob_elitismo = 0.2;


    _prob_select_cross = 2.5;
    _prob_select_mutation = 3;
    _prob_select_elitismo = 1.5;
}

// COSÌ HO CONTROLLO SIA SULLA PROBABILITÀ DI INCROCIARE/MUTARE SIA SU CHI INCROCIARE E CHI MUTARE
// DI BASE VORREI INCROCIARE CROMOSOMI "FORTI" E MUTARE QUELLI "DEBOLI"
void Population :: evolve(double cross, double mutation, double pass, double select_cross, double select_mutation) {
    // Evolve al generazione sia sendo metodi di cross che secondo mutazioni. Favorire il cross
    // CREO UNA GENERAZIONE MOMENTANEA CHE POI SOSTITUIRÀ LA CORRENTE
    // Non sto a creare proprio una classe Popolazione, solo il filed di cromosomi
    arma::field<Chromosome> next_gen(_number_chromosomes);
    int index = 0;
    //if (averageFitnessBestHalf.second / averageFitnessBestHalf.first < 1e-10) {mutation += 0.075; cross -= 0.075;}

    while (index < _number_chromosomes) {

        // SE SI DECIDE DI ACCOPPIARE ALLORA PREDO DUE POSTI IN next_gen() figlio e figlia, index +=2
        double rnd_cross = rnd_pop.Rannyu();
        // SE SI DECIDE DI MUTARE ALLORA PREDO UN SOLO POSTO IN next_gen() genitore sigle mutato
        double rnd_mutation = rnd_pop.Rannyu();
        double rnd_pass = rnd_pop.Rannyu();

        if (rnd_cross < cross && index < _number_chromosomes - 1) {
            Chromosome mother = select(select_cross);
            Chromosome father = select(select_cross);
            next_gen(index) = mother;
            next_gen(index + 1) = father;
            crossover_cut(mother, father, next_gen(index), next_gen(index + 1));
            index += 2;
        }if (rnd_mutation < mutation && index < _number_chromosomes) {
            Chromosome guinea_pig = select(select_mutation);
            guinea_pig.mutate(rnd_pop);
            next_gen(index) = guinea_pig;
            index++;
        }
        if (rnd_pass < pass && index < _number_chromosomes) {
            // fallback: copia un cromosoma esistente senza mutazioni, ELITISMO
            next_gen(index) = select(1.5); // è uno forte che sopravvive all'ambiente e passa alla generazione successiva
            index++;
        }
        // cerr << "index: " << index << endl;
        //cerr << "numero: " << _number_chromosomes << endl;
    }
    _chromosomes = next_gen;
}

void Population :: sort() {
    double min;
    int min_index;

    for (int i = 0; i < _number_chromosomes; i++) {
        min = _chromosomes(i).evaluateFit();
        min_index = i;

        for (int j = i + 1; j < _number_chromosomes; j++) {
            if (_chromosomes(j).evaluateFit() < min) {
                min = _chromosomes(j).evaluateFit();
                min_index = j;
            }
        }
        // Swap manuale
        Chromosome temp = _chromosomes(i);
        _chromosomes(i) = _chromosomes(min_index);
        _chromosomes(min_index) = temp;

        //swap(_chromosomes(i), _chromosomes(min_index));
    }

    for (int i = 0; i < _number_chromosomes; i++) {
        if (!_chromosomes(i).check()) cerr << "ERRORE" << endl;
    }
}

// COn p > 1 pesco più in cima, ho j più piccoli, prendo i più forti
Chromosome Population :: select(double p) {
    int j = static_cast<int>((_number_chromosomes - 1.) * pow(rnd_pop.Rannyu(), p)) + 1;
    return _chromosomes(j);
}

void Population :: crossover_cut(Chromosome mother, Chromosome father, Chromosome &son, Chromosome &daughter) {

    int N = father.get_genes().size();
    int cut_point = static_cast<int>(rnd_pop.Rannyu(1, N - 1));

    for (uword i = 0; i < cut_point; ++i) {
        son.get_genes()(i) = father.get_genes()(i);
        daughter.get_genes()(i) = mother.get_genes()(i);
    }

    int index = 0;
    for (int k = 0; k < N; k++) {
        bool already_here = false;

        for (int j = 0; j < cut_point; j++) {
            if (mother.get_genes()(k).index == son.get_genes()(j).index) already_here = true;
        }

        if (!already_here) {
            son.get_genes()(cut_point + index) = mother.get_genes()(k);
            index++;
        }
    }

    index = 0;
    for (int k = 0; k < N; k++) {
        bool already_here = false;

        for (int j = 0; j < cut_point; j++) {
            if (father.get_genes()(k).index == daughter.get_genes()(j).index) already_here = true;
        }

        if (!already_here) {
            daughter.get_genes()(cut_point + index) = father.get_genes()(k);
            index++;
        }
    }

}

double Population :: getBestFit() {
    return _chromosomes(0).evaluateFit();
}

Chromosome Population :: getBest() {
    return _chromosomes(0);
}

pair <double, double> Population :: getAverageFit() {
    double somma = 0.0;
    double varianza = 0.0;
    int size = _chromosomes.size();

    for (int i = 0; i < size; i++) {
        double fit = _chromosomes(i).evaluateFit();
        somma += fit;
    }

    double media = somma / size;

    for (int i = 0; i < size; i++) {
        varianza += (_chromosomes(i).evaluateFit() - media) * (_chromosomes(i).evaluateFit() - media);
    }

    varianza /= size;

    double err = sqrt(varianza); 

    averageFitness.first = media;
    averageFitness.second = err;

    return make_pair(media, err);
}

pair <double, double> Population :: getAverageFitBestHalf() {
    double somma = 0.0;
    double somma_2 = 0.0;
    double varianza = 0.0;
    int size = _chromosomes.size() / 2;

    for (int i = 0; i < size; i++) {
        double fit = _chromosomes(i).evaluateFit();
        somma += fit;
        somma_2 += fit * fit;
    }

    double media = somma / size;

    for (int i = 0; i < size; i++) {
        varianza += (_chromosomes(i).evaluateFit() - media) * (_chromosomes(i).evaluateFit() - media);
    }

    varianza /= size;

    double err = sqrt(varianza); 

    averageFitnessBestHalf.first = media;
    averageFitnessBestHalf.second = err;

    return make_pair(media, err);
}


void Population::evolve() {
    // Variabili statiche per tracking storico
    static double prev_best = std::numeric_limits<double>::max();
    static int stagnation_count = 0;
    static int generation = 0;
    
    // Calcolo metriche attuali
    auto fit = getAverageFitBestHalf();
    double avg = fit.first;
    double err = fit.second;
    double best_current = getBestFit();
    double diversity_factor = (avg > 0) ? err / avg : 0.5; // Protezione divisione per zero
    
    // Controllo miglioramento
    double improvement = (prev_best - best_current) / std::max(prev_best, 1.0);
    if (improvement > 0.001) {
        stagnation_count = 0; // Reset se c'è miglioramento
    } else {
        stagnation_count++;
    }
    
    // Adattamento graduale e bilanciato
    double lr = 0.02; // Learning rate molto più conservativo
    
    // Strategia principale: bilanciare diversità e sfruttamento
    if (diversity_factor < 0.1 || stagnation_count > 50) {
        // Popolazione troppo convergente o stagnante -> più esplorazione
        _prob_mutation = std::min(_prob_mutation + lr * 2, 0.6);
        _prob_cross = std::max(_prob_cross - lr, 0.3);
        _prob_elitismo = std::max(_prob_elitismo - lr, 0.05);
    } else if (diversity_factor > 0.4) {
        // Troppa diversità -> più sfruttamento
        _prob_mutation = std::max(_prob_mutation - lr, 0.1);
        _prob_cross = std::min(_prob_cross + lr, 0.7);
        _prob_elitismo = std::min(_prob_elitismo + lr * 0.5, 0.3);
    }
    
    // Normalizzazione per evitare che la somma sia troppo alta
    double total = _prob_mutation + _prob_cross + _prob_elitismo;
    if (total > 1.2) {
        double factor = 1.0 / total;
        _prob_mutation *= factor;
        _prob_cross *= factor;
        _prob_elitismo *= factor;
    }
    
    // Debug periodico
    if (generation % 200 == 0) {
        std::cout << "Gen " << generation << ": diversity=" << diversity_factor 
                  << ", stagnation=" << stagnation_count 
                  << ", probs=[" << _prob_mutation << "," << _prob_cross << "," << _prob_elitismo << "]" << std::endl;
    }
    
    // Generazione nuova popolazione - METODO SEQUENZIALE GARANTITO
    arma::field<Chromosome> next_gen(_number_chromosomes);
    int filled = 0;
    
    // Prima: garantisci alcuni crossover (sono importanti)
    int target_crossovers = static_cast<int>(_prob_cross * _number_chromosomes / 2) * 2; // Pari
    for (int i = 0; i < target_crossovers && filled < _number_chromosomes - 1; i += 2) {
        Chromosome mother = select(_prob_select_cross);
        Chromosome father = select(_prob_select_cross);
        next_gen(filled) = mother;
        next_gen(filled + 1) = father;
        crossover_cut(mother, father, next_gen(filled), next_gen(filled + 1));
        filled += 2;
    }
    
    // Secondo: mutazioni
    int target_mutations = static_cast<int>(_prob_mutation * _number_chromosomes);
    for (int i = 0; i < target_mutations && filled < _number_chromosomes; i++) {
        Chromosome mutant = select(_prob_select_mutation);
        mutant.mutate(rnd_pop);
        next_gen(filled) = mutant;
        filled++;
    }
    
    // Terzo: riempi il resto con elitismo
    while (filled < _number_chromosomes) {
        next_gen(filled) = select(_prob_select_elitismo);
        filled++;
    }
    
    // Sostituisco la popolazione
    _chromosomes = next_gen;
    prev_best = best_current;
    generation++;
}
// VEDERE QUANOT MANCA AL VALORE PER RAGGIUNGERE IL SUO MASSIMO CHE NON VOLGIO SUPERI E DIVIDERE FRATTO 2
// COSÌ NON LO RAGGIUNGERÀ MAI MA CI SI AVVICINA SEMPRE

/*
void Population::evolve() {
    // Calcolo media e deviazione della metà migliore
    auto fit = getAverageFitBestHalf();
    double avg = fit.first;
    double err = fit.second;
    double factor = err / avg; // basso -> popolazione simile, alto -> diversità
    
    // Parametri di adattamento
    double learning_rate = 0.5; // quanto velocemente si adatta
    double diversity_threshold = 0.2; // soglia per considerare "bassa diversità"
    
    // Adattamento reattivo basato sulle condizioni attuali
    if (factor < diversity_threshold) {
        // Popolazione troppo convergente -> aumenta esplorazione
        _prob_mutation += learning_rate * (diversity_threshold - factor);
        _prob_cross -= learning_rate * 0.5 * (diversity_threshold - factor);
        _prob_elitismo -= learning_rate * (diversity_threshold - factor);
    } else {
        // Alta diversità -> stabilizza e sfrutta
        _prob_mutation -= learning_rate * (factor - diversity_threshold) * 0.5;
        _prob_cross += learning_rate * (factor - diversity_threshold);
        _prob_elitismo += learning_rate * (factor - diversity_threshold) * 0.3;
    }
    
    // Adattamento delle probabilità di selezione basato sulla performance
    static double prev_avg_fitness = avg; // mantiene memoria della generazione precedente
    double fitness_improvement = (avg - prev_avg_fitness) / std::max(prev_avg_fitness, 1e-10);
    
    if (fitness_improvement > 0.01) {
        // Buon miglioramento -> rinforza strategie attuali
        _prob_select_cross *= (1.0 + learning_rate * fitness_improvement);
        _prob_select_mutation *= (1.0 + learning_rate * fitness_improvement * 0.5);
    } else if (fitness_improvement < -0.01) {
        // Peggioramento -> aumenta esplorazione
        _prob_select_mutation *= (1.0 + learning_rate * std::abs(fitness_improvement));
        _prob_select_cross *= (1.0 - learning_rate * std::abs(fitness_improvement) * 0.3);
    }
    
    prev_avg_fitness = avg;
    
    // Clamp dei valori per restare in intervalli sensati
    _prob_mutation = std::min(0.8, std::max(0.05, _prob_mutation));
    _prob_cross = std::min(0.9, std::max(0.1, _prob_cross));
    _prob_elitismo = std::min(0.5, std::max(0.05, _prob_elitismo));
    
    _prob_select_cross = std::min(3.0, std::max(0.5, _prob_select_cross));
    _prob_select_mutation = std::min(3.0, std::max(0.5, _prob_select_mutation));
    _prob_select_elitismo = 1.5; // mantieni fisso per stabilità
    
    // Creazione della nuova generazione (resto del codice invariato)
    arma::field<Chromosome> next_gen(_number_chromosomes);
    int index = 0;
    
    while (index < _number_chromosomes) {
        double rnd_cross = rnd_pop.Rannyu();
        double rnd_mutation = rnd_pop.Rannyu();
        double rnd_pass = rnd_pop.Rannyu();
        
        // Crossover tra cromosomi forti
        if (rnd_cross < _prob_cross && index < _number_chromosomes - 1) {
            Chromosome mother = select(_prob_select_cross);
            Chromosome father = select(_prob_select_cross);
            next_gen(index) = mother;
            next_gen(index + 1) = father;
            crossover_cut(mother, father, next_gen(index), next_gen(index + 1));
            index += 2;
        }
        
        // Mutazione di cromosomi deboli
        if (rnd_mutation < _prob_mutation && index < _number_chromosomes) {
            Chromosome guinea_pig = select(_prob_select_mutation);
            guinea_pig.mutate(rnd_pop);
            next_gen(index) = guinea_pig;
            index++;
        }
        
        // Elitismo: copia diretta dei migliori
        if (rnd_pass < _prob_elitismo && index < _number_chromosomes) {
            next_gen(index) = select(_prob_select_elitismo);
            index++;
        }
    }
    
    // Sostituisco la popolazione corrente con la nuova
    _chromosomes = next_gen;
}
*/