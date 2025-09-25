#include "chromosome.h"

Chromosome :: Chromosome(int n, string filename) {

    _genes_number = n;
    _genes.set_size(_genes_number);

    ifstream infile(filename);
    if (!infile.is_open()) {
        throw invalid_argument("Errore apertura file input");
    }
    int i = 0;
    double val1, val2;
    while (infile >> val1 >> val2) {
        _genes(i).index = i;
        _genes(i).x = val1;
        _genes(i).y = val2;
        i++;
    }

    infile.close();

}

// DA MIGLIORARE, NON È PER NULLA EFFICACE
bool Chromosome :: check() {
    bool number = true;
    int n = 0;
    for (int i = 0; i < _genes.size(); i++) {
        n = 0;
        for (int j = 0; j < _genes.size(); j++) {
            if (_genes(j).index == i) n++;
        }
        if (n == 0 || n > 1) {
                number = false; 
                cerr << "ERROR: Indice non presente o ripetuto" << std::endl;
            }
    }
    return number;
}

// FATTA DA AI
/*
bool Chromosome :: check() {
    // Crea un vettore per tracciare gli indici trovati
    std::vector<bool> found(_genes_number, false);
    bool valid = true;
    
    for (int i = 0; i < _genes.size(); i++) {
        int current_index = _genes(i).index;
        
        // Controlla se l'indice è nel range valido
        if (current_index < 0 || current_index >= _genes_number) {
            cerr << "ERROR: Indice " << current_index << " fuori range" << endl;
            valid = false;
            continue;
        }
        
        // Controlla se l'indice è già stato trovato
        if (found[current_index]) {
            cerr << "ERROR: Indice " << current_index << " duplicato" << endl;
            valid = false;
        }
        
        found[current_index] = true;
    }
    
    // Controlla se ci sono indici mancanti
    for (int i = 0; i < _genes_number; i++) {
        if (!found[i]) {
            cerr << "ERROR: Indice " << i << " mancante" << endl;
            valid = false;
        }
    }
    
    return valid;
}
*/

// Calcola la lunghezza totale del percorso rappresentato dal cromosoma
double Chromosome::evaluateFit() {
    bool L_1 = true;
    double sum = 0.0;
    // Somma delle distanze euclidee tra ogni coppia di geni consecutivi
    for (int i = 0; i < _genes.size() - 1; i++) {
        if (L_1) sum += abs((_genes(i).x - _genes(i+1).x)) + abs((_genes(i).y - _genes(i+1).y));
        else sum += pow((_genes(i).x - _genes(i+1).x), 2) + pow((_genes(i).y - _genes(i+1).y), 2);
        
    }
    // Aggiunge la distanza tra l'ultimo e il primo gene per chiudere il ciclo
    if (L_1) sum += abs((_genes(_genes.size()-1).x - _genes(0).x)) + abs((_genes(_genes.size()-1).y - _genes(0).y));
    else sum += pow((_genes(_genes.size()-1).x - _genes(0).x), 2) + pow((_genes(_genes.size()-1).y - _genes(0).y), 2);
    
    return sum;
}

void Chromosome::printPath(string filename){
    ofstream off_file;
    off_file.open(filename);
    if (off_file.is_open()){

        for (int i = 0; i < _genes_number; i++) {
            off_file << _genes(i).x << " " << _genes(i).y << endl;
        }
        off_file << _genes(0).x << " " << _genes(0).y << endl;

    } else cerr << "PROBLEM: Unable to open path file" << endl;
    off_file.close();
    return;
}

void Chromosome::printString() {
    // cerr << "Numero di geni: " << _genes_number << endl;
    for (int i = 0; i < _genes_number; i++) {
            cerr << _genes(i).index << " ";
        }
    cerr << endl;
}

// Scambia due geni in posizioni casuali, escludendo la prima posizione (0)
void Chromosome::mutate_Permutation(Random &rnd) {
    int a = static_cast<int>(rnd.Rannyu(1, _genes_number)); // a ∈ [1, _genes_number - 1]
    int b = static_cast<int>(rnd.Rannyu(1, _genes_number)); // b ∈ [1, _genes_number - 1]
    std::swap(_genes[a], _genes[b]);
}

void Chromosome::mutate_Shift(Random &rnd) {
    uword size = _genes_number;

    if (size < 3) { // sicurezza: serve almeno 3 geni per fare uno shift
        cerr << "Cromosoma troppo piccolo per shift\n";
        return;
    }

    // Scegli posizione di inizio del blocco (dopo la prima città)
    uword start = static_cast<uword>(rnd.Rannyu(1, size - 1)); // [1, size-2]

    // Lunghezza del blocco m (assicurata dentro il cromosoma)
    uword max_block = size - start;
    if (max_block == 0) return; // sicurezza
    uword m = static_cast<uword>(rnd.Rannyu(1, max_block)); // [1, max_block]

    cerr << "Shift: start=" << start << ", m=" << m << endl;

    // Shift massimo possibile senza superare la fine del cromosoma
    uword max_shift = size - start - m;
    if (max_shift == 0) return; // nessun posto per shift
    uword n = static_cast<uword>(rnd.Rannyu(1, max_shift + 1)); // [1, max_shift]

    cerr << "Shift: n=" << n << endl;

    // Copia i geni in un nuovo array temporaneo
    arma::field<Gene> new_genes(size);

    // Copia i geni prima del blocco invariati
    for (uword i = 0; i < start; ++i)
        new_genes(i) = _genes(i);

    // Copia i geni dopo il blocco inizialmente (dopo start + m)
    for (uword i = start + m; i < size; ++i)
        new_genes(i) = _genes(i);

    // Sposta il blocco all’interno del cromosoma
    for (uword i = 0; i < m; ++i)
        new_genes[start + n + i] = _genes[start + i];

    _genes = new_genes;
}

void Chromosome::mutate_Blocks(Random &rnd) {
    uword size = _genes_number;
    //cerr << "NUMERO DI GENI: " << size << endl;

    // Limitiamo m a < N/4 per sicurezza
    uword m = static_cast<uword>(rnd.Rannyu(1, size / 4.));
    //cerr << "LUNGHEZZA BLOCCO: " << m << endl;

    // Primo blocco (escludendo la prima città)
    uword start1 = static_cast<uword>(rnd.Rannyu(1, size - m));

    // Costruisco un vettore di posizioni valide per start2
    std::vector<uword> valid_positions;
    for (uword s = 1; s <= size - m; ++s) {
        if (s + m <= start1 || s >= start1 + m) { // non sovrappone a start1
            valid_positions.push_back(s);
        }
    }

    // Se non ci sono posizioni valide, esco
    if (valid_positions.empty()) {
        cerr << "Impossibile trovare un secondo blocco non sovrapposto\n";
        return;
    }

    // Scelgo start2 casualmente tra le posizioni valide
    uword idx = static_cast<uword>(rnd.Rannyu(0, valid_positions.size()));
    uword start2 = valid_positions[idx];

    //cerr << "Primo blocco in: " << start1 << ", Secondo blocco in: " << start2 << endl;

    // Salva i due blocchi
    arma::field<Gene> block1(m);
    arma::field<Gene> block2(m);
    for (uword i = 0; i < m; ++i) {
        block1(i) = _genes(start1 + i);
        block2(i) = _genes(start2 + i);
    }

    // Scambia i blocchi
    for (uword i = 0; i < m; ++i) {
        _genes(start1 + i) = block2(i);
        _genes(start2 + i) = block1(i);
    }
}

// Inverte un blocco di dimensione m a partire da una posizione casuale, con wrapping
void Chromosome::mutate_Inversion(Random &rnd) {
    uword size = _genes_number - 1; // Esclude il primo gene

    uword m = static_cast<uword>(rnd.Rannyu(1, size));     // Lunghezza della finestra
    uword start = static_cast<uword>(rnd.Rannyu(1, size)); // Posizione iniziale del blocco

    // Copia il blocco da invertire
    field<Gene> temp(m);
    for (uword i = 0; i < m; ++i) {
        temp(i) = _genes(1 + ((start - 1 + i) % size));
    }

    // Scrive il blocco invertito nella posizione originale
    for (uword i = 0; i < m; ++i) {
        _genes(1 + ((start - 1 + i) % size)) = temp(m - 1 - i);
    }
}

void Chromosome::mutate(Random &rnd) {
    const double p_perm = 0.5;
    const double p_shift = 0.0; // disabilitato
    const double p_blocks = 0.5;
    const double p_inversion = 0.5;

    // Decidi quante mutazioni fare (può essere più di una)
    double p = 1.5;
    int num_mutations = static_cast<int>(4. * pow(rnd.Rannyu(), p)) + 1;

    for (int i = 0; i < num_mutations; ++i) {
        if (rnd.Rannyu() < p_perm) {
            mutate_Permutation(rnd);
        }
        if (rnd.Rannyu() < p_shift) {
            // mutate_Shift(rnd);
        }
        if (rnd.Rannyu() < p_blocks) {
            mutate_Blocks(rnd);
        }
        if (rnd.Rannyu() < p_inversion) {
            mutate_Inversion(rnd);
        }
    }
}