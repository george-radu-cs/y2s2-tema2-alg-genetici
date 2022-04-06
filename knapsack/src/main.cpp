#include "GeneticAlgorithm.h"

int main() {
    GeneticAlgorithm ga = GeneticAlgorithm("../data/settings.txt", "../data/evolution.txt", true);
    ga.generatePopulation();
    ga.generateObjects();
    ga.print_objects();
    ga.runGeneticAlghorithm();

    return 0;
}
