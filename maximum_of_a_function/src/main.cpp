#include "GeneticAlgorithm.h"

int main() {
    GeneticAlgorithm ga = GeneticAlgorithm("../data/settings.txt", "../data/evolution.txt", false);
    ga.generatePopulation();
    ga.runGeneticAlghorithm();

    return 0;
}
