#include "Chromosome.h"
#include "RandomGenerator.h"

Chromosome::Chromosome(const int number_of_genes)
        : number_of_genes(number_of_genes) {
    this->generateRandomGenome();
}

Chromosome::Chromosome(const std::vector<bool> &genes)
        : number_of_genes((int) genes.size()) {
    setGenome(genes);
}

Chromosome::Chromosome(const Chromosome &chromosome) {
    this->setNumberOfGenes(chromosome.getNumberOfGenes());
    this->setGenome(chromosome.getGenome());
}

int Chromosome::getNumberOfGenes() const {
    return this->number_of_genes;
}

const std::vector<bool> &Chromosome::getGenome() const {
    return this->genes;
}

void Chromosome::setGenome(const std::vector<bool> &genes) {
    size_t genomeSize = genes.size();
    this->genes.resize(genomeSize);
    for (size_t i = 0; i < genomeSize; ++i) {
        this->genes[i] = genes[i];
    }
}

void Chromosome::generateRandomGenome() {
    // clear the current genome, in case if exists, generate a random genome
    this->genes.clear();
    this->genes.resize(this->number_of_genes);
    for (size_t i = 0; i < this->number_of_genes; ++i) {
        genes[i] = RandomGenerator::getOneOrZero();
    }
}

std::string Chromosome::getStringRepresentationOfGenome() const {
    std::string genomeAsString;
    for (const auto &gene: this->genes) {
        genomeAsString += gene ? "1" : "0";
    }
    return genomeAsString;
}

void Chromosome::setNumberOfGenes(int number_of_genes) {
    this->number_of_genes = number_of_genes;
}
