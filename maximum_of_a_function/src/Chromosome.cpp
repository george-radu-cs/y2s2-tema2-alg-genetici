#include "Chromosome.h"
#include "RandomGenerator.h"

Chromosome::Chromosome(const std::pair<long double, long double> &domain, long double precision)
        : domain(domain), precision(precision) {
    this->calculateNumberOfGenes();
    this->generateRandomGenome();
    this->computeValue();
}

Chromosome::Chromosome(const std::pair<long double, long double> &domain, long double precision,
                       const std::vector<bool> &genes)
        : domain(domain), precision(precision), number_of_genes((int) genes.size()) {
    setGenome(genes);
}

Chromosome::Chromosome(const Chromosome &chromosome) {
    this->setDomain(chromosome.getDomain());
    this->setPrecision(chromosome.getPrecision());
    this->setNumberOfGenes(chromosome.getNumberOfGenes());
    this->setGenome(chromosome.getGenome());
}

void Chromosome::calculateNumberOfGenes() {
    // meshing the interval [a,b] in (b-a) * 10**precision subintervals (elements)
    long double number_of_genes = log((this->domain.second - this->domain.first) * pow(10, this->precision));
    this->number_of_genes = ceil(abs(number_of_genes));
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
    this->computeValue();
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

void Chromosome::computeValue() {
    // linear translation from the binary codification of the genome into a number from the domain
    // gene_value will be initially the value of the least significant gene of the chromosome
    long double gene_value = (this->domain.second - this->domain.first) / (pow(2, this->number_of_genes) - 1);
    // the value of the chromosome will be the value `a` from the domain plus the number codified in the genome
    long double value{this->domain.first};
    for (int i = this->number_of_genes - 1; i >= 0; --i) {
        if (this->genes[i]) { // the gene at position i is marked / the bit is flagged
            value += gene_value;
        }
        gene_value *= 2; // compute the value for the next gene based on binary codification
    }
    // try to round to the precision required despite floating point representation
    long double value_rounded_to_precision = (value * pow(10, this->precision)) / pow(10, this->precision);
    this->setValue(value_rounded_to_precision);
}

void Chromosome::setValue(long double value) {
    this->value = value;
}

long double Chromosome::getValue() const {
    return this->value;
}

long double Chromosome::getPrecision() const {
    return this->precision;
}

void Chromosome::setPrecision(long double precision) {
    this->precision = precision;
}

std::pair<long double, long double> Chromosome::getDomain() const {
    return this->domain;
}

void Chromosome::setDomain(const std::pair<long double, long double> &domain) {
    this->domain.first = domain.first;
    this->domain.second = domain.second;
}

void Chromosome::setNumberOfGenes(int number_of_genes) {
    this->number_of_genes = number_of_genes;
}
