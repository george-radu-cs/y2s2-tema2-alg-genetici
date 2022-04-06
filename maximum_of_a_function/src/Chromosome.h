#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <vector>
#include <string>
#include <cmath>

class Chromosome {
private:
    int number_of_genes; // chromosome number_of_genes
    std::pair<long double, long double> domain;
    long double precision;
    // from cpp reference page the compiler will optimize the space for this vector, using only one bit per bool
    std::vector<bool> genes;
    long double value;

public:
    Chromosome(const std::pair<long double, long double> &domain, long double precision);

    Chromosome(const std::pair<long double, long double> &domain, long double precision,
               const std::vector<bool> &genes);

    Chromosome(const Chromosome&chromosome);

    void calculateNumberOfGenes();

    int getNumberOfGenes() const;

    const std::vector<bool> &getGenome() const;

    void setGenome(const std::vector<bool> &genes);

    std::pair<long double, long double> getDomain() const;

    void setDomain(const std::pair<long double, long double> &domain);

    long double getPrecision() const;

    void setPrecision(long double precision);

    void generateRandomGenome();

    std::string getStringRepresentationOfGenome() const;

    // used to compute calculate the value of the chromosome and save it's value for later usages
    // O(m) where m = len(genome)
    void computeValue();

    void setValue(long double value);

    long double getValue() const;

    void setNumberOfGenes(int number_of_genes);
};


#endif // CHROMOSOME_H