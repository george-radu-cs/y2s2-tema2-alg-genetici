#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <vector>
#include <string>
#include <cmath>

class Chromosome {
private:
    int number_of_genes; // chromosome number_of_genes
    // from cpp reference page the compiler will optimize the space for this vector, using only one bit per bool
    std::vector<bool> genes;

public:
    Chromosome(const int number_of_genes);

    Chromosome(const std::vector<bool> &genes);

    Chromosome(const Chromosome &chromosome);

    int getNumberOfGenes() const;

    const std::vector<bool> &getGenome() const;

    void setGenome(const std::vector<bool> &genes);

    void generateRandomGenome();

    std::string getStringRepresentationOfGenome() const;

    void setNumberOfGenes(int number_of_genes);
};


#endif // CHROMOSOME_H