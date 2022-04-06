#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <iomanip>
#include <fstream>
#include <string>
#include <limits>
#include <sstream>
#include <iostream>
#include "Chromosome.h"

#define OUTPUT_PRECISION 20

class GeneticAlgorithm {
private:
    // settings for the algorithm
    int population_size{};
    std::pair<long double, long double> domain_of_definition;
    std::tuple<long double, long double, long double> function_params;
    int precision{}; // used to discretize the interval
    long double recombination_probability{}; // crossover
    long double mutation_probability{};
    int steps{};
    bool keep_elitist{};
    bool show_steps{false};
    std::ofstream file_output;
    int points_of_crossing{1};

    std::vector<Chromosome> population;
    std::vector<long double> selection_probability;
    std::vector<long double> selection_probability_interval;
    std::vector<Chromosome> selected_population;
    std::vector<size_t> population_selected_for_crossing;

    GeneticAlgorithm(const GeneticAlgorithm &);

    GeneticAlgorithm &operator=(const GeneticAlgorithm &);

    void setShowSteps(bool show_steps);

    void computeNextGeneration();

    long double fitnessFunction(const Chromosome &chromosome) const;

    long double findMaximumValueOfFitnessFunction() const;

    size_t findElitistChromosomeIndexFromPopulation() const;

    long double
    computeSelectionProbabilityForChromosome(const Chromosome &chromosome, const long double &sumValueFunction) const;

    long double computeSumFunctionForPopulation() const;

    void computeSelectionProbabilitiesForPopulation();

    void computeSelectionProbabilityInterval();

    size_t findChromosomeIndexFromSelectionProbabilityInterval(long double value) const;

    size_t findChromosomeIndexFromSelectionProbabilityIntervalBinarySearch(long double value, size_t left, size_t right,
                                                                           size_t length) const;

    void selectNextPopulation();

    void computeCrossingParticipationFromPopulation();

    static void
    crossGenomesOfTwoChromosomes(Chromosome &chromosome1, Chromosome &chromosome2, const size_t point_of_crossing);

    static void
    crossGenomesOfTwoChromosomesWithTwoPointsOfCrossingSwitchMiddle(Chromosome &chromosome1, Chromosome &chromosome2,
                                                                    const size_t point_of_crossing_left,
                                                                    const size_t point_of_crossing_right);

    void crossPopulation();

    void mutatePopulation();

    void saveNextGeneration();

    void printInitialPopulation() const;

    void printSelectionProbability() const;

    void printSelectionProbabilityInterval() const;

    void printSelectedPopulationInfo() const;

    void printPopulationInfoAfterCrossing() const;

    void printPopulationInfoAfterMutation() const;

public:
    GeneticAlgorithm(const std::string &settings_filename, std::string const &output_filename,
                     bool keep_elitist = false);

    ~GeneticAlgorithm();

    void generatePopulation();

    void runGeneticAlghorithm();
};


#endif // GENETIC_ALGORITHM_H