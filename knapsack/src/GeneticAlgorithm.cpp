#include "GeneticAlgorithm.h"
#include "RandomGenerator.h"

GeneticAlgorithm::GeneticAlgorithm(const std::string &settings_filename, const std::string &output_filename,
                                   bool keep_elitist) : file_output(output_filename) {
    // read the settings from the file and save them
    std::ifstream fin(settings_filename);
    fin >> this->population_size;

    fin >> this->knapsack_weight;

    fin >> this->number_of_objects; // number of objects
    this->population_size = number_of_objects;

    fin >> this->object_max_value;

    // set the output precision for real numbers
    std::cout << std::setprecision(OUTPUT_PRECISION);

    fin >> this->recombination_probability;

    fin >> this->mutation_probability;

    fin >> this->steps;

    fin >> this->points_of_crossing;

    fin.close(); // close the settings file

    this->keep_elitist = keep_elitist;

    // because of a bug when trying to replace the std::cout with ofstream I decided to redirect the output
    std::cout.rdbuf(file_output.rdbuf()); //redirect std::cout to output_file
}

GeneticAlgorithm::~GeneticAlgorithm() {
    this->file_output.close(); // close the output file
}

void GeneticAlgorithm::setShowSteps(bool show_steps) {
    this->show_steps = show_steps;
}

void GeneticAlgorithm::generateObjects() {
    this->objects_values.clear();
    this->objects_weights.clear();
    int value, weight;
    for (size_t i = 0; i < this->population_size; ++i) {
        value = RandomGenerator::getRandomInt(0, this->object_max_value);
        weight = RandomGenerator::getRandomInt(0, this->knapsack_weight / 10);
        this->objects_values.emplace_back(value);
        this->objects_weights.emplace_back(weight);
    }
}

void GeneticAlgorithm::generatePopulation() {
    // create all the chromosomes
    this->population.clear();
    for (size_t i = 0; i < this->population_size; ++i) {
        this->population.emplace_back(Chromosome(this->number_of_objects));
    }
}

void GeneticAlgorithm::computeNextGeneration() {
    // selection step - choose based on probabilities the next population
    this->computeSelectionProbabilitiesForPopulation();
    this->computeSelectionProbabilityInterval();
    if (show_steps) {
        this->printSelectionProbability();
        this->printSelectionProbabilityInterval();
    }
    this->selectNextPopulation();
    if (show_steps) {
        this->printSelectedPopulationInfo();
    }

    // crossover step - some chromosomes will create children which will replace them
    this->computeCrossingParticipationFromPopulation();
    this->crossPopulation();
    if (show_steps) {
        this->printPopulationInfoAfterCrossing();
    }

    // mutation step - allow mutations of the genes of the population
    this->mutatePopulation();
    if (show_steps) {
        this->printPopulationInfoAfterMutation();
    }

    this->saveNextGeneration();
}

void GeneticAlgorithm::runGeneticAlghorithm() {
    // for the first step show the steps of the algorithm
    this->setShowSteps(true);
    this->printInitialPopulation();
    this->computeNextGeneration(); // compute the first generation

    // for the rest of the steps show only the evolution of the maximum
    this->setShowSteps(false);
    std::cout << "\nEvolution of maximum:\n";
    std::cout << this->findMaximumValueOfFitnessFunction() << "\n";
    for (size_t step = 1; step < this->steps; ++step) {
        this->computeNextGeneration();
        std::cout << this->findMaximumValueOfFitnessFunction() << "\n";

        auto elitist_index = this->findElitistChromosomeIndexFromPopulation();
        std::cout << this->population[elitist_index].getStringRepresentationOfGenome() << '\n';
    }
}

long double GeneticAlgorithm::fitnessFunction(const Chromosome &chromosome) const {
    long double chromosome_value{0}, chromosome_weight{0};
    auto chromosome_genome = chromosome.getGenome();
    for (size_t i = 0; i < this->number_of_objects; ++i) {
        if (chromosome_genome[i]) {
            chromosome_value += this->objects_values[i];
            chromosome_weight += this->objects_weights[i];
        }
    }
    return chromosome_weight <= this->knapsack_weight ? chromosome_value : 0;
}

long double GeneticAlgorithm::findMaximumValueOfFitnessFunction() const {
    long double maximum_value = -1 * std::numeric_limits<long double>::max();
    long double chromosome_value;
    for (const auto &chromsome: this->population) {
        chromosome_value = this->fitnessFunction(chromsome);
        if (maximum_value < chromosome_value) {
            maximum_value = chromosome_value;
        }
    }
    return maximum_value;
}

size_t GeneticAlgorithm::findElitistChromosomeIndexFromPopulation() const {
    size_t elitist_chromosome_index{0};
    long double maximum_value = -1 * std::numeric_limits<long double>::max();
    long double chromosome_value;
    for (size_t i = 0; i < this->population_size; ++i) {
        chromosome_value = this->fitnessFunction(this->population[i]);
        if (maximum_value < chromosome_value) {
            maximum_value = chromosome_value;
            elitist_chromosome_index = i;
        }
    }
    return elitist_chromosome_index;
}

long double GeneticAlgorithm::computeSelectionProbabilityForChromosome(const Chromosome &chromosome,
                                                                       const long double &sumValueFunction) const {
    return this->fitnessFunction(chromosome) / sumValueFunction;
}

long double GeneticAlgorithm::computeSumFunctionForPopulation() const {
    long double sum_function_population{0};
    for (const auto &chromosome: this->population) {
        sum_function_population += this->fitnessFunction(chromosome);
    }
    return sum_function_population;
}

void GeneticAlgorithm::computeSelectionProbabilitiesForPopulation() {
    long double sum_function_population = this->computeSumFunctionForPopulation();
    for (const auto &chromosome: this->population) {
        long double chromosome_selection_probability =
                this->computeSelectionProbabilityForChromosome(chromosome, sum_function_population);
        this->selection_probability.push_back(chromosome_selection_probability);
    }
}

void GeneticAlgorithm::computeSelectionProbabilityInterval() {
    long double value{0};
    for (const auto &sp: this->selection_probability) {
        this->selection_probability_interval.push_back(value);
        value += sp;
    }
    this->selection_probability_interval.push_back(1);
}

size_t GeneticAlgorithm::findChromosomeIndexFromSelectionProbabilityInterval(long double value) const {
    // each chromosome from the population has an interval [a,b) based on it's probability, given a random value chosen
    // if the value belongs to an interval I for the chromosome c then we will select chromosome c
    // the selection is made based of the index of the chromosome from the population vector
    size_t index{0};
    for (size_t i = 1; i <= this->population_size; ++i) {
        if (value < this->selection_probability_interval[i]) {
            index = i - 1;
            break;
        }
    }
    return index;
}

size_t
GeneticAlgorithm::findChromosomeIndexFromSelectionProbabilityIntervalBinarySearch(long double value, size_t left,
                                                                                  size_t right, size_t length) const {
    // each chromosome from the population has an interval [a,b) based on it's probability, given a random value chosen
    // if the value belongs to an interval I for the chromosome c then we will select chromosome c
    // the chromosome c is at index a, if it's interval is [a, b)
    // the selection is made based of the index of the chromosome from the population vector

    // find the middle for the binary search
    size_t middle = (right - left) / 2 + left;

    // if value is on the left side of the array
    if (value < this->selection_probability_interval[middle]) {
        if (middle > 0) { // check for boundaries
            // check if value is in the interval [ value[middle-1] , value[middle] )
            if (this->selection_probability_interval[middle - 1] < value) {
                return middle - 1;
            }
            return this->findChromosomeIndexFromSelectionProbabilityIntervalBinarySearch(value, left, middle, length);
        }
        return middle;
    } else { // value is on the right side of the array
        if (middle + 1 < length) { // check for boundaries
            // check if value is in the interval [ value[middle] , value[middle+1] )
            if (value < this->selection_probability_interval[middle + 1]) {
                return middle;
            }
            return this->findChromosomeIndexFromSelectionProbabilityIntervalBinarySearch(value, middle, right, length);
        }
        return middle;
    }
}

void GeneticAlgorithm::selectNextPopulation() {
    // we will select at each step a chromosome from the current population based on its probability to be selected
    // we will mimic this probability using a random generator from the interval [0,1] and the vector of intervals
    // for each chromosome computed before; at each step we will compute a real number from [0,1] and find its
    // associated chromosome and save it for the next population
    long double generated_value_from_interval;
    size_t selected_chromosome_index;
    // if we choose to keep the elitist chromosome than we will select one less chromosome
    size_t i = this->keep_elitist ? 1 : 0;
    if (keep_elitist) {
        size_t elitist_chromosome_index = this->findElitistChromosomeIndexFromPopulation();
        this->selected_population.emplace_back(this->population[elitist_chromosome_index]);
    }
    for (; i < this->population_size; ++i) {
        generated_value_from_interval = RandomGenerator::getRandomNumberBetweenZeroAndOne();
        selected_chromosome_index = this->findChromosomeIndexFromSelectionProbabilityIntervalBinarySearch(
                generated_value_from_interval, 0, this->selection_probability_interval.size(),
                this->selection_probability_interval.size());
        this->selected_population.emplace_back(this->population[selected_chromosome_index]);
        if (this->show_steps) {
            std::cout << "\nu=" << generated_value_from_interval << " choose chromosome "
                      << std::setw(ceil(log10(this->population_size))) << selected_chromosome_index + 1;
        }
    }
}

void GeneticAlgorithm::computeCrossingParticipationFromPopulation() {
    if (show_steps) {
        std::cout << "\nCrossing probability: " << this->recombination_probability << "\n";
    }
    // from the new population, some chromosomes will participate in crossover, we will decide which ones based on
    // a probability already decided in the settings config
    // to model this probability, if the real number generated appears in the interval [0, recombination_probability]
    // then the chromosome will appear in the crossover step, otherwise will go directly to the next step
    long double u;
    bool is_selected;
    // if we choose to keep the elitist chromosome than we will select one less chromosome
    size_t i = this->keep_elitist ? 1 : 0;
    for (; i < this->population_size; ++i) {
        u = RandomGenerator::getRandomNumberBetweenZeroAndOne();
        is_selected = u <= this->recombination_probability;
        if (is_selected) {
            this->population_selected_for_crossing.push_back(i);
        }
        if (show_steps) {
            std::cout << i + 1 << ":\t" << this->selected_population[i].getStringRepresentationOfGenome() << "\tu="
                      << u;
            if (is_selected) {
                std::cout << "<" << this->recombination_probability << " he participates";
            }
            std::cout << "\n";
        }
    }
    if (show_steps) {
        std::cout << "\n";
    }
}

void GeneticAlgorithm::crossGenomesOfTwoChromosomes(Chromosome &chromosome1, Chromosome &chromosome2,
                                                    const size_t point_of_crossing) {
    // given 2 parents chromosomes we will create 2 children which will replace the parents in the next population
    // we will randomly select a point of crossing from their genome and reverse their prefixes of genomes
    // the prefix will start from gene 0 and end before the point of crossing chosen
    std::vector<bool> chromosome1_genome = chromosome1.getGenome();
    std::vector<bool> chromosome2_genome = chromosome2.getGenome();
    bool temporary_gene;
    for (size_t i = 0; i < point_of_crossing; ++i) {
        temporary_gene = chromosome1_genome[i];
        chromosome1_genome[i] = chromosome2_genome[i];
        chromosome2_genome[i] = temporary_gene;
    }
    chromosome1.setGenome(chromosome1_genome);
    chromosome2.setGenome(chromosome2_genome);
}

void GeneticAlgorithm::crossGenomesOfTwoChromosomesWithTwoPointsOfCrossingSwitchMiddle(Chromosome &chromosome1,
                                                                                       Chromosome &chromosome2,
                                                                                       const size_t point_of_crossing_left,
                                                                                       const size_t point_of_crossing_right) {
    // given 2 parents chromosomes we will create 2 children which will replace the parents in the next population
    // given 2 points of crossing switch the genes in this interval left-right
    std::vector<bool> chromosome1_genome = chromosome1.getGenome();
    std::vector<bool> chromosome2_genome = chromosome2.getGenome();
    bool temporary_gene;
    // the point of crossing left isn't necessary < than point of crossing right
    size_t left = point_of_crossing_left < point_of_crossing_right ? point_of_crossing_left : point_of_crossing_right;
    size_t right = point_of_crossing_left < point_of_crossing_right ? point_of_crossing_right : point_of_crossing_left;
    for (size_t i = left + 1; i < right; ++i) {
        temporary_gene = chromosome1_genome[i];
        chromosome1_genome[i] = chromosome2_genome[i];
        chromosome2_genome[i] = temporary_gene;
    }
    chromosome1.setGenome(chromosome1_genome);
    chromosome2.setGenome(chromosome2_genome);
}

void GeneticAlgorithm::crossPopulation() {
    // if the algorithm chosen an odd number of genes for crossing, then the last one will remain unchanged
    // we'll group the genes 2 by 2
    size_t number_of_groups = this->population_selected_for_crossing.size() / 2;
    size_t chromosome_index1, chromosome_index2;
    size_t point_of_crossing, point_of_crossing_left, point_of_crossing_right;
    for (size_t group_index = 0; group_index < number_of_groups; ++group_index) {
        chromosome_index1 = this->population_selected_for_crossing[2 * group_index];
        chromosome_index2 = this->population_selected_for_crossing[2 * group_index + 1];
        if (this->points_of_crossing == 1) {
            point_of_crossing = RandomGenerator::getRandomInt(0, this->population[0].getNumberOfGenes() - 1);
            if (show_steps) {
                std::cout << "Recombination between chromosome " << 2 * group_index << " and chromosome "
                          << 2 * group_index + 1 << ":\n";
                std::cout << this->selected_population[chromosome_index1].getStringRepresentationOfGenome() << " "
                          << this->selected_population[chromosome_index2].getStringRepresentationOfGenome()
                          << " point of crossing " << point_of_crossing << "\n";
            }
            crossGenomesOfTwoChromosomes(this->selected_population[chromosome_index1],
                                         this->selected_population[chromosome_index2], point_of_crossing);
        } else if (this->points_of_crossing == 2) {
            point_of_crossing_left = RandomGenerator::getRandomInt(0, this->population[0].getNumberOfGenes() - 1);
            point_of_crossing_right = RandomGenerator::getRandomInt(0, this->population[0].getNumberOfGenes() - 1);
            if (show_steps) {
                std::cout << "Recombination between chromosome " << 2 * group_index << " and chromosome "
                          << 2 * group_index + 1 << ":\n";
                std::cout << this->selected_population[chromosome_index1].getStringRepresentationOfGenome() << " "
                          << this->selected_population[chromosome_index2].getStringRepresentationOfGenome()
                          << " left point of crossing " << point_of_crossing_left
                          << " right point of crossing " << point_of_crossing_right << "\n";
            }
            crossGenomesOfTwoChromosomesWithTwoPointsOfCrossingSwitchMiddle(
                    this->selected_population[chromosome_index1],
                    this->selected_population[chromosome_index2],
                    point_of_crossing_left, point_of_crossing_right);
        }
        if (show_steps) {
            std::cout << "Result " << this->selected_population[chromosome_index1].getStringRepresentationOfGenome()
                      << " " << this->selected_population[chromosome_index2].getStringRepresentationOfGenome()
                      << "\n";
        }
    }
}

void GeneticAlgorithm::mutatePopulation() {
    // in the mutation step each gene from the genome of each chromosome from the proposed next generation can mutate
    // and change its value
    std::vector<bool> chromosome_genome;
    size_t number_of_genes = this->selected_population[0].getNumberOfGenes();
    long double u;
    bool chromosome_was_mutated;
    bool population_was_mutated{false};
    if (show_steps) {
        std::cout << "The following chromosomes have mutated:\n";
    }
    size_t chromosome_index = this->keep_elitist ? 1 : 0;
    for (; chromosome_index < this->population_size; ++chromosome_index) {
        chromosome_was_mutated = false;
        chromosome_genome = this->selected_population[chromosome_index].getGenome();
        for (size_t i = 0; i < number_of_genes; ++i) {
            u = RandomGenerator::getRandomNumberBetweenZeroAndOne();
            if (u < this->mutation_probability) {
                // not all genes in our representation have the same importance in the chromosome value, so we should
                // choose a different probability for each gene, so the significant genes have a lower probability
                // to mutate and create a huge difference between the original one and the mutated one
                chromosome_genome[i] = !chromosome_genome[i];
                chromosome_was_mutated = true;
                population_was_mutated = true;
            }
        }
        if (chromosome_was_mutated) {
            this->selected_population[chromosome_index].setGenome(chromosome_genome);
            if (show_steps) {
                std::cout << chromosome_index + 1 << "\n";
            }
        }
    }
    if (show_steps && !population_was_mutated) {
        std::cout << "No chromosome mutated\n";
    }
}

void GeneticAlgorithm::saveNextGeneration() {
    // remove the old population from the memory
    this->population.clear();
    // copy the computed generation as the new current population
    for (const auto &chromosome: this->selected_population) {
        this->population.emplace_back(Chromosome(chromosome));
    }
    // clear auxiliary space used by algorithm
    this->selection_probability.clear();
    this->selection_probability_interval.clear();
    this->selected_population.clear();
    this->population_selected_for_crossing.clear();
}

void GeneticAlgorithm::printInitialPopulation() const {
    std::cout << "Initial population:\n";
    std::streamsize ss = std::cout.precision();
    for (size_t i = 0; i < this->population_size; ++i) {
        std::cout << std::setw(ceil(log10(this->population_size))) << i + 1 << ": "
                  << this->population[i].getStringRepresentationOfGenome() << " f="
                  << this->fitnessFunction(this->population[i]) << '\n';
    }
}

void GeneticAlgorithm::printSelectionProbability() const {
    std::cout << "\nSelection probabilities:\n";
    for (size_t i = 0; i < this->population_size; ++i) {
        std::cout << "chromosome " << std::setw(ceil(log10(this->population_size))) << i + 1 << " probability "
                  << this->selection_probability[i] << '\n';
    }
}

void GeneticAlgorithm::printSelectionProbabilityInterval() const {
    std::cout << "\nSelection probabilities interval\n";
    for (const auto &value: this->selection_probability_interval) {
        std::cout << value << " ";
    }
}

void GeneticAlgorithm::printSelectedPopulationInfo() const {
    std::cout << "\nAfter selection\n";
    std::streamsize ss = std::cout.precision();
    for (size_t i = 0; i < this->population_size; ++i) {
        std::cout << std::setw(ceil(log10(this->population_size))) << i << ": "
                  << this->selected_population[i].getStringRepresentationOfGenome() << " f= "
                  << this->fitnessFunction(this->selected_population[i]) << '\n';
    }
}

void GeneticAlgorithm::printPopulationInfoAfterCrossing() const {
    std::cout << "After crossing\n";
    std::streamsize ss = std::cout.precision();
    for (size_t i = 0; i < this->population_size; ++i) {
        std::cout << std::setw(ceil(log10(this->population_size))) << i + 1 << ": "
                  << this->selected_population[i].getStringRepresentationOfGenome() << " f="
                  << this->fitnessFunction(this->selected_population[i]) << '\n';
    }
}

void GeneticAlgorithm::printPopulationInfoAfterMutation() const {
    std::cout << "\nAfter mutation\n";
    std::streamsize ss = std::cout.precision();
    for (size_t i = 0; i < this->population_size; ++i) {
        std::cout << std::setw(ceil(log10(this->population_size))) << i + 1 << ": "
                  << this->selected_population[i].getStringRepresentationOfGenome() << " f= "
                  << this->fitnessFunction(this->selected_population[i]) << '\n';
    }
}
