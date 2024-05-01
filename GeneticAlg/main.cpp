#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <random>
#include <numeric>

constexpr int POPULATION_SIZE = 1000;
constexpr int MAX_GENERATIONS = 10000;
constexpr double MUTATION_RATE = 0.7;
constexpr double MIN_VAL = 0.0;
constexpr double MAX_VAL = 1.0;

// Use random number engine for better random number generation
std::random_device rd;
std::mt19937 gen(rd());

double randVal(double min_val_, double max_val_)
{
    std::uniform_real_distribution<double> dis(min_val_, max_val_);
    return dis(gen);
}

// The problem to optimize
double foo(const double & x, const double & y, const double & z)
{
    return fabs(6 * pow(x, 3) + 9 * pow(y, 2) + 90 * z - 25);
}

// 2. Fitness Evaluation - Evaluate the fitness of each solution in the population based on some predefined criteria
double fitnessEval(const double & x, const double & y, const double & z)
{
    double ans = foo(x, y, z);
    return (ans == 0) ? INT_MAX : fabs(1 / ans);
}

// A potential solution
struct Solution
{
    double x, y, z;
	double fitness;

    Solution() = default;

    Solution(const double & x_, const double & y_, const double & z_)
        : x(x_)
        , y(y_)
        , z(z_)
        , fitness(fitnessEval(x_, y_, z_)) {}

    bool operator>(const Solution & other) const
    {
        return fitness > other.fitness;
    }

    friend std::ostream & operator<<(std::ostream & os, const Solution & obj)
    {
        os << "x = " << obj.x << ", "
           << "y = " << obj.y << ", "
           << "z = " << obj.z << ", "
           << "Fitness = " << obj.fitness;
        return os;
    }
};

// 1. Initialization - Generate an initial population of potential solutions randomly
void initPopulation(std::vector<Solution> & population)
{
	for (int i = 0; i < POPULATION_SIZE; ++i)
	{
        population[i] = Solution(randVal(MIN_VAL, MAX_VAL), randVal(MIN_VAL, MAX_VAL), randVal(MIN_VAL, MAX_VAL));
	}
}

// 3. Selection - Select individuals from the current population based on their fitness, favoring better individuals
void selection(std::vector<Solution> & population)
{
    population.resize((size_t)(std::ceil(population.size() * 0.6)));
}

// 4. Recombination (Crossover) - Create new solutions by combining genetic material from selected individuals
Solution crossover(const Solution & parent1, const Solution & parent2)
{
    Solution child((parent1.x + parent2.x) / 2.0, (parent1.y + parent2.y) / 2.0, (parent1.z + parent2.z) / 2.0);
    return child;
}

// 5. Mutation - Introduce random changes in the new solutions to maintain genetic diversity and explore new regions of the solution space
void mutate(Solution & child)
{
    if (randVal(0, 1) < MUTATION_RATE)
    {
        child.x *= randVal(0.99, 1.01);
        child.y *= randVal(0.99, 1.01);
        child.z *= randVal(0.99, 1.01);
        child.fitness = fitnessEval(child.x, child.y, child.z);
    }
}

std::vector<Solution> generation(std::vector<Solution> & population)
{
    selection(population);
    std::vector<Solution> newPopulation;

    std::uniform_int_distribution<int> dis(0, population.size() - 1);
    for (int i = 0; i < population.size(); i += 2)
    {
        Solution child = crossover(population[dis(gen)], population[dis(gen)]);
        mutate(child);
        newPopulation.push_back(child); // 6. Replacement - Replace the old population with the new population of solutions
    }

    return newPopulation;
}

void geneticAlgorithm()
{
    std::vector<Solution> population(POPULATION_SIZE);

    initPopulation(population);
    std::sort(population.begin(), population.end(), std::greater<Solution>());
    std::cout << "Generation " << 0 << " best solution: " << std::endl << population[0] << std::endl;

    for (int i = 1; i <= MAX_GENERATIONS; ++i)
    {
        population = generation(population);
        std::sort(population.begin(), population.end(), std::greater<Solution>());
        std::cout << "Generation " << i << " best solution: " << std::endl << population[0] << std::endl;

        // 7. Termination - Check if termination conditions are met (satisfactory solution found or maximum number of generations reached)
        if (population[0].fitness > 999)
        {
            break;
        }
    }

    std::cout << "foo(x, y, z) = " << foo(population[0].x, population[0].y, population[0].z);

    int a;
    std::cin >> a;
}

int main()
{
    geneticAlgorithm();
    return 0;
}
