#ifndef _GENETIC_ALGORITHM_HPP_INCLUDED_SFOJ3W94889YKLSFDJHXVCMSAJKL4389SFDJKSFAJKLSFAFIOOSFIDWREOFF
#define _GENETIC_ALGORITHM_HPP_INCLUDED_SFOJ3W94889YKLSFDJHXVCMSAJKL4389SFDJKSFAJKLSFAFIOOSFIDWREOFF

#include <memory>
#include <vector>
#include <cstddef>

#include <singleton.hpp>

namespace ga
{
    template<
                class Chromosome_Representation_Type,
                class Fitness_Type,
                class Fitness_Comparision_Type
            >
    struct chromosome
    {
        std::shared_ptr<Chromosome_Representation_Type>  ch;
        Fitness_Type                                     fit;

        friend bool operator < ( const chromosome& lhs, const chromosome& rhs )
        {
            return Fitness_Comparision_Type()(lhs, rhs); 
        }
    };//struct chromosome

    template<
                class Chromosome_Representation_Type,
                class Problem_Representation_Type,
                class Fitness_Type,
                class Fitness_Comparision_Type,
                class Chromosome_Mutation_Method,
                class Chromosome_Xover_Method,
                class Chromosome_Random_Generation_Method,
                class Evaluation_Method,
                class Runtime_Manager_Method
            >
    struct genetic_algorithm
    {
        typedef chromosome<Chromosome_Representation_Type, Fitness_Type, Fitness_Comparision_Type> chromosome_type;
        typedef std::vector<chromosome_type> population_type; 
        
        population_type current_population;

        Problem_Representation_Type const 
        operator()( const std::size_t runtime_n, const std::size_t individual_n ) 
        {
            





                


        }

    };//struct genetic_algorithm 

#include "genetic_algorithm.tcc"

}//namespace ga

#endif//_GENETIC_ALGORITHM_HPP_INCLUDED_SFOJ3W94889YKLSFDJHXVCMSAJKL4389SFDJKSFAJKLSFAFIOOSFIDWREOFF

