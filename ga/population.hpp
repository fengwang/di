#ifndef _POPULATION_HPP_INCLUDED_SDOIFHJ398TYUHAFLKHJASOFDIUH34E8Y7GHASFDIH48T37YSFIUDH4387YASFIUHAS
#define _POPULATION_HPP_INCLUDED_SDOIFHJ398TYUHAFLKHJASOFDIUH34E8Y7GHASFDIH48T37YSFIUDH4387YASFIUHAS

#include <vector>
#include <memory>

namespace ga
{

    // the data structure should looks like
    // Container<Chromosome1, Chromosome2, ..., ChromosomeN>
    //
    // an example:
    //              typedef chromosome<feng::matrix<double>> chromosome_type;
    //              typedef population<std::vector, chromosome_type> population_type;
    //
    template<class Chromosome_Type>
    struct population
    {
        typedef std::shared_ptr<Chromosome_Type> chromosome_pointer_type;

        typedef std::vector<chromosome_pointer_type> chromosome_pool_type;

        chromosome_pool_type pool;      //the individuals living in the current pool
 
        chromosome_pointer_type elite;  //the chromosome fit best
    };//struct population

}//namespace ga

#endif//_POPULATION_HPP_INCLUDED_SDOIFHJ398TYUHAFLKHJASOFDIUH34E8Y7GHASFDIH48T37YSFIUDH4387YASFIUHAS

