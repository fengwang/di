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
    template<template<class, class> class Container_Type, class Chromosome_Type>
    struct population
    {
        Container_Type<Chromosome_Type, std::allocator<Chromosome_Type>> individual_pool;




    };//struct population

}//namespace ga

#endif//_POPULATION_HPP_INCLUDED_SDOIFHJ398TYUHAFLKHJASOFDIUH34E8Y7GHASFDIH48T37YSFIUDH4387YASFIUHAS

