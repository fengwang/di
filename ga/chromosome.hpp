#ifndef _CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU
#define _CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU

#include <matrix.hpp>

#include <vector>

namespace ga
{
    template< typename Chromosome_Dominance_Type, 
              typename Fitness_Type = double, 
              typename Probability_type = double >
    struct chromosome
    {
        typedef Chromosome_Dominance_Type value_type;

        Chromosome_Dominance_Type chrom;
        Fitness_Type fit;
        Probability_type prob;
    };//struct chromosome 

}//namespace ga

#endif//_CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU

