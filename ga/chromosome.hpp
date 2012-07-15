#ifndef _CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU
#define _CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU

#include <matrix.hpp>

#include <vector>

namespace ga
{
    template < typename Chromosome_Dominance_Type, typename Fitness_Type = double >
    struct chromosome
    {
        typedef Chromosome_Dominance_Type value_type;

        Chromosome_Dominance_Type chrom; //a kind of data structure representing the value of the chromosome
        Fitness_Type fit; //the greater the better, no normalization required
    };//struct chromosome

}//namespace ga

#endif//_CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU

