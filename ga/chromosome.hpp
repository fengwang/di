#ifndef _CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU
#define _CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU

namespace ga
{
    // example:
    // std::shared_ptr<chromosome<feng::matrix<uint64>>> sfoijeoijsfda;

    template < typename Chromosome_Dominance_Type, typename Fitness_Type = double >
    struct chromosome
    {
        typedef Chromosome_Dominance_Type value_type;

        Chromosome_Dominance_Type chrom; //a kind of data structure representing the value of the chromosome
        Fitness_Type fit; //the greater the better, no normalization required

        bool modification_after_evaluation_flag; // true if this chromosome has been modified after last evalution
    };//struct chromosome

    template< typename Chromosome_Dominance_Type, typename Fitness_Type >
    bool operator < ( const chromosome<Chromosome_Dominance_Type, Fitness_Type>& lhs,
                      const chromosome<Chromosome_Dominance_Type, Fitness_Type>& rhs )
    {
        return lhs.fit < rhs.fit;
    }

}//namespace ga

#endif//_CHROMOSOME_HPP_INCLUDED_SOFIJW98RUASFOIUHEW8RUFHASFDLKUHASFIOUHVSDKJHBVCJHBSAKFIUYHASF8IGYU

