#ifndef _EVALUATION_MANAGER_HPP_INCLUDED_SOJIWORIASFOIJ4E98YASFKHJAWF9EYUH8RRRUHSAFILUHTI9YUHSAFHUUI
#define _EVALUATION_MANAGER_HPP_INCLUDED_SOJIWORIASFOIJ4E98YASFKHJAWF9EYUH8RRRUHSAFILUHTI9YUHSAFHUUI

#include <chromosome.hpp>

namespace ga
{

    template<typename Function_Argument_Type>
    struct evaluator
    {
        //caches the arguement type to save runtime
        // e.g.
        //      chromosome<feng::matrix<uint64_t>> ch;
        //      symmetric_matrix_to_real_complex_matrix_translator tr( 0.1, 0.3 );
        //      ...func
        //      evaluator<feng::matrix<std::complex<double>>> ev;
        // 
        //      auto fit = ev( ch, tr, func );
        //
        Function_Argument_Type arg;

        template<typename Chromosome_Dominance_Type, typename Fitness_Type, 
                 typename Chromosome_Translator, typename Function_Type>
        void
        operator()(const chromosome<Chromsome_Dominance_Type, Fitness_Type>& ch, 
                   const Chromosome_Translator& tr,
                   const Function_Type& f )
        {
            if ( ch.modification_after_evaluation_flag )
            {
                tr( ch.chrom, arg );
                ch.fit = f(arg);
                ch.modification_after_evaluation_flag = false; 
            }
        }

        template< typename Chromosome_Itor, typename Chromosome_Translator, typename Function_Type >
        void operator()( Chromosome_Itor ch_begin, Chromosome_Itor ch_end,
                         const Chromosome_Translator& tr, 
                         const Function_Type& f )
        {
            while ( ch_begin != ch_end )
            {
                operator()( *ch_begin++, tr, f );
            }
        }





    };//struct evaluator

    struct evaluator_without_cache
    {
            
    };

}//namespace ga

#endif//_EVALUATION_MANAGER_HPP_INCLUDED_SOJIWORIASFOIJ4E98YASFKHJAWF9EYUH8RRRUHSAFILUHTI9YUHSAFHUUI

