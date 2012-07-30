#ifndef _RANDOM_MATRIX_INITIALIZE_HPP_INCLUDED_FSOIJ498LGKJS8IJH4KLNFKLJAKLJHSUI9498YFSDJKLFSU84FSFF
#define _RANDOM_MATRIX_INITIALIZE_HPP_INCLUDED_FSOIJ498LGKJS8IJH4KLNFKLJAKLJHSUI9498YFSDJKLFSU84FSFF

#include <matrix.hpp>
#include <vg.hpp>
#include <chromosome.hpp>

#include <cstdint>
#include <limits>
#include <memory>

namespace ga
{
   
    template<typename Uint_Type = std::uint64_t>
    struct random_initialize
    {
        typedef Uint_Type uint_type;

        uint_type operator()() const 
        {
            auto& vg = feng::singleton<vg::variate_generator<long double>>::instance();

            long double var = vg() * ( std::numeric_limits<uint_type>::max() - std::numeric_limits<uint_type>::min() ) + std::numeric_limits<uint_type>::min();

            return var;

        }
    };//struct random_initialize 

    struct symmetric_random_matrix_initialize
    {
        template<typename T, std::size_t D, typename A>
        void operator()(feng::matrix<T,D,A>& m ) const 
        {
            assert( m.row() == m.col() );

            for ( std::size_t i = 1; i != m.row(); ++i )
            {
                feng::for_each( m.upper_diag_begin(i), m.upper_diag_end(i), 
                                [](T& v){v = random_initialize<T>()();} ); 
            
                std::copy( m.upper_diag_begin(i), m.upper_diag_end(i), 
                           m.lower_diag_begin(i) );
            }
        }
    };

#if 0
    // usage:
    //          std::size_t n = 16;
    //          auto m = symmetric_random_matrix_initialize()(n);
    template<typename Uint_Type = std::uint64_t>
    struct symmetric_random_matrix_initialize
    {
        typedef Uint_Type uint_type;
        
        feng::matrix<uint_type> const 
        operator()( const std::size_t n ) const 
        {
            feng::matrix<uint_type> m(n, n);

            for ( std::size_t i = 1; i != n; ++i )
            {
                std::generate( m.upper_diag_begin(i), m.upper_diag_end(i), random_initialize<uint_type>() );
                std::copy( m.upper_diag_begin(i), m.upper_diag_end(i), m.lower_diag_begin(i) );
            }

            return m;
        }
        template<typename T, std::size_t D, typename A>
        void operator()(feng::matrix<T,D,A>& m ) const 
        {
            feng::for_each( m.begin(), m.end(), [](T& v){v = random_initialize<T>()();} ); 
        }

    };//struct random_matrix_initialize 

    // generate a random matrix for a chromosome as initial status
    template<typename Uint_Type = std::uint64_t, typename Fit_Type = double>
    struct chromosome_random_initialize
    {
        typedef feng::matrix<Uint_Type> matrix_type;
        typedef Fit_Type fit_type;
        typedef chromosome<matrix_type, fit_type> chromosome_type;

        std::shared_ptr<chromosome_type> 
        operator()(std::size_t n) const 
        {
            auto ans = std::make_shared<chromosome_type>();
            (*ans).chrom = symmetric_random_matrix_initialize<Uint_Type>()(n);
            (*ans).modification_after_evaluation_flag = true;
            return ans;
        }
    
    };//chromosomei_random_initialize 

#endif








}//namespace ga

#endif//_RANDOM_MATRIX_INITIALIZE_HPP_INCLUDED_FSOIJ498LGKJS8IJH4KLNFKLJAKLJHSUI9498YFSDJKLFSU84FSFF

