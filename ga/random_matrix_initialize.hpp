#ifndef _RANDOM_MATRIX_INITIALIZE_HPP_INCLUDED_FSOIJ498LGKJS8IJH4KLNFKLJAKLJHSUI9498YFSDJKLFSU84FSFF
#define _RANDOM_MATRIX_INITIALIZE_HPP_INCLUDED_FSOIJ498LGKJS8IJH4KLNFKLJAKLJHSUI9498YFSDJKLFSU84FSFF

#include <matrix.hpp>
#include <vg.hpp>

#include <cstdint>
#include <limits>

namespace ga
{
   
    struct random_initialize
    {
        typedef std::uint64_t uint_type;

        uint_type operator()() const 
        {
            auto& vg = feng::singleton<vg::variate_generator<long double>>::instance();

            long double var = vg() * ( std::numeric_limits<uint_type>::max() - std::numeric_limits<uint_type>::min() ) + std::numeric_limits<uint_type>::min();

            return var;

        }
    };//struct random_initialize 

    // usage:
    //          std::size_t n = 16;
    //          auto m = symmetric_random_matrix_initialize()(n);
    struct symmetric_random_matrix_initialize
    {
        typedef std::uint64_t uint_type;
        
        feng::matrix<uint_type> const 
        operator()( const std::size_t n ) const 
        {
            feng::matrix<uint_type> m(n, n);

            for ( std::size_t i = 1; i != n; ++i )
            {
                std::generate( m.upper_diag_begin(i), m.upper_diag_end(i), random_initialize() );
                std::copy( m.upper_diag_begin(i), m.upper_diag_end(i), m.lower_diag_begin(i) );
            }

            return m;
        
        }

    };//struct random_matrix_initialize 


}//namespace ga

#endif//_RANDOM_MATRIX_INITIALIZE_HPP_INCLUDED_FSOIJ498LGKJS8IJH4KLNFKLJAKLJHSUI9498YFSDJKLFSU84FSFF

