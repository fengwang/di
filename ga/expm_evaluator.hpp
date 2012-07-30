#ifndef _EXPM_EVALUATOR_HPP_INCLUDED_SDFOIU4398YUASFKJH4978YASFIUHYU7HASFIUHT9YU8SFUISFUIWREUI998473
#define _EXPM_EVALUATOR_HPP_INCLUDED_SDFOIU4398YUASFKJH4978YASFIUHYU7HASFIUHT9YU8SFUISFUIWREUI998473

#include <simulation.hpp>

namespace feng
{

template<typename T=double>
struct expm_evaluator
{
    typedef T value_type;
    typedef feng::construct_a<value_type>                host_type;
    typedef typename host_type::complex_type             complex_type;
    typedef typename host_type::matrix_type              matrix_type;
    typedef typename host_type::complex_matrix_type      complex_matrix_type;
    typedef typename host_type::vector_type              vector_type;

    matrix_type offset;
    vector_type gxy2;
    complex_matrix_type A_real;
    matrix_type Is_real;
    //cache for acceleration
    complex_matrix_type A_guess;
    matrix_type Is_guess;
    complex_matrix_type A_expm_cache;

    expm_evaluator():   offset( host_type().make_gxy2_offset() ),
                        gxy2( host_type().make_gxy2() ),
                        A_real( host_type().make_new_a() )
    {
        Is_real.resize( A_real.row(), offset.col() );
        for ( size_t i = 0; i != offset.col(); ++i )
        //for ( size_t i = 0; i != offset.col() && i != A_real.row(); ++i )
        {
            auto I = host_type(). make_new_i_with_offset(A_real, 7.879, offset.col_begin(i));
            auto const mid = I.row() >> 1;
            std::copy( I.col_begin(mid), I.col_end(mid), Is_real.col_begin(i) );
        }    

    }

    value_type operator()( const complex_matrix_type& A )
    {
        assert( A.row() == A_real.row() );
        assert( A.col() == A_real.col() );
        A_guess = A;
        feng::for_each( A_guess.diag_begin(), A_guess.diag_end(), gxy2.begin(), []( complex_type&c, const value_type& v ){ c = complex_type(-v, 0); } );

        Is_guess.resize( A_guess.row(), offset.col() );

        std::size_t const n = A.row();

        for ( size_t i = 0; i != offset.col(); ++i )
        {
            host_type().make_new_i_with_offset( A_expm_cache, A_guess, 7.879, offset.col_begin(i), Is_guess.col_begin(i) );
        }    

        value_type diff = 0;

        feng::for_each( Is_guess.begin(), Is_guess.end(), Is_real.begin(), [&diff](value_type const ig, value_type const ir ) { auto const d=ig-ir; diff += d*d; } );

        return diff;
    }

};
}//namespace feng

#endif//_EXPM_EVALUATOR_HPP_INCLUDED_SDFOIU4398YUASFKJH4978YASFIUHYU7HASFIUHT9YU8SFUISFUIWREUI998473

