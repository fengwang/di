#include <simulation.hpp>

#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>

#include <ctime>

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
        for ( size_t i = 0; i != offset.col() && i != A_real.row(); ++i )
        {
            auto I = host_type(). make_new_i_with_offset(A_real, 7.879, offset.col_begin(i));
            auto const mid = I.row() >> 1;
            std::copy( I.col_begin(mid), I.col_end(mid), Is_real.col_begin(i) );
        }    

    }

    value_type operator()( const complex_matrix_type& A )
    {
        A_guess = A;
        feng::for_each( A_guess.diag_begin(), A_guess.diag_end(), gxy2.begin(), []( complex_type&c, const value_type& v ){ c = complex_type(-v, 0); } );

        Is_guess.resize( A_guess.row(), offset.col() );

        std::size_t const n = A.row();

        for ( size_t i = 0; (i != offset.col()) && (i != n); ++i )
        {
            host_type().make_new_i_with_offset( A_expm_cache, A_guess, 7.879, offset.col_begin(i), Is_guess.col_begin(i) );
        }    

        value_type diff = 0;

        feng::for_each( Is_guess.begin(), Is_guess.end(), Is_real.begin(), [&diff](value_type const ig, value_type const ir ) { auto const d=ig-ir; diff += d*d; } );

        return diff;
    }

};

int main()
{
#if 1
    auto A = feng::construct_a<double>().make_new_a();
    feng::for_each( A.begin(), A.end(), [](std::complex<double>&d){ if(std::abs(d.real())<1.8e-8) d.real(0); if(std::abs(d.imag())<1.0e-8) d.imag(0); } );

    auto const offset = feng::construct_a<double>().make_gxy2_offset();

    feng::construct_a<double> :: matrix_type Is( A.row(), offset.col() );


    std::cout << "\nA=\n" << A;

    for ( size_t i = 0; i != offset.col(); ++i )
    {
        auto I = feng::construct_a<double>(). make_new_i_with_offset(A,7.879, offset.col_begin(i));
        feng::for_each( I.begin(), I.end(), [](double&d){ if (std::abs(d)<=1.0e-8) d=0; } );
        auto const mid = I.row() >> 1;
        std::copy( I.col_begin(mid), I.col_end(mid), Is.col_begin(i) );
    }

    auto II = feng::pow( Is, 2 );

    std::cout << "\ntotal energy is " << std::accumulate( II.begin(), II.end(), double(0) );

    std::cout << "\nthe columns of current I is " << II.col();
    
    std::cout << "\n current I is \n" << Is;

#endif

#if 0
    expm_evaluator<double> ee;
    auto A_ = ee.A_real;
    std::fill( A_.diag_begin(), A_.diag_end(), 0.0 );
    std::cout << ee( A );
#endif

    return 0;









}

