#ifndef _MATRIX_BINARY_TO_REAL_HPP_INCLUDED_DSFOJI4398UASF98UYH34TY89OEF8HUAWIFUYH34ETY78HSFDKUHFSDU
#define _MATRIX_BINARY_TO_REAL_HPP_INCLUDED_DSFOJI4398UASF98UYH34TY89OEF8HUAWIFUYH34ETY78HSFDKUHFSDU

#include <matrix.hpp>
#include <complex>

#include <cassert>
#include <limits>

namespace ga
{
    //          decode a binary matrix to a real complex matrix 
    // usage:
    //          matrix<uint64_t> m1;
    //          matrix<std::complex<double>> m2;
    //          double low = -0.6;
    //          double hig = 0.6;
    //          binary_matrix_to_real_complex_matrix()( m1, m2, low, hig );
    struct binary_matrix_to_real_complex_matrix_symmetric
    {
        // convert a binary matrix to a real complex matrix, the Block wave potential matrix  Ug
        // the diagonal value should all be zero, and for symmetric crystals, the imag part is zero
        //template<typename Binary_Matrix, typename Real_Complex_Matrix, typename Real_Type=double>
        //void operator()( const Binary_Matrix& bm, Real_Complex_Matrix& rm, const Real_Type lower=0.0, const Real_Type upper=0.06 ) const
        template<typename T1, std::size_t D1, typename A1, typename T2, std::size_t D2, typename A2, typename Real_Type=double>
        void operator()( const feng::matrix<T1,D1,A1>& bm, feng::matrix<T2,D2,A2>& rm, const Real_Type lower=-0.06, const Real_Type upper=0.06 ) const
        {
            typedef T1 int_type;
            typedef T2 complex_type;

            assert( bm.row() == bm.col() ); 
            auto const n = bm.row();

            rm.resize( n, n );

            auto const f = [lower,upper]( int_type V ) -> Real_Type
            {
                Real_Type const Upper = std::numeric_limits<int_type>::max();
                Real_Type const Lower = std::numeric_limits<int_type>::min();
                Real_Type const Diff = Upper - Lower;
                Real_Type const diff = upper - lower;
                Real_Type const d = V - Lower;
                Real_Type const ans = d * diff / Diff + lower;
                return ans;
            };

            for ( std::size_t i = 1; i != n; ++i )
            {
                feng::for_each( bm.upper_diag_begin(i), bm.upper_diag_end(i), rm.upper_diag_begin(i), 
                                [f](int_type const v, complex_type& c) { c.real(f(v)); c.imag(0);} 
                              );
                std::copy( rm.upper_diag_begin(i), rm.upper_diag_end(i), rm.lower_diag_begin(i) );
            }

            std::fill( rm.diag_begin(), rm.diag_end(), complex_type(0,0) );
        }

        template<typename Binary_Matrix, typename Real_Type=double>
        feng::matrix<std::complex<Real_Type>> const 
        operator()( const Binary_Matrix& bm, const Real_Type lower=-0.06, const Real_Type upper=0.06 ) const 
        {
            feng::matrix<std::complex<Real_Type>> rm( bm.row(), bm.col() );
            binary_matrix_to_real_complex_matrix_symmetric()( bm, rm, lower, upper );
            return rm;
        }

    };//struct binary_matrix_to_real_complex_matrix 

#if 0
    template<typename Type = double>
    struct symmetric_binary_matrix_to_real_complex_matrix_translator
    {
        typedef Type value_type;
        value_type lower_boundary;
        value_type upper_boundary;

        template<typename Binary_Matrix, typename Real_Complex_Matrix>
        void operator()( const Binary_Matrix& bm, Real_Complex_Matrix& rm ) const
        {
            return binary_matrix_to_real_complex_matrix_symmetric()(bm, rm, lower_boundary, upper_boundary);
        }

        symmetric_binary_matrix_to_real_complex_matrix_translator( const value_type lower_boundary_=-0.06, const value_type upper_boundary_=0.06 ) :
            lower_boundary(lower_boundary_), upper_boundary(upper_boundary_)
        {}

    };
#endif
}//namespace ga

#endif//_MATRIX_BINARY_TO_REAL_HPP_INCLUDED_DSFOJI4398UASF98UYH34TY89OEF8HUAWIFUYH34ETY78HSFDKUHFSDU

