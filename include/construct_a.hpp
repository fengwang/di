#ifndef _CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF
#define _CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

#include <array.hpp>
#include <lattice.hpp>
#include <matrix.hpp>

#include <vector>
#include <cstddef>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <complex>
#include <iterator>

#include <iostream>

namespace feng
{
template<typename T = double>
struct construct_a
{
    typedef std::size_t                     size_type;
    typedef T                               value_type;
    typedef std::complex<value_type>        complex_type;
    typedef array<value_type>               array_type;
    typedef matrix<value_type>              matrix_type;
    typedef matrix<complex_type>            complex_matrix_type;
    typedef matrix<array_type>              array_matrix_type;
    typedef construct_a                     self_type;
    typedef lattice<value_type>             lattice_type;
    typedef std::vector<value_type>         vector_type;
    typedef std::vector<complex_type>       complex_vector_type;
    typedef std::vector<array_type>         array_vector_type;

    array_type      kt;
    array_type      zone;
    array_type      gx;
    value_type      v0;
    lattice_type    la;

    construct_a( const self_type& ) = default;
    construct_a( self_type && ) = default;
    self_type& operator = ( const self_type& ) = default;
    self_type& operator = ( self_type && ) = default;

    construct_a( const array_type& kt_, const array_type& zone_,
                 const array_type& gx_, const value_type  v0_,
                 const lattice_type& la_ ) :
        kt( kt_ ), zone( zone_ ), gx( gx_ ), v0( v0_ ), la( la_ )
    {}

    construct_a()
    {
        kt      = array_type {0, 0, 0};
        zone    = array_type {1, 1, 0};
        gx      = array_type {0, 0, 1};
        v0      = 200.0;
        la      = lattice_type {5.43, 5.43, 5.43, 1.5707963268, 1.5707963268, 1.5707963268 };
    }

    const array_type make_gyvec() const
    {
        return cross_product( zone, gx );
    }

    const matrix_type make_matrix() const
    {
        const auto a = la.a;
        const auto b = la.b;
        const auto c = la.c;
        const auto ca = std::cos( la.alpha );
        const auto cb = std::cos( la.beta );
        const auto cg = std::cos( la.gamma );
        const auto sg = std::sin( la.gamma );
        matrix_type A( 3, 3 );
        A[0][0] = a;
        A[0][1] = 0;
        A[0][2] = 0;
        A[1][0] = b * cg;
        A[1][1] = b * sg;
        A[1][2] = 0;
        A[2][0] = c * cb;
        A[2][1] = c * ( ca - cb * cg ) / sg;
        A[2][2] = c * std::sqrt( T( 1 ) - ca * ca - cb * cb - cg * cg + T( 2 ) * ca * cb * cg ) / sg;
        return A;
    }

    const matrix_type make_gaussian_electron( const array_type& s, const value_type v ) const
    {
        matrix_type ans{ 1, s.size(), s.begin(), s.end() };
        auto const factor = v / value_type( 511 ) + value_type( 1 );
        std::for_each( ans.begin(), ans.end(), [factor]( value_type & s_ )
        {
            const value_type A[] = {0.3626, 0.9737, 2.7209, 1.7660};
            const value_type B[] = {0.4281, 3.5770, 19.3905, 64.3334};
            value_type tmp[4];
            std::fill( tmp, tmp + 4, value_type() );
            auto const ss = s_ * s_;
            std::transform( B, B + 4, tmp, [ss]( value_type b )
            {
                return std::exp( -b * ss );
            } );
            s_ = factor * std::inner_product( A, A + 4, tmp, value_type() );
        }
                     );
        return repmat( ans, 8, 1 );
    }

    const matrix_type make_gaussian_electron( const matrix_type& s, const value_type v ) const
    {
        auto ans = s;
        auto const factor = v / value_type( 511 ) + value_type( 1 );
        std::for_each( ans.begin(), ans.end(), [factor]( value_type & s_ )
        {
            const value_type A[] = {0.3626, 0.9737, 2.7209, 1.7660};
            const value_type B[] = {0.4281, 3.5770, 19.3905, 64.3334};
            value_type tmp[4];
            std::fill( tmp, tmp + 4, value_type() );
            auto const ss = s_ * s_;
            std::transform( B, B + 4, tmp, [ss]( value_type b )
            {
                return std::exp( -b * ss );
            } );
            s_ = factor * std::inner_product( A, A + 4, tmp, value_type() );
        }
                     );
        return repmat( ans, 8, 1 );
    }

    const complex_matrix_type make_ug( const matrix_type& G, const matrix_type& A, const matrix_type& D ) const
    {
        assert( G.col() == 3 );
        assert( A.col() == 3 );
        assert( D.col() == 1 );
        assert( A.row() == D.row() );
        auto const M = make_matrix();
        auto const S = G * ( M.inverse() );
        matrix_type s( 1, S.row() );

        for ( size_type i = 0; i < S.row(); ++ i )
        {
            s[0][i] = value_type( 0.5 ) * std::sqrt( std::inner_product( S.row_begin( i ), S.row_end( i ), S.row_begin( i ), value_type( 0 ) ) );
        }

        auto const piomega =  3.141592553590 * feng::inner_product( array_type( M[0][0], M[1][0], M[2][0] ), feng::cross_product( array_type( M[0][1], M[1][1], M[2][1] ), array_type( M[0][2], M[1][2], M[2][2] ) ) );
        auto const atomcellfacte = make_gaussian_electron( s, v0 );
        const complex_matrix_type dwss = D * feng::pow( s, value_type( 2 ) );
        const complex_matrix_type piag = A * G.transpose();
        auto fact = feng::exp( - dwss - piag * complex_type( 0, 6.2831853071796 ) );
        std::transform( fact.begin(), fact.end(), atomcellfacte.begin(), fact.begin(), [piomega]( const complex_type f,  const value_type a )
        {
            return f * a / piomega;
        } );
        complex_matrix_type Ug( fact.col(), 1 );

        for ( size_type i = 0; i < fact.col(); ++i )
        {
            Ug[i][0] = std::accumulate( fact.col_begin( i ), fact.col_end( i ), complex_type() );
            if ( std::abs(Ug[i][0].real()) < 1.0e-5 ) Ug[i][0].real(0);
            if ( std::abs(Ug[i][0].imag()) < 1.0e-5 ) Ug[i][0].imag(0);
        }

        return Ug;
    }

    const matrix_type make_atomic_positions() const
    {
        const value_type pos_x[] = { -0.125, -0.125, 0.375, 0.375, 0.125, 0.125, 0.625, 0.625};
        const value_type pos_y[] = { -0.125, 0.375, -0.125, 0.375, 0.125, 0.625, 0.125, 0.625};
        const value_type pos_z[] = { -0.25,  0.25,   0.25, -0.25,   0.0,   0.5,   0.5, 0.0};
        matrix_type m( 8, 3 );
        std::copy( pos_x, pos_x + 8, m.col_begin( 0 ) );
        std::copy( pos_y, pos_y + 8, m.col_begin( 1 ) );
        std::copy( pos_z, pos_z + 8, m.col_begin( 2 ) );
        std::for_each( m.col_begin( 2 ), m.col_end( 2 ), []( value_type & v )
        {
            v += 0.125;
        } );
        return m;
    }

    struct complex_size
    {
        complex_type c;
        size_type    s;
        complex_size( const complex_type& c_ = complex_type( 0, 0 ), const size_type s_ = size_type( 0 ) ) : c( c_ ), s( s_ ) {}
    };//complex_size

    const array_vector_type make_beam_vector() const
    {
        auto const mx = 7;// make_beamx_max();
        auto const my = 3;// make_beamy_max();
        auto const gy = make_gyvec();
        array_vector_type beams {1};

        for ( int y = 0; y <= my; ++y )
            for ( int x = 0; x <= mx; ++x )
            {
                if ( x )
                {
                    array_type const b0 = ( x * gx ) + ( y * gy );
                    beams.push_back( b0 );
                    beams.push_back( -b0 );
                }

                if ( y )
                {
                    array_type const b0 = - ( x * gx ) + ( y * gy );
                    beams.push_back( b0 );
                    beams.push_back( -b0 );
                }
            }

        return beams;
    }

    const array_vector_type make_gd( const array_vector_type& beams ) const
    {
        array_vector_type gd {beams.size() };
        auto const nm = beams.size() >> 1;

        for ( size_type i = 1; i <= nm; ++i )
        {
            gd[nm + i] = beams[i + i - 1];
            gd[nm - i] = beams[i + i];
        }

        return gd;
    }

    const array_matrix_type make_gm( const array_vector_type& gd ) const
    {
        auto const N = gd.size();
        array_matrix_type gr( N, N );
        array_matrix_type gc( N, N );

        for ( size_type i = 0; i != N; ++i )
        {
            std::fill( gr.row_begin( i ), gr.row_end( i ), gd[i] );
            std::fill( gc.col_begin( i ), gc.col_end( i ), gd[i] );
        }

        auto gm = gr - gc;
        //std::copy( gd.begin(), gd.end(), gm.diag_begin() );
        std::fill( gm.diag_begin(), gm.diag_end(), gm[0][1] );
        return gm;
    }

    const array_vector_type make_unique_beams( const array_matrix_type& Gm ) const
    {
        array_vector_type gm( Gm.begin(), Gm.end() );
        std::sort( gm.begin(), gm.end() );
        gm.resize( std::distance( gm.begin(), std::unique( gm.begin(), gm.end() ) ) );
        return gm;
    }

    const complex_matrix_type make_ug( const array_vector_type& g ) const
    {
        matrix_type bArray( g.size(), 3 );
        for ( size_type i = 0; i < g.size(); ++i )
            std::copy( g[i].begin(), g[i].end(), bArray.row_begin(i) );

        auto const A = make_atomic_positions();
        auto const D = matrix_type {8, 1, 0.43};
        return make_ug( bArray, A, D );
    }

    const complex_matrix_type make_ug_se( const matrix_type& bArray ) const
    {
        auto const A = make_atomic_positions();
        auto const D = matrix_type {8, 1, 0.43};
        return make_ug( bArray, A, D );
    }

    const complex_matrix_type operator()() const
    {
        auto const Gm = make_gm(make_gd(make_beam_vector()));
        auto const Ub = make_unique_beams(Gm);
        auto const Ug = make_ug(Ub);
        auto const R  = make_matrix().inverse();
        complex_matrix_type A( Ug.size(), Ug.size() );
        for ( size_type h = 0; h != Ug.size(); ++h )
        {
            auto const vr = Gm[h][h] * R;
            A[h][h] = std::inner_product( vr.begin(), vr.end(), vr.begin(), value_type(0) );

            for ( size_type k = 0; k != Ug.size(); ++k )
            {
                if ( h == k ) continue;
                auto itor = std::lower_bound( Ub.begin(), Ub.end(), vr );
                assert( itor != Ub.end() );
                A[h][k] = *( Ug.begin() + std::distance( Ub.begin(), std::lower_bound( Ub.begin(), Ub.end(), vr ) ) );


            }
        }


    }

};//sturct construct_a

}//namespace feng

#endif//_CONSTRUCT_A_HPP_INCLUDED_SFOIDDDDFFSF8I3498FDAKJSSSSDFSFFFFFFFFFFFFFFFFFFFFFFFFF

